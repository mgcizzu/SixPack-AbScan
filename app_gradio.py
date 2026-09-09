#!/usr/bin/env python3
"""Gradio app for SixPack-AbScan."""

from __future__ import annotations

import atexit
import os
import shutil
import socket
import traceback
from collections.abc import Generator
from datetime import UTC, datetime
from pathlib import Path

import gradio as gr
import pandas as pd

from analysis_modes import (
    ANALYSIS_MODE_CHOICES,
    AnalysisMode,
    OutputPolicy,
    output_policy_for,
    parse_analysis_mode,
)
from fasta_input import (
    MAX_FASTA_SOURCE_BYTES,
    InputPreparationError,
    format_bytes,
    prepare_fasta,
)
from protected_catalogue import (
    ProtectedCatalogueError,
    get_configured_catalogue,
    get_configured_catalogue_choices,
    get_configured_catalogue_paths,
    scan_protected_catalogue_with_progress,
)
from sixpack_abscan import (
    DEFAULT_GENETIC_CODE_TABLE,
    GENETIC_CODE_TABLE_NAMES,
    genetic_code_table_name,
    normalize_epitope,
    read_epitope_table,
    scan_epitopes_with_progress,
    write_six_frame_fasta_with_progress,
)

RUNS_DIR = Path("runs")
_SESSION_RUN_DIRS: set[Path] = set()
NUCLEOTIDE_MODE = "Nucleotide FASTA (will be 6-frame translated automatically)"
PROTEIN_MODE = "Protein FASTA (precomputed proteome)"
FASTA_FILE_TYPES = [".fa", ".fasta", ".fna", ".faa", ".fas", ".gz"]
GENETIC_CODE_CHOICES = [
    (f"{table_id} — {name}", table_id)
    for table_id, name in GENETIC_CODE_TABLE_NAMES.items()
]
APP_CSS = """
.gradio-container {
    font-size: 18px;
}

.gradio-container h1 {
    font-size: 2.4rem;
}

.gradio-container h2,
.gradio-container h3 {
    font-size: 1.5rem;
}

.gradio-container label,
.gradio-container .prose,
.gradio-container .gr-markdown,
.gradio-container .gr-button,
.gradio-container input,
.gradio-container textarea,
.gradio-container table,
.gradio-container select,
.gradio-container .wrap,
.gradio-container .message {
    font-size: 1.05rem;
}
"""
APP_HEAD = """
<meta name="description" content="SixPack-AbScan predicts crossreactivity of monoclonal antibodies on non-target species by searching epitope sequences against nucleotide or protein FASTA inputs.">
<meta property="og:title" content="SixPack-AbScan">
<meta property="og:description" content="Predict crossreactivity of monoclonal antibodies on non-target species.">
<meta property="og:type" content="website">
<meta property="og:image" content="https://images.unsplash.com/photo-1707863081130-7048e715688f?auto=format&fit=crop&fm=jpg&ixid=M3wxMjA3fDB8MHxwaG90by1wYWdlfHx8fGVufDB8fHx8fA%3D%3D&q=60&w=1200">
<meta name="twitter:card" content="summary">
<meta name="twitter:title" content="SixPack-AbScan">
<meta name="twitter:description" content="Predict crossreactivity of monoclonal antibodies on non-target species.">
<meta name="twitter:image" content="https://images.unsplash.com/photo-1707863081130-7048e715688f?auto=format&fit=crop&fm=jpg&ixid=M3wxMjA3fDB8MHxwaG90by1wYWdlfHx8fGVufDB8fHx8fA%3D%3D&q=60&w=1200">
"""


def _find_free_port(start: int = 7860, end: int = 7870) -> int:
    for port in range(start, end + 1):
        with socket.socket(socket.AF_INET, socket.SOCK_STREAM) as sock:
            sock.setsockopt(socket.SOL_SOCKET, socket.SO_REUSEADDR, 1)
            try:
                sock.bind(("0.0.0.0", port))
            except OSError:
                continue
            return port
    raise RuntimeError(f"No free port found in range {start}-{end}.")


def _blocked_catalogue_paths() -> list[str]:
    """Deny Gradio file-route access to every configured private catalogue."""

    return get_configured_catalogue_paths()


def _cleanup_previous_runs() -> None:
    if not RUNS_DIR.exists():
        return
    for entry in RUNS_DIR.iterdir():
        if entry.is_dir() and entry.name.startswith("run_"):
            shutil.rmtree(entry, ignore_errors=True)


def _cleanup_session_runs() -> None:
    for run_dir in list(_SESSION_RUN_DIRS):
        shutil.rmtree(run_dir, ignore_errors=True)


if os.getenv("CLEANUP_RUNS_ON_START", "1") == "1":
    _cleanup_previous_runs()
if os.getenv("CLEANUP_RUNS_ON_EXIT", "1") == "1":
    atexit.register(_cleanup_session_runs)


def _run_scan(
    analysis_mode: str,
    catalogue_key: str | None,
    input_mode: str,
    genetic_code_table: int,
    fasta_file: str | None,
    ncbi_url: str | None,
    epitope_file: str | None,
    epitope_column: str | None,
    epitope_separator: str,
    progress: gr.Progress = gr.Progress(),  # noqa: B008 - Gradio dependency injection
) -> Generator[tuple, None, None]:
    try:
        mode = parse_analysis_mode(analysis_mode)
    except ValueError as exc:
        raise gr.Error(str(exc)) from exc
    output_policy = output_policy_for(mode)

    epitope_df: pd.DataFrame | None = None
    catalogue = None
    if output_policy is OutputPolicy.FULL:
        if not epitope_file:
            raise gr.Error("Please upload an epitope file (CSV/TSV/XLSX).")
        if not epitope_column:
            raise gr.Error("Please select an epitope column from the dropdown.")
        epitope_path = Path(epitope_file)
    else:
        try:
            catalogue = get_configured_catalogue(catalogue_key)
        except ProtectedCatalogueError as exc:
            raise gr.Error(str(exc)) from exc

    if not fasta_file and not (ncbi_url or "").strip():
        sequence_type = "nucleotide" if input_mode == NUCLEOTIDE_MODE else "protein"
        raise gr.Error(
            f"Please upload a {sequence_type} FASTA file or provide an NCBI URL."
        )

    genetic_code_name: str | None = None
    if input_mode == NUCLEOTIDE_MODE:
        try:
            genetic_code_table = int(genetic_code_table)
            genetic_code_name = genetic_code_table_name(genetic_code_table)
        except (TypeError, ValueError) as exc:
            raise gr.Error(str(exc)) from exc

    def report_progress(stage: str, completed: int | None, total: int | None) -> None:
        if completed is not None and total:
            progress((completed, total), desc=stage)
        elif completed is not None:
            progress(0, desc=f"{stage} ({format_bytes(completed)})")
        else:
            progress(0, desc=stage)

    try:
        prepared_fasta = prepare_fasta(
            uploaded_path=fasta_file,
            ncbi_url=ncbi_url,
            progress=report_progress,
        )
    except InputPreparationError as exc:
        raise gr.Error(str(exc)) from exc

    nucleotide_path = prepared_fasta.path if input_mode == NUCLEOTIDE_MODE else None
    protein_path = prepared_fasta.path if input_mode == PROTEIN_MODE else None

    try:
        run_id = datetime.now(tz=UTC).strftime("%Y%m%d_%H%M%S_%f")
        output_dir = RUNS_DIR / f"run_{run_id}"
        _SESSION_RUN_DIRS.add(output_dir)

        if output_policy is OutputPolicy.FULL:
            epitope_df = read_epitope_table(epitope_path, epitope_separator)
            if epitope_column not in epitope_df.columns:
                raise gr.Error(
                    f"Selected column '{epitope_column}' is not in the epitope file."
                )
            scan_item_count = (
                epitope_df[epitope_column]
                .dropna()
                .astype(str)
                .str.strip()
                .ne("")
                .sum()
            )
        else:
            assert catalogue is not None
            scan_item_count = len(catalogue.entries)

        empty_df = pd.DataFrame()
        translated_output: Path | None = None
        source_note = " (decompressed from gzip)" if prepared_fasta.was_gzip else ""
        safe_source_name = prepared_fasta.source_name.replace("`", "'")
        genetic_code_summary = ""
        if input_mode == NUCLEOTIDE_MODE:
            assert genetic_code_name is not None
            genetic_code_summary = (
                f"- Genetic code: `{genetic_code_table} — {genetic_code_name}`\n"
            )
            seq_count = prepared_fasta.record_count
            yield (
                (
                    "Computing 6-frame translation, please be patient.\n\n"
                    "This can take up to 5 minutes for large datasets.\n\n"
                    f"- FASTA source: `{safe_source_name}`{source_note}\n"
                    f"- Input nucleotide sequences: `{seq_count}`\n"
                    f"{genetic_code_summary}"
                    + (
                        f"- Epitopes to scan: `{int(scan_item_count)}`"
                        if output_policy is OutputPolicy.FULL
                        else (
                            "- Catalogue antibodies to scan: "
                            f"`{scan_item_count}`"
                        )
                    )
                ),
                empty_df,
                empty_df,
                None,
                None,
                None,
            )
            output_dir.mkdir(parents=True, exist_ok=True)
            translated_output = output_dir / "output6frame.fasta"
            assert nucleotide_path is not None
            for translated_count, total_count in write_six_frame_fasta_with_progress(
                nucleotide_path,
                translated_output,
                genetic_code_table,
            ):
                progress(
                    (translated_count, total_count),
                    desc="Computing six-frame translation",
                )
        else:
            seq_count = prepared_fasta.record_count

        output_dir.mkdir(parents=True, exist_ok=True)
        protein_to_scan = translated_output or protein_path
        assert protein_to_scan is not None
        scan_record_count = seq_count * 6 if translated_output else seq_count

        unique_epitopes: list[str] = []
        if output_policy is OutputPolicy.FULL:
            assert epitope_df is not None
            assert epitope_column is not None
            epitope_df = epitope_df.copy()
            epitope_df["epitope_query"] = epitope_df[epitope_column].apply(
                normalize_epitope
            )
            epitope_df = epitope_df.dropna(subset=["epitope_query"])
            unique_epitopes = sorted(set(epitope_df["epitope_query"].tolist()))

        scan_intro = (
            "Scanning translated/protein sequences for epitope matches.\n\n"
            if output_policy is OutputPolicy.FULL
            else "Scanning the protected antibody catalogue.\n\n"
        )
        yield (
            (
                scan_intro
                + f"- FASTA source: `{safe_source_name}`{source_note}\n"
                f"- Protein sequences to scan: `{scan_record_count}`\n"
                f"{genetic_code_summary}"
                + (
                    f"- Epitopes to scan: `{len(unique_epitopes)}`"
                    if output_policy is OutputPolicy.FULL
                    else f"- Catalogue antibodies to scan: `{scan_item_count}`"
                )
            ),
            empty_df,
            empty_df,
            None,
            None,
            None,
        )

        if output_policy is OutputPolicy.FULL:
            hits_df = pd.DataFrame(
                columns=["epitope_query", "target_id", "target_description"]
            )
            scan_gen = scan_epitopes_with_progress(
                unique_epitopes,
                protein_to_scan,
                total_records=scan_record_count,
            )
            while True:
                try:
                    scanned_records, total_records, _hits_so_far = next(scan_gen)
                    progress(
                        (scanned_records, total_records),
                        desc="Scanning protein sequences",
                    )
                except StopIteration as stop:
                    hits_df = stop.value
                    break

            hits_path = output_dir / "epitope_hits.csv"
            if hits_df.empty:
                hits_df = pd.DataFrame(
                    columns=["epitope_query", "target_id", "target_description"]
                )
            hits_df.to_csv(hits_path, index=False)

            assert epitope_df is not None
            merged = epitope_df.merge(hits_df, on="epitope_query", how="inner")
            matched_path = output_dir / "matched_epitope_rows.csv"
            merged.to_csv(matched_path, index=False)

            hits_df = pd.read_csv(hits_path)
            matched_df = pd.read_csv(matched_path)
            summary = (
                f"Run complete.\n\n"
                f"- Output directory: `{output_dir}`\n"
                f"- Unique epitopes scanned: `{len(unique_epitopes)}`\n"
                f"{genetic_code_summary}"
                f"- Total hits: `{len(hits_df)}`\n"
                f"- Matched metadata rows: `{len(matched_df)}`"
            )
            progress(1, desc="Run complete")
            yield (
                summary,
                hits_df,
                matched_df,
                str(hits_path),
                str(matched_path),
                str(translated_output) if translated_output else None,
            )
        else:
            assert catalogue is not None
            restricted_df = pd.DataFrame()
            target_match_count = 0
            scan_gen = scan_protected_catalogue_with_progress(
                catalogue,
                protein_to_scan,
                total_records=scan_record_count,
            )
            while True:
                try:
                    scanned_records, total_records, target_match_count = next(scan_gen)
                    progress(
                        (scanned_records, total_records),
                        desc="Scanning protein sequences",
                    )
                except StopIteration as stop:
                    restricted_df = stop.value
                    break

            restricted_path = output_dir / "catalogue_matches.csv"
            restricted_df.to_csv(restricted_path, index=False)
            summary = (
                "Protected catalogue scan complete.\n\n"
                f"{genetic_code_summary}"
                f"- Positive antibody–target matches: `{target_match_count}`"
            )
            progress(1, desc="Run complete")
            yield (
                summary,
                restricted_df,
                empty_df,
                str(restricted_path),
                None,
                str(translated_output) if translated_output else None,
            )
    finally:
        prepared_fasta.cleanup()


def _upload_received(uploaded_path: str | None) -> str:
    """Report when Gradio's native browser-to-server upload has completed."""

    if not uploaded_path:
        return ""
    path = Path(uploaded_path)
    try:
        size = path.stat().st_size
    except OSError:
        return "Upload could not be read by the server."
    return (
        f"Upload complete: `{path.name}` ({format_bytes(size)}). "
        "It will be decompressed and validated when the scan starts."
    )


def _load_epitope_columns(epitope_file: str | None, epitope_separator: str):
    if not epitope_file:
        return gr.update(choices=[], value=None, interactive=False)

    epitope_path = Path(epitope_file)
    try:
        try:
            headers = list(pd.read_excel(epitope_path, nrows=0).columns)
        except Exception:
            headers = list(
                pd.read_csv(epitope_path, sep=epitope_separator, nrows=0).columns
            )
    except Exception:
        return gr.update(choices=[], value=None, interactive=False)

    if not headers:
        return gr.update(choices=[], value=None, interactive=False)

    default_column = (
        "epitope_specificity" if "epitope_specificity" in headers else headers[0]
    )
    return gr.update(choices=headers, value=default_column, interactive=True)


def _set_analysis_mode(analysis_mode: str):
    """Show only the inputs and result surfaces permitted for the mode."""

    try:
        mode = parse_analysis_mode(analysis_mode)
    except ValueError:
        mode = AnalysisMode.PROTECTED_CATALOGUE
    user_mode = mode is AnalysisMode.USER_SUPPLIED
    return (
        gr.update(visible=user_mode),
        gr.update(visible=not user_mode),
        gr.update(
            label="Epitope hits" if user_mode else "Protected catalogue matches",
            value=None,
        ),
        gr.update(visible=user_mode, value=None),
        gr.update(
            label=(
                "Download: epitope_hits.csv"
                if user_mode
                else "Download: catalogue_matches.csv"
            ),
            value=None,
        ),
        gr.update(visible=user_mode, value=None),
        gr.update(visible=True, value=None),
    )


def _update_genetic_code_visibility(input_mode: str):
    return gr.update(visible=input_mode == NUCLEOTIDE_MODE)


def build_app() -> gr.Blocks:
    try:
        catalogue_choices = get_configured_catalogue_choices()
    except ProtectedCatalogueError:
        catalogue_choices = []

    with gr.Blocks(title="SixPack-AbScan", css=APP_CSS, head=APP_HEAD) as app:
        gr.Markdown(
            """
# SixPack-AbScan

Interactive epitope matching for antibody cross-reactivity prediction.

<span style="font-size: 0.9em; color: gray;">
For the screening of large target sequence files (>2-3 Gb) we recommend the use of the command line version of SixPack-AbScan, or a local deployment of the app.
</span>
"""
        )

        gr.Markdown("### Analysis mode")

        analysis_mode = gr.Radio(
            choices=ANALYSIS_MODE_CHOICES,
            value=AnalysisMode.USER_SUPPLIED.value,
            label="Choose how antibody epitopes are supplied",
        )

        with gr.Group(visible=True) as user_epitope_inputs:
            gr.Markdown("### Antibody information")

            with gr.Row():
                epitope_file = gr.File(
                    label="Upload here your file with the list of epitopes to search",
                    file_count="single",
                    type="filepath",
                )

            with gr.Row():
                epitope_column = gr.Dropdown(
                    label="Epitope column",
                    choices=[],
                    value=None,
                    interactive=False,
                )
                epitope_separator = gr.Textbox(
                    label="CSV/TSV separator",
                    value=";",
                )

        with gr.Group(visible=False) as protected_catalogue_information:
            gr.Markdown(
                "### Protected antibody catalogue\n"
                "Choose the pre-loaded commercial antibody catalogue. "
                "The exact epitope mapping information is confidential, hence not "
                "available for download, but is available to the app and used in the "
                "back-end for cross-reactivity prediction."
            )
            catalogue_key = gr.Dropdown(
                choices=catalogue_choices,
                value=catalogue_choices[0][1] if catalogue_choices else None,
                label="Choose an antibody catalogue from the drop-down menu",
                interactive=bool(catalogue_choices),
            )

        gr.Markdown("### Crossreactivity prediction")

        with gr.Row():
            input_mode = gr.Radio(
                choices=[
                    NUCLEOTIDE_MODE,
                    PROTEIN_MODE,
                ],
                value=NUCLEOTIDE_MODE,
                label="On which file type you want to perform the search?",
            )

            genetic_code_table = gr.Dropdown(
                choices=GENETIC_CODE_CHOICES,
                value=DEFAULT_GENETIC_CODE_TABLE,
                label="Genetic code (NCBI translation table)",
                info="Used only when translating nucleotide FASTA input.",
            )

        gr.Markdown(
            "Provide **one** sequence source for the selected mode: upload a FASTA/"
            "FASTA.GZ file, or paste a direct HTTPS file URL on an NCBI host. "
            "A progress bar will show while an upload is in progress."
        )

        fasta_file = gr.File(
            label="Upload sequence FASTA or FASTA.GZ",
            file_count="single",
            file_types=FASTA_FILE_TYPES,
            type="filepath",
        )
        fasta_upload_status = gr.Markdown()

        ncbi_url = gr.Textbox(
            label="Or use a direct NCBI FASTA file URL",
            placeholder="https://ftp.ncbi.nlm.nih.gov/.../genomic.fna.gz",
            info="HTTPS only; the URL and every redirect must remain on ncbi.nlm.nih.gov.",
        )

        run_button = gr.Button("Run Scan", variant="primary")

        summary = gr.Markdown()
        hits_table = gr.Dataframe(label="Epitope hits", interactive=False)
        matched_table = gr.Dataframe(
            label="Matched epitope metadata rows", interactive=False
        )

        with gr.Row():
            hits_download = gr.File(label="Download: epitope_hits.csv")
            matched_download = gr.File(label="Download: matched_epitope_rows.csv")
            translated_download = gr.File(
                label="Download: output6frame.fasta (if generated)"
            )

        run_button.click(
            fn=_run_scan,
            inputs=[
                analysis_mode,
                catalogue_key,
                input_mode,
                genetic_code_table,
                fasta_file,
                ncbi_url,
                epitope_file,
                epitope_column,
                epitope_separator,
            ],
            outputs=[
                summary,
                hits_table,
                matched_table,
                hits_download,
                matched_download,
                translated_download,
            ],
            show_progress="minimal",
            show_progress_on=summary,
        )

        analysis_mode.change(
            fn=_set_analysis_mode,
            inputs=[analysis_mode],
            outputs=[
                user_epitope_inputs,
                protected_catalogue_information,
                hits_table,
                matched_table,
                hits_download,
                matched_download,
                translated_download,
            ],
            show_progress="hidden",
        )

        fasta_file.upload(
            fn=_upload_received,
            inputs=[fasta_file],
            outputs=[fasta_upload_status],
            show_progress="hidden",
        )

        epitope_file.change(
            fn=_load_epitope_columns,
            inputs=[epitope_file, epitope_separator],
            outputs=[epitope_column],
            show_progress="hidden",
        )
        epitope_separator.change(
            fn=_load_epitope_columns,
            inputs=[epitope_file, epitope_separator],
            outputs=[epitope_column],
            show_progress="hidden",
        )
        input_mode.change(
            fn=_update_genetic_code_visibility,
            inputs=[input_mode],
            outputs=[genetic_code_table],
            show_progress="hidden",
        )

        gr.Markdown(
            "## How to cite this app.\n"
            "Please cite this web as: Grillo, 2026. SixPack-AbScan: "
            "predict crossreactivity of monoclonal antibodies on non-target species. "
            "DOI: https://doi.org/10.82595/scilifelab.4c3a-tp57"
        )

    app.queue()
    return app


def main() -> None:
    env_port = os.getenv("PORT") or os.getenv("GRADIO_SERVER_PORT")
    server_port = int(env_port) if env_port else _find_free_port()
    server_name = os.getenv("HOST", "0.0.0.0")
    root_path = os.getenv("ROOT_PATH")
    print(
        "Starting SixPack-AbScan",
        flush=True,
    )
    print(
        f"Runtime config: HOST={server_name} PORT={server_port} ROOT_PATH={root_path!r}",
        flush=True,
    )
    print(
        f"Working directory: {Path.cwd()}",
        flush=True,
    )
    print(
        f"Runs directory: {RUNS_DIR.resolve()}",
        flush=True,
    )
    try:
        build_app().launch(
            server_name=server_name,
            server_port=server_port,
            root_path=root_path,
            max_file_size=MAX_FASTA_SOURCE_BYTES,
            blocked_paths=_blocked_catalogue_paths(),
        )
    except Exception:
        print("Application startup failed:", flush=True)
        traceback.print_exc()
        raise


if __name__ == "__main__":
    main()
