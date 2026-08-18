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

from fasta_input import (
    MAX_FASTA_SOURCE_BYTES,
    InputPreparationError,
    format_bytes,
    prepare_fasta,
)
from sixpack_abscan import (
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
    input_mode: str,
    fasta_file: str | None,
    ncbi_url: str | None,
    epitope_file: str | None,
    epitope_column: str | None,
    epitope_separator: str,
    progress: gr.Progress = gr.Progress(),  # noqa: B008 - Gradio dependency injection
) -> Generator[tuple, None, None]:
    if not epitope_file:
        raise gr.Error("Please upload an epitope file (CSV/TSV/XLSX).")
    if not epitope_column:
        raise gr.Error("Please select an epitope column from the dropdown.")

    epitope_path = Path(epitope_file)

    if not fasta_file and not (ncbi_url or "").strip():
        sequence_type = "nucleotide" if input_mode == NUCLEOTIDE_MODE else "protein"
        raise gr.Error(
            f"Please upload a {sequence_type} FASTA file or provide an NCBI URL."
        )

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

        epitope_df = read_epitope_table(epitope_path, epitope_separator)
        if epitope_column not in epitope_df.columns:
            raise gr.Error(
                f"Selected column '{epitope_column}' is not in the epitope file."
            )
        epitope_count = (
            epitope_df[epitope_column].dropna().astype(str).str.strip().ne("").sum()
            if epitope_column in epitope_df.columns
            else 0
        )

        empty_df = pd.DataFrame()
        translated_output: Path | None = None
        source_note = " (decompressed from gzip)" if prepared_fasta.was_gzip else ""
        safe_source_name = prepared_fasta.source_name.replace("`", "'")
        if input_mode == NUCLEOTIDE_MODE:
            seq_count = prepared_fasta.record_count
            yield (
                (
                    "Computing 6-frame translation, please be patient.\n\n"
                    "This can take up to 5 minutes for large datasets.\n\n"
                    f"- FASTA source: `{safe_source_name}`{source_note}\n"
                    f"- Input nucleotide sequences: `{seq_count}`\n"
                    f"- Epitopes to scan: `{int(epitope_count)}`"
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
                nucleotide_path, translated_output
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

        epitope_df = epitope_df.copy()
        epitope_df["epitope_query"] = epitope_df[epitope_column].apply(
            normalize_epitope
        )
        epitope_df = epitope_df.dropna(subset=["epitope_query"])
        unique_epitopes = sorted(set(epitope_df["epitope_query"].tolist()))

        yield (
            (
                "Scanning translated/protein sequences for epitope matches.\n\n"
                f"- FASTA source: `{safe_source_name}`{source_note}\n"
                f"- Protein sequences to scan: `{scan_record_count}`\n"
                f"- Epitopes to scan: `{len(unique_epitopes)}`"
            ),
            empty_df,
            empty_df,
            None,
            None,
            None,
        )

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

        merged = epitope_df.merge(hits_df, on="epitope_query", how="inner")
        matched_path = output_dir / "matched_epitope_rows.csv"
        merged.to_csv(matched_path, index=False)
        translated_path = translated_output

        hits_df = pd.read_csv(hits_path)
        matched_df = pd.read_csv(matched_path)

        summary = (
            f"Run complete.\n\n"
            f"- Output directory: `{output_dir}`\n"
            f"- Unique epitopes scanned: `{len(unique_epitopes)}`\n"
            f"- Total hits: `{len(hits_df)}`\n"
            f"- Matched metadata rows: `{len(matched_df)}`"
        )

        translated_download = str(translated_path) if translated_path else None
        progress(1, desc="Run complete")
        yield (
            summary,
            hits_df,
            matched_df,
            str(hits_path),
            str(matched_path),
            translated_download,
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


def build_app() -> gr.Blocks:
    with gr.Blocks(title="SixPack-AbScan", css=APP_CSS, head=APP_HEAD) as app:
        gr.Markdown(
            "# SixPack-AbScan\n"
            "Interactive epitope matching for antibody cross-reactivity prediction."
        )

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

        gr.Markdown(
            "Provide **one** sequence source for the selected mode: upload a FASTA/"
            "FASTA.GZ file, or paste a direct HTTPS file URL on an NCBI host. "
            "Gradio shows transfer progress while an upload is in progress."
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
                input_mode,
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
        )
    except Exception:
        print("Application startup failed:", flush=True)
        traceback.print_exc()
        raise


if __name__ == "__main__":
    main()
