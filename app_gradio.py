#!/usr/bin/env python3
"""Gradio app for SixPack-AbScan."""

from __future__ import annotations

import atexit
from datetime import datetime
import os
from pathlib import Path
import socket
from typing import Generator
import shutil

import gradio as gr
import pandas as pd
from Bio import SeqIO

from sixpack_abscan import (
    normalize_epitope,
    read_epitope_table,
    scan_epitopes_with_progress,
    write_six_frame_fasta_with_progress,
)

RUNS_DIR = Path("runs")
_SESSION_RUN_DIRS: set[Path] = set()


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
    nucleotide_fasta: str | None,
    protein_fasta: str | None,
    epitope_file: str | None,
    epitope_column: str | None,
    epitope_separator: str,
) -> Generator[tuple, None, None]:
    if not epitope_file:
        raise gr.Error("Please upload an epitope file (CSV/TSV/XLSX).")
    if not epitope_column:
        raise gr.Error("Please select an epitope column from the dropdown.")

    nucleotide_path = Path(nucleotide_fasta) if nucleotide_fasta else None
    protein_path = Path(protein_fasta) if protein_fasta else None
    epitope_path = Path(epitope_file)

    if input_mode == "Nucleotide FASTA (will be 6-frame translated automatically)":
        if not nucleotide_path:
            raise gr.Error("Please upload a nucleotide FASTA file.")
        protein_path = None
    else:
        if not protein_path:
            raise gr.Error("Please upload a precomputed protein FASTA file.")
        nucleotide_path = None

    run_id = datetime.now().strftime("%Y%m%d_%H%M%S")
    output_dir = Path("runs") / f"run_{run_id}"
    _SESSION_RUN_DIRS.add(output_dir)

    epitope_df = read_epitope_table(epitope_path, epitope_separator)
    if epitope_column not in epitope_df.columns:
        raise gr.Error(f"Selected column '{epitope_column}' is not in the epitope file.")
    epitope_count = (
        epitope_df[epitope_column].dropna().astype(str).str.strip().ne("").sum()
        if epitope_column in epitope_df.columns
        else 0
    )

    empty_df = pd.DataFrame()
    translated_output: Path | None = None
    if input_mode == "Nucleotide FASTA (will be 6-frame translated automatically)":
        seq_count = sum(1 for _ in SeqIO.parse(str(nucleotide_path), "fasta"))
        yield (
            "Computing 6-frame translation, please be patient.\n\n"
            "This can take up to 5 minutes for large datasets.\n\n"
            f"- Input nucleotide sequences: `{seq_count}`\n"
            f"- Epitopes to scan: `{int(epitope_count)}`",
            empty_df,
            empty_df,
            None,
            None,
            None,
        )
        output_dir.mkdir(parents=True, exist_ok=True)
        translated_output = output_dir / "output6frame.fasta"
        for translated_count, total_count in write_six_frame_fasta_with_progress(
            nucleotide_path, translated_output
        ):
            if (
                translated_count == 1
                or translated_count == total_count
                or translated_count % 25 == 0
            ):
                yield (
                    "Computing 6-frame translation, please be patient.\n\n"
                    "This can take up to 5 minutes for large datasets.\n\n"
                    f"- Translated sequences: `{translated_count}/{total_count}`\n"
                    f"- Epitopes to scan: `{int(epitope_count)}`",
                    empty_df,
                    empty_df,
                    None,
                    None,
                    None,
                )
    else:
        seq_count = sum(1 for _ in SeqIO.parse(str(protein_path), "fasta"))
        yield (
            "Scanning protein FASTA for epitope matches, please wait.\n\n"
            f"- Input protein sequences: `{seq_count}`\n"
            f"- Epitopes to scan: `{int(epitope_count)}`",
            empty_df,
            empty_df,
            None,
            None,
            None,
        )

    output_dir.mkdir(parents=True, exist_ok=True)
    protein_to_scan = translated_output or protein_path
    assert protein_to_scan is not None

    epitope_df = epitope_df.copy()
    epitope_df["epitope_query"] = epitope_df[epitope_column].apply(normalize_epitope)
    epitope_df = epitope_df.dropna(subset=["epitope_query"])
    unique_epitopes = sorted(set(epitope_df["epitope_query"].tolist()))

    total_hits_so_far = 0
    hits_df = pd.DataFrame(columns=["epitope_query", "target_id", "target_description"])
    scan_gen = scan_epitopes_with_progress(unique_epitopes, protein_to_scan)
    while True:
        try:
            scanned, total, total_hits_so_far = next(scan_gen)
            if scanned == 1 or scanned == total or scanned % 10 == 0:
                yield (
                    "Scanning translated/protein sequences for epitope matches.\n\n"
                    f"- Scanned epitopes: `{scanned}/{total}`\n"
                    f"- Hits found so far: `{total_hits_so_far}`",
                    empty_df,
                    empty_df,
                    None,
                    None,
                    None,
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
    yield (
        summary,
        hits_df,
        matched_df,
        str(hits_path),
        str(matched_path),
        translated_download,
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
        "epitope_specificity"
        if "epitope_specificity" in headers
        else headers[0]
    )
    return gr.update(choices=headers, value=default_column, interactive=True)


def build_app() -> gr.Blocks:
    with gr.Blocks(title="SixPack-AbScan") as app:
        gr.Markdown(
            "# SixPack-AbScan\n"
            "Interactive epitope matching for antibody cross-reactivity prediction."
        )

        with gr.Row():
            input_mode = gr.Radio(
                choices=[
                    "Nucleotide FASTA (will be 6-frame translated automatically)",
                    "Protein FASTA (precomputed proteome)",
                ],
                value="Nucleotide FASTA (will be 6-frame translated automatically)",
                label="What do you want to predict cross reactivity on?",
            )

        with gr.Row():
            nucleotide_fasta = gr.File(
                label="Nucleotide FASTA",
                file_count="single",
                type="filepath",
            )
            protein_fasta = gr.File(
                label="Protein FASTA",
                file_count="single",
                type="filepath",
            )

        with gr.Row():
            epitope_file = gr.File(
                label="Epitope file (CSV/TSV/XLSX)",
                file_count="single",
                type="filepath",
            )
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

        run_button = gr.Button("Run Scan", variant="primary")

        summary = gr.Markdown()
        hits_table = gr.Dataframe(label="Epitope hits", interactive=False)
        matched_table = gr.Dataframe(label="Matched epitope metadata rows", interactive=False)

        with gr.Row():
            hits_download = gr.File(label="Download: epitope_hits.csv")
            matched_download = gr.File(label="Download: matched_epitope_rows.csv")
            translated_download = gr.File(label="Download: output6frame.fasta (if generated)")

        run_button.click(
            fn=_run_scan,
            inputs=[
                input_mode,
                nucleotide_fasta,
                protein_fasta,
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
            show_progress="hidden",
        )

        epitope_file.change(
            fn=_load_epitope_columns,
            inputs=[epitope_file, epitope_separator],
            outputs=[epitope_column],
        )
        epitope_separator.change(
            fn=_load_epitope_columns,
            inputs=[epitope_file, epitope_separator],
            outputs=[epitope_column],
        )

    app.queue()
    return app


if __name__ == "__main__":
    env_port = os.getenv("PORT") or os.getenv("GRADIO_SERVER_PORT")
    server_port = int(env_port) if env_port else _find_free_port()
    server_name = os.getenv("HOST", "0.0.0.0")
    root_path = os.getenv("ROOT_PATH")
    build_app().launch(
        server_name=server_name,
        server_port=server_port,
        root_path=root_path,
    )
