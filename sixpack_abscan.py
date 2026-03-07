#!/usr/bin/env python3
"""SixPack-AbScan CLI.

Predict antibody cross-reactivity by scanning epitope peptides against:
1) a six-frame translation of a nucleotide FASTA, or
2) a precomputed protein FASTA.
"""

from __future__ import annotations

import argparse
from pathlib import Path
from typing import Generator, Iterable

import pandas as pd
from Bio import SeqIO


GENETIC_CODE = {
    "ATA": "I",
    "ATC": "I",
    "ATT": "I",
    "ATG": "M",
    "ACA": "T",
    "ACC": "T",
    "ACG": "T",
    "ACT": "T",
    "AAC": "N",
    "AAT": "N",
    "AAA": "K",
    "AAG": "K",
    "AGC": "S",
    "AGT": "S",
    "AGA": "R",
    "AGG": "R",
    "CTA": "L",
    "CTC": "L",
    "CTG": "L",
    "CTT": "L",
    "CCA": "P",
    "CCC": "P",
    "CCG": "P",
    "CCT": "P",
    "CAC": "H",
    "CAT": "H",
    "CAA": "Q",
    "CAG": "Q",
    "CGA": "R",
    "CGC": "R",
    "CGG": "R",
    "CGT": "R",
    "GTA": "V",
    "GTC": "V",
    "GTG": "V",
    "GTT": "V",
    "GCA": "A",
    "GCC": "A",
    "GCG": "A",
    "GCT": "A",
    "GAC": "D",
    "GAT": "D",
    "GAA": "E",
    "GAG": "E",
    "GGA": "G",
    "GGC": "G",
    "GGG": "G",
    "GGT": "G",
    "TCA": "S",
    "TCC": "S",
    "TCG": "S",
    "TCT": "S",
    "TTC": "F",
    "TTT": "F",
    "TTA": "L",
    "TTG": "L",
    "TAC": "Y",
    "TAT": "Y",
    "TAA": "_",
    "TAG": "_",
    "TGC": "C",
    "TGT": "C",
    "TGA": "_",
    "TGG": "W",
}

BASE_COMPLEMENT = str.maketrans({"A": "T", "C": "G", "G": "C", "T": "A", "N": "N"})


def reverse_complement(sequence: str) -> str:
    return sequence.upper().translate(BASE_COMPLEMENT)[::-1]


def translate_frame(sequence: str, offset: int) -> str:
    translated = []
    seq = sequence[offset:]
    for i in range(0, len(seq) - 2, 3):
        codon = seq[i : i + 3]
        translated.append(GENETIC_CODE.get(codon, "X"))
    return "".join(translated)


def six_frame_translation(sequence: str) -> list[str]:
    seq = sequence.upper()
    rc = reverse_complement(seq)
    return [
        translate_frame(seq, 0),
        translate_frame(seq, 1),
        translate_frame(seq, 2),
        translate_frame(rc, 0),
        translate_frame(rc, 1),
        translate_frame(rc, 2),
    ]


def wrap_fasta(seq: str, width: int = 80) -> Iterable[str]:
    for i in range(0, len(seq), width):
        yield seq[i : i + width]


def write_six_frame_fasta(input_fasta: Path, output_fasta: Path) -> None:
    for _ in write_six_frame_fasta_with_progress(input_fasta, output_fasta):
        pass


def write_six_frame_fasta_with_progress(
    input_fasta: Path, output_fasta: Path
) -> Iterable[tuple[int, int]]:
    total_records = sum(1 for _ in SeqIO.parse(str(input_fasta), "fasta"))
    written = 0
    with output_fasta.open("w", encoding="utf-8") as out_handle:
        for record in SeqIO.parse(str(input_fasta), "fasta"):
            frames = six_frame_translation(str(record.seq))
            for idx, frame_seq in enumerate(frames, start=1):
                out_handle.write(f">{record.id}|frame{idx}\n")
                for line in wrap_fasta(frame_seq):
                    out_handle.write(f"{line}\n")
            written += 1
            yield written, total_records


def read_epitope_table(path: Path, sep: str) -> pd.DataFrame:
    suffix = path.suffix.lower()
    if suffix in {".xlsx", ".xls"}:
        return pd.read_excel(path)
    if suffix in {".csv", ".tsv", ".txt"}:
        return pd.read_csv(path, sep=sep)

    # Uploaded files (e.g. via Gradio temp storage) may lose the original extension.
    # Try Excel first, then fallback to delimited text.
    try:
        return pd.read_excel(path)
    except Exception:
        return pd.read_csv(path, sep=sep)


def normalize_epitope(value: object) -> str | None:
    if pd.isna(value):
        return None
    cleaned = str(value).strip().upper().replace(" ", "")
    return cleaned or None


def scan_epitopes(epitopes: list[str], protein_fasta: Path) -> pd.DataFrame:
    hits: list[dict[str, str]] = []
    for record in SeqIO.parse(str(protein_fasta), "fasta"):
        target_id = record.id
        target_description = record.description
        sequence = str(record.seq).upper()
        for epitope in epitopes:
            if epitope in sequence:
                hits.append(
                    {
                        "epitope_query": epitope,
                        "target_id": target_id,
                        "target_description": target_description,
                    }
                )
    return pd.DataFrame(hits)


def scan_epitopes_with_progress(
    epitopes: list[str], protein_fasta: Path
) -> Generator[tuple[int, int, int], None, pd.DataFrame]:
    """Scan epitopes and yield progress as (done, total, hits_so_far)."""
    hits: list[dict[str, str]] = []
    total = len(epitopes)
    for done, epitope in enumerate(epitopes, start=1):
        for record in SeqIO.parse(str(protein_fasta), "fasta"):
            target_id = record.id
            target_description = record.description
            sequence = str(record.seq).upper()
            if epitope in sequence:
                hits.append(
                    {
                        "epitope_query": epitope,
                        "target_id": target_id,
                        "target_description": target_description,
                    }
                )
        yield done, total, len(hits)
    return pd.DataFrame(hits)


def run_abscan(
    *,
    output_dir: Path,
    epitope_file: Path,
    epitope_column: str = "epitope_specificity",
    epitope_separator: str = ";",
    input_nucleotide_fasta: Path | None = None,
    input_protein_fasta: Path | None = None,
    translated_fasta_name: str = "output6frame.fasta",
    hits_output_name: str = "epitope_hits.csv",
    matched_metadata_output_name: str = "matched_epitope_rows.csv",
) -> dict[str, Path | int]:
    if not input_nucleotide_fasta and not input_protein_fasta:
        raise ValueError(
            "Provide either input_nucleotide_fasta or input_protein_fasta."
        )

    output_dir.mkdir(parents=True, exist_ok=True)
    translated_output = output_dir / translated_fasta_name

    protein_fasta: Path
    if input_protein_fasta:
        protein_fasta = input_protein_fasta
        print(f"Using precomputed protein FASTA: {protein_fasta}")
    else:
        assert input_nucleotide_fasta is not None
        print(f"Generating six-frame translation from: {input_nucleotide_fasta}")
        write_six_frame_fasta(input_nucleotide_fasta, translated_output)
        protein_fasta = translated_output
        print(f"Six-frame FASTA written to: {translated_output}")

    epitope_df = read_epitope_table(epitope_file, epitope_separator)
    if epitope_column not in epitope_df.columns:
        raise KeyError(f"Column '{epitope_column}' was not found in {epitope_file}.")

    epitope_df = epitope_df.copy()
    epitope_df["epitope_query"] = epitope_df[epitope_column].apply(normalize_epitope)
    epitope_df = epitope_df.dropna(subset=["epitope_query"])

    unique_epitopes = sorted(set(epitope_df["epitope_query"].tolist()))
    print(f"Loaded {len(unique_epitopes)} unique epitopes for scanning.")

    hits_df = scan_epitopes(unique_epitopes, protein_fasta)
    hits_output = output_dir / hits_output_name
    if hits_df.empty:
        hits_df = pd.DataFrame(
            columns=["epitope_query", "target_id", "target_description"]
        )
    hits_df.to_csv(hits_output, index=False)
    print(f"Total hits: {len(hits_df)}")
    print(f"Saved epitope hits to: {hits_output}")

    merged = epitope_df.merge(hits_df, on="epitope_query", how="inner")
    matched_output = output_dir / matched_metadata_output_name
    merged.to_csv(matched_output, index=False)
    print(f"Saved matched epitope metadata rows to: {matched_output}")

    return {
        "hits_output": hits_output,
        "matched_output": matched_output,
        "translated_output": translated_output if not input_protein_fasta else None,
        "num_unique_epitopes": len(unique_epitopes),
        "num_hits": len(hits_df),
        "num_matched_rows": len(merged),
    }


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(
        description="Predict antibody cross-reactivity using epitope exact matches.",
    )
    parser.add_argument(
        "--input-nucleotide-fasta",
        type=Path,
        help="Nucleotide FASTA used to generate six-frame translations.",
    )
    parser.add_argument(
        "--input-protein-fasta",
        type=Path,
        help="Precomputed protein FASTA to scan directly (skips translation).",
    )
    parser.add_argument(
        "--epitope-file",
        type=Path,
        required=True,
        help="CSV/TSV/XLSX file containing epitope sequences.",
    )
    parser.add_argument(
        "--epitope-column",
        default="epitope_specificity",
        help="Column name in epitope file containing peptide sequences.",
    )
    parser.add_argument(
        "--epitope-separator",
        default=";",
        help="Delimiter used when reading CSV/TSV epitope files.",
    )
    parser.add_argument(
        "--output-dir",
        type=Path,
        required=True,
        help="Directory for all generated outputs.",
    )
    parser.add_argument(
        "--translated-fasta-name",
        default="output6frame.fasta",
        help="Output FASTA name for six-frame translations.",
    )
    parser.add_argument(
        "--hits-output-name",
        default="epitope_hits.csv",
        help="Output CSV name containing all epitope-target hits.",
    )
    parser.add_argument(
        "--matched-metadata-output-name",
        default="matched_epitope_rows.csv",
        help="Output CSV name containing matched rows from epitope metadata.",
    )
    args = parser.parse_args()
    if not args.input_nucleotide_fasta and not args.input_protein_fasta:
        parser.error(
            "Provide either --input-nucleotide-fasta or --input-protein-fasta."
        )
    return args


def main() -> None:
    args = parse_args()
    run_abscan(
        output_dir=args.output_dir,
        epitope_file=args.epitope_file,
        epitope_column=args.epitope_column,
        epitope_separator=args.epitope_separator,
        input_nucleotide_fasta=args.input_nucleotide_fasta,
        input_protein_fasta=args.input_protein_fasta,
        translated_fasta_name=args.translated_fasta_name,
        hits_output_name=args.hits_output_name,
        matched_metadata_output_name=args.matched_metadata_output_name,
    )


if __name__ == "__main__":
    main()
