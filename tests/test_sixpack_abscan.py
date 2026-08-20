from __future__ import annotations

import tempfile
import unittest
from pathlib import Path

import pandas as pd

from sixpack_abscan import (
    scan_epitopes,
    scan_epitopes_with_progress,
    six_frame_translation,
    translate_frame,
    write_six_frame_fasta,
)


def consume_progress_scan(
    epitopes: list[str], fasta_path: Path, *, total_records: int | None = None
) -> tuple[pd.DataFrame, list[tuple[int, int, int]]]:
    generator = scan_epitopes_with_progress(
        epitopes,
        fasta_path,
        total_records=total_records,
    )
    updates: list[tuple[int, int, int]] = []
    while True:
        try:
            updates.append(next(generator))
        except StopIteration as stopped:
            return stopped.value, updates


class EpitopeScanTests(unittest.TestCase):
    def test_record_first_progress_scan_returns_the_same_hits(self) -> None:
        fasta_text = (
            ">record-one first description\n"
            "acdacd\n"
            ">duplicate-id second description\n"
            "MNPQRST\n"
            ">duplicate-id third description\n"
            "XXACDYYMNP\n"
        )
        epitopes = ["ACD", "MNP", "NOTPRESENT"]

        with tempfile.TemporaryDirectory() as directory:
            fasta_path = Path(directory) / "proteins.faa"
            fasta_path.write_text(fasta_text, encoding="utf-8")

            expected = scan_epitopes(epitopes, fasta_path)
            actual, updates = consume_progress_scan(
                epitopes,
                fasta_path,
                total_records=3,
            )

        columns = ["epitope_query", "target_id", "target_description"]
        expected = expected.sort_values(columns).reset_index(drop=True)
        actual = actual.sort_values(columns).reset_index(drop=True)
        pd.testing.assert_frame_equal(actual, expected)
        self.assertEqual(updates[-1], (3, 3, len(actual)))

    def test_progress_scan_can_count_records_when_total_is_unknown(self) -> None:
        with tempfile.TemporaryDirectory() as directory:
            fasta_path = Path(directory) / "proteins.faa"
            fasta_path.write_text(
                ">one\nAAAA\n>two\nBBBB\n",
                encoding="utf-8",
            )

            hits, updates = consume_progress_scan(["AA"], fasta_path)

        self.assertEqual(len(hits), 1)
        self.assertEqual(updates[-1], (2, 2, 1))


class SixFrameTranslationTests(unittest.TestCase):
    def test_standard_code_remains_the_default(self) -> None:
        sequence = "ATATGAAGA"

        self.assertEqual(translate_frame(sequence, 0), "I_R")
        self.assertEqual(translate_frame(sequence, 0, 1), "I_R")

    def test_vertebrate_mitochondrial_code_changes_expected_codons(self) -> None:
        sequence = "ATATGAAGA"

        self.assertEqual(translate_frame(sequence, 0, 2), "MW_")
        self.assertEqual(six_frame_translation(sequence, 2)[0], "MW_")

    def test_alternative_start_codons_are_not_forced_to_methionine(self) -> None:
        self.assertEqual(translate_frame("GTG", 0, 11), "V")

    def test_writer_uses_the_selected_genetic_code(self) -> None:
        with tempfile.TemporaryDirectory() as directory:
            directory_path = Path(directory)
            input_path = directory_path / "input.fna"
            output_path = directory_path / "translated.faa"
            input_path.write_text(">sequence\nATATGAAGA\n", encoding="utf-8")

            write_six_frame_fasta(input_path, output_path, 2)

            output_lines = output_path.read_text(encoding="utf-8").splitlines()

        self.assertEqual(output_lines[:2], [">sequence|frame1", "MW_"])

    def test_unknown_genetic_code_is_rejected(self) -> None:
        with self.assertRaisesRegex(ValueError, "Unknown NCBI genetic code table"):
            translate_frame("ATG", 0, 999)


if __name__ == "__main__":
    unittest.main()
