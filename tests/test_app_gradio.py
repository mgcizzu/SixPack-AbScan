from __future__ import annotations

import os
import tempfile
import unittest
from pathlib import Path
from unittest.mock import patch

os.environ["CLEANUP_RUNS_ON_START"] = "0"
os.environ["CLEANUP_RUNS_ON_EXIT"] = "0"

import app_gradio  # noqa: E402
from analysis_modes import AnalysisMode  # noqa: E402
from app_gradio import NUCLEOTIDE_MODE, _run_scan, build_app  # noqa: E402


class GradioAppConfigurationTests(unittest.TestCase):
    @classmethod
    def setUpClass(cls) -> None:
        cls.config = build_app().get_config_file()

    def test_has_one_shared_sequence_upload(self) -> None:
        sequence_uploads = [
            component
            for component in self.config["components"]
            if component.get("type") == "file"
            and component.get("props", {}).get("label")
            == "Upload sequence FASTA or FASTA.GZ"
        ]

        self.assertEqual(len(sequence_uploads), 1)

        run_dependency = next(
            dependency
            for dependency in self.config["dependencies"]
            if dependency.get("api_name") == "_run_scan"
        )
        self.assertEqual(len(run_dependency["inputs"]), 9)

    def test_has_two_analysis_modes(self) -> None:
        analysis_modes = [
            component
            for component in self.config["components"]
            if component.get("type") == "radio"
            and component.get("props", {}).get("label")
            == "Choose how antibody epitopes are supplied"
        ]

        self.assertEqual(len(analysis_modes), 1)
        self.assertEqual(
            analysis_modes[0]["props"]["value"], AnalysisMode.USER_SUPPLIED.value
        )

        visibility_dependency = next(
            dependency
            for dependency in self.config["dependencies"]
            if dependency.get("api_name") == "_set_analysis_mode"
        )
        self.assertEqual(visibility_dependency["show_progress"], "hidden")

    def test_has_nucleotide_genetic_code_dropdown(self) -> None:
        genetic_code_dropdowns = [
            component
            for component in self.config["components"]
            if component.get("type") == "dropdown"
            and component.get("props", {}).get("label")
            == "Genetic code (NCBI translation table)"
        ]

        self.assertEqual(len(genetic_code_dropdowns), 1)
        self.assertEqual(genetic_code_dropdowns[0]["props"]["value"], 1)

        visibility_dependency = next(
            dependency
            for dependency in self.config["dependencies"]
            if dependency.get("api_name") == "_update_genetic_code_visibility"
        )
        self.assertEqual(visibility_dependency["show_progress"], "hidden")

    def test_uses_only_explicit_processing_progress(self) -> None:
        dependencies = self.config["dependencies"]
        run_dependency = next(
            dependency
            for dependency in dependencies
            if dependency.get("api_name") == "_run_scan"
        )
        upload_dependency = next(
            dependency
            for dependency in dependencies
            if dependency.get("api_name") == "_upload_received"
        )
        column_dependencies = [
            dependency
            for dependency in dependencies
            if str(dependency.get("api_name", "")).startswith(
                "_load_epitope_columns"
            )
        ]

        self.assertEqual(run_dependency["show_progress"], "minimal")
        summary_id = run_dependency["outputs"][0]
        self.assertEqual(run_dependency["show_progress_on"], [summary_id])
        self.assertEqual(upload_dependency["show_progress"], "hidden")
        self.assertEqual(len(column_dependencies), 2)
        self.assertTrue(
            all(
                dependency["show_progress"] == "hidden"
                for dependency in column_dependencies
            )
        )


class ScanProgressTests(unittest.TestCase):
    def test_translation_updates_progress_without_refreshing_summary(self) -> None:
        with tempfile.TemporaryDirectory() as directory:
            directory_path = Path(directory)
            fasta_path = directory_path / "input.fna"
            fasta_path.write_text(">sequence\nATGGCCATG\n", encoding="utf-8")
            epitope_path = directory_path / "epitopes.csv"
            epitope_path.write_text(
                "epitope_specificity\nMA\n", encoding="utf-8"
            )
            progress_descriptions: list[str | None] = []

            def capture_progress(_value: object, *, desc: str | None = None) -> None:
                progress_descriptions.append(desc)

            with (
                patch.object(app_gradio, "RUNS_DIR", directory_path / "runs"),
                patch.object(app_gradio, "_SESSION_RUN_DIRS", set()),
            ):
                updates = list(
                    _run_scan(
                        AnalysisMode.USER_SUPPLIED.value,
                        None,
                        NUCLEOTIDE_MODE,
                        1,
                        str(fasta_path),
                        None,
                        str(epitope_path),
                        "epitope_specificity",
                        ";",
                        progress=capture_progress,
                    )
                )

            summaries = [update[0] for update in updates]
            self.assertTrue(
                any(
                    description == "Computing six-frame translation"
                    for description in progress_descriptions
                )
            )
            self.assertFalse(
                any("Translated sequences" in summary for summary in summaries)
            )

    def test_selected_genetic_code_is_used_for_translation(self) -> None:
        with tempfile.TemporaryDirectory() as directory:
            directory_path = Path(directory)
            fasta_path = directory_path / "input.fna"
            fasta_path.write_text(">sequence\nATATGAAGA\n", encoding="utf-8")
            epitope_path = directory_path / "epitopes.csv"
            epitope_path.write_text(
                "epitope_specificity\nMW\n", encoding="utf-8"
            )

            with (
                patch.object(app_gradio, "RUNS_DIR", directory_path / "runs"),
                patch.object(app_gradio, "_SESSION_RUN_DIRS", set()),
            ):
                updates = list(
                    _run_scan(
                        AnalysisMode.USER_SUPPLIED.value,
                        None,
                        NUCLEOTIDE_MODE,
                        2,
                        str(fasta_path),
                        None,
                        str(epitope_path),
                        "epitope_specificity",
                        ";",
                        progress=lambda *_args, **_kwargs: None,
                    )
                )

        final_summary, final_hits = updates[-1][0], updates[-1][1]
        self.assertIn("2 — Vertebrate Mitochondrial", final_summary)
        self.assertEqual(final_hits["epitope_query"].tolist(), ["MW"])


if __name__ == "__main__":
    unittest.main()
