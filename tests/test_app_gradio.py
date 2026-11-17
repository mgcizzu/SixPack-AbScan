from __future__ import annotations

import csv
from collections import Counter
import os
import tempfile
import unittest
from pathlib import Path
from unittest.mock import patch

os.environ["CLEANUP_RUNS_ON_START"] = "0"
os.environ["CLEANUP_RUNS_ON_EXIT"] = "0"

import app_gradio  # noqa: E402
from analysis_modes import AnalysisMode  # noqa: E402
from app_gradio import (  # noqa: E402
    DEMO_EPITOPE_FILE,
    DEMO_EXPECTED_HITS,
    DMEL_PROTEOME_URL,
    NUCLEOTIDE_MODE,
    PROTEIN_MODE,
    PUBLIC_DSHB_CATALOGUE,
    _load_dmel_demo,
    _load_epitope_columns,
    _run_scan,
    build_app,
)


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
        self.assertEqual(
            analysis_modes[0]["props"]["choices"],
            [
                ("Upload my epitope list", AnalysisMode.USER_SUPPLIED.value),
                (
                    "Scan a pre-loaded commercial antibody catalogue",
                    AnalysisMode.PROTECTED_CATALOGUE.value,
                ),
            ],
        )

        visibility_dependency = next(
            dependency
            for dependency in self.config["dependencies"]
            if dependency.get("api_name") == "_set_analysis_mode"
        )
        self.assertEqual(visibility_dependency["show_progress"], "hidden")

    def test_uses_requested_catalogue_copy(self) -> None:
        catalogue_dropdown = next(
            component
            for component in self.config["components"]
            if component.get("type") == "dropdown"
            and component.get("props", {}).get("label")
            == "Choose an antibody catalogue from the drop-down menu"
        )
        self.assertNotIn("info", catalogue_dropdown["props"])

        markdown_values = [
            component.get("props", {}).get("value", "")
            for component in self.config["components"]
            if component.get("type") == "markdown"
        ]
        self.assertTrue(
            any(
                "The exact epitope mapping information is confidential, hence not "
                "available for download, but is available to the app and used in the "
                "back-end for cross-reactivity prediction."
                in value
                for value in markdown_values
            )
        )

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

    def test_includes_help_tutorial_examples_and_license(self) -> None:
        markdown_values = [
            component.get("props", {}).get("value", "")
            for component in self.config["components"]
            if component.get("type") == "markdown"
        ]
        combined_markdown = "\n".join(markdown_values)

        self.assertIn("What the application does", combined_markdown)
        self.assertIn(
            "Worked example: DSHB antibodies against Drosophila proteins",
            combined_markdown,
        )
        self.assertIn("MIT License", combined_markdown)
        self.assertIn("free of charge for non-commercial use", combined_markdown)
        self.assertIn("Load demo data", combined_markdown)
        self.assertIn("13 isoform-level rows", combined_markdown)
        self.assertIn("Optionally, download the test files", combined_markdown)
        self.assertNotIn(
            "not experimentally validated cross-reactivity", combined_markdown
        )
        self.assertNotIn("Load quick demo", combined_markdown)
        self.assertNotIn("Try the Drosophila demo", combined_markdown)

        demo_buttons = [
            component
            for component in self.config["components"]
            if component.get("type") == "button"
            and component.get("props", {}).get("value") == "Load demo data"
        ]
        self.assertEqual(len(demo_buttons), 1)

        download_labels = {
            component.get("props", {}).get("label")
            for component in self.config["components"]
            if component.get("type") == "downloadbutton"
        }
        self.assertIn("Download demo DSHB subset", download_labels)
        self.assertIn("Download expected Dmel hits", download_labels)
        self.assertNotIn("Download complete DSHB catalogue", download_labels)

        dependency_names = {
            dependency.get("api_name")
            for dependency in self.config["dependencies"]
        }
        self.assertIn("_load_dmel_demo", dependency_names)


class BundledExampleTests(unittest.TestCase):
    def test_public_documentation_uses_current_demo_controls(self) -> None:
        documentation = "\n".join(
            path.read_text(encoding="utf-8")
            for path in (
                Path("README.md"),
                Path("docs/user-guide.md"),
                Path("docs/tutorial.md"),
            )
        )

        self.assertIn("Load demo data", documentation)
        self.assertNotIn("Load Dmel demo", documentation)
        self.assertNotIn("Try the Drosophila demo", documentation)
        self.assertNotIn(
            "provides the complete public 49-record DSHB catalogue",
            documentation,
        )

    def test_example_files_are_present(self) -> None:
        self.assertTrue(PUBLIC_DSHB_CATALOGUE.is_file())
        self.assertTrue(DEMO_EPITOPE_FILE.is_file())
        self.assertTrue(DEMO_EXPECTED_HITS.is_file())

    def test_public_catalogue_has_the_documented_records(self) -> None:
        with PUBLIC_DSHB_CATALOGUE.open(newline="", encoding="utf-8") as handle:
            rows = list(csv.DictReader(handle))

        self.assertEqual(len(rows), 49)
        self.assertEqual(
            set(rows[0]),
            {
                "antibody_id",
                "manufacturer",
                "catalog_number",
                "antibody_name",
                "epitope_sequence",
            },
        )

        column_update = _load_epitope_columns(
            str(PUBLIC_DSHB_CATALOGUE),
            ",",
        )
        self.assertEqual(column_update["value"], "epitope_sequence")

    def test_demo_loader_selects_the_full_proteome(self) -> None:
        update = _load_dmel_demo()

        self.assertEqual(update[0], AnalysisMode.USER_SUPPLIED.value)
        self.assertEqual(update[8], PROTEIN_MODE)
        self.assertIsNone(update[10])
        self.assertEqual(update[12], DMEL_PROTEOME_URL)
        self.assertEqual(update[13], str(DEMO_EPITOPE_FILE))
        self.assertEqual(update[15], ",")

    def test_expected_proteome_result_has_13_documented_hits(self) -> None:
        with DEMO_EXPECTED_HITS.open(newline="", encoding="utf-8") as handle:
            rows = list(csv.DictReader(handle))

        self.assertEqual(len(rows), 13)
        self.assertEqual(
            Counter(row["epitope_query"] for row in rows),
            {
                "ELGLNEAQIKI": 5,
                "DPKELRRWLREMENKLESSPTITELTLFGNAELQRHLAVHS": 4,
                "FDIGLSRSMGLGVPPH": 4,
            },
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
