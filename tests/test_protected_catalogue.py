from __future__ import annotations

import contextlib
import io
import os
import tempfile
import unittest
from pathlib import Path
from unittest.mock import patch

import pandas as pd

import app_gradio
import protected_catalogue
from analysis_modes import AnalysisMode, OutputPolicy, output_policy_for
from protected_catalogue import (
    RESTRICTED_RESULT_COLUMNS,
    ProtectedCatalogueValidationError,
    get_configured_catalogue_choices,
    load_protected_catalogue,
    scan_protected_catalogue_with_progress,
)

SECRET_MATCH = "SIINFEKL"
SECRET_MISS = "PEPTIDER"


def _write_catalogue(path: Path) -> None:
    path.write_text(
        "antibody_id,manufacturer,catalog_number,antibody_name,epitope_sequence\n"
        f"internal-a,ExampleCo,AB-100,Antibody A,{SECRET_MATCH}\n"
        f"internal-b,ExampleCo,AB-200,Antibody B,{SECRET_MISS}\n",
        encoding="utf-8",
    )


def _finish_generator(generator):
    updates = []
    while True:
        try:
            updates.append(next(generator))
        except StopIteration as stopped:
            return updates, stopped.value


class _ProgressRecorder:
    def __init__(self) -> None:
        self.calls = []

    def __call__(self, *args, **kwargs) -> None:
        self.calls.append((args, kwargs))


class AnalysisModeTests(unittest.TestCase):
    def test_each_mode_has_an_explicit_disclosure_policy(self) -> None:
        self.assertIs(output_policy_for(AnalysisMode.USER_SUPPLIED), OutputPolicy.FULL)
        self.assertIs(
            output_policy_for(AnalysisMode.PROTECTED_CATALOGUE),
            OutputPolicy.RESTRICTED,
        )

    def test_protected_ui_hides_detailed_results_but_keeps_translation(self) -> None:
        updates = app_gradio._set_analysis_mode(AnalysisMode.PROTECTED_CATALOGUE.value)

        self.assertFalse(updates[0]["visible"])
        self.assertTrue(updates[1]["visible"])
        self.assertFalse(updates[3]["visible"])
        self.assertFalse(updates[5]["visible"])
        self.assertTrue(updates[6]["visible"])

    def test_private_catalogue_is_blocked_from_gradio_file_serving(self) -> None:
        with tempfile.TemporaryDirectory() as directory:
            private_path = Path(directory) / "catalogue.csv"
            _write_catalogue(private_path)
            with patch.dict(
                os.environ,
                {"PROTECTED_CATALOGUE_PATH": str(private_path)},
                clear=False,
            ):
                self.assertEqual(
                    app_gradio._blocked_catalogue_paths(),
                    [str(private_path.resolve())],
                )

    def test_catalogue_dropdown_uses_public_labels_and_opaque_keys(self) -> None:
        with tempfile.TemporaryDirectory() as directory:
            directory_path = Path(directory)
            _write_catalogue(directory_path / "Atlas_Antibodies.csv")
            _write_catalogue(directory_path / "Example_Catalogue.csv")
            with patch.dict(
                os.environ,
                {"PROTECTED_CATALOGUE_DIR": str(directory_path)},
                clear=False,
            ):
                choices = get_configured_catalogue_choices()

        self.assertEqual(
            [label for label, _key in choices],
            ["Atlas Antibodies", "Example Catalogue"],
        )
        self.assertTrue(all("/" not in key and ".csv" not in key for _, key in choices))

    def test_browser_config_does_not_expose_private_catalogue_data(self) -> None:
        with tempfile.TemporaryDirectory() as directory:
            private_path = Path(directory) / "private-catalogue.csv"
            _write_catalogue(private_path)
            with patch.dict(
                os.environ,
                {
                    "PROTECTED_CATALOGUE_PATH": str(private_path),
                    "PROTECTED_CATALOGUE_DIR": "",
                },
                clear=False,
            ):
                browser_config = repr(app_gradio.build_app().get_config_file())

        self.assertNotIn(str(private_path), browser_config)
        self.assertNotIn(SECRET_MATCH, browser_config)
        self.assertNotIn(SECRET_MISS, browser_config)
        self.assertNotIn("internal-a", browser_config)

    def test_serve_project_volume_is_used_without_environment_configuration(
        self,
    ) -> None:
        with tempfile.TemporaryDirectory() as directory:
            catalogue_directory = Path(directory) / "catalogues"
            catalogue_directory.mkdir()
            _write_catalogue(catalogue_directory / "Atlas_Antibodies.csv")
            with (
                patch.dict(os.environ, {}, clear=False),
                patch.object(
                    protected_catalogue,
                    "DEFAULT_CATALOGUE_DIRECTORY",
                    catalogue_directory,
                ),
            ):
                os.environ.pop("PROTECTED_CATALOGUE_DIR", None)
                os.environ.pop("PROTECTED_CATALOGUE_PATH", None)
                choices = get_configured_catalogue_choices()

        self.assertEqual([label for label, _key in choices], ["Atlas Antibodies"])


class ProtectedCatalogueTests(unittest.TestCase):
    def test_single_pass_scanner_returns_only_positive_restricted_hits(self) -> None:
        with tempfile.TemporaryDirectory() as directory:
            directory_path = Path(directory)
            catalogue_path = directory_path / "catalogue.csv"
            protein_path = directory_path / "target.faa"
            catalogue_path.write_text(
                "antibody_id,manufacturer,catalog_number,antibody_name,"
                "epitope_sequence\n"
                f"internal-a,ExampleCo,AB-100,Antibody A,{SECRET_MATCH}\n"
                f"internal-b,ExampleCo,AB-200,Antibody B,{SECRET_MISS}\n"
                "internal-c,ExampleCo,AB-050,Antibody C,INFEK\n",
                encoding="utf-8",
            )
            protein_path.write_text(
                ">isoform-1 Protein isoform 1\nAASIINFEKLYY\n"
                ">isoform-2 Protein isoform 2\nQQSIINFEKLZZ\n"
                ">protein-3 Unmatched protein\nNOTAMATCH\n",
                encoding="utf-8",
            )

            catalogue = load_protected_catalogue(catalogue_path)
            updates, results = _finish_generator(
                scan_protected_catalogue_with_progress(
                    catalogue,
                    protein_path,
                    total_records=3,
                )
            )

        self.assertEqual(list(results.columns), list(RESTRICTED_RESULT_COLUMNS))
        self.assertNotIn("perfect_match", results.columns)
        self.assertEqual(
            results["catalog_number"].tolist(),
            ["AB-050", "AB-050", "AB-100", "AB-100"],
        )
        self.assertEqual(
            results["target_id"].tolist(),
            ["isoform-1", "isoform-2", "isoform-1", "isoform-2"],
        )
        self.assertEqual(
            results["target_description"].tolist(),
            [
                "isoform-1 Protein isoform 1",
                "isoform-2 Protein isoform 2",
                "isoform-1 Protein isoform 1",
                "isoform-2 Protein isoform 2",
            ],
        )
        self.assertNotIn("AB-200", results.to_csv(index=False))
        self.assertEqual(updates[-1], (3, 3, 4))
        serialized = results.to_csv(index=False)
        self.assertNotIn(SECRET_MATCH, serialized)
        self.assertNotIn(SECRET_MISS, serialized)
        self.assertNotIn("internal-a", serialized)

    def test_invalid_private_value_is_not_repeated_in_error(self) -> None:
        with tempfile.TemporaryDirectory() as directory:
            path = Path(directory) / "catalogue.csv"
            invalid_secret = "DO-NOT-LEAK"
            path.write_text(
                "antibody_id,manufacturer,catalog_number,antibody_name,"
                "epitope_sequence\n"
                f"internal-a,ExampleCo,AB-100,Antibody A,{invalid_secret}\n",
                encoding="utf-8",
            )

            with self.assertRaises(ProtectedCatalogueValidationError) as raised:
                load_protected_catalogue(path)

        self.assertNotIn(invalid_secret, str(raised.exception))

    def test_protected_app_run_does_not_publish_or_persist_secrets(self) -> None:
        with tempfile.TemporaryDirectory() as directory:
            directory_path = Path(directory)
            catalogue_path = directory_path / "private-catalogue.csv"
            protein_path = directory_path / "target.faa"
            runs_path = directory_path / "runs"
            _write_catalogue(catalogue_path)
            protein_path.write_text(">protein-1\nAASIINFEKLYY\n", encoding="utf-8")
            progress = _ProgressRecorder()
            captured_stdout = io.StringIO()
            captured_stderr = io.StringIO()

            with (
                patch.dict(
                    os.environ,
                    {"PROTECTED_CATALOGUE_PATH": str(catalogue_path)},
                    clear=False,
                ),
                patch.object(app_gradio, "RUNS_DIR", runs_path),
                patch.object(app_gradio, "_SESSION_RUN_DIRS", set()),
                contextlib.redirect_stdout(captured_stdout),
                contextlib.redirect_stderr(captured_stderr),
            ):
                events = list(
                    app_gradio._run_scan(
                        AnalysisMode.PROTECTED_CATALOGUE.value,
                        None,
                        app_gradio.PROTEIN_MODE,
                        1,
                        str(protein_path),
                        None,
                        None,
                        None,
                        ";",
                        progress=progress,
                    )
                )

            summary, results, matched, result_file, matched_file, translated = events[
                -1
            ]
            self.assertIsInstance(results, pd.DataFrame)
            self.assertTrue(matched.empty)
            self.assertIsNone(matched_file)
            self.assertIsNone(translated)
            self.assertEqual(Path(result_file).name, "catalogue_matches.csv")
            self.assertEqual(
                sorted(path.name for path in runs_path.rglob("*") if path.is_file()),
                ["catalogue_matches.csv"],
            )

            published = "\n".join(
                [
                    str(summary),
                    results.to_csv(index=False),
                    Path(result_file).read_text(encoding="utf-8"),
                    captured_stdout.getvalue(),
                    captured_stderr.getvalue(),
                    repr(progress.calls),
                    repr(events),
                ]
            )
            self.assertNotIn(SECRET_MATCH, published)
            self.assertNotIn(SECRET_MISS, published)
            self.assertNotIn("internal-a", published)
            self.assertNotIn("epitope_query", published)
            self.assertNotIn("AB-200", published)
            self.assertIn("target_id,target_description", published)
            self.assertIn("protein-1", published)
            self.assertIn("Positive antibody–target matches: `1`", summary)

    def test_protected_nucleotide_run_keeps_genetic_code_and_translation(self) -> None:
        with tempfile.TemporaryDirectory() as directory:
            directory_path = Path(directory)
            catalogue_path = directory_path / "private-catalogue.csv"
            nucleotide_path = directory_path / "target.fna"
            runs_path = directory_path / "runs"
            catalogue_path.write_text(
                "antibody_id,manufacturer,catalog_number,antibody_name,"
                "epitope_sequence\n"
                "internal-mw,ExampleCo,AB-MW,Antibody MW,MW\n",
                encoding="utf-8",
            )
            nucleotide_path.write_text(
                ">target\nATATGAAGA\n",
                encoding="utf-8",
            )

            with (
                patch.dict(
                    os.environ,
                    {"PROTECTED_CATALOGUE_PATH": str(catalogue_path)},
                    clear=False,
                ),
                patch.object(app_gradio, "RUNS_DIR", runs_path),
                patch.object(app_gradio, "_SESSION_RUN_DIRS", set()),
            ):
                final = list(
                    app_gradio._run_scan(
                        AnalysisMode.PROTECTED_CATALOGUE.value,
                        None,
                        app_gradio.NUCLEOTIDE_MODE,
                        2,
                        str(nucleotide_path),
                        None,
                        None,
                        None,
                        ";",
                        progress=_ProgressRecorder(),
                    )
                )[-1]

            summary, results = final[0], final[1]
            self.assertIn("2 — Vertebrate Mitochondrial", summary)
            self.assertNotIn("perfect_match", results.columns)
            self.assertTrue(results["target_id"].iloc[0].startswith("target|frame"))
            self.assertEqual(Path(final[-1]).name, "output6frame.fasta")
            self.assertEqual(
                sorted(path.name for path in runs_path.rglob("*") if path.is_file()),
                ["catalogue_matches.csv", "output6frame.fasta"],
            )


if __name__ == "__main__":
    unittest.main()
