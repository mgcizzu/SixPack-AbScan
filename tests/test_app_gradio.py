from __future__ import annotations

import os
import unittest

os.environ["CLEANUP_RUNS_ON_START"] = "0"
os.environ["CLEANUP_RUNS_ON_EXIT"] = "0"

from app_gradio import build_app  # noqa: E402


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
        self.assertEqual(len(run_dependency["inputs"]), 6)

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


if __name__ == "__main__":
    unittest.main()
