"""Backend-only catalogue loading and restricted exact-match scanning.

The private sequence is used only while the server scans the target. Returned
rows, exceptions, progress values, and catalogue result files never contain it.
"""

from __future__ import annotations

import csv
import hashlib
import os
import re
from collections.abc import Generator
from dataclasses import dataclass
from functools import lru_cache
from pathlib import Path

import pandas as pd
from Bio import SeqIO

CATALOGUE_PATH_ENV = "PROTECTED_CATALOGUE_PATH"
CATALOGUE_DIRECTORY_ENV = "PROTECTED_CATALOGUE_DIR"
CATALOGUE_LABEL_ENV = "PROTECTED_CATALOGUE_LABEL"
DEFAULT_CATALOGUE_DIRECTORY = Path("/srv/project_vol/catalogues")
PRIVATE_COLUMNS = (
    "antibody_id",
    "manufacturer",
    "catalog_number",
    "antibody_name",
    "epitope_sequence",
)
RESTRICTED_RESULT_COLUMNS = (
    "manufacturer",
    "catalog_number",
    "antibody_name",
    "target_id",
    "target_description",
)
MAX_SCAN_PROGRESS_UPDATES = 100
_AMINO_ACID_SEQUENCE = re.compile(r"[A-Z]+")
_SPREADSHEET_FORMULA_PREFIXES = ("=", "+", "-", "@")


class ProtectedCatalogueError(RuntimeError):
    """Base class for safe, user-displayable catalogue errors."""


class ProtectedCatalogueConfigurationError(ProtectedCatalogueError):
    """Raised when the server has no usable catalogue configured."""


class ProtectedCatalogueValidationError(ProtectedCatalogueError):
    """Raised when private catalogue structure or values are invalid."""


@dataclass(frozen=True, slots=True)
class ProtectedCatalogueEntry:
    """One private query plus the metadata permitted in results."""

    antibody_id: str
    manufacturer: str
    catalog_number: str
    antibody_name: str
    epitope_sequence: str


@dataclass(frozen=True, slots=True)
class ProtectedCatalogue:
    """An immutable catalogue loaded only by the backend."""

    entries: tuple[ProtectedCatalogueEntry, ...]


@dataclass(frozen=True, slots=True)
class ProtectedCatalogueOption:
    """A server-side catalogue exposed to the UI only by label and opaque key."""

    key: str
    label: str
    path: Path


def _safe_metadata(value: object, *, field: str, row_number: int) -> str:
    cleaned = str(value or "").strip()
    if not cleaned:
        raise ProtectedCatalogueValidationError(
            f"Protected catalogue row {row_number} has an empty '{field}' field."
        )
    if "\n" in cleaned or "\r" in cleaned:
        raise ProtectedCatalogueValidationError(
            f"Protected catalogue row {row_number} has invalid public metadata."
        )
    if cleaned.startswith(_SPREADSHEET_FORMULA_PREFIXES):
        cleaned = f"'{cleaned}"
    return cleaned


def _private_sequence(value: object, *, row_number: int) -> str:
    sequence = str(value or "").strip().upper().replace(" ", "")
    if not sequence or not _AMINO_ACID_SEQUENCE.fullmatch(sequence):
        raise ProtectedCatalogueValidationError(
            f"Protected catalogue row {row_number} has an invalid epitope sequence."
        )
    return sequence


def _safe_target_metadata(value: object) -> str:
    """Make user-supplied FASTA metadata safe for CSV display/download."""

    cleaned = str(value or "").strip().replace("\r", " ").replace("\n", " ")
    if cleaned.startswith(_SPREADSHEET_FORMULA_PREFIXES):
        cleaned = f"'{cleaned}"
    return cleaned


def load_protected_catalogue(path: Path) -> ProtectedCatalogue:
    """Load a private CSV/TSV without repeating values in error messages."""

    try:
        delimiter = "\t" if path.suffix.lower() == ".tsv" else ","
        with path.open("r", encoding="utf-8-sig", newline="") as handle:
            reader = csv.DictReader(handle, delimiter=delimiter)
            fieldnames = set(reader.fieldnames or ())
            missing = [column for column in PRIVATE_COLUMNS if column not in fieldnames]
            if missing:
                raise ProtectedCatalogueValidationError(
                    "Protected catalogue is missing required columns: "
                    + ", ".join(missing)
                )

            entries: list[ProtectedCatalogueEntry] = []
            seen_ids: set[str] = set()
            for row_number, row in enumerate(reader, start=2):
                antibody_id = _safe_metadata(
                    row.get("antibody_id"),
                    field="antibody_id",
                    row_number=row_number,
                )
                if antibody_id in seen_ids:
                    raise ProtectedCatalogueValidationError(
                        f"Protected catalogue row {row_number} repeats an antibody ID."
                    )
                seen_ids.add(antibody_id)
                entries.append(
                    ProtectedCatalogueEntry(
                        antibody_id=antibody_id,
                        manufacturer=_safe_metadata(
                            row.get("manufacturer"),
                            field="manufacturer",
                            row_number=row_number,
                        ),
                        catalog_number=_safe_metadata(
                            row.get("catalog_number"),
                            field="catalog_number",
                            row_number=row_number,
                        ),
                        antibody_name=_safe_metadata(
                            row.get("antibody_name"),
                            field="antibody_name",
                            row_number=row_number,
                        ),
                        epitope_sequence=_private_sequence(
                            row.get("epitope_sequence"), row_number=row_number
                        ),
                    )
                )
    except ProtectedCatalogueError:
        raise
    except (OSError, UnicodeError, csv.Error) as exc:
        raise ProtectedCatalogueConfigurationError(
            "The protected catalogue could not be loaded by the server."
        ) from exc

    if not entries:
        raise ProtectedCatalogueValidationError(
            "The protected catalogue contains no antibody records."
        )
    return ProtectedCatalogue(entries=tuple(entries))


@lru_cache(maxsize=8)
def _load_cached_catalogue(resolved_path: str) -> ProtectedCatalogue:
    return load_protected_catalogue(Path(resolved_path))


def _catalogue_key(path: Path) -> str:
    return hashlib.sha256(str(path).encode("utf-8")).hexdigest()[:16]


def _catalogue_label(path: Path, explicit_label: str | None = None) -> str:
    label = (explicit_label or path.stem.replace("_", " ").replace("-", " ")).strip()
    if not label or "\n" in label or "\r" in label:
        raise ProtectedCatalogueConfigurationError(
            "A protected catalogue has an invalid public label."
        )
    return label if explicit_label else label.title()


def get_configured_catalogue_options() -> tuple[ProtectedCatalogueOption, ...]:
    """Resolve configured private files to public labels and opaque UI keys."""

    candidates: list[tuple[Path, str | None]] = []
    configured_directory = (os.getenv(CATALOGUE_DIRECTORY_ENV) or "").strip()
    directory_path = (
        Path(configured_directory).expanduser()
        if configured_directory
        else DEFAULT_CATALOGUE_DIRECTORY
    )
    if configured_directory or directory_path.is_dir():
        try:
            directory = directory_path.resolve(strict=True)
            if not directory.is_dir():
                raise OSError
            candidates.extend(
                (path.resolve(strict=True), None)
                for path in sorted(
                    directory.iterdir(), key=lambda item: item.name.lower()
                )
                if path.is_file() and path.suffix.lower() in {".csv", ".tsv"}
            )
        except OSError as exc:
            raise ProtectedCatalogueConfigurationError(
                "The protected catalogue directory could not be loaded by the server."
            ) from exc

    configured_path = (os.getenv(CATALOGUE_PATH_ENV) or "").strip()
    if configured_path:
        try:
            path = Path(configured_path).expanduser().resolve(strict=True)
            if not path.is_file() or path.suffix.lower() not in {".csv", ".tsv"}:
                raise OSError
            candidates.append(
                (path, (os.getenv(CATALOGUE_LABEL_ENV) or "").strip() or None)
            )
        except OSError as exc:
            raise ProtectedCatalogueConfigurationError(
                "The protected catalogue could not be loaded by the server."
            ) from exc

    unique_candidates: dict[Path, str | None] = {}
    for path, explicit_label in candidates:
        unique_candidates[path] = explicit_label or unique_candidates.get(path)

    options = tuple(
        ProtectedCatalogueOption(
            key=_catalogue_key(path),
            label=_catalogue_label(path, explicit_label),
            path=path,
        )
        for path, explicit_label in unique_candidates.items()
    )
    labels = [option.label.casefold() for option in options]
    if len(labels) != len(set(labels)):
        raise ProtectedCatalogueConfigurationError(
            "Protected catalogue public labels must be unique."
        )
    return options


def get_configured_catalogue_choices() -> list[tuple[str, str]]:
    """Return public dropdown labels paired with opaque selection keys."""

    return [(option.label, option.key) for option in get_configured_catalogue_options()]


def get_configured_catalogue_paths() -> list[str]:
    """Return every private path that Gradio must block from file serving."""

    return [str(option.path) for option in get_configured_catalogue_options()]


def get_configured_catalogue(catalogue_key: str | None = None) -> ProtectedCatalogue:
    """Return the selected runtime catalogue, cached for the server process."""

    options = get_configured_catalogue_options()
    if not options:
        raise ProtectedCatalogueConfigurationError(
            "Protected catalogue mode is not configured on this server."
        )
    if catalogue_key is None and len(options) == 1:
        selected = options[0]
    else:
        selected = next(
            (option for option in options if option.key == catalogue_key), None
        )
    if selected is None:
        raise ProtectedCatalogueConfigurationError(
            "Please select an available protected antibody catalogue."
        )
    return _load_cached_catalogue(str(selected.path))


def scan_protected_catalogue_with_progress(
    catalogue: ProtectedCatalogue,
    protein_fasta: Path,
    *,
    total_records: int | None = None,
) -> Generator[tuple[int, int, int], None, pd.DataFrame]:
    """Return one restricted row per positive antibody-target record match."""

    if total_records is None:
        with protein_fasta.open("r", encoding="utf-8") as count_handle:
            total_records = sum(1 for _ in SeqIO.parse(count_handle, "fasta"))

    restricted_rows: list[dict[str, str]] = []
    update_interval = max(1, total_records // MAX_SCAN_PROGRESS_UPDATES)
    processed_records = 0
    last_reported = 0

    with protein_fasta.open("r", encoding="utf-8") as protein_handle:
        for processed_records, record in enumerate(
            SeqIO.parse(protein_handle, "fasta"), start=1
        ):
            sequence = str(record.seq).upper()
            target_id = _safe_target_metadata(record.id)
            target_description = _safe_target_metadata(record.description)
            for entry in catalogue.entries:
                if entry.epitope_sequence in sequence:
                    restricted_rows.append(
                        {
                            "manufacturer": entry.manufacturer,
                            "catalog_number": entry.catalog_number,
                            "antibody_name": entry.antibody_name,
                            "target_id": target_id,
                            "target_description": target_description,
                        }
                    )

            if (
                processed_records == 1
                or processed_records == total_records
                or processed_records % update_interval == 0
            ):
                last_reported = processed_records
                yield processed_records, total_records, len(restricted_rows)

    if processed_records and last_reported != processed_records:
        yield processed_records, total_records, len(restricted_rows)

    return (
        pd.DataFrame(restricted_rows, columns=RESTRICTED_RESULT_COLUMNS)
        .sort_values(
            ["catalog_number", "antibody_name", "target_id", "target_description"],
            kind="stable",
        )
        .reset_index(drop=True)
    )
