"""Securely prepare uploaded or remote FASTA inputs for local processing."""

from __future__ import annotations

import gzip
import ipaddress
import os
import socket
import tempfile
from collections.abc import Callable
from dataclasses import dataclass, field
from pathlib import Path
from typing import Self
from urllib.error import HTTPError, URLError
from urllib.parse import unquote, urljoin, urlsplit
from urllib.request import HTTPRedirectHandler, Request, build_opener

MIB = 1024 * 1024
GIB = 1024 * MIB
DOWNLOAD_CHUNK_SIZE = MIB
DECOMPRESSION_CHUNK_SIZE = MIB
VALIDATION_PROGRESS_INTERVAL = 4 * MIB
REDIRECT_STATUSES = {301, 302, 303, 307, 308}
ALLOWED_NCBI_DOMAIN = "ncbi.nlm.nih.gov"
ALLOWED_SEQUENCE_CHARACTERS = frozenset(
    b"ABCDEFGHIJKLMNOPQRSTUVWXYZabcdefghijklmnopqrstuvwxyz*-?.-"
)


def _positive_int_env(name: str, default: int) -> int:
    raw_value = os.getenv(name)
    if raw_value is None:
        return default
    try:
        value = int(raw_value)
    except ValueError as exc:
        raise RuntimeError(f"{name} must be a positive integer.") from exc
    if value <= 0:
        raise RuntimeError(f"{name} must be a positive integer.")
    return value


MAX_FASTA_SOURCE_BYTES = _positive_int_env("MAX_FASTA_SOURCE_BYTES", 2 * GIB)
MAX_FASTA_DECOMPRESSED_BYTES = _positive_int_env(
    "MAX_FASTA_DECOMPRESSED_BYTES", 8 * GIB
)
NCBI_DOWNLOAD_TIMEOUT_SECONDS = _positive_int_env("NCBI_DOWNLOAD_TIMEOUT_SECONDS", 30)
MAX_NCBI_REDIRECTS = _positive_int_env("MAX_NCBI_REDIRECTS", 5)

ProgressCallback = Callable[[str, int | None, int | None], None]


class InputPreparationError(ValueError):
    """Raised when a FASTA source cannot be acquired or safely prepared."""


@dataclass
class PreparedFasta:
    """A validated, uncompressed FASTA path and its owned temporary storage."""

    path: Path
    source_name: str
    record_count: int
    was_gzip: bool
    _temporary_directory: tempfile.TemporaryDirectory[str] | None = field(
        default=None, repr=False
    )

    def cleanup(self) -> None:
        if self._temporary_directory is not None:
            self._temporary_directory.cleanup()
            self._temporary_directory = None

    def __enter__(self) -> Self:
        return self

    def __exit__(self, *_args: object) -> None:
        self.cleanup()


class _NoRedirectHandler(HTTPRedirectHandler):
    """Expose redirects to the caller so every destination can be validated."""

    def redirect_request(self, *_args: object, **_kwargs: object) -> None:
        return None


def _report(
    progress: ProgressCallback | None,
    stage: str,
    completed: int | None = None,
    total: int | None = None,
) -> None:
    if progress is not None:
        progress(stage, completed, total)


def format_bytes(byte_count: int) -> str:
    value = float(byte_count)
    for unit in ("B", "KiB", "MiB", "GiB", "TiB"):
        if value < 1024 or unit == "TiB":
            precision = 0 if unit == "B" else 1
            return f"{value:.{precision}f} {unit}"
        value /= 1024
    return f"{value:.1f} TiB"


def validate_ncbi_url(url: str) -> str:
    """Return a stripped NCBI HTTPS URL or raise a user-facing error."""

    normalized_url = url.strip()
    if not normalized_url:
        raise InputPreparationError("Please provide an NCBI file URL.")
    if len(normalized_url) > 4096:
        raise InputPreparationError("The NCBI URL is too long.")

    parsed = urlsplit(normalized_url)
    if parsed.scheme.lower() != "https":
        raise InputPreparationError("NCBI file URLs must use HTTPS.")
    if parsed.username is not None or parsed.password is not None:
        raise InputPreparationError("NCBI file URLs cannot include credentials.")
    if parsed.fragment:
        raise InputPreparationError("NCBI file URLs cannot include a fragment.")
    if not parsed.hostname:
        raise InputPreparationError("The NCBI URL is missing a hostname.")

    try:
        hostname = parsed.hostname.encode("idna").decode("ascii").lower().rstrip(".")
        port = parsed.port
    except (UnicodeError, ValueError) as exc:
        raise InputPreparationError(
            "The NCBI URL has an invalid hostname or port."
        ) from exc

    if not (
        hostname == ALLOWED_NCBI_DOMAIN or hostname.endswith(f".{ALLOWED_NCBI_DOMAIN}")
    ):
        raise InputPreparationError(
            "Only HTTPS URLs hosted by ncbi.nlm.nih.gov are allowed."
        )
    if port not in (None, 443):
        raise InputPreparationError("NCBI file URLs may only use HTTPS port 443.")

    return normalized_url


def _assert_public_ncbi_resolution(
    url: str, resolver: Callable[..., object] = socket.getaddrinfo
) -> None:
    """Reject private/special DNS answers before opening an allowed NCBI URL."""

    hostname = urlsplit(url).hostname
    assert hostname is not None
    try:
        address_info = resolver(
            hostname,
            443,
            family=socket.AF_UNSPEC,
            type=socket.SOCK_STREAM,
        )
    except OSError as exc:
        raise InputPreparationError(
            f"Could not resolve the NCBI host '{hostname}'."
        ) from exc

    addresses: set[ipaddress.IPv4Address | ipaddress.IPv6Address] = set()
    for item in address_info:  # type: ignore[union-attr]
        try:
            addresses.add(ipaddress.ip_address(item[4][0]))
        except (IndexError, TypeError, ValueError):
            continue

    if not addresses or any(not address.is_global for address in addresses):
        raise InputPreparationError(
            "The NCBI hostname did not resolve exclusively to public IP addresses."
        )


def _response_status(response: object) -> int:
    status = getattr(response, "status", None)
    if status is None:
        status = response.getcode()  # type: ignore[attr-defined]
    return int(status)


def _open_url(opener: object, request: Request, timeout: int) -> object:
    try:
        return opener.open(request, timeout=timeout)  # type: ignore[attr-defined]
    except HTTPError as exc:
        if exc.code in REDIRECT_STATUSES:
            return exc
        raise InputPreparationError(
            f"NCBI returned HTTP {exc.code} while downloading the FASTA file."
        ) from exc
    except (TimeoutError, URLError, OSError) as exc:
        raise InputPreparationError(
            "The NCBI FASTA download failed. Please check the URL and try again."
        ) from exc


def _source_name_from_url(url: str) -> str:
    name = Path(unquote(urlsplit(url).path)).name
    return name[:255] if name else "ncbi_fasta"


def download_ncbi_file(
    url: str,
    destination: Path,
    *,
    progress: ProgressCallback | None = None,
    max_bytes: int = MAX_FASTA_SOURCE_BYTES,
    timeout: int = NCBI_DOWNLOAD_TIMEOUT_SECONDS,
    max_redirects: int = MAX_NCBI_REDIRECTS,
    opener: object | None = None,
    resolver: Callable[..., object] = socket.getaddrinfo,
) -> str:
    """Download an allowlisted NCBI URL without permitting unsafe redirects."""

    if max_bytes <= 0:
        raise ValueError("max_bytes must be positive.")
    current_url = validate_ncbi_url(url)
    url_opener = opener or build_opener(_NoRedirectHandler())

    try:
        for redirect_count in range(max_redirects + 1):
            current_url = validate_ncbi_url(current_url)
            _assert_public_ncbi_resolution(current_url, resolver)
            request = Request(
                current_url,
                headers={
                    "Accept": "application/octet-stream, text/plain;q=0.9, */*;q=0.1",
                    "User-Agent": "SixPack-AbScan/1.0",
                },
                method="GET",
            )
            response = _open_url(url_opener, request, timeout)
            with response:  # type: ignore[attr-defined]
                status = _response_status(response)
                headers = response.headers  # type: ignore[attr-defined]
                if status in REDIRECT_STATUSES:
                    location = headers.get("Location")
                    if not location:
                        raise InputPreparationError(
                            "NCBI returned a redirect without a destination."
                        )
                    if redirect_count >= max_redirects:
                        raise InputPreparationError(
                            "The NCBI URL exceeded the redirect limit."
                        )
                    current_url = validate_ncbi_url(urljoin(current_url, location))
                    continue
                if status < 200 or status >= 300:
                    raise InputPreparationError(
                        f"NCBI returned HTTP {status} while downloading the FASTA file."
                    )

                content_length_header = headers.get("Content-Length")
                content_length: int | None = None
                if content_length_header:
                    try:
                        content_length = int(content_length_header)
                    except ValueError as exc:
                        raise InputPreparationError(
                            "NCBI returned an invalid Content-Length header."
                        ) from exc
                    if content_length < 0:
                        raise InputPreparationError(
                            "NCBI returned an invalid Content-Length header."
                        )
                    if content_length > max_bytes:
                        raise InputPreparationError(
                            "The NCBI file is larger than the configured download limit "
                            f"({format_bytes(max_bytes)})."
                        )

                downloaded = 0
                _report(progress, "Downloading NCBI FASTA", 0, content_length)
                destination.parent.mkdir(parents=True, exist_ok=True)
                with destination.open("wb") as output_handle:
                    while True:
                        chunk = response.read(DOWNLOAD_CHUNK_SIZE)  # type: ignore[attr-defined]
                        if not chunk:
                            break
                        downloaded += len(chunk)
                        if downloaded > max_bytes:
                            raise InputPreparationError(
                                "The NCBI download exceeded the configured limit "
                                f"({format_bytes(max_bytes)})."
                            )
                        output_handle.write(chunk)
                        _report(
                            progress,
                            "Downloading NCBI FASTA",
                            downloaded,
                            content_length,
                        )

                if downloaded == 0:
                    raise InputPreparationError("NCBI returned an empty file.")
                return _source_name_from_url(current_url)
    except BaseException:
        destination.unlink(missing_ok=True)
        raise

    raise InputPreparationError("The NCBI URL exceeded the redirect limit.")


def _is_gzip(path: Path) -> bool:
    with path.open("rb") as handle:
        return handle.read(2) == b"\x1f\x8b"


def _decompress_gzip(
    source: Path,
    destination: Path,
    *,
    progress: ProgressCallback | None,
    max_decompressed_bytes: int,
) -> None:
    compressed_size = source.stat().st_size
    decompressed_size = 0
    _report(progress, "Decompressing FASTA", 0, compressed_size)
    try:
        with (
            source.open("rb") as raw_handle,
            gzip.GzipFile(fileobj=raw_handle, mode="rb") as gzip_handle,
            destination.open("wb") as output_handle,
        ):
            while True:
                chunk = gzip_handle.read(DECOMPRESSION_CHUNK_SIZE)
                if not chunk:
                    break
                decompressed_size += len(chunk)
                if decompressed_size > max_decompressed_bytes:
                    raise InputPreparationError(
                        "The decompressed FASTA exceeded the configured limit "
                        f"({format_bytes(max_decompressed_bytes)})."
                    )
                output_handle.write(chunk)
                _report(
                    progress,
                    "Decompressing FASTA",
                    min(raw_handle.tell(), compressed_size),
                    compressed_size,
                )
    except (gzip.BadGzipFile, EOFError, OSError) as exc:
        destination.unlink(missing_ok=True)
        raise InputPreparationError("The gzip FASTA is corrupt or incomplete.") from exc

    if decompressed_size == 0:
        destination.unlink(missing_ok=True)
        raise InputPreparationError("The gzip FASTA decompressed to an empty file.")


def validate_fasta(path: Path, *, progress: ProgressCallback | None = None) -> int:
    """Validate basic FASTA structure and return the number of records."""

    file_size = path.stat().st_size
    if file_size == 0:
        raise InputPreparationError("The FASTA file is empty.")

    record_count = 0
    current_record_has_sequence = False
    processed = 0
    next_progress_update = 0
    _report(progress, "Validating FASTA", 0, file_size)

    with path.open("rb") as handle:
        for line_number, raw_line in enumerate(handle, start=1):
            processed += len(raw_line)
            stripped = raw_line.strip()
            if not stripped:
                continue
            if b"\x00" in stripped:
                raise InputPreparationError(
                    f"The input is not a text FASTA file (line {line_number})."
                )

            if stripped.startswith(b">"):
                if record_count and not current_record_has_sequence:
                    raise InputPreparationError(
                        f"FASTA record {record_count} has no sequence."
                    )
                if not stripped[1:].strip():
                    raise InputPreparationError(
                        f"FASTA header on line {line_number} is empty."
                    )
                try:
                    stripped[1:].decode("utf-8")
                except UnicodeDecodeError as exc:
                    raise InputPreparationError(
                        f"FASTA header on line {line_number} is not valid UTF-8 text."
                    ) from exc
                record_count += 1
                current_record_has_sequence = False
            else:
                if record_count == 0:
                    raise InputPreparationError(
                        "The input is not FASTA: the first non-empty line must start with '>'."
                    )
                sequence = b"".join(stripped.split())
                if not sequence or any(
                    character not in ALLOWED_SEQUENCE_CHARACTERS
                    for character in sequence
                ):
                    raise InputPreparationError(
                        f"FASTA sequence line {line_number} contains invalid characters."
                    )
                current_record_has_sequence = True

            if processed >= next_progress_update:
                _report(progress, "Validating FASTA", processed, file_size)
                next_progress_update = processed + VALIDATION_PROGRESS_INTERVAL

    if record_count == 0:
        raise InputPreparationError("The input does not contain any FASTA records.")
    if not current_record_has_sequence:
        raise InputPreparationError(f"FASTA record {record_count} has no sequence.")

    _report(progress, "Validating FASTA", file_size, file_size)
    return record_count


def prepare_fasta(
    *,
    uploaded_path: str | Path | None = None,
    ncbi_url: str | None = None,
    progress: ProgressCallback | None = None,
    max_source_bytes: int = MAX_FASTA_SOURCE_BYTES,
    max_decompressed_bytes: int = MAX_FASTA_DECOMPRESSED_BYTES,
    download_timeout: int = NCBI_DOWNLOAD_TIMEOUT_SECONDS,
    max_redirects: int = MAX_NCBI_REDIRECTS,
    opener: object | None = None,
    resolver: Callable[..., object] = socket.getaddrinfo,
) -> PreparedFasta:
    """Acquire, decompress if needed, and validate one FASTA source."""

    remote_url = ncbi_url.strip() if ncbi_url else ""
    if bool(uploaded_path) == bool(remote_url):
        raise InputPreparationError(
            "Provide exactly one FASTA source: an uploaded file or an NCBI URL."
        )
    if max_source_bytes <= 0 or max_decompressed_bytes <= 0:
        raise ValueError("FASTA size limits must be positive.")

    temporary_directory: tempfile.TemporaryDirectory[str] | None = None
    try:
        if uploaded_path:
            source_path = Path(uploaded_path)
            if not source_path.is_file():
                raise InputPreparationError("The uploaded FASTA file is unavailable.")
            source_size = source_path.stat().st_size
            if source_size > max_source_bytes:
                raise InputPreparationError(
                    "The uploaded FASTA is larger than the configured limit "
                    f"({format_bytes(max_source_bytes)})."
                )
            if source_size == 0:
                raise InputPreparationError("The uploaded FASTA file is empty.")
            source_name = source_path.name
            _report(progress, "Upload received by server", source_size, source_size)
        else:
            temporary_directory = tempfile.TemporaryDirectory(
                prefix="sixpack_abscan_fasta_"
            )
            source_path = Path(temporary_directory.name) / "downloaded_input"
            assert remote_url
            source_name = download_ncbi_file(
                remote_url,
                source_path,
                progress=progress,
                max_bytes=max_source_bytes,
                timeout=download_timeout,
                max_redirects=max_redirects,
                opener=opener,
                resolver=resolver,
            )

        is_gzip = _is_gzip(source_path)
        if source_name.lower().endswith(".gz") and not is_gzip:
            raise InputPreparationError(
                "The file name ends in .gz, but its contents are not valid gzip data."
            )

        prepared_path = source_path
        if is_gzip:
            if temporary_directory is None:
                temporary_directory = tempfile.TemporaryDirectory(
                    prefix="sixpack_abscan_fasta_"
                )
            prepared_path = Path(temporary_directory.name) / "prepared.fasta"
            _decompress_gzip(
                source_path,
                prepared_path,
                progress=progress,
                max_decompressed_bytes=max_decompressed_bytes,
            )

        record_count = validate_fasta(prepared_path, progress=progress)
        return PreparedFasta(
            path=prepared_path,
            source_name=source_name,
            record_count=record_count,
            was_gzip=is_gzip,
            _temporary_directory=temporary_directory,
        )
    except BaseException:
        if temporary_directory is not None:
            temporary_directory.cleanup()
        raise
