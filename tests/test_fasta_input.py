from __future__ import annotations

import gzip
import socket
import tempfile
import unittest
from io import BytesIO
from pathlib import Path
from typing import Self

from fasta_input import (
    InputPreparationError,
    download_ncbi_file,
    prepare_fasta,
    validate_ncbi_url,
)

FASTA = b">record-one\nACGTN\n>record-two description\nTTAA\n"
NCBI_URL = "https://ftp.ncbi.nlm.nih.gov/genomes/example.fna.gz"


def public_ncbi_resolver(*_args: object, **_kwargs: object) -> list[tuple]:
    return [
        (
            socket.AF_INET,
            socket.SOCK_STREAM,
            socket.IPPROTO_TCP,
            "",
            ("130.14.29.120", 443),
        )
    ]


class FakeResponse:
    def __init__(
        self,
        body: bytes = b"",
        *,
        status: int = 200,
        headers: dict[str, str] | None = None,
    ) -> None:
        self.status = status
        self.headers = headers or {}
        self._body = BytesIO(body)
        self.read_calls = 0

    def read(self, size: int = -1) -> bytes:
        self.read_calls += 1
        return self._body.read(size)

    def getcode(self) -> int:
        return self.status

    def __enter__(self) -> Self:
        return self

    def __exit__(self, *_args: object) -> None:
        return None


class FakeOpener:
    def __init__(self, *responses: FakeResponse) -> None:
        self.responses = list(responses)
        self.opened_urls: list[str] = []

    def open(self, request: object, *, timeout: int) -> FakeResponse:
        del timeout
        self.opened_urls.append(request.full_url)  # type: ignore[attr-defined]
        if not self.responses:
            raise AssertionError("No fake response remains")
        return self.responses.pop(0)


class PrepareFastaTests(unittest.TestCase):
    def test_plain_upload_is_validated_without_copying(self) -> None:
        with tempfile.TemporaryDirectory() as directory:
            source = Path(directory) / "input.fasta"
            source.write_bytes(FASTA)

            with prepare_fasta(uploaded_path=source) as prepared:
                self.assertEqual(prepared.path, source)
                self.assertEqual(prepared.record_count, 2)
                self.assertFalse(prepared.was_gzip)

    def test_gzip_upload_is_decompressed_to_temporary_fasta(self) -> None:
        with tempfile.TemporaryDirectory() as directory:
            source = Path(directory) / "input.fna.gz"
            source.write_bytes(gzip.compress(FASTA))

            prepared = prepare_fasta(uploaded_path=source)
            prepared_path = prepared.path
            self.assertNotEqual(prepared_path, source)
            self.assertEqual(prepared_path.read_bytes(), FASTA)
            self.assertEqual(prepared.record_count, 2)
            self.assertTrue(prepared.was_gzip)

            prepared.cleanup()
            self.assertFalse(prepared_path.exists())

    def test_decompression_limit_rejects_gzip_bomb(self) -> None:
        with tempfile.TemporaryDirectory() as directory:
            source = Path(directory) / "large.fasta.gz"
            source.write_bytes(gzip.compress(b">record\n" + b"A" * 5000 + b"\n"))

            with self.assertRaisesRegex(
                InputPreparationError, "decompressed FASTA exceeded"
            ):
                prepare_fasta(
                    uploaded_path=source,
                    max_source_bytes=10_000,
                    max_decompressed_bytes=100,
                )

    def test_url_gzip_download_uses_the_same_preparation_path(self) -> None:
        compressed = gzip.compress(FASTA)
        opener = FakeOpener(
            FakeResponse(
                compressed,
                headers={"Content-Length": str(len(compressed))},
            )
        )

        with prepare_fasta(
            ncbi_url=NCBI_URL,
            opener=opener,
            resolver=public_ncbi_resolver,
            max_source_bytes=10_000,
            max_decompressed_bytes=10_000,
        ) as prepared:
            self.assertEqual(prepared.path.read_bytes(), FASTA)
            self.assertEqual(prepared.record_count, 2)
            self.assertTrue(prepared.was_gzip)

        self.assertEqual(opener.opened_urls, [NCBI_URL])

    def test_download_progress_reports_bytes(self) -> None:
        response = FakeResponse(FASTA)
        opener = FakeOpener(response)
        updates: list[tuple[str, int | None, int | None]] = []

        with tempfile.TemporaryDirectory() as directory:
            destination = Path(directory) / "download"
            download_ncbi_file(
                "https://www.ncbi.nlm.nih.gov/data/example.fasta",
                destination,
                opener=opener,
                resolver=public_ncbi_resolver,
                max_bytes=10_000,
                progress=lambda stage, completed, total: updates.append(
                    (stage, completed, total)
                ),
            )

        self.assertEqual(updates[0], ("Downloading NCBI FASTA", 0, None))
        self.assertEqual(updates[-1], ("Downloading NCBI FASTA", len(FASTA), None))


class NcbiUrlSecurityTests(unittest.TestCase):
    def test_only_ncbi_https_urls_are_allowed(self) -> None:
        self.assertEqual(validate_ncbi_url(NCBI_URL), NCBI_URL)
        self.assertEqual(
            validate_ncbi_url("https://www.ncbi.nlm.nih.gov/file.fa"),
            "https://www.ncbi.nlm.nih.gov/file.fa",
        )

        rejected = [
            "http://ftp.ncbi.nlm.nih.gov/file.fa",
            "https://example.com/file.fa",
            "https://ftp.ncbi.nlm.nih.gov.evil.example/file.fa",
            "https://user@ftp.ncbi.nlm.nih.gov/file.fa",
            "https://ftp.ncbi.nlm.nih.gov:8443/file.fa",
            "https://127.0.0.1/file.fa",
        ]
        for url in rejected:
            with self.subTest(url=url), self.assertRaises(InputPreparationError):
                validate_ncbi_url(url)

    def test_redirect_to_non_ncbi_host_is_rejected_before_request(self) -> None:
        opener = FakeOpener(
            FakeResponse(status=302, headers={"Location": "https://example.com/a.fa"})
        )

        with tempfile.TemporaryDirectory() as directory:
            destination = Path(directory) / "download"
            with self.assertRaisesRegex(InputPreparationError, "Only HTTPS URLs"):
                download_ncbi_file(
                    NCBI_URL,
                    destination,
                    opener=opener,
                    resolver=public_ncbi_resolver,
                )
            self.assertFalse(destination.exists())

        self.assertEqual(opener.opened_urls, [NCBI_URL])

    def test_allowed_ncbi_redirect_is_revalidated_and_followed(self) -> None:
        redirected_url = "https://download.ncbi.nlm.nih.gov/data/example.fna"
        opener = FakeOpener(
            FakeResponse(status=302, headers={"Location": redirected_url}),
            FakeResponse(FASTA, headers={"Content-Length": str(len(FASTA))}),
        )

        with tempfile.TemporaryDirectory() as directory:
            destination = Path(directory) / "download"
            name = download_ncbi_file(
                NCBI_URL,
                destination,
                opener=opener,
                resolver=public_ncbi_resolver,
                max_bytes=10_000,
            )
            self.assertEqual(destination.read_bytes(), FASTA)
            self.assertEqual(name, "example.fna")

        self.assertEqual(opener.opened_urls, [NCBI_URL, redirected_url])

    def test_content_length_and_stream_are_both_size_limited(self) -> None:
        declared_too_large = FakeResponse(FASTA, headers={"Content-Length": "10001"})
        streamed_too_large = FakeResponse(FASTA)

        for response in (declared_too_large, streamed_too_large):
            with self.subTest(headers=response.headers):
                opener = FakeOpener(response)
                with tempfile.TemporaryDirectory() as directory:
                    destination = Path(directory) / "download"
                    with self.assertRaisesRegex(InputPreparationError, "limit"):
                        download_ncbi_file(
                            NCBI_URL,
                            destination,
                            opener=opener,
                            resolver=public_ncbi_resolver,
                            max_bytes=10,
                        )
                    self.assertFalse(destination.exists())

        self.assertEqual(declared_too_large.read_calls, 0)
        self.assertGreater(streamed_too_large.read_calls, 0)

    def test_private_dns_answer_is_rejected_before_request(self) -> None:
        opener = FakeOpener(FakeResponse(FASTA))

        def private_resolver(*_args: object, **_kwargs: object) -> list[tuple]:
            return [
                (
                    socket.AF_INET,
                    socket.SOCK_STREAM,
                    socket.IPPROTO_TCP,
                    "",
                    ("127.0.0.1", 443),
                )
            ]

        with (
            tempfile.TemporaryDirectory() as directory,
            self.assertRaisesRegex(InputPreparationError, "public IP"),
        ):
            download_ncbi_file(
                NCBI_URL,
                Path(directory) / "download",
                opener=opener,
                resolver=private_resolver,
            )

        self.assertEqual(opener.opened_urls, [])


if __name__ == "__main__":
    unittest.main()
