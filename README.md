# SixPack-AbScan

A command-line tool to predict antibody cross-reactivity on non-target species by exact epitope matching.

## What It Does

Given either:
- a nucleotide FASTA (genome/transcriptome), or
- a precomputed protein FASTA,

and a table of epitope peptide sequences, the tool:
1. Optionally generates a 6-frame translation (3 forward + 3 reverse-complement frames).
2. Scans translated/protein sequences for exact epitope matches.
3. Writes hit tables for downstream analysis.

## Limitations

- Only exact peptide matches are reported (no mismatches/partial matches).
- Only linear sequence matches are considered (no structural/conformational modeling).
- Results depend on the input assembly/transcriptome quality and completeness.

## Installation

```bash
python -m venv .venv
source .venv/bin/activate
pip install -r requirements.txt
```

## Usage

### 1) Start from nucleotide FASTA (build 6-frame translation first)

```bash
python sixpack_abscan.py \
  --input-nucleotide-fasta /path/to/transcriptome_or_genome.fasta \
  --epitope-file /path/to/epitopes.csv \
  --epitope-column epitope_specificity \
  --epitope-separator ';' \
  --output-dir /path/to/output
```

### 2) Use precomputed protein FASTA directly

```bash
python sixpack_abscan.py \
  --input-protein-fasta /path/to/proteome.fasta \
  --epitope-file /path/to/epitopes.xlsx \
  --epitope-column epitope_specificity \
  --output-dir /path/to/output
```

## Interactive App (Gradio)

Run:

```bash
python app_gradio.py
```

Then open the local URL printed in the terminal (typically `http://127.0.0.1:7860`).

In the app:
1. Choose input mode (`Nucleotide FASTA` or `Protein FASTA`).
2. Upload a plain/gzipped sequence FASTA and an epitope table, or paste a direct
   `https://*.ncbi.nlm.nih.gov/...` FASTA URL instead of uploading the sequence.
3. Set epitope column/separator if needed.
4. Run the scan, follow the preparation/translation/search progress, inspect the
   tables, and download output files.

Gradio displays its native browser-to-server transfer status during uploads. Once
an upload reaches the server, or while an NCBI file is downloaded, the app reports
progress for acquisition, gzip decompression, FASTA validation, translation, and
epitope scanning. Gzipped inputs are expanded into temporary local FASTA files
because the existing search pipeline expects a normal, seekable path; those
temporary files are removed after the run.

Remote inputs are restricted to HTTPS URLs on `ncbi.nlm.nih.gov` or its
subdomains. Redirect destinations are revalidated, DNS answers must be public,
and compressed/downloaded and decompressed sizes are bounded. The defaults can be
adjusted with these environment variables (values are bytes unless noted):

- `MAX_FASTA_SOURCE_BYTES` (default: 2 GiB)
- `MAX_FASTA_DECOMPRESSED_BYTES` (default: 8 GiB)
- `NCBI_DOWNLOAD_TIMEOUT_SECONDS` (default: 30)
- `MAX_NCBI_REDIRECTS` (default: 5)

The source-size limit is also passed to Gradio as its maximum upload size.

## Tests

Run the focused Gradio configuration, input-preparation, and URL security tests
with:

```bash
python -m unittest discover -v
```

## SciLifeLab Serve Preparation (Steps 0-4)

This repository now includes everything needed up to publish automation:

1. Gradio app entrypoint
   - `app_gradio.py`
   - Uses container-friendly launch settings (`HOST`, `PORT`, optional `ROOT_PATH`).
2. Containerization files
   - `Dockerfile`
   - `.dockerignore`
3. Local container test (optional)

```bash
docker build -t sixpack-abscan:local .
docker run --rm -p 7860:7860 sixpack-abscan:local
```

4. Automated image publish to GHCR
   - Workflow: `.github/workflows/publish-ghcr.yml`
   - Triggers on pushes to `main`, version tags (`v*`), and manual dispatch.
   - Publishes image tags to `ghcr.io/<owner>/<repo>`.

## Deploy To SciLifeLab Serve

After pushing this repository to GitHub and letting the GHCR workflow publish an image, create the app in Serve with:

- Application type: `Gradio`
- Container image: `ghcr.io/<owner>/<repo>:latest`
- Port: `7860`
- Health path: `/`

This app already supports the typical Serve runtime settings:

- `HOST=0.0.0.0`
- `PORT=7860`
- optional `ROOT_PATH` for reverse-proxy mounting

The container runs as a non-root user and stores temporary uploads and run outputs in writable app-owned directories.

## Outputs

Written to `--output-dir`:
- `output6frame.fasta` (only when `--input-nucleotide-fasta` is used)
- `epitope_hits.csv`: one row per matched `(epitope_query, target_id)`
- `matched_epitope_rows.csv`: original epitope metadata rows joined with matching targets

## Notes

- CSV/TSV epitope tables are read with `--epitope-separator` (default `;`).
- Excel input (`.xlsx`/`.xls`) is supported automatically.
- The default epitope column is `epitope_specificity`.

## Legacy Notebook

The original workflow notebook is still included:
- `6frame_translation and epitope search.ipynb`

The CLI is recommended for scripted/reproducible runs, and Gradio is recommended for interactive use.
