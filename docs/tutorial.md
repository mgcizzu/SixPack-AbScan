# Tutorial: DSHB antibodies against the Drosophila proteome

This tutorial demonstrates the complete web workflow using three public DSHB
antibody records and the version-pinned *Drosophila melanogaster* NCBI RefSeq
proteome. It covers loading the example, running an exact-match scan,
interpreting the result, and downloading the output.

## Example inputs

- [`dshb_demo_subset.csv`](../examples/dshb_demo_subset.csv) contains three
  records selected from the public DSHB catalogue.
- The target is the complete NCBI RefSeq proteome for assembly
  `GCF_000001215.4`, approximately 9.8 MB compressed:
  [`GCF_000001215.4_Release_6_plus_ISO1_MT_protein.faa.gz`](https://ftp.ncbi.nlm.nih.gov/genomes/all/GCF/000/001/215/GCF_000001215.4_Release_6_plus_ISO1_MT/GCF_000001215.4_Release_6_plus_ISO1_MT_protein.faa.gz)
- [`expected_dmel_proteome_hits.csv`](../examples/expected_dmel_proteome_hits.csv)
  contains the expected 13-row result.

## Run the demo in the web app

1. Open **Tutorial & example downloads** and select **Load demo data**.
2. Confirm that the form now shows:
   - Analysis mode: **Upload my epitope list**
   - Input type: **Protein FASTA**
   - NCBI URL ending in `GCF_000001215.4_Release_6_plus_ISO1_MT_protein.faa.gz`
   - Epitope column: `epitope_sequence`
   - CSV/TSV separator: `,`
3. Select **Run Scan**.
4. Wait while the app downloads, validates, and decompresses the proteome.
5. Inspect the two result tables and download the CSV outputs.

No manual download or upload is required.

## Expected result

The scan returns 13 rows in `epitope_hits.csv`:

| DSHB antibody | Exact epitope | Matching protein records |
| --- | --- | ---: |
| `4D9 anti-engrailed/invected` | `ELGLNEAQIKI` | 5 |
| `KLAR-N 9E1B` | `DPKELRRWLREMENKLESSPTITELTLFGNAELQRHLAVHS` | 4 |
| `Seven-up 2D3` | `FDIGLSRSMGLGVPPH` | 4 |

The multiple rows represent annotated protein isoforms or the related
engrailed/invected proteins. The exact accessions and descriptions are listed in
[`expected_dmel_proteome_hits.csv`](../examples/expected_dmel_proteome_hits.csv).

`matched_epitope_rows.csv` also contains 13 rows. It joins each hit back to its
DSHB catalogue number, antibody name, manufacturer, and original epitope
sequence.

## Optional manual-upload files

No download or manual upload is required for the preconfigured demo. To test
the upload workflow manually, use the two optional downloads in **Tutorial &
example downloads**:

- `dshb_demo_subset.csv`: upload this as the epitope table.
- `expected_dmel_proteome_hits.csv`: use this to compare the scan result.

For the uploaded DSHB subset, use:

- Epitope column: `epitope_sequence`
- Separator: `,`

The three-record subset keeps the tutorial result compact and interpretable.

## Command-line equivalent

Download and decompress the same version-pinned proteome, then run:

```bash
python sixpack_abscan.py \
  --input-protein-fasta GCF_000001215.4_Release_6_plus_ISO1_MT_protein.faa \
  --epitope-file examples/dshb_demo_subset.csv \
  --epitope-column epitope_sequence \
  --epitope-separator ',' \
  --output-dir tutorial_output
```

Compare `tutorial_output/epitope_hits.csv` with
[`examples/expected_dmel_proteome_hits.csv`](../examples/expected_dmel_proteome_hits.csv).

## Interpretation

These results establish exact linear sequence identity only. They do not
establish that an antibody binds a target experimentally. Epitope accessibility,
protein expression, sample preparation, conformational effects, and other
biological factors remain outside the scope of the analysis.
