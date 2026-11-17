# SixPack-AbScan user guide

SixPack-AbScan predicts possible antibody cross-reactivity by finding **exact,
linear peptide matches** between antibody epitopes and target sequences. It can
scan a protein FASTA directly or translate every nucleotide FASTA record in all
six reading frames before scanning.

## Choose an analysis mode

- **Upload my epitope list**: upload your own epitope table. The app returns the
  matched peptide, target identifier, target description, and matched input
  metadata.
- **Scan a pre-loaded commercial antibody catalogue**: choose a catalogue held
  on the server. To protect licensed epitope data, the app returns only positive
  antibody-to-target matches and does not reveal epitope sequences or positions.

## Prepare the inputs

### Sequence input

Provide exactly one sequence source:

1. Upload an uncompressed or gzipped FASTA (`.fa`, `.fasta`, `.fna`, `.faa`,
   `.fas`, or `.gz`), or
2. Paste a direct HTTPS URL to a FASTA file hosted on `ncbi.nlm.nih.gov`.

Select **Nucleotide FASTA** when the input is DNA. The app translates the three
forward and three reverse-complement frames. Select the NCBI genetic code that
matches the organism or organelle; table 1 is the standard default.

Select **Protein FASTA** for an already translated proteome or protein set. No
translation is performed in this mode.

### Epitope table

For an uploaded epitope list, use CSV, TSV, or Excel format. The table may
contain any metadata columns, but one column must contain the peptide sequences.
After upload, select that column in **Epitope column**. For CSV/TSV files, enter
the delimiter used by the file; the bundled DSHB example uses a comma (`,`).

Epitope values are stripped of spaces and converted to uppercase. Blank cells
are ignored. A query is reported when its complete normalized sequence occurs
in a target sequence.

## Run and interpret a scan

1. Select an analysis mode and input type.
2. Provide the required files or NCBI URL.
3. Check the epitope column, separator, and (for DNA) genetic code.
4. Select **Run Scan** and wait for the preparation, translation, and scan
   progress to complete.
5. Inspect the result tables and download the generated files.

The uploaded-list workflow produces:

- `epitope_hits.csv`: one row for each matched epitope and target FASTA record.
- `matched_epitope_rows.csv`: matching input metadata joined to the target.
- `output6frame.fasta`: the translated sequences, only for nucleotide input.

The protected-catalogue workflow produces `catalogue_matches.csv`, containing
only permitted catalogue metadata and the matched target identifiers.

Multiple rows can be biologically meaningful: the same peptide may occur in
several proteins, isoforms, or reading frames. A missing match means only that
the exact linear peptide was not present in the supplied sequences.

## Examples and tutorial

The **Tutorial & example downloads** panel provides one preconfigured example.
**Load demo data** loads a curated three-record DSHB subset and fills a direct,
version-pinned NCBI URL for the complete `GCF_000001215.4` reference proteome.
The scan produces 13 isoform-level hits.

After loading the example, select **Run Scan**. The
[worked tutorial](tutorial.md) documents every setting and expected result.
No download is required for the preconfigured demo. For optional manual testing,
the panel provides the curated three-record DSHB subset and the expected 13-row
hit table.

## Method and limitations

SixPack-AbScan is a sequence-screening tool, not an experimental validation of
antibody specificity. It does not model substitutions, gapped or partial
matches, post-translational modifications, three-dimensional/conformational
epitopes, epitope accessibility, or expression. Results also depend on the
quality and completeness of the input assembly or proteome. Candidate matches
should therefore be evaluated using biological context and, where appropriate,
experimental validation.

## Privacy and large inputs

Uploaded files are processed for the current run and temporary run directories
are cleaned by the configured deployment lifecycle. Do not upload confidential
data unless the operator's data-handling policy permits it. For sequence files
larger than roughly 2–3 GB, use the command-line program or a local deployment.

## License and citation

The software is distributed under the [MIT License](../LICENSE). The public web
service is available free of charge for non-commercial use. The MIT license also
permits reuse, modification, distribution, and commercial use of the software,
subject to its copyright and notice conditions.

The public DSHB catalogue and derived example subset are research data, not
software, and are not covered by the MIT license. Their provenance and reuse
notice are documented in the [example data README](../examples/README.md).

When using SixPack-AbScan in published work, please cite:

> Grillo, M. (2026). SixPack-AbScan: predict crossreactivity of monoclonal
> antibodies on non-target species. https://doi.org/10.82595/scilifelab.4c3a-tp57
