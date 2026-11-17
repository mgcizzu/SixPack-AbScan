# SixPack-AbScan example data

These files support the in-app Drosophila proteome demo and the worked tutorial.

## Demo files

- `dshb_demo_subset.csv` contains three DSHB antibodies with distinct,
  interpretable exact matches in the Drosophila reference proteome.
- `expected_dmel_proteome_hits.csv` is the 13-row result expected from scanning
  the complete version-pinned proteome with `epitope_sequence` as the epitope
  column and comma as the separator.

The demo downloads the complete version-pinned NCBI RefSeq proteome for assembly
`GCF_000001215.4`:
[NCBI](https://ftp.ncbi.nlm.nih.gov/genomes/all/GCF/000/001/215/GCF_000001215.4_Release_6_plus_ISO1_MT/GCF_000001215.4_Release_6_plus_ISO1_MT_protein.faa.gz).

## Complete public DSHB catalogue

`dshb_updated_2026-09-10.csv` is the complete 49-record catalogue supplied for
public use with SixPack-AbScan. `dshb_demo_subset.csv` is a three-record subset
of that file; no epitope sequence was modified.

The catalogue data originate from the
[Developmental Studies Hybridoma Bank (DSHB)](https://dshb.biology.uiowa.edu/).
They are distributed here as research data and are not covered by the software's
MIT license. Users should retain the source attribution and consult the source
provider's applicable terms when redistributing or reusing the data.

## Scope

The examples demonstrate software behavior, not experimentally validated
cross-reactivity. Exact sequence identity does not establish antibody binding.
