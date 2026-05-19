# Find transcripts covering given splice junction breakpoints

This function identifies transcripts that cover specified junction
breakpoints.

## Usage

``` r
bp_junction_transcript_overlap(bp_junctions = NULL, ref_data = NULL)
```

## Arguments

- bp_junctions:

  A data frame with columns:

  - BP_CHROM: Chromosome of the breakpoint junction

  - BP_POSITION: Chromosome position of the breakpoint junction

- ref_data:

  PCGR reference data bundle (list)

## Value

A data frame with columns:

- CHROM: Chromosome of the breakpoint junction

- BP_POSITION: Position of the breakpoint junction

- ENSEMBL_TRANSCRIPT_ID: Ensembl transcript ID covering the splice
  junction

- ENSEMBL_GENE_ID: Ensembl gene ID

- GENE_BIOTYPE: Biotype of the transcript

- TRANSCRIPT_START: Start position of the transcript

- TRANSCRIPT_END: End position of the transcript
