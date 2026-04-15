# Plot allele-specific copy number segments

Function that plots allele-specific copy number segments (minor + total
allele copies)

## Usage

``` r
plot_cna_segments(
  chrom_coordinates = NULL,
  cna_segment = NULL,
  cna_gene = NULL
)
```

## Arguments

- chrom_coordinates:

  data frame with assembly-specific chromosome coordinate data (length
  etc)

- cna_segment:

  data frame with annotated copy number segments

- cna_gene:

  data frame with gene-level copy number data
