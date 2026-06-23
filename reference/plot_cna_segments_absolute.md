# Plot allele-specific copy number segments (absolute copies)

Function that plots allele-specific copy number segments (minor + total
allele copies)

## Usage

``` r
plot_cna_segments_absolute(
  chrom_coordinates = NULL,
  cna_segment = NULL,
  cna_gene = NULL,
  tumor_ploidy = 2,
  amp_threshold_effective = 5,
  threshold_mode = "absolute",
  gain_threshold_effective = 3,
  del_threshold_effective = 1,
  color_palette = pcgrr::color_palette
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

- tumor_ploidy:

  numeric tumor ploidy used as the neutral baseline. A dotted reference
  line is always drawn at this value. Default 2.

- amp_threshold_effective:

  numeric effective amplification threshold in absolute copy number
  units. A dotted reference line is drawn at this value when
  `threshold_mode` is "absolute" or "combined". Default 5.

- threshold_mode:

  character thresholding mode applied to all CNA tiers: "absolute",
  "relative", or "combined". The amplification threshold line is shown
  only when mode is "absolute" or "combined". Default "absolute".

- gain_threshold_effective:

  numeric effective gain threshold in absolute copy number units. A
  dotted reference line is drawn at this value. Default 3.

- del_threshold_effective:

  numeric effective heterozygous deletion threshold in absolute copy
  number units. A dotted reference line is drawn at this value. Default
  1.
