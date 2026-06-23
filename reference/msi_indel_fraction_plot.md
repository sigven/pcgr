# Function that plots the indel fraction for a given sample and contrasts this with the distribution for MSI-H/MSS samples from TCGA

Function that plots the indel fraction for a given sample and contrasts
this with the distribution for MSI-H/MSS samples from TCGA

## Usage

``` r
msi_indel_fraction_plot(
  tcga_msi_dataset,
  indel_fraction,
  color_palette = pcgrr::color_palette
)
```

## Arguments

- tcga_msi_dataset:

  underlying dataset from TCGA used for development of statistical
  classifier

- indel_fraction:

  fraction of indels of all mutations (SNVs + indels)

## Value

p
