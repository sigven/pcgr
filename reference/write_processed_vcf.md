# Function that writes a VCF intended for mutational signature analysis

Function that writes a VCF intended for mutational signature analysis

## Usage

``` r
write_processed_vcf(
  calls,
  sample_name = NULL,
  output_directory = NULL,
  vcf_fname = NULL,
  snv_only = TRUE
)
```

## Arguments

- calls:

  data frame with calls

- sample_name:

  sample name

- output_directory:

  Output directory for output file

- vcf_fname:

  filename for VCF

- snv_only:

  logical, if TRUE only SNVs are written to VCF
