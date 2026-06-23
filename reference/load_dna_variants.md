# Function that reads and validates CNA or SNV/InDel TSV files file from PCGR/CPSR pre-report (Python) pipeline

Function that reads and validates CNA or SNV/InDel TSV files file from
PCGR/CPSR pre-report (Python) pipeline

## Usage

``` r
load_dna_variants(
  fname = NA,
  cols = NULL,
  ref_data = NULL,
  settings = NULL,
  vartype = "snv_indel",
  primary_site = "Any",
  retained_info_tags = "None",
  variant_origin = "Somatic"
)
```

## Arguments

- fname:

  Path to raw input file with DNA aberrations (PCGR/CPSR)

- cols:

  column type definitions of raw input file

- ref_data:

  PCGR reference data object

- settings:

  PCGR run/configuration settings

- vartype:

  type of DNA aberrations ('snv_indel','cna')

- primary_site:

  primary site of tumor

- retained_info_tags:

  VCF INFO tags to be retained in output (SNVs/InDels)

- variant_origin:

  Germline/Somatic
