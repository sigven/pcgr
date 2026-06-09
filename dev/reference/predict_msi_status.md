# Function that predicts MSI status based on fraction of indels among calls

Function that predicts MSI status based on fraction of indels among
calls

## Usage

``` r
predict_msi_status(
  variant_set,
  ref_data,
  msi_prediction_model,
  msi_prediction_dataset,
  msi_feature_variance_table = NULL,
  target_size_mb,
  n_calls = NA_integer_,
  sample_name = "Test"
)
```

## Arguments

- variant_set:

  data frame with somatic mutations/indels

- ref_data:

  PCGR reference data object

- msi_prediction_model:

  statistical model for MSI prediction

- msi_prediction_dataset:

  underlying dataset from TCGA used for development of statistical
  classifier

- msi_feature_variance_table:

  lookup table (from model training) mapping mutation-count bins to
  expected feature SD; used to qualify prediction confidence. NULL if
  not available.

- target_size_mb:

  size of targeted genomic region (coding)

- n_calls:

  number of variants passed to the classifier

- sample_name:

  name of sample

## Value

msi_data
