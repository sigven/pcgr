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
  target_size_mb,
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

- target_size_mb:

  size of targeted genomic region (coding)

- sample_name:

  name of sample

## Value

msi_data
