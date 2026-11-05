# Function that initiates report element with HRD results

The HRD report element holds the result of an external HRD assay/tool
('external', `--input_hrd`). The genomic scar score estimated by PCGR
from the copy number segments (`--estimate_hrd`) is kept in the settings
(`settings$conf$somatic_cna`); 'internal_eval' indicates whether it is
available.

## Usage

``` r
init_hrd_content()
```
