# Number of Mitelman DB records of a fusion in a given tumor site, parsed from the annotated Mitelman evidence ('**Site** (*n = 3*) - ...')

Number of Mitelman DB records of a fusion in a given tumor site, parsed
from the annotated Mitelman evidence ('**Site** (*n = 3*) - ...')

## Usage

``` r
mitdb_site_count(evidence, site)
```

## Arguments

- evidence:

  Mitelman evidence per fusion (MITDB_EVIDENCE)

- site:

  primary tumor site (e.g. 'Colon/Rectum')

## Value

integer vector (0 when the site is not recorded)
