# Internal evidence tiers of prognostic and diagnostic biomarker evidence shown in the reports (HTML/PDF sections, key findings)

Prognostic and diagnostic evidence is not tiered in the reports
(AMP/ASCO/CAP tiers refer to therapeutic evidence only); the internal
tiers are used to select which evidence is shown. With a primary site,
evidence matching the primary site is shown (internal tiers 1/2:
strong/weak evidence). Without a primary site (site 'Any'), evidence
from any tumor type is considered, for which weak evidence (levels
C/D/E) always gets internal tier 3 - which is then also shown.

## Usage

``` r
progn_diagn_tiers(report)
```

## Arguments

- report:

  PCGR report object

## Value

integer vector of tiers
