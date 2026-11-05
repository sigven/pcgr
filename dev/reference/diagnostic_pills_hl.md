# Diagnostic hits to highlight as pills in the key findings

Only when a primary tumor site is specified: alterations (SNVs/InDels,
copy number aberrations, RNA fusions) with *strong* diagnostic evidence
(CIViC levels A/B, OncoKB Dx1) matching that site. With a site set, the
diagnostic evidence is restricted to the site (see
assign_amp_asco_cap_tiers()), so the internal evidence tier 1 is exactly
strong, site-matching evidence. Nothing is returned for site 'Any' (code
0), where evidence from any tumor type is considered - such findings are
only listed among the supporting findings.

## Usage

``` r
diagnostic_pills_hl(report, max_pills = 3)
```

## Arguments

- report:

  PCGR report object

- max_pills:

  maximum number of pills (remaining hits stay in the supporting
  findings)

## Value

list with 'pills' (plain-text pill labels, e.g. 'Diagnostic: EWSR1::FLI1
fusion (Ewing sarcoma)') and 'labels' (the corresponding markdown labels
of tiered_alterations_hl(), used to avoid listing the same alterations
again among the supporting findings)
