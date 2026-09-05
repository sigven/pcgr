# Pad segments narrower than a minimum width for genome-wide plot legibility

Widens (symmetrically, around the segment midpoint) any segment whose
genome-coordinate width is below `min_width_bp`, so it remains visible
as a rendered mark on a genome-wide plot. Padding is clamped to the
segment's own chromosome boundaries (`genome_start`/`genome_end`) so it
cannot visually bleed into a neighboring chromosome. Purely cosmetic:
does not touch SEGMENT_START/SEGMENT_END or any other reported column -
only the SegmentStart/SegmentEnd columns used for plot x-coordinates.

## Usage

``` r
widen_short_segments_for_plot(df, min_width_bp = 1e+07)
```

## Arguments

- df:

  data frame with SegmentStart, SegmentEnd, genome_start, genome_end

- min_width_bp:

  numeric minimum plotted width in genome-coordinate bp. Default 1e7 (10
  Mb, ~0.3% of the human genome).

## Value

df with SegmentStart/SegmentEnd widened where needed
