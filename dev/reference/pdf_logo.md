# Logo of an annotation source (CIViC, CGI, OncoKB, cancerhotspots.org) as inline Typst markup - images copied with the templates to the rendering directory (typst/logo\_\*.png); other sources as plain text

Logo of an annotation source (CIViC, CGI, OncoKB, cancerhotspots.org) as
inline Typst markup - images copied with the templates to the rendering
directory (typst/logo\_\*.png); other sources as plain text

## Usage

``` r
pdf_logo(
  id,
  height = "1.2em",
  sources = c("civic", "cgi", "oncokb", "cancerhotspots")
)
```

## Arguments

- id:

  source identifier(s) (e.g. 'civic', 'cgi', 'oncokb', 'cancerhotspots';
  NA as empty string)

- height:

  logo height (Typst length)

- sources:

  source identifiers with a logo image (typst/logo\_.png) in the
  rendering directory

## Value

character vector with Typst markup
