# Gene symbols as links to the Open Targets Platform target page (by Ensembl gene id) in a PDF report table cell (raw Typst markup; the column must not be escaped by format_tt()). Symbols without a valid Ensembl gene id are shown as plain (escaped) text.

Gene symbols as links to the Open Targets Platform target page (by
Ensembl gene id) in a PDF report table cell (raw Typst markup; the
column must not be escaped by format_tt()). Symbols without a valid
Ensembl gene id are shown as plain (escaped) text.

## Usage

``` r
pdf_gene_link(symbol, ensembl_gene_id = NULL)
```

## Arguments

- symbol:

  gene symbol(s)

- ensembl_gene_id:

  Ensembl gene id(s) (e.g. ENSG00000133703)

## Value

character vector with Typst markup
