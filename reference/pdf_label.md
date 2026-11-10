# Label above a table in a PDF report section (e.g. 'Therapeutic sensitivity') - a sticky block (pcgr-label(), defined in typst/\_preamble.typ), kept on the same page as the table that follows rather than left alone at the bottom of a page

Label above a table in a PDF report section (e.g. 'Therapeutic
sensitivity') - a sticky block (pcgr-label(), defined in
typst/\_preamble.typ), kept on the same page as the table that follows
rather than left alone at the bottom of a page

## Usage

``` r
pdf_label(label, suffix = NULL)
```

## Arguments

- label:

  label text, in bold (plain text - Typst markup characters are escaped)

- suffix:

  optional text following the label, not in bold (plain text)

## Value

character, to be printed with cat() in a chunk with output 'asis'
