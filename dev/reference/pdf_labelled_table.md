# Table with a label above it (pdf_label()) for a PDF report section, as raw Typst. A sticky label only stays with the start of a breakable table - when the first row(s) do not fit at the bottom of a page, the table moves on and the label is left behind. Tables that fit on a page (at most max_lines lines of text) are therefore kept together with their label in one unbreakable block; longer tables may break across pages

Table with a label above it (pdf_label()) for a PDF report section, as
raw Typst. A sticky label only stays with the start of a breakable
table - when the first row(s) do not fit at the bottom of a page, the
table moves on and the label is left behind. Tables that fit on a page
(at most max_lines lines of text) are therefore kept together with their
label in one unbreakable block; longer tables may break across pages

## Usage

``` r
pdf_labelled_table(label, tt, n_lines, max_lines = 35)
```

## Arguments

- label:

  label text, in bold (plain text)

- tt:

  tinytable object

- n_lines:

  number of text lines of the table body (e.g. rows, plus extra lines of
  multi-line cells)

- max_lines:

  maximum number of lines for a table kept together

## Value

character, to be printed with cat() in a chunk with output 'asis'
