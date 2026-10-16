# Note box in a PDF report section (raw Typst block)

Note box in a PDF report section (raw Typst block)

## Usage

``` r
pdf_note(
  text,
  title = NULL,
  kind = c("note", "warning", "empty", "placeholder")
)
```

## Arguments

- text:

  note text - markdown for 'note' and 'warning' (rendered by Quarto),
  plain text for 'empty' and 'placeholder' (Typst markup characters are
  escaped)

- title:

  short title in bold (default: 'Note', or 'No findings' for 'empty';
  not used for 'placeholder')

- kind:

  'note' (callout, main report color), 'warning' (callout, amber),
  'empty' (grey box for a section without findings), or 'placeholder'
  (content not yet implemented). 'note' and 'placeholder' use Typst
  definitions of the PCGR report preamble (pcgr-main, pcgr-stripe,
  pcgr-placeholder()); 'warning' and 'empty' are self-contained

## Value

character, to be printed with cat() in a chunk with output 'asis'
