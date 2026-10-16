# Evaluate an expression (report rendering) with a UTF-8 locale

Quarto/knitr run R in a subprocess that inherits the locale; with a
non-UTF-8 locale (e.g. 'C', as in a conda environment without LANG set),
non-ASCII characters of the report text (e.g. dashes, middle dots) are
written as escape codes ('\<U+00B7\>'). If the current locale is not
UTF-8, LC_ALL is temporarily set to an available UTF-8 locale.

## Usage

``` r
with_utf8_locale(expr)
```

## Arguments

- expr:

  expression to evaluate
