# Strip HTML tags

Remove HTML tags and comments from text. From
https://github.com/yihui/xfun/blob/ccee26/R/string.R#L329.

## Usage

``` r
strip_html(x)
```

## Arguments

- x:

  A character vector.

## Value

A character vector with HTML tags and comments stripped off.

## Examples

``` r
strip_html('<a href="#">Hello <!-- comment -->world!</a>')
#> [1] "Hello world!"
```
