# Reference sources of the expression similarity analysis, as display labels (run configuration / report settings)

Reference sources of the expression similarity analysis, as display
labels (run configuration / report settings)

## Usage

``` r
expression_similarity_sources(conf)
```

## Arguments

- conf:

  PCGR run configuration (settings - conf); the sources are the names of
  conf - expression - similarity_db ('–expression_sim_db')

## Value

character vector of source labels (in the order TCGA, DepMap,
Treehouse), empty if the similarity analysis is not run
