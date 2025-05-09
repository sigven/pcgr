% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/mutational_signatures.R
\name{plot_signature_contributions}
\alias{plot_signature_contributions}
\title{Function that makes plots of mutational signature contributions
in a given sample (both ggplot and plotly)}
\usage{
plot_signature_contributions(
  signature_contributions = NULL,
  per_signature = TRUE
)
}
\arguments{
\item{signature_contributions}{A list with two data frames: 'per_group' and
'per_signature'.}

\item{per_signature}{Logical. If TRUE, the plot will show the contribution
per signature. If FALSE, the plot will show the contribution per group.}
}
\description{
Function that makes plots of mutational signature contributions
in a given sample (both ggplot and plotly)
}
