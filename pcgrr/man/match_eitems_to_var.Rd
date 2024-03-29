% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/biomarkers.R
\name{match_eitems_to_var}
\alias{match_eitems_to_var}
\title{Function that matches variants to evidence items}
\usage{
match_eitems_to_var(
  sample_calls,
  db = "civic",
  colset = NULL,
  eitems = NULL,
  region_marker = T
)
}
\arguments{
\item{sample_calls}{data frame with variant calls}

\item{db}{database with evidence items ('civic','cgi')}

\item{colset}{character vector with column names to pull out from sample_calls}

\item{eitems}{raw list of evidence items}

\item{region_marker}{logical indication if region biomarkers are to be matched or not}
}
\description{
Function that matches variants to evidence items
}
