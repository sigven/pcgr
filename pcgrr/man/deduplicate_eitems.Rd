% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/biomarkers.R
\name{deduplicate_eitems}
\alias{deduplicate_eitems}
\title{Function that removes redundancy in variant evidence items (i.e. if
a variant is assicated with evidence at the codon level, evidence
at the exon/gene level is ignored)}
\usage{
deduplicate_eitems(
  var_eitems = NULL,
  target_type = "exact",
  target_other = c("codon", "exon", "gene")
)
}
\arguments{
\item{var_eitems}{data frame with variant evidence items}

\item{target_type}{which resolution level should be used as the
"best" level ('exact' or 'codon)}

\item{target_other}{resolution levels for other evidence items
that should be ignored if evidence is found at the target_type level}
}
\description{
Function that removes redundancy in variant evidence items (i.e. if
a variant is assicated with evidence at the codon level, evidence
at the exon/gene level is ignored)
}
