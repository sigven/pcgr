% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/biomarkers.R
\name{load_eitems}
\alias{load_eitems}
\title{Function that loads specific set of clinical variant evidence items (CIViC + CGI) based
on given parameters (mutation type, variant origin, tumor type etc)}
\usage{
load_eitems(
  eitems_raw = NULL,
  ontology = NULL,
  alteration_types = c("MUT"),
  origin = "Somatic",
  tumor_type_specificity = NULL,
  tumor_type = NULL
)
}
\arguments{
\item{eitems_raw}{complete set of clinical variant evidence items}

\item{ontology}{phenotype ontology data frame}

\item{alteration_types}{types of alteration ('MUT', 'CNA', 'MUT_LOF')}

\item{origin}{variant origin ('Somatic','Germline')}

\item{tumor_type_specificity}{tumor type specificity ('any', 'specific')}

\item{tumor_type}{primary tumor site}
}
\value{
eitems variant evidence items
}
\description{
Function that loads specific set of clinical variant evidence items (CIViC + CGI) based
on given parameters (mutation type, variant origin, tumor type etc)
}
