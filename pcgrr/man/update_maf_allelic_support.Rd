% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/utils.R
\name{update_maf_allelic_support}
\alias{update_maf_allelic_support}
\title{Function that updates a MAF file (produced by vcf2maf) with allelic support data}
\usage{
update_maf_allelic_support(calls, maf_fname_tmp, maf_fname, delete_raw = T)
}
\arguments{
\item{calls}{data frame with variant calls}

\item{maf_fname_tmp}{Filename for MAF produced by vcf2maf}

\item{maf_fname}{Filename for MAF that includes variant allelic support (final MAF)}

\item{delete_raw}{logical indicating if initial MAF should be deleted}
}
\description{
Function that updates a MAF file (produced by vcf2maf) with allelic support data
}
