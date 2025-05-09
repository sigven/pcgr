% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/utils.R
\name{exclude_non_chrom_variants}
\alias{exclude_non_chrom_variants}
\title{Function that excludes genomic aberrations from non-nuclear chromosomes}
\usage{
exclude_non_chrom_variants(vcf_df, chrom_var = "CHROM")
}
\arguments{
\item{vcf_df}{data frame}

\item{chrom_var}{variable name of chromosome in data frame}
}
\value{
vcf_df data frame with mutations from nuclear chromosomes only
}
\description{
Function that excludes genomic aberrations from non-nuclear chromosomes
}
