% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/acmg.R
\name{assign_acmg_tiers}
\alias{assign_acmg_tiers}
\title{Function that assigns tier classifications to somatic CNA segments and
SNVs/InDels, based on the presence of biomarker evidence found in
the variant set}
\usage{
assign_acmg_tiers(
  vartype = "snv_indel",
  primary_site = "Any",
  variants_df = NULL,
  biomarker_items = NULL
)
}
\arguments{
\item{vartype}{variant type ('snv_indel' or 'cna')}

\item{primary_site}{primary tumor site}

\item{variants_df}{data frame with variants (SNVs/InDels or CNAs)}

\item{biomarker_items}{data frame with biomarker evidence items}
}
\description{
Function that assigns tier classifications to somatic CNA segments and
SNVs/InDels, based on the presence of biomarker evidence found in
the variant set
}
