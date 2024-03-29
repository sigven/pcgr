% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/main.R
\name{generate_report_data_snv_indel}
\alias{generate_report_data_snv_indel}
\title{Function that generates tiered variant sets for SNVs/InDels}
\usage{
generate_report_data_snv_indel(
  sample_calls,
  pcgr_data,
  sample_name,
  config,
  callset = "somatic calls",
  biomarker_mapping_stringency = 1,
  tier_model = "pcgr_acmg"
)
}
\arguments{
\item{sample_calls}{variant calls subject to mutational signature analysis}

\item{pcgr_data}{object with PCGR annotation data}

\item{sample_name}{sample identifier}

\item{config}{Object with PCGR configuration parameters}

\item{callset}{type of calls}

\item{biomarker_mapping_stringency}{quality level for biomarkers}

\item{tier_model}{tier model (pcgr_acmg)}
}
\value{
pcg_report_data data frame with all report elements
}
\description{
Function that generates tiered variant sets for SNVs/InDels
}
