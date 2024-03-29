% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/main.R
\name{generate_pcgr_report}
\alias{generate_pcgr_report}
\title{Function that generates all contents of the cancer genome report (PCGR)}
\usage{
generate_pcgr_report(
  project_directory = NULL,
  pcgr_data = NULL,
  config = NULL,
  tier_model = "pcgr_acmg"
)
}
\arguments{
\item{project_directory}{name of project directory}

\item{pcgr_data}{List object with multiple PCGR data bundle annotations}

\item{config}{Object with PCGR configuration parameters}

\item{tier_model}{Variant tier model}
}
\description{
Function that generates all contents of the cancer genome report (PCGR)
}
