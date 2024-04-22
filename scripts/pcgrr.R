#!/usr/bin/env Rscript

options(warn=-1)
.libPaths(R.home("library")) # use conda R pkgs, not e.g. user's local installation

suppressWarnings(suppressPackageStartupMessages(library(pcgrr)))
suppressWarnings(suppressPackageStartupMessages(library(log4r)))
suppressWarnings(suppressPackageStartupMessages(library(argparse)))
suppressWarnings(suppressPackageStartupMessages(library(CNAqc)))

args <- commandArgs(trailingOnly=TRUE)

yaml_fname <- as.character(args[1])

my_log4r_layout <- function(level, ...) {
  paste0(format(Sys.time()), " - pcgr-report-generation - ",
         level, " - ", ..., "\n", collapse = "")
}

log4r_logger <-
  log4r::logger(
    threshold = "INFO", appenders = log4r::console_appender(my_log4r_layout))

# this gets passed on to all the log4r_* functions inside the pkg
options("PCGRR_LOG4R_LOGGER" = log4r_logger)

## Generate report content
pcg_report <- pcgrr::generate_report(
  yaml_fname = yaml_fname
)

# Write result files (HTML, xlsx, TSV)
if (!is.null(pcg_report)) {
  pcgrr::write_report_quarto_html(report = pcg_report)
  pcgrr::write_report_excel(report = pcg_report)
  pcgrr::write_report_tsv(report = pcg_report, variant_type = 'snv_indel')
  pcgrr::write_report_tsv(report = pcg_report, variant_type = 'cna_gene')
}
