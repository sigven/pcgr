#!/usr/bin/env Rscript

options(warn=-1)
.libPaths(R.home("library")) # use conda R pkgs, not e.g. user's local installation

suppressWarnings(suppressPackageStartupMessages(library(pcgrr)))
suppressWarnings(suppressPackageStartupMessages(library(log4r)))
suppressWarnings(suppressPackageStartupMessages(library(argparse)))

args <- commandArgs(trailingOnly=TRUE)

## YAML file produced by PCGR Python workflow
## - settings and paths to reference data and annotated input sample files
yaml_fname <- as.character(args[1])
quarto_evars_path <- as.character(args[2])
pcgrr::export_quarto_evars(quarto_evars_path)

my_log4r_layout <- function(level, ...) {
  paste0(format(Sys.time()), " - pcgr-report-generation - ",
         level, " - ", ..., "\n", collapse = "")
}

log4r_logger <-
  log4r::logger(
    threshold = "INFO", appenders = log4r::console_appender(my_log4r_layout))

## this gets passed on to all the log4r_* functions inside the pkg
options("PCGRR_LOG4R_LOGGER" = log4r_logger)

## Generate report content
pcg_report <- pcgrr::generate_report(
  yaml_fname = yaml_fname
)

## Write report contents to output files (HTML, XLSX, TSV)
if (!is.null(pcg_report)) {
  if(pcg_report$settings$conf$other$no_html == FALSE){
    pcgrr::write_report_quarto_html(report = pcg_report)
  }
  else{
    pcgrr::log4r_info("Skipping HTML report generation (option '--no_html' set to TRUE)")
  }

  pcgrr::write_report_excel(report = pcg_report)
  pcgrr::write_report_tsv(report = pcg_report, output_type = 'snv_indel')
  pcgrr::write_report_tsv(report = pcg_report, output_type = 'snv_indel_unfiltered')
  pcgrr::write_report_tsv(report = pcg_report, output_type = 'cna_gene')
  pcgrr::write_report_tsv(report = pcg_report, output_type = 'msigs')
}
