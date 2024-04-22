#!/usr/bin/env Rscript

options(warn=-1)
.libPaths(R.home("library")) # use conda R pkgs, not e.g. user's local installation

args <- commandArgs(trailingOnly=TRUE)

suppressWarnings(suppressPackageStartupMessages(library(cpsr)))
suppressWarnings(suppressPackageStartupMessages(library(log4r)))

## YAML file produced by CPSR Python workflow
## - settings and paths to reference data and annotated input sample files
yaml_fname <- as.character(args[1])

my_log4r_layout <- function(level, ...) {
  paste0(format(Sys.time()), " - cpsr-report-generation - ",
         level, " - ", ..., "\n", collapse = "")
}

log4r_logger <-
  log4r::logger(
    threshold = "INFO", appenders = log4r::console_appender(my_log4r_layout))

# this gets passed on to all the log4r_* functions inside the pkg
options("PCGRR_LOG4R_LOGGER" = log4r_logger)

## Generate report content
cps_report <- cpsr::generate_cpsr_report(
  yaml_fname = yaml_fname
)

## Write report contents to output files (HTML, XLSX, TSV)
if(!is.null(cps_report)){
  cpsr::write_cpsr_output(
      cps_report,
      output_format = 'tsv')
  cpsr::write_cpsr_output(
    cps_report,
    output_format = 'xlsx')
  cpsr::write_cpsr_output(
     cps_report,
     output_format = 'html')
}
