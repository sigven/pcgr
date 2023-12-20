#!/usr/bin/env Rscript

.libPaths(R.home("library")) # use conda R pkgs, not e.g. user's local installation

suppressWarnings(suppressPackageStartupMessages(library(argparse)))
suppressWarnings(suppressPackageStartupMessages(library(pcgrr)))
suppressWarnings(suppressPackageStartupMessages(library(stringr)))


# my_log4r_layout <- function(level, ...) {
#   paste0(format(Sys.time()), " - pcgr-report-generation - ",
#          level, " - ", ..., "\n", collapse = "")
# }

# log4r_logger <- log4r::logger(threshold = "INFO",
#                               appenders = log4r::console_appender(my_log4r_layout))

# # this gets passed on to all the log4r_* functions inside the pkg
# options("PCGRR_LOG4R_LOGGER" = log4r_logger)


# pcg_report <- NULL

# defaultW <- getOption("warn")
# options(warn = -1)

# # ## Generate report object
# pcg_report <-
#   pcgrr::generate_pcgr_report(
#     project_directory  = pcgr_config[['required_args']][['output_dir']],
#     pcgr_data = pcgr_data,
#     config = pcgr_config,
#     tier_model = 'pcgr_acmg')

# options(warn = defaultW)


# # ## Write report and result files
# if (!is.null(pcg_report)) {

#   pcgrr::write_report_output(
#     pcg_report,
#     pcgr_config,
#     output_format = 'snv_tsv')
#   pcgrr::write_report_output(
#     pcg_report,
#     pcgr_config,
#     output_format = 'msigs_tsv')
#   if (pcgr_config[['assay_props']][['vcf_tumor_only']] == T){
#     pcgrr::write_report_output(
#       pcg_report,
#       pcgr_config,
#       output_format = 'snv_tsv_unfiltered')
#   }
#   pcgrr::write_report_output(
#     pcg_report,
#     pcgr_config,
#     output_format = 'cna_tsv')
#   pcgrr::write_report_output(
#     pcg_report,
#     pcgr_config,
#     output_format = 'html')
#   pcgrr::write_report_output(
#     pcg_report,
#     pcgr_config,
#     output_format = 'html',
#     flexdb = T)
#   pcgrr::write_report_output(
#     pcg_report,
#     pcgr_config,
#     output_format = 'json')
# }
