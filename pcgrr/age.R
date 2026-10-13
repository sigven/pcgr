.libPaths(c(Sys.getenv("SCRATCH_LIB"), .libPaths()))
suppressPackageStartupMessages(library(pcgrr))
options("PCGRR_LOG4R_LOGGER" = log4r::logger(threshold = "WARN", appenders = log4r::console_appender()))
rep <- readRDS("BLCA_report.rds"); rep$settings$conf$debug <- FALSE
rep$settings$conf$sample_properties$age <- 37L
pcgrr::write_report_pdf(rep)
