#!/usr/bin/env Rscript

options(warn=-1)
args = commandArgs(trailingOnly=TRUE)

suppressWarnings(suppressPackageStartupMessages(library(pcgrr)))
suppressWarnings(suppressPackageStartupMessages(library(magrittr)))
suppressWarnings(suppressPackageStartupMessages(library(BSgenome.Hsapiens.UCSC.hg19)))
suppressWarnings(suppressPackageStartupMessages(library(deconstructSigs)))
suppressWarnings(suppressPackageStartupMessages(library(randomForest)))
suppressWarnings(suppressPackageStartupMessages(library(caret)))
suppressWarnings(suppressPackageStartupMessages(library(RcppTOML)))



dir <- as.character(args[1])
query_vcf <- as.character(args[2])
query_cnv <- as.character(args[3])
sample_name <- as.character(args[4])
configuration_file <- as.character(args[5])
version <- as.character(args[6])

rlogging::SetTimeStampFormat(ts.format="%Y-%m-%d %H:%M:%S ")
rlogging::SetLogFile(NULL)

load('/data/data/rda/pcgr_data.rda')

eval_tier1 <- FALSE
eval_tier2 <- FALSE
eval_tier3 <- FALSE
eval_tier4 <- FALSE
eval_tier5 <- FALSE
eval_signature_report <- FALSE
eval_missing_signature_data <- FALSE
eval_cna_segments <- FALSE
eval_cna_plot <- FALSE
eval_cna_loss <- FALSE
eval_cna_gain <- FALSE
eval_cna_biomarker <- FALSE
eval_msi_report <- FALSE
eval_missing_msi_data <- FALSE
eval_tumor_only <- FALSE


pcgr_config <- NULL
default_configuration_file <- '/data/data/pcgr_configuration_default.toml'
if(file.exists(default_configuration_file)){
	pcgr_config <- RcppTOML::parseTOML(default_configuration_file, fromFile = T)
}
user_config <- RcppTOML::parseTOML(configuration_file, fromFile = T)

## overwrite default config
for(section in names(pcgr_config)){
  if(!is.null(user_config[[section]])){
    for(element in names(pcgr_config[[section]])){
      if(!is.null(user_config[[section]][[element]])){
        pcgr_config[[section]][[element]] <- user_config[[section]][[element]]
      }
    }
  }
}


pcgrr::generate_report(project_directory = dir, query_vcf = query_vcf, cna_segments_tsv = query_cnv, pcgr_data = pcgr_data, pcgr_config = pcgr_config, sample_name = sample_name, pcgr_version = version)
