#!/usr/bin/env Rscript

options(warn=-1)
args = commandArgs(trailingOnly=TRUE)

suppressWarnings(suppressPackageStartupMessages(library(pcgrr)))
suppressWarnings(suppressPackageStartupMessages(library(magrittr)))
suppressWarnings(suppressPackageStartupMessages(library(BSgenome.Hsapiens.UCSC.hg19)))
suppressWarnings(suppressPackageStartupMessages(library(deconstructSigs)))
suppressWarnings(suppressPackageStartupMessages(library(randomForest)))
suppressWarnings(suppressPackageStartupMessages(library(caret)))


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

pcgrr::generate_report(project_directory = dir, query_vcf = query_vcf, cna_segments_tsv = query_cnv, pcgr_data = pcgr_data, sample_name = sample_name, configuration_file = configuration_file, pcgr_version = version)
