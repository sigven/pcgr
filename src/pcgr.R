#!/usr/bin/env Rscript

options(warn=-1)
args = commandArgs(trailingOnly=TRUE)

suppressWarnings(suppressPackageStartupMessages(library(pcgrr2)))
suppressWarnings(suppressPackageStartupMessages(library(magrittr)))
suppressWarnings(suppressPackageStartupMessages(library(BSgenome.Hsapiens.UCSC.hg19)))
suppressWarnings(suppressPackageStartupMessages(library(deconstructSigs)))

dir <- as.character(args[1])
query_vcf <- as.character(args[2])
query_cnv <- as.character(args[3])
sample_name <- as.character(args[4])
logR_threshold_amplification <- as.numeric(args[5])
logR_threshold_homozygous_deletion <- as.numeric(args[6])

rlogging::SetTimeStampFormat(ts.format="%Y-%m-%d %H:%M:%S ")
rlogging::SetLogFile(NULL)

load('/data/data/rda/pcgr_data.rda')

pcgrr2::generate_pcg_report(project_directory = dir, query_vcf = query_vcf, cna_segments_tsv = query_cnv, sample_name = sample_name, logR_threshold_amplification = logR_threshold_amplification, logR_threshold_homozygous_deletion = logR_threshold_homozygous_deletion, print_biomarkers = TRUE, print_tier_variants = TRUE, print_mutational_signatures = TRUE, print_cna_segments = TRUE, print_maf = TRUE, print_html_report = TRUE)
