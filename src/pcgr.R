#!/usr/bin/env Rscript

options(warn=-1)
args = commandArgs(trailingOnly=TRUE)

suppressWarnings(suppressPackageStartupMessages(library(pcgrr2)))
suppressWarnings(suppressPackageStartupMessages(library(magrittr)))
suppressWarnings(suppressPackageStartupMessages(library(BSgenome.Hsapiens.UCSC.hg19)))
suppressWarnings(suppressPackageStartupMessages(library(deconstructSigs)))
suppressWarnings(suppressPackageStartupMessages(library(randomForest)))
suppressWarnings(suppressPackageStartupMessages(library(caret)))


dir <- as.character(args[1])
query_vcf <- as.character(args[2])
query_cnv <- as.character(args[3])
sample_name <- as.character(args[4])
logR_gain <- as.numeric(args[5])
logR_homdel <- as.numeric(args[6])
pred_MSI <- as.integer(args[7])
ident_msigs <- as.integer(args[8])
signatures_limit <- as.integer(args[9])
signature_normalization <- as.character(args[10])
show_noncoding_variants <- as.integer(args[11])
tumor_dp_tag <- as.character(args[12])
tumor_af_tag <- as.character(args[13])
normal_dp_tag <- as.character(args[14])
normal_af_tag <- as.character(args[15])
call_conf_tag <- as.character(args[16])

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
eval_cna_loss <- FALSE
eval_cna_gain <- FALSE
eval_cna_biomarker <- FALSE
eval_msi_report <- FALSE
eval_missing_msi_data <- FALSE
predict_MSI <- FALSE
show_noncoding <- FALSE
if(pred_MSI == 1){
	predict_MSI <- TRUE
}
identify_msigs <- FALSE
if(ident_msigs == 1){
	identify_msigs <- TRUE
}
if(show_noncoding_variants == 1){
	show_noncoding <- TRUE
}


pcgrr2::generate_pcg_report(project_directory = dir, query_vcf = query_vcf, cna_segments_tsv = query_cnv, logR_gain = logR_gain, logR_homdel = logR_homdel, pcgr_data = pcgr_data, sample_name = sample_name, signature_normalization_method = signature_normalization, signatures_limit = signatures_limit, print_biomarkers = TRUE, print_tier_variants = TRUE, print_mutational_signatures = TRUE, print_cna_segments = TRUE, print_maf = TRUE, print_html_report = TRUE, show_noncoding = show_noncoding, predict_MSI = predict_MSI, identify_msigs = identify_msigs, tumor_dp_tag = tumor_dp_tag, tumor_af_tag = tumor_af_tag, normal_dp_tag = normal_dp_tag, normal_af_tag = normal_af_tag, call_conf_tag = call_conf_tag)
