#!/usr/bin/env Rscript

options(warn=-1)
args = commandArgs(trailingOnly=TRUE)

suppressWarnings(suppressPackageStartupMessages(library(pcgrr)))
suppressWarnings(suppressPackageStartupMessages(library(magrittr)))
suppressWarnings(suppressPackageStartupMessages(library(BSgenome.Hsapiens.UCSC.hg19)))
suppressWarnings(suppressPackageStartupMessages(library(BSgenome.Hsapiens.UCSC.hg38)))
suppressWarnings(suppressPackageStartupMessages(library(deconstructSigs)))
suppressWarnings(suppressPackageStartupMessages(library(randomForest)))
suppressWarnings(suppressPackageStartupMessages(library(caret)))
suppressWarnings(suppressPackageStartupMessages(library(RcppTOML)))

dir <- as.character(args[1])
query_vcf2tsv <- as.character(args[2])
query_cnv <- as.character(args[3])
sample_name <- as.character(args[4])
configuration_file <- as.character(args[5])
version <- as.character(args[6])
genome_assembly <- as.character(args[7])
data_dir <- as.character(args[8])

rlogging::SetTimeStampFormat(ts.format="%Y-%m-%d %H:%M:%S ")
rlogging::SetLogFile(NULL)

load(paste0(data_dir,'/data/',genome_assembly,'/rda/pcgr_data.rda'))

pcgr_config <- NULL
default_configuration_file <- paste0(data_dir,'/data/',genome_assembly,'/pcgr_configuration_somatic_default.toml')
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

if(pcgr_config$tier_model$tier_model == 'pcgr'){
  pcgrr::generate_report(dir, query_vcf2tsv, pcgr_data, pcgr_config, sample_name, query_cnv, version, genome_assembly = genome_assembly)
}else{
  pcgrr::generate_report_acmg(dir, query_vcf2tsv, pcgr_data, pcgr_config, sample_name, query_cnv, version, genome_assembly = genome_assembly)
}
