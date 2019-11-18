#!/usr/bin/env Rscript

options(warn=-1)
args = commandArgs(trailingOnly=TRUE)

suppressWarnings(suppressPackageStartupMessages(library(plyr)))
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
query_cna_plot <- as.character(args[9])
purity_estimate <- as.character(args[10])
ploidy_estimate <- as.character(args[11])
t_only <- as.integer(args[12])
target_size_mb <- as.numeric(args[13])
tumor_type <- as.character(args[14])

rlogging::SetTimeStampFormat(ts.format="%Y-%m-%d %H:%M:%S ")
rlogging::SetLogFile(NULL)

pcgr_data <- readRDS(paste0(data_dir,'/data/',genome_assembly,'/rds/pcgr_data.rds'))

pcgr_data[['assembly']][['seqinfo']] <- 
   GenomeInfoDb::Seqinfo(seqnames = GenomeInfoDb::seqlevels(GenomeInfoDb::seqinfo(BSgenome.Hsapiens.UCSC.hg38)), seqlengths = GenomeInfoDb::seqlengths(GenomeInfoDb::seqinfo(BSgenome.Hsapiens.UCSC.hg38)), genome = 'hg38')
pcgr_data[['assembly']][['bsg']] <- BSgenome.Hsapiens.UCSC.hg38
if(genome_assembly == 'grch37'){
  pcgr_data[['assembly']][['bsg']] <- BSgenome.Hsapiens.UCSC.hg19
  pcgr_data[['assembly']][['seqinfo']] <- 
     GenomeInfoDb::Seqinfo(seqnames = GenomeInfoDb::seqlevels(GenomeInfoDb::seqinfo(BSgenome.Hsapiens.UCSC.hg19)), seqlengths = GenomeInfoDb::seqlengths(GenomeInfoDb::seqinfo(BSgenome.Hsapiens.UCSC.hg19)), genome = 'hg19')
}

pcgr_config <- NULL
default_configuration_file <- paste0(data_dir,'/data/',genome_assembly,'/pcgr_configuration_default.toml')
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

## append optional estimates of tumor purity and ploidy provided by user
pcgr_config[['tumor_properties']] <- list()
pcgr_config[['tumor_properties']][['tumor_purity']] <- purity_estimate
pcgr_config[['tumor_properties']][['tumor_ploidy']] <- ploidy_estimate
pcgr_config[['mutational_burden']][['target_size_mb']] <- target_size_mb
pcgr_config[['tumor_type']] <- list()
pcgr_config[['tumor_type']][['type']] <- tumor_type
if(tumor_type == "Cancer_NOS"){
  pcgr_config[['tumor_type']][['type']] <- ""
}
pcg_report <- pcgrr::generate_report(dir, query_vcf2tsv, pcgr_data, pcgr_version = version, pcgr_config, sample_name = sample_name, cna_segments_tsv = query_cnv, 
				     cna_plot = query_cna_plot, tier_model = 'pcgr_acmg', tumor_only = t_only)
pcgrr::write_report(dir, pcg_report, sample_name, genome_assembly, tier_model = 'pcgr_acmg', format = 'html')
pcgrr::write_report(dir, pcg_report, sample_name, genome_assembly, tier_model = 'pcgr_acmg', format = 'json')

