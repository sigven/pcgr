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
sample_name <- as.character(args[3])
configuration_file <- as.character(args[4])
version <- as.character(args[5])
genome_assembly <- as.character(args[6])
virtual_panel_id <- as.integer(args[7])
diagnostic_grade_only <- as.integer(args[8])
data_dir <- as.character(args[9])

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

cpsr_config <- NULL
default_configuration_file <- paste0(data_dir,'/data/',genome_assembly,'/cpsr_configuration_default.toml')
if(file.exists(default_configuration_file)){
	cpsr_config <- RcppTOML::parseTOML(default_configuration_file, fromFile = T)
}
user_config <- RcppTOML::parseTOML(configuration_file, fromFile = T)

## overwrite default config
for(section in names(cpsr_config)){
  if(!is.null(user_config[[section]])){
    for(element in names(cpsr_config[[section]])){
      if(!is.null(user_config[[section]][[element]])){
        cpsr_config[[section]][[element]] <- user_config[[section]][[element]]
      }
    }
  }
}

cps_report <- pcgrr::generate_predisposition_report(dir, query_vcf2tsv, pcgr_data, cpsr_config, virtual_panel_id, diagnostic_grade_only, sample_name)
pcgrr::write_report(dir, cps_report, sample_name, genome_assembly, tier_model = "cpsr", format = 'html')
pcgrr::write_report(dir, cps_report, sample_name, genome_assembly, tier_model = "cpsr", format = 'json')

