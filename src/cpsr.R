#!/usr/bin/env Rscript

options(warn=-1)
.libPaths(R.home("library"))

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
pcgr_version <- as.character(args[5])
cpsr_version <- as.character(args[6])
genome_assembly <- as.character(args[7])
virtual_panel_id <- as.integer(args[8])
custom_bed <- as.character(args[9])
maf_upper_threshold <- as.numeric(args[10])
diagnostic_grade_only <- as.integer(args[11])
data_dir <- as.character(args[12])
ignore_noncoding <- as.integer(args[13])
clinvar_ignore_noncancer <- as.integer(args[14])
secondary_findings <- as.integer(args[15])
gwas_findings <- as.integer(args[16])
classify_all <- as.integer(args[17])

rlogging::SetTimeStampFormat(ts.format="%Y-%m-%d %H:%M:%S ")
rlogging::SetLogFile(NULL)

pcgr_data <- readRDS(paste0(data_dir,'/data/',genome_assembly,'/rds/pcgr_data.rds'))

pcgr_data[['assembly']][['seqinfo']] <- 
   GenomeInfoDb::Seqinfo(seqnames = GenomeInfoDb::seqlevels(GenomeInfoDb::seqinfo(BSgenome.Hsapiens.UCSC.hg38)), 
                         seqlengths = GenomeInfoDb::seqlengths(GenomeInfoDb::seqinfo(BSgenome.Hsapiens.UCSC.hg38)), genome = 'hg38')
pcgr_data[['assembly']][['bsg']] <- BSgenome.Hsapiens.UCSC.hg38
if(genome_assembly == 'grch37'){
  pcgr_data[['assembly']][['bsg']] <- BSgenome.Hsapiens.UCSC.hg19
  pcgr_data[['assembly']][['seqinfo']] <- 
     GenomeInfoDb::Seqinfo(seqnames = GenomeInfoDb::seqlevels(GenomeInfoDb::seqinfo(BSgenome.Hsapiens.UCSC.hg19)), 
                           seqlengths = GenomeInfoDb::seqlengths(GenomeInfoDb::seqinfo(BSgenome.Hsapiens.UCSC.hg19)), genome = 'hg19')
}

## Load default CPSR configurations
cpsr_config <- NULL
default_configuration_file <- paste0(data_dir,'/data/',genome_assembly,'/cpsr_configuration_default.toml')
if(file.exists(default_configuration_file)){
	cpsr_config <- RcppTOML::parseTOML(default_configuration_file, fromFile = T)
}
user_config <- RcppTOML::parseTOML(configuration_file, fromFile = T)

## Override configurations with user-defined settings
for(section in names(cpsr_config)){
  if(!is.null(user_config[[section]])){
    for(element in names(cpsr_config[[section]])){
      if(!is.null(user_config[[section]][[element]])){
        cpsr_config[[section]][[element]] <- user_config[[section]][[element]]
      }
    }
  }
}

cpsr_config[['ignore_noncoding']] <- as.logical(ignore_noncoding)
cpsr_config[['secondary_findings']] <- as.logical(secondary_findings)
cpsr_config[['gwas_findings']] <- as.logical(gwas_findings)
cpsr_config[['classify_all']] <- as.logical(classify_all)
cpsr_config[['maf_upper_threshold']] <- maf_upper_threshold
cpsr_config[['diagnostic_grade_only']] <- as.logical(diagnostic_grade_only)
cpsr_config[['clinvar_ignore_noncancer']] <- as.logical(clinvar_ignore_noncancer)

## Generate report content
cps_report <- pcgrr::generate_cpsr_report(dir, query_vcf2tsv, custom_bed, pcgr_data, pcgr_version, 
                                          cpsr_version, cpsr_config, virtual_panel_id, sample_name)

## Write report contents to output files
if(!is.null(cps_report)){
  pcgrr::write_report_output(dir, cps_report, sample_name, genome_assembly, tier_model = 'cpsr', output_format = 'snv_tsv')
  pcgrr::write_report_output(dir, cps_report, sample_name, genome_assembly, tier_model = "cpsr", output_format = 'json')
  #if(cps_report[['content']][['snv_indel']][['max_dt_rows']] < 5000){
  pcgrr::write_report_output(dir, cps_report, sample_name, genome_assembly, tier_model = "cpsr", output_format = 'html')
  #}
  #else{
  #  rlogging::message("Dataset is too big for client-side DataTables (i.e. >= 5000 pr. tier) - stand-alone HTML report file will not be produced")
  #}

}

