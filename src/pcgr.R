#!/usr/bin/env Rscript

options(warn=-1)
.libPaths(R.home("library"))

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
assay <- as.character(args[12])
t_only <- as.integer(args[13])
estimate_tmb <- as.integer(args[14])
tmb_algorithm <- as.character(args[15])
estimate_msi <- as.integer(args[16])
estimate_signatures <- as.integer(args[17])
target_size_mb <- as.numeric(args[18])
logr_homdel <- as.numeric(args[19])
logr_gain <- as.numeric(args[20])
cna_overlap_pct <- as.numeric(args[21])
min_snv_signatures <- as.integer(args[22])
all_reference_signatures <- as.integer(args[23])
tumor_af_min <- as.numeric(args[24])
tumor_dp_min <- as.numeric(args[25])
control_af_max <- as.numeric(args[26])
control_dp_min <- as.numeric(args[27])
cell_line <- as.integer(args[28])
include_trials <- as.integer(args[29])
tumor_site <- as.character(args[30])


tumor_site <- stringr::str_replace_all(tumor_site,"_"," ")
tumor_site <- stringr::str_replace_all(tumor_site,"@","/")

rlogging::SetTimeStampFormat(ts.format="%Y-%m-%d %H:%M:%S ")
rlogging::SetLogFile(NULL)

rlogging::message(paste0("Tumor primary site: ",tumor_site))

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

## Read default PCGR configurations
pcgr_config <- NULL
default_configuration_file <- paste0(data_dir,'/data/',genome_assembly,'/pcgr_configuration_default.toml')
if(file.exists(default_configuration_file)){
	pcgr_config <- RcppTOML::parseTOML(default_configuration_file, fromFile = T)
}
user_config <- RcppTOML::parseTOML(configuration_file, fromFile = T)

## Override with user-defined settings
for(section in names(pcgr_config)){
  if(!is.null(user_config[[section]])){
    for(element in names(pcgr_config[[section]])){
      if(!is.null(user_config[[section]][[element]])){
        pcgr_config[[section]][[element]] <- user_config[[section]][[element]]
      }
    }
  }
}

if(query_cnv == "None"){
  query_cnv <- NULL
}

## Append settings to PCGR configuration object (fed as arguments to main Python script (pcgr.py))

## Tumor properties (t_props): purity, ploidy, tumor type
pcgr_config[['t_props']] <- list()
pcgr_config[['t_props']][['tumor_purity']] <- purity_estimate
pcgr_config[['t_props']][['tumor_ploidy']] <- ploidy_estimate
pcgr_config[['t_props']][['tumor_type']] <- tumor_site
if(tumor_site == "Any"){
  pcgr_config[['t_props']][['tumor_type']] <- "Cancer, NOS"
}

## Sequencing assay properties (VCF)
## Target (WES/WGS/TARGETED), mode (tumor-control/tumor-only), coding target size
pcgr_config[['assay_props']] <- list()
pcgr_config[['assay_props']][['mode']] <- 'Tumor-Control'
pcgr_config[['assay_props']][['vcf_tumor_only']] <- FALSE
pcgr_config[['assay_props']][['target_size_mb']] <- target_size_mb
pcgr_config[['assay_props']][['type']] <- assay

if(t_only == 1){
  pcgr_config[['assay_props']][['vcf_tumor_only']] <- TRUE
  pcgr_config[['assay_props']][['mode']] <- 'Tumor-Only'
  if(cell_line == 1){
      pcgr_config[['assay_props']][['mode']] <- 'Cell line (Tumor-Only)'
  }
}

## Clinical trials 
pcgr_config[['clinicaltrials']] <- list()
pcgr_config[['clinicaltrials']][['run']] <- as.logical(include_trials)

if(pcgr_config[['t_props']][['tumor_type']] == "Cancer, NOS"){
  rlogging::message(paste0("Clinical trials will not be included in the report when primary site is not specified - skipping"))
  pcgr_config[['clinicaltrials']][['run']] <- F

}

## Analyses to be performed and included in report (with options for mutational signatures)
pcgr_config[['tmb']] <- list()
pcgr_config[['tmb']][['run']] <- as.logical(estimate_tmb)
pcgr_config[['tmb']][['algorithm']] <- as.character(tmb_algorithm)
pcgr_config[['msi']] <- list()
pcgr_config[['msi']][['run']] <- as.logical(estimate_msi)
pcgr_config[['msigs']] <- list()
pcgr_config[['msigs']][['run']] <- as.logical(estimate_signatures)
pcgr_config[['msigs']][['all_reference_signatures']] <- as.logical(all_reference_signatures)
pcgr_config[['msigs']][['mutation_limit']] <- min_snv_signatures

## Copy number aberration (CNA) settings
pcgr_config[['cna']] <- list()
pcgr_config[['cna']][['log_r_homdel']] <- logr_homdel
pcgr_config[['cna']][['log_r_gain']] <- logr_gain
pcgr_config[['cna']][['cna_overlap_pct']] <- cna_overlap_pct

## Allelic support settings (VCF)
pcgr_config[['allelic_support']][['tumor_af_min']] <- tumor_af_min
pcgr_config[['allelic_support']][['tumor_dp_min']] <- tumor_dp_min
pcgr_config[['allelic_support']][['control_dp_min']] <- control_dp_min
pcgr_config[['allelic_support']][['control_af_max']] <- control_af_max

## Generate report object
pcg_report <- pcgrr::generate_pcgr_report(dir, query_vcf2tsv,pcgr_data, 
                                     pcgr_version = version, pcgr_config, sample_name = sample_name, 
                                     cna_segments_tsv = query_cnv, cna_plot = query_cna_plot, 
                                     tier_model = 'pcgr_acmg')

## Write report and result files
if(!is.null(pcg_report)){
  pcgrr::write_report_output(dir, pcg_report, sample_name, genome_assembly, output_format = 'json')
  pcgrr::write_report_output(dir, pcg_report, sample_name, genome_assembly, output_format = 'snv_tsv')
  pcgrr::write_report_output(dir, pcg_report, sample_name, genome_assembly, output_format = 'msigs_tsv')
  if(pcgr_config[['assay_props']][['vcf_tumor_only']] == T){
    pcgrr::write_report_output(dir, pcg_report, sample_name, genome_assembly, output_format = 'snv_tsv_unfiltered')
  }
  pcgrr::write_report_output(dir, pcg_report, sample_name, genome_assembly, output_format = 'cna_tsv')
  pcgrr::write_report_output(dir, pcg_report, sample_name, genome_assembly, output_format = 'html')
  pcgrr::write_report_output(dir, pcg_report, sample_name, genome_assembly, output_format = 'html', flexdb = T)

}

