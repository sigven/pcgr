#!/usr/bin/env Rscript

options(warn=-1)
.libPaths(R.home("library"))

args = commandArgs(trailingOnly=TRUE)

suppressWarnings(suppressPackageStartupMessages(library(plyr)))
suppressWarnings(suppressPackageStartupMessages(library(pcgrr)))
suppressWarnings(suppressPackageStartupMessages(library(magrittr)))
suppressWarnings(suppressPackageStartupMessages(library(BSgenome.Hsapiens.UCSC.hg19)))
suppressWarnings(suppressPackageStartupMessages(library(BSgenome.Hsapiens.UCSC.hg38)))
suppressWarnings(suppressPackageStartupMessages(library(randomForest)))
suppressWarnings(suppressPackageStartupMessages(library(caret)))

dir <- as.character(args[1])

## Append settings to PCGR configuration object processed in R package
pcgr_config <- list()

pcgr_config[['required_args']] <- list()
pcgr_config[['required_args']][['query_vcf2tsv']] <- as.character(args[2])
pcgr_config[['required_args']][['query_cna']] <- as.character(args[3])
pcgr_config[['required_args']][['query_rna_fusion']] <- as.character(args[4])
pcgr_config[['required_args']][['query_rna_expression']] <- as.character(args[5])
pcgr_config[['required_args']][['cpsr_report']] <- as.character(args[6])
pcgr_config[['required_args']][['sample_name']] <- as.character(args[7])
pcgr_config[['required_args']][['pcgr_version']] <- as.character(args[8])
pcgr_config[['required_args']][['genome_assembly']] <- as.character(args[9])
pcgr_config[['required_args']][['output_dir']] <- dir
pcgr_config[['required_args']][['data_dir']] <- as.character(args[10])

if(pcgr_config[['required_args']][['query_cna']] == "None"){
  pcgr_config[['required_args']][['query_cna']] <- NULL
}
if(pcgr_config[['required_args']][['query_rna_fusion']] == "None"){
  pcgr_config[['required_args']][['query_rna_fusion']] <- NULL
}
if(pcgr_config[['required_args']][['query_rna_expression']] == "None"){
  pcgr_config[['required_args']][['query_rna_expression']] <- NULL
}
if(pcgr_config[['required_args']][['cpsr_report']] == "None"){
  pcgr_config[['required_args']][['cpsr_report']] <- NULL
}


## Tumor properties (t_props): purity, ploidy, tumor type
pcgr_config[['t_props']] <- list()
pcgr_config[['t_props']][['tumor_purity']] <- as.numeric(args[11])
pcgr_config[['t_props']][['tumor_ploidy']] <- as.numeric(args[12])
pcgr_config[['t_props']][['tumor_type']] <- stringr::str_replace_all(
  stringr::str_replace_all(
    as.character(args[13]),"_"," "),
  "@","/")
if(pcgr_config[['t_props']][['tumor_type']] == "Any"){
  pcgr_config[['t_props']][['tumor_type']] <- "Cancer, NOS"
}

## Sequencing assay properties (VCF)
## Target (WES/WGS/TARGETED), mode (tumor-control/tumor-only), coding target size
pcgr_config[['assay_props']] <- list()
pcgr_config[['assay_props']][['mode']] <- 'Tumor-Control'
pcgr_config[['assay_props']][['vcf_tumor_only']] <- FALSE
pcgr_config[['assay_props']][['target_size_mb']] <- as.numeric(args[14])
pcgr_config[['assay_props']][['type']] <- as.character(args[15])

if(as.integer(args[16]) == 1){
  pcgr_config[['assay_props']][['vcf_tumor_only']] <- TRUE
  pcgr_config[['assay_props']][['mode']] <- 'Tumor-Only'
  if(as.integer(args[17]) == 1){
    pcgr_config[['assay_props']][['mode']] <- 'Cell line (Tumor-Only)'
  }
}

pcgr_config[['tumor_only']] <- list()

arg_counter <- 18
for(maf_pop in c('maf_onekg_afr','maf_onekg_amr','maf_onekg_eas',
                 'maf_onekg_eur','maf_onekg_sas','maf_onekg_global')){
  pcgr_config[['tumor_only']][[maf_pop]] <- 
    as.numeric(args[arg_counter])
  arg_counter <- arg_counter + 1
}

for(maf_pop in c('maf_gnomad_afr','maf_gnomad_amr','maf_gnomad_asj',
                 'maf_gnomad_eas','maf_gnomad_fin','maf_gnomad_nfe',
                 'maf_gnomad_other','maf_gnomad_sas','maf_gnomad_global')){
  pcgr_config[['tumor_only']][[maf_pop]] <- 
    as.numeric(args[arg_counter])
  arg_counter <- arg_counter + 1
  
}

for(maf_filter in c('exclude_pon','exclude_likely_het_germline',
                    'exclude_likely_hom_germline','exclude_dbsnp_nonsomatic',
                    'exclude_nonexonic')){
  pcgr_config[['tumor_only']][[maf_filter]] <- 
    as.logical(as.integer(args[arg_counter]))
  arg_counter <- arg_counter + 1
  
}

## Settings for analyses to be performed and included in report
## TMB
pcgr_config[['tmb']] <- list()
pcgr_config[['tmb']][['run']] <- as.logical(as.integer(args[arg_counter]))
arg_counter <- arg_counter + 1
pcgr_config[['tmb']][['algorithm']] <- as.character(args[arg_counter])
arg_counter <- arg_counter + 1

## MSI
pcgr_config[['msi']] <- list()
pcgr_config[['msi']][['run']] <- as.logical(as.integer(args[arg_counter]))
arg_counter <- arg_counter + 1

## Mutational signatures
pcgr_config[['msigs']] <- list()
pcgr_config[['msigs']][['run']] <- as.logical(as.integer(args[arg_counter]))
arg_counter <- arg_counter + 1

pcgr_config[['msigs']][['mutation_limit']] <- as.numeric(args[arg_counter])
arg_counter <- arg_counter + 1

pcgr_config[['msigs']][['all_reference_signatures']] <- as.logical(as.integer(args[arg_counter]))
arg_counter <- arg_counter + 1

pcgr_config[['msigs']][['include_artefact_signatures']] <- as.logical(as.integer(args[arg_counter]))
arg_counter <- arg_counter + 1

## CNA arguments
## Copy number aberration (CNA)
pcgr_config[['cna']] <- list()
for(cna_config in c('log_r_homdel','log_r_gain','cna_overlap_pct')){
  pcgr_config[['cna']][[cna_config]] <- as.numeric(args[arg_counter])
  arg_counter <- arg_counter + 1
}

## Allelic support settings (VCF) - max/min support
for(as_setting in c('tumor_af_min','tumor_dp_min',
                    'control_dp_min','control_af_max')){
  pcgr_config[['allelic_support']][[as_setting]] <- as.numeric(args[arg_counter])
  arg_counter <- arg_counter + 1
}

## Allelic support settings (VCF) - VCF INFO tags
for(as_setting in c('tumor_af_tag','tumor_dp_tag',
                    'control_af_tag','control_dp_tag',
                    'call_conf_tag')){
  pcgr_config[['allelic_support']][[as_setting]] <- as.character(args[arg_counter])
  arg_counter <- arg_counter + 1
}

## Other arguments
## Other settings (VEP, vcf2maf, visual theme, VCF validation, custom tags)
pcgr_config[['clinicaltrials']] <- list()
pcgr_config[['clinicaltrials']][['run']] <- as.logical(as.integer(args[arg_counter]))
arg_counter <- arg_counter + 1
pcgr_config[['other']] <- list()
pcgr_config[['other']][['vep_n_forks']] <- as.integer(args[arg_counter])
arg_counter <- arg_counter + 1
pcgr_config[['other']][['vep_buffer_size']] <- as.integer(args[arg_counter])
arg_counter <- arg_counter + 1
pcgr_config[['other']][['vep_no_intergenic']] <- as.logical(as.integer(args[arg_counter]))
arg_counter <- arg_counter + 1
pcgr_config[['other']][['vep_pick_order']] <- as.character(args[arg_counter])
arg_counter <- arg_counter + 1
pcgr_config[['other']][['vep_regulatory']] <- as.logical(as.integer(args[arg_counter]))
arg_counter <- arg_counter + 1
pcgr_config[['other']][['vcf2maf']] <- as.logical(as.integer(args[arg_counter]))
arg_counter <- arg_counter + 1
pcgr_config[['other']][['list_noncoding']] <- as.logical(as.integer(args[arg_counter]))
arg_counter <- arg_counter + 1
#pcgr_config[['preserved_info_tags']] <- list()
pcgr_config[['preserved_info_tags']] <- as.character(args[arg_counter])
arg_counter <- arg_counter + 1
pcgr_config[['visual']] <- list()
pcgr_config[['visual']][['report_theme']] <- as.character(args[arg_counter])
arg_counter <- arg_counter + 1
pcgr_config[['visual']][['nonfloating_toc']] <- as.logical(as.integer(args[arg_counter]))
arg_counter <- arg_counter + 1
pcgr_config[['other']][['vcf_no_validation']] <- as.logical(as.integer(args[arg_counter]))

saveRDS(pcgr_config, file=paste0(dir,"/", pcgr_config[['required_args']][['sample_name']],".pcgr_config.rds"))

pcgr_data <- readRDS(paste0(pcgr_config[['required_args']][['data_dir']],
                            '/data/',
                            pcgr_config[['required_args']][['genome_assembly']],
                            '/rds/pcgr_data.rds'))

pcgr_data[['assembly']][['seqinfo']] <- 
   GenomeInfoDb::Seqinfo(seqnames = GenomeInfoDb::seqlevels(GenomeInfoDb::seqinfo(BSgenome.Hsapiens.UCSC.hg38)), 
                         seqlengths = GenomeInfoDb::seqlengths(GenomeInfoDb::seqinfo(BSgenome.Hsapiens.UCSC.hg38)), 
                         genome = 'hg38')
pcgr_data[['assembly']][['bsg']] <- BSgenome.Hsapiens.UCSC.hg38
if(pcgr_config[['required_args']][['genome_assembly']] == 'grch37'){
  pcgr_data[['assembly']][['bsg']] <- BSgenome.Hsapiens.UCSC.hg19
  pcgr_data[['assembly']][['seqinfo']] <- 
     GenomeInfoDb::Seqinfo(
       seqnames = GenomeInfoDb::seqlevels(GenomeInfoDb::seqinfo(BSgenome.Hsapiens.UCSC.hg19)), 
       seqlengths = GenomeInfoDb::seqlengths(GenomeInfoDb::seqinfo(BSgenome.Hsapiens.UCSC.hg19)), genome = 'hg19')
}

if(pcgr_config[['other']][['vep_regulatory']] == F){
  for(e in c('tier4_display','tier5_display','all','tsv')){
    pcgr_data[['annotation_tags']][[e]] <-
      pcgr_data[['annotation_tags']][[e]][
        pcgr_data[['annotation_tags']][[e]] != "REGULATORY_ANNOTATION"]
  }
}

my_log4r_layout <- function(level, ...) {
  paste0(format(Sys.time()), " - pcgr-report-generation - ", level, " - ", ..., "\n", collapse = "")
}

log4r_logger <- log4r::logger(threshold = "INFO", appenders = log4r::console_appender(my_log4r_layout))

## Clinical trials 
if(pcgr_config[['t_props']][['tumor_type']] == "Cancer, NOS"){
  pcgrr:::log4r_info(paste0("Clinical trials will not be included in the report when primary site is not specified - skipping"))
  pcgr_config[['clinicaltrials']][['run']] <- F
}

pcgrr::log4r_info(paste0("Tumor primary site: ",pcgr_config[['t_props']][['tumor_type']]))

pcg_report <- NULL

# ## Generate report object
pcg_report <- 
  pcgrr::generate_pcgr_report(
    project_directory  = pcgr_config[['required_args']][['output_dir']], 
    pcgr_data = pcgr_data, 
    config = pcgr_config, 
    tier_model = 'pcgr_acmg')

# ## Write report and result files
if(!is.null(pcg_report)){
  pcgrr::write_report_output(
    pcg_report, 
    pcgr_config,
    output_format = 'json')
  pcgrr::write_report_output(
    pcg_report, 
    pcgr_config,
    output_format = 'snv_tsv')
  pcgrr::write_report_output(
    pcg_report, 
    pcgr_config,
    output_format = 'msigs_tsv')
  if(pcgr_config[['assay_props']][['vcf_tumor_only']] == T){
    pcgrr::write_report_output(
      pcg_report, 
      pcgr_config,
      output_format = 'snv_tsv_unfiltered')
  }
  pcgrr::write_report_output(
    pcg_report, 
    pcgr_config,
    output_format = 'cna_tsv')
  pcgrr::write_report_output(
    pcg_report, 
    pcgr_config,
    output_format = 'html')
  pcgrr::write_report_output(
    pcg_report, 
    pcgr_config,
    output_format = 'html', 
    flexdb = T)

}

