#!/usr/bin/env Rscript

options(warn=-1)
.libPaths(R.home("library"))

args = commandArgs(trailingOnly=TRUE)

suppressWarnings(suppressPackageStartupMessages(library(pcgrr)))
suppressWarnings(suppressPackageStartupMessages(library(magrittr)))
suppressWarnings(suppressPackageStartupMessages(library(BSgenome.Hsapiens.UCSC.hg19)))
suppressWarnings(suppressPackageStartupMessages(library(BSgenome.Hsapiens.UCSC.hg38)))

dir <- as.character(args[1])

cpsr_config <- list()
cpsr_config[['required_args']] <- list()
cpsr_config[['required_args']][['query_vcf2tsv']] <- as.character(args[2])
cpsr_config[['required_args']][['sample_name']] <- as.character(args[3])
cpsr_config[['required_args']][['pcgr_version']] <- as.character(args[4])
cpsr_config[['required_args']][['cpsr_version']] <- as.character(args[5])
cpsr_config[['required_args']][['genome_assembly']] <- as.character(args[6])
cpsr_config[['required_args']][['data_dir']] <- as.character(args[7])
cpsr_config[['required_args']][['output_dir']] <- as.character(dir)
cpsr_config[['required_args']][['virtual_panel_id']] <- as.character(args[8])
cpsr_config[['preserved_info_tags']] <- as.character(args[9])

cpsr_config[['custom_panel']] <- list()
cpsr_config[['custom_panel']][['bed_fname']] <- as.character(args[10])
cpsr_config[['custom_panel']][['name']] <- as.character(args[11])
cpsr_config[['custom_panel']][['version']] <- '1.0'
cpsr_config[['custom_panel']][['url']] <- ''

cpsr_config[['assay_props']] <- list()
cpsr_config[['assay_props']][['type']] <- "Germline"

cpsr_config[['visual']] <- list()
cpsr_config[['visual']][['report_theme']] <- as.character(args[12])
cpsr_config[['visual']][['table_display']] <- as.character(args[13])
cpsr_config[['visual']][['nonfloating_toc']] <- as.logical(as.integer(args[14]))
cpsr_config[['gwas']] <- list()
cpsr_config[['gwas']][['run']] <- as.logical(as.integer(args[15]))
cpsr_config[['gwas']][['p_value_min']] <- as.numeric(args[16])
cpsr_config[['popgen']] <- list()
cpsr_config[['popgen']][['pop_gnomad']] <- as.character(args[17])
cpsr_config[['popgen']][['maf_upper_threshold']] <- as.numeric(args[18])

cpsr_config[['other']] <- list()
cpsr_config[['other']][['vep_pick_order']] <- as.character(args[19])
cpsr_config[['other']][['vep_n_forks']] <- as.integer(args[20])
cpsr_config[['other']][['vep_buffer_size']] <- as.integer(args[21])
cpsr_config[['other']][['vep_skip_intergenic']] <- as.logical(as.integer(args[22]))
cpsr_config[['other']][['vep_regulatory']] <- as.logical(as.integer(args[23]))

arg_counter <- 24
for (opt in c('secondary_findings','classify_all','ignore_noncoding',
              'clinvar_ignore_noncancer','diagnostic_grade_only')){
  cpsr_config[opt] <- as.logical(as.integer(args[arg_counter]))
  arg_counter <- arg_counter + 1
}

saveRDS(cpsr_config, file=paste0(dir,"/", cpsr_config[['required_args']][['sample_name']],".cpsr_config.rds"))

pcgr_data <- readRDS(paste0(cpsr_config[['required_args']][['data_dir']],
                            '/data/',
                            cpsr_config[['required_args']][['genome_assembly']],
                            '/rds/pcgr_data.rds'))

pcgr_data[['assembly']][['seqinfo']] <- 
   GenomeInfoDb::Seqinfo(seqnames = GenomeInfoDb::seqlevels(GenomeInfoDb::seqinfo(BSgenome.Hsapiens.UCSC.hg38)), 
                         seqlengths = GenomeInfoDb::seqlengths(GenomeInfoDb::seqinfo(BSgenome.Hsapiens.UCSC.hg38)), genome = 'hg38')
pcgr_data[['assembly']][['bsg']] <- BSgenome.Hsapiens.UCSC.hg38
if(cpsr_config[['required_args']][['genome_assembly']] == 'grch37'){
  pcgr_data[['assembly']][['bsg']] <- BSgenome.Hsapiens.UCSC.hg19
  pcgr_data[['assembly']][['seqinfo']] <- 
     GenomeInfoDb::Seqinfo(seqnames = GenomeInfoDb::seqlevels(GenomeInfoDb::seqinfo(BSgenome.Hsapiens.UCSC.hg19)), 
                           seqlengths = GenomeInfoDb::seqlengths(GenomeInfoDb::seqinfo(BSgenome.Hsapiens.UCSC.hg19)), genome = 'hg19')
}

my_log4r_layout <- function(level, ...) {
  paste0(format(Sys.time()), " - cpsr-report-generation - ", level, " - ", ..., "\n", collapse = "")
}

log4r_logger <- log4r::logger(threshold = "INFO", appenders = log4r::console_appender(my_log4r_layout))


## Generate report content
cps_report <- pcgrr::generate_cpsr_report(
  project_directory = cpsr_config[['required_args']][['output_dir']], 
  pcgr_data, 
  cpsr_config
)

# ## Write report contents to output files
if(!is.null(cps_report)){
  pcgrr::write_report_output(
    cps_report, 
    cpsr_config,
    tier_model = 'cpsr', 
    output_format = 'snv_tsv')
  pcgrr::write_report_output(
    cps_report, 
    cpsr_config,
    tier_model = "cpsr", 
    output_format = 'json')
  pcgrr::write_report_output(
    cps_report, 
    cpsr_config,
    tier_model = "cpsr", 
    output_format = 'html')
}

