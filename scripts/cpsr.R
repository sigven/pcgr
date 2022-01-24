#!/usr/bin/env Rscript

options(warn=-1)
.libPaths(R.home("library")) # use conda R pkgs, not e.g. user's local installation

args = commandArgs(trailingOnly=TRUE)

suppressWarnings(suppressPackageStartupMessages(library(pcgrr)))
suppressWarnings(suppressPackageStartupMessages(library(GenomeInfoDb)))
suppressWarnings(suppressPackageStartupMessages(library(log4r)))

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
cpsr_config[['other']][['vep_gencode_all']] <- as.logical(args[22])
cpsr_config[['other']][['vep_skip_intergenic']] <- as.logical(as.integer(args[23]))
cpsr_config[['other']][['vep_regulatory']] <- as.logical(as.integer(args[24]))

arg_counter <- 25
for (opt in c('secondary_findings','classify_all','ignore_noncoding',
              'clinvar_ignore_noncancer','diagnostic_grade_only')){
  cpsr_config[opt] <- as.logical(as.integer(args[arg_counter]))
  arg_counter <- arg_counter + 1
}

saveRDS(cpsr_config,
        file = file.path(dir, cpsr_config[['required_args']][['sample_name']],".cpsr_config.rds"))

pcgr_data <- readRDS(
  file.path(
    cpsr_config[['required_args']][['data_dir']],
    'data', cpsr_config[['required_args']][['genome_assembly']],
    'rds', 'pcgr_data.rds')
  )

## temporary type fix
pcgr_data$biomarkers$cgi$ACTIONABILITY_SCORE <-
  as.numeric(pcgr_data$biomarkers$cgi$ACTIONABILITY_SCORE)
for(col in c('VARIANT_TYPE','DRUG_INTERACTION_TYPE','GDNA')){
  pcgr_data$biomarkers$cgi[,col] <- as.character(
    pcgr_data$biomarkers$cgi[,col]
  )
}

# set up genome assembly
genome_assembly <- cpsr_config[['required_args']][['genome_assembly']]
bsgenome_obj <- pcgrr::get_genome_obj(genome_assembly)
genome_grch2hg <- c("grch38" = "hg38", "grch37" = "hg19")
pcgr_data[['assembly']][['seqinfo']] <-
  GenomeInfoDb::Seqinfo(
    seqnames = GenomeInfoDb::seqlevels(GenomeInfoDb::seqinfo(bsgenome_obj)),
    seqlengths = GenomeInfoDb::seqlengths(GenomeInfoDb::seqinfo(bsgenome_obj)),
    genome = genome_grch2hg[genome_assembly])
pcgr_data[['assembly']][['bsg']] <- bsgenome_obj

my_log4r_layout <- function(level, ...) {
  paste0(format(Sys.time()), " - cpsr-report-generation - ", level, " - ", ..., "\n", collapse = "")
}

log4r_logger <- log4r::logger(threshold = "INFO", appenders = log4r::console_appender(my_log4r_layout))

# this gets passed on to all the log4r_* functions inside the pkg
options("PCGRR_LOG4R_LOGGER" = log4r_logger)

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
