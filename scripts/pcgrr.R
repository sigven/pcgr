#!/usr/bin/env Rscript

.libPaths(R.home("library")) # use conda R pkgs, not e.g. user's local installation

suppressWarnings(suppressPackageStartupMessages(library(argparse)))
suppressWarnings(suppressPackageStartupMessages(library(pcgrr)))
suppressWarnings(suppressPackageStartupMessages(library(GenomeInfoDb)))
suppressWarnings(suppressPackageStartupMessages(library(stringr)))

##---- Argument Parsing ----##
p <- argparse::ArgumentParser(description='PCGR HTML generation step', prog='pcgrr')

# required args
required_args <- c(
  'output_dir', 'query_vcf2tsv', 'query_cna', 'query_rna_fusion',
  'query_rna_expression', 'cpsr_report', 'sample_name',
  'pcgr_version', 'genome_assembly', 'data_dir')
p$add_argument('output_dir')           #  1
p$add_argument('query_vcf2tsv')        #  2
p$add_argument('query_cna')            #  3
p$add_argument('query_rna_fusion')     #  4
p$add_argument('query_rna_expression') #  5
p$add_argument('cpsr_report')          #  6
p$add_argument('sample_name')          #  7
p$add_argument('pcgr_version')         #  8
p$add_argument('genome_assembly')      #  9
p$add_argument('data_dir')             # 10
# tumor props
t_props_args <- c('tumor_purity', 'tumor_ploidy', 'tumor_type')
p$add_argument('tumor_purity')         # 11
p$add_argument('tumor_ploidy')         # 12
p$add_argument('tumor_type')           # 13
# assay props
assay_props_args <- c('target_size_mb', 'type')
p$add_argument('target_size_mb', type='double')     # 14
p$add_argument('type')                              # 15
# used to populate 'vcf_tumor_only' and 'mode' in the 'assay_props' list
p$add_argument('vcf_tumor_only1')                   # 16
p$add_argument('mode1')                             # 17
# tumor only args
tumor_only_args <- list(
  onekg = c('maf_onekg_afr', 'maf_onekg_amr', 'maf_onekg_eas',
            'maf_onekg_eur', 'maf_onekg_sas', 'maf_onekg_global'),
  gnomad = c('maf_gnomad_afr','maf_gnomad_amr','maf_gnomad_asj',
             'maf_gnomad_eas','maf_gnomad_fin','maf_gnomad_nfe',
             'maf_gnomad_other','maf_gnomad_sas','maf_gnomad_global'),
  filter = c('exclude_pon', 'exclude_likely_het_germline',
             'exclude_likely_hom_germline', 'exclude_dbsnp_nonsomatic',
             'exclude_nonexonic'))
p$add_argument('maf_onekg_afr', type='double')      # 18
p$add_argument('maf_onekg_amr', type='double')      # 19
p$add_argument('maf_onekg_eas', type='double')      # 20
p$add_argument('maf_onekg_eur', type='double')      # 21
p$add_argument('maf_onekg_sas', type='double')      # 22
p$add_argument('maf_onekg_global', type='double')   # 23
p$add_argument('maf_gnomad_afr', type='double')     # 24
p$add_argument('maf_gnomad_amr', type='double')     # 25
p$add_argument('maf_gnomad_asj', type='double')     # 26
p$add_argument('maf_gnomad_eas', type='double')     # 27
p$add_argument('maf_gnomad_fin', type='double')     # 28
p$add_argument('maf_gnomad_nfe', type='double')     # 29
p$add_argument('maf_gnomad_other', type='double')   # 30
p$add_argument('maf_gnomad_sas', type='double')     # 31
p$add_argument('maf_gnomad_global', type='double')  # 32
p$add_argument('exclude_pon', type='integer')                 # 33
p$add_argument('exclude_likely_het_germline', type='integer') # 34
p$add_argument('exclude_likely_hom_germline', type='integer') # 35
p$add_argument('exclude_dbsnp_nonsomatic', type='integer')    # 36
p$add_argument('exclude_nonexonic', type='integer')           # 37
# tmb args
p$add_argument('tmb_run', type='integer') # 38
p$add_argument('tmb_algo')                # 39
# msi/msigs args
p$add_argument('msi_run', type='integer') # 40
p$add_argument('msigs_run', type='integer')           # 41
p$add_argument('msigs_mut_lim', type='double')        # 42
p$add_argument('msigs_all_ref_sigs', type='integer')  # 43
p$add_argument('msigs_incl_art_sigs', type='integer') # 44
# cna args
cna_args <- c('log_r_homdel', 'log_r_gain', 'cna_overlap_pct')
p$add_argument('log_r_homdel', type='double')      # 45
p$add_argument('log_r_gain', type='double')        # 46
p$add_argument('cna_overlap_pct', type='double')   # 47
# allelic support args
allelic_support_args <- c(
  'tumor_af_min', 'tumor_dp_min',
  'control_dp_min', 'control_af_max',
  'tumor_af_tag', 'tumor_dp_tag',
  'control_af_tag', 'control_dp_tag',
  'call_conf_tag')
p$add_argument('tumor_af_min', type='double')   # 48
p$add_argument('tumor_dp_min', type='double')   # 49
p$add_argument('control_dp_min', type='double') # 50
p$add_argument('control_af_max', type='double') # 51
p$add_argument('tumor_af_tag')    # 52
p$add_argument('tumor_dp_tag')    # 53
p$add_argument('control_af_tag')  # 54
p$add_argument('control_dp_tag')  # 55
p$add_argument('call_conf_tag')   # 56
# clinicaltrials
p$add_argument('clinicaltrials_run', type='integer')  # 57
# other
p$add_argument('vep_n_forks', type='integer')         # 58
p$add_argument('vep_buffer_size', type='integer')     # 59
p$add_argument('vep_no_intergenic', type='integer')   # 60
p$add_argument('vep_pick_order')                      # 61
p$add_argument('vep_regulatory', type='integer')      # 62
p$add_argument('vcf2maf', type='integer')             # 63
p$add_argument('list_noncoding', type='integer')      # 64
# preserved_info_tags
p$add_argument('preserved_info_tags')                 # 65
# visual
p$add_argument('report_theme')                        # 66
p$add_argument('nonfloating_toc', type='integer')     # 67
# other
p$add_argument('vcf_no_validation', type='integer')   # 68

args <- p$parse_args()

# main pcgr_config list (processed later)
pcgr_config <- list(
  required_args = args[required_args],
  t_props = args[t_props_args],
  assay_props = args[assay_props_args],
  tumor_only = args[unlist(tumor_only_args)],
  tmb = list(run = as.logical(args[['tmb_run']]), algorithm = args[['tmb_algo']]),
  msi = list(run = as.logical(args[['msi_run']])),
  msigs = list(run = as.logical(args[['msigs_run']]),
               mutation_limit = args[['msigs_mut_lim']],
               all_reference_signatures = as.logical(args[['msigs_all_ref_sigs']]),
               include_artefact_signatures = as.logical(args[['msigs_incl_art_sigs']])
  ),
  cna = args[cna_args],
  allelic_support = args[allelic_support_args],
  clinicaltrials = list(run = as.logical(args[['clinicaltrials_run']])),
  other = list(vep_n_forks = args[['vep_n_forks']],
               vep_buffer_size = args[['vep_buffer_size']],
               vep_no_intergenic = as.logical(args[['vep_no_intergenic']]),
               vep_pick_order = args[['vep_pick_order']],
               vep_regulatory = as.logical(args[['vep_regulatory']]),
               vcf2maf = as.logical(args[['vcf2maf']]),
               list_noncoding = as.logical(args[['list_noncoding']]),
               vcf_no_validation = as.logical(args[['vcf_no_validation']])
  ),
  preserved_info_tags = args[['preserved_info_tags']],
  visual = list(report_theme = args[['report_theme']],
                nonfloating_toc = as.logical(args[['nonfloating_toc']])
  )
)

##---- Argument Processing ----##
# required args
if (pcgr_config[['required_args']][['query_cna']] == "None"){
  pcgr_config[['required_args']][['query_cna']] <- NULL
}
if (pcgr_config[['required_args']][['query_rna_fusion']] == "None"){
  pcgr_config[['required_args']][['query_rna_fusion']] <- NULL
}
if (pcgr_config[['required_args']][['query_rna_expression']] == "None"){
  pcgr_config[['required_args']][['query_rna_expression']] <- NULL
}
if (pcgr_config[['required_args']][['cpsr_report']] == "None"){
  pcgr_config[['required_args']][['cpsr_report']] <- NULL
}

# tumor props
# Handle case when 'NA'
purity <- pcgr_config[['t_props']][['tumor_purity']]
ploidy <- pcgr_config[['t_props']][['tumor_ploidy']]
pcgr_config[['t_props']][['tumor_purity']] <-
  ifelse(purity == 'NA', NA_real_, as.numeric(purity))
pcgr_config[['t_props']][['tumor_ploidy']] <-
  ifelse(ploidy == 'NA', NA_real_, as.numeric(ploidy))
pcgr_config[['t_props']][['tumor_type']] <- stringr::str_replace_all(
  stringr::str_replace_all(
    args[['tumor_type']], "_", " "),
  "@", "/")
if (pcgr_config[['t_props']][['tumor_type']] == "Any") {
  pcgr_config[['t_props']][['tumor_type']] <- "Cancer, NOS"
}

### Sequencing assay properties (VCF)
### Target (WES/WGS/TARGETED), mode (tumor-control/tumor-only), coding target size
pcgr_config[['assay_props']][['mode']] <- 'Tumor-Control'
pcgr_config[['assay_props']][['vcf_tumor_only']] <- FALSE

if (as.integer(args[['vcf_tumor_only1']]) == 1) {
  pcgr_config[['assay_props']][['vcf_tumor_only']] <- TRUE
  pcgr_config[['assay_props']][['mode']] <- 'Tumor-Only'
  if (as.integer(args[['mode1']]) == 1){
    pcgr_config[['assay_props']][['mode']] <- 'Cell line (Tumor-Only)'
  }
}

# tumor_only maf filter
for (mf in tumor_only_args[['filter']]) {
  pcgr_config[['tumor_only']][[mf]] <- as.logical(args[[mf]])
}

pcgr_config_rds <- file.path(
  pcgr_config[['required_args']][['output_dir']],
  paste0(pcgr_config[['required_args']][['sample_name']],
         ".pcgr_config.rds"))
saveRDS(pcgr_config, file = pcgr_config_rds)
### Arg processing END

pcgr_data <- readRDS(
  file.path(pcgr_config[['required_args']][['data_dir']],
            'data',
            pcgr_config[['required_args']][['genome_assembly']],
            'rds/pcgr_data.rds'))

# set up genome assembly
genome_assembly <- pcgr_config[['required_args']][['genome_assembly']]
bsgenome_obj <- pcgrr::get_genome_obj(genome_assembly)
genome_grch2hg <- c("grch38" = "hg38", "grch37" = "hg19")
pcgr_data[['assembly']][['seqinfo']] <-
  GenomeInfoDb::Seqinfo(
    seqnames = GenomeInfoDb::seqlevels(GenomeInfoDb::seqinfo(bsgenome_obj)),
    seqlengths = GenomeInfoDb::seqlengths(GenomeInfoDb::seqinfo(bsgenome_obj)),
    genome = genome_grch2hg[genome_assembly])
pcgr_data[['assembly']][['bsg']] <- bsgenome_obj

if (pcgr_config[['other']][['vep_regulatory']] == F){
  for (e in c('tier4_display','tier5_display','all','tsv')){
    pcgr_data[['annotation_tags']][[e]] <-
      pcgr_data[['annotation_tags']][[e]][
        pcgr_data[['annotation_tags']][[e]] != "REGULATORY_ANNOTATION"]
  }
}

my_log4r_layout <- function(level, ...) {
  paste0(format(Sys.time()), " - pcgr-report-generation - ",
         level, " - ", ..., "\n", collapse = "")
}

log4r_logger <- log4r::logger(threshold = "INFO",
                              appenders = log4r::console_appender(my_log4r_layout))

# this gets passed on to all the log4r_* functions inside the pkg
options("PCGRR_LOG4R_LOGGER" = log4r_logger)

## Clinical trials
if (pcgr_config[['t_props']][['tumor_type']] == "Cancer, NOS"){
  pcgrr:::log4r_info(paste0("Clinical trials will not be included in the report when primary site is not specified - skipping"))
  pcgr_config[['clinicaltrials']][['run']] <- F
}

pcgrr:::log4r_info(paste0("Tumor primary site: ", pcgr_config[['t_props']][['tumor_type']]))

pcg_report <- NULL

defaultW <- getOption("warn")
options(warn = -1)



# ## Generate report object
pcg_report <-
  pcgrr::generate_pcgr_report(
    project_directory  = pcgr_config[['required_args']][['output_dir']],
    pcgr_data = pcgr_data,
    config = pcgr_config,
    tier_model = 'pcgr_acmg')

options(warn = defaultW)


# ## Write report and result files
if (!is.null(pcg_report)) {
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
  if (pcgr_config[['assay_props']][['vcf_tumor_only']] == T){
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
