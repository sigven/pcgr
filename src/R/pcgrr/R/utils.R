
check_common_colnames <- function(df1 = NULL, df2 = NULL, cnames = NULL){
  invisible(assertthat::assert_that(is.data.frame(df1), msg = "Object df1 is not of type data.frame"))
  invisible(assertthat::assert_that(is.data.frame(df2), msg = "Object df2 is not of type data.frame"))

  existing_common_columns <- F
  colnames_in_df1 <- unique(cnames %in% colnames(df1))
  colnames_in_df2 <- unique(cnames %in% colnames(df2))

  if(length(colnames_in_df1) == 1 & length(colnames_in_df2) == 1){
    if(colnames_in_df1 == T & colnames_in_df2 == T){
      existing_common_columns <- T
    }
  }
  return(existing_common_columns)

}

remove_cols_from_df <- function(df, cnames = NULL){

  invisible(assertthat::assert_that(is.data.frame(df), msg = "Object 'df' is not of type data.frame"))
  if(!is.null(cnames)){
    for(c in cnames){
      col_name <- as.character(c)
      if(col_name %in% colnames(df)){
        df[,col_name] <- NULL
      }
    }
  }
  return(df)

}


#' Function that plots a histogram of the the variant allelic support (tumor) - grouped by tiers
#'
#' @param tier_df data frame with somatic mutations
#' @param bin_size size of bins for allelic frequency
#' @return p geom_histogram plot from ggplot2
#'
#'
tier_af_distribution <- function(tier_df, bin_size = 0.1) {

  af_bin_df <- data.frame()
  assertable::assert_colnames(tier_df, c("AF_TUMOR"), only_colnames = F, quiet = T)
  i <- 1
  num_bins <- as.integer(1 / bin_size)
  bin_start <- 0
  while (i <= num_bins) {
    bin_end <- bin_start + bin_size
    bin_name <- as.character(paste0(bin_start, " - ", bin_end))
    j <- 1
    while (j <= 4) {
      TIER <- paste0("TIER ", j)
      df <- data.frame(bin_name = bin_name,
                       bin_start = bin_start,
                       bin_end = bin_end,
                       bin = as.integer(i),
                       TIER = TIER, stringsAsFactors = F)
      af_bin_df <- rbind(af_bin_df, df)
      j <- j + 1
    }
    TIER <- "NONCODING"
    df <- data.frame(bin_name = bin_name,
                     bin_start = bin_start,
                     bin_end = bin_end,
                     bin = as.integer(i),
                     TIER = TIER, stringsAsFactors = F)
    af_bin_df <- rbind(af_bin_df, df)
    bin_start <- bin_end
    i <- i + 1
  }

  tier_df_trans <- transform(tier_df, bin = cut(AF_TUMOR, breaks = seq(from = 0, to = 1, by = bin_size),
                                                right = F, include.lowest = T, labels = F))
  tier_df_trans_bin <- as.data.frame(dplyr::group_by(tier_df_trans, TIER, bin) %>%
                                       dplyr::summarise(Count = dplyr::n()))
  af_bin_df <- af_bin_df %>%
    dplyr::left_join(tier_df_trans_bin, by = c("bin", "TIER")) %>%
    dplyr::mutate(Count = dplyr::if_else(is.na(Count), as.numeric(0), as.numeric(Count)))

  return(af_bin_df)

}

#' Function that check for valid mutations/chromosomes in input
#'
#' @param vcf_data_df data frame
#' @param chromosome_column name of columns for which string replace is to be performed
#' @param pattern pattern to replace
#' @return vcf_data_df_valid data frame with valid mutations
#'
get_valid_chromosomes <- function(vcf_data_df, chromosome_column = "CHROM", bsg = BSgenome.Hsapiens.UCSC.hg19) {
  assertthat::assert_that(is.data.frame(vcf_data_df), msg = paste0("Argument 'vcf_data_df' must be of type data.frame, not ", class(vcf_data_df)))
  assertable::assert_colnames(vcf_data_df, chromosome_column, only_colnames = F, quiet = T)
  vcf_data_df_valid <- vcf_data_df
  vcf_data_df_valid[, chromosome_column] <- factor(vcf_data_df_valid[, chromosome_column])
  levels(vcf_data_df_valid[, chromosome_column]) <- sub("^([0-9XY])", "chr\\1", levels(vcf_data_df_valid[, chromosome_column]))
  levels(vcf_data_df_valid[, chromosome_column]) <- sub("^MT", "chrM", levels(vcf_data_df_valid[, chromosome_column]))
  levels(vcf_data_df_valid[, chromosome_column]) <- sub("^(GL[0-9]+).[0-9]", "chrUn_\\L\\1", levels(vcf_data_df_valid[, chromosome_column]), perl = T)
  unknown.regions <- levels(vcf_data_df_valid[, chromosome_column])[which(!(levels(vcf_data_df_valid[, chromosome_column]) %in% GenomeInfoDb::seqnames(bsg)))]
  if (length(unknown.regions) > 0) {
    unknown.regions <- paste(unknown.regions, collapse = ",\ ")
    rlogging::warning(paste("Check chr names -- not all match BSgenome.Hsapiens object:\n", unknown.regions, sep = " "))
    vcf_data_df_valid <- vcf_data_df_valid[vcf_data_df_valid[, chromosome_column] %in% GenomeInfoDb::seqnames(bsg), ]
  }
  vcf_data_df_valid[, chromosome_column] <- as.character(vcf_data_df_valid[, chromosome_column])
  return(vcf_data_df_valid)

}

#' Function that excludes genomic aberrations from non-nuclear chromosomes
#'
#' @param vcf_df data frame
#' @param chrom_var variable name of chromosome in data frame
#' @return vcf_df data frame with mutations from nuclear chromosomes only
#'
get_ordinary_chromosomes <- function(vcf_df, chrom_var = "CHROM") {
  invisible(assertthat::assert_that(is.data.frame(vcf_df), msg = "Argument 'vcf_df' must be of type data.frame"))
  assertable::assert_colnames(vcf_df, chrom_var, only_colnames = F, quiet = T)
  vcf_df <- vcf_df %>%
    dplyr::mutate(!!rlang::sym(chrom_var) := as.character(!!rlang::sym(chrom_var)))
  n_before_exclusion <- nrow(vcf_df)
  nuc_chromosomes_df <- data.frame(c(as.character(seq(1:22)), "X", "Y"), stringsAsFactors = F)
  colnames(nuc_chromosomes_df) <- c(chrom_var)
  vcf_df <- dplyr::semi_join(vcf_df, nuc_chromosomes_df, by = chrom_var)
  n_after_exclusion <- nrow(vcf_df)
  rlogging::message(paste0("Excluding ", n_before_exclusion - n_after_exclusion, " variants from non-nuclear chromosomes/scaffolds"))
  return(vcf_df)

}

#' Function that orders genomic aberrations according to order of chromosomes and chromosomal position
#'
#' @param vcf_df data frame
#' @param chrom_var variable name of chromosome in data frame
#' @param pos_var variable name for chromosomal position
#' @return vcf_df data frame with ordered mutations
#'
order_variants <- function(vcf_df, chrom_var = "CHROM", pos_var = "POS") {
  stopifnot(is.data.frame(vcf_df) & chrom_var %in% colnames(vcf_df) & pos_var %in% colnames(vcf_df))
  if (nrow(vcf_df) == 0)return(vcf_df)
  vcf_df %>%
    dplyr::mutate(!!rlang::sym(chrom_var) := factor(!!rlang::sym(chrom_var), ordered = T, levels = c(as.character(seq(1:22)), "X", "Y"))) %>%
    dplyr::arrange(!!rlang::sym(chrom_var), !!rlang::sym(pos_var)) %>%
    dplyr::mutate(!!rlang::sym(chrom_var) := as.character(!!rlang::sym(chrom_var)))
}


#' Function that sorts chromosomal segments according to chromosome and chromosomal start/end position
#'
#' @param df data frame with chromosome and start + end segment
#' @param chromosome_column name of column for chromosome name is sigven
#' @param start_segment name of column that indicates start of chromosomal segment
#' @param end_segment name of column that indicates end of chromosomal segment
#' @return df_final data frame with sorted chromosomal segments
#'
#'
sort_chromosomal_segments <- function(df, chromosome_column = "CHROM",
                                      start_segment = "POS", end_segment = "POS") {

  invisible(assertthat::assert_that(!is.null(df), msg = "Argument 'df' must be a non-NULL object"))
  invisible(assertthat::assert_that(is.data.frame(df), msg = paste0("Argument 'df' must be of type data.frame, not ", class(df))))
  assertable::assert_colnames(df, c(chromosome_column, start_segment, end_segment), only_colnames = F, quiet = T)
  if (nrow(df) == 0) {
    return(df)
  }


  df[, start_segment] <- as.integer(df[, start_segment])
  df[, end_segment] <- as.integer(df[, end_segment])
  df_sorted <- df

  chr_prefix <- FALSE
  chromosome_names <- unique(df[, chromosome_column])
  for (m in chromosome_names) {
    if (startsWith(m, "chr")) {
      chr_prefix <- TRUE
    }
  }

  chrOrder <- c(as.character(paste0("chr", c(1:22))), "chrX", "chrY")
  if (chr_prefix == FALSE) {
    chrOrder <- c(as.character(c(1:22)), "X", "Y")
  }
  df_sorted[, chromosome_column] <- factor(df_sorted[, chromosome_column], levels = chrOrder)
  df_sorted <- df_sorted[order(df_sorted[, chromosome_column]), ]

  df_final <- NULL
  for (chrom in chrOrder) {
    if (nrow(df_sorted[!is.na(df_sorted[, chromosome_column]) & df_sorted[, chromosome_column] == chrom, ]) > 0) {
      chrom_regions <- df_sorted[df_sorted[, chromosome_column] == chrom, ]
      chrom_regions_sorted <- chrom_regions[with(chrom_regions, order(chrom_regions[, start_segment], chrom_regions[, end_segment])), ]
      df_final <- rbind(df_final, chrom_regions_sorted)
    }
  }
  return(df_final)
}


#' Function that performs stringr::str_replace on strings of multiple string columns of a dataframe
#'
#' @param df data frame
#' @param strings name of columns for which string replace is to be performed
#' @param pattern pattern to replace
#' @param replacement string to replace
#' @param replace_all logical - replace all occurrences
#' @return df
#'
#'
df_string_replace <- function(df, strings, pattern, replacement, replace_all = F) {
  stopifnot(is.data.frame(df))
  for (column_name in strings) {
    if (column_name %in% colnames(df)) {
      if (replace_all == F) {
        df[, column_name] <- stringr::str_replace(df[, column_name],
                                                  pattern = pattern,
                                                  replacement = replacement)
      }else{
        df[, column_name] <- stringr::str_replace_all(df[, column_name],
                                                      pattern = pattern,
                                                      replacement = replacement)
      }
    }
  }
  return(df)
}



#' Function that generate snv/indel + coding/noncoding stats for a given variant set
#'
#' @param calls data frame with variants in predisposition_genes
#' @param name type of variant group
#'
variant_stats_report <- function(calls, name = "variant_statistic") {

  call_stats <- list()
  call_stats[[name]] <- list()
  for (n in c("n", "n_snv", "n_indel", "n_coding", "n_noncoding")) {
    call_stats[[name]][[n]] <- 0
  }

  call_stats[[name]][["n"]] <- calls %>% nrow()
  call_stats[[name]][["n_snv"]] <- calls %>% dplyr::filter(VARIANT_CLASS == "SNV") %>% nrow()
  call_stats[[name]][["n_indel"]] <- calls %>% dplyr::filter(VARIANT_CLASS != "SNV") %>% nrow()
  call_stats[[name]][["n_coding"]] <- calls %>% dplyr::filter(CODING_STATUS == "coding") %>% nrow()
  call_stats[[name]][["n_noncoding"]] <- calls %>% dplyr::filter(CODING_STATUS != "coding") %>% nrow()

  return(call_stats)
}


#' Function that generates all contents of the cancer genome report (PCGR)
#'
#' @param project_directory name of project directory
#' @param query_vcf2tsv name of gzipped TSV file (vcf2tsv) with annotated query SNVs/InDels
#' @param pcgr_data List object with multiple PCGR data bundle annotations
#' @param pcgr_version PCGR software version
#' @param config Object with PCGR configuration parameters
#' @param sample_name sample identifier
#' @param cna_segments_tsv name of CNA segments file (tab-separated values)
#' @param cna_plot Path to PNG image with CNA plot
#' @param tier_model Variant tier model (pcgr_acmg/cpsr)
#'

generate_pcgr_report <- function(project_directory, query_vcf2tsv,
                            pcgr_data, pcgr_version = "dev",
                            config = NULL, sample_name = "SampleX",
                            cna_segments_tsv = NULL, cna_plot = NULL, tier_model = "pcgr_acmg") {

  rlogging::message(paste0("Initializing PCGR report - sample ", sample_name))

  ## Assert argument values
  invisible(assertthat::assert_that(!is.null(config) & is.list(config),
                                    msg = "Argument 'config' must be a non-NULL list object"))
  if(!is.null(cna_segments_tsv)){
    invisible(assertthat::assert_that(file.exists(cna_segments_tsv),
                                      msg = paste0("Filename provided for argument 'cna_segments_tsv' (",
                                                   cna_segments_tsv,") does not exist")))
    invisible(assertthat::assert_that(file.size(cna_segments_tsv) > 0,
                                      msg = paste0("File provided for argument 'cna_segments_tsv' (",
                                                   cna_segments_tsv,") has a filesize of zero")))
  }
  if(!is.null(query_vcf2tsv)){
    invisible(assertthat::assert_that(file.exists(query_vcf2tsv),
                                      msg = paste0("Filename provided for argument 'query_vcf2tsv' (",
                                                   query_vcf2tsv,") does not exist")))
    invisible(assertthat::assert_that(file.size(query_vcf2tsv) > 0,
                                      msg = paste0("File provided for argument 'query_vcf2tsv' (",
                                                   query_vcf2tsv,") has a filesize of zero")))

  }
  pcg_report <- pcgrr::init_report(config, sample_name, class = NULL,
                                   pcgr_data = pcgr_data, pcgr_version = pcgr_version)

  ## Set output and temporary file names
  sample_fname_pattern <-
    paste(sample_name, tier_model, pcgr_data[["assembly"]][["grch_name"]], sep = ".")

  fnames <- list()
  fnames[["maf"]] <-
    file.path(project_directory, paste0(sample_fname_pattern, ".maf"))
  fnames[["maf_tmp"]] <-
    file.path(project_directory, paste0(sample_fname_pattern, ".tmp.maf"))
  fnames[['vcf_mp']] <-
    file.path(project_directory, paste0(sample_fname_pattern, ".mp_input.vcf"))


  ## Retrieve relevant clinical trials for the tumor type in question

  if(config[['clinicaltrials']][['run']] == T){
    pcg_report_trials <- pcgrr::generate_report_data_trials(pcgr_data = pcgr_data,
                                                            sample_name = sample_name,
                                                            config = config)
    ## Update genome report with trial data
    pcg_report <-
      pcgrr::update_report(pcg_report, pcg_report_trials,
                           analysis_element = "clinicaltrials")
  }



  ## Load sample calls (check integrity, filter variants for depth/allelic fraction, append links etc)
  sample_calls <-
    pcgrr::get_calls(query_vcf2tsv, pcgr_data, sample_name,
                     config, oncotree = pcg_report[["metadata"]][["phenotype_ontology"]][["oncotree_query"]],
                     filenames = fnames)


  if (nrow(sample_calls) > 0) {

    ## Update MAF file with allelic support data
    #pcgrr::update_maf_allelic_support(sample_calls, fnames[["maf_tmp"]], fnames[["maf"]])

    ## Perform analyses in tumor-only mode
    if (pcg_report[['metadata']][['config']][['assay_properties']][["vcf_tumor_only"]] == TRUE) {
      pcg_report_tumor_only <- pcgrr::generate_report_data_tumor_only(sample_calls, sample_name, config)

      ## Generate data for SNVs/InDels
      ## -
      pcg_report_snv_indel_filtered <-
        pcgrr::generate_report_data_snv_indel(pcg_report_tumor_only[["variant_set"]][["filtered"]],
                                              pcgr_data, sample_name,
                                              config, callset = "germline-filtered callset",
                                              tier_model = tier_model)

      pcg_report_tumor_only[["upset_data"]] <-
        pcgrr::make_upset_plot_data(pcg_report_tumor_only$variant_set$unfiltered, config)
      num_upset_sources <- 0
      for (c in colnames(pcg_report_tumor_only[["upset_data"]])) {
        if (c != "VAR_ID") {
          if (sum(pcg_report_tumor_only[["upset_data"]][, c]) > 0) {
            num_upset_sources <- num_upset_sources + 1
          }
        }
      }
      if (num_upset_sources >= 2) {
        pcg_report_tumor_only[["upset_plot_valid"]] <- TRUE
      }

      ## Update genome report with SNV/InDels (display, tiers etc)
      pcg_report <-
        pcgrr::update_report(pcg_report, pcg_report_snv_indel_filtered,
                            analysis_element = "snv_indel")
      pcg_report <-
        pcgrr::update_report(pcg_report, pcg_report_tumor_only,
                            analysis_element = "tumor_only")

      ## Generate data for rainfall plot (SNVs)
      pcg_report_rainfall <-
        pcgrr::generate_report_data_rainfall(pcg_report$content$snv_indel$variant_set$tsv,
                                             build = pcg_report$metadata$genome_assembly)
      ## Update genome report
      pcg_report <-
        pcgrr::update_report(pcg_report, pcg_report_rainfall, analysis_element = "rainfall")

    }else{
      ## Generate report data for SNVs/InDels
      pcg_report_snv_indel <-
        pcgrr::generate_report_data_snv_indel(sample_calls, pcgr_data, sample_name,
                                              config, tier_model = tier_model)

      ## Update genome report
      pcg_report <- pcgrr::update_report(pcg_report, pcg_report_snv_indel,
                                             analysis_element = "snv_indel")
    }

    ## Estimate contribution of mutational signatures
    if (pcg_report[['metadata']][['config']][["msigs"]][["run"]] == T) {

      pcgrr::write_processed_vcf(calls = pcg_report$content$snv_indel$variant_set$tsv,
                                 sample_name = sample_name,
                                 output_directory = project_directory,
                                 vcf_fname = fnames[['vcf_mp']])

      pcg_report_signatures <-
        pcgrr::generate_report_data_signatures_mp(vcf_fname = paste0(fnames[['vcf_mp']],".gz"), pcgr_data,
                                                  sample_name, config,
                                                  type_specific =
                                                    !config[['msigs']][['all_reference_signatures']])

      ## Update genome report with signature info
      pcg_report <- pcgrr::update_report(pcg_report, pcg_report_signatures,
                                             analysis_element = "m_signature_mp")

      ## Generate report data for rainfall plot
      pcg_report_rainfall <-
        pcgrr::generate_report_data_rainfall(variant_set = pcg_report$content$snv_indel$variant_set$tsv,
                                             build = pcg_report$metadata$genome_assembly)

      ## Update genome report
      pcg_report <-
        pcgrr::update_report(pcg_report, pcg_report_rainfall, analysis_element = "rainfall")

      ## Generate report data for kataegis events (for WES/WGS runs)
      if(stringr::str_detect(pcg_report[['metadata']][['config']][['assay_properties']][['type']],"WGS|WES")){
        pcg_report_kataegis <-
          pcgrr::generate_report_data_kataegis(variant_set = pcg_report$content$snv_indel$variant_set$tsv,
                                               sample_name = sample_name,
                                               build = pcg_report$metadata$genome_assembly)
        ## Update genome report
        pcg_report <- pcgrr::update_report(pcg_report, pcg_report_kataegis, analysis_element = "kataegis")
      }
    }

    ## If assay is Tumor-Control and WES/WGS - perform MSI prediction
    if (pcg_report[['metadata']][['config']][['msi']][['run']] == T &
        stringr::str_detect(pcg_report[['metadata']][['config']][['assay_properties']][['type']],"WGS|WES") &
        pcg_report[['metadata']][['config']][['assay_properties']][["vcf_tumor_only"]] == FALSE) {
      pcg_report_msi <-
        pcgrr::generate_report_data_msi(pcg_report[["content"]][["snv_indel"]][["variant_set"]][["tsv"]],
                                        pcgr_data, sample_name, config)

      ## Update genome report with MSI info
      pcg_report <- pcgrr::update_report(pcg_report, pcg_report_msi,
                                             analysis_element = "msi")
    }

    ## Generate report contents for analysis of mutational burden (TMB)
    if (pcg_report[['metadata']][["config"]][["tmb"]][["run"]] == T) {
      pcg_report_tmb <-
        pcgrr::generate_report_data_tmb(pcg_report[["content"]][["snv_indel"]][["variant_set"]][["tsv"]],
                                        pcgr_data, sample_name, config)

      ## Update genome report with TMB info
      pcg_report <- pcgrr::update_report(pcg_report, pcg_report_tmb, analysis_element = "tmb")
    }
  }else{
    pcg_report[["content"]][["snv_indel"]][["zero"]] <- TRUE
    pcg_report[["metadata"]][["config"]][["other"]][["list_noncoding"]] <- FALSE
  }

  if (!is.null(cna_segments_tsv)) {
    pcg_report_cna <-
      pcgrr::generate_report_data_cna(cna_segments_tsv, pcgr_data, sample_name,
                                      config, transcript_overlap_pct = config[["cna"]][["cna_overlap_pct"]])
    pcg_report <-
      pcgrr::update_report(pcg_report, pcg_report_cna, analysis_element = "cna")

  }

  pcg_report_value_box <- pcgrr::generate_report_data_value_box(pcg_report, pcgr_data, sample_name, config)
  pcg_report <- pcgrr::update_report(pcg_report, pcg_report_value_box, analysis_element = "value_box")

  for (elem in c("tier1", "tier2", "tier3", "tier4")) {
    stat <- paste0("n_", elem)
    pcg_report[["content"]][["snv_indel"]][["variant_statistic"]][[stat]] <-
      nrow(pcg_report[["content"]][["snv_indel"]][["variant_set"]][[elem]])
    pcg_report[["content"]][["snv_indel"]][["variant_set"]][[elem]] <- NULL
  }
  pcg_report[["content"]][["snv_indel"]][["variant_set"]][["noncoding"]] <- NULL
  pcg_report[["content"]][["snv_indel"]][["variant_set"]][["coding"]] <- NULL
  pcg_report[["content"]][["snv_indel"]][["variant_set"]][["all"]] <- NULL
  pcg_report[["content"]][["cna"]][["variant_set"]][["cna_print"]] <- NULL
  pcg_report[["metadata"]][["phenotype_ontology"]] <- list()

  if (!is.null(cna_plot) && cna_plot != "None") {
    pcg_report[["content"]][["cna_plot"]][["png"]] <- cna_plot
    pcg_report[["content"]][["cna_plot"]][["eval"]] <- TRUE
  }
  return(pcg_report)
}

#' Function that writes contents of PCGR object to various output formats (Rmarkdown/flexdashboard HTML reports, JSON, tab-separated etc)
#'
#' @param project_directory working directory
#' @param report List object with all report data (PCGR/CPSR)
#' @param sample_name sample name
#' @param genome_assembly genome assembly (grch37/grch38)
#' @param tier_model type of tier model
#' @param output_format contents/file format of output (html/json/tsv/cna_tsv etc)
#' @param flexdb logical indicating if HTML output should be dashboard

write_report_output <- function(project_directory, report, sample_name, genome_assembly,
                         tier_model = "pcgr_acmg", output_format = "html", flexdb = F) {

  sample_fname_pattern <-
    paste(sample_name, tier_model, genome_assembly, sep = ".")

  fnames <- list()
  fnames[["snv_tsv_unfiltered"]] <-
    file.path(project_directory, paste0(sample_fname_pattern, ".snvs_indels.tiers.unfiltered.tsv"))
  fnames[["msigs_tsv"]] <-
    file.path(project_directory, paste0(sample_fname_pattern, ".mutational_signatures.tsv"))
  fnames[["snv_tsv"]] <-
    file.path(project_directory, paste0(sample_fname_pattern, ".snvs_indels.tiers.tsv"))
  fnames[["cna_tsv"]] <-
    file.path(project_directory, paste0(sample_fname_pattern, ".cna_segments.tsv"))
  fnames[["maf"]] <-
    file.path(project_directory, paste0(sample_fname_pattern, ".maf"))
  fnames[["maf_tmp"]] <-
    file.path(project_directory, paste0(sample_fname_pattern, ".tmp.maf"))
  fnames[["json"]] <-
    file.path(project_directory, paste0(sample_fname_pattern, ".json"))
  fnames[['html']] <-
    file.path(project_directory, paste0(sample_fname_pattern, ".html"))
  if(flexdb == T){
    fnames[['html']] <-
      file.path(project_directory, paste0(sample_fname_pattern, ".flexdb.html"))
  }

  disclaimer <- "disclaimer.md"
  report_theme <- report[["metadata"]][["config"]][["visual"]][["report_theme"]]


  sequencing_assay <- 'Unknown'
  sequencing_design <- 'Tumor-Control'
  css_fname <- system.file("templates","pcgr_flexdb_tumor_control.css", package = "pcgrr")

  if (tier_model == "cpsr") {
    sequencing_design <- 'Germline'
    disclaimer <- "disclaimer_predisposition.md"
  }else{
    sequencing_assay <- report[['metadata']][['config']][['assay_properties']][['type']]
    if(report[['metadata']][['config']][['assay_properties']][['vcf_tumor_only']] == T){
      sequencing_design <- 'Tumor-Only'
      css_fname <- system.file("templates","pcgr_flexdb_tumor_only.css", package = "pcgrr")
    }
  }

  if (output_format == "html") {

    markdown_input <- system.file("templates", "pcgr_rmarkdown_report.Rmd", package = "pcgrr")
    if (flexdb == T) {
      rlogging::message("------")
      rlogging::message("Writing HTML file (.html) with report contents - flexdashboard")
      navbar_items <- list()
      navbar_items[[1]] <- list('title' = paste0('<b>',sample_name,'</b> | <i>',
                                                 report[['metadata']][['config']][['tumor_properties']][['tumor_type']],
                                                 '</i> | ', sequencing_design,' | ',sequencing_assay),
                                href = "", target = "_blank", align = "right")
      navbar_items[[2]] <- list('icon' = "fa-github",
                                href = "https://github.com/sigven/pcgr", target = "_blank",
                                align = "right")

      markdown_input <- system.file("templates", "pcgr_flexdb_report.Rmd", package = "pcgrr")
      rmarkdown::render(markdown_input,
                        output_format = flexdashboard::flex_dashboard(orientation = "rows",
                                                                      theme = 'cosmo',
                                                      #vertical_layout = "fill",
                                                      css = css_fname,
                                                      navbar = navbar_items),
                        output_file = fnames[["html"]],
                        output_dir = project_directory,
                        clean = T,
                        intermediates_dir = project_directory, quiet = T)
    }else{
      rlogging::message("------")
      rlogging::message("Writing HTML file (.html) with report contents - rmarkdown")
      css_fname <- NULL
      if (tier_model == "cpsr") {
        markdown_input <- system.file("templates", "cpsr_rmarkdown_report.Rmd", package = "pcgrr")
        css_fname <- system.file("templates","cpsr.css", package = "pcgrr")
      }
      rmarkdown::render(markdown_input,
                        output_format =
                          rmarkdown::html_document(theme = report_theme,
                                                   toc = T, toc_depth = 3,
                                                   toc_float = T,
                                                   number_sections = F,
                                                   css = css_fname,
                                                   includes = rmarkdown::includes(after_body = disclaimer)),
                        output_file = fnames[["html"]],
                        output_dir = project_directory,
                        clean = T,
                        intermediates_dir = project_directory, quiet = T)
    }
  }
  if(output_format == "json"){
    if (!is.null(report[["cna_plot"]][["png"]])) {
      report[["cna_plot"]][["png"]] <- NULL
    }
    if (!is.null(report[["tmb"]][["tcga_tmb"]])) {
      report[["tmb"]][["tcga_tmb"]] <- NULL
    }
    rlogging::message("------")
    rlogging::message("Writing JSON file (.json) with report contents")
    pcgr_json <- jsonlite::toJSON(report, pretty = T, na = "string", null = "null", force = T)
    write(pcgr_json, fnames[["json"]])
    gzip_command <- paste0("gzip -f ", fnames[["json"]])
    system(gzip_command, intern = F)
  }

  if(output_format == "snv_tsv" | output_format == "snv_tsv_unfiltered"){
    output_format_slim = stringr::str_replace(output_format,"snv_","")
    if (NROW(report[["content"]][["snv_indel"]][["variant_set"]][[output_format_slim]]) > 0) {
      rlogging::message("------")
      if(tier_model == "pcgr_acmg"){
        rlogging::message(paste0("Writing SNV/InDel tab-separated output file with PCGR annotations - ('",
                               output_format_slim,"')"))
      }
      if(tier_model == "cpsr"){
        rlogging::message(paste0("Writing SNV/InDel tab-separated output file with CPSR annotations - ('",
                                 output_format_slim,"')"))
      }
      write.table(report[["content"]][["snv_indel"]][["variant_set"]][[output_format_slim]],
                  file = fnames[[output_format]], sep = "\t", col.names = T,
                  row.names = F, quote = F)
    }
  }

  if(output_format == "msigs_tsv"){
    if (NROW(report[["content"]][["m_signature_mp"]][["result"]][['tsv']]) > 0) {
      rlogging::message("------")
      rlogging::message("Writing tab-separated output file with details of contributing mutational signatures - ('tsv')")
      write.table(report[["content"]][["m_signature_mp"]][["result"]][['tsv']],
                  file = fnames[[output_format]], sep = "\t", col.names = T,
                  row.names = F, quote = F)
    }
  }

  if(output_format == "cna_tsv"){
    if (NROW(report[["content"]][["cna"]][["variant_set"]][['tsv']]) > 0) {
      rlogging::message("------")
      rlogging::message("Writing CNA tab-separated output file with PCGR annotations (.tsv.gz)")
      write.table(report[["content"]][["cna"]][["variant_set"]][['tsv']],
                  file = fnames[['cna_tsv']], sep = "\t", col.names = T,
                  row.names = F, quote = F)
      gzip_command <- paste0("gzip -f ", fnames[['cna_tsv']])
      system(gzip_command, intern = F)
    }
  }

}


#' Function that generates a data frame with basic biomarker annotations from tier1 variants
#'
#' @param tier1_variants df with tier 1 variants
#' @param sample_id Sample identifier
#'
#' @return tsv_variants data frame with all tier 1 biomarkers for tab-separated output
#'
generate_biomarker_tsv <- function(clinical_evidence_items, sample_name = "test", tier = "tier 1") {

  bm_tags <- c("CLINICAL_SIGNIFICANCE", "EVIDENCE_LEVEL", "EVIDENCE_TYPE", "EVIDENCE_DIRECTION", "CANCER_TYPE",
               "DESCRIPTION", "BIOMARKER_MAPPING", "THERAPEUTIC_CONTEXT", "RATING", "CITATION")
  all_biomarker_tags <- c(c("GENOMIC_CHANGE", "GENOME_VERSION", "VCF_SAMPLE_ID", "SYMBOL", "CONSEQUENCE"), bm_tags)
  tsv_biomarkers <- NULL
  if (nrow(clinical_evidence_items) > 0) {
    biomarker_tsv <- clinical_evidence_items
    biomarker_tsv$VCF_SAMPLE_ID <- sample_name
    biomarker_tsv <- biomarker_tsv %>% dplyr::select(dplyr::one_of(all_biomarker_tags))
    biomarker_tsv$TIER <- "TIER 1"
    biomarker_tsv$TIER_DESCRIPTION <- "Clinical biomarker - Predictive/prognostic/diagnostic"

    tmp2 <- as.data.frame(
      biomarker_tsv %>%
        dplyr::rowwise() %>%
        dplyr::mutate(CITATION2 = paste(unlist(stringr::str_replace_all(stringr::str_match_all(CITATION, ">.+<"), "^>|<$", "")), collapse = ";")))

    biomarker_tsv$CITATION <- tmp2$CITATION2
    tsv_biomarkers <- biomarker_tsv %>% dplyr::distinct()
  }
  return(tsv_biomarkers)
}

#' Function that generates dense and tiered annotated variant datasets
#' @param variant_set List with tiered variants
#' @param pcgr_data List of data frames with PCGR data annotations
#' @param sample_name Sample identifier
#'
#' @return tsv_variants data frame with tier-annotated list of variants for tab-separated output
#'
generate_tier_tsv <- function(variant_set, pcgr_data, config, sample_name = "test") {

  tags <- NULL
  if (!is.null(config[["custom_tags"]])) {
    if (config[["custom_tags"]][["custom_tags"]] != "") {
      tags <- stringr::str_split(config[["custom_tags"]][["custom_tags"]], pattern = ",")[[1]]
    }
  }
  rlogging::message("Generating tiered set of result variants for output in tab-separated values (TSV) file")
  tsv_variants <- NULL
  for (tier in c("tier1", "tier2", "tier3", "tier4", "noncoding")) {
    if (nrow(variant_set[[tier]]) > 0) {
      tierset <- variant_set[[tier]]
      tierset$VCF_SAMPLE_ID <- sample_name
      tsv_columns <- pcgr_data[["annotation_tags"]][["tsv"]]
      if (!is.null(tags)) {
        for (t in tags) {
          t <- stringr::str_trim(t)
          if (t %in% colnames(tierset)) {
            tsv_columns <- c(tsv_columns, t)
          }
        }
      }

      if (tier == "tier1") {
        tierset$TIER_DESCRIPTION <- "Variants of strong clinical significance"
        tierset$TIER <- "TIER 1"
      }
      if (tier == "tier2") {
        tierset$TIER_DESCRIPTION <- "Variants of potential clinical significance"
        tierset$TIER <- "TIER 2"
      }
      if (tier == "tier3") {
        tierset$TIER_DESCRIPTION <- "Variants of uncertain significance"
        tierset$TIER <- "TIER 3"
      }
      if (tier == "tier4") {
        tierset$TIER_DESCRIPTION <- "Other coding mutation"
        tierset$TIER <- "TIER 4"
      }
      if (tier == "noncoding") {
        tierset$TIER_DESCRIPTION <- "Noncoding mutation"
        tierset$TIER <- "NONCODING"
      }
      tierset <- tierset %>%
        dplyr::select(dplyr::one_of(tsv_columns)) %>%
        dplyr::distinct()

      tsv_variants <- dplyr::bind_rows(tsv_variants, tierset)
    }
  }
  tsv_variants$GENE_NAME <- unlist(lapply(stringr::str_match_all(tsv_variants$GENE_NAME, ">.+<"), paste, collapse = ","))
  tsv_variants$GENE_NAME <- stringr::str_replace_all(tsv_variants$GENE_NAME, ">|<", "")
  tsv_variants$CLINVAR <- unlist(lapply(stringr::str_match_all(tsv_variants$CLINVAR, ">.+<"), paste, collapse = ","))
  tsv_variants$CLINVAR <- stringr::str_replace_all(tsv_variants$CLINVAR, ">|<", "")
  tsv_variants$PROTEIN_DOMAIN <- unlist(lapply(stringr::str_match_all(tsv_variants$PROTEIN_DOMAIN, ">.+<"), paste, collapse = ","))
  tsv_variants$PROTEIN_DOMAIN <- stringr::str_replace_all(tsv_variants$PROTEIN_DOMAIN, ">|<", "")
  tsv_variants$TCGA_FREQUENCY <- stringr::str_replace_all(tsv_variants$TCGA_FREQUENCY, "<a href='https://portal.gdc.cancer.gov/projects/TCGA-[A-Z]{1,}' target=\"_blank\">|</a>", "")
  tsv_variants <- tsv_variants %>% dplyr::distinct()

  return(tsv_variants)
}


#' Function that generates tiered variant sets for SNVs/InDels
#'
#' @param sample_calls variant calls subject to mutational signature analysis
#' @param pcgr_data object with PCGR annotation data
#' @param sample_name sample identifier
#' @param pcgr_config Object with PCGR configuration parameters
#' @param biomarker_mapping_stringency quality level for biomarkers
#' @param callset type of calls
#' @param tier_model tier model (pcgr_acmg)
#'
#' @return pcg_report_data data frame with all report elements
#'
generate_report_data_snv_indel <- function(sample_calls, pcgr_data, sample_name, config, callset = "somatic calls",
                                           biomarker_mapping_stringency = 1, tier_model = "pcgr_acmg") {

  rlogging::message("------")
  rlogging::message(paste0("Generating data for tiered cancer genome report - ",
                           callset, " tier model '", tier_model, "'"))

  pcg_report_snv_indel <- pcgrr::init_report(config, sample_name, class = "snv_indel")
  pcg_report_snv_indel[["eval"]] <- TRUE
  pcg_report_snv_indel[["variant_set"]][["all"]] <- sample_calls

  ## Get basic variant statistics (type, coding status)
  call_stats <- pcgrr::variant_stats_report(sample_calls, name = "variant_statistic")
  for(stat in c('n','n_snv','n_indel','n_coding','n_noncoding')){
    pcg_report_snv_indel[['variant_statistic']][[stat]] <- call_stats[['variant_statistic']][[stat]]
  }
  rlogging::message(paste0("Number of protein-coding variants: ",
                           pcg_report_snv_indel[["variant_statistic"]][["n_coding"]]))

  ## Add custom annotation tags to lists of tags to display
  if (!is.null(config[["custom_tags"]])) {
    if (config[["custom_tags"]][["custom_tags"]] != "") {
      tags <- stringr::str_split(config[["custom_tags"]][["custom_tags"]],
                                 pattern = ",")[[1]]
      for (t in tags) {
        t <- stringr::str_trim(t)
        if (t %in% colnames(sample_calls)) {
          pcgr_data[["annotation_tags"]][["all"]] <- c(pcgr_data[["annotation_tags"]][["all"]], t)
        }
      }
    }
  }

  if (pcg_report_snv_indel[["variant_statistic"]][["n"]] > 0) {

    tumor_type <- config[['tumor_properties']][['tumor_type']]

    ## Get all clinical evidence items that overlap query set (NOT tumor-type specific)
    biomarker_hits_snv_indels_any <-
      pcgrr::get_clin_assocs_snv_indel(pcg_report_snv_indel[["variant_set"]][["all"]], pcgr_data,
                                       tumor_type, tumor_type_specificity = "any",
                                       mapping_stringency = biomarker_mapping_stringency)

    ## Assign putative TIER 2 variant set
    pcg_report_snv_indel[["clin_eitem"]][["any_ttype"]] <-
      biomarker_hits_snv_indels_any$clin_eitem
    pcg_report_snv_indel[["variant_set"]][["tier2"]] <-
      biomarker_hits_snv_indels_any$variant_set

    ## Get all clinical evidence items that overlap query set (if tumor type is specified)
    if(tumor_type != "Cancer, NOS"){
      biomarker_hits_snv_indels_specific <-
        pcgrr::get_clin_assocs_snv_indel(pcg_report_snv_indel[["variant_set"]][["all"]], pcgr_data,
                                         tumor_type, tumor_type_specificity = "specific",
                                         mapping_stringency = biomarker_mapping_stringency)

      ## Assign putative TIER 1 variant set
      pcg_report_snv_indel[["clin_eitem"]][["specific_ttype"]] <-
        biomarker_hits_snv_indels_specific$clin_eitem
      pcg_report_snv_indel[["variant_set"]][["tier1"]] <-
        biomarker_hits_snv_indels_specific$variant_set
    }

    ## Remove potential overlap/redundancies and assign final TIER1/TIER2 classification
    pcg_report_snv_indel <- pcgrr::assign_tier1_tier2_acmg(pcg_report_snv_indel)
    tier12 <- rbind(pcg_report_snv_indel[["variant_display"]][["tier1"]],
                    pcg_report_snv_indel[["variant_display"]][["tier2"]])

    ## Determine TIER 3 variant set: coding mutations in oncogenes/tumor suppressors/cancer census genes
    pcg_report_snv_indel[["variant_set"]][["tier3"]] <-
      dplyr::select(pcg_report_snv_indel[["variant_set"]][["all"]],
                    dplyr::one_of(pcgr_data[["annotation_tags"]][["all"]])) %>%
      dplyr::filter(CODING_STATUS == "coding") %>%
      dplyr::filter(ONCOGENE == TRUE | TUMOR_SUPPRESSOR == TRUE)
    if (nrow(tier12) > 0 & nrow(pcg_report_snv_indel[["variant_set"]][["tier3"]]) > 0) {
      pcg_report_snv_indel[["variant_set"]][["tier3"]] <-
        dplyr::anti_join(pcg_report_snv_indel[["variant_set"]][["tier3"]],
                         tier12, by = c("GENOMIC_CHANGE"))
    }
    tier123 <- tier12
    if (nrow(pcg_report_snv_indel[["variant_set"]][["tier3"]]) > 0) {
      pcg_report_snv_indel[["variant_set"]][["tier3"]] <-
        pcg_report_snv_indel[["variant_set"]][["tier3"]] %>%
        dplyr::arrange(desc(OPENTARGETS_RANK), desc(ONCOSCORE))
      tier123 <- tier12 %>%
        dplyr::bind_rows(dplyr::select(pcg_report_snv_indel[["variant_set"]][["tier3"]],
                                       GENOMIC_CHANGE)) %>%
        dplyr::distinct()
      pcg_report_snv_indel[["variant_display"]][["tier3"]][["proto_oncogene"]] <-
        dplyr::select(pcg_report_snv_indel[["variant_set"]][["tier3"]],
                      dplyr::one_of(pcgr_data[["annotation_tags"]][["tier3_display"]])) %>%
        dplyr::filter(ONCOGENE == TRUE & (is.na(TUMOR_SUPPRESSOR) | TUMOR_SUPPRESSOR == FALSE))
      pcg_report_snv_indel[["variant_display"]][["tier3"]][["tumor_suppressor"]] <-
        dplyr::select(pcg_report_snv_indel[["variant_set"]][["tier3"]],
                      dplyr::one_of(pcgr_data[["annotation_tags"]][["tier3_display"]])) %>%
        dplyr::filter(TUMOR_SUPPRESSOR == TRUE)
    }

    ## Determine TIER 4: Other coding mutations
    pcg_report_snv_indel[["variant_set"]][["tier4"]] <-
      dplyr::select(pcg_report_snv_indel[["variant_set"]][["all"]],
                    dplyr::one_of(pcgr_data[["annotation_tags"]][["all"]])) %>%
      dplyr::filter(CODING_STATUS == "coding")
    if (nrow(tier123) > 0 & nrow(pcg_report_snv_indel[["variant_set"]][["tier4"]]) > 0) {
      pcg_report_snv_indel[["variant_set"]][["tier4"]] <-
        dplyr::anti_join(pcg_report_snv_indel[["variant_set"]][["tier4"]],
                         tier123, by = c("GENOMIC_CHANGE"))
    }
    if (nrow(pcg_report_snv_indel[["variant_set"]][["tier4"]]) > 0) {
      pcg_report_snv_indel[["variant_set"]][["tier4"]] <-
        pcg_report_snv_indel[["variant_set"]][["tier4"]] %>%
        dplyr::arrange(desc(OPENTARGETS_RANK), desc(ONCOSCORE))
      pcg_report_snv_indel[["variant_display"]][["tier4"]] <-
        dplyr::select(pcg_report_snv_indel[["variant_set"]][["tier4"]],
                      dplyr::one_of(pcgr_data[["annotation_tags"]][["tier4_display"]]))
    }

    ## Determine non-coding mutation set
    pcg_report_snv_indel[["variant_set"]][["noncoding"]] <-
      dplyr::select(pcg_report_snv_indel[["variant_set"]][["all"]],
                    dplyr::one_of(pcgr_data[["annotation_tags"]][["all"]])) %>%
      dplyr::filter(CODING_STATUS == "noncoding")
    if (nrow(pcg_report_snv_indel[["variant_set"]][["noncoding"]]) > 0) {
      if (nrow(tier123) > 0) {
        pcg_report_snv_indel[["variant_set"]][["noncoding"]] <-
          dplyr::anti_join(pcg_report_snv_indel[["variant_set"]][["noncoding"]], tier123,
                           by = c("GENOMIC_CHANGE"))
      }
      pcg_report_snv_indel[["variant_set"]][["noncoding"]] <-
        pcg_report_snv_indel[["variant_set"]][["noncoding"]] %>%
        dplyr::arrange(desc(OPENTARGETS_RANK), desc(ONCOSCORE))
      pcg_report_snv_indel[["variant_display"]][["noncoding"]] <-
        dplyr::select(pcg_report_snv_indel[["variant_set"]][["noncoding"]],
                      dplyr::one_of(pcgr_data[["annotation_tags"]][["tier5_display"]]))
    }

    ## Make TSV content with variant set
    pcg_report_snv_indel[["variant_statistic"]][["n_noncoding"]] <-
      pcg_report_snv_indel[["variant_set"]][["noncoding"]] %>% nrow()
    pcg_report_snv_indel[["variant_set"]][["tsv"]] <-
      pcgrr::generate_tier_tsv(pcg_report_snv_indel[["variant_set"]],
                               pcgr_data = pcgr_data, config, sample_name = sample_name)

  }

  return(pcg_report_snv_indel)

}


append_ucsc_segment_link <- function(var_df, hgname = "hg38", chrom = NULL, start = NULL, end = NULL) {
  ucsc_browser_prefix <- paste0("http://genome.ucsc.edu/cgi-bin/hgTracks?db=", hgname, "&position=")
  if (!is.null(chrom) & !is.null(start) & !is.null(end) &
      chrom %in% colnames(var_df) & start %in% colnames(var_df) &
      end %in% colnames(var_df)) {
    var_df <- var_df %>%
      dplyr::mutate(SEGMENT_LINK = paste0("<a href='", paste0(ucsc_browser_prefix,
                                                              paste0(!!rlang::sym(chrom),
                                                                     ":", !!rlang::sym(start),
                                                                     "-", !!rlang::sym(end)),
                                                              "' target=\"_blank\">",
                                                              paste0(!!rlang::sym(chrom), ":",
                                                                     !!rlang::sym(start), "-",
                                                                     !!rlang::sym(end)), "</a>")))
  }else{
    var_df$SEGMENT_LINK <- NA
  }
  return(var_df)


}

#' Function that adds HTML links to different genetic variant identifiers
#'
#' @param var_df data frame with variants
#' @param linktype type of link
#' @param pcgr_data PCGR data structure
#'
append_tcga_var_link <- function(var_df, pcgr_data = NULL,
                                 linktype = "dbsource"){

  rlogging::message("Adding annotation links - TCGA")

  if (any(grepl(paste0("^TCGA_FREQUENCY$"), names(var_df))) &
      any(grepl(paste0("^VAR_ID$"), names(var_df))) &
      !is.null(pcgr_data)) {
    var_df_unique_slim <- dplyr::select(var_df, VAR_ID, TCGA_FREQUENCY) %>%
      dplyr::filter(!is.na(TCGA_FREQUENCY)) %>%
      dplyr::distinct()
    if (nrow(var_df_unique_slim) > 0) {
      var_df_unique_slim_melted <- var_df_unique_slim %>%
        tidyr::separate_rows(TCGA_FREQUENCY, sep = ",") %>%
        tidyr::separate(TCGA_FREQUENCY, c("tumor", "percentage", "affected", "cohort_size"),
                        sep = "\\|", convert = T) %>%
        dplyr::left_join(pcgr_data[["tcga"]][["projects"]], by = "tumor") %>%
        dplyr::arrange(VAR_ID, desc(percentage))
      if (linktype == "dbsource") {
        var_df_unique_slim_melted <- var_df_unique_slim_melted %>%
          dplyr::mutate(tmp_assoc = paste0("<a href='https://portal.gdc.cancer.gov/projects/TCGA-",
                                           tumor, "' target=\"_blank\">", name, "</a>: ",
                                           percentage, "% (", affected, "/", cohort_size, ")"))
      }

      var_df_links <- dplyr::group_by(var_df_unique_slim_melted, VAR_ID) %>%
        dplyr::summarise(TCGALINK = unlist(paste(tmp_assoc, collapse = ", ")), .groups = "drop") %>%
        dplyr::select(VAR_ID, TCGALINK) %>%
        magrittr::set_colnames(c('VAR_ID','TCGA_FREQUENCY'))
      var_df <- dplyr::rename(var_df, TCGA_FREQUENCY_RAW = TCGA_FREQUENCY)
      var_df <- dplyr::left_join(var_df, var_df_links, by = c("VAR_ID" = "VAR_ID"))
    }
  }
  return(var_df)
}

append_dbnsfp_var_link <- function(var_df, linktype = "dbsource") {

  rlogging::message("Adding annotation links - dbNSFP")

  #assertable::assert_colnames(var_df, colnames = c("EFFECT_PREDICTIONS"), only_colnames = F, quiet = T)
  if (any(grepl(paste0("EFFECT_PREDICTIONS"), names(var_df)))) {
    var_df <- var_df %>%
      dplyr::mutate(PREDICTED_EFFECT = EFFECT_PREDICTIONS)
    i <- 1
    while (i <= nrow(pcgrr::effect_prediction_algos)) {
      str_to_replace <- paste0(pcgrr::effect_prediction_algos[i, "algorithm"], ":")
      replacement_str <-
        paste0("<a href='", pcgrr::effect_prediction_algos[i, "url"], "' target='_blank'>",
               pcgrr::effect_prediction_algos[i, "display_name"], "</a>:")
      algorithm_display <- paste0(pcgrr::effect_prediction_algos[i, "display_name"], ":")
      var_df <- var_df %>%
        dplyr::mutate(PREDICTED_EFFECT = stringr::str_replace(PREDICTED_EFFECT, str_to_replace, replacement_str))
      i <- i + 1
    }
  }
  else{
    var_df$PREDICTED_EFFECT <- NA
  }
  return(var_df)

}

#' Function that adds HTML links to different genetic variant identifiers
#'
#' @param var_df data frame with variants
#' @param linktype type of link
#' @param pcgr_data PCGR data structure
#'
append_drug_var_link <- function(var_df, pcgr_data = NULL,
                             linktype = "dbsource"){

  rlogging::message("Adding annotation links - targeted cancer drugs")
  if (any(grepl(paste0("^CHEMBL_COMPOUND_ID$"), names(var_df))) &
      any(grepl(paste0("^SYMBOL$"), names(var_df))) &
      any(grepl(paste0("^VAR_ID$"), names(var_df))) &
      !is.null(pcgr_data)) {
    var_df_unique_slim <- dplyr::select(var_df, VAR_ID, SYMBOL, CHEMBL_COMPOUND_ID) %>%
      dplyr::filter(!is.na(CHEMBL_COMPOUND_ID)) %>%
      dplyr::distinct()
    if (nrow(var_df_unique_slim) > 0) {
      var_df_unique_slim_melted <- var_df_unique_slim %>%
        tidyr::separate_rows(CHEMBL_COMPOUND_ID, sep = "&")
      chembl_drugs <-
        dplyr::select(pcgr_data[["antineopharma"]][["antineopharma"]],
                      molecule_chembl_id, symbol, nci_concept_display_name) %>%
        dplyr::arrange(symbol) %>%
        dplyr::distinct()
      var_df_unique_slim_melted <- var_df_unique_slim_melted %>%
        dplyr::left_join(chembl_drugs, by = c("CHEMBL_COMPOUND_ID" = "molecule_chembl_id",
                                              "SYMBOL" = "symbol")) %>%
        dplyr::filter(!is.na(nci_concept_display_name)) %>%
        dplyr::distinct()
      if (nrow(var_df_unique_slim_melted) > 0) {
        if (linktype == "dbsource") {
          var_df_unique_slim_melted <- var_df_unique_slim_melted %>%
            dplyr::mutate(tmp_antineopharma =
                            paste0("<a href='https://www.targetvalidation.org/summary?drug=",
                                   CHEMBL_COMPOUND_ID, "' target=\"_blank\">",
                                   nci_concept_display_name, "</a>"))
        }
        var_df_unique_slim_melted_terms <- dplyr::select(var_df_unique_slim_melted, VAR_ID, nci_concept_display_name)
        var_df_terms <- dplyr::group_by(var_df_unique_slim_melted_terms, VAR_ID) %>%
          dplyr::summarise(CHEMBL_COMPOUND_TERMS = paste(nci_concept_display_name, collapse = ", "))
        var_df_links <- dplyr::group_by(var_df_unique_slim_melted, VAR_ID) %>%
          dplyr::summarise(ANTINEOPHARMALINK = unlist(paste(tmp_antineopharma, collapse = ", "))) %>%
          dplyr::select(VAR_ID, ANTINEOPHARMALINK) %>%
          dplyr::distinct()
        var_df <- dplyr::left_join(var_df, var_df_links, by = c("VAR_ID" = "VAR_ID"))
        var_df <- dplyr::left_join(var_df, var_df_terms, by = c("VAR_ID" = "VAR_ID"))
      }else{
        var_df$ANTINEOPHARMALINK <- NA
        var_df$CHEMBL_COMPOUND_TERMS <- NA
      }
    }
    else{
      var_df$ANTINEOPHARMALINK <- NA
      var_df$CHEMBL_COMPOUND_TERMS <- NA
    }
  }
  else{
    cat("WARNING: Could not generate links with targeted compounds - no ANTINEOPHARMA info provided in annotated VCF", sep = "\n")
    var_df$ANTINEOPHARMALINK <- NA
    var_df$CHEMBL_COMPOUND_TERMS <- NA
  }
  return(var_df)
}

#' Function that adds HTML links to different genetic variant identifiers
#'
#' @param var_df data frame with variants
#' @param linktype type of link
#' @param pcgr_data PCGR data structure
#' @param oncotree Oncotree data frame
#'
append_otargets_pheno_link <- function(var_df, pcgr_data = NULL, oncotree = NULL,
                                       linktype = "dbsource"){

  rlogging::message("Adding annotation links - gene-cancer type associations (OpenTargets Platform)")

  if (any(grepl(paste0("^OPENTARGETS_DISEASE_ASSOCS$"), names(var_df))) &
      any(grepl(paste0("^ENSEMBL_GENE_ID$"), names(var_df))) &
      !is.null(pcgr_data) & nrow(oncotree) > 0) {
    var_df_unique_slim <- dplyr::select(var_df, ENSEMBL_GENE_ID,
                                        OPENTARGETS_DISEASE_ASSOCS) %>%
      dplyr::filter(!is.na(OPENTARGETS_DISEASE_ASSOCS)) %>%
      dplyr::distinct()
    associations_found <- 0
    if (nrow(var_df_unique_slim) > 0) {
      var_df_unique_slim_melted <- as.data.frame(
        var_df_unique_slim %>%
          tidyr::separate_rows(OPENTARGETS_DISEASE_ASSOCS, sep = "&") %>%
          tidyr::separate(OPENTARGETS_DISEASE_ASSOCS, c("cui", "efo_id", "ot_is_direct", "ot_score"), sep = ":", remove = T) %>%
          dplyr::inner_join(dplyr::select(oncotree, cui), by = c("cui" = "cui")) %>%
          dplyr::left_join(pcgr_data[["phenotype_ontology"]][["umls"]], by = c("cui" = "cui")) %>%
          dplyr::mutate(ot_score = as.numeric(ot_score))
      )

      if (nrow(var_df_unique_slim_melted) > 0) {
        associations_found <- 1
        var_df_unique_slim_melted <- as.data.frame(
          var_df_unique_slim_melted %>%
            dplyr::group_by(ENSEMBL_GENE_ID, efo_id, cui_name) %>%
            dplyr::summarise(score = max(ot_score, na.rm = T)) %>%
            dplyr::distinct() %>%
            dplyr::arrange(desc(score))
        )

        if (linktype == "dbsource") {
          var_df_unique_slim_melted <- var_df_unique_slim_melted %>%
            dplyr::mutate(tmp_assoc = paste0("<a href='https://www.targetvalidation.org/evidence/",
                                             ENSEMBL_GENE_ID, "/", efo_id, "' target=\"_blank\">",
                                             cui_name, "</a>"))
        }

        var_df_unique_slim_melted_terms <- dplyr::select(var_df_unique_slim_melted, ENSEMBL_GENE_ID, cui_name)
        var_df_terms <- dplyr::group_by(var_df_unique_slim_melted_terms, ENSEMBL_GENE_ID) %>%
          dplyr::summarise(OT_DISEASE_TERMS = paste(cui_name, collapse = "&"))
        var_df_links <- dplyr::group_by(var_df_unique_slim_melted, ENSEMBL_GENE_ID) %>%
          dplyr::summarise(OT_DISEASE_LINK = unlist(paste(tmp_assoc, collapse = ", ")),
                           OPENTARGETS_RANK = max(score))
        var_df_links <- dplyr::select(var_df_links, ENSEMBL_GENE_ID, OT_DISEASE_LINK, OPENTARGETS_RANK)
        var_df <- dplyr::left_join(var_df, var_df_links, by = c("ENSEMBL_GENE_ID" = "ENSEMBL_GENE_ID"))
        var_df <- dplyr::left_join(var_df, var_df_terms, by = c("ENSEMBL_GENE_ID" = "ENSEMBL_GENE_ID"))
        var_df <- var_df %>% dplyr::mutate(OPENTARGETS_RANK = dplyr::if_else(is.na(OPENTARGETS_RANK),
                                                                             as.numeric(0), OPENTARGETS_RANK))

      }else{
        var_df$OT_DISEASE_LINK <- NA
        var_df$OT_DISEASE_TERMS <- NA
        var_df$OPENTARGETS_RANK <- 0
      }
    }else{
      if (associations_found == 0) {
        var_df$OT_DISEASE_LINK <- NA
        var_df$OT_DISEASE_TERMS <- NA
        var_df$OPENTARGETS_RANK <- 0
      }
    }
  }else{
    rlogging::warning("Could not generate Open Targets association links - no Open Targets annotations provided in annotated VCF")
    var_df$OT_DISEASE_LINK <- NA
    var_df$OT_DISEASE_TERMS <- NA
    var_df$OPENTARGETS_RANK <- 0
  }
  return(var_df)
}

append_disgenet_var_link <- function(var_df, linktype = "dbsource",
                                      pcgr_data = NULL) {

  rlogging::message("Adding annotation links - gene-cancer type associations (DisGenet)")

  if (any(grepl(paste0("^DISGENET_CUI$"), names(var_df))) &
      any(grepl(paste0("^SYMBOL$"), names(var_df))) &
      !is.null(pcgr_data)) {
    var_df_unique_slim <- dplyr::select(var_df, SYMBOL, DISGENET_CUI) %>%
      dplyr::filter(!is.na(DISGENET_CUI)) %>%
      dplyr::distinct()
    if (nrow(var_df_unique_slim) > 0) {
      var_df_unique_slim_melted <- var_df_unique_slim %>%
        tidyr::separate_rows(DISGENET_CUI, sep = "&") %>%
        dplyr::left_join(pcgr_data[["phenotype_ontology"]][["umls"]], by = c("DISGENET_CUI" = "cui"))
      if (linktype == "dbsource") {
        var_df_unique_slim_melted <- var_df_unique_slim_melted %>%
          dplyr::mutate(tmp_assoc = paste0("<a href='https://www.ncbi.nlm.nih.gov/medgen/",
                                           DISGENET_CUI, "' target=\"_blank\">", cui_name, "</a>"))
      }

      var_df_unique_slim_melted_terms <- dplyr::select(var_df_unique_slim_melted, SYMBOL, cui_name)
      var_df_terms <- dplyr::group_by(var_df_unique_slim_melted_terms, SYMBOL) %>%
        dplyr::summarise(DISGENET_TERMS = paste(cui_name, collapse = "&"))
      var_df_links <- dplyr::group_by(var_df_unique_slim_melted, SYMBOL) %>%
        dplyr::summarise(DISGENET_LINK = unlist(paste(tmp_assoc, collapse = ", ")))
      var_df_links <- dplyr::select(var_df_links, SYMBOL, DISGENET_LINK)
      var_df <- dplyr::left_join(var_df, var_df_links, by = c("SYMBOL" = "SYMBOL"))
      var_df <- dplyr::left_join(var_df, var_df_terms, by = c("SYMBOL" = "SYMBOL"))
    }
    else{
      var_df$DISGENET_LINK <- NA
      var_df$DISGENET_TERMS <- NA
    }
  }
  else{
    cat("WARNING: Could not generate cancer gene association links - no Disgenet cancer associations provided in annotated VCF", sep = "\n")
    var_df$DISGENET_LINK <- NA
    var_df$DISGENET_TERMS <- NA
  }
  return(var_df)

}

#' Function that adds read support (depth, allelic fraction) for tumor and normal and filters according to settings
#'
#' @param vcf_df data frame with variants
#' @param config list with workflow configuration values
#' @param precision number of significant digits for allelic fraction estimation
#'
#' @return vcf_df
#'
append_read_support <- function(vcf_df, config = NULL, precision = 3) {

  invisible(assertthat::assert_that(!is.null(vcf_df), msg = "Argument 'vcf_df' cannot not be NULL"))
  invisible(assertthat::assert_that(is.data.frame(vcf_df), msg = paste0("Argument 'vcf_df' must be of type 'data.frame'")))
  invisible(assertthat::assert_that(!is.null(config), msg = "Argument 'config' cannot not be NULL"))
  invisible(assertthat::assert_that(class(config)[2] == "list",
                                    msg = paste0("Argument 'config' must be of type list, not ", class(config)[2])))
  if (is.null(config$allelic_support))return(vcf_df)

  for (v in c("DP_TUMOR", "AF_TUMOR", "DP_CONTROL", "AF_CONTROL", "CALL_CONFIDENCE")) {
    vcf_df[v] <- NA
  }
  found_tumor_tag <- 0
  for (tag_name in names(config$allelic_support)) {
    if (config$allelic_support[[tag_name]] != "" &
        tag_name != "tumor_dp_min" &
        tag_name != "tumor_af_min" &
        tag_name != "control_dp_min" &
        tag_name != "control_af_max") {
      config$allelic_support[[tag_name]] <- stringr::str_replace_all(config$allelic_support[[tag_name]], "-", ".")
      if (config$allelic_support[[tag_name]] %in% colnames(vcf_df)) {
        if (tag_name == "control_af_tag") {
          vcf_df[, "AF_CONTROL"] <- round(as.numeric(vcf_df[, config$allelic_support[[tag_name]]]), digits = precision)
        }
        if (tag_name == "control_dp_tag") {
          vcf_df[, "DP_CONTROL"] <- as.integer(vcf_df[, config$allelic_support[[tag_name]]])
        }
        if (tag_name == "tumor_af_tag") {
          found_tumor_tag <- 1
          vcf_df[, "AF_TUMOR"] <- round(as.numeric(vcf_df[, config$allelic_support[[tag_name]]]), digits = precision)
        }
        if (tag_name == "tumor_dp_tag") {
          found_tumor_tag <- 1
          vcf_df[, "DP_TUMOR"] <- as.integer(vcf_df[, config$allelic_support[[tag_name]]])
        }
        if (tag_name == "call_conf_tag") {
          vcf_df[, "CALL_CONFIDENCE"] <- as.character(vcf_df[, config$allelic_support[[tag_name]]])
        }
      }
    }
  }
  return(vcf_df)
}

#' Function that adds SwissProt feature descriptions based on keys coming from PCGR pipeline
#'
#' @param vcf_data_df Data frame of sample variants from VCF
#' @param feature_descriptions Data frame with descriptions of functional features (active sites, disulphide bonds)
#' @return vcf_data_df
#'
#'
append_pfeature_descriptions <- function(vcf_data_df, feature_descriptions) {
  rlogging::message("Extending annotation descriptions related to UniprotKB/SwissProt protein features")

  invisible(assertthat::assert_that(!is.null(vcf_data_df), msg = "Argument 'vcf_data_df' cannot not be NULL"))
  invisible(assertthat::assert_that(is.data.frame(vcf_data_df), msg = paste0("Argument 'vcf_data_df' must be of type 'data.frame'")))
  assertable::assert_colnames(vcf_data_df, c("UNIPROT_FEATURE", "VAR_ID", "UNIPROT_ID", "CONSEQUENCE"), only_colnames = F, quiet = T)
  assertable::assert_colnames(feature_descriptions, c("UNIPROT_FEATURE", "PF", "UNIPROT_ID"), only_colnames = T, quiet = T)

  #if ("UNIPROT_FEATURE" %in% colnames(vcf_data_df) & "VAR_ID" %in% colnames(vcf_data_df) & "UNIPROT_ID" %in% colnames(vcf_data_df)) {
  feature_df <- dplyr::select(vcf_data_df, UNIPROT_FEATURE, VAR_ID, CONSEQUENCE, UNIPROT_ID) %>%
    dplyr::filter(!is.na(UNIPROT_FEATURE) & !is.na(UNIPROT_ID)) %>%
    dplyr::distinct()
  if (nrow(feature_df) == 0) {
    vcf_data_df$PROTEIN_FEATURE <- NA
    return(vcf_data_df)
  }
  feature_df <- as.data.frame(
    feature_df %>%
      tidyr::separate_rows(UNIPROT_FEATURE, sep = "&") %>%
      tidyr::separate_rows(UNIPROT_FEATURE, sep = ",") %>%
      dplyr::left_join(
        dplyr::select(feature_descriptions, UNIPROT_FEATURE, UNIPROT_ID, PF), by = c("UNIPROT_FEATURE", "UNIPROT_ID")
      ) %>%
      dplyr::filter(!is.na(PF)) %>%
      dplyr::group_by(VAR_ID, CONSEQUENCE) %>%
      dplyr::summarise(PROTEIN_FEATURE = paste(PF, collapse = ", ")) %>%
      dplyr::ungroup() %>%
      dplyr::mutate(PROTEIN_FEATURE = dplyr::if_else(PROTEIN_FEATURE == "NA",
                                                     as.character(NA),
                                                     as.character(PROTEIN_FEATURE))) %>%
      dplyr::mutate(PROTEIN_FEATURE = dplyr::if_else(!is.na(PROTEIN_FEATURE) & !stringr::str_detect(CONSEQUENCE, "^(synonymous|missense|stop|start)"),
                                                     as.character(NA),
                                                     as.character(PROTEIN_FEATURE))) %>%
      dplyr::filter(!is.na(PROTEIN_FEATURE))
  )

  if (nrow(feature_df) > 0) {
    vcf_data_df <- dplyr::select(vcf_data_df, -UNIPROT_FEATURE) %>%
      dplyr::left_join(feature_df, by = c("VAR_ID", "CONSEQUENCE"))
  }
  else{
    vcf_data_df$PROTEIN_FEATURE <- NA
  }

  return(vcf_data_df)

}

#' Function that adds GWAS citation/phenotype to GWAS hit found through PCGR annotation
#'
#' @param vcf_data_df Data frame of sample variants from VCF
#' @param pcgr_data PCGR data structure
#' @param p_value_threshold Required p-value to report associations from GWAS catalog
#'
#' @return vcf_data_df
#'
#'
append_gwas_citation_phenotype <- function(vcf_data_df, gwas_citations_phenotypes, p_value_threshold = 1e-6) {
  gwas_citations_phenotypes <- gwas_citations_phenotypes %>%
    dplyr::filter(p_value_num <= p_value_threshold)
  rlogging::message("Adding citations/phenotypes underlying GWAS hits (NHGRI-EBI GWAS Catalog)")

  invisible(assertthat::assert_that(!is.null(vcf_data_df), msg = "Argument 'vcf_data_df' cannot not be NULL"))
  invisible(assertthat::assert_that(is.data.frame(vcf_data_df), msg = paste0("Argument 'vcf_data_df' must be of type 'data.frame'")))
  invisible(assertthat::assert_that(!is.null(gwas_citations_phenotypes), msg = "Argument 'gwas_citations_phenotypes' cannot not be NULL"))
  invisible(assertthat::assert_that(is.data.frame(gwas_citations_phenotypes), msg = paste0("Argument 'gwas_citations_phenotypes' must be of type 'data.frame'")))
  assertable::assert_colnames(vcf_data_df, c("GWAS_HIT", "VAR_ID"), only_colnames = F, quiet = T)
  assertable::assert_colnames(gwas_citations_phenotypes, c("gwas_key", "GWAS_PH", "GWAS_CIT"), only_colnames = F, quiet = T)

  if ("GWAS_HIT" %in% colnames(vcf_data_df) & "VAR_ID" %in% colnames(vcf_data_df)) {
    feature_df <- dplyr::select(vcf_data_df, GWAS_HIT, VAR_ID) %>%
      dplyr::filter(!is.na(GWAS_HIT)) %>%
      dplyr::distinct()
    if (nrow(feature_df) == 0) {
      vcf_data_df$GWAS_CITATION <- NA
      vcf_data_df$GWAS_PHENOTYPE <- NA
      return(vcf_data_df)
    }
    feature_df <- as.data.frame(feature_df %>%
                                  tidyr::separate_rows(GWAS_HIT, sep = ",") %>%
                                  tidyr::separate(GWAS_HIT,
                                                  into = c("rsid", "risk_allele", "pmid",
                                                           "tagsnp", "p_value", "efo_id", "do_id"),
                                                  sep = "\\|", remove = F) %>%
                                  dplyr::mutate(gwas_key = paste(rsid, efo_id, pmid, sep = "_")) %>%
                                  dplyr::left_join(dplyr::select(gwas_citations_phenotypes,
                                                                 gwas_key, GWAS_PH, GWAS_CIT),
                                                   by = c("gwas_key")) %>%
                                  dplyr::filter(!is.na(GWAS_CIT)) %>%
                                  dplyr::group_by(VAR_ID) %>%
                                  dplyr::summarise(GWAS_PHENOTYPE = paste(unique(GWAS_PH), collapse = "; "),
                                                   GWAS_CITATION = paste(unique(GWAS_CIT), collapse = "; ")) %>%
                                  dplyr::mutate(GWAS_CITATION = dplyr::if_else(GWAS_CITATION == "NA",
                                                                               as.character(NA),
                                                                               as.character(GWAS_CITATION))) %>%
                                  dplyr::mutate(GWAS_PHENOTYPE = dplyr::if_else(GWAS_PHENOTYPE == "NA",
                                                                                as.character(NA),
                                                                                as.character(GWAS_PHENOTYPE))) %>%
                                  dplyr::filter(!is.na(GWAS_CITATION) & !is.na(GWAS_PHENOTYPE))
    )
    if (nrow(feature_df) > 0) {
      vcf_data_df <- vcf_data_df %>%
        dplyr::left_join(feature_df, by = c("VAR_ID" = "VAR_ID"))
    }else{
      vcf_data_df$GWAS_CITATION <- NA
      vcf_data_df$GWAS_PHENOTYPE <- NA
    }
  }else{
    vcf_data_df$GWAS_CITATION <- NA
    vcf_data_df$GWAS_PHENOTYPE <- NA
  }

  return(vcf_data_df)

}


#' Function that filters variant set on (depth, allelic fraction) for tumor and normal and filters according to settings
#'
#' @param vcf_df data frame with variants
#' @param config list with workflow configuration values
#' @param precision number of significant digits for allelic fraction estimation
#'
#' @return vcf_df
#'
filter_read_support <- function(vcf_df, config = NULL, precision = 3) {

  rlogging::message("Filtering tumor variants based on allelic depth/fraction (min_dp_tumor=", config$allelic_support$tumor_dp_min, ", min_af_tumor=", config$allelic_support$tumor_af_min, ")")
  rlogging::message("Filtering tumor variants based on allelic depth/fraction (min_dp_control=", config$allelic_support$control_dp_min, ", max_af_control=", config$allelic_support$control_af_max, ")")
  n_before_dp_af_filtering <- nrow(vcf_df)
  if (!any(is.na(vcf_df$DP_TUMOR))) {
    vcf_df <- dplyr::filter(vcf_df, DP_TUMOR >= config$allelic_support$tumor_dp_min)
  }
  if (!any(is.na(vcf_df$AF_TUMOR))) {
    vcf_df <- dplyr::filter(vcf_df, AF_TUMOR >= config$allelic_support$tumor_af_min)
  }
  if (!any(is.na(vcf_df$AF_CONTROL))) {
    vcf_df <- dplyr::filter(vcf_df, AF_CONTROL <= config$allelic_support$control_af_max)
  }
  if (!any(is.na(vcf_df$DP_CONTROL))) {
    vcf_df <- dplyr::filter(vcf_df, DP_CONTROL >= config$allelic_support$control_dp_min)
  }
  n_removed <- n_before_dp_af_filtering - nrow(vcf_df)
  percentage <- round(as.numeric((n_removed / n_before_dp_af_filtering) * 100), digits = 2)
  rlogging::message(paste0("Removed ", n_removed, " tumor variants (", percentage, "%) based on thresholds for allelic depth/fraction"))

  return(vcf_df)
}


#' A function that generates a HTML link for selected identifiers (DBSNP, COSMIC, CLINVAR, ENTREZ)
#'
#' @param df data frame
#' @param vardb type of database
#' @param group_by_var variable used for grouping (VAR_ID)
#' @param url_prefix url prefix for link generation
#' @param link_key_var variable used in url for linking
#' @param link_display_var variable used in url for display
#' @return df_annotation_links
#'

generate_annotation_link <- function(df, vardb = "DBSNP", group_by_var = "VAR_ID",
                                     url_prefix = "http://www.ncbi.nlm.nih.gov/projects/SNP/snp_ref.cgi?rs=",
                                     link_key_var = "DBSNPRSID", link_display_var = "DBSNPRSID",
                                     genome_assembly = "grch37") {

  df_annotation_links <- data.frame()
  invisible(assertthat::assert_that(is.data.frame(df), msg = paste0("Object df must be of type data.frame")))
  assertable::assert_colnames(df, c(group_by_var, link_key_var, link_display_var), only_colnames = F, quiet = T)

  if (group_by_var %in% colnames(df) & link_key_var %in% colnames(df) & link_display_var %in% colnames(df)) {
    #selected_vars <- dplyr::quos(unique(c(group_by_var, link_key_var, link_display_var)))
    selected_vars <- unique(c(group_by_var, link_key_var, link_display_var))
    tmp_df <- df[,selected_vars]
    tmp_df <- tmp_df[!is.na(tmp_df[,link_key_var]) & !is.na(tmp_df[,link_display_var]),] %>%
      dplyr::distinct()

    # tmp_df <- df %>%
    #   dplyr::select(!!enquo(selected_vars)) %>%
    #   dplyr::filter(!is.na(!!enquo(link_key_var)) & !is.na(!!enquo(link_display_var))) %>%
    #   dplyr::distinct()
    if (nrow(tmp_df) > 0) {
      df_annotation_links <- data.frame()
      if(vardb == "DBSNP"){
        df_annotation_links <- as.data.frame(
          tmp_df %>%
            tidyr::separate_rows(DBSNPRSID, sep ="&|,") %>%
            dplyr::mutate(tmp_link = paste0("<a href='",url_prefix,DBSNPRSID, "'target=\"_blank\">",DBSNPRSID,"</a>"))
        )
      }

      if(vardb == "CLINVAR"){
        df_annotation_links <- as.data.frame(
          tmp_df %>%
            tidyr::separate_rows(CLINVAR_MSID, sep ="&|,") %>%
            dplyr::mutate(tmp_link = paste0("<a href='",url_prefix,CLINVAR_MSID, "'target=\"_blank\">",CLINVAR_TRAITS_ALL,"</a>"))
        )
      }

      if(vardb == "GENE_NAME"){
        df_annotation_links <- as.data.frame(
          tmp_df %>%
            tidyr::separate_rows(ENTREZ_ID, sep ="&|,") %>%
            dplyr::mutate(tmp_link = paste0("<a href='",url_prefix,ENTREZ_ID, "'target=\"_blank\">",GENENAME,"</a>"))
        )
      }

      if(vardb == "PROTEIN_DOMAIN"){
        df_annotation_links <- as.data.frame(
          tmp_df %>%
            tidyr::separate_rows(PFAM_DOMAIN, sep ="&|,") %>%
            dplyr::mutate(tmp_link = paste0("<a href='",url_prefix,PFAM_DOMAIN, "'target=\"_blank\">",PFAM_DOMAIN_NAME,"</a>"))
        )
      }

      if(vardb == "COSMIC"){
        df_annotation_links <- as.data.frame(
          tmp_df %>%
            tidyr::separate_rows(COSMIC_MUTATION_ID, sep ="&|,") %>%
            dplyr::mutate(tmp_link = paste0("<a href='",url_prefix,COSMIC_MUTATION_ID, "'target=\"_blank\">",COSMIC_MUTATION_ID,"</a>"))
        )
      }

      if(vardb == "NCBI_REFSEQ"){
        df_annotation_links <- as.data.frame(
          tmp_df %>%
            tidyr::separate_rows(REFSEQ_MRNA, sep ="&|,") %>%
            dplyr::mutate(tmp_link = paste0("<a href='",url_prefix,REFSEQ_MRNA, "'target=\"_blank\">",REFSEQ_MRNA,"</a>"))
        )
      }

      if(nrow(df_annotation_links) > 0){
        df_annotation_links <- as.data.frame(
          df_annotation_links %>%
            dplyr::group_by(VAR_ID) %>%
            dplyr::summarise(link = unlist(paste(tmp_link, collapse=", ")), .groups = "drop")
        )
      }




      # tmp_df_melted <- tmp_df %>%
      #   tidyr::separate_rows(!!enquo(link_key_var), sep = "&|,") %>%
      #   dplyr::mutate(tmp_link = paste0("<a href='", url_prefix, !!enquo(link_key_var), "' target=\"_blank\">", !!rlang::sym(link_display_var), "</a>"))

      # df_annotation_links <- as.data.frame(
      #   tmp_df_melted %>%
      #     dplyr::group_by(!! rlang::sym(group_by_var)) %>%
      #     dplyr::summarise(link = unlist(paste(tmp_link, collapse = ", ")))
      #   )
    }
  }
  else{
    cat("WARNING: Could not generate HTML URL links - missing url_key or grouping varable in df", sep = "\n")
  }
  return(df_annotation_links)
}


#' Function that converts a list to a data frame
#'
#' @param listfordf list to be converted
#'
#' @return df
#'
list_to_df <- function(listfordf) {
  if (!is.list(listfordf)) stop("it should be a list")

  df <- list(list.element = listfordf)
  class(df) <- c("tbl_df", "data.frame")
  attr(df, "row.names") <- .set_row_names(length(listfordf))

  if (!is.null(names(listfordf))) {
    df$name <- names(listfordf)
  }

  return(as.data.frame(df))
}




#' Function that assigns genotype (het/hom) from VCF GT tag
#'
#' @param vcf_df
#'
#' @return vcf_df
#'
determine_genotype <- function(vcf_df) {

  invisible(assertthat::assert_that(is.data.frame(vcf_df), msg = paste0("Argument 'vcf_df' must be of type data.frame, not", class(vcf_df))))
  vcf_df$GENOTYPE <- "NA"
  if ("GT" %in% colnames(vcf_df)) {
    vcf_df <- vcf_df %>%
      dplyr::mutate(GENOTYPE = dplyr::if_else(GT %in% pcgrr::heterozygous_states, "heterozygous", as.character("ND"))) %>%
      dplyr::mutate(GENOTYPE = dplyr::if_else(GT %in% pcgrr::homozygous_states, "homozygous", as.character(GENOTYPE)))
  }
  return(vcf_df)
}


data_integrity_check <- function(vcf_df, pcgr_data, workflow = "pcgr") {
  rlogging::message("Verifying data integrity of input callset")

  stopifnot(is.data.frame(vcf_df) & !is.null(pcgr_data))
  stopifnot(!is.null(pcgr_data[["annotation_tags"]][["vcf_cpsr"]]) &
              !is.null(pcgr_data[["annotation_tags"]][["vcf_pcgr"]]))

  vars_required <- pcgr_data[["annotation_tags"]][["vcf_cpsr"]]
  if (workflow == "pcgr") {
    vars_required <- pcgr_data[["annotation_tags"]][["vcf_pcgr"]]
    vars_required <- vars_required[!vars_required$tag == "PANEL_OF_NORMALS", ]
  }

  vcf_required_vars <- dplyr::bind_rows(data.frame(tag = "CHROM", stringsAsFactors = F),
                             data.frame(tag = "REF", stringsAsFactors = F),
                             data.frame(tag = "ALT", stringsAsFactors = F),
                             data.frame(tag = "POS", stringsAsFactors = F),
                             data.frame(tag = "QUAL", stringsAsFactors = F),
                             data.frame(tag = "FILTER", stringsAsFactors = F))
  vcf_required_vars$tag <- as.character(vcf_required_vars$tag)
  vars_required <- dplyr::bind_rows(vars_required, vcf_required_vars)
  i <- 1
  while (i <= nrow(vars_required)) {
    tag <- vars_required[i, "tag"]
    if (!(tag %in% colnames(vcf_df))) {
      if (workflow == "pcgr") {
        rlogging::stop(paste0("Missing required variable (", tag, ") in annotated TSV file from PCGR workflow - quitting"))
      }else{
        rlogging::stop(paste0("Missing required variable (", tag, ") in annotated TSV file from CPSR workflow - quitting"))
      }
    }
    i <- i + 1
  }



  return(vcf_df)
}

recode_logical_vars <- function(vcf_data_df){
  for (v in c("ONCOGENE", "TUMOR_SUPPRESSOR", "NETWORK_CG", "LAST_EXON", "LAST_INTRON", "PANEL_OF_NORMALS",
              "NULL_VARIANT", "SPLICE_DONOR_RELEVANT", "WINMASKER_HIT", "SIMPLEREPEATS_HIT")) {
    if (v %in% colnames(vcf_data_df)) {
      vcf_data_df[, v] <- as.logical(dplyr::recode(vcf_data_df[, v], True = TRUE, False = FALSE))
    }
  }
  return(vcf_data_df)

}

append_annotation_links <- function(vcf_data_df, pcgr_data = NULL, skip = NULL){
  i <- 1
  while (i <= nrow(pcgrr::variant_db_url)) {

    name <- pcgrr::variant_db_url[i, ]$name
    if(!name %in% skip){
      rlogging::message(paste0("Adding annotation links - ", name))
      group_by_var <- pcgrr::variant_db_url[i, ]$group_by_var
      url_prefix <- pcgrr::variant_db_url[i, ]$url_prefix
      link_key_var <- pcgrr::variant_db_url[i, ]$link_key_var
      link_display_var <- pcgrr::variant_db_url[i, ]$link_display_var
      if (!(name %in% colnames(vcf_data_df))) {
        annotation_links <-
          pcgrr::generate_annotation_link(vcf_data_df, vardb = name, group_by_var = group_by_var,
                                          url_prefix = url_prefix,
                                          link_key_var = link_key_var, link_display_var = link_display_var,
                                          genome_assembly = pcgr_data[["assembly"]][["grch_name"]])
        if (nrow(annotation_links) > 0) {
          vcf_data_df <- vcf_data_df %>%
            dplyr::left_join(dplyr::rename(annotation_links, !!rlang::sym(name) := link), by = c("VAR_ID"))
        }else{
          vcf_data_df[, name] <- NA
        }
      }
    }
    i <- i + 1
  }
  return(vcf_data_df)
}

#' Function that reads a fully annotated VCF from PCGR VEP/vcfanno pipeline
#'
#' @param tsv_gz_file Bgzipped VCF file
#' @param pcgr_data object with PCGR annotation data
#' @param sample_name sample identifier
#' @param config Object with PCGR configuration parameters
#' @param oncotree data frame with phenotype terms from UMLS that is relevant for the sample
#' @param cpsr logical indicating CPSR workflow
#'
#' @return vcf_data_df
#'
get_calls <- function(tsv_gz_file, pcgr_data, sample_name,
                      config, oncotree = NULL, cpsr = F, n_lines_skip = 1,
                      filenames = NULL) {

  ## check that arguments are valid
  invisible(assertthat::assert_that(!is.null(sample_name), msg = "Argument 'sample_name' cannot not be NULL"))
  invisible(assertthat::assert_that(class(sample_name) == "character",
                                    msg = paste0("Argument 'sample_name' must be of type 'character', not ", class(sample_name))))
  invisible(assertthat::assert_that(!is.null(tsv_gz_file), msg = "Argument 'tsv_gz_file' cannot not be NULL"))
  invisible(assertthat::assert_that(file.exists(tsv_gz_file),
                                    msg = paste0("File ", tsv_gz_file, " does not exist")))
  invisible(assertthat::assert_that(nat.utils::is.gzip(tsv_gz_file) == T,
                                    msg = paste0("File ", tsv_gz_file, " is not a gzipped TSV file")))
  invisible(assertthat::assert_that(!is.null(pcgr_data), msg = "Argument 'pcgr_data' cannot not be NULL"))
  invisible(assertthat::assert_that(class(pcgr_data) == "list",
                                    msg = paste0("Argument 'pcgr_data' must be of type list, not ", class(pcgr_data))))
  invisible(assertthat::assert_that(!is.null(config), msg = "Argument 'config' cannot not be NULL"))
  invisible(assertthat::assert_that(class(config)[2] == "list",
                                    msg = paste0("Argument 'config' must be of type list, not ", class(config)[2])))
  invisible(assertthat::assert_that(!is.null(oncotree), msg = "Argument 'oncotree' cannot not be NULL"))

  ## read VCF data
  vcf_data_df <- read.table(gzfile(tsv_gz_file), skip = n_lines_skip,
                            sep = "\t", header = T, stringsAsFactors = F,
                            quote = "", comment.char = "", na.strings = c("."))
  if (nrow(vcf_data_df) == 0) return(vcf_data_df)

  ## Fix erroneous parsing of cases where all reference/alternate alleles are "T",
  ## converted to logical with read.table
  if(is.logical(vcf_data_df$REF)){
    vcf_data_df$REF <- as.character("T")
  }
  if(is.logical(vcf_data_df$ALT)){
    vcf_data_df$ALT <- as.character("T")
  }


  wflow <- "pcgr"
  if (cpsr == T) {
    wflow <- "cpsr"
  }

  ## Annotation variables added by get_calls()
  pcgr_columns <- c("GENOME_VERSION", "PROTEIN_CHANGE",
                    "CONSEQUENCE", "GENOMIC_CHANGE",
                    "VAR_ID", "OPENTARGETS_ASSOCIATIONS",
                    "DOCM_DISEASE", "DOCM_LITERATURE",
                    "CLINVAR", "CLINVAR_TRAITS_ALL",
                    "GENE_NAME", "GENENAME", "CANCERGENE_SUPPORT",
                    "OPENTARGETS_RANK", "TARGETED_DRUGS",
                    "CANCER_ASSOCIATIONS", "DBSNP",
                    "COSMIC", "PROTEIN_DOMAIN",
                    "CLINVAR_PHENOTYPE", "NCBI_REFSEQ",
                    "AF_TUMOR", "DP_TUMOR",
                    "AF_CONTROL", "DP_CONTROL",
                    "CALL_CONFIDENCE", "PFAM_DOMAIN_NAME",
                    "GWAS_CITATION", "GWAS_PHENOTYPE")
  vcf_data_df <- vcf_data_df[, !(colnames(vcf_data_df) %in% pcgr_columns)]

  ## Perform data integrity check of input
  vcf_data_df <- vcf_data_df %>%
    pcgrr::data_integrity_check(pcgr_data = pcgr_data, workflow = wflow) %>%
    dplyr::mutate(GENOMIC_CHANGE = paste0(CHROM, ":g.", POS, REF, ">", ALT)) %>%
    dplyr::mutate(VAR_ID = paste(CHROM, POS, REF, ALT, sep = "_")) %>%
    pcgrr::detect_vcf_sample_name(cpsr = cpsr, sample_name = sample_name) %>%
    pcgrr::get_ordinary_chromosomes(chrom_var = "CHROM")
  if (nrow(vcf_data_df) == 0) return(vcf_data_df)

  vcf_data_df <- vcf_data_df %>%
    pcgrr::order_variants() %>%
    dplyr::rename(CONSEQUENCE = Consequence)

  ##
  if (cpsr == T) {
    if ("LoF" %in% colnames(vcf_data_df)) {
      vcf_data_df <- vcf_data_df %>%
        dplyr::mutate(LOSS_OF_FUNCTION = dplyr::if_else(!is.na(LoF) & LoF == "HC",
                                                        TRUE, FALSE, FALSE)) %>%
        ## Ignore LoF predictions for missense variants (bug in LofTee?)
        dplyr::mutate(LOSS_OF_FUNCTION = dplyr::if_else(!is.na(CONSEQUENCE) &
                                                          CONSEQUENCE == "missense_variant" &
                                                          LOSS_OF_FUNCTION == T, FALSE, LOSS_OF_FUNCTION))
    }
    if (!is.null(config[["gwas"]][["p_value_min"]])) {
      vcf_data_df <- vcf_data_df %>%
        pcgrr::append_gwas_citation_phenotype(gwas_citations_phenotypes = pcgr_data[["gwas"]][["citations_phenotypes"]],
                                           p_value_threshold = config[["gwas"]][["p_value_min"]])
    }
    vcf_data_df <- vcf_data_df %>%
      pcgrr::determine_genotype()

  }else{
    if (!("PANEL_OF_NORMALS" %in% colnames(vcf_data_df))) {
      vcf_data_df <- vcf_data_df %>%
        dplyr::mutate(PANEL_OF_NORMALS = "False")
    }
    vcf_data_df <- vcf_data_df %>%
      pcgrr::append_read_support(config = config)

    if("maf_tmp" %in% names(filenames) & "maf" %in% names(filenames)){
      pcgrr::update_maf_allelic_support(vcf_data_df, filenames[["maf_tmp"]], filenames[["maf"]])
    }

   vcf_data_df <- vcf_data_df %>%
      pcgrr::filter_read_support(config = config) %>%
      dplyr::mutate(PUTATIVE_DRIVER_MUTATION = dplyr::if_else(!is.na(PUTATIVE_DRIVER_MUTATION),
                                                              TRUE, FALSE))
  }

  if (nrow(vcf_data_df) == 0) return(vcf_data_df)

  vcf_data_df <- vcf_data_df %>%
    dplyr::mutate(SYMBOL = SYMBOL_ENTREZ) %>%
    dplyr::mutate(CLINVAR_CONFLICTED = dplyr::case_when(CLINVAR_CONFLICTED == "1" ~ TRUE,
                                                        CLINVAR_CONFLICTED == "0" ~ FALSE,
                                                        FALSE ~ as.logical(CLINVAR_CONFLICTED))) %>%
    dplyr::mutate(GENOME_VERSION = pcgr_data[["assembly"]][["grch_name"]], PROTEIN_CHANGE = HGVSp) %>%
    dplyr::mutate(PROTEIN_CHANGE = dplyr::if_else(stringr::str_detect(PROTEIN_CHANGE, ":"),
                                                  stringr::str_split_fixed(PROTEIN_CHANGE, pattern = ":", 2)[, 2],
                                                  as.character(PROTEIN_CHANGE))) %>%
    dplyr::mutate(PROTEIN_CHANGE = dplyr::if_else(stringr::str_detect(PROTEIN_CHANGE, "^ENSP"),
                                                  as.character(NA),
                                                  as.character(PROTEIN_CHANGE))) %>%
    dplyr::mutate(CLINVAR_MSID = as.integer(CLINVAR_MSID)) %>%
    dplyr::mutate(EXON = dplyr::if_else(!is.na(EXON) & stringr::str_detect(EXON,"^([0-9]{1,}/[0-9]{1,})$"),
                                        as.integer(stringr::str_split_fixed(EXON,"/",2)[,1]),
                                        as.integer(NA))) %>%
    pcgrr::append_pfeature_descriptions(feature_descriptions = pcgr_data[["protein_features"]][["swissprot"]]) %>%
    dplyr::mutate(PFAM_DOMAIN = as.character(stringr::str_replace(PFAM_DOMAIN, "Pfam_domain:", ""))) %>%
    dplyr::left_join(dplyr::select(pcgr_data[["protein_domains"]][["pfam"]], pfam_id, pfam_name), by = c("PFAM_DOMAIN" = "pfam_id")) %>%
    dplyr::rename(PFAM_DOMAIN_NAME = pfam_name) %>%
    dplyr::left_join(pcgr_data[["biomarkers"]][["docm"]], by = c("VAR_ID")) %>%
    dplyr::mutate(ENTREZ_ID = as.character(ENTREZ_ID)) %>%
    dplyr::mutate(Gene = as.character(Gene)) %>%
    pcgrr::recode_logical_vars()



  rlogging::message(paste0("Number of PASS variants: ", nrow(vcf_data_df)))
  if (any(grepl(paste0("VARIANT_CLASS$"), names(vcf_data_df)))) {
    n_snvs <-
      vcf_data_df %>% dplyr::filter(!is.na(VARIANT_CLASS) & VARIANT_CLASS == "SNV") %>% nrow()
    n_deletions <- vcf_data_df %>%
      dplyr::filter(!is.na(VARIANT_CLASS) & (VARIANT_CLASS == "deletion" | VARIANT_CLASS == "indel")) %>%
      nrow()
    n_insertions <- vcf_data_df %>%
      dplyr::filter(!is.na(VARIANT_CLASS) & VARIANT_CLASS == "insertion") %>%
      nrow()
    n_substitutions <- vcf_data_df %>%
      dplyr::filter(!is.na(VARIANT_CLASS) & VARIANT_CLASS == "substitution") %>%
      nrow()
    rlogging::message(paste0("Number of SNVs: ", n_snvs))
    rlogging::message(paste0("Number of deletions: ", n_deletions))
    rlogging::message(paste0("Number of insertions: ", n_insertions))
    rlogging::message(paste0("Number of block substitutions: ", n_substitutions))
  }



  if (nrow(vcf_data_df) == 0)return(vcf_data_df)

  vcf_data_df_1 <-
    dplyr::left_join(dplyr::filter(vcf_data_df, !is.na(ENTREZ_ID)),
                     dplyr::filter(dplyr::select(pcgr_data[["gene_xref"]][["gencode"]],
                                                 ENTREZ_ID, Gene, GENENAME, CANCERGENE_SUPPORT),
                                   !is.na(ENTREZ_ID)),
                     by = c("ENTREZ_ID", "Gene"))
  vcf_data_df_2 <- dplyr::left_join(dplyr::filter(vcf_data_df, is.na(ENTREZ_ID)),
                                    dplyr::select(pcgr_data[["gene_xref"]][["gencode"]],
                                                  Gene, GENENAME, CANCERGENE_SUPPORT),
                                    by = c("Gene"))
  vcf_data_df <- rbind(vcf_data_df_1, vcf_data_df_2) %>%
    pcgrr::order_variants()

  rlogging::message("Extending annotation descriptions related to KEGG pathways")
  vcf_data_df <- dplyr::left_join(vcf_data_df, pcgr_data[["kegg"]][["pathway_links"]],
                                  by = c("ENTREZ_ID" = "gene_id")) %>%
    dplyr::rename(KEGG_PATHWAY = kegg_pathway_urls)

  clinvar <- dplyr::select(pcgr_data[["clinvar"]][["variants"]],
                           CLINVAR_TRAITS_ALL, CLINVAR_MSID, VAR_ID) %>%
    dplyr::filter(!is.na(CLINVAR_MSID) & !is.na(VAR_ID))
  if ("CLINVAR_MSID" %in% colnames(vcf_data_df) & "VAR_ID" %in% colnames(vcf_data_df)) {
    rlogging::message("Extending annotation descriptions related to ClinVar")
    vcf_data_df <- dplyr::left_join(vcf_data_df, clinvar, by = c("CLINVAR_MSID", "VAR_ID")) %>%
      dplyr::mutate(CLINVAR_MSID = as.character(CLINVAR_MSID))
  }

  vcf_data_df <-
    pcgrr::df_string_replace(vcf_data_df, strings = c("CONSEQUENCE", "REFSEQ_MRNA"),
                             pattern = "&", replacement = ", ",
                             replace_all = T)
  vcf_data_df <-
    pcgrr::df_string_replace(vcf_data_df,
                             strings = c("VEP_ALL_CSQ", "DOCM_DISEASE",
                                         "MUTATION_HOTSPOT_CANCERTYPE",
                                         "ICGC_PCAWG_OCCURRENCE"),
                             pattern = ",", replacement = ", ", replace_all = T)


  if ("EFFECT_PREDICTIONS" %in% colnames(vcf_data_df)) {
    vcf_data_df <- vcf_data_df %>%
      dplyr::mutate(EFFECT_PREDICTIONS = stringr::str_replace_all(EFFECT_PREDICTIONS, "\\.&|\\.$", "NA&"))  %>%
      dplyr::mutate(EFFECT_PREDICTIONS = stringr::str_replace_all(EFFECT_PREDICTIONS, "&$", "")) %>%
      dplyr::mutate(EFFECT_PREDICTIONS = stringr::str_replace_all(EFFECT_PREDICTIONS, "&", ", ")) %>%
      pcgrr::append_dbnsfp_var_link()
  }


  if ("TCGA_FREQUENCY" %in% colnames(vcf_data_df)) {
    vcf_data_df <- vcf_data_df %>%
      pcgrr::append_tcga_var_link(pcgr_data = pcgr_data)
  }

  vcf_data_df <- vcf_data_df %>%
    pcgrr::append_annotation_links(pcgr_data = pcgr_data)

  if (!("TARGETED_DRUGS" %in% colnames(vcf_data_df))) {
    vcf_data_df <- vcf_data_df %>%
      pcgrr::append_drug_var_link(pcgr_data = pcgr_data) %>%
      dplyr::rename(TARGETED_DRUGS = ANTINEOPHARMALINK)
  }


  if (!("CANCER_ASSOCIATIONS" %in% colnames(vcf_data_df))) {
    vcf_data_df <- vcf_data_df %>%
      pcgrr::append_disgenet_var_link(pcgr_data = pcgr_data) %>%
      dplyr::rename(CANCER_ASSOCIATIONS = DISGENET_LINK)
  }

  if (!("OPENTARGETS_ASSOCIATIONS" %in% colnames(vcf_data_df))) {
    vcf_data_df <- vcf_data_df %>%
      pcgrr::append_otargets_pheno_link(pcgr_data = pcgr_data, oncotree = oncotree) %>%
      dplyr::rename(OPENTARGETS_ASSOCIATIONS = OT_DISEASE_LINK)
  }


  return(vcf_data_df)

}



write_processed_vcf <- function(calls, sample_name = NULL,
                                output_directory = NULL,
                                vcf_fname = NULL){

  rlogging::message("Writing VCF file with input calls for signature analysis")

  header_lines <- c("##fileformat=VCFv4.2",
                    "##INFO=<ID=SAMPLE_ID,Number=1,Type=String,Description=\"Sample identifier\">",
                    "#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO")

  vcf_df <- calls %>%
    dplyr::select(CHROM, POS, REF, ALT) %>%
    dplyr::mutate(QUAL = ".", FILTER = "PASS", ID = ".",
                  INFO = paste0("SAMPLE_ID=",sample_name)) %>%
    dplyr::distinct()

  sample_vcf_content_fname <- file.path(output_directory,
                                        paste0(sample_name, ".",
                                               sample(100000,1),".vcf_content.tsv"))
  write(header_lines,file = vcf_fname,sep="\n")

  sample_vcf <- vcf_df[,c("CHROM","POS","ID","REF","ALT","QUAL","FILTER","INFO")] %>%
    dplyr::mutate(REF = as.character(REF)) %>%
    dplyr::mutate(ALT = as.character(ALT))

  options(scipen = 999)
  write.table(sample_vcf, file=sample_vcf_content_fname,sep="\t",col.names = F,quote=F, row.names = F)

  #system(paste0("cat ", sample_vcf_content_fname," | egrep -v \"^[XYM]\" | sort -k1,1n -k2,2n -k4,4 -k5,5 | awk '{printf(\"%s\\n\",$0);}' >> ",vcf_fname))
  #system(paste0("cat ", sample_vcf_content_fname," | egrep \"^[XYM]\" | sort -k1,1 -k2,2n -k4,4 -k5,5 | awk '{printf(\"%s\\n\",$0);}' >> ",vcf_fname))

  system(paste0("cat ", sample_vcf_content_fname," | egrep -v \"^[XYM]\" | sort -k1,1n -k2,2n -k4,4 -k5,5  >> ",vcf_fname))
  system(paste0("cat ", sample_vcf_content_fname," | egrep \"^[XYM]\" | sort -k1,1 -k2,2n -k4,4 -k5,5 >> ",vcf_fname))

  system(paste0("bgzip -f ", vcf_fname))
  system(paste0("tabix -p vcf ", vcf_fname,".gz"))

  system(paste0('rm -f ', sample_vcf_content_fname))
}


#' A function that detects whether the sample name in variant data frame is unique, throws an error if
#' multiple sample names are present for the CPSR workflow
#
#' @param df VCF data frame
#' @param cpsr logical indicating CPSR workflow
#' @return df Vranges object
#'
#'
detect_vcf_sample_name <- function(df, sample_name = NULL, cpsr = FALSE) {
  stopifnot(is.data.frame(df) & !is.null(sample_name))
  if ("VCF_SAMPLE_ID" %in% colnames(df)) {
    unique_sample_names <- unique(df$VCF_SAMPLE_ID)
    rlogging::message(paste0("Found the following VCF sample names: ", paste(unique_sample_names, collapse = ", ")))

    if (length(unique_sample_names) > 1 & cpsr == T) {
      rlogging::stop("Found more than one sample name - VCF with somatic calls? Expecting single sample germline VCF for CPSR")
    }
  }
  df <- df %>% dplyr::mutate(VCF_SAMPLE_ID = sample_name)
  return(df)
}

parse_transvar_file <- function(transvar_output_fname) {

  transvar_output_raw <- as.data.frame(
    readr::read_tsv(transvar_output_fname, col_names = T) %>%
    janitor::clean_names() %>%
    tidyr::separate(coordinates_g_dna_c_dna_protein, c("gdna", "cdna", "hgvsp"), sep = "/") %>%
    dplyr::rename(transcript_id = transcript, transvar_id = input) %>%
    dplyr::mutate(transcript_id = stringr::str_replace(transcript_id, " \\(protein_coding\\)", "")) %>%
    dplyr::filter(stringr::str_detect(gdna, ">|del[A-Z]{1,}ins[A-Z]{1,}")) %>%
    dplyr::distinct()
  )

  transvar_gdna_snv <- transvar_output_raw %>%
    dplyr::select(gdna, transcript_id, transvar_id) %>%
    dplyr::distinct()

  transvar_gdna_mnv <- dplyr::select(transvar_output_raw, transvar_id, transcript_id, info) %>%
    dplyr::mutate(info = stringr::str_replace_all(info, "CSQN=\\S+;reference_codon=\\S+;candidate_codons=((A|C|G|T){3},?){1,};", "")) %>%
    dplyr::mutate(info = stringr::str_replace_all(info, "source=GENCODE$", "NA")) %>%
    dplyr::mutate(info = stringr::str_replace_all(info, "candidate_mnv_variants=|candidate_snv_variants=", "")) %>%
    dplyr::mutate(info = stringr::str_replace_all(info, ";?aliases=ENSP[0-9]{1,}\\.[0-9]{1,}", "")) %>%
    dplyr::mutate(info = stringr::str_replace_all(info, ";source=GENCODE$", "")) %>%
    dplyr::mutate(info = stringr::str_replace_all(info, ";", ",")) %>%
    dplyr::filter(!is.na(info) & info != "NA" & nchar(info) > 0) %>%
    dplyr::rename(gdna = info) %>%
    tidyr::separate_rows(gdna, sep = ",") %>%
    dplyr::filter(gdna != "NA")

  transvar_output <- as.data.frame(
      dplyr::bind_rows(transvar_gdna_snv, transvar_gdna_mnv) %>%
      dplyr::group_by(gdna, transvar_id) %>%
      dplyr::summarise(transvar_transcript_id = paste(transcript_id, collapse = ",")))

  ##insertion/deletions
  transvar_output_indels <- transvar_output %>% dplyr::filter(stringr::str_detect(gdna, "del|ins"))
  transvar_output_indels$chr_start <-
    stringr::str_split_fixed(stringr::str_replace(transvar_output_indels$gdna, "chr.{1,}:g\\.", ""), "_", n = 2)[, 1]
  alleles <-
    stringr::str_split_fixed(stringr::str_replace(transvar_output_indels$gdna, "chr([0-9]{1,}|X|Y):g\\.[0-9]{1,}_[0-9]{1,}del", ""), "ins", 2)
  transvar_output_indels$refbase <- alleles[, 1]
  transvar_output_indels$altbase <- alleles[, 2]
  transvar_output_indels$chr_stop <- as.integer(transvar_output_indels$chr_start) + nchar(transvar_output_indels$altbase) - 1
  transvar_output_indels$chromosome <- stringr::str_replace(stringr::str_split_fixed(transvar_output_indels$gdna, ":", 2)[, 1], "chr", "")

  ##snvs
  transvar_output_snvs <- transvar_output %>% dplyr::filter(!stringr::str_detect(gdna, "del|ins"))
  transvar_output_snvs$chr_start <-
    stringr::str_replace(stringr::str_replace(transvar_output_snvs$gdna, "chr([0-9]{1,}|X|Y):g\\.", ""), "(A|G|C|T)>(A|G|C|T)$", "")
  transvar_output_snvs$chr_stop <- transvar_output_snvs$chr_start
  alleles <- stringr::str_split_fixed(stringr::str_replace(transvar_output_snvs$gdna, "chr([0-9]{1,}|X|Y):g\\.[0-9]{1,}", ""), ">", 2)
  transvar_output_snvs$refbase <- alleles[, 1]
  transvar_output_snvs$altbase <- alleles[, 2]
  transvar_output_snvs$chromosome <- stringr::str_replace(stringr::str_split_fixed(transvar_output_snvs$gdna, ":", 2)[, 1], "chr", "")


  all_transvar <- rbind(transvar_output_indels, transvar_output_snvs)
  all_transvar <- dplyr::filter(all_transvar, !stringr::str_detect(gdna, "Error|no_valid_transcript_found"))

  return(all_transvar)

}

#' #' A function that reads a VCF file into a VRanges object
#' #'
#' #' This function runs the VariantAnnotation::readVcfAsVranges
#' #'
#' #' @param vcf_file VCF file
#' #' @param genotype_fields Character vector with genotype fields to be imported into the VRanges object
#' #' @param samples samples to include
#' #' @return vr Vranges object
#' #' @examples
#' #' read_vcf_vranges('chr1.vcf.gz',c('GT','DP','AD'))
#' #'
#' #'
#' read_vcf_vranges <- function(vcf_file, genotype_fields = NULL,
#'                              samples = NULL, bsg = BSgenome.Hsapiens.UCSC.hg19) {
#'
#'   info_tags <- rownames(VariantAnnotation::info(VariantAnnotation::scanVcfHeader(vcf_file)))
#'   vcfparam <- VariantAnnotation::ScanVcfParam(fixed = c("ALT", "FILTER"), info = info_tags)
#'   if (!is.null(genotype_fields)) {
#'     vcfparam <- VariantAnnotation::ScanVcfParam(fixed = c("ALT", "FILTER"),
#'                                                 geno = genotype_fields, info = info_tags)
#'     if (!is.null(samples)) {
#'       vcfparam <- VariantAnnotation::ScanVcfParam(fixed = c("ALT", "FILTER"),
#'                                                   geno = genotype_fields, info = info_tags,
#'                                                   samples = samples)
#'     }
#'   }
#'   vr <- VariantAnnotation::readVcfAsVRanges(vcf_file, "hg19", param = vcfparam)
#'   seqlevels(vr) <- as.character(GenomeInfoDb::mapSeqlevels(unique(seqlevels(vr)), "UCSC", best.only = TRUE, drop = TRUE))
#'   seqlengths(vr) <- seqlengths(BSgenome.Hsapiens.UCSC.hg19)[seqlevels(bsg) %in% unique(seqlevels(vr))]
#'   strand(vr) <- "+"
#'
#'
#'   return(vr)
#' }


update_maf_allelic_support <- function(calls, maf_fname_tmp,
                                       maf_fname, delete_raw = T) {

  assertable::assert_colnames(calls,c("CHROM","POS","DP_TUMOR","AF_TUMOR",
                                      "DP_CONTROL","AF_CONTROL","VARIANT_CLASS"),
                              only_colnames = F, quiet = T)
  rlogging::message("Updating MAF file with information regarding variant sequencing depths and allelic support")
  maf_data <- NULL
  if (file.exists(maf_fname_tmp) & file.info(maf_fname_tmp)$size > 0) {
    maf_data <- read.table(maf_fname_tmp, skip = 1, header = T,
                           sep = "\t", na.strings = c(""), stringsAsFactors = F,
                           quote = "", comment.char = "") %>%
      dplyr::mutate(Chromosome = as.character(Chromosome))

  allele_columns <- c("Reference_Allele","Tumor_Seq_Allele1","Tumor_Seq_Allele2",
                     "Match_Norm_Seq_Allele1","Match_Norm_Seq_Allele2", "Tumor_Validation_Allele1",
                     "Tumor_Validation_Allele2","Match_Norm_Validation_Allele1",
                     "Match_Norm_Validation_Allele2")
  for(c in allele_columns){
    if(!is.null(maf_data[,c])){
      if(is.logical(maf_data[,c])){
        maf_data[,c] <- as.character("T")
      }
    }
  }

  calls_maf <- calls %>%
    dplyr::select(CHROM, POS, DP_TUMOR, AF_TUMOR, DP_CONTROL, AF_CONTROL, VARIANT_CLASS) %>%
    dplyr::rename(Chromosome = CHROM, Start_Position = POS) %>%
    dplyr::mutate(Chromosome = as.character(Chromosome)) %>%
    dplyr::mutate(Start_Position = dplyr::if_else(VARIANT_CLASS == "deletion",
                                                  Start_Position + 1,
                                                  as.double(Start_Position)))

    if (!is.null(maf_data)) {
      if (!any(is.na(calls_maf$DP_TUMOR)) & !any(is.na(calls_maf$AF_TUMOR))) {
        calls_maf$t_depth_estimate <- calls_maf$DP_TUMOR
        calls_maf$t_ref_count_estimate <- calls_maf$DP_TUMOR - round(calls_maf$AF_TUMOR * calls_maf$DP_TUMOR, digits = 0)
        calls_maf$t_alt_count_estimate <- calls_maf$DP_TUMOR - calls_maf$t_ref_count_estimate
        calls_maf <- calls_maf %>% dplyr::select(-c(DP_TUMOR, AF_TUMOR)) %>%
          dplyr::mutate(Chromosome = as.character(Chromosome))

        maf_data <- maf_data %>%
          dplyr::mutate(Chromosome = as.character(Chromosome)) %>%
          dplyr::left_join(dplyr::select(calls_maf,-c(DP_CONTROL,AF_CONTROL)),
                           by=c("Chromosome","Start_Position","VARIANT_CLASS")) %>%
          dplyr::mutate(t_depth = t_depth_estimate, t_ref_count = t_ref_count_estimate,
                        t_alt_count = t_alt_count_estimate) %>%
          dplyr::mutate(t_depth_estimate = NULL, t_ref_count_estimate = NULL, t_alt_count_estimate = NULL)

      }
      if (!any(is.na(calls_maf$DP_CONTROL)) & !any(is.na(calls_maf$AF_CONTROL))) {
        calls_maf$n_depth_estimate <- calls_maf$DP_CONTROL
        calls_maf$n_ref_count_estimate <- calls_maf$DP_CONTROL - round(calls_maf$AF_CONTROL * calls_maf$DP_CONTROL, digits = 0)
        calls_maf$n_alt_count_estimate <- calls_maf$DP_CONTROL - calls_maf$n_ref_count_estimate
        calls_maf <- calls_maf %>%
          dplyr::select(-c(DP_CONTROL, AF_CONTROL)) %>%
          dplyr::mutate(Chromosome = as.character(Chromosome))


        maf_data <- maf_data %>%
          dplyr::left_join(calls_maf, by = c("Chromosome", "Start_Position", "VARIANT_CLASS")) %>%
          dplyr::mutate(n_depth = n_depth_estimate, n_ref_count = n_ref_count_estimate, n_alt_count = n_alt_count_estimate) %>%
          dplyr::mutate(n_depth_estimate = NULL, n_ref_count_estimate = NULL, n_alt_count_estimate = NULL)

      }

      write("#version 2.4", file = maf_fname, sep = "\n")
      write.table(maf_data, file = maf_fname, col.names = T, append = T,
                  row.names = F, na = "", quote = F, sep = "\t")
      if (delete_raw == T) {
        system(paste0("rm -f ",maf_fname_tmp))
      }

    }
  }

}

targeted_drugs_pr_ttype <- function(ttype, pcgr_data,
                                    ignore_antimetabolites = T,
                                    ignore_on_label_early_phase = T,
                                    inhibitors_only = T,
                                    ignore_channel_blocker_openers = T){

  rlogging::message(paste0("Retrieving targeted drugs (on-label and off-label) for indications relevant for tumor type - ", ttype))
  assertthat::assert_that(!is.null(ttype))

  if(ttype == "Cancer, NOS"){
    ttype = "Any"
  }

  assertthat::assert_that(!is.null(pcgr_data[['antineopharma']][['drug_pr_site']][[ttype]]))

  site_candidates <- pcgr_data[['antineopharma']][['drug_pr_site']][[ttype]]

  drug_candidates <- list()
  drug_candidates[['on_label']] <-  site_candidates[['on_label']][['late_phase']]
  drug_candidates[['on_label_early_phase']] <-  site_candidates[['on_label']][['early_phase']]
  drug_candidates[['off_label']] <-  site_candidates[['off_label']]

  if(ignore_antimetabolites == T){
    if(nrow(drug_candidates[['on_label']]) > 0){
      drug_candidates[['on_label']] <- drug_candidates[['on_label']] %>%
        dplyr::filter(antimetabolite == F)
    }
    if(nrow(drug_candidates[['on_label_early_phase']]) > 0){
      drug_candidates[['on_label_early_phase']] <- drug_candidates[['on_label_early_phase']] %>%
        dplyr::filter(antimetabolite == F)
    }
    if(nrow(drug_candidates[['off_label']]) > 0){
      drug_candidates[['off_label']] <- drug_candidates[['off_label']] %>%
        dplyr::filter(antimetabolite == F)
    }
  }
  if(inhibitors_only == T){
    if(nrow(drug_candidates[['on_label']]) > 0){
      drug_candidates[['on_label']] <- drug_candidates[['on_label']] %>%
        dplyr::filter(moa_inhibition == T)
    }
    if(nrow(drug_candidates[['on_label_early_phase']]) > 0){
      drug_candidates[['on_label_early_phase']] <- drug_candidates[['on_label_early_phase']] %>%
        dplyr::filter(moa_inhibition == T)
    }
    if(nrow(drug_candidates[['off_label']]) > 0){
      drug_candidates[['off_label']] <- drug_candidates[['off_label']] %>%
        dplyr::filter(moa_inhibition == T)
    }
  }

  if(ignore_channel_blocker_openers == T){
    if(nrow(drug_candidates[['on_label']]) > 0){
      drug_candidates[['on_label']] <- drug_candidates[['on_label']] %>%
        dplyr::filter(channel_blocker_opener == F)
    }
    if(nrow(drug_candidates[['on_label_early_phase']]) > 0){
      drug_candidates[['on_label_early_phase']] <- drug_candidates[['on_label_early_phase']] %>%
        dplyr::filter(channel_blocker_opener == F)
    }
    if(nrow(drug_candidates[['off_label']]) > 0){
      drug_candidates[['off_label']] <- drug_candidates[['off_label']] %>%
        dplyr::filter(channel_blocker_opener == F)
    }
  }


  if(nrow(drug_candidates[['on_label']]) > 0){
    drug_candidates[['on_label']] <- as.data.frame(
      drug_candidates[['on_label']] %>%
        dplyr::select(symbol, nci_concept_display_name, drug_link, max_phase, drug_indication_label) %>%
        tidyr::separate_rows(drug_indication_label, sep="\\|") %>%
        dplyr::arrange(symbol, desc(max_phase)) %>%
        dplyr::distinct() %>%
        dplyr::group_by(symbol, nci_concept_display_name, drug_link) %>%
        dplyr::summarise(drug_indication_label = paste(unique(drug_indication_label), collapse=", ")) %>%
        dplyr::ungroup() %>%
        dplyr::mutate(drug_indication_label = paste0(drug_indication_label," (<b>",nci_concept_display_name,"</b>)")) %>%
        dplyr::group_by(symbol) %>%
        dplyr::summarise(DRUGS_ON_LABEL = paste0(drug_link, collapse=", "),
                         DRUGS_ON_LABEL_INDICATIONS = paste0(drug_indication_label, collapse="; "))
    )
  }
  if(nrow(drug_candidates[['on_label_early_phase']]) > 0){
    drug_candidates[['on_label_early_phase']] <- as.data.frame(
      drug_candidates[['on_label_early_phase']] %>%
        dplyr::select(symbol, nci_concept_display_name, drug_link, max_phase, drug_indication_label) %>%
        tidyr::separate_rows(drug_indication_label, sep="\\|") %>%
        dplyr::arrange(symbol, desc(max_phase)) %>%
        dplyr::distinct() %>%
        dplyr::group_by(symbol, nci_concept_display_name, drug_link) %>%
        dplyr::summarise(drug_indication_label = paste(unique(drug_indication_label), collapse=", ")) %>%
        dplyr::ungroup() %>%
        dplyr::mutate(drug_indication_label = paste0(drug_indication_label," (<b>",nci_concept_display_name,"</b>)")) %>%
        dplyr::group_by(symbol) %>%
        dplyr::summarise(DRUGS_ON_LABEL_EARLY_PHASE = paste0(drug_link, collapse=", "),
                         DRUGS_ON_LABEL_EARLY_PHASE_INDICATIONS = paste0(drug_indication_label, collapse="; "))
    )
  }

  if(nrow(drug_candidates[['off_label']]) > 0){
    drug_candidates[['off_label']] <- as.data.frame(
      drug_candidates[['off_label']] %>%
        dplyr::select(symbol, nci_concept_display_name, drug_link, max_phase, drug_indication_label) %>%
        tidyr::separate_rows(drug_indication_label, sep="\\|") %>%
        dplyr::arrange(symbol, desc(max_phase), nci_concept_display_name) %>%
        dplyr::distinct() %>%
        dplyr::group_by(symbol, nci_concept_display_name, drug_link) %>%
        dplyr::summarise(drug_indication_label = paste(unique(drug_indication_label), collapse=", ")) %>%
        dplyr::ungroup() %>%
        dplyr::mutate(drug_indication_label = paste0(drug_indication_label," (<b>",nci_concept_display_name,"</b>)")) %>%
        dplyr::group_by(symbol) %>%
        dplyr::summarise(DRUGS_OFF_LABEL = paste0(drug_link, collapse=", "),
                         DRUGS_OFF_LABEL_INDICATIONS = paste0(drug_indication_label, collapse="; "))
    )
  }

  all_candidates <- dplyr::full_join(drug_candidates[['off_label']],
                                     drug_candidates[['on_label_early_phase']],
                                     by = "symbol")
  if(nrow(drug_candidates[['on_label']]) > 0){
    all_candidates <- dplyr::full_join(all_candidates,
                                       drug_candidates[['on_label']],
                                       by = "symbol")
  }else{
    all_candidates$DRUGS_ON_LABEL <- NA
    all_candidates$DRUGS_ON_LABEL_INDICATIONS <- NA
  }

  all_drug_target_candidates <- all_candidates %>%
    dplyr::rename(SYMBOL = symbol) %>%
    dplyr::arrange(SYMBOL)

  if(ignore_on_label_early_phase == T){
    all_drug_target_candidates <- all_drug_target_candidates %>%
      dplyr::select(-c(DRUGS_ON_LABEL_EARLY_PHASE,
                       DRUGS_ON_LABEL_EARLY_PHASE_INDICATIONS)) %>%
      dplyr::filter(!is.na(DRUGS_ON_LABEL) | !is.na(DRUGS_OFF_LABEL))
  }




  return(all_drug_target_candidates)


}


