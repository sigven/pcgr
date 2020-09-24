
#' Function that removes copy number segments that go beyond chrosomal lengths for the given assembly
#'
#' @param cna_segments_df genomic segments with CNA
#' @param genome_assembly genome assembly (grch37/grch38)
#' @param bsg BSgenome object
#'
get_valid_chromosome_segments <- function(cna_segments_df,
                                          genome_assembly = "grch37",
                                          bsg = BSgenome.Hsapiens.UCSC.hg19) {

  invisible(assertthat::assert_that(!is.null(cna_segments_df) & is.data.frame(cna_segments_df),
                                    msg = "Argument 'cna_segments_df' must be a valid data.frame() object"))
  assertable::assert_colnames(cna_segments_df, c("chromosome", "segment_end"), only_colnames = F, quiet = T)
  assertable::assert_coltypes(cna_segments_df, list(chromosome = character(), segment_end = integer()), quiet = T)

  chromosome_lengths <- data.frame(chromosome = head(names(seqlengths(bsg)), 24),
                                   chrom_length = head(seqlengths(bsg), 24),
                                   stringsAsFactors = F, row.names = NULL)
  cna_segments_df <- as.data.frame(
    cna_segments_df %>%
      dplyr::left_join(chromosome_lengths, by = c("chromosome")) %>%
      dplyr::rowwise() %>%
      dplyr::mutate(segment_error = segment_end > chrom_length)
    )
  if (nrow(dplyr::filter(cna_segments_df, segment_error == T)) > 0) {
    n_removed <- nrow(dplyr::filter(cna_segments_df, segment_error == T))
    rlogging::warning("Skipping ", n_removed, " copy number segments that span beyond chromosomal lengths for ",
                      genome_assembly, " (make sure chromosomal segments are consistent with assembly)")
  }
  cna_segments_df <- cna_segments_df %>%
    dplyr::filter(segment_error == F) %>%
    dplyr::select(-c(segment_error, chrom_length))

  return(cna_segments_df)
}

#' Function that gets the chromosome bands of copy number segments
#'
#' @param cna_df genomic ranges object with copy number aberrations
#' @param pcgr_data pcgr data bundle object
#'
get_cna_cytoband <- function(cna_df, pcgr_data = NULL) {

  cna_gr <-
    GenomicRanges::makeGRangesFromDataFrame(cna_df, keep.extra.columns = T,
                                            seqinfo = pcgr_data[["assembly"]][["seqinfo"]],
                                            seqnames.field = "chromosome",
                                            start.field = "segment_start",
                                            end.field = "segment_end",
                                            ignore.strand = T, starts.in.df.are.0based = T)

  cytoband_gr <- pcgr_data[["genomic_ranges"]][["cytoband"]]

  invisible(assertthat::assert_that("focalCNAthreshold" %in% names(mcols(cytoband_gr))))
  cyto_hits <- GenomicRanges::findOverlaps(cna_gr, cytoband_gr, type = "any", select = "all")
  ranges <- cytoband_gr[S4Vectors::subjectHits(cyto_hits)]
  mcols(ranges) <- c(mcols(ranges), mcols(cna_gr[S4Vectors::queryHits(cyto_hits)]))
  cyto_df <- as.data.frame(mcols(ranges)) %>%
    dplyr::mutate(segment_start = start(ranges(cna_gr[S4Vectors::queryHits(cyto_hits)]))) %>%
    dplyr::mutate(segment_end = end(ranges(cna_gr[S4Vectors::queryHits(cyto_hits)]))) %>%
    dplyr::mutate(segment_length = width(ranges(cna_gr[S4Vectors::queryHits(cyto_hits)])))

  cyto_stats <- as.data.frame(
    cyto_df %>%
      dplyr::group_by(SEGMENT_ID, segment_length) %>%
      dplyr::summarise(CYTOBAND = paste(name, collapse = ", "),
                       chromosome_arm = paste(unique(arm), collapse = ","),
                       focalCNAthresholds = paste(unique(focalCNAthreshold), collapse = ","),
                       .groups = "drop") %>%
      dplyr::ungroup() %>%
      dplyr::mutate(focalCNAthresholds = dplyr::if_else(stringr::str_detect(focalCNAthresholds, ","),
                                                        as.character(NA),
                                                        as.character(focalCNAthresholds))) %>%
      dplyr::mutate(focalCNAthresholds = dplyr::if_else(!is.na(focalCNAthresholds),
                                                        as.numeric(focalCNAthresholds),
                                                        as.numeric(NA))) %>%
      dplyr::rowwise() %>%
      dplyr::mutate(broad_cnv_event = segment_length > focalCNAthresholds) %>%
      dplyr::mutate(EVENT_TYPE = "broad") %>%
      dplyr::mutate(EVENT_TYPE = dplyr::if_else(!is.na(broad_cnv_event) & broad_cnv_event == F,
                                                "focal", EVENT_TYPE)) %>%
      dplyr::mutate(EVENT_TYPE = dplyr::if_else(is.na(broad_cnv_event),
                                                as.character("broad"), EVENT_TYPE)) %>%
      dplyr::mutate(CYTOBAND = stringr::str_replace(CYTOBAND, ", (\\S+, ){0,}", " - ")) %>%
      tidyr::separate(SEGMENT_ID, sep = ":", into = c("chrom", "start_stop"), remove = F) %>%
      dplyr::mutate(CYTOBAND = paste0(chrom, ":", CYTOBAND)) %>%
      dplyr::select(SEGMENT_ID, CYTOBAND, EVENT_TYPE)

  )
  cna_df <- cna_df %>%
    dplyr::left_join(cyto_stats, by = "SEGMENT_ID")

  return(cna_df)

}

#' Function that annotates CNV segment files
#'
#' @param cna_segments_tsv CNV file name with chromosomal log(2)-ratio segments
#' @param pcgr_data object with PCGR annotation data
#' @param sample_name sample identifier
#' @param pcgr_config Object with PCGR configuration parameters
#' @param transcript_overlap_pct required aberration overlap fraction (percent) for reported transcripts (default 100 percent)
#'
generate_report_data_cna <- function(cna_segments_tsv, pcgr_data, sample_name, pcgr_config,
                                     transcript_overlap_pct = 100) {

  invisible(assertthat::assert_that(file.exists(cna_segments_tsv),
                                    msg = paste0("File 'cna_segments_tsv' (", cna_segments_tsv, ") does not exist")))
  pcg_report_cna <- pcgrr::init_report(pcgr_config, sample_name, class = "cna")
  logR_homdel <- pcgr_config[["cna"]][["logR_homdel"]]
  logR_gain <- pcgr_config[["cna"]][["logR_gain"]]
  tumor_type <- pcgr_config[['tumor_properties']][['tumor_type']]
  MEGABASE <- 1000000

  rlogging::message("------")
  rlogging::message(paste0("Generating report data for copy number segment file ", cna_segments_tsv))

  ## READ INPUT FILE, VALIDATE INPUT CHROMOSOMES AND SEGMENTS, ADD CYTOBAND INFO
  cna_df <- read.table(file = cna_segments_tsv, header = T,
                           stringsAsFactors = F, sep = "\t",
                           comment.char = "", quote = "") %>%
    dplyr::rename(chromosome = Chromosome, LogR = Segment_Mean,
                  segment_start = Start, segment_end = End) %>%
    dplyr::distinct() %>%
    dplyr::select(chromosome, LogR, segment_start, segment_end) %>%
    dplyr::mutate(chromosome = stringr::str_replace(chromosome, "^chr", "")) %>%
    pcgrr::get_valid_chromosomes(chromosome_column = "chromosome",
                                 bsg = pcgr_data[["assembly"]][["bsg"]]) %>%
    pcgrr::get_valid_chromosome_segments(genome_assembly = pcgr_data[["assembly"]][["grch_name"]],
                                         bsg = pcgr_data[["assembly"]][["bsg"]]) %>%
    dplyr::filter(!is.na(LogR)) %>%
    dplyr::mutate(LogR = round(as.numeric(LogR), digits = 3)) %>%
    dplyr::mutate(SEGMENT_ID = paste0(chromosome, ":", segment_start, "-", segment_end)) %>%
    pcgrr::get_cna_cytoband(pcgr_data = pcgr_data) %>%
    dplyr::mutate(SAMPLE_ID = sample_name) %>%
    pcgrr::append_ucsc_segment_link(hgname = pcgr_data[["assembly"]][["hg_name"]],
                                    chrom = "chromosome",
                                    start = "segment_start",
                                    end = "segment_end") %>%
    dplyr::mutate(SEGMENT_LENGTH_MB = round((as.numeric((segment_end - segment_start) / MEGABASE)),
                                            digits = 5)) %>%
    dplyr::rename(SEGMENT = SEGMENT_LINK, LOG_R = LogR)

  ## MAKE SIMPLE SEGMENTS DATA FRAME FOR FILTERING IN REPORT
  cna_segments <- cna_df %>%
    dplyr::select(SEGMENT, SEGMENT_LENGTH_MB, CYTOBAND, LOG_R, EVENT_TYPE) %>%
    dplyr::distinct()

  #### FIND AND APPEND GENCODE TRANSCRIPTS THAT OVERLAP
  cna_transcript_df <- pcgrr::get_cna_overlapping_transcripts(cna_df, pcgr_data = pcgr_data)

  #### GENERATE DATAFRAME OF UNIQUE TRANSCRIPT-CNA SEGMENTS FOR OUTPUT TSV
  cna_transcript_df_print <- cna_transcript_df %>%
    dplyr::select(chrom, segment_start, segment_end, SEGMENT_ID, SEGMENT_LENGTH_MB,
                  EVENT_TYPE, CYTOBAND, LOG_R, SAMPLE_ID, ensembl_gene_id,
                  symbol, ensembl_transcript_id, transcript_start,
                  transcript_end, transcript_overlap_percent, name, biotype,
                  tumor_suppressor, oncogene, intogen_drivers, chembl_compound_id,
                  gencode_tag, gencode_release) %>%
    magrittr::set_colnames(tolower(names(.)))

  avg_transcript_overlap <- as.data.frame(
    cna_transcript_df %>%
      dplyr::filter(biotype == "protein_coding") %>%
      dplyr::group_by(SEGMENT_ID, symbol) %>%
      dplyr::summarise(MEAN_TRANSCRIPT_CNA_OVERLAP = mean(transcript_overlap_percent),
                       TRANSCRIPTS = paste0(ensembl_transcript_id, collapse = ", ")) %>%
      dplyr::rename(SYMBOL = symbol) %>%
      dplyr::mutate(MEAN_TRANSCRIPT_CNA_OVERLAP = round(MEAN_TRANSCRIPT_CNA_OVERLAP, digits = 2))
  )

  cna_transcript_df <- dplyr::select(cna_transcript_df, -ensembl_transcript_id) %>%
    dplyr::filter(biotype == "protein_coding") %>%
    dplyr::distinct() %>%
    dplyr::mutate(VAR_ID = as.character(rep(1:nrow(.)))) %>%
    magrittr::set_colnames(toupper(names(.))) %>%
    dplyr::select(VAR_ID, SEGMENT_ID, SYMBOL, ONCOGENE,
                  ONCOGENE_EVIDENCE, TUMOR_SUPPRESSOR,
                  TUMOR_SUPPRESSOR_EVIDENCE, CANCERGENE_SUPPORT,
                  ENTREZGENE, CHROM, NAME, EVENT_TYPE,
                  SEGMENT_LENGTH_MB, SEGMENT,
                  TRANSCRIPT_OVERLAP_PERCENT, LOG_R) %>%
    dplyr::mutate(ENTREZ_ID = as.character(ENTREZGENE)) %>%
    dplyr::rename(GENENAME = NAME,
                  TRANSCRIPT_OVERLAP = TRANSCRIPT_OVERLAP_PERCENT,
                  CHROMOSOME = CHROM) %>%
    dplyr::left_join(pcgr_data[["kegg"]][["pathway_links"]], by = c("ENTREZ_ID" = "gene_id")) %>%
    dplyr::rename(KEGG_PATHWAY = kegg_pathway_urls)

  ## Get gene annotation links
  entrezgene_annotation_links <-
    pcgrr::generate_annotation_link(cna_transcript_df,
                                    vardb = "GENE_NAME",
                                    group_by_var = "VAR_ID",
                                    link_key_var = "ENTREZ_ID",
                                    link_display_var = "GENENAME",
                                    url_prefix = "http://www.ncbi.nlm.nih.gov/gene/")

  cna_transcript_df <- cna_transcript_df %>%
    dplyr::left_join(dplyr::rename(entrezgene_annotation_links, GENE_NAME = link),
                     by = c("VAR_ID")) %>%
    dplyr::select(SEGMENT_ID, CHROMOSOME, SYMBOL, GENE_NAME, KEGG_PATHWAY,
                  TUMOR_SUPPRESSOR, TUMOR_SUPPRESSOR_EVIDENCE, ONCOGENE,
                  ONCOGENE_EVIDENCE, CANCERGENE_SUPPORT, SEGMENT_LENGTH_MB,
                  SEGMENT, EVENT_TYPE, LOG_R) %>%
    dplyr::distinct() %>%
    dplyr::left_join(avg_transcript_overlap, by = c("SEGMENT_ID", "SYMBOL"))


  n_cna_loss <- dplyr::filter(cna_segments, LOG_R <= logR_homdel) %>% nrow()
  n_cna_gain <- dplyr::filter(cna_segments, LOG_R >= logR_gain) %>% nrow()
  cna_segments_filtered <- cna_segments %>%
    dplyr::filter(LOG_R >= logR_gain | LOG_R <= logR_homdel) %>%
    dplyr::arrange(desc(LOG_R))
  rlogging::message(paste0("Detected ", nrow(cna_segments_filtered),
                           " segments subject to amplification/deletion (",
                           n_cna_loss, " deletions, ", n_cna_gain,
                           " gains according to user-defined log(2) ratio thresholds)"))


  ## Get aberration sets related to tumor suppressor genes/oncogenes/drug targets
  onco_ts_sets <- pcgrr::get_oncogene_tsgene_target_sets(cna_transcript_df,
                                                  logR_homdel = logR_homdel,
                                                  logR_gain = logR_gain,
                                                  tumor_type = tumor_type,
                                                  pcgr_data = pcgr_data)

  ## Get all clinical evidence items that are related to
  ## tumor suppressor genes/oncogenes/drug targets (NOT tumor-type specific)
  biomarker_hits_cna_any <-
    pcgrr::get_clin_assocs_cna(onco_ts_sets, pcgr_data, tumor_type,
                                         tumor_type_specificity = "any")
  pcg_report_cna[["clin_eitem"]][["any_ttype"]] <-
    biomarker_hits_cna_any[["clin_eitem"]]
  pcg_report_cna[["variant_set"]][["tier2"]] <-
    biomarker_hits_cna_any$variant_set

  ## Get all clinical evidence items that overlap query set (if tumor type is specified)
  if(tumor_type != "Cancer, NOS"){
    biomarker_hits_cna_specific <-
      pcgrr::get_clin_assocs_cna(onco_ts_sets, pcgr_data,
                                tumor_type, tumor_type_specificity = "specific")

    ## Assign putative TIER 1 variant set
    pcg_report_cna[["clin_eitem"]][["specific_ttype"]] <-
      biomarker_hits_cna_specific$clin_eitem
    pcg_report_cna[["variant_set"]][["tier1"]] <-
      biomarker_hits_cna_specific$variant_set
  }

  pcg_report_cna[["eval"]] <- T
  pcg_report_cna[["variant_set"]][["tsv"]] <- cna_transcript_df_print
  pcg_report_cna[["variant_statistic"]][["n_cna_gain"]] <- n_cna_gain
  pcg_report_cna[["variant_statistic"]][["n_cna_loss"]] <- n_cna_loss
  pcg_report_cna[["variant_display"]][["segment"]] <- cna_segments_filtered
  pcg_report_cna[["variant_display"]][["oncogene_gain"]] <- onco_ts_sets[["oncogene_gain"]]
  pcg_report_cna[["variant_display"]][["tsgene_loss"]] <- onco_ts_sets[["tsgene_loss"]]
  pcg_report_cna[["variant_display"]][["other_target"]] <- onco_ts_sets[["other_target"]]


  pcg_report_cna <- pcgrr::assign_tier1_tier2_acmg_cna(pcg_report_cna)
  #pcg_report_cna[['clin_eitem']][['other_ttype']] <- tier_1_2_biomarkers[['clin_eitem']][['other_ttype']]
  #pcg_report_cna[["variant_display"]][["tier1"]] <- tier_1_2_biomarkers[['tier1']]
  #pcg_report_cna[["variant_display"]][["tier2"]] <- tier_1_2_biomarkers[['tier2']]

  return(pcg_report_cna)
}

get_oncogene_tsgene_target_sets <- function(cna_df, transcript_overlap_pct = 100,
                                     logR_gain = 0.8, logR_homdel = -0.8, tumor_type = "Any",
                                     pcgr_data = NULL){

  invisible(assertthat::assert_that(!is.null(pcgr_data)))
  invisible(assertthat::assert_that(is.data.frame(cna_df), msg = "Argument 'cna_df' must be of type data.frame"))
  assertable::assert_colnames(cna_df, c("ONCOGENE", "TUMOR_SUPPRESSOR", "MEAN_TRANSCRIPT_CNA_OVERLAP", "LOG_R",
                                        "SYMBOL", "KEGG_PATHWAY", "SEGMENT_ID", "CHROMOSOME",
                                        "GENE_NAME","EVENT_TYPE"), only_colnames = F, quiet = T)

  onco_ts_sets <- list()
  onco_ts_sets[["oncogene_gain"]] <- data.frame()
  onco_ts_sets[["oncogene_gain"]] <- dplyr::filter(cna_df, ONCOGENE == T & TUMOR_SUPPRESSOR == F &
                                                     MEAN_TRANSCRIPT_CNA_OVERLAP >= transcript_overlap_pct &
                                                     LOG_R >= logR_gain)
  onco_ts_sets[["tsgene_loss"]] <- data.frame()
  onco_ts_sets[["tsgene_loss"]] <- dplyr::filter(cna_df, TUMOR_SUPPRESSOR == T &
                                                   MEAN_TRANSCRIPT_CNA_OVERLAP >= transcript_overlap_pct  &
                                                   LOG_R <= logR_homdel)

  onco_ts_sets[["other_target"]] <- data.frame()
  onco_ts_sets[["other_target"]] <- dplyr::filter(cna_df, TUMOR_SUPPRESSOR == F & ONCOGENE == F &
                                                   MEAN_TRANSCRIPT_CNA_OVERLAP >= transcript_overlap_pct  &
                                                   LOG_R >= logR_gain)

  drug_target_site <- pcgrr::targeted_drugs_pr_ttype(tumor_type, pcgr_data,
                                                     ignore_antimetabolites = T,
                                                     inhibitors_only = T,
                                                     ignore_on_label_early_phase = T,
                                                     ignore_channel_blocker_openers = T)

  if(nrow(onco_ts_sets[['other_target']]) > 0){
    onco_ts_sets[['other_target']] <- onco_ts_sets[['other_target']] %>%
      dplyr::left_join(drug_target_site, by = "SYMBOL") %>%
      dplyr::filter(!is.na(DRUGS_ON_LABEL) | !is.na(DRUGS_OFF_LABEL)) %>%
      dplyr::select(-c(TUMOR_SUPPRESSOR, ONCOGENE,
                       TUMOR_SUPPRESSOR_EVIDENCE, ONCOGENE_EVIDENCE)) %>%
      dplyr::distinct() %>%
      dplyr::select(CHROMOSOME, SYMBOL, GENE_NAME, SEGMENT, EVENT_TYPE,
                    DRUGS_ON_LABEL, DRUGS_OFF_LABEL,
                    SEGMENT_LENGTH_MB, LOG_R, SEGMENT_ID,
                    DRUGS_ON_LABEL_INDICATIONS,
                    DRUGS_OFF_LABEL_INDICATIONS, MEAN_TRANSCRIPT_CNA_OVERLAP,
                    KEGG_PATHWAY, TRANSCRIPTS) %>%
      dplyr::mutate(MEAN_TRANSCRIPT_CNA_OVERLAP = paste0(MEAN_TRANSCRIPT_CNA_OVERLAP, "%")) %>%
      dplyr::mutate(CNA_TYPE = "gain")

    rlogging::message(paste0("Detected ", nrow(onco_ts_sets[['other_target']]),
                             " drug targets to amplification (log(2) ratio >= ",
                             logR_gain, "): ", paste0(unique(onco_ts_sets[['other_target']]$SYMBOL),
                                                      collapse = ", ")))

  }

  for (t in c("oncogene_gain", "tsgene_loss")) {
    if (nrow(onco_ts_sets[[t]]) > 0) {
      onco_ts_sets[[t]] <- onco_ts_sets[[t]] %>%
        dplyr::select(-c(TUMOR_SUPPRESSOR, ONCOGENE)) %>%
        dplyr::distinct() %>%
        dplyr::select(CHROMOSOME, SYMBOL, GENE_NAME, SEGMENT, LOG_R, EVENT_TYPE,
                      SEGMENT_LENGTH_MB,  MEAN_TRANSCRIPT_CNA_OVERLAP,
                      KEGG_PATHWAY, TRANSCRIPTS, SEGMENT_ID) %>%
        dplyr::mutate(MEAN_TRANSCRIPT_CNA_OVERLAP = paste0(MEAN_TRANSCRIPT_CNA_OVERLAP, "%")) %>%
        dplyr::mutate(CNA_TYPE = dplyr::if_else(t == "oncogene_gain" | t == "other_target", "gain", "loss"))
      if (t == "oncogene_gain") {
        rlogging::message(paste0("Detected ", nrow(onco_ts_sets[[t]]),
                                 " proto-oncogene(s) subject to amplification (log(2) ratio >= ",
                                 logR_gain, "): ", paste0(unique(onco_ts_sets[[t]]$SYMBOL), collapse = ", ")))
        onco_ts_sets[[t]] <- onco_ts_sets[[t]] %>%
          dplyr::left_join(drug_target_site, by = "SYMBOL") %>%
          dplyr::select(CHROMOSOME, SYMBOL, GENE_NAME, SEGMENT, LOG_R, EVENT_TYPE,
                        DRUGS_ON_LABEL, DRUGS_OFF_LABEL, dplyr::everything()) %>%
          dplyr::arrange(DRUGS_ON_LABEL, DRUGS_OFF_LABEL)
      }else{
          rlogging::message(paste0("Detected ", nrow(onco_ts_sets[[t]]),
                                 " tumor suppressor gene(s) subject to homozygous deletions (log(2) ratio <= ",
                                 logR_homdel, "): ", paste0(unique(onco_ts_sets[[t]]$SYMBOL), collapse = ", ")))
      }
    }else{
      if (t == "tsgene_loss") {
        rlogging::message(paste0("Detected 0 tumor suppressor genes subject to homozygous deletion (log(2) ratio <= ", logR_homdel))
      }else{
        if(t == "other_target"){
          rlogging::message(paste0("Detected 0 other drug targets subject to amplification (log(2) ratio >= ", logR_gain))
        }else{
          rlogging::message(paste0("Detected 0 proto-oncogenes subject to amplification (log(2) ratio >= ", logR_gain))
        }
      }
    }
  }

  return(onco_ts_sets)

}

get_cna_overlapping_transcripts <- function(cna_df, pcgr_data){

  cna_gr <-
    GenomicRanges::makeGRangesFromDataFrame(cna_df, keep.extra.columns = T, seqinfo = pcgr_data[["assembly"]][["seqinfo"]],
                                            seqnames.field = "chromosome", start.field = "segment_start", end.field = "segment_end",
                                            ignore.strand = T, starts.in.df.are.0based = T)

  hits <- GenomicRanges::findOverlaps(cna_gr, pcgr_data[["genomic_ranges"]][["gencode_genes"]],
                                      type = "any", select = "all")
  ranges <- pcgr_data[["genomic_ranges"]][["gencode_genes"]][S4Vectors::subjectHits(hits)]
  mcols(ranges) <- c(mcols(ranges), mcols(cna_gr[S4Vectors::queryHits(hits)]))

  cna_transcript_df <-
    as.data.frame(
      as.data.frame(mcols(ranges)) %>%
        dplyr::mutate(segment_start = as.integer(start(ranges(cna_gr[S4Vectors::queryHits(hits)])))) %>%
        dplyr::mutate(segment_end = as.integer(end(ranges(cna_gr[S4Vectors::queryHits(hits)])))) %>%
        dplyr::mutate(transcript_start = start(ranges)) %>%
        dplyr::mutate(transcript_end = end(ranges)) %>%
        dplyr::mutate(chrom = as.character(seqnames(ranges))) %>%
        dplyr::rowwise() %>%
        dplyr::mutate(transcript_overlap_percent =
                        round(as.numeric((min(transcript_end, segment_end) - max(segment_start, transcript_start)) /
                                           (transcript_end - transcript_start)) * 100, digits = 2))
    ) %>%
    pcgrr::sort_chromosomal_segments(chromosome_column = "chrom", start_segment = "segment_start",
                                     end_segment = "segment_end")

  return(cna_transcript_df)
}



