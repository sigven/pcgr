# .library(BSgenome.Hsapiens.UCSC.hg19)
# library(BSgenome.Hsapiens.UCSC.hg38)
# library(magrittr)
# library(deconstructSigs)
# library(data.table)
# library(ggplot2)
# library(plotly)
# library(rlogging)
# library(configr)
# library(RcppTOML)

#' Function that plots a histogram of the the variant allelic support (tumor) - grouped by tiers
#'
#' @param tier_df data frame with somatic mutations
#' @param bin_size size of bins for allelic frequency
#' @return p geom_histogram plot from ggplot2
#'
#'
tier_af_distribution <- function(tier_df, bin_size = 0.1){
  af_bin_df <- data.frame()
  i <- 1
  num_bins <- as.integer(1 / bin_size)
  bin_start <- 0
  while (i <= num_bins){
    bin_end <- bin_start + bin_size
    bin_name <- as.character(paste0(bin_start, " - ", bin_end))
    j <- 1
    while (j <= 4){
      TIER <- paste0("TIER ", j)
      df <- data.frame(bin_name = bin_name, bin_start = bin_start, bin_end = bin_end, bin = as.integer(i), TIER = TIER, stringsAsFactors = F)
      af_bin_df <- rbind(af_bin_df, df)
      j <- j + 1
    }
    TIER <- "NONCODING"
    df <- data.frame(bin_name = bin_name, bin_start = bin_start, bin_end = bin_end, bin = as.integer(i), TIER = TIER, stringsAsFactors = F)
    af_bin_df <- rbind(af_bin_df, df)
    bin_start <- bin_end
    i <- i + 1
  }

  tier_df_trans <- transform(tier_df, bin = cut(AF_TUMOR, breaks = seq(from = 0, to = 1, by = bin_size), right = F, include.lowest = T, labels = F))
  tier_df_trans_bin <- as.data.frame(dplyr::group_by(tier_df_trans, TIER, bin) %>% dplyr::summarise(Count = n()))

  af_bin_df <- dplyr::left_join(af_bin_df, tier_df_trans_bin, by = c("bin", "TIER"))
  if (nrow(af_bin_df[is.na(af_bin_df$Count), ]) > 0){
    af_bin_df[is.na(af_bin_df$Count), ]$Count <- 0
  }

  p <- ggplot2::ggplot(data = af_bin_df) + ggplot2::geom_bar(mapping = ggplot2::aes(x = bin_name, y = Count, fill = TIER), stat = "identity") +
    ggplot2::scale_color_brewer(palette = "Dark2") +
    ggplot2::scale_fill_brewer(palette = "Dark2") +
    ggplot2::theme_classic() +
    ggplot2::ylab("Number of variants") +
    ggplot2::xlab("Variant allelic fraction - tumor") +
    ggplot2::theme(legend.title = ggplot2::element_blank(),
                   #axis.text.x = element_blank(),
                   axis.text.x = ggplot2::element_text(angle = 45, family = "Helvetica", size = 12, vjust = -0.1),
                   axis.title.x = ggplot2::element_text(family = "Helvetica", size = 12, vjust = -2),
                   axis.text.y = ggplot2::element_text(family = "Helvetica", size = 12),
                   axis.title.y = ggplot2::element_text(family = "Helvetica", size = 12, vjust = 1.5),
                   plot.margin = (grid::unit(c(0.5, 2, 2, 0.5), "cm")),
                   legend.text = ggplot2::element_text(family = "Helvetica", size = 12))

  return(p)

}

#' Function that check for valid mutations/chromosomes in input
#'
#' @param vcf_data_df data frame
#' @param chromosome_column name of columns for which string replace is to be performed
#' @param pattern pattern to replace
#' @return vcf_data_df_valid data frame with valid mutations
#'
get_valid_chromosomes <- function(vcf_data_df, chromosome_column = "CHROM", bsg = BSgenome.Hsapiens.UCSC.hg19){
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
df_string_replace <- function(df, strings, pattern, replacement, replace_all = F){

  for (column_name in strings){
    if (column_name %in% colnames(df)){
      if (replace_all == F){
        df[, column_name] <- stringr::str_replace(df[, column_name], pattern = pattern, replacement = replacement)
      }else{
        df[, column_name] <- stringr::str_replace_all(df[, column_name], pattern = pattern, replacement = replacement)
      }
    }
  }
  return(df)
}

#' Function that transforms a tier-structured variant data frame into a MAF-like data frame (for input to 2020plus, MutSigCV)
#'
#' @param tier_df data frame with somatic mutations
#' @return maf_df
#'
#'

tier_to_maf <- function(tier_df){
  maf_df <- dplyr::select(tier_df, SYMBOL, GENOMIC_CHANGE, VCF_SAMPLE_ID, CONSEQUENCE, VARIANT_CLASS, PROTEIN_CHANGE)
  maf_df$Hugo_Symbol <- maf_df$SYMBOL
  locus_info <- tidyr::separate(dplyr::select(maf_df, GENOMIC_CHANGE), GENOMIC_CHANGE, c("Chromosome", "pos_alleles"), sep = ":", convert = T)
  #maf_df$Chromosome <- stringr::str_replace(locus_info$chrom, pattern = "chr", replacement = '')
  maf_df$Chromosome <- locus_info$Chromosome
  locus_info$pos_alleles <- stringr::str_replace(locus_info$pos_alleles, pattern = "g\\.", replacement = "")
  maf_df$Reference_Allele <- stringr::str_split_fixed(stringr::str_replace_all(locus_info$pos_alleles, "[0-9]{1,}", ""), pattern = ">", n = 2)[, 1]
  maf_df$Tumor_Seq_Allele2 <- stringr::str_split_fixed(stringr::str_replace_all(locus_info$pos_alleles, "[0-9]{1,}", ""), pattern = ">", n = 2)[, 2]
  maf_df$Start_Position <- as.integer(stringr::str_replace_all(locus_info$pos_alleles, "[A-Z]{1,}>[A-Z]{1,}", ""))
  maf_df$End_Position <- maf_df$Start_Position + nchar(maf_df$Reference_Allele) - 1
  maf_df$Tumor_Sample_Barcode <- maf_df$VCF_SAMPLE_ID
  maf_df$NCBI_Build <- "GRCh37"
  maf_df$Amino_Acid_Change <- maf_df$PROTEIN_CHANGE

  maf_df$Variant_Classification <- character(nrow(maf_df))
  maf_df <- dplyr::filter(maf_df, !is.na(CONSEQUENCE))
  if(nrow(maf_df[stringr::str_detect(maf_df$CONSEQUENCE, "^(splice_acceptor_variant|splice_donor_variant|transcript_ablation|exon_loss_variant)"), ]) > 0){
    maf_df[stringr::str_detect(maf_df$CONSEQUENCE, "^(splice_acceptor_variant|splice_donor_variant|transcript_ablation|exon_loss_variant)"), ]$Variant_Classification <- "Splice_Site"
  }
  if (nrow(maf_df[stringr::str_detect(maf_df$CONSEQUENCE, "stop_gained"), ]) > 0){
    maf_df[stringr::str_detect(maf_df$CONSEQUENCE, "stop_gained"), ]$Variant_Classification <- "Nonsense_Mutation"
  }
  if (nrow(maf_df[stringr::str_detect(maf_df$CONSEQUENCE, "frameshift_variant") & maf_df$VARIANT_CLASS == "deletion", ]) > 0){
    maf_df[stringr::str_detect(maf_df$CONSEQUENCE, "frameshift_variant") & maf_df$VARIANT_CLASS == "deletion", ]$Variant_Classification <- "Frame_Shift_Del"
  }
  if (nrow(maf_df[stringr::str_detect(maf_df$CONSEQUENCE, "frameshift_variant") & maf_df$VARIANT_CLASS == "insertion", ]) > 0){
    maf_df[stringr::str_detect(maf_df$CONSEQUENCE, "frameshift_variant") & maf_df$VARIANT_CLASS == "insertion", ]$Variant_Classification <- "Frame_Shift_Ins"
  }
  if (nrow(maf_df[stringr::str_detect(maf_df$CONSEQUENCE, "stop_lost"), ]) > 0){
    maf_df[stringr::str_detect(maf_df$CONSEQUENCE, "stop_lost"), ]$Variant_Classification <- "Nonstop_Mutation"
  }
  if (nrow(maf_df[stringr::str_detect(maf_df$CONSEQUENCE, "initiator_codon_variant|start_lost"), ])){
    maf_df[stringr::str_detect(maf_df$CONSEQUENCE, "initiator_codon_variant|start_lost"), ]$Variant_Classification <- "Translation_Start_Site"
  }
  if (nrow(maf_df[stringr::str_detect(maf_df$CONSEQUENCE, "inframe_insertion") & maf_df$VARIANT_CLASS == "insertion", ]) > 0){
    maf_df[stringr::str_detect(maf_df$CONSEQUENCE, "inframe_insertion") & maf_df$VARIANT_CLASS == "insertion", ]$Variant_Classification <- "In_Frame_Ins"
  }
  if (nrow(maf_df[stringr::str_detect(maf_df$CONSEQUENCE, "inframe_deletion") & maf_df$VARIANT_CLASS == "deletion", ]) > 0){
    maf_df[stringr::str_detect(maf_df$CONSEQUENCE, "inframe_deletion") & maf_df$VARIANT_CLASS == "deletion", ]$Variant_Classification <- "In_Frame_Del"
  }
  if (nrow(maf_df[stringr::str_detect(maf_df$CONSEQUENCE, "^(missense_variant|coding_sequence_variant|conservative_missense_variant|rare_amino_acid_variant)"), ]) > 0){
    maf_df[stringr::str_detect(maf_df$CONSEQUENCE, "^(missense_variant|coding_sequence_variant|conservative_missense_variant|rare_amino_acid_variant)"), ]$Variant_Classification <- "Missense_Mutation"
  }
  if (nrow(maf_df[stringr::str_detect(maf_df$CONSEQUENCE, "3_prime_UTR_variant"), ]) > 0){
    maf_df[stringr::str_detect(maf_df$CONSEQUENCE, "3_prime_UTR_variant"), ]$Variant_Classification <- "3'UTR"
  }
  if (nrow(maf_df[stringr::str_detect(maf_df$CONSEQUENCE, "5_prime_UTR_variant"), ]) > 0){
    maf_df[stringr::str_detect(maf_df$CONSEQUENCE, "5_prime_UTR_variant"), ]$Variant_Classification <- "5'UTR"
  }
  if (nrow(maf_df[stringr::str_detect(maf_df$CONSEQUENCE, "^(TF_binding_site_variant|regulatory_region_variant|regulatory_region|intergenic_variant|intergenic_region)"), ]) > 0){
    maf_df[stringr::str_detect(maf_df$CONSEQUENCE, "^(TF_binding_site_variant|regulatory_region_variant|regulatory_region|intergenic_variant|intergenic_region)"), ]$Variant_Classification <- "IGR"
  }
  if (nrow(maf_df[stringr::str_detect(maf_df$CONSEQUENCE, "^(incomplete_terminal_codon_variant|synonymous_variant|stop_retained_variant|NMD_transcript_variant)"), ]) > 0){
    maf_df[stringr::str_detect(maf_df$CONSEQUENCE, "^(incomplete_terminal_codon_variant|synonymous_variant|stop_retained_variant|NMD_transcript_variant)"), ]$Variant_Classification <- "Silent"
  }
  if (nrow(maf_df[stringr::str_detect(maf_df$CONSEQUENCE, "^(mature_miRNA_variant|exon_variant|non_coding_exon_variant|non_coding_transcript_exon_variant|non_coding_transcript_variant|nc_transcript_variant)"), ]) > 0){
    maf_df[stringr::str_detect(maf_df$CONSEQUENCE, "^(mature_miRNA_variant|exon_variant|non_coding_exon_variant|non_coding_transcript_exon_variant|non_coding_transcript_variant|nc_transcript_variant)"), ]$Variant_Classification <- "RNA"
  }
  if (nrow(maf_df[stringr::str_detect(maf_df$CONSEQUENCE, "^(splice_region_variant, intron_variant|transcript_amplification|intron_variant|INTRAGENIC|intragenic_variant)"), ]) > 0){
    maf_df[stringr::str_detect(maf_df$CONSEQUENCE, "^(splice_region_variant, intron_variant|transcript_amplification|intron_variant|INTRAGENIC|intragenic_variant)"), ]$Variant_Classification <- "Intron"
  }
  if (nrow(maf_df[stringr::str_detect(maf_df$CONSEQUENCE, "downstream_gene_variant"), ]) > 0){
    maf_df[stringr::str_detect(maf_df$CONSEQUENCE, "downstream_gene_variant"), ]$Variant_Classification <- "3'Flank"
  }
  if (nrow(maf_df[stringr::str_detect(maf_df$CONSEQUENCE, "upstream_gene_variant"), ])){
    maf_df[stringr::str_detect(maf_df$CONSEQUENCE, "upstream_gene_variant"), ]$Variant_Classification <- "5'Flank"
  }

  maf_df <- dplyr::select(maf_df, Hugo_Symbol, Chromosome, NCBI_Build, Start_Position, End_Position, Reference_Allele, Tumor_Seq_Allele2, Tumor_Sample_Barcode, Variant_Classification, Amino_Acid_Change)
  chrom_order <- c("1", "2", "3", "4", "5", "6", "7", "8", "9", "10", "11", "12", "13", "14", "15", "16", "17", "18", "19", "20", "21", "22", "X", "Y")
  maf_df$Chromosome <- factor(maf_df$Chromosome, levels = chrom_order)
  maf_df <- dplyr::arrange(maf_df, Chromosome, Start_Position, End_Position)
  return(maf_df)
}


#' Function that generates cancer genome report - Tier model pcgr.0
#'
#' @param project_directory name of project directory
#' @param query_vcf2tsv name of gzipped TSV file (vcf2tsv) with annotated query SNVs/InDels
#' @param pcgr_data List of data frames with PCGR data annotations
#' @param pcg_config Object with PCGR configuration parameters
#' @param cna_segments_tsv name of CNA segments file (tab-separated values)
#' @param sample_name sample identifier
#' @param pcgr_version PCGR software version
#' @param genome_assembly human genome assembly version
#'

generate_report <- function(project_directory, query_vcf2tsv, pcgr_data, pcgr_config = NULL, sample_name = "SampleX",
                            cna_segments_tsv = NULL, pcgr_version = "0.6.0", genome_assembly = "grch37"){

  pcg_report <- pcgrr::init_pcg_report(pcgr_config, sample_name, pcgr_version, genome_assembly, class = NULL, pcgr_data = pcgr_data)

  genome_seq <- BSgenome.Hsapiens.UCSC.hg38
  assembly <- "hg38"
  if (genome_assembly == "grch37"){
    genome_seq <- BSgenome.Hsapiens.UCSC.hg19
    assembly <- "hg19"
  }

  fnames <- list()
  fnames[["tsv_unfiltered"]] <- paste0(project_directory, "/", sample_name, ".pcgr.", genome_assembly, ".snvs_indels.tiers.unfiltered.tsv")
  fnames[["tsv"]] <- paste0(project_directory, "/", sample_name, ".pcgr.", genome_assembly, ".snvs_indels.tiers.tsv")
  fnames[["cna_print"]] <- paste0(project_directory, "/", sample_name, ".pcgr.", genome_assembly, ".cna_segments.tsv")
  fnames[["maf"]] <- paste0(project_directory, "/", sample_name, ".pcgr.", genome_assembly, ".maf")
  fnames[["json"]] <- paste0(project_directory, "/", sample_name, ".pcgr.", genome_assembly, ".json")

  if (query_vcf2tsv != "None.gz"){
    if (!file.exists(query_vcf2tsv) | file.size(query_vcf2tsv) == 0){
      rlogging::warning(paste0("File ", query_vcf2tsv, " does not exist or has zero size"))
    }
    else{
      if (!is.null(pcgr_config) & query_vcf2tsv != "None.gz"){
        sample_calls <- pcgrr::get_calls(query_vcf2tsv, pcgr_data, pcgr_version, sample_name, pcgr_config, genome_seq, assembly)
        pcg_report_seqmode <- pcgrr::init_pcg_report(pcgr_config, sample_name, pcgr_version, genome_assembly, class = "sequencing_mode")
        pcg_report_seqmode[["eval"]] <- TRUE
        if (nrow(sample_calls) > 0){
          if (pcgr_config[["tumor_only"]][["vcf_tumor_only"]] == TRUE){
            pcg_report_seqmode[["mode"]] <- "Tumor-only (no matching control)"
            pcg_report_seqmode[["tumor_only"]] <- TRUE
            pcg_report_tumor_only <- pcgrr::generate_report_data_tumor_only(sample_calls, pcgr_data, pcgr_version,
                                                                            sample_name, pcgr_config, genome_seq, genome_assembly = assembly)
            # pcg_report_snv_indel_unfiltered <- pcgrr::generate_report_data_snv_indel_pcgr(pcg_report_tumor_only[['variant_set']][['unfiltered']],
            #                                                                               pcgr_data, pcgr_version, sample_name, pcgr_config, genome_seq,
            #                                                                               assembly, callset = "unfiltered callset")
            pcg_report_snv_indel_filtered <- pcgrr::generate_report_data_snv_indel_pcgr(pcg_report_tumor_only[["variant_set"]][["filtered"]],
                                                                                        pcgr_data, pcgr_version, sample_name, pcgr_config, genome_seq,
                                                                                        assembly, callset = "germline-filtered callset")

            pcg_report_tumor_only[["variant_set"]] <- NULL
            pcg_report <- pcgrr::update_pcg_report(pcg_report, pcg_report_snv_indel_filtered, analysis_element = "snv_indel")
            pcg_report <- pcgrr::update_pcg_report(pcg_report, pcg_report_tumor_only, analysis_element = "tumor_only")
          }else{
            pcg_report_snv_indel <- pcgrr::generate_report_data_snv_indel_pcgr(sample_calls, pcgr_data, pcgr_version, sample_name, pcgr_config, genome_seq, assembly)
            pcg_report <- pcgrr::update_pcg_report(pcg_report, pcg_report_snv_indel, analysis_element = "snv_indel")
            if (pcgr_config[["mutational_signatures"]][["mutsignatures"]] == T){
              pcg_report_signatures <- pcgrr::generate_report_data_signatures(sample_calls, pcgr_data, pcgr_version, sample_name, pcgr_config, genome_seq, assembly)
              pcg_report <- pcgrr::update_pcg_report(pcg_report, pcg_report_signatures, analysis_element = "m_signature")
            }
            if (pcgr_config[["msi"]][["msi"]] == T){
              pcg_report_msi <- pcgrr::generate_report_data_msi(sample_calls, pcgr_data, pcgr_version, sample_name, pcgr_config, genome_seq, assembly)
              pcg_report <- pcgrr::update_pcg_report(pcg_report, pcg_report_msi, analysis_element = "msi")
            }
            if (pcgr_config[["mutational_burden"]][["mutational_burden"]] == T){
              pcg_report_tmb <- pcgrr::generate_report_data_tmb(sample_calls, pcgr_data, pcgr_version, sample_name, pcgr_config, genome_seq, assembly)
              pcg_report <- pcgrr::update_pcg_report(pcg_report, pcg_report_tmb, analysis_element = "tmb")
            }
          }
          pcg_report <- pcgrr::update_pcg_report(pcg_report, pcg_report_seqmode, analysis_element = "sequencing_mode")
        }else{
          pcg_report$snv_indel$zero <- TRUE
          pcg_report[["pcgr_config"]][["other"]][["list_noncoding"]] <- FALSE
          pcg_report <- pcgrr::update_pcg_report(pcg_report, pcg_report_seqmode, analysis_element = "sequencing_mode")
        }
      }
    }
  }
  else{
    pcg_report[["pcgr_config"]][["other"]][["list_noncoding"]] <- F
  }

  if (!is.null(cna_segments_tsv)){
    if (file.exists(cna_segments_tsv)){
      pcg_report_cna <- pcgrr::generate_report_data_cna(cna_segments_tsv, pcgr_data, pcgr_version, sample_name, pcgr_config, genome_seq, genome_assembly, transcript_overlap_pct = pcgr_config[["cna"]][["cna_overlap_pct"]])
      pcg_report <- pcgrr::update_pcg_report(pcg_report, pcg_report_cna, analysis_element = "cna")
    }
  }
  for (fname_key in c("maf", "tsv", "tsv_unfiltered", "cna_print")){
    if (fname_key == "cna_print"){
      if (!is.null(pcg_report[["cna"]][["variant_set"]][[fname_key]])){
        if (nrow(pcg_report[["cna"]][["variant_set"]][[fname_key]]) > 0){
          write.table(pcg_report[["cna"]][["variant_set"]][[fname_key]], file = fnames[[fname_key]], sep = "\t", col.names = T, row.names = F, quote = F)
          gzip_command <- paste0("gzip -f ", fnames[[fname_key]])
          system(gzip_command, intern = F)
        }
      }
    }
    else{
      if (!is.null(pcg_report[["snv_indel"]][["variant_set"]][[fname_key]])){
        if (nrow(pcg_report[["snv_indel"]][["variant_set"]][[fname_key]]) > 0){
          write.table(pcg_report[["snv_indel"]][["variant_set"]][[fname_key]], file = fnames[[fname_key]], sep = "\t", col.names = T, row.names = F, quote = F)
        }
      }
    }
  }

  pcg_report_value_box <- pcgrr::generate_report_data_value_box(pcg_report, pcgr_data, pcgr_version, sample_name, pcgr_config, genome_seq, genome_assembly)
  pcg_report <- pcgrr::update_pcg_report(pcg_report, pcg_report_value_box, analysis_element = "value_box")

  for (elem in c("tier1", "tier2", "tier3", "tier4")){
    stat <- paste0("n_", elem)
    pcg_report[["snv_indel"]][["variant_statistic"]][[stat]] <- nrow(pcg_report[["snv_indel"]][["variant_set"]][[elem]])
    pcg_report[["snv_indel"]][["variant_set"]][[elem]] <- NULL
  }
  pcg_report[["snv_indel"]][["variant_set"]][["noncoding"]] <- NULL
  pcg_report[["snv_indel"]][["variant_set"]][["coding"]] <- NULL
  pcg_report[["snv_indel"]][["variant_set"]][["all"]] <- NULL
  pcg_report[['cna']][['variant_set']][['cna_print']] <- NULL

  pcgr_json <- jsonlite::toJSON(pcg_report, pretty = T, na = "string", null = "null")
  write(pcgr_json, fnames[["json"]])
  gzip_command <- paste0("gzip -f ", fnames[['json']])
  system(gzip_command, intern = F)

  rmarkdown::render(system.file("templates", "report.Rmd", package = "pcgrr"), output_format = rmarkdown::html_document(theme = pcg_report[["pcgr_config"]][["visual"]][["report_theme"]], toc = T, toc_depth = 3, toc_float = T, number_sections = F, includes = rmarkdown::includes(after_body = "disclaimer.md")), output_file = paste0(sample_name, ".pcgr.", genome_assembly, ".html"), output_dir = project_directory, clean = T, intermediates_dir = project_directory, quiet = T)


}

#' Function that generates a data frame with basic biomarker annotations from tier1 variants
#'
#' @param tier1_variants df with tier 1 variants
#' @param sample_id Sample identifier
#'
#' @return tsv_variants data frame with all tier 1 biomarkers for tab-separated output
#'
generate_biomarker_tsv <- function(clinical_evidence_items, sample_name = "test", tier = "tier 1"){

  bm_tags <- c("CLINICAL_SIGNIFICANCE", "EVIDENCE_LEVEL", "EVIDENCE_TYPE", "EVIDENCE_DIRECTION", "CANCER_TYPE",
               "DESCRIPTION", "BIOMARKER_MAPPING", "THERAPEUTIC_CONTEXT", "RATING", "CITATION")
  all_biomarker_tags <- c(c("GENOMIC_CHANGE", "GENOME_VERSION", "VCF_SAMPLE_ID", "SYMBOL", "CONSEQUENCE"),bm_tags)
  tsv_biomarkers <- NULL
  if(nrow(clinical_evidence_items) > 0){
    biomarker_tsv <- clinical_evidence_items
    biomarker_tsv$VCF_SAMPLE_ID <- sample_name
    biomarker_tsv <- biomarker_tsv %>% dplyr::select(dplyr::one_of(all_biomarker_tags))
    biomarker_tsv$TIER <- "TIER 1"
    biomarker_tsv$TIER_DESCRIPTION <- "Clinical biomarker - Predictive/prognostic/diagnostic"

    tmp2 <- as.data.frame(biomarker_tsv %>%
                            dplyr::rowwise() %>%
                            dplyr::mutate(CITATION2 = paste(unlist(stringr::str_replace_all(stringr::str_match_all(CITATION,">.+<"),"^>|<$", "")),collapse =";")))

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
generate_tier_tsv <- function(variant_set, pcgr_data, pcgr_config, sample_name = "test"){

  tags <- NULL
  if(!is.null(pcgr_config[['custom_tags']])){
    if(pcgr_config[['custom_tags']][['custom_tags']] != ""){
      tags <- stringr::str_split(pcgr_config[['custom_tags']][['custom_tags']],pattern = ",")[[1]]
    }
  }
  rlogging::message("Generating tiered set of result variants for output in tab-separated values (TSV) file")
  tsv_variants <- NULL
  for(tier in c("tier1", "tier2", "tier3", "tier4", "noncoding")){
    if(nrow(variant_set[[tier]]) > 0){
      tierset <- variant_set[[tier]]
      tierset$VCF_SAMPLE_ID <- sample_name
      tsv_columns <- pcgr_data$pcgr_tsv_tiered_columns
      if (!is.null(tags)){
        for(t in tags){
          t <- stringr::str_trim(t)
          if(t %in% colnames(tierset)){
            tsv_columns <- c(tsv_columns,t)
          }
        }
      }

      if(tier == "tier1"){
        tierset$TIER_DESCRIPTION <- "Clinical biomarker - Predictive/prognostic/diagnostic"
        if(pcgr_config[["tier_model"]][["tier_model"]] == "pcgr_acmg"){
          tierset$TIER_DESCRIPTION <- "Variants of strong clinical significance"
        }
        tierset$TIER <- "TIER 1"
      }
      if(tier == "tier2"){
        tierset$TIER_DESCRIPTION <- "Other cancer mutation hotspot/predicted driver mutation/curated cancer-associated mutation"
        if(pcgr_config[["tier_model"]][["tier_model"]] == "pcgr_acmg"){
          tierset$TIER_DESCRIPTION <- "Variants of potential clinical significance"
        }
        tierset$TIER <- "TIER 2"
      }
      if(tier == "tier3"){
        tierset$TIER_DESCRIPTION <- "Other proto-oncogene/tumor suppressor mutation"
        if(pcgr_config[["tier_model"]][["tier_model"]] == "pcgr_acmg"){
          tierset$TIER_DESCRIPTION <- "Variants of uncertain significance"
        }
        tierset$TIER <- "TIER 3"
      }
      if(tier == "tier4"){
        tierset$TIER_DESCRIPTION <- "Other coding mutation"
        tierset$TIER <- "TIER 4"
      }
      if(tier == "noncoding"){
        tierset$TIER_DESCRIPTION <- "Noncoding mutation"
        tierset$TIER <- "NONCODING"
      }
      tierset <- tierset %>% dplyr::select(dplyr::one_of(tsv_columns)) %>% dplyr::distinct()

      tsv_variants <- dplyr::bind_rows(tsv_variants, tierset)
    }
  }
  tsv_variants$GENE_NAME <- unlist(lapply(stringr::str_match_all(tsv_variants$GENE_NAME,">.+<"),paste,collapse=","))
  tsv_variants$GENE_NAME <- stringr::str_replace_all(tsv_variants$GENE_NAME,">|<", "")
  tsv_variants$CLINVAR <- unlist(lapply(stringr::str_match_all(tsv_variants$CLINVAR,">.+<"),paste,collapse=","))
  tsv_variants$CLINVAR <- stringr::str_replace_all(tsv_variants$CLINVAR,">|<", "")
  tsv_variants$PROTEIN_DOMAIN <- unlist(lapply(stringr::str_match_all(tsv_variants$PROTEIN_DOMAIN,">.+<"),paste,collapse=","))
  tsv_variants$PROTEIN_DOMAIN <- stringr::str_replace_all(tsv_variants$PROTEIN_DOMAIN,">|<", "")
  tsv_variants$TCGA_FREQUENCY <- stringr::str_replace_all(tsv_variants$TCGA_FREQUENCY,"<a href='https://portal.gdc.cancer.gov/projects/TCGA-[A-Z]{1,}' target=\"_blank\">|</a>","")
  tsv_variants <- tsv_variants %>% dplyr::distinct()

  return(tsv_variants)
}


#' Function that tiered variant sets for SNVs/InDels
#'
#' @param sample_calls variant calls subject to mutational signature analysis
#' @param pcgr_data object with PCGR annotation data
#' @param pcgr_version PCGR software version
#' @param sample_name sample identifier
#' @param pcgr_config Object with PCGR configuration parameters
#' @param genome_seq BSgenome object
#' @param genome_assembly human genome assembly version
#' @param biomarker_mapping_stringency quality level for biomarkers
#' @param callset type of calls
#'
#' @return pcg_report_data data frame with all report elements
#'
generate_report_data_snv_indel_pcgr <- function(sample_calls, pcgr_data, pcgr_version, sample_name, pcgr_config,
                                           genome_seq, genome_assembly ,callset = "somatic calls", biomarker_mapping_stringency = 1){

  rlogging::message("------")
  rlogging::message(paste0("Generating data for tiered cancer genome report - ", callset, " tier model '", pcgr_config[["tier_model"]][["tier_model"]],"'"))

  pcg_report_snv_indel <- pcgrr::init_pcg_report(pcgr_config, sample_name, pcgr_version , genome_assembly, class = "snv_indel")
  pcg_report_snv_indel[["eval"]] <- TRUE
  pcg_report_snv_indel[["variant_set"]][["all"]] <- sample_calls
  pcg_report_snv_indel[["variant_statistic"]][["n"]] <- sample_calls %>% nrow()
  pcg_report_snv_indel[["variant_statistic"]][["n_snv"]] <- sample_calls %>% dplyr::filter(VARIANT_CLASS == "SNV") %>% nrow()
  pcg_report_snv_indel[["variant_statistic"]][["n_indel"]] <- sample_calls %>% dplyr::filter(VARIANT_CLASS != "SNV") %>% nrow()
  pcg_report_snv_indel[["variant_statistic"]][["n_coding"]] <- sample_calls %>% dplyr::filter(CODING_STATUS == "coding") %>% nrow()
  #pcg_report_snv_indel[["variant_statistic"]][["n_noncoding"]] <- sample_calls %>% dplyr::filter(CODING_STATUS == "noncoding") %>% nrow()
  rlogging::message(paste0("Number of protein-coding variants: ", pcg_report_snv_indel[["variant_statistic"]][["n_coding"]]))
  #rlogging::message(paste0("Number of noncoding/silent variants: ", pcg_report_snv_indel[["variant_statistic"]][["n_noncoding"]]))

  if(!is.null(pcgr_config[['custom_tags']])){
    if(pcgr_config[['custom_tags']][['custom_tags']] != ""){
      tags <- stringr::str_split(pcgr_config[['custom_tags']][['custom_tags']],pattern = ",")[[1]]
      for(t in tags){
        t <- stringr::str_trim(t)
        if(t %in% colnames(sample_calls)){
          pcgr_data$pcgr_all_annotation_columns <- c(pcgr_data$pcgr_all_annotation_columns,t)
        }
      }
    }
  }

  if(pcg_report_snv_indel[["variant_statistic"]][["n"]] > 0){

    ## Analyze Tier1: actionable mutations and variants of clinical significance (diagnosis/prognosis etc)
    biomarker_hits_snv_indels_any <- pcgrr::get_clinical_associations_snv_indel(pcg_report_snv_indel[["variant_set"]][["all"]],pcgr_data, pcgr_config,tumor_type_specificity = "any_tumortype", biomarker_mapping_stringency = biomarker_mapping_stringency)

    if(nrow(biomarker_hits_snv_indels_any$variant_set) > 0){
      pcg_report_snv_indel[["clinical_evidence_item"]][["any_tumortype"]] <- biomarker_hits_snv_indels_any$clinical_evidence_item
      pcg_report_snv_indel[["variant_display"]][["tier1"]] <- biomarker_hits_snv_indels_any$variant_set %>% dplyr::select(GENOMIC_CHANGE) %>% dplyr::distinct()
      pcg_report_snv_indel[["variant_set"]][["tier1"]] <- biomarker_hits_snv_indels_any$variant_set
    }

    ## Analyze Tier 2: curated mutations, cancer mutation hotspots and predicted driver mutations
    pcg_report_snv_indel[["variant_set"]][["tier2"]] <- dplyr::select(pcg_report_snv_indel[["variant_set"]][["all"]], dplyr::one_of(pcgr_data$pcgr_all_annotation_columns)) %>% dplyr::filter(CODING_STATUS == "coding") %>% dplyr::filter(INTOGEN_DRIVER_MUT == TRUE | !is.na(MUTATION_HOTSPOT) | !is.na(DOCM_DISEASE))
    if(nrow(pcg_report_snv_indel[["variant_set"]][["tier1"]]) > 0){
      pcg_report_snv_indel[["variant_set"]][["tier2"]] <- dplyr::anti_join(pcg_report_snv_indel[["variant_set"]][["tier2"]], pcg_report_snv_indel[["variant_display"]][["tier1"]], by=c("GENOMIC_CHANGE"))
    }
    tier12 <- pcg_report_snv_indel[["variant_display"]][["tier1"]]
    if(nrow(pcg_report_snv_indel[["variant_set"]][["tier2"]]) > 0){
      pcg_report_snv_indel[["variant_set"]][["tier2"]] <- pcg_report_snv_indel[["variant_set"]][["tier2"]] %>% dplyr::arrange(desc(ONCOSCORE))
      tier12 <- rbind(pcg_report_snv_indel[["variant_display"]][["tier1"]],dplyr::select(pcg_report_snv_indel[["variant_set"]][["tier2"]],GENOMIC_CHANGE)) %>% dplyr::distinct()
      all_tier2_variants_display <- dplyr::select(pcg_report_snv_indel[["variant_set"]][["tier2"]], dplyr::one_of(pcgr_data$tier2_tags_display))
      pcg_report_snv_indel[["variant_display"]][["tier2"]][["hotspot"]] <- all_tier2_variants_display %>% dplyr::filter(!is.na(MUTATION_HOTSPOT))
      pcg_report_snv_indel[["variant_display"]][["tier2"]][["curated_mutation"]] <- all_tier2_variants_display %>% dplyr::filter(is.na(MUTATION_HOTSPOT) & !is.na(DOCM_DISEASE))
      pcg_report_snv_indel[["variant_display"]][["tier2"]][["predicted_driver"]] <- all_tier2_variants_display %>% dplyr::filter(is.na(MUTATION_HOTSPOT) & is.na(DOCM_DISEASE) & INTOGEN_DRIVER_MUT == TRUE)
    }

    ## Analyze Tier 3: coding mutations in oncogenes/tumor suppressors/cancer census genes
    pcg_report_snv_indel[["variant_set"]][["tier3"]] <- dplyr::select(pcg_report_snv_indel[["variant_set"]][["all"]], dplyr::one_of(pcgr_data$pcgr_all_annotation_columns)) %>% dplyr::filter(CODING_STATUS == "coding") %>% dplyr::filter(ONCOGENE == TRUE | TUMOR_SUPPRESSOR == TRUE)
    if(nrow(tier12) > 0 & nrow(pcg_report_snv_indel[["variant_set"]][["tier3"]]) > 0){
      pcg_report_snv_indel[["variant_set"]][["tier3"]] <- dplyr::anti_join(pcg_report_snv_indel[["variant_set"]][["tier3"]],tier12, by=c("GENOMIC_CHANGE"))
    }
    tier123 <- tier12
    if(nrow(pcg_report_snv_indel[["variant_set"]][["tier3"]]) > 0){
      pcg_report_snv_indel[["variant_set"]][["tier3"]] <- pcg_report_snv_indel[["variant_set"]][["tier3"]] %>% dplyr::arrange(desc(ONCOSCORE))
      tier123 <- rbind(tier12,dplyr::select(pcg_report_snv_indel[["variant_set"]][["tier3"]],GENOMIC_CHANGE)) %>% dplyr::distinct()
      pcg_report_snv_indel[["variant_display"]][["tier3"]] <- dplyr::select(pcg_report_snv_indel[["variant_set"]][["tier3"]], dplyr::one_of(pcgr_data$tier3_tags_display))
    }

    ## Analyze Tier 4: Other coding mutations
    pcg_report_snv_indel[["variant_set"]][["tier4"]] <- dplyr::select(pcg_report_snv_indel[["variant_set"]][["all"]], dplyr::one_of(pcgr_data$pcgr_all_annotation_columns)) %>% dplyr::filter(CODING_STATUS == "coding")
    if(nrow(tier123) > 0 & nrow(pcg_report_snv_indel[["variant_set"]][["tier4"]]) > 0){
      pcg_report_snv_indel[["variant_set"]][["tier4"]] <- dplyr::anti_join(pcg_report_snv_indel[["variant_set"]][["tier4"]],tier123, by=c("GENOMIC_CHANGE"))
    }
    if(nrow(pcg_report_snv_indel[["variant_set"]][["tier4"]]) > 0){
      pcg_report_snv_indel[["variant_set"]][["tier4"]] <- pcg_report_snv_indel[["variant_set"]][["tier4"]] %>% dplyr::arrange(desc(ONCOSCORE))
      pcg_report_snv_indel[["variant_display"]][["tier4"]] <- dplyr::select(pcg_report_snv_indel[["variant_set"]][["tier4"]], dplyr::one_of(pcgr_data$tier4_tags_display))
    }

    ## Analyze noncoding mutations
    pcg_report_snv_indel[["variant_set"]][["noncoding"]] <- dplyr::select(pcg_report_snv_indel[["variant_set"]][["all"]], dplyr::one_of(pcgr_data$pcgr_all_annotation_columns)) %>% dplyr::filter(CODING_STATUS == "noncoding")
    if(nrow(pcg_report_snv_indel[["variant_set"]][["noncoding"]]) > 0){
      if(nrow(tier123) > 0){
        pcg_report_snv_indel[["variant_set"]][["noncoding"]] <- dplyr::anti_join(pcg_report_snv_indel[["variant_set"]][["noncoding"]],tier123, by=c("GENOMIC_CHANGE"))
      }
      pcg_report_snv_indel[["variant_set"]][["noncoding"]] <- pcg_report_snv_indel[["variant_set"]][["noncoding"]] %>% dplyr::arrange(desc(ONCOSCORE))
      pcg_report_snv_indel[["variant_display"]][["noncoding"]] <- dplyr::select(pcg_report_snv_indel[["variant_set"]][["noncoding"]], dplyr::one_of(pcgr_data$tier5_tags_display))
    }

    pcg_report_snv_indel[["variant_statistic"]][["n_noncoding"]] <- pcg_report_snv_indel[["variant_set"]][["noncoding"]] %>% nrow()
    pcg_report_snv_indel[["variant_set"]][["tsv"]] <- pcgrr::generate_tier_tsv(pcg_report_snv_indel[["variant_set"]], pcgr_data = pcgr_data, pcgr_config, sample_name = sample_name)
    pcg_report_snv_indel[["variant_set"]][["maf"]] <- pcgrr::tier_to_maf(pcg_report_snv_indel[["variant_set"]][["tsv"]])
    #pcg_report_snv_indel[["variant_set"]][["biomarker"]] <- pcgrr::generate_biomarker_tsv(evidence_items_any, sample_name = sample_name)
  }

  rlogging::message("------")
  return(pcg_report_snv_indel)

}

#' Function that adds HTML links to different genetic variant identifiers
#'
#' @param var_df data frame with variants
#' @param vardb type of variant database
#' @param linktype type of link
#'
#' @return var_df
#'
annotate_variant_link <- function(var_df, vardb = "DBSNP", linktype = "dbsource", pcgr_data = NULL){

  if(vardb == "DBNSFP"){
    if(any(grepl(paste0("EFFECT_PREDICTIONS"),names(var_df)))){
      var_df$PREDICTED_EFFECT <- var_df$EFFECT_PREDICTIONS
      var_df$PREDICTED_EFFECT <- stringr::str_replace(var_df$PREDICTED_EFFECT,"metalr:", "<a href='https://www.ncbi.nlm.nih.gov/pubmed/25552646' target=\"_blank\">Ensembl-LogisticRegression</a>:")
      var_df$PREDICTED_EFFECT <- stringr::str_replace(var_df$PREDICTED_EFFECT,"metasvm:", "<a href='https://www.ncbi.nlm.nih.gov/pubmed/25552646' target=\"_blank\">Ensembl-SVM</a>:")
      var_df$PREDICTED_EFFECT <- stringr::str_replace(var_df$PREDICTED_EFFECT,"mutationassessor:", "<a href='http://mutationassessor.org' target=\"_blank\">MutationAssessor</a>:")
      var_df$PREDICTED_EFFECT <- stringr::str_replace(var_df$PREDICTED_EFFECT,"mutationtaster:", "<a href='http://www.mutationtaster.org' target=\"_blank\">MutationTaster</a>:")
      var_df$PREDICTED_EFFECT <- stringr::str_replace(var_df$PREDICTED_EFFECT,"fathmm:", "<a href='http://fathmm.biocompute.org.uk' target=\"_blank\">FATHMM</a>:")
      var_df$PREDICTED_EFFECT <- stringr::str_replace(var_df$PREDICTED_EFFECT,"fathmm_mkl_coding:", "<a href='http://fathmm.biocompute.org.uk/fathmmMKL.htm' target=\"_blank\">FATHMM-mkl</a>:")
      var_df$PREDICTED_EFFECT <- stringr::str_replace(var_df$PREDICTED_EFFECT,"sift:", "<a href='http://provean.jcvi.org/index.php' target=\"_blank\">SIFT</a>:")
      var_df$PREDICTED_EFFECT <- stringr::str_replace(var_df$PREDICTED_EFFECT,"lrt:", "<a href='http://www.genetics.wustl.edu/jflab/lrt_query.html' target=\"_blank\">LRT</a>:")
      var_df$PREDICTED_EFFECT <- stringr::str_replace(var_df$PREDICTED_EFFECT,"provean:", "<a href='http://provean.jcvi.org/index.php' target=\"_blank\">PROVEAN</a>:")
      var_df$PREDICTED_EFFECT <- stringr::str_replace(var_df$PREDICTED_EFFECT,"mutpred:", "<a href='http://mutpred.mutdb.org' target=\"_blank\">MutPred</a>:")
      var_df$PREDICTED_EFFECT <- stringr::str_replace(var_df$PREDICTED_EFFECT,"m-cap:", "<a href='http://bejerano.stanford.edu/MCAP/' target=\"_blank\">M-CAP</a>:")
      var_df$PREDICTED_EFFECT <- stringr::str_replace(var_df$PREDICTED_EFFECT,"splice_site_rf:", "<a href='http://nar.oxfordjournals.org/content/42/22/13534' target=\"_blank\">Splice site effect (Random forest)</a>:")
      var_df$PREDICTED_EFFECT <- stringr::str_replace(var_df$PREDICTED_EFFECT,"splice_site_ada:", "<a href='http://nar.oxfordjournals.org/content/42/22/13534' target=\"_blank\">Splice site effect (Adaptive boosting)</a>:")
      var_df$PREDICTED_EFFECT <- stringr::str_replace(var_df$PREDICTED_EFFECT,"gerp_rs:", "<a href='http://mendel.stanford.edu/SidowLab/downloads/gerp/' target=\"_blank\">GERP++ RS score</a>:")
    }
    else{
      var_df$PREDICTED_EFFECT <- NA
    }
  }

  if(vardb == "DGIDB"){
    if(any(grepl(paste0("^CHEMBL_COMPOUND_ID$"),names(var_df))) & any(grepl(paste0("^SYMBOL$"),names(var_df))) & any(grepl(paste0("^VAR_ID$"),names(var_df))) & !is.null(pcgr_data)){
      var_df_unique_slim <- dplyr::select(var_df, VAR_ID, SYMBOL, CHEMBL_COMPOUND_ID) %>% dplyr::filter(!is.na(CHEMBL_COMPOUND_ID)) %>% dplyr::distinct()
      if(nrow(var_df_unique_slim) > 0){
        var_df_unique_slim_melted <- var_df_unique_slim %>% tidyr::separate_rows(CHEMBL_COMPOUND_ID,sep="&")
        chembl_drugs <- dplyr::select(pcgr_data$dgidb,chembl_compound_id,genesymbol,drug_name,score) %>% dplyr::arrange(genesymbol,desc(score)) %>% dplyr::distinct()
        var_df_unique_slim_melted <- dplyr::left_join(var_df_unique_slim_melted, chembl_drugs, by=c("CHEMBL_COMPOUND_ID" = "chembl_compound_id", "SYMBOL" = "genesymbol")) %>% dplyr::distinct()
        if(linktype == "dbsource"){
          var_df_unique_slim_melted <- var_df_unique_slim_melted %>% dplyr::mutate(tmp_dgidb = paste0("<a href='https://www.ebi.ac.uk/chembl/compound/inspect/", CHEMBL_COMPOUND_ID,"' target=\"_blank\">", drug_name,"</a>"))
        }
        var_df_unique_slim_melted_terms <- dplyr::select(var_df_unique_slim_melted, VAR_ID, drug_name)
        var_df_terms <- dplyr::group_by(var_df_unique_slim_melted_terms, VAR_ID) %>% dplyr::summarise(CHEMBL_COMPOUND_TERMS = paste(drug_name,collapse = "&"))
        var_df_links <- dplyr::group_by(var_df_unique_slim_melted, VAR_ID) %>% dplyr::summarise(DGIDBLINK = unlist(paste(tmp_dgidb, collapse = ", ")))
        var_df_links <- dplyr::select(var_df_links, VAR_ID, DGIDBLINK) %>% dplyr::distinct()
        var_df <- dplyr::left_join(var_df, var_df_links,by=c("VAR_ID" = "VAR_ID"))
        var_df <- dplyr::left_join(var_df, var_df_terms,by=c("VAR_ID" = "VAR_ID"))

      }
      else{
        var_df$DGIDBLINK <- NA
        var_df$CHEMBL_COMPOUND_TERMS <- NA
      }
    }
    else{
      cat("WARNING: Could not generate DGIdb links - no DGIDB info provided in annotated VCF", sep="\n")
      var_df$DGIDBLINK <- NA
      var_df$CHEMBL_COMPOUND_TERMS <- NA
    }
  }

  if(vardb == "DISGENET"){
    if(any(grepl(paste0("^DISGENET_CUI$"),names(var_df))) & any(grepl(paste0("^VAR_ID$"),names(var_df))) & !is.null(pcgr_data)){
      var_df_unique_slim <- dplyr::select(var_df, VAR_ID, DISGENET_CUI) %>% dplyr::filter(!is.na(DISGENET_CUI)) %>% dplyr::distinct()
      if(nrow(var_df_unique_slim) > 0){
        var_df_unique_slim_melted <- var_df_unique_slim %>% tidyr::separate_rows(DISGENET_CUI,sep="&")
        var_df_unique_slim_melted <- dplyr::left_join(var_df_unique_slim_melted, pcgr_data$medgen_map, by=c("DISGENET_CUI" = "cui"))
        if(linktype == "dbsource"){
          var_df_unique_slim_melted <- var_df_unique_slim_melted %>% dplyr::mutate(tmp_assoc = paste0("<a href='https://www.ncbi.nlm.nih.gov/medgen/", DISGENET_CUI,"' target=\"_blank\">", cui_name,"</a>"))
        }

        var_df_unique_slim_melted_terms <- dplyr::select(var_df_unique_slim_melted, VAR_ID, cui_name)
        var_df_terms <- dplyr::group_by(var_df_unique_slim_melted_terms, VAR_ID) %>% dplyr::summarise(DISGENET_TERMS = paste(cui_name,collapse = "&"))
        var_df_links <- dplyr::group_by(var_df_unique_slim_melted, VAR_ID) %>% dplyr::summarise(DISGENET_LINK = unlist(paste(tmp_assoc, collapse = ", ")))
        var_df_links <- dplyr::select(var_df_links, VAR_ID, DISGENET_LINK)
        var_df <- dplyr::left_join(var_df, var_df_links,by=c("VAR_ID" = "VAR_ID"))
        var_df <- dplyr::left_join(var_df, var_df_terms,by=c("VAR_ID" = "VAR_ID"))
      }
      else{
        var_df$DISGENET_LINK <- NA
        var_df$DISGENET_TERMS <- NA
      }
    }
    else{
      cat("WARNING: Could not generate cancer gene association links - no Disgenet cancer associations provided in annotated VCF", sep="\n")
      var_df$DISGENET_LINK <- NA
      var_df$DISGENET_TERMS <- NA
    }
  }

  if(vardb == "TCGA"){
    if(any(grepl(paste0("^TCGA_FREQUENCY$"),names(var_df))) & any(grepl(paste0("^VAR_ID$"), names(var_df))) & !is.null(pcgr_data)){
      var_df_unique_slim <- dplyr::select(var_df, VAR_ID, TCGA_FREQUENCY) %>% dplyr::filter(!is.na(TCGA_FREQUENCY)) %>% dplyr::distinct()
      if(nrow(var_df_unique_slim) > 0){
        var_df_unique_slim_melted <- var_df_unique_slim %>% tidyr::separate_rows(TCGA_FREQUENCY, sep = ",")
        var_df_unique_slim_melted <- tidyr::separate(var_df_unique_slim_melted, TCGA_FREQUENCY, c("tumor", "percentage", "affected", "cohort_size"), sep = "\\|", convert = T)
        var_df_unique_slim_melted <- dplyr::left_join(var_df_unique_slim_melted, pcgr_data$tcga_projects, by = "tumor")
        var_df_unique_slim_melted <- dplyr::arrange(var_df_unique_slim_melted, VAR_ID, desc(percentage))
        if(linktype == "dbsource"){
          var_df_unique_slim_melted <- var_df_unique_slim_melted %>% dplyr::mutate(tmp_assoc = paste0("<a href='https://portal.gdc.cancer.gov/projects/TCGA-", tumor,"' target=\"_blank\">", name, "</a>: ", percentage, "% (", affected, "/", cohort_size, ")"))
        }

        var_df_links <- dplyr::group_by(var_df_unique_slim_melted, VAR_ID) %>% dplyr::summarise(TCGALINK = unlist(paste(tmp_assoc, collapse = ", ")))
        var_df_links <- dplyr::select(var_df_links, VAR_ID, TCGALINK)
        names(var_df_links) <- c("VAR_ID", "TCGA_FREQUENCY")
        var_df <- dplyr::rename(var_df, TCGA_FREQUENCY_RAW = TCGA_FREQUENCY)
        var_df <- dplyr::left_join(var_df, var_df_links, by = c("VAR_ID" = "VAR_ID"))
      }
    }
  }


  return(var_df)

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

generate_annotation_link <- function(df, vardb = "DBSNP", group_by_var = "VAR_ID", url_prefix = "http://www.ncbi.nlm.nih.gov/projects/SNP/snp_ref.cgi?rs=", link_key_var = "DBSNPRSID", link_display_var = "DBSNPRSID"){

  df_annotation_links <- data.frame()
  if (group_by_var %in% colnames(df) & link_key_var %in% colnames(df) & link_display_var %in% colnames(df)){
    selected_vars <- dplyr::quos(unique(c(group_by_var, link_key_var, link_display_var)))
    tmp_df <- df %>% dplyr::select(!!! selected_vars) %>% dplyr::filter(!is.na(!!rlang::sym(link_key_var)) & !is.na(!!rlang::sym(link_display_var))) %>% dplyr::distinct()
    if (nrow(tmp_df) > 0){
      tmp_df_melted <- tmp_df %>% tidyr::separate_rows(!!rlang::sym(link_key_var), sep = "&|,")
      if (vardb == "COSMIC"){
        tmp_df_melted <- tmp_df_melted %>%
          dplyr::mutate(tmp_link = paste0("<a href='", url_prefix, stringr::str_replace(!!rlang::sym(link_key_var), "COSM", ""), "' target=\"_blank\">", !!rlang::sym(link_display_var), "</a>"))
      }
      else{
        tmp_df_melted <- tmp_df_melted %>%
          dplyr::mutate(tmp_link = paste0("<a href='", url_prefix, !!rlang::sym(link_key_var), "' target=\"_blank\">", !!rlang::sym(link_display_var), "</a>"))
      }
      df_annotation_links <- as.data.frame(tmp_df_melted %>% dplyr::group_by(!! rlang::sym(group_by_var)) %>% dplyr::summarise(link = unlist(paste(tmp_link, collapse = ", "))))
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
list_to_df <- function(listfordf){
  if (!is.list(listfordf)) stop("it should be a list")

  df <- list(list.element = listfordf)
  class(df) <- c("tbl_df", "data.frame")
  attr(df, "row.names") <- .set_row_names(length(listfordf))

  if (!is.null(names(listfordf))) {
    df$name <- names(listfordf)
  }

  return(as.data.frame(df))
}

#' Function that adds read support (depth, allelic fraction) for tumor and normal and filters according to settings
#'
#' @param vcf_data_df data frame with variants
#' @param pcgr_config list with configuration values
#'
#' @return vcf_data_df
#'
add_filter_read_support <- function(vcf_data_df, pcgr_config){
  for (v in c("DP_TUMOR", "AF_TUMOR", "DP_NORMAL", "AF_NORMAL", "CALL_CONFIDENCE")){
    vcf_data_df[v] <- NA
  }
  found_tumor_tag <- 0
  for (tag_name in names(pcgr_config$allelic_support)){
    if (pcgr_config$allelic_support[[tag_name]] != "" & tag_name != "tumor_dp_min" & tag_name != "tumor_af_min" & tag_name != "normal_dp_min" & tag_name != "normal_af_max"){
      pcgr_config$allelic_support[[tag_name]] <- stringr::str_replace_all(pcgr_config$allelic_support[[tag_name]], "-", ".")
      if (pcgr_config$allelic_support[[tag_name]] %in% colnames(vcf_data_df)){
        if (tag_name == "normal_af_tag"){
          vcf_data_df[, "AF_NORMAL"] <- round(as.numeric(vcf_data_df[, pcgr_config$allelic_support[[tag_name]]]), digits = 3)
        }
        if (tag_name == "normal_dp_tag"){
          vcf_data_df[, "DP_NORMAL"] <- as.integer(vcf_data_df[, pcgr_config$allelic_support[[tag_name]]])
        }
        if (tag_name == "tumor_af_tag"){
          found_tumor_tag <- 1
          vcf_data_df[, "AF_TUMOR"] <- round(as.numeric(vcf_data_df[, pcgr_config$allelic_support[[tag_name]]]), digits = 3)
        }
        if (tag_name == "tumor_dp_tag"){
          found_tumor_tag <- 1
          vcf_data_df[, "DP_TUMOR"] <- as.integer(vcf_data_df[, pcgr_config$allelic_support[[tag_name]]])
        }
        if (tag_name == "call_conf_tag"){
          vcf_data_df[, "CALL_CONFIDENCE"] <- as.character(vcf_data_df[, pcgr_config$allelic_support[[tag_name]]])
        }
      }
    }
  }

  if (found_tumor_tag == 1){
    rlogging::message("Filtering tumor variants based on allelic depth/fraction (min_dp_tumor=", pcgr_config$allelic_support$tumor_dp_min, ", min_af_tumor=", pcgr_config$allelic_support$tumor_af_min, ")")
    rlogging::message("Filtering tumor variants based on allelic depth/fraction (min_dp_normal=", pcgr_config$allelic_support$normal_dp_min, ", max_af_normal=", pcgr_config$allelic_support$normal_af_max, ")")
    n_before_dp_af_filtering <- nrow(vcf_data_df)
    if (!any(is.na(vcf_data_df$DP_TUMOR))){
      vcf_data_df <- dplyr::filter(vcf_data_df, DP_TUMOR >= pcgr_config$allelic_support$tumor_dp_min)
    }
    if (!any(is.na(vcf_data_df$AF_TUMOR))){
      vcf_data_df <- dplyr::filter(vcf_data_df, AF_TUMOR >= pcgr_config$allelic_support$tumor_af_min)
    }
    if (!any(is.na(vcf_data_df$AF_NORMAL))){
      vcf_data_df <- dplyr::filter(vcf_data_df, AF_NORMAL <= pcgr_config$allelic_support$normal_af_max)
    }
    if (!any(is.na(vcf_data_df$DP_NORMAL))){
      vcf_data_df <- dplyr::filter(vcf_data_df, DP_NORMAL >= pcgr_config$allelic_support$normal_dp_min)
    }
    n_removed <- n_before_dp_af_filtering - nrow(vcf_data_df)
    percentage <- round(as.numeric((n_removed/n_before_dp_af_filtering) * 100), digits = 2)
    rlogging::message(paste0("Removed ", n_removed, " tumor variants (", percentage, "%) based on allelic depth/fraction"))
  }

  return(vcf_data_df)
}


#' Function that adds SwissProt feature descriptions based on keys coming from pcgr pipeline
#'
#' @param vcf_data_df Data frame of sample variants from VCF
#'
#' @return vcf_data_df
#'
#'
add_swissprot_feature_descriptions <- function(vcf_data_df, pcgr_data){
  rlogging::message("Extending annotation descriptions related to UniprotKB/SwissProt protein features")
  swissprot_features <- pcgr_data$swissprot_features
  swissprot_features$UNIPROT_FEATURE <- swissprot_features$key
  swissprot_features$PF <- paste(swissprot_features$feat_type, swissprot_features$description, tidyr::separate(swissprot_features, key, into = c("uniprot_id", "feat_type", "amino_acids"), sep = ":")$amino_acids, sep = ":")
  swissprot_features$UNIPROT_ID <- tidyr::separate(swissprot_features, key, into = c("UNIPROT_ID", "feat_type", "amino_acids"), sep = ":")$UNIPROT_ID
  swissprot_features <- dplyr::select(swissprot_features, UNIPROT_FEATURE, PF, UNIPROT_ID)

  if ("UNIPROT_FEATURE" %in% colnames(vcf_data_df) & "VAR_ID" %in% colnames(vcf_data_df) & "UNIPROT_ID" %in% colnames(vcf_data_df)){
    feature_df <- dplyr::select(vcf_data_df, UNIPROT_FEATURE, VAR_ID, UNIPROT_ID) %>% dplyr::filter(!is.na(UNIPROT_FEATURE) & !is.na(UNIPROT_ID)) %>% dplyr::distinct()
    if (nrow(feature_df) == 0){
      vcf_data_df$PROTEIN_FEATURE <- NA
      return(vcf_data_df)
    }
    feature_df <- feature_df %>% tidyr::separate_rows(UNIPROT_FEATURE, sep = "&") %>% tidyr::separate_rows(UNIPROT_FEATURE, sep = ",")
    feature_df <- as.data.frame(dplyr::left_join(feature_df, dplyr::select(swissprot_features, UNIPROT_FEATURE, UNIPROT_ID, PF), by = c("UNIPROT_FEATURE","UNIPROT_ID")))
    feature_df <- as.data.frame(dplyr::filter(feature_df, !is.na(PF)) %>% dplyr::group_by(VAR_ID) %>% dplyr::summarise(PROTEIN_FEATURE = paste(PF, collapse = ", ")))
    if (nrow(feature_df) > 0){
      if (nrow(feature_df[!is.na(feature_df$PROTEIN_FEATURE) & feature_df$PROTEIN_FEATURE == "NA", ]) > 0){
        feature_df[!is.na(feature_df$PROTEIN_FEATURE) & feature_df$PROTEIN_FEATURE == "NA", ]$PROTEIN_FEATURE <- NA
      }
      vcf_data_df <- dplyr::left_join(dplyr::select(vcf_data_df, -UNIPROT_FEATURE), feature_df, by = c("VAR_ID" = "VAR_ID"))
      if("CONSEQUENCE" %in% colnames(vcf_data_df)){
        if(nrow(vcf_data_df[!is.na(vcf_data_df$PROTEIN_FEATURE) & !stringr::str_detect(vcf_data_df$CONSEQUENCE,"^(synonymous|missense|stop|start)"),]) > 0){
          vcf_data_df[!is.na(vcf_data_df$PROTEIN_FEATURE) & !stringr::str_detect(vcf_data_df$CONSEQUENCE,"^(synonymous|missense|stop|start)"),]$PROTEIN_FEATURE <- NA
        }
      }
    }
    else{
      vcf_data_df$PROTEIN_FEATURE <- NA
    }
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
add_gwas_citation_phenotype <- function(vcf_data_df, pcgr_data, p_value_threshold = 1e-6){
  gwas_citations_phenotypes <- pcgr_data$gwas_citations_phenotypes

  gwas_citations_phenotypes <- dplyr::filter(gwas_citations_phenotypes, p_value <= p_value_threshold)
  gwas_citations_phenotypes$gwas_key <- paste(gwas_citations_phenotypes$rsid,gwas_citations_phenotypes$efo_id,gwas_citations_phenotypes$pmid,sep="_")
  gwas_citations_phenotypes <- gwas_citations_phenotypes %>%
    dplyr::mutate(GWAS_CIT = paste0(stringr::str_to_title(efo_name),", ",link," (association p-value = ",p_value,")")) %>%
    dplyr::mutate(GWAS_PH = stringr::str_to_title(efo_name)) %>%
    dplyr::filter(!is.na(GWAS_PH))
  rlogging::message("Adding citations/phenotypes underlying GWAS hits (NHGRI-EBI GWAS Catalog)")

  if ("GWAS_HIT" %in% colnames(vcf_data_df) & "VAR_ID" %in% colnames(vcf_data_df)){
    feature_df <- dplyr::select(vcf_data_df, GWAS_HIT, VAR_ID) %>% dplyr::filter(!is.na(GWAS_HIT)) %>% dplyr::distinct()
    if (nrow(feature_df) == 0){
      vcf_data_df$GWAS_CITATION <- NA
      vcf_data_df$GWAS_PHENOTYPE <- NA
      return(vcf_data_df)
    }
    feature_df <- feature_df %>% tidyr::separate_rows(GWAS_HIT, sep = ",")
    feature_df <- tidyr::separate(feature_df, GWAS_HIT, into = c("rsid", "pmid", "tagsnp","p_value","efo_id","do_id"), sep = "\\|", remove = F) %>%
      dplyr::mutate(gwas_key = paste(rsid,efo_id,pmid, sep="_"))
    feature_df <- as.data.frame(dplyr::left_join(feature_df, dplyr::select(gwas_citations_phenotypes, gwas_key, GWAS_PH, GWAS_CIT), by = c("gwas_key")))
    feature_df <- as.data.frame(feature_df %>% dplyr::filter(!is.na(GWAS_PH)) %>% dplyr::group_by(VAR_ID) %>% dplyr::summarise(GWAS_PHENOTYPE = paste(unique(GWAS_PH), collapse = "; "), GWAS_CITATION = paste(unique(GWAS_CIT), collapse = "; ")))
    if (nrow(feature_df) > 0){
      for(c in c('GWAS_CITATION','GWAS_PHENOTYPE')){
        if (nrow(feature_df[!is.na(feature_df[,c]) & feature_df[,c] == "NA", ]) > 0){
          feature_df[!is.na(feature_df[,c]) & feature_df[,c] == "NA", ][,c] <- NA
        }
      }
      vcf_data_df <- dplyr::left_join(vcf_data_df, feature_df, by = c("VAR_ID" = "VAR_ID"))
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

#' Function that assigns genotype (het/hom) from VCF GT tag
#'
#' @param vcf_data_df
#'
#' @return vcf_data_df
#'
determine_genotype <- function(vcf_data_df){

  heterozygous_states <- c()
  ref_allele_index <- 0
  while(ref_allele_index < 20){
    alt_allele_index <- ref_allele_index + 1
    while(alt_allele_index <= 20){
      phased_gt_1 <- paste0(ref_allele_index,'|',alt_allele_index)
      phased_gt_2 <- paste0(alt_allele_index,'|',ref_allele_index)
      unphased_gt_1 <- paste0(ref_allele_index,'/',alt_allele_index)
      unphased_gt_2 <- paste0(alt_allele_index,'/',ref_allele_index)
      heterozygous_states <- c(heterozygous_states,phased_gt_1,phased_gt_2,unphased_gt_1,unphased_gt_2)
      alt_allele_index <- alt_allele_index + 1
    }
    ref_allele_index <- ref_allele_index + 1
  }
  homozygous_states <- c()
  hom_allele_index <- 1
  while(hom_allele_index <= 10){
    phased_gt <- paste0(hom_allele_index,'|',hom_allele_index)
    unphased_gt <- paste0(hom_allele_index,'/',hom_allele_index)
    homozygous_states <- c(homozygous_states,phased_gt,unphased_gt)
    hom_allele_index <- hom_allele_index + 1
  }

  if("GT" %in% colnames(vcf_data_df)){
    vcf_data_df$GENOTYPE <- "ND"
    vcf_data_df[!is.na(vcf_data_df$GT) & vcf_data_df$GT %in% heterozygous_states,]$GENOTYPE <- 'heterozygous'
    vcf_data_df[!is.na(vcf_data_df$GT) & vcf_data_df$GT %in% homozygous_states,]$GENOTYPE <- 'homozygous'
  }
  else{
    vcf_data_df$GENOTYPE <- "NA"
  }
  return(vcf_data_df)
}
#' Function that reads a fully annotated VCF from PCGR VEP/vcfanno pipeline
#'
#' @param tsv_gz_file Bgzipped VCF file
#' @param pcgr_data object with PCGR annotation data
#' @param pcgr_version PCGR software version
#' @param sample_name sample identifier
#' @param pcgr_config Object with PCGR configuration parameters
#' @param genome_seq BSgenome object
#' @param genome_assembly human genome assembly version
#'
#' @return vcf_data_df
#'
get_calls <- function(tsv_gz_file, pcgr_data, pcgr_version, sample_name, pcgr_config, genome_seq = BSgenome.Hsapiens.UCSC.hg19, genome_assembly = "hg19"){

  vcf_data_df <- read.table(gzfile(tsv_gz_file), skip = 1, sep = "\t", header = T, stringsAsFactors = F, quote = "", comment.char = "", na.strings = c("."))
  if (nrow(vcf_data_df) == 0){
    return(vcf_data_df)
  }

  vcf_data_df$CHROM <- as.character(vcf_data_df$CHROM)
  n_before_exclusion <- nrow(vcf_data_df)
  nuclear_chromosomes <- c("1", "2", "3", "4", "5", "6", "7", "8", "9", "10", "11", "12", "13", "14", "15", "16", "17", "18", "19", "20", "21", "22", "X", "Y")
  nuclear_chromosomes_df <- data.frame("CHROM" = nuclear_chromosomes, stringsAsFactors = F)
  vcf_data_df <- dplyr::semi_join(vcf_data_df, nuclear_chromosomes_df, by = "CHROM")
  n_after_exclusion <- nrow(vcf_data_df)
  rlogging::message(paste0("Excluding ", n_before_exclusion - n_after_exclusion, " variants from non-nuclear chromosomes/scaffolds"))
  if (nrow(vcf_data_df) == 0){
    return(vcf_data_df)
  }
  vcf_data_df <- pcgrr::order_variants(vcf_data_df)

  hg_version <- "grch38"
  if (genome_assembly == "hg19"){
    hg_version <- "grch37"
  }

  pcgr_columns <- c("GENOME_VERSION", "PROTEIN_CHANGE", "CONSEQUENCE", "CODING_STATUS", "GENOMIC_CHANGE", "VAR_ID",
                    "DOCM_DISEASE", "DOCM_LITERATURE", "CLINVAR", "CLINVAR_TRAITS_ALL", "GENE_NAME", "GENENAME",
                    "TARGETED_DRUGS", "CANCER_ASSOCIATIONS", "DBSNP", "COSMIC", "PROTEIN_DOMAIN","GENOTYPE","CLINVAR_PHENOTYPE",
                    "AF_TUMOR","DP_TUMOR","AF_NORMAL","DP_NORMAL","CALL_CONFIDENCE","PFAM_DOMAIN_NAME","GWAS_CITATION","GWAS_PHENOTYPE")
  vcf_data_df <- vcf_data_df[, !(colnames(vcf_data_df) %in% pcgr_columns)]


  vcf_data_df <- dplyr::mutate(vcf_data_df, GENOME_VERSION = hg_version, PROTEIN_CHANGE = HGVSp) %>%
    dplyr::rename(CONSEQUENCE = Consequence)

  if (nrow(vcf_data_df[!is.na(vcf_data_df$PROTEIN_CHANGE) & stringr::str_detect(vcf_data_df$PROTEIN_CHANGE, ":"), ]) > 0){
    vcf_data_df[!is.na(vcf_data_df$PROTEIN_CHANGE) & stringr::str_detect(vcf_data_df$PROTEIN_CHANGE, ":"), ]$PROTEIN_CHANGE <-
      stringr::str_split_fixed(vcf_data_df[!is.na(vcf_data_df$PROTEIN_CHANGE) & stringr::str_detect(vcf_data_df$PROTEIN_CHANGE, ":"), ]$PROTEIN_CHANGE, pattern = ":", 2)[, 2]
  }
  if (nrow(vcf_data_df[!is.na(vcf_data_df$PROTEIN_CHANGE) & stringr::str_detect(vcf_data_df$PROTEIN_CHANGE, "^ENSP"), ]) > 0){
    vcf_data_df[!is.na(vcf_data_df$PROTEIN_CHANGE) & stringr::str_detect(vcf_data_df$PROTEIN_CHANGE, "^ENSP"), ]$PROTEIN_CHANGE <- NA
  }

  vcf_data_df$CODING_STATUS <- "noncoding"
  coding_consequence_pattern <- "^(stop_|start_lost|frameshift_|missense_variant|splice_donor|splice_acceptor|inframe_)"
  if (nrow(vcf_data_df[!is.na(vcf_data_df$CONSEQUENCE) & stringr::str_detect(vcf_data_df$CONSEQUENCE, coding_consequence_pattern), ]) > 0){
    vcf_data_df[!is.na(vcf_data_df$CONSEQUENCE) & stringr::str_detect(vcf_data_df$CONSEQUENCE, coding_consequence_pattern), ]$CODING_STATUS <- "coding"
  }
  for (v in c("ONCOGENE", "TUMOR_SUPPRESSOR", "NETWORK_CG")){
    vcf_data_df[, v] <- as.logical(dplyr::recode(vcf_data_df[, v], True = TRUE, False = FALSE))
  }

  if("ENTREZ_ID" %in% colnames(vcf_data_df$ENTREZ_ID)){
    vcf_data_df$ENTREZ_ID <- as.character(vcf_data_df$ENTREZ_ID)
  }

  if("INTOGEN_DRIVER_MUT" %in% colnames(vcf_data_df)){
    if (nrow(vcf_data_df[!is.na(vcf_data_df$INTOGEN_DRIVER_MUT), ]) > 0){
      vcf_data_df[!is.na(vcf_data_df$INTOGEN_DRIVER_MUT), ]$INTOGEN_DRIVER_MUT <- TRUE
    }
    if (nrow(vcf_data_df[is.na(vcf_data_df$INTOGEN_DRIVER_MUT), ]) > 0){
      vcf_data_df[is.na(vcf_data_df$INTOGEN_DRIVER_MUT), ]$INTOGEN_DRIVER_MUT <- FALSE
    }
  }

  if (!is.null(sample_name) & nrow(vcf_data_df) > 0){
    vcf_data_df$VCF_SAMPLE_ID <- sample_name
  }

  vcf_data_df <- pcgrr::add_filter_read_support(vcf_data_df, pcgr_config)
  vcf_data_df <- pcgrr::determine_genotype(vcf_data_df)

  rlogging::message(paste0("Number of PASS variants: ", nrow(vcf_data_df)))
  if (any(grepl(paste0("VARIANT_CLASS$"), names(vcf_data_df)))){
    n_snvs <- nrow(vcf_data_df[!is.na(vcf_data_df$VARIANT_CLASS) & vcf_data_df$VARIANT_CLASS == "SNV", ])
    n_deletions <- nrow(vcf_data_df[!is.na(vcf_data_df$VARIANT_CLASS) & vcf_data_df$VARIANT_CLASS == "deletion", ])
    n_insertions <- nrow(vcf_data_df[!is.na(vcf_data_df$VARIANT_CLASS) & vcf_data_df$VARIANT_CLASS == "insertion", ])
    n_substitutions <- nrow(vcf_data_df[!is.na(vcf_data_df$VARIANT_CLASS) & vcf_data_df$VARIANT_CLASS == "substitution", ])
    rlogging::message(paste0("Number of SNVs: ", n_snvs))
    rlogging::message(paste0("Number of deletions: ", n_deletions))
    rlogging::message(paste0("Number of insertions: ", n_insertions))
    rlogging::message(paste0("Number of block substitutions: ", n_substitutions))
  }

  if (nrow(vcf_data_df) == 0){
    return(vcf_data_df)
  }

  vcf_data_df$GENOMIC_CHANGE <- paste0(vcf_data_df$CHROM, ":g.", vcf_data_df$POS, vcf_data_df$REF, ">", vcf_data_df$ALT)
  vcf_data_df$VAR_ID <- paste(vcf_data_df$CHROM, vcf_data_df$POS, vcf_data_df$REF, vcf_data_df$ALT, sep = "_")
  for (col in c("PFAM_DOMAIN", "ENTREZ_ID", "UNIPROT_FEATURE", "Gene")){
    vcf_data_df[, col] <- as.character(vcf_data_df[, col])
  }

  if("CLINVAR_MSID" %in% colnames(vcf_data_df)){
    vcf_data_df$CLINVAR_MSID <- as.integer(vcf_data_df$CLINVAR_MSID)
  }

  if ("LoF" %in% colnames(vcf_data_df)){
    vcf_data_df$LOSS_OF_FUNCTION <- FALSE
    if(nrow(vcf_data_df[!is.na(vcf_data_df$LoF),]) > 0){
      vcf_data_df[!is.na(vcf_data_df$LoF) & vcf_data_df$LoF == "HC", ]$LOSS_OF_FUNCTION <- TRUE
    }
    ## Ignore LoF predictions for missense variants (bug in LofTee?)
    if(nrow(vcf_data_df[!is.na(vcf_data_df$CONSEQUENCE) & vcf_data_df$CONSEQUENCE == 'missense_variant' & vcf_data_df$LOSS_OF_FUNCTION == T,]) > 0){
      vcf_data_df[!is.na(vcf_data_df$CONSEQUENCE) & vcf_data_df$CONSEQUENCE == 'missense_variant' & vcf_data_df$LOSS_OF_FUNCTION == T,]$LOSS_OF_FUNCTION <- FALSE
    }
  }

  if ("EXON" %in% colnames(vcf_data_df)){
    vcf_data_df$EXON <- as.integer(stringr::str_split_fixed(vcf_data_df$EXON, "/", 2)[, 1])
  }
  vcf_data_df <- pcgrr::add_swissprot_feature_descriptions(vcf_data_df, pcgr_data = pcgr_data)
  if("GWAS_HIT" %in% colnames(vcf_data_df)){
    vcf_data_df <- pcgrr::add_gwas_citation_phenotype(vcf_data_df, pcgr_data = pcgr_data, p_value_threshold = pcgr_config$gwas$p_value_min)
  }
  vcf_data_df <- dplyr::left_join(vcf_data_df, dplyr::select(pcgr_data$pfam_domains, pfam_id, pfam_name), by = c("PFAM_DOMAIN" = "pfam_id")) %>%
    dplyr::rename(PFAM_DOMAIN_NAME = pfam_name)

  rlogging::message("Extending annotation descriptions related to Database of Curated Mutations (DoCM)")
  vcf_data_df <- dplyr::left_join(vcf_data_df, pcgr_data$docm_literature, by = c("VAR_ID"))

  vcf_data_df_1 <- dplyr::left_join(dplyr::filter(vcf_data_df, !is.na(ENTREZ_ID)),
                                    dplyr::filter(dplyr::select(pcgr_data$gene_xref, ENTREZ_ID, Gene, GENENAME), !is.na(ENTREZ_ID)), by = c("ENTREZ_ID", "Gene"))
  vcf_data_df_2 <- dplyr::left_join(dplyr::filter(vcf_data_df, is.na(ENTREZ_ID)), dplyr::select(pcgr_data$gene_xref, Gene, GENENAME), by = c("Gene"))
  vcf_data_df <- rbind(vcf_data_df_1, vcf_data_df_2)
  vcf_data_df <- pcgrr::order_variants(vcf_data_df)

  rlogging::message("Extending annotation descriptions related to KEGG pathways")
  vcf_data_df <- dplyr::left_join(vcf_data_df, pcgr_data$kegg_gene_pathway_links, by = c("ENTREZ_ID" = "gene_id")) %>%
    dplyr::rename(KEGG_PATHWAY = kegg_pathway_urls)

  clinvar <- dplyr::select(pcgr_data$clinvar, CLINVAR_TRAITS_ALL, CLINVAR_MSID, var_id) %>%
    dplyr::rename(VAR_ID = var_id) %>%
    dplyr::mutate(VAR_ID = stringr::str_replace(VAR_ID, "chr", ""))
  if ("CLINVAR_MSID" %in% colnames(vcf_data_df)){
    rlogging::message("Extending annotation descriptions related to ClinVar")
    vcf_data_df <- dplyr::left_join(vcf_data_df, clinvar, by = c("CLINVAR_MSID", "VAR_ID"))
  }

  vcf_data_df <- pcgrr::df_string_replace(vcf_data_df, strings = c("CONSEQUENCE"), pattern = "&", replacement = ", ", replace_all = T)
  vcf_data_df <- pcgrr::df_string_replace(vcf_data_df, strings = c("VEP_ALL_CONSEQUENCE", "DOCM_DISEASE","MUTATION_HOTSPOT_CANCERTYPE","ICGC_PCAWG_OCCURRENCE"), pattern = ",", replacement = ", ", replace_all = T)

  if ("EFFECT_PREDICTIONS" %in% colnames(vcf_data_df)){
    vcf_data_df$EFFECT_PREDICTIONS <- stringr::str_replace_all(vcf_data_df$EFFECT_PREDICTIONS, "\\.&|\\.$", "NA&")
    vcf_data_df$EFFECT_PREDICTIONS <- stringr::str_replace_all(vcf_data_df$EFFECT_PREDICTIONS, "&$", "")
    vcf_data_df$EFFECT_PREDICTIONS <- stringr::str_replace_all(vcf_data_df$EFFECT_PREDICTIONS, "&", ", ")
    vcf_data_df <- pcgrr::annotate_variant_link(vcf_data_df, vardb = "DBNSFP")
  }

  if ("TCGA_FREQUENCY" %in% colnames(vcf_data_df)){
    vcf_data_df <- pcgrr::annotate_variant_link(vcf_data_df, vardb = "TCGA", pcgr_data = pcgr_data)
  }

  if ("ONCOSCORE" %in% colnames(vcf_data_df)){
    if (nrow(vcf_data_df[is.na(vcf_data_df$ONCOSCORE), ]) > 0){
      vcf_data_df[is.na(vcf_data_df$ONCOSCORE), ]$ONCOSCORE <- 0
    }
  }

  if (!("DBSNP" %in% colnames(vcf_data_df))){
    dbsnp_annotation_links <- pcgrr::generate_annotation_link(vcf_data_df,
                                                     group_by_var = "VAR_ID",
                                                     url_prefix = "http://www.ncbi.nlm.nih.gov/projects/SNP/snp_ref.cgi?rs=",
                                                     link_key_var = "DBSNPRSID",
                                                     link_display_var = "DBSNPRSID")
    if (nrow(dbsnp_annotation_links) > 0){
      vcf_data_df <- dplyr::left_join(vcf_data_df, dplyr::rename(dbsnp_annotation_links, DBSNP = link), by = c("VAR_ID"))
    }else{
      vcf_data_df$DBSNP <- NA
    }
  }
  if (!("CLINVAR" %in% colnames(vcf_data_df))){
    clinvar_annotation_links <- pcgrr::generate_annotation_link(vcf_data_df,
                                                       group_by_var = "VAR_ID",
                                                       url_prefix = "http://www.ncbi.nlm.nih.gov/clinvar/variation/",
                                                       link_key_var = "CLINVAR_MSID",
                                                       link_display_var = "CLINVAR_TRAITS_ALL")
    if (nrow(clinvar_annotation_links) > 0){
      vcf_data_df <- dplyr::left_join(vcf_data_df, dplyr::rename(clinvar_annotation_links, CLINVAR = link), by = c("VAR_ID"))
    }else{
      vcf_data_df$CLINVAR <- NA
    }
  }
  if (!("GENE_NAME" %in% colnames(vcf_data_df))){
    entrezgene_annotation_links <- pcgrr::generate_annotation_link(vcf_data_df,
                                                          group_by_var = "VAR_ID",
                                                          url_prefix = "http://www.ncbi.nlm.nih.gov/gene/",
                                                          link_key_var = "ENTREZ_ID",
                                                          link_display_var = "GENENAME")
    if (nrow(entrezgene_annotation_links) > 0){
      vcf_data_df <- dplyr::left_join(vcf_data_df, dplyr::rename(entrezgene_annotation_links, GENE_NAME = link), by = c("VAR_ID"))
    }else{
      vcf_data_df$GENE_NAME <- NA
    }
  }
  if (!("PROTEIN_DOMAIN" %in% colnames(vcf_data_df))){
    pfam_annotation_links <- pcgrr::generate_annotation_link(vcf_data_df,
                                                                   group_by_var = "VAR_ID",
                                                                   url_prefix = "http://pfam.xfam.org/family/",
                                                                   link_key_var = "PFAM_DOMAIN",
                                                                   link_display_var = "PFAM_DOMAIN_NAME")
    if (nrow(pfam_annotation_links) > 0){
      vcf_data_df <- dplyr::left_join(vcf_data_df, dplyr::rename(pfam_annotation_links, PROTEIN_DOMAIN = link), by = c("VAR_ID"))
    }else{
      vcf_data_df$PROTEIN_DOMAIN <- NA
    }
  }
  if (!("COSMIC" %in% colnames(vcf_data_df))){
    cosmic_annotation_links <- pcgrr::generate_annotation_link(vcf_data_df,
                                                               vardb = "COSMIC",
                                                      group_by_var = "VAR_ID",
                                                      url_prefix = "http://cancer.sanger.ac.uk/cosmic/mutation/overview?id=",
                                                      link_key_var = "COSMIC_MUTATION_ID",
                                                      link_display_var = "COSMIC_MUTATION_ID")
    if (nrow(cosmic_annotation_links) > 0){
      vcf_data_df <- dplyr::left_join(vcf_data_df, dplyr::rename(cosmic_annotation_links, COSMIC = link), by = c("VAR_ID"))
    }else{
      vcf_data_df$COSMIC <- NA
    }
  }

  if (!("TARGETED_DRUGS" %in% colnames(vcf_data_df))){
    vcf_data_df <- pcgrr::annotate_variant_link(vcf_data_df, vardb = "DGIDB", pcgr_data = pcgr_data)
    vcf_data_df <- dplyr::rename(vcf_data_df, TARGETED_DRUGS = DGIDBLINK)
  }
  if (!("CANCER_ASSOCIATIONS" %in% colnames(vcf_data_df))){
    vcf_data_df <- pcgrr::annotate_variant_link(vcf_data_df, vardb = "DISGENET", pcgr_data = pcgr_data)
    vcf_data_df <- dplyr::rename(vcf_data_df, CANCER_ASSOCIATIONS = DISGENET_LINK)
  }

  return(vcf_data_df)

}

order_variants <- function(df){
  if ("CHROM" %in% colnames(df) & "POS" %in% colnames(df)){
    df$CHROM <- ordered(df$CHROM, levels = c("1", "2", "3", "4", "5", "6", "7", "8", "9", "10", "11", "12", "13", "14", "15", "16", "17", "18", "19", "20", "21", "22", "X", "Y"))
    df <- dplyr::arrange(df, CHROM, POS)
    df$CHROM <- as.character(df$CHROM)
  }
  return(df)
}

parse_transvar_file <- function(transvar_output_fname){

  transvar_output_raw <- as.data.frame(
    readr::read_tsv(transvar_output_fname,col_names = T) %>%
    janitor::clean_names() %>%
    tidyr::separate(coordinates_g_dna_c_dna_protein,c('gdna','cdna','hgvsp'),sep="/") %>%
    dplyr::rename(transcript_id = transcript, transvar_id = input) %>%
    dplyr::mutate(transcript_id = stringr::str_replace(transcript_id," \\(protein_coding\\)","")) %>%
    dplyr::filter(stringr::str_detect(gdna,">|del[A-Z]{1,}ins[A-Z]{1,}")) %>%
    dplyr::distinct()
  )

  transvar_gdna_snv <- transvar_output_raw %>%
    dplyr::select(gdna, transcript_id, transvar_id) %>%
    dplyr::distinct()

  transvar_gdna_mnv <- dplyr::select(transvar_output_raw, transvar_id,transcript_id,info) %>%
    dplyr::mutate(info = stringr::str_replace_all(info,"^CSQN=\\S+;reference_codon=\\S+;candidate_codons=((A|C|G|T){3},?){1,};","")) %>%
    dplyr::mutate(info = stringr::str_replace_all(info,"^source=GENCODE$","NA")) %>%
    dplyr::mutate(info = stringr::str_replace_all(info,"candidate_mnv_variants=|candidate_snv_variants=","")) %>%
    dplyr::mutate(info = stringr::str_replace_all(info,";?aliases=ENSP[0-9]{1,}\\.[0-9]{1,}","")) %>%
    dplyr::mutate(info = stringr::str_replace_all(info,";source=GENCODE$","")) %>%
    dplyr::mutate(info = stringr::str_replace_all(info,";",",")) %>%
    dplyr::filter(info != 'NA' & nchar(info) > 0) %>%
    dplyr::rename(gdna = info) %>%
    tidyr::separate_rows(gdna,sep=",")

  transvar_output <- as.data.frame(
      dplyr::bind_rows(transvar_gdna_snv, transvar_gdna_mnv) %>%
      dplyr::group_by(gdna,transvar_id) %>%
      dplyr::summarise(transvar_transcript_id = paste(transcript_id, collapse=",")))

  ##insertion/deletions
  transvar_output_indels <- transvar_output %>% dplyr::filter(stringr::str_detect(gdna,"del|ins"))
  transvar_output_indels$chr_start <- stringr::str_split_fixed(stringr::str_replace(transvar_output_indels$gdna,"chr.{1,}:g\\.",''),"_",n=2)[,1]
  alleles <- stringr::str_split_fixed(stringr::str_replace(transvar_output_indels$gdna,'chr([0-9]{1,}|X|Y):g\\.[0-9]{1,}_[0-9]{1,}del',''),'ins',2)
  transvar_output_indels$refbase <- alleles[,1]
  transvar_output_indels$altbase <- alleles[,2]
  transvar_output_indels$chr_stop <- as.integer(transvar_output_indels$chr_start) + nchar(transvar_output_indels$altbase) - 1
  transvar_output_indels$chromosome <- stringr::str_replace(stringr::str_split_fixed(transvar_output_indels$gdna,':',2)[,1],'chr','')

  ##snvs
  transvar_output_snvs <- transvar_output %>% dplyr::filter(!stringr::str_detect(gdna,"del|ins"))
  transvar_output_snvs$chr_start <- stringr::str_replace(stringr::str_replace(transvar_output_snvs$gdna,'chr([0-9]{1,}|X|Y):g\\.',''),'(A|G|C|T)>(A|G|C|T)$','')
  transvar_output_snvs$chr_stop <- transvar_output_snvs$chr_start
  alleles <- stringr::str_split_fixed(stringr::str_replace(transvar_output_snvs$gdna,'chr([0-9]{1,}|X|Y):g\\.[0-9]{1,}',''),'>',2)
  transvar_output_snvs$refbase <- alleles[,1]
  transvar_output_snvs$altbase <- alleles[,2]
  transvar_output_snvs$chromosome <- stringr::str_replace(stringr::str_split_fixed(transvar_output_snvs$gdna,':',2)[,1],'chr','')


  all_transvar <- rbind(transvar_output_indels,transvar_output_snvs)
  all_transvar <- dplyr::filter(all_transvar, !stringr::str_detect(gdna,"Error|no_valid_transcript_found"))

  return(all_transvar)

}
