#' Function that generates predisposition_report - CPSR
#'
#' @param project_directory name of project/output directory
#' @param pcgr_data Object with underlying knowledge base (PCGR/CPSR)
#' @param cpsr_config Object with CPSR configuration parameters
#' (data directory, virtual panels, sample names, genome assembly etc)
#' @export

generate_cpsr_report <- function(project_directory = NULL,
                                 pcgr_data = NULL,
                                 cpsr_config = NULL){


  invisible(assertthat::assert_that(
    !is.null(cpsr_config),
    msg = "Object 'cpsr_config' cannot be NULL"))

  query_vcf2tsv <- cpsr_config[['required_args']][['query_vcf2tsv']]
  virtual_panel_id <- cpsr_config[['required_args']][['virtual_panel_id']]
  sample_name <- cpsr_config[['required_args']][['sample_name']]
  custom_bed <- cpsr_config[['custom_panel']][['bed_fname']]

  cps_report <-
    pcgrr::init_report(
      cpsr_config,
      class = NULL,
      pcgr_data = pcgr_data,
      type = "germline",
      virtual_panel_id = virtual_panel_id,
      custom_bed = custom_bed)
  if (is.null(cps_report$metadata$gene_panel)) {
    return(NULL)
  }

  cpsr_tier_display_cols <-
    c("SYMBOL",
      "CLINVAR_PHENOTYPE",
      "CONSEQUENCE",
      "PROTEIN_CHANGE",
      "GENOTYPE",
      "GENE_NAME",
      "PROTEIN_DOMAIN",
      "HGVSp",
      "HGVSc",
      "NCBI_REFSEQ",
      "CDS_CHANGE",
      "REFSEQ_MRNA",
      "MUTATION_HOTSPOT",
      "RMSK_HIT",
      "PROTEIN_FEATURE",
      "PREDICTED_EFFECT",
      "miRNA_TARGET_HIT",
      "miRNA_TARGET_HIT_PREDICTION",
      "TF_BINDING_SITE_VARIANT",
      "TF_BINDING_SITE_VARIANT_INFO",
      "GERP_SCORE",
      "LOSS_OF_FUNCTION",
      "DBSNP",
      "CLINVAR",
      "CLINVAR_CLASSIFICATION",
      "CLINVAR_REVIEW_STATUS_STARS",
      "CLINVAR_CONFLICTED",
      "CLINVAR_VARIANT_ORIGIN",
      "CPSR_CLASSIFICATION_SOURCE",
      "CPSR_CLASSIFICATION",
      "CPSR_PATHOGENICITY_SCORE",
      "CPSR_CLASSIFICATION_DOC",
      "CPSR_CLASSIFICATION_CODE",
      "ONCOGENE",
      "TUMOR_SUPPRESSOR",
      "GLOBAL_AF_GNOMAD",
      cps_report[["metadata"]][["config"]][["popgen"]][["vcftag_gnomad"]],
      "GENOMIC_CHANGE", "GENOME_VERSION")

  ## define tags/variables to display in data tables (secondary findings)
  cpsr_sf_display_cols <-
    c("SYMBOL",
      "CONSEQUENCE",
      "CLINVAR_CLASSIFICATION",
      "CLINVAR_PHENOTYPE",
      "PROTEIN_CHANGE",
      "GENOTYPE",
      "GENE_NAME",
      "PROTEIN_DOMAIN",
      "HGVSp",
      "HGVSc",
      "NCBI_REFSEQ",
      "CDS_CHANGE",
      "REFSEQ_MRNA",
      "MUTATION_HOTSPOT",
      "RMSK_HIT",
      "PROTEIN_FEATURE",
      "PREDICTED_EFFECT",
      "LOSS_OF_FUNCTION",
      "DBSNP",
      "CLINVAR",
      "CLINVAR_REVIEW_STATUS_STARS",
      "CLINVAR_CONFLICTED",
      "CLINVAR_VARIANT_ORIGIN",
      "GWAS_CITATION",
      "ONCOGENE",
      "TUMOR_SUPPRESSOR",
      "GLOBAL_AF_GNOMAD",
      cps_report[["metadata"]][["config"]][["popgen"]][["vcftag_gnomad"]],
      "GENOMIC_CHANGE", "GENOME_VERSION")

  ## define tags/variables to display in data tables (GWAS findings)
  cpsr_gwas_display_cols <-
    c("SYMBOL",
      "CONSEQUENCE",
      "GWAS_CITATION",
      "PROTEIN_CHANGE",
      "GENOTYPE",
      "LOSS_OF_FUNCTION",
      "PROTEIN_CHANGE",
      "GENE_NAME",
      "GWAS_PHENOTYPE",
      "PROTEIN_DOMAIN",
      "HGVSp",
      "HGVSc",
      "NCBI_REFSEQ",
      "CDS_CHANGE",
      "CODING_STATUS",
      "miRNA_TARGET_HIT",
      "miRNA_TARGET_HIT_PREDICTION",
      "TF_BINDING_SITE_VARIANT",
      "TF_BINDING_SITE_VARIANT_INFO",
      "GERP_SCORE",
      "REFSEQ_MRNA",
      "PROTEIN_FEATURE",
      "PREDICTED_EFFECT",
      "DBSNP",
      "GLOBAL_AF_GNOMAD",
      "GENOMIC_CHANGE",
      "GENOME_VERSION")

  ## define tags/variables to display in output TSV
  cpsr_tsv_cols <-
    c("GENOMIC_CHANGE",
      "VAR_ID",
      "GENOTYPE",
      "GENOME_VERSION",
      "VCF_SAMPLE_ID",
      "VARIANT_CLASS",
      "CODING_STATUS",
      "SYMBOL",
      "GENE_NAME",
      "CCDS",
      "ENTREZ_ID",
      "UNIPROT_ID",
      "ENSEMBL_GENE_ID",
      "ENSEMBL_TRANSCRIPT_ID",
      "REFSEQ_MRNA",
      "ONCOGENE",
      "TUMOR_SUPPRESSOR",
      "AMINO_ACID_START",
      "AMINO_ACID_END",
      "CONSEQUENCE",
      "VEP_ALL_CSQ",
      "CIVIC_ID",
      "CIVIC_ID_SEGMENT",
      "PROTEIN_CHANGE",
      "PROTEIN_DOMAIN",
      "HGVSp",
      "HGVSc",
      "LAST_EXON",
      "EXON",
      "CODING_STATUS",
      "EXON_POSITION",
      "INTRON_POSITION",
      "CDS_CHANGE",
      "CANCER_PHENOTYPE",
      "MUTATION_HOTSPOT",
      "RMSK_HIT",
      "PROTEIN_FEATURE",
      "EFFECT_PREDICTIONS",
      "LOSS_OF_FUNCTION",
      "DBMTS",
      "miRNA_TARGET_HIT",
      "miRNA_TARGET_HIT_PREDICTION",
      "REGULATORY_ANNOTATION",
      "TF_BINDING_SITE_VARIANT",
      "TF_BINDING_SITE_VARIANT_INFO",
      "GERP_SCORE",
      "DBSNP",
      "CLINVAR_CLASSIFICATION",
      "CLINVAR_MSID",
      "CLINVAR_VARIANT_ORIGIN",
      "CLINVAR_CONFLICTED",
      "CLINVAR_PHENOTYPE",
      "CLINVAR_REVIEW_STATUS_STARS",
      "N_INSILICO_CALLED",
      "N_INSILICO_DAMAGING",
      "N_INSILICO_TOLERATED",
      "N_INSILICO_SPLICING_NEUTRAL",
      "N_INSILICO_SPLICING_AFFECTED",
      "GLOBAL_AF_GNOMAD",
      cpsr_config[["popgen"]][["vcftag_gnomad"]])

  if (query_vcf2tsv != "None.gz") {
    if (!file.exists(query_vcf2tsv) | file.size(query_vcf2tsv) == 0) {
      log4r_warn(
        paste0("File ",
               query_vcf2tsv,
               " does not exist or has zero size"))

    }else{
      if (!is.null(cpsr_config) & query_vcf2tsv != "None.gz") {

        log4r_info("Retrieving ALL variants: GWAS tag SNPs / variants in ACMG secondary findings list / target predisposition genes")
        ## read calls
        calls <-
          pcgrr::get_calls(
            query_vcf2tsv,
            pcgr_data,
            sample_name,
            cpsr_config,
            oncotree =
              cps_report[["metadata"]][["phenotype_ontology"]][["oncotree_query"]],
            cpsr = TRUE)
        if (nrow(calls) == 0) {
          log4r_warn(paste0("There are zero calls in input file ",
                            "- no report will be produced"))
          return(NULL)
        }

        ## remove any existing VCF INFO column that coincide with variables
        ## established by CPSR
        cpsr_generated_cols <-
          c("CPSR_CLASSIFICATION_SOURCE",
            "CPSR_CLASSIFICATION",
            "CPSR_PATHOGENICITY_SCORE",
            "CPSR_CLASSIFICATION_DOC",
            "CPSR_CLASSIFICATION_CODE",
            "FINAL_CLASSIFICATION",
            "N_INSILICO_CALLED",
            "N_INSILICO_DAMAGING",
            "N_INSILICO_TOLERATED",
            "N_INSILICO_SPLICING_NEUTRAL",
            "N_INSILICO_SPLICING_AFFECTED",
            pcgrr::cpsr_acmg$evidence_codes$cpsr_evidence_code)

        calls <- calls[, !(colnames(calls) %in% cpsr_generated_cols)]

        log4r_info("------")
        log4r_info(
          paste0("Considering variants in the targeted predisposition genes: ",
                 paste(unique(sort(cps_report$metadata$gene_panel$genes$symbol)),
                       collapse=", ")))

        if (cpsr_config[["ignore_noncoding"]] == T) {
          n_noncoding_vars <- calls %>%
            dplyr::filter(.data$CODING_STATUS == "noncoding") %>%
            NROW()
          log4r_info(
            paste0(
              "Excluding n = ",
              n_noncoding_vars,
              " variants for classification (option --ignore_noncoding)"
              )
          )
          calls <- calls %>%
            dplyr::filter(.data$CODING_STATUS == "coding")
          if (NROW(calls) == 0) {
            log4r_warn(paste0(
              "There are zero remaining protein-coding ",
              "calls in input file - no report will be produced")
              )
            return(NULL)
          }

        }
        ## get overall call statistics
        call_stats <-
          pcgrr::variant_stats_report(calls, name = "v_stat")

        calls <- dplyr::mutate(calls, SYMBOL = as.character(.data$SYMBOL))
        cpg_calls <-
          dplyr::inner_join(calls,
                            cps_report[["metadata"]][["gene_panel"]][["genes"]],
                            by = c("SYMBOL" = "symbol"))
        cpg_call_stats <-
          pcgrr::variant_stats_report(
            cpg_calls, name = "v_stat_cpg")

        log4r_info(
          paste0("Total number of variants in target cancer predisposition genes (for TIER output): ",
                 cpg_call_stats[["v_stat_cpg"]][["n"]]))
        log4r_info(
          paste0("Number of coding variants in target cancer predisposition genes (for TIER output): ",
                 cpg_call_stats[["v_stat_cpg"]][["n_coding"]]))
        log4r_info(
          paste0(
            "Number of non-coding variants in cancer predisposition genes (for TIER output): ",
            cpg_call_stats[["v_stat_cpg"]][["n_noncoding"]]))


        secondary_calls <-
          pcgrr::retrieve_secondary_calls(
            calls,
            umls_map = pcgr_data$phenotype_ontology$umls)
        secondary_call_stats <-
          pcgrr::variant_stats_report(secondary_calls,
                                      name = "v_stat_secondary")

        if (nrow(cpg_calls) == 0) {
          cps_report$content$snv_indel$v_stat_cpg <-
            cpg_call_stats$v_stat_cpg
          return(cps_report)
        }

        ## Assign ACMG evidence codes (>30 criteria) to variants
        cpg_calls <-
          pcgrr::assign_pathogenicity_evidence(
            cpg_calls, cpsr_config, pcgr_data) %>%
          pcgrr::determine_pathogenicity_classification() %>%
          pcgrr::detect_cancer_traits_clinvar(
            oncotree =
              cps_report[["metadata"]][["phenotype_ontology"]][["oncotree"]],
            umls_map =
              pcgr_data[["phenotype_ontology"]][["umls"]])

        if (cpsr_config[['clinvar_ignore_noncancer']] == T &
           "VAR_ID" %in% colnames(cpg_calls) &
           "CLINVAR_MSID" %in% colnames(cpg_calls) &
           "CANCER_PHENOTYPE" %in% colnames(cpg_calls)) {
          n_clinvar_noncancer <- cpg_calls %>%
            dplyr::filter(!is.na(.data$CLINVAR_MSID)) %>%
            dplyr::filter(.data$CANCER_PHENOTYPE == 0) %>%
            NROW()

          log4r_info(
            paste0("ClinVar variants related to non-cancer conditions excluded",
            " from report: ", n_clinvar_noncancer)
          )

          if (n_clinvar_noncancer > 0) {
            cpg_calls_exclude <- cpg_calls %>%
              dplyr::filter(!is.na(.data$CLINVAR_MSID)) %>%
              dplyr::filter(.data$CANCER_PHENOTYPE == 0) %>%
              dplyr::select(.data$VAR_ID)

            cpg_calls <- cpg_calls %>%
              dplyr::anti_join(cpg_calls_exclude, by = "VAR_ID")

            log4r_info(
              paste0("Variants remaining after exclusion of non-cancer related",
              " ClinVar variants: ", NROW(cpg_calls)))

          }
        }

        if (nrow(cpg_calls) == 0) {
          cps_report$content$snv_indel$v_stat_cpg <-
            cpg_call_stats$v_stat_cpg
          return(cps_report)
        }

        ## Assign calls to tiers (ClinVar calls + CPSR classification
        ## for novel, non-ClinVar variants)
        snv_indel_report <-
          pcgrr::assign_variant_tiers(
            cpg_calls,
            config = cps_report[["metadata"]][["config"]],
            cpsr_display_cols = cpsr_tier_display_cols,
            cpsr_tsv_cols = cpsr_tsv_cols)

        snv_indel_report$v_stat <-
          call_stats$v_stat
        snv_indel_report$v_stat_cpg <-
          cpg_call_stats$v_stat_cpg
        snv_indel_report$v_stat_secondary <-
          secondary_call_stats$v_stat_secondary

        ## Retrieve and assign potential biomarkers to report

        ## load all clinical evidence items related to germline variants in
        ## cancer (mutations and loss-of-function mutations)
        eitems_mut_germline <- pcgrr::load_all_eitems(
          eitems_raw = pcgr_data$biomarkers,
          alteration_type = "MUT",
          origin = "Germline")

        eitems_mut_lof_germline <- pcgrr::load_all_eitems(
          eitems_raw = pcgr_data$biomarkers,
          alteration_type = "MUT_LOF",
          origin = "Germline")

        eitems_all <- dplyr::bind_rows(eitems_mut_germline,
                                       eitems_mut_lof_germline) %>%
          pcgrr::remove_cols_from_df(cnames = c("DISEASE_ONTOLOGY_ID",
                                                "VARIANT_ORIGIN",
                                                "SOURCE_DB"))

        snv_indel_report$clin_eitem <-
          pcgrr::get_germline_biomarkers(
            sample_calls = snv_indel_report[["variant_set"]][["tsv"]],
            colset = colnames(snv_indel_report[["variant_set"]][["tsv"]]),
            eitems = eitems_all)

        cps_report <-
          pcgrr::update_report(cps_report, report_data = snv_indel_report)

        cps_report[["content"]][["snv_indel"]][["variant_set"]][["tsv"]] <-
          cps_report[["content"]][["snv_indel"]][["variant_set"]][["tsv"]] %>%
          pcgrr::remove_cols_from_df(
            cnames = c("CIVIC_ID",
                       "CIVIC_ID_SEGMENT",
                       "AMINO_ACID_END",
                       "AMINO_ACID_START",
                       "EXON"))

        gene_hits <- paste(
          unique(
            cps_report[["content"]][["snv_indel"]][["variant_set"]][["tsv"]]$SYMBOL),
          collapse = ", ")
        log4r_info(paste0("Variants were found in the following cancer ",
                          "predisposition genes: ", gene_hits))

        ## secondary findings
        if (cpsr_config[["secondary_findings"]] == TRUE) {
          if (nrow(secondary_calls) > 0) {
            log4r_info(paste0("Assignment of other variants in genes ",
                              "recommended for reporting as secondary ",
                              "findings (ACMG SF v3.0)"))
            cps_report[["content"]][["snv_indel"]][["disp"]][["secondary"]] <-
              secondary_calls %>%
              dplyr::arrange(.data$LOSS_OF_FUNCTION, .data$CODING_STATUS) %>%
              dplyr::select(dplyr::one_of(cpsr_sf_display_cols))
            log4r_info(paste0(
              "Number of pathogenic variants in the ACMG secondary findings list - other ",
              "genes of clinical significance: ",
              cps_report[["content"]][["snv_indel"]][["v_stat_secondary"]][["n_coding"]]))
          }
        }

        cps_report[["content"]][["snv_indel"]][["eval"]] <- TRUE

        if (cpsr_config[["gwas"]][["run"]] == TRUE) {
          log4r_info(paste0("Assignment of other variants to hits ",
                            "from genome-wide association studies"))

          ## Assign GWAS hits to cps_report object
          cps_report[["content"]][["snv_indel"]][["disp"]][["gwas"]] <-
            dplyr::filter(calls, !is.na(.data$GWAS_HIT) & !is.na(.data$GWAS_CITATION))
          if (nrow(cps_report[["content"]][["snv_indel"]][["disp"]][["gwas"]]) > 0) {
            if (nrow(cps_report[["content"]][["snv_indel"]][["variant_set"]][["tsv"]]) > 0) {

              ## Omit variants that is already present in TIER 1-5
              cps_report[["content"]][["snv_indel"]][["disp"]][["gwas"]] <-
                dplyr::anti_join(
                  cps_report[["content"]][["snv_indel"]][["disp"]][["gwas"]],
                  cps_report[["content"]][["snv_indel"]][["variant_set"]][["tsv"]],
                  by = c("GENOMIC_CHANGE"))
            }
            ## Select variables to include for GWAS hits and arrange results
            if (nrow(cps_report[["content"]][["snv_indel"]][["disp"]][["gwas"]]) > 0) {
              cps_report[["content"]][["snv_indel"]][["disp"]][["gwas"]] <-
                cps_report[["content"]][["snv_indel"]][["disp"]][["gwas"]] %>%
                dplyr::select(dplyr::one_of(cpsr_gwas_display_cols)) %>%
                dplyr::arrange(.data$LOSS_OF_FUNCTION, .data$CODING_STATUS)
            }
          }
        }
      }
    }
  }

  cps_report[["content"]][["snv_indel"]][['max_dt_rows']] <-
    get_max_rows_pr_datatable(cps_report)

  return(cps_report)
}


#' Function that gets the maximum number of rows across different
#' tier data frames in CPSR report
#'
#' @param cps_report CPSR report structure with tier data frames
#'
#' @return max_row_nr maximum number of rows
#'
#' @export
get_max_rows_pr_datatable <- function(cps_report) {

  max_row_nr <- 0
  if (!is.null(cps_report[["content"]][["snv_indel"]][["disp"]])) {
    for(c in c("class1", "class2", "class3",
               "class4", "class5", "sf", "gwas")) {
      if (NROW(cps_report[["content"]][["snv_indel"]][["disp"]][[c]]) == 0) {
        next
      }
      t1 <- cps_report[["content"]][["snv_indel"]][["disp"]][[c]]
      if (c == "sf" | c == "gwas") {
        if (nrow(t1) > max_row_nr) {
          max_row_nr <- nrow(t1)
        }
      }else{
        t1 <- cps_report[["content"]][["snv_indel"]][["disp"]][[c]]
        num_rows_clinvar <- t1 %>%
          dplyr::filter(.data$CPSR_CLASSIFICATION_SOURCE == "ClinVar") %>%
          nrow()
        num_rows_other <- t1 %>%
          dplyr::filter(.data$CPSR_CLASSIFICATION_SOURCE == "Other") %>%
          nrow()
        if (num_rows_other > max_row_nr) {
          max_row_nr <- num_rows_other
        }
        if (num_rows_clinvar > max_row_nr) {
          max_row_nr <- num_rows_clinvar
        }
      }
    }
  }
  return(max_row_nr)
}

#' Function that counts insilico predictions of variant effects
#' (i.e. damaging/tolerated) from dbNSFP
#'
#' @param cpg_calls sample calls with dbNSFP annotations
#'
#' @return cpg_calls
#'
#' @export
get_insilico_prediction_statistics <- function(cpg_calls) {

  insilico_pathogenicity_pred_algos <-
    c("DBNSFP_SIFT", "DBNSFP_PROVEAN",
      "DBNSFP_META_RNN", "DBNSFP_FATHMM",
      "DBNSFP_MUTATIONTASTER", "DBNSFP_DEOGEN2",
      "DBNSFP_PRIMATEAI", "DBNSFP_MUTATIONASSESSOR",
      "DBNSFP_FATHMM_MKL", "DBNSFP_M_CAP",
      "DBNSFP_LIST_S2", "DBNSFP_BAYESDEL_ADDAF",
      "DBNSFP_SPLICE_SITE_ADA", "DBNSFP_SPLICE_SITE_RF")
  for (v in c("CALLED", "DAMAGING", "TOLERATED",
              "SPLICING_NEUTRAL", "SPLICING_AFFECTED")) {
    cpg_calls[, paste0("N_INSILICO_", v)] <- 0
  }

  for (algo in insilico_pathogenicity_pred_algos) {
    if (algo %in% colnames(cpg_calls)) {
      cpg_calls[!is.na(cpg_calls[, algo]) &
                  cpg_calls[, algo] != ".", "N_INSILICO_CALLED"] <-
        cpg_calls[!is.na(cpg_calls[, algo]) &
                    cpg_calls[, algo] != ".", "N_INSILICO_CALLED"] + 1
      cpg_calls[!is.na(cpg_calls[, algo]) &
                  (cpg_calls[, algo] == "D" |
                     cpg_calls[, algo] == "PD"),
                "N_INSILICO_DAMAGING"] <-
        cpg_calls[!is.na(cpg_calls[, algo]) &
                    (cpg_calls[, algo] == "D" |
                       cpg_calls[, algo] == "PD"),
                  "N_INSILICO_DAMAGING"] + 1
      cpg_calls[!is.na(cpg_calls[, algo]) &
                  cpg_calls[, algo] == "T",
                "N_INSILICO_TOLERATED"] <-
        cpg_calls[!is.na(cpg_calls[, algo]) &
                    cpg_calls[, algo] == "T",
                  "N_INSILICO_TOLERATED"] + 1
      cpg_calls[!is.na(cpg_calls[, algo]) &
                  cpg_calls[, algo] == "AS",
                "N_INSILICO_SPLICING_AFFECTED"] <-
        cpg_calls[!is.na(cpg_calls[, algo]) &
                    cpg_calls[, algo] == "AS",
                  "N_INSILICO_SPLICING_AFFECTED"] + 1
      cpg_calls[!is.na(cpg_calls[, algo]) &
                  cpg_calls[, algo] == "SN",
                "N_INSILICO_SPLICING_NEUTRAL"] <-
        cpg_calls[!is.na(cpg_calls[, algo]) &
                    cpg_calls[, algo] == "SN",
                  "N_INSILICO_SPLICING_NEUTRAL"] + 1
    }
  }
  return(cpg_calls)
}

#' Function that assigns variant pathogenicity evidence based on ACMG guidelines
#'
#' @param cpg_calls sample calls with dbnsfp annotations
#' @param cpsr_config cpsr configuration object
#' @param pcgr_data pcgr data object
#'
#' @return cpg_calls
#'
#' @export
assign_pathogenicity_evidence <- function(cpg_calls, cpsr_config, pcgr_data) {

  gad_population <- toupper(cpsr_config[["popgen"]][["pop_gnomad"]])
  gad_AN_tag <- paste0("NON_CANCER_AN_", gad_population)
  gad_AF_tag <- paste0("NON_CANCER_AF_", gad_population)
  gad_NHOMALT_tag <- paste0("NON_CANCER_NHOMALT_", gad_population)
  gad_AC_tag <- paste0("NON_CANCER_AC_", gad_population)

  #pathogenic_range_ac <- 20
  pathogenic_range_af <- pcgrr::cpsr_acmg[["pathogenic_range_gnomad"]][["af"]]
  min_an <- pcgrr::cpsr_acmg[["pathogenic_range_gnomad"]][["min_an"]]

  acmg_ev_codes <-
    c("ACMG_BA1_AD",
      ## Very high MAF (> 0.5% in gnomAD non-cancer pop subset) -
      ## min AN = 12,000, - Dominant mechanism of disease
      "ACMG_BS1_1_AD",
      ## High MAF (> 0.1% in gnomAD non-cancer pop subset) -
      ## min AN = 12,000 - Dominant mechanism of disease
      "ACMG_BS1_2_AD",
      ## Somewhat high MAF (> 0.005% in gnomAD non-cancer pop subset) -
      ## min AN = 12,000 - Dominant mechanism of disease
      "ACMG_BA1_AR",
      ## Very high MAF (> 1% in gnomAD non-cancer pop subset) -
      ## min AN = 12,000 - Recessive mechanism of disease
      "ACMG_BS1_1_AR",
      ## High MAF (> 0.3% in gnomAD non-cancer pop subset) -
      ## min AN = 12,000 - Recessive mechanism of disease
      "ACMG_BS1_2_AR",
      ## Somewhat high MAF (> 0.005% in gnomAD non-cancer pop subset) -
      ## min AN = 12,000 - Recessive mechanism of disease
      #"ACMG_BS2_1",
      ## 1 homozygote in gnomAD non-cancer pop subset -
      ## severe, early onset, highly penetrant
      #"ACMG_BS2_2",
      ## 2 homozygotes in gnomAD non-cancer pop subset -
      ## severe, early onset, highly penetrant
      #"ACMG_BS2_3",
      ## 2 homozygotes in gnomAD non-cancer pop subset -
      ## moderate, early onset, variably penetrant
      "ACMG_PM2_1",
      ## Allele count within pathogenic range (MAF < 0.005% in the
      ## population-specific non-cancer gnomAD subset, min AN = 12,000)
      "ACMG_PM2_2",
      ## Alternate allele absent in the population-specific
      ## non-cancer gnomAD subset
      "ACMG_PVS1_1",
      ## Null variant - predicted as LoF by LOFTEE - within pathogenic range
      ## - LoF established for gene
      "ACMG_PVS1_2",
      ## Null variant - not predicted as LoF by LOFTEE -
      ## within pathogenic range - LoF established for gene
      "ACMG_PVS1_3",
      ## Null variant - predicted as LoF by LOFTEE - within pathogenic range -
      ## LoF not established for gene
      "ACMG_PVS1_4",
      ## Null variant - not predicted as LoF by LOFTEE --
      ## within pathogenic range - LoF not established for gene
      "ACMG_PVS1_5",
      ## start lost - within pathogenic range - Lof established for gene
      "ACMG_PVS1_6",
      ## start lost - within pathogenic range - LoF not established for gene
      "ACMG_PVS1_7",
      ## donor/acceptor variant - predicted as LoF by LOFTEE -
      ## within pathogenic range
      ## - not last intron - LoF established for gene
      "ACMG_PVS1_8",
      ## donor/acceptor variant - last intron - within pathogenic range -
      ## LoF established for gene
      "ACMG_PVS1_9",
      ## donor/acceptor variant - not last intron - within pathogenic range
      ## - LoF not established for gene
      "ACMG_PVS1_10",
      ## donor variant at located at the +3, +4 or +5 position of the intron -
      ## within the pathogenic range (i.e. MAF < 0.005% in gnOMAD))
      "ACMG_PS1",
      ## Same amino acid change as a previously established pathogenic
      ## variant (ClinVar) regardless of nucleotide change
      "ACMG_PP2",
      ## Missense variant in a gene that has a relatively low rate of
      ## benign missense variation (<20%) and
      ## where missense variants are a common mechanism of disease
      ## (>50% of high-confidence pathogenic variants (ClinVar))
      "ACMG_PM4",
      ## Protein length changes due to inframe indels or nonstop variant
      ## in non-repetitive regions of genes
      ## that harbor variants with a dominant mode of inheritance.
      "ACMG_PPC1",
      ## Protein length changes due to inframe indels or nonstop variant
      ## in non-repetitive regions of genes
      ## that harbor variants with a recessive mode of inheritance.
      "ACMG_PM5",
      ## Novel missense change at an amino acid residue where a different
      ## missense change determined to be pathogenic
      ## has been seen before (ClinVar)
      "ACMG_PP3",
      ## Multiple lines of computational evidence support a
      ## deleterious effect on the gene or gene product
      ## (conservation, evolutionary, splicing impact, etc. - from dbNSFP
      "ACMG_BP4",
      ## Multiple lines of computational evidence support a benign
      ## effect on the gene or gene product
      ## (conservation, evolutionary, splicing impact, etc. - from dbNSFP
      "ACMG_BMC1",
      ## Peptide change is at the same location of a
      ## known benign change (ClinVar)
      "ACMG_BSC1",
      ## Peptide change is reported as benign (ClinVar),
      "ACMG_BP3",
      ## Variants in promoter or untranslated regions
      "ACMG_BP7",
      ## Silent/intronic variant outside of the splice site consensus
      "ACMG_BP1")
       ## Missense variant in a gene for which primarily truncating
       ##variants are known to cause disease (ClinVar)


    path_columns <-
      c(acmg_ev_codes, "N_INSILICO_CALLED",
        "N_INSILICO_DAMAGING", "N_INSILICO_TOLERATED",
        "N_INSILICO_SPLICING_NEUTRAL",
        "N_INSILICO_SPLICING_AFFECTED", "codon_prefix",
        "clinvar_pathogenic_codon",
        "clinvar_pathogenic", "clinvar_benign",
        "clinvar_benign_codon", "hotspot_symbol",
        "hotspot_codon", "hotspot_pvalue",
        "MOD", "PATH_TRUNCATION_RATE", "BENIGN_MISSENSE_RATE")
  cpg_calls <- cpg_calls[, !(colnames(cpg_calls) %in% path_columns)]


  cpg_calls <- cpg_calls %>%
    dplyr::mutate(
      cpsr_gene_moi =
        dplyr::if_else(stringr::str_detect(.data$CANCER_PREDISPOSITION_MOI,
                                           "AD|AD/AR"),
                       "AD", as.character(NA), as.character(NA))) %>%
    dplyr::mutate(
      cpsr_gene_moi =
        dplyr::if_else(!stringr::str_detect(.data$CANCER_PREDISPOSITION_MOI,
                                            "AD|Mosaic"),
                       "AR", .data$cpsr_gene_moi, .data$cpsr_gene_moi))

  predisposition_gene_info <-
    dplyr::select(pcgr_data[["predisposition"]][["genes"]], .data$symbol,
                  .data$mechanism_of_disease, .data$path_truncation_rate,
                  .data$benign_missense_rate) %>%
    dplyr::rename(SYMBOL = .data$symbol,
                  MOD = .data$mechanism_of_disease,
                  PATH_TRUNCATION_RATE = .data$path_truncation_rate,
                  BENIGN_MISSENSE_RATE = .data$benign_missense_rate)

  cpg_calls <-
    dplyr::left_join(cpg_calls,
                     predisposition_gene_info, by = c("SYMBOL" = "SYMBOL"))

  ## Assign logical ACMG evidence indicators
  #
  #
  # ACMG_PP3 - Multiple lines (>=5) of insilico evidence support a
  #             deleterious effect on the gene or gene product
  ##           (conservation, evolutionary, splicing impact, etc.)
  # ACMG_BP4 - Multiple lines (>=5) of insilico evidence support a benign effect.
  #
  # Computational evidence for deleterious/benign effect is taken from
  # invidual algorithm predictions in dbNSFP: SIFT,Provean,MutationTaster,
  # MutationAssessor,M_CAP,MutPred,FATHMM,FATHMM-mkl,DBNSFP_RNN,dbscSNV_RF,
  # dbscSNV_AdaBoost
  # Default scheme (from default TOML file):
  # 1) Damaging: Among all possible protein variant effect predictions, at
  #              least six algorithms must have made a call,
  #              with at least 5 predicted as damaging/D
  #              (possibly_damaging/PD), and at most two
  #              predicted as tolerated/T (PP3)
  #       - at most 1 prediction for a splicing neutral effect
  #    Exception: if both splice site predictions indicate damaging effects;
  #    ignore other criteria
  # 2) Tolerated: Among all possible protein variant effect predictions, at
  #    least six algorithms must have made a call,
  #    with at least 5 predicted as tolerated, and at most one
  #    predicted as damaging (BP4)
  #    - 0 predictions of splice site affected

  dbnsfp_min_majority <- 7
  dbnsfp_max_minority <- 3
  dbnsfp_min_called <- dbnsfp_min_majority

  cpg_calls <- pcgrr::get_insilico_prediction_statistics(cpg_calls)

  cpg_calls <- cpg_calls %>%
    dplyr::mutate(
      ACMG_PP3 =
        dplyr::if_else(.data$N_INSILICO_CALLED >= dbnsfp_min_called &
                         .data$N_INSILICO_DAMAGING >= dbnsfp_min_majority &
                         .data$N_INSILICO_TOLERATED <= dbnsfp_max_minority &
                         .data$N_INSILICO_SPLICING_NEUTRAL <= 1, TRUE,
                       FALSE, FALSE)) %>%
    dplyr::mutate(
      ACMG_BP4 = dplyr::if_else(.data$N_INSILICO_CALLED >= dbnsfp_min_called &
                                  .data$N_INSILICO_TOLERATED >= dbnsfp_min_majority &
                                  .data$N_INSILICO_DAMAGING <= dbnsfp_max_minority &
                                  .data$N_INSILICO_SPLICING_AFFECTED == 0, TRUE,
                                FALSE, FALSE)) %>%
    dplyr::mutate(ACMG_PP3 = dplyr::case_when(
      .data$N_INSILICO_SPLICING_AFFECTED == 2 ~ TRUE, TRUE ~ as.logical(.data$ACMG_PP3)))

  ## Assign logical ACMG evidence indicators based on population frequency
  ## data in non-cancer samples from gnomAD (Dominant vs. recessive
  ## modes of inheritance)
  # 'ACMG_BA1_AD'   - Very high MAF (> 0.5% in gnomAD non-cancer pop subset) -
  #                   min AN = 12,000 - Dominant mechanism of disease
  # 'ACMG_BS1_1_AD' - High MAF (> 0.1% in gnomAD non-cancer pop subset) -
  #                   min AN = 12,000 - Dominant mechanism of disease
  # 'ACMG_BS1_2_AD' - Somewhat high MAF (> 0.005% in gnomAD non-cancer pop
  #                   subset) - Dominant mechanism of disease
  # 'ACMG_BA1_AR'   - Very high MAF (> 1% in gnomAD non-cancer pop subset) -
  #                   min AN = 12,000 - Recessive mechanism of disease
  # 'ACMG_BS1_1_AR' - High MAF (> 0.3% in gnomAD non-cancer pop subset) -
  #                   min AN = 12,000 - Recessive mechanism of disease
  # 'ACMG_BS1_2_AR' - Somewhat high MAF (> 0.005% in gnomAD non-cancer pop
  #                   subset) - Recessive mechanism of disease
  # 'ACMG_PM2_1'    - Allele count within pathogenic range (MAF <= 0.005%
  #                   in the population-specific non-cancer gnomAD subset,
  #                   min AN = 12,000)
  # 'ACMG_PM2_2'    - Alternate allele absent in the population-specific
  #                   non-cancer gnomAD subset
  if (gad_AN_tag %in% colnames(cpg_calls) &
      gad_AC_tag %in% colnames(cpg_calls) &
      gad_NHOMALT_tag %in% colnames(cpg_calls)) {

    cpg_calls <- cpg_calls %>%
      dplyr::mutate(
        gad_af =
          dplyr::if_else(!!rlang::sym(gad_AN_tag) >= min_an,
                         as.numeric(!!rlang::sym(gad_AC_tag) /
                                      !!rlang::sym(gad_AN_tag)),
                         as.double(NA), as.double(NA))) %>%
      dplyr::mutate(
        ACMG_PM2_1 =
          dplyr::if_else(!!rlang::sym(gad_AN_tag) >= min_an &
                           !is.na(!!rlang::sym(gad_AC_tag)) &
                           .data$gad_af <= pathogenic_range_af,
                         TRUE, FALSE, FALSE)) %>%
      dplyr::mutate(
        ACMG_PM2_2 = dplyr::if_else(is.na(!!rlang::sym(gad_AC_tag)),
                                    TRUE, FALSE, FALSE)) %>%
      dplyr::mutate(
        ACMG_BA1_AD = dplyr::if_else(.data$ACMG_PM2_2 == FALSE &
                                       .data$gad_af >= 0.005 &
                                       .data$cpsr_gene_moi == "AD",
                                     TRUE, FALSE, FALSE)) %>%
      dplyr::mutate(
        ACMG_BS1_1_AD = dplyr::if_else(.data$ACMG_BA1_AD == FALSE &
                                         .data$ACMG_PM2_2 == FALSE &
                                         .data$gad_af >= 0.001 &
                                         .data$cpsr_gene_moi == "AD",
                                       TRUE, FALSE, FALSE)) %>%
      dplyr::mutate(
        ACMG_BS1_2_AD = dplyr::if_else(.data$ACMG_BS1_1_AD == FALSE &
                                         .data$ACMG_BA1_AD == FALSE &
                                         .data$ACMG_PM2_2 == FALSE &
                                         .data$gad_af > pathogenic_range_af &
                                         .data$cpsr_gene_moi == "AD",
                                       TRUE, FALSE, FALSE)) %>%
      dplyr::mutate(
        ACMG_BA1_AR = dplyr::if_else(.data$ACMG_PM2_2 == FALSE &
                                       .data$gad_af >= 0.01 &
                                       (.data$cpsr_gene_moi == "AR" |
                                          is.na(.data$cpsr_gene_moi)),
                                     TRUE, FALSE, FALSE)) %>%
      dplyr::mutate(
        ACMG_BS1_1_AR = dplyr::if_else(.data$ACMG_BA1_AR == FALSE &
                                         .data$ACMG_PM2_2 == FALSE &
                                         .data$gad_af >= 0.003 &
                                         (.data$cpsr_gene_moi == "AR" |
                                            is.na(.data$cpsr_gene_moi)),
                                       TRUE, FALSE, FALSE)) %>%
      dplyr::mutate(
        ACMG_BS1_2_AR = dplyr::if_else(.data$ACMG_BA1_AR == FALSE &
                                         .data$ACMG_BS1_1_AR == FALSE &
                                         .data$ACMG_PM2_2 == FALSE &
                                         .data$gad_af > pathogenic_range_af &
                                         (.data$cpsr_gene_moi == "AR" |
                                            is.na(.data$cpsr_gene_moi)),
                                       TRUE, FALSE, FALSE))
  }

  ## Assign logical ACMG evidence indicators on NULL variants in known
  # predisposition genes (LoF established as mechanism of disease or not,
  # presumed loss of mRNA/protein (LOFTEE) or not)
  #
  # 'ACMG_PVS1_1' - Null variant (frameshift, nonsense) -
  #  predicted as LoF by LOFTEE - within pathogenic range - LoF established
  # 'ACMG_PVS1_2' - Null variant (frameshift, nonsense) -
  # not predicted as LoF by LOFTEE - within pathogenic range - LoF established
  # 'ACMG_PVS1_3' - Null variant (frameshift, nonsense) -
  # predicted as LoF by LOFTEE - within pathogenic range - LoF not established
  # 'ACMG_PVS1_4' - Null variant (frameshift, nonsense) -
  # not predicted as LoF by LOFTEE -- within pathogenic range - LoF not
  # established for gene
  # 'ACMG_PVS1_5' - start lost - within pathogenic range - Lof established
  # 'ACMG_PVS1_6' - start lost - within pathogenic range - LoF not established
  # 'ACMG_PVS1_7' - splice acceptor/donor variant - predicted as LoF by LOFTEE
  # - not last intron - within pathogenic range - Lof established
  # 'ACMG_PVS1_8' - splice acceptor/donor variant - predicted as LoF by LOFTEE
  # - last intron - within pathogenic range - Lof established
  # 'ACMG_PVS1_9' - splice acceptor/donor variant - predicted as LoF by LOFTEE
  # - not last intron - within pathogenic range - Lof established
  # 'ACMG_PVS1_10' - splice variant involving a donor at +3A/G, +4A or +5G -
  # predicted as damaging by insilico predictions - within pathogenic range

  cpg_calls <- cpg_calls %>%
    dplyr::mutate(
      ACMG_PVS1_1 =
        dplyr::if_else(.data$NULL_VARIANT == T & .data$LOSS_OF_FUNCTION == T &
                         .data$MOD == "LoF" &
                         (.data$ACMG_PM2_1 == TRUE | .data$ACMG_PM2_2 == TRUE),
                       TRUE, FALSE, FALSE)) %>%
    dplyr::mutate(
      ACMG_PVS1_3 =
        dplyr::if_else(.data$NULL_VARIANT == T & .data$LOSS_OF_FUNCTION == T &
                         (is.na(.data$MOD) | .data$MOD != "LoF") &
                         (.data$ACMG_PM2_1 == TRUE | .data$ACMG_PM2_2 == TRUE),
                       TRUE, FALSE, FALSE)) %>%
    dplyr::mutate(
      ACMG_PVS1_2 =
        dplyr::if_else(.data$NULL_VARIANT == T & .data$LOSS_OF_FUNCTION == F &
                         .data$MOD == "LoF" &
                         (.data$ACMG_PM2_1 == TRUE | .data$ACMG_PM2_2 == TRUE),
                       TRUE, FALSE, FALSE)) %>%
    dplyr::mutate(
      ACMG_PVS1_4 =
        dplyr::if_else(.data$NULL_VARIANT == T & .data$LOSS_OF_FUNCTION == F &
                         (is.na(.data$MOD) | .data$MOD != "LoF") &
                         (.data$ACMG_PM2_1 == TRUE | .data$ACMG_PM2_2 == TRUE),
                       TRUE, FALSE, FALSE)) %>%
    dplyr::mutate(
      ACMG_PVS1_5 =
        dplyr::if_else(.data$CONSEQUENCE == "start_lost" & .data$MOD == "LoF" &
                         (.data$ACMG_PM2_1 == TRUE | .data$ACMG_PM2_2 == TRUE),
                       TRUE, FALSE, FALSE)) %>%
    dplyr::mutate(
      ACMG_PVS1_6 =
        dplyr::if_else(.data$CONSEQUENCE == "start_lost" & (is.na(.data$MOD) |
                                                        .data$MOD != "LoF") &
                         (.data$ACMG_PM2_1 == TRUE | .data$ACMG_PM2_2 == TRUE),
                       TRUE, FALSE, FALSE)) %>%
    dplyr::mutate(
      ACMG_PVS1_7 =
        dplyr::if_else(.data$LOSS_OF_FUNCTION == T &
                         stringr::str_detect(.data$CONSEQUENCE, "_donor|_acceptor") &
                         .data$LAST_INTRON == F & .data$MOD == "LoF" &
                         (.data$ACMG_PM2_1 == TRUE | .data$ACMG_PM2_2 == TRUE),
                       TRUE, FALSE, FALSE)) %>%
    dplyr::mutate(
      ACMG_PVS1_8 =
        dplyr::if_else(.data$LOSS_OF_FUNCTION == T &
                         stringr::str_detect(.data$CONSEQUENCE, "_donor|_acceptor") &
                         .data$LAST_INTRON == T & .data$MOD == "LoF" &
                         (.data$ACMG_PM2_1 == TRUE | .data$ACMG_PM2_2 == TRUE),
                       TRUE, FALSE, FALSE)) %>%
    dplyr::mutate(
      ACMG_PVS1_9 =
        dplyr::if_else(.data$LOSS_OF_FUNCTION == T &
                         stringr::str_detect(.data$CONSEQUENCE, "_donor|_acceptor") &
                         .data$LAST_INTRON == F & (is.na(.data$MOD) | .data$MOD != "LoF") &
                         (.data$ACMG_PM2_1 == TRUE | .data$ACMG_PM2_2 == TRUE),
                       TRUE, FALSE, FALSE)) %>%
    dplyr::mutate(
      ACMG_PVS1_10 =
        dplyr::if_else(.data$SPLICE_DONOR_RELEVANT == T & .data$ACMG_PP3 == TRUE &
                         (.data$ACMG_PM2_1 == TRUE | .data$ACMG_PM2_2 == TRUE),
                       TRUE, FALSE,  FALSE))



  # Assign logical ACMG evidence indicators
  # # TODO - BA1 -  exceptions for high population germline frequency
  #  (gnomAD) - HFE/SERPINA1

  ## Assign logical ACMG evidence indicator
  # PM4 - Protein length changes (in non-repetitive regions) due to
  # inframe indels or nonstop variant of genes that harbor variants with
  # a dominant mode of inheritance
  #
  # PPC1 - Protein length changes (in non-repetitive regions) due to
  # inframe indels or nonstop variant of genes that harbor variants with a
  # recessive mode of inheritance (and unknown MOI) - PPC1
  if ("RMSK_HIT" %in% colnames(cpg_calls)) {
    cpg_calls <- cpg_calls %>%
      dplyr::mutate(
        ACMG_PM4 =
          dplyr::if_else(
            stringr::str_detect(
              .data$CONSEQUENCE, "stop_lost|inframe_deletion|inframe_insertion") &
              is.na(.data$RMSK_HIT) & .data$cpsr_gene_moi == "AD",
            TRUE, FALSE, FALSE)) %>%
      dplyr::mutate(
        ACMG_PPC1 =
          dplyr::if_else(
            stringr::str_detect(
              .data$CONSEQUENCE, "stop_lost|inframe_deletion|inframe_insertion") &
              is.na(.data$RMSK_HIT) & (.data$cpsr_gene_moi == "AR" | is.na(.data$cpsr_gene_moi)),
            TRUE, FALSE, FALSE))
  }

  ## Assign logical ACMG evidence indicator
  # ACMG_PP2 - Missense variant in a gene that has a relatively low rate
  # of benign missense variation and where missense variants are a
  # common mechanism of disease
  cpg_calls <- cpg_calls %>%
    dplyr::mutate(
      ACMG_PP2 =
        dplyr::if_else(
          (is.na(.data$BENIGN_MISSENSE_RATE) | .data$BENIGN_MISSENSE_RATE <= 0.1) &
            (is.na(.data$PATH_TRUNCATION_RATE) | .data$PATH_TRUNCATION_RATE < 0.5) &
            stringr::str_detect(.data$CONSEQUENCE, "^missense_variant"),
          TRUE, FALSE, FALSE))

  ## Assign a logical ACMG evidence indicator
  # ACMG_BP1 - Missense variant in a gene for which primarily truncating
  # variants (> 90%, as given in Maxwell et al.) are known to cause disease
  cpg_calls <- cpg_calls %>%
    dplyr::mutate(
      ACMG_BP1 =
        dplyr::if_else(.data$PATH_TRUNCATION_RATE > 0.90 &
                         stringr::str_detect(.data$CONSEQUENCE, "^missense_variant"),
                       TRUE, FALSE, FALSE))

  ## Assign a logical ACMG evidence indicator
  # ACMG_BP7 - Silent/intronic variant outside of the splice site consensus
  cpg_calls <- cpg_calls %>%
    dplyr::mutate(
      ACMG_BP7 =
        dplyr::if_else((
          (.data$INTRON_POSITION < 0 & .data$INTRON_POSITION < -3) |
            (.data$INTRON_POSITION > 0 & .data$INTRON_POSITION > 6) |
            (.data$EXON_POSITION < 0 & .data$EXON_POSITION < -2) |
            (.data$EXON_POSITION > 0 & .data$EXON_POSITION > 1)) &
            stringr::str_detect(
              .data$CONSEQUENCE,
              "^(synonymous_variant|intron_variant|splice_region_variant)"),
          TRUE, FALSE, FALSE))

  ## Assign a logical ACMG evidence indicator
  # ACMG_BP3 - Variants in promoter or untranslated regions
  cpg_calls <- cpg_calls %>%
    dplyr::mutate(
      ACMG_BP3 =
        dplyr::if_else(
          stringr::str_detect(
            .data$CONSEQUENCE,
            "^(downstream|upstream|5_prime_UTR_variant|3_prime_UTR_variant)"),
          TRUE, FALSE, FALSE))


  ## Assign logical ACMG evidence indicators
  # ACMG_PS1 - coinciding with known pathogenic missense variants
  # (yet with different nucleotide change)
  # ACMG_PM5 - occurs at the same codon as a known pathogenic missense variant
  # ACMG_BSC1 - coinciding with known benign missense variants
  # ACMG_BMC1 - occurs at the same codon as a known benign missense variant

  cpg_calls$codon_prefix <- NA
  if (nrow(cpg_calls[!is.na(cpg_calls$CONSEQUENCE) &
                     cpg_calls$CONSEQUENCE == "missense_variant", ]) > 0) {
    cpg_calls[cpg_calls$CONSEQUENCE == "missense_variant", ]$codon_prefix <-
      stringr::str_match(cpg_calls[cpg_calls$CONSEQUENCE == "missense_variant", ]$HGVSp_short, "p\\.[A-Z]{1}[0-9]{1,}")
  }

  if (nrow(cpg_calls[!is.na(cpg_calls$codon_prefix), ]) > 0) {
    cpg_calls_pathogenic_codon <-
      dplyr::left_join(
        dplyr::filter(dplyr::select(cpg_calls, .data$VAR_ID, .data$codon_prefix, .data$SYMBOL),
                      !is.na(.data$codon_prefix)),
        pcgr_data[["clinvar"]][["cpg_loci"]][["hi"]][["P"]][["codon"]],
        by = c("codon_prefix" = "codon_prefix", "SYMBOL" = "symbol")) %>%
      dplyr::filter(.data$clinvar_pathogenic_codon == T)
    cpg_calls <- cpg_calls %>%
      dplyr::left_join(dplyr::select(cpg_calls_pathogenic_codon,
                                     .data$VAR_ID, .data$clinvar_pathogenic_codon),
                       by = c("VAR_ID"))

    cpg_calls_benign_codon <-
      dplyr::left_join(
        dplyr::filter(dplyr::select(cpg_calls, .data$VAR_ID, .data$codon_prefix, .data$SYMBOL),
                      !is.na(.data$codon_prefix)),
        pcgr_data[["clinvar"]][["cpg_loci"]][["hi"]][["B"]][["codon"]],
        by = c("codon_prefix" = "codon_prefix", "SYMBOL" = "symbol")) %>%
      dplyr::filter(.data$clinvar_benign_codon == T)

    cpg_calls <- cpg_calls %>%
      dplyr::left_join(
        dplyr::select(cpg_calls_benign_codon, .data$VAR_ID, .data$clinvar_benign_codon),
        by = c("VAR_ID")) %>%
      dplyr::mutate(ACMG_PM5 = dplyr::if_else(
        .data$clinvar_pathogenic_codon == TRUE, TRUE, FALSE, FALSE)) %>%
      dplyr::mutate(ACMG_BMC1 = dplyr::if_else(
        .data$clinvar_benign_codon == TRUE, TRUE, FALSE, FALSE))
  }else{
    cpg_calls$ACMG_PM5 <- FALSE
    cpg_calls$ACMG_BMC1 <- FALSE
    cpg_calls$clinvar_pathogenic_codon <- NA
    cpg_calls$clinvar_benign_codon <- NA
  }


  if (nrow(cpg_calls[!is.na(cpg_calls$HGVSp_short), ]) > 0) {
    cpg_calls_pathogenic_hgvsp <-
      dplyr::left_join(
        dplyr::filter(
          dplyr::select(cpg_calls, .data$VAR_ID, .data$HGVSp_short, .data$SYMBOL),
          !is.na(.data$HGVSp_short)),
        pcgr_data[["clinvar"]][["cpg_loci"]][["hi"]][["P"]][["peptide_change"]],
        by = c("HGVSp_short" = "hgvs_p", "SYMBOL" = "symbol")) %>%
      dplyr::filter(.data$clinvar_pathogenic == T)

    cpg_calls <- cpg_calls %>%
      dplyr::left_join(
        dplyr::select(cpg_calls_pathogenic_hgvsp,
                      .data$VAR_ID, .data$clinvar_pathogenic), by = c("VAR_ID"))

    cpg_calls_benign_hgvsp <-
      dplyr::left_join(
        dplyr::filter(
          dplyr::select(cpg_calls, .data$VAR_ID, .data$HGVSp_short, .data$SYMBOL),
          !is.na(.data$HGVSp_short)),
        pcgr_data[["clinvar"]][["cpg_loci"]][["hi"]][["B"]][["peptide_change"]],
        by = c("HGVSp_short" = "hgvs_p", "SYMBOL" = "symbol")) %>%
      dplyr::filter(.data$clinvar_benign == T)

    cpg_calls <- cpg_calls %>%
      dplyr::left_join(
        dplyr::select(cpg_calls_benign_hgvsp, .data$VAR_ID, .data$clinvar_benign),
                       by = c("VAR_ID")) %>%
      dplyr::mutate(ACMG_PS1 =
                      dplyr::if_else(.data$clinvar_pathogenic == TRUE,
                                     TRUE, FALSE, FALSE)) %>%
      dplyr::mutate(ACMG_BSC1 =
                      dplyr::if_else(.data$clinvar_benign == TRUE,
                                     TRUE, FALSE, FALSE))
  }else{
    cpg_calls$ACMG_PS1 <- FALSE
    cpg_calls$ACMG_BSC1 <- FALSE
    cpg_calls$clinvar_pathogenic <- NA
    cpg_calls$clinvar_benign <- NA
  }

  ## if previously found coinciding with pathogenic variant (ACMG_PS1),
  # set ACMG_PM5 to false
  cpg_calls <- cpg_calls %>%
    dplyr::mutate(
      ACMG_PM5 =
        dplyr::case_when(.data$ACMG_PM5 == T & .data$ACMG_PS1 == T ~ FALSE,
                         TRUE ~ as.logical(.data$ACMG_PM5))) %>%
    ## if previously found coinciding with benign variant (ACMG_BSC1),
    ##  set ACMG_BMC1 to false
    dplyr::mutate(
      ACMG_BMC1 =
        dplyr::case_when(.data$ACMG_BMC1 == T & .data$ACMG_BSC1 == T ~ FALSE,
                         TRUE ~ as.logical(.data$ACMG_BMC1))) %>%
    dplyr::select(
      -c(.data$clinvar_pathogenic_codon, .data$clinvar_benign_codon,
         .data$clinvar_pathogenic,
         .data$clinvar_benign, .data$cpsr_gene_moi, .data$gad_af))

  ##Assign logical ACMG level
  # PM1 - missense variant in a somatic mutation hotspot as
  # determined by cancerhotspots.org (v2)
  cpg_calls <- cpg_calls %>%
    tidyr::separate(
      .data$MUTATION_HOTSPOT,
      c("hotspot_symbol", "hotspot_codon", "hotspot_pvalue"),
      sep = "\\|", remove = F, extra = "drop") %>%
    dplyr::mutate(
      hotspot_codon =
        dplyr::if_else(!is.na(.data$hotspot_codon),
                       paste0("p.", .data$hotspot_codon),
                       as.character(NA))) %>%
    dplyr::mutate(
      ACMG_PM1 =
        dplyr::if_else(!is.na(.data$hotspot_codon) &
                         !is.na(.data$hotspot_symbol) &
                         !is.na(.data$codon_prefix) &
                         .data$SYMBOL == .data$hotspot_symbol &
                         .data$hotspot_codon == .data$codon_prefix,
                       TRUE, FALSE))

  return(cpg_calls)
}


#' Function that assigns final pathogenicity classification (B, LB, VUS, P, LP)
#' based on accumulated scores from different ACMG criteria and pre-defined
#' cutoffs (calibrated against ClinVar)
#'
#' @param cpg_calls data frame with variant calls in predisposition genes
#'
#' @return cpg_calls data frame with pathogenicity classification appended
#'
#' @export
determine_pathogenicity_classification <- function(cpg_calls) {

  evidence_codes <- pcgrr::cpsr_acmg[["evidence_codes"]]

  path_cols <- c("CPSR_CLASSIFICATION", "CPSR_CLASSIFICATION_DOC",
                 "CPSR_CLASSIFICATION_CODE",
                "cpsr_score_pathogenic", "cpsr_score_benign")
  cpg_calls <- cpg_calls[, !(colnames(cpg_calls) %in% path_cols)]

  cpg_calls$CPSR_CLASSIFICATION <- "VUS"
  cpg_calls$CPSR_CLASSIFICATION_DOC <- ""
  cpg_calls$CPSR_CLASSIFICATION_CODE <- ""
  cpg_calls$cpsr_score_pathogenic <- 0
  cpg_calls$cpsr_score_benign <- 0

  i <- 1
  while (i <= nrow(evidence_codes)) {
    category <- evidence_codes[i, ]$category
    pole <- evidence_codes[i, ]$pathogenicity_pole
    description <- evidence_codes[i, ]$description
    cpsr_evidence_code <- evidence_codes[i, ]$cpsr_evidence_code
    score <- evidence_codes[i, ]$path_score
    if (cpsr_evidence_code %in% colnames(cpg_calls)) {
      cpg_calls <- cpg_calls %>%
        dplyr::mutate(
          cpsr_score_benign = .data$cpsr_score_benign +
            dplyr::if_else(pole == "B" & !!rlang::sym(cpsr_evidence_code) == T,
                           score, 0)) %>%
        dplyr::mutate(
          cpsr_score_pathogenic = .data$cpsr_score_pathogenic +
            dplyr::if_else(pole == "P" & !!rlang::sym(cpsr_evidence_code) == T,
                           score, 0)) %>%
        dplyr::mutate(
          CPSR_CLASSIFICATION_DOC =
            paste0(.data$CPSR_CLASSIFICATION_DOC,
                   dplyr::if_else(!!rlang::sym(cpsr_evidence_code) == T,
                                  paste0("- ", description), ""),
                   sep = "<br>")) %>%
        dplyr::mutate(
          CPSR_CLASSIFICATION_CODE =
            paste0(.data$CPSR_CLASSIFICATION_CODE,
                   dplyr::if_else(!!rlang::sym(cpsr_evidence_code) == T,
                                  cpsr_evidence_code, ""), sep = "|"))
    }
    i <- i + 1
  }

  lb_upper_limit <- -1.5
  lb_lower_limit <- -4.5
  b_upper_limit <- -5.0
  vus_lower_limit <- -1.0
  vus_upper_limit <- 2.0
  lp_lower_limit <- 2.5
  lp_upper_limit <- 4.5
  p_lower_limit <- 5.0


  cpg_calls <- cpg_calls %>%
    dplyr::mutate(
      CPSR_CLASSIFICATION_CODE =
        stringr::str_replace_all(
          stringr::str_replace_all(
            .data$CPSR_CLASSIFICATION_CODE,
            "(\\|{2,})", "|"),
          "(^\\|)|(\\|$)", "")
      ) %>%
    dplyr::mutate(
      CPSR_CLASSIFICATION_DOC =
        stringr::str_replace_all(
          stringr::str_replace_all(
            .data$CPSR_CLASSIFICATION_DOC,
            "(<br>){2,}", "<br>"), "(^(<br>))|((<br>)$)", "")) %>%

    ## Adjust scores in cases where critera are acting as a
    ## prerequisite for other criteria
    dplyr::mutate(
      cpsr_score_pathogenic =
        dplyr::if_else(
          stringr::str_detect(.data$CPSR_CLASSIFICATION_CODE, "ACMG_PVS") &
            stringr::str_detect(.data$CPSR_CLASSIFICATION_CODE, "ACMG_PM2_2"),
          .data$cpsr_score_pathogenic - 1, .data$cpsr_score_pathogenic)) %>%
    dplyr::mutate(
      cpsr_score_pathogenic =
        dplyr::if_else(
          stringr::str_detect(.data$CPSR_CLASSIFICATION_CODE, "ACMG_PVS") &
            stringr::str_detect(.data$CPSR_CLASSIFICATION_CODE, "ACMG_PM2_1"),
          .data$cpsr_score_pathogenic - 0.5, .data$cpsr_score_pathogenic)) %>%
    dplyr::mutate(
      cpsr_score_pathogenic =
        dplyr::if_else(
          stringr::str_detect(.data$CPSR_CLASSIFICATION_CODE, "ACMG_PVS1_10") &
            stringr::str_detect(.data$CPSR_CLASSIFICATION_CODE, "ACMG_PP3"),
          .data$cpsr_score_pathogenic - 0.5, .data$cpsr_score_pathogenic)) %>%

    ## Add scores accumulated with benign criteria and pathogenic criteria
    dplyr::mutate(
      CPSR_PATHOGENICITY_SCORE =
        dplyr::if_else(.data$cpsr_score_benign == 0,
                       .data$cpsr_score_pathogenic,
                       .data$cpsr_score_benign)) %>%
    dplyr::mutate(
      CPSR_PATHOGENICITY_SCORE =
        dplyr::if_else(.data$cpsr_score_benign < 0 &
                         .data$cpsr_score_pathogenic > 0,
                       .data$cpsr_score_benign + .data$cpsr_score_pathogenic,
                       .data$CPSR_PATHOGENICITY_SCORE)) %>%


    dplyr::mutate(
      CPSR_CLASSIFICATION =
        dplyr::case_when(
          .data$CPSR_PATHOGENICITY_SCORE <= lb_upper_limit &
            .data$CPSR_PATHOGENICITY_SCORE >= lb_lower_limit ~ "Likely_Benign",
          .data$CPSR_PATHOGENICITY_SCORE <= b_upper_limit ~ "Benign",
          .data$CPSR_PATHOGENICITY_SCORE <= vus_upper_limit &
            .data$CPSR_PATHOGENICITY_SCORE >= vus_lower_limit ~ "VUS",
          .data$CPSR_PATHOGENICITY_SCORE >= p_lower_limit ~ "Pathogenic",
          .data$CPSR_PATHOGENICITY_SCORE >= lp_lower_limit &
            .data$CPSR_PATHOGENICITY_SCORE <= lp_upper_limit ~ "Likely_Pathogenic",
          TRUE ~ as.character("VUS"))) %>%
    dplyr::select(-c(.data$cpsr_score_benign, .data$cpsr_score_pathogenic))

  # dplyr::mutate(
  #   CPSR_CLASSIFICATION =
  #     dplyr::case_when(
  #       CPSR_PATHOGENICITY_SCORE <= -1.5 &
  #         CPSR_PATHOGENICITY_SCORE >= -4.5 ~ "Likely_Benign",
  #       CPSR_PATHOGENICITY_SCORE <= -5 ~ "Benign",
  #       CPSR_PATHOGENICITY_SCORE <= 1.5 &
  #         CPSR_PATHOGENICITY_SCORE >= -1.0 ~ "VUS",
  #       CPSR_PATHOGENICITY_SCORE >= 5 ~ "Pathogenic",
  #       CPSR_PATHOGENICITY_SCORE >= 2.0 &
  #         CPSR_PATHOGENICITY_SCORE <= 4.5 ~ "Likely_Pathogenic",
  #       TRUE ~ as.character("VUS"))) %>%
  #   dplyr::select(-c(cpsr_score_benign, cpsr_score_pathogenic))

  return(cpg_calls)

}

#' Function that assign variants to different tiers for
#' prioritization of germline variants
#'
#' @param cpg_calls data frame with variants in predisposition_genes
#' @param config CPSR configuration object with run settings
#' @param cpsr_display_cols list of variables for display in report
#' @param cpsr_tsv_cols list of variables for display in report
#' @export
assign_variant_tiers <-
  function(cpg_calls,
           config = NULL,
           cpsr_display_cols = NULL,
           cpsr_tsv_cols = NULL) {


  #dot_args <- list(...)

  evidence_codes <- pcgrr::cpsr_acmg[['evidence_codes']] %>%
    dplyr::filter(.data$cpsr_evidence_code != "ACMG_BS2_1" &
                    .data$cpsr_evidence_code != "ACMG_BS2_2" &
                    .data$cpsr_evidence_code != "ACMG_BS2_3")
  cpsr_tsv_cols <-
    c(cpsr_tsv_cols,
      evidence_codes$cpsr_evidence_code,
      c("FINAL_CLASSIFICATION",
        "CPSR_CLASSIFICATION",
        "CPSR_PATHOGENICITY_SCORE",
        "CPSR_CLASSIFICATION_CODE",
        "CPSR_CLASSIFICATION_DOC",
        "CPSR_CLASSIFICATION_SOURCE"))

  ## Add custom annotation tags to lists of tags to display
  if (!is.null(config[["preserved_info_tags"]])) {
    if (config[["preserved_info_tags"]] != "None") {
      tags <- stringr::str_split(config[["preserved_info_tags"]],
                                 pattern = ",")[[1]]
      for (t in tags) {
        if(nchar(t) == 0){
          next
        }
        t <- stringr::str_trim(t)
        if (t %in% colnames(cpg_calls)) {
          cpsr_tsv_cols <- c(cpsr_tsv_cols, t)
        }else{
          log4r_warn(
            paste0("Could NOT detect the following tag in query VCF: ", tag))
        }
      }
    }
  }


  log4r_info(paste0("Generating tiered set of result variants for ",
                    "output in tab-separated values (TSV) file"))

  snv_indel_report <- pcgrr::init_report(config = config,
                                         type = "germline",
                                         class = "snv_indel")

  #predispose_tags <- cpsr_tsv_cols

  snv_indel_report[["variant_set"]][["class5"]] <- cpg_calls %>%
    dplyr::filter(!is.na(.data$CLINVAR_CLASSIFICATION) &
                    .data$CLINVAR_CLASSIFICATION == "Pathogenic") %>%
    dplyr::mutate(CPSR_CLASSIFICATION_SOURCE = "ClinVar")

  snv_indel_report[["variant_set"]][["class4"]] <- cpg_calls %>%
    dplyr::filter(!is.na(.data$CLINVAR_CLASSIFICATION) &
                    .data$CLINVAR_CLASSIFICATION == "Likely_Pathogenic") %>%
    dplyr::mutate(CPSR_CLASSIFICATION_SOURCE = "ClinVar")

  snv_indel_report[["variant_set"]][["class3"]] <- cpg_calls %>%
    dplyr::filter(!is.na(.data$CLINVAR_CLASSIFICATION) &
                    .data$CLINVAR_CLASSIFICATION == "VUS") %>%
    dplyr::mutate(CPSR_CLASSIFICATION_SOURCE = "ClinVar")

  snv_indel_report[["variant_set"]][["class2"]] <- cpg_calls %>%
    dplyr::filter(!is.na(.data$CLINVAR_CLASSIFICATION) &
                    .data$CLINVAR_CLASSIFICATION == "Likely_Benign") %>%
    dplyr::mutate(CPSR_CLASSIFICATION_SOURCE = "ClinVar")

  snv_indel_report[["variant_set"]][["class1"]] <- cpg_calls %>%
    dplyr::filter(!is.na(.data$CLINVAR_CLASSIFICATION) &
                    .data$CLINVAR_CLASSIFICATION == "Benign") %>%
    dplyr::mutate(CPSR_CLASSIFICATION_SOURCE = "ClinVar")

  ## identify remaining calls not registered in ClinVar
  all_clinvar_calls <- data.frame()
  for (c in c("class1", "class2", "class3", "class4", "class5")) {
    all_clinvar_calls <- all_clinvar_calls %>%
      dplyr::bind_rows(
        dplyr::select(
          snv_indel_report[["variant_set"]][[c]],
          .data$VAR_ID)
        )
  }
  cpg_calls_non_clinvar <- cpg_calls %>%
    dplyr::anti_join(all_clinvar_calls, by = c("VAR_ID"))

  n_nonclinvar <- nrow(cpg_calls_non_clinvar)

  cpg_calls_non_clinvar <- cpg_calls_non_clinvar %>%
    dplyr::filter(is.na(.data$GLOBAL_AF_GNOMAD) |
                    .data$GLOBAL_AF_GNOMAD <= config[["popgen"]][["maf_upper_threshold"]])
  n_maf_filtered <- n_nonclinvar - nrow(cpg_calls_non_clinvar)
  log4r_info(
    paste0("Ignoring n = ", n_maf_filtered,
           " unclassified variants with a global MAF frequency above ",
           config[["popgen"]][["maf_upper_threshold"]]))

  n_after_maf_filtering <- nrow(cpg_calls_non_clinvar)

  non_clinvar_calls <- list()
  non_clinvar_calls[["class5"]] <- cpg_calls_non_clinvar %>%
    dplyr::filter(.data$CPSR_CLASSIFICATION == "Pathogenic") %>%
    dplyr::mutate(CPSR_CLASSIFICATION_SOURCE = "Other")

  non_clinvar_calls[["class4"]] <- cpg_calls_non_clinvar %>%
    dplyr::filter(.data$CPSR_CLASSIFICATION == "Likely_Pathogenic") %>%
    dplyr::mutate(CPSR_CLASSIFICATION_SOURCE = "Other")

  non_clinvar_calls[["class3"]] <- cpg_calls_non_clinvar %>%
    dplyr::filter(.data$CPSR_CLASSIFICATION == "VUS") %>%
    dplyr::mutate(CPSR_CLASSIFICATION_SOURCE = "Other")

  non_clinvar_calls[["class2"]] <- cpg_calls_non_clinvar %>%
    dplyr::filter(.data$CPSR_CLASSIFICATION == "Likely_Benign") %>%
    dplyr::mutate(CPSR_CLASSIFICATION_SOURCE = "Other")

  non_clinvar_calls[["class1"]] <- cpg_calls_non_clinvar %>%
    dplyr::filter(.data$CPSR_CLASSIFICATION == "Benign") %>%
    dplyr::mutate(CPSR_CLASSIFICATION_SOURCE = "Other")

  for (c in c("class1", "class2", "class3", "class4", "class5")) {

    log4r_info(paste0("Merging ClinVar-classified variants and CPSR-classified (novel) variants - ",c))
    snv_indel_report[["variant_set"]][[c]] <-
      dplyr::bind_rows(non_clinvar_calls[[c]],
                       snv_indel_report[["variant_set"]][[c]])


    if(nrow(snv_indel_report[["variant_set"]][[c]]) == 0){
      log4r_info(paste0("Zero variants found - ", c))
    }

    ## set FINAL_CLASSIFICATION col
    snv_indel_report[["variant_set"]][[c]] <-
      snv_indel_report[["variant_set"]][[c]] %>%
      dplyr::mutate(
        FINAL_CLASSIFICATION = dplyr::case_when(
          !is.na(.data$CLINVAR_CLASSIFICATION) ~ as.character(.data$CLINVAR_CLASSIFICATION),
          is.na(.data$CLINVAR_CLASSIFICATION) ~ as.character(.data$CPSR_CLASSIFICATION),
          TRUE ~ as.character(NA)
        )
      )


    ## If not 'classify_all' is turned on,
    ## remove CPSR classifications for existing ClinVar classifications
    if (config[["classify_all"]] == F) {
      snv_indel_report[["variant_set"]][[c]] <-
        snv_indel_report[["variant_set"]][[c]] %>%
        dplyr::mutate(
          CPSR_CLASSIFICATION =
            dplyr::if_else(!is.na(.data$CLINVAR_CLASSIFICATION),
                           "", as.character(.data$CPSR_CLASSIFICATION))) %>%
        dplyr::mutate(
          CPSR_CLASSIFICATION_DOC =
            dplyr::if_else(!is.na(.data$CLINVAR_CLASSIFICATION),
                           "", as.character(.data$CPSR_CLASSIFICATION_DOC))) %>%
        dplyr::mutate(
          CPSR_CLASSIFICATION_CODE =
            dplyr::if_else(!is.na(.data$CLINVAR_CLASSIFICATION),
                           "", as.character(.data$CPSR_CLASSIFICATION_CODE))) %>%
        dplyr::mutate(
          CPSR_PATHOGENICITY_SCORE =
            dplyr::if_else(!is.na(.data$CLINVAR_CLASSIFICATION),
                           as.numeric(NA),
                           as.numeric(.data$CPSR_PATHOGENICITY_SCORE)))
    }
    snv_indel_report[["variant_set"]][[c]] <-
      snv_indel_report[["variant_set"]][[c]] %>%
      dplyr::arrange(.data$CPSR_CLASSIFICATION_SOURCE,
                     dplyr::desc(.data$CANCER_PHENOTYPE),
                     dplyr::desc(.data$CPSR_PATHOGENICITY_SCORE))
      #dplyr::select(-CANCER_PHENOTYPE)

    if(config[['visual']][['table_display']] == 'full'){
      cols_in_tsv_not_display_cols <-
        setdiff(cpsr_tsv_cols, cpsr_display_cols)

      cpsr_display_cols <-
        c(cpsr_display_cols,
          cols_in_tsv_not_display_cols)

      cpsr_display_cols <-
        cpsr_display_cols[!cpsr_display_cols %in% c("AMINO_ACID_START","AMINO_ACID_END","EXON","CIVIC_ID","CIVIC_ID_SEGMENT")]
    }

    snv_indel_report[["disp"]][[c]] <-
      dplyr::select(snv_indel_report[["variant_set"]][[c]],
                    dplyr::one_of(cpsr_display_cols))
    snv_indel_report[["variant_set"]][[c]] <-
      dplyr::select(snv_indel_report[["variant_set"]][[c]],
                    dplyr::one_of(cpsr_tsv_cols))

    snv_indel_report[["variant_set"]][[c]]$DBSNP <-
      unlist(lapply(
        stringr::str_match_all(snv_indel_report[["variant_set"]][[c]]$DBSNP,
                               ">rs[0-9]{1,}<"), paste, collapse = ","))
    snv_indel_report[["variant_set"]][[c]]$DBSNP <-
      stringr::str_replace_all(
        snv_indel_report[["variant_set"]][[c]]$DBSNP, ">|<", "")
    snv_indel_report[["variant_set"]][[c]]$GENE_NAME <-
      unlist(lapply(
        stringr::str_match_all(
          snv_indel_report[["variant_set"]][[c]]$GENE_NAME, ">.+<"),
        paste, collapse = ","))
    snv_indel_report[["variant_set"]][[c]]$GENE_NAME <-
      stringr::str_replace_all(
        snv_indel_report[["variant_set"]][[c]]$GENE_NAME, ">|<", "")
    snv_indel_report[["variant_set"]][[c]]$PROTEIN_DOMAIN <-
      unlist(lapply(
        stringr::str_match_all(
          snv_indel_report[["variant_set"]][[c]]$PROTEIN_DOMAIN, ">.+<"),
        paste, collapse = ","))
    snv_indel_report[["variant_set"]][[c]]$PROTEIN_DOMAIN <-
      stringr::str_replace_all(
        snv_indel_report[["variant_set"]][[c]]$PROTEIN_DOMAIN, ">|<", "")
    snv_indel_report[["variant_set"]][[c]]$CPSR_CLASSIFICATION_DOC <-
      stringr::str_replace_all(
        snv_indel_report[["variant_set"]][[c]]$CPSR_CLASSIFICATION_DOC,
        "<br>-", ",")
    snv_indel_report[["variant_set"]][[c]]$CPSR_CLASSIFICATION_DOC <-
      stringr::str_replace_all(
        snv_indel_report[["variant_set"]][[c]]$CPSR_CLASSIFICATION_DOC,
        "^, ", "")

    snv_indel_report[["variant_set"]][[c]] <-
      snv_indel_report[["variant_set"]][[c]] %>%
      dplyr::select(c("GENOMIC_CHANGE",
                      "VAR_ID",
                      "GENOTYPE",
                      "CPSR_CLASSIFICATION_SOURCE",
                      "GENOME_VERSION",
                      "VCF_SAMPLE_ID",
                      "VARIANT_CLASS",
                      "CODING_STATUS",
                      "SYMBOL",
                      "GENE_NAME",
                      "CCDS",
                      "ENTREZ_ID",
                      "UNIPROT_ID",
                      "ENSEMBL_GENE_ID",
                      "ENSEMBL_TRANSCRIPT_ID",
                      "REFSEQ_MRNA",
                      "ONCOGENE",
                      "TUMOR_SUPPRESSOR",
                      "CONSEQUENCE",
                      "VEP_ALL_CSQ",
                      "REGULATORY_ANNOTATION",
                      "PROTEIN_CHANGE",
                      "PROTEIN_DOMAIN",
                      "DBSNP",
                      "HGVSp",
                      "HGVSc",
                      "LAST_EXON",
                      "EXON_POSITION",
                      "INTRON_POSITION",
                      "CDS_CHANGE",
                      "MUTATION_HOTSPOT",
                      "RMSK_HIT",
                      "PROTEIN_FEATURE",
                      "EFFECT_PREDICTIONS",
                      "LOSS_OF_FUNCTION",
                      "DBSNP",
                      "CANCER_PHENOTYPE",
                      "CLINVAR_CLASSIFICATION",
                      "CLINVAR_MSID",
                      "CLINVAR_VARIANT_ORIGIN",
                      "CLINVAR_CONFLICTED",
                      "CLINVAR_PHENOTYPE",
                      "CLINVAR_REVIEW_STATUS_STARS"),
                    dplyr::everything())

    for (col in colnames(snv_indel_report[["variant_set"]][[c]])) {
      if (nrow(snv_indel_report[["variant_set"]][[c]][!is.na(
        snv_indel_report[["variant_set"]][[c]][, col]) &
        snv_indel_report[["variant_set"]][[c]][, col] == "", ]) > 0) {
        snv_indel_report[["variant_set"]][[c]][!is.na(
          snv_indel_report[["variant_set"]][[c]][, col]) &
            snv_indel_report[["variant_set"]][[c]][, col] == "", col] <- NA
      }
    }

    population_tags <- unique(
      c("GLOBAL_AF_GNOMAD", config[["popgen"]][["vcftag_gnomad"]]))
    for (tag in population_tags) {
      if (tag %in% colnames(snv_indel_report[["disp"]][[c]])) {
        if (nrow(snv_indel_report[["disp"]][[c]][is.na(
          snv_indel_report[["disp"]][[c]][, tag]), ]) > 0) {
          snv_indel_report[["disp"]][[c]][is.na(
            snv_indel_report[["disp"]][[c]][, tag]), tag] <- 0.00
        }
      }
    }
  }

  snv_indel_report[["variant_set"]][["tsv"]] <-
    dplyr::bind_rows(snv_indel_report[["variant_set"]][["class5"]],
                     snv_indel_report[["variant_set"]][["class4"]],
                     snv_indel_report[["variant_set"]][["class3"]],
                     snv_indel_report[["variant_set"]][["class2"]],
                     snv_indel_report[["variant_set"]][["class1"]])


  return(snv_indel_report)
}

#' Function that retrieves variants in genes recommended for secondary
#' (secondary) findings
#'
#' @param calls data frame with variants in genes recommended for
#' secondary findings reporting
#' @param umls_map data frame with UMLS phenotype terms
#'
#' @export
retrieve_secondary_calls <- function(calls, umls_map) {

  secondary_calls <- calls %>%
    dplyr::filter(!is.na(.data$CANCER_PREDISPOSITION_SOURCE) &
                    .data$CANCER_PREDISPOSITION_SOURCE == "ACMG_SF30" &
                    !is.na(.data$CLINVAR_CLASSIFICATION)) %>%
    dplyr::filter(.data$CLINVAR_CLASSIFICATION == "Pathogenic" |
                    .data$CLINVAR_CLASSIFICATION == "Likely_Pathogenic") %>%
    ## only LOF for TTN
    dplyr::filter(.data$SYMBOL != "TTN" | (.data$SYMBOL == "TTN" & .data$LOSS_OF_FUNCTION == T))


  ## AR genes
  min_two_variants_required <-
    data.frame(SYMBOL = 'MUTYH', stringsAsFactors = F) %>%
    dplyr::bind_rows(
      data.frame(SYMBOL = 'RPE65', stringsAsFactors = F),
      data.frame(SYMBOL = 'GAA', stringsAsFactors = F),
      data.frame(SYMBOL = 'BTD', stringsAsFactors = F),
      data.frame(SYMBOL = 'ATP7B', stringsAsFactors = F)
    )

  if (nrow(secondary_calls) > 0) {

    ## Skip MUTYH/RPE65/GAA/BTD/ATP7B if they only occur with one variant
    genes_lacking_twohit_evidence <- secondary_calls %>%
      dplyr::group_by(.data$SYMBOL) %>%
      dplyr::summarise(n = dplyr::n(), .groups = "drop") %>%
      dplyr::filter(.data$n == 1)

    if(nrow(genes_lacking_twohit_evidence) > 0){
      genes_lacking_twohit_evidence <- genes_lacking_twohit_evidence %>%
        dplyr::inner_join(min_two_variants_required)

      if(nrow(genes_lacking_twohit_evidence) > 0){
        secondary_calls <- secondary_calls %>%
          dplyr::anti_join(genes_lacking_twohit_evidence, by = "SYMBOL")
      }
    }

    secondary_calls_per_trait <-
      tidyr::separate_rows(secondary_calls, .data$CLINVAR_UMLS_CUI, sep = ",") %>%
      dplyr::select(.data$VAR_ID, .data$CLINVAR_UMLS_CUI) %>%
      dplyr::left_join(umls_map, by = c("CLINVAR_UMLS_CUI" = "cui")) %>%
      dplyr::rename(CLINVAR_PHENOTYPE = .data$cui_name) %>%
      dplyr::group_by(.data$VAR_ID) %>%
      dplyr::summarise(CLINVAR_PHENOTYPE =
                         paste(unique(.data$CLINVAR_PHENOTYPE), collapse="; ")) %>%
      dplyr::distinct()

    secondary_calls <-
      dplyr::inner_join(secondary_calls,
                        dplyr::select(secondary_calls_per_trait,
                                      .data$VAR_ID, .data$CLINVAR_PHENOTYPE),
                        by = c("VAR_ID" = "VAR_ID"))

    log4r_info(paste0(
      "Found n = ",
      nrow(secondary_calls),
      " variants in genes recommended for return as ",
      "secondary findings from clinical sequencing"))

  }

  return(secondary_calls)

}

#' Function that retrieves variants in cancer predisposition genes linked
#' to cancer-related conditions according to ClinVar
#'
#' @param cpg_calls data frame with variants in predisposition_genes
#' @param oncotree data frame with hereditary cancer phenotypes from OncoTree
#' @param umls_map data frame with UMLS phenotype entries
#'
#' @export
detect_cancer_traits_clinvar <- function(cpg_calls, oncotree, umls_map) {
  if (nrow(cpg_calls) > 0 & "CLINVAR_UMLS_CUI" %in% colnames(cpg_calls)
      & "VAR_ID" %in% colnames(cpg_calls)) {

    oncotree <- oncotree %>%
      dplyr::select(.data$cui) %>%
      dplyr::mutate(CANCER_PHENOTYPE = 1) %>%
      dplyr::distinct()

    n_clinvar <- cpg_calls %>% dplyr::filter(!is.na(.data$CLINVAR_UMLS_CUI)) %>%
      nrow()

    if (n_clinvar > 0) {
      cpg_calls_traits <- as.data.frame(
        tidyr::separate_rows(cpg_calls, .data$CLINVAR_UMLS_CUI, sep = ",") %>%
        dplyr::select(.data$VAR_ID, .data$CLINVAR_UMLS_CUI) %>%
        dplyr::left_join(umls_map, by = c("CLINVAR_UMLS_CUI" = "cui")) %>%
        dplyr::distinct() %>%
        dplyr::filter(!is.na(.data$cui_name)) %>%
        dplyr::left_join(oncotree, by = c("CLINVAR_UMLS_CUI" = "cui")) %>%
        dplyr::mutate(
          CANCER_PHENOTYPE = dplyr::if_else(is.na(.data$CANCER_PHENOTYPE),
                                            as.integer(0),
                                            as.integer(.data$CANCER_PHENOTYPE))) %>%
        dplyr::mutate(
          CANCER_PHENOTYPE =
            dplyr::if_else(stringr::str_detect(tolower(.data$cui_name),
                                               pcgrr::cancer_phenotypes_regex),
                           as.integer(1),
                           as.integer(.data$CANCER_PHENOTYPE))) %>%
        dplyr::group_by(.data$VAR_ID) %>%
        dplyr::summarise(
          CLINVAR_PHENOTYPE = paste(unique(.data$cui_name), collapse = "; "),
          CANCER_PHENOTYPE = max(.data$CANCER_PHENOTYPE),
          .groups = "drop") %>%
          dplyr::mutate(
            CANCER_PHENOTYPE =
              dplyr::if_else(
                stringr::str_detect(
                  .data$CLINVAR_PHENOTYPE,
                  "^(not specified; not provided|not specified|not provided)"),
                as.integer(1),
                as.integer(.data$CANCER_PHENOTYPE)))

      )

      cpg_calls <- cpg_calls %>%
        dplyr::left_join(cpg_calls_traits, by = "VAR_ID")
    }else{
      cpg_calls$CLINVAR_PHENOTYPE <- NA
      cpg_calls$CANCER_PHENOTYPE <- NA
    }


  }
  return(cpg_calls)
}

#' Function that makes a piechart showing the number of variants at
#' each significance level
#'
#' @param variants_tsv data frame with variants in predisposition_genes
#' @param plot_type ClinVar or Other
#'
#' @export
summary_donut_chart <- function(variants_tsv, plot_type = "ClinVar") {

  title <- "ClinVar variants"
  p <- NULL

  if (nrow(variants_tsv) > 0) {

    set_clinvar <- variants_tsv %>%
      dplyr::filter(!is.na(.data$CLINVAR_CLASSIFICATION) &
                      !(.data$CLINVAR_CLASSIFICATION == "NA"))
    set_other <- variants_tsv %>%
      dplyr::filter(nchar(.data$CPSR_CLASSIFICATION) > 0 &
                      (is.na(.data$CLINVAR_CLASSIFICATION) |
                         .data$CLINVAR_CLASSIFICATION == "NA"))

    if ((plot_type == "ClinVar" & nrow(set_clinvar) > 0) |
        (plot_type != "ClinVar" & nrow(set_other) > 0)) {

      m <- data.frame()

      if (plot_type == "ClinVar") {
        if (nrow(set_clinvar) > 0) {
          t <- paste0("n = ", nrow(set_clinvar))
          title <- bquote("ClinVar variants, "~bold(.(t)))
          m <- as.data.frame(set_clinvar %>%
            dplyr::group_by(.data$CLINVAR_CLASSIFICATION) %>%
            dplyr::summarise(n = dplyr::n()) %>%
            dplyr::rename(level = .data$CLINVAR_CLASSIFICATION)) %>%
            dplyr::mutate(
              level =
                factor(
                  .data$level, levels =
                    pcgrr::color_palette[["pathogenicity"]][["levels"]])) %>%
            dplyr::arrange(.data$level) %>%
            dplyr::mutate(prop = as.numeric(.data$n / sum(.data$n))) %>%
            dplyr::mutate(lab.ypos = cumsum(.data$prop) - 0.5 * .data$prop) %>%
            dplyr::mutate(n = as.character(.data$n))
        }
      }else{
        if (nrow(set_other) > 0) {
          t <- paste0("n = ", nrow(set_other))
          title <- bquote("Other variants, CPSR-classified, "~bold(.(t)))
          m <- as.data.frame(set_other %>%
            dplyr::group_by(.data$CPSR_CLASSIFICATION) %>%
            dplyr::summarise(n = dplyr::n()) %>%
            dplyr::rename(level = .data$CPSR_CLASSIFICATION)) %>%
            dplyr::mutate(
              level = factor(
                .data$level,
                levels = pcgrr::color_palette[["pathogenicity"]][["levels"]])) %>%
            dplyr::arrange(.data$level) %>%
            dplyr::mutate(prop = as.numeric(.data$n / sum(.data$n))) %>%
            dplyr::mutate(lab.ypos = cumsum(.data$prop) - 0.5 * .data$prop) %>%
            dplyr::mutate(n = as.character(.data$n))

        }
      }


      p <- ggplot2::ggplot(m, ggplot2::aes(x = 2, y = .data$prop, fill = .data$level)) +
        ggplot2::geom_bar(stat = "identity", color = "white") +
        ggplot2::coord_polar(theta = "y", start = 0) +
        ggplot2::geom_text(ggplot2::aes(y = 1 - .data$lab.ypos, label = .data$n),
                           color = "white", family = "Helvetica", size = 6) +
        ggplot2::scale_fill_manual(
          values = pcgrr::color_palette[["pathogenicity"]][["values"]],
          labels = pcgrr::color_palette[["pathogenicity"]][["levels"]],
          drop = F) +
        ggplot2::theme_void() +
        ggplot2::xlim(0.5, 2.5) +
        ggplot2::ggtitle(title) +
        ggplot2::theme(plot.title =
                         ggplot2::element_text(family = "Helvetica",
                                               size = 16, vjust = -1,
                                               hjust = 0.5),
                       legend.title = ggplot2::element_blank(),
                       plot.margin = ggplot2::margin(0, 0, 0, 0, "cm"),
                       legend.text = ggplot2::element_text(
                         family = "Helvetica", size = 14))
    }


  }
  return(p)

}


#' Function that makes a HTML display of virtual gene panel
#'
#' @param gene_df data frame with genes targeted in virtual panel
#'
#'
#' @export
virtual_panel_display_html <- function(gene_df) {

  i <- 1
  gene_df <- gene_df %>%
    dplyr::arrange(dplyr::desc(.data$confidence_level), symbol)

  html_string <- "<div id=\"container\">"
  while(i <= nrow(gene_df)) {
    confidence_level <- gene_df[i,"confidence_level"]
    css_class <- "exploratory"
    if (confidence_level == 3) {
      css_class <- "green"
    }
    if (confidence_level == 2) {
      css_class <- "amber"
    }
    if (confidence_level == 1) {
      css_class <- "red"
    }
    if (confidence_level == -1) {
      css_class <- "custom"
    }
    if (confidence_level == 0) {
      css_class <- "nolist"
    }
    if (confidence_level == 5) {
      css_class <- "app_combo"
    }
    symbol <- gene_df[i, "symbol"]
    name <- gene_df[i, "genename"]
    entrezgene <- gene_df[i, "entrezgene"]
    panel_id <- gene_df[i, "panel_id"]

    gene_url <- paste0("https://www.ncbi.nlm.nih.gov/gene/", entrezgene)
    if (!is.na(panel_id)) {
      gene_url <- paste0("https://panelapp.genomicsengland.co.uk/panels/",
                         panel_id, "/", symbol)
    }

    entry_string <- paste0("  <div class=\"", css_class, "\"><a href=\"",
                           gene_url, "\" target=\"_blank\" title=\"",
                           name, "\">", symbol, "</a></div>")
    html_string <- paste0(html_string, entry_string)
    if (i %% 7 == 0) {
      html_string <- paste0(html_string, "</div>  <div id=\"container\">")
    }
    i <- i + 1
  }
  html_string <- paste0(html_string, "</div>")
  return(html_string)
}


#' Function that retrieves clinical evidence items (CIVIC, CGI) for
#'somatic cancer variants
#'
#' @param sample_calls data frame with germline variant callset from
#' query sample
#' @param colset vector with column names to display for each report element
#' @param eitems all clinical evidence items linking germline variants with
#' impact on therapeutic response/diagnosis/prognosis etc
#'
#' @return list
#' @export
get_germline_biomarkers <- function(sample_calls,
                                    colset = NULL,
                                    eitems = NULL) {

  log4r_info(paste0("Matching variant set with existing genomic ",
                    "biomarkers from CIViC (germline)"))

  clin_eitems_list <- list()
  for (type in c("diagnostic", "prognostic", "predictive", "predisposing")) {
    clin_eitems_list[[type]] <- data.frame()
  }

  all_var_evidence_items <- data.frame()
  var_eitems <- list()
  for (m in c("codon", "exon", "gene", "exact")) {
    var_eitems[[m]] <- data.frame()
  }

  sample_calls_p_lp <- sample_calls %>%
    dplyr::filter(!is.na(.data$CPSR_CLASSIFICATION) |
                  !is.na(.data$CLINVAR_CLASSIFICATION)) %>%
    dplyr::filter(stringr::str_detect(.data$CPSR_CLASSIFICATION, "Pathogenic") |
                  stringr::str_detect(.data$CLINVAR_CLASSIFICATION, "Pathogenic"))


  ## match clinical evidence items against
  ## query variants (non-regional - exact), civic + cgi
  var_eitems[["exact"]] <-
    pcgrr::match_eitems_to_var(
      sample_calls,
      db = "civic",
      colset = colset,
      eitems = eitems,
      region_marker = F)


  ## For evidence items reported at "regional level" - codon/exon/gene
  ## consider only pathogenic/likely pathogenic variants in the sample
  ## as candidates for match
  var_eitems_regional <- data.frame()
  if(nrow(sample_calls_p_lp) > 0){
    var_eitems_regional <-
      pcgrr::match_eitems_to_var(
        sample_calls_p_lp,
        db = "civic",
        colset = colset,
        eitems = eitems,
        region_marker = T)
  }

  ## for regional biomarkers - perform additional quality checks
  ## (making sure variants are of correct consequence,
  ## at the correct codon/exon etc), and is loss-of-function if
  ## this is specified
  for (m in c("codon", "exon", "gene")) {
    if (nrow(var_eitems_regional) > 0) {
      var_eitems[[m]] <-
        pcgrr::qc_var_eitems(var_eitems = var_eitems_regional,
                             marker_type = m)
    }
  }

  var_eitems <- pcgrr::deduplicate_eitems(var_eitems = var_eitems,
                                          target_type = "exact",
                                          target_other =
                                            c("codon", "exon", "gene"))

  var_eitems <- pcgrr::deduplicate_eitems(var_eitems = var_eitems,
                                          target_type = "codon",
                                          target_other =
                                            c("exon", "gene"))

  ## log the types and number of clinical
  ## evidence items found (exact / codon / exon)
  pcgrr::log_var_eitem_stats(var_eitems = var_eitems,
                             target_type = "exact")
  pcgrr::log_var_eitem_stats(var_eitems = var_eitems,
                             target_type = "codon")
  pcgrr::log_var_eitem_stats(var_eitems = var_eitems,
                             target_type = "exon")
  pcgrr::log_var_eitem_stats(var_eitems = var_eitems,
                             target_type = "gene")

  all_var_evidence_items <- all_var_evidence_items %>%
    dplyr::bind_rows(var_eitems[["exact"]]) %>%
    dplyr::bind_rows(var_eitems[["codon"]]) %>%
    dplyr::bind_rows(var_eitems[["exon"]]) %>%
    dplyr::bind_rows(var_eitems[["gene"]])


  ## Organize all variants in a list object 'clin_items', organized through
  ## evidence type (diagnostic|prognostic|predictive|predisposing)

  if (nrow(all_var_evidence_items) > 0) {
    for (type in c("prognostic", "diagnostic", "predictive", "predisposing")) {
      clin_eitems_list[[type]] <- all_var_evidence_items %>%
        dplyr::filter(.data$EVIDENCE_TYPE == stringr::str_to_title(.data$type)) %>%
        dplyr::arrange(.data$EVIDENCE_LEVEL, .data$RATING) %>%
        dplyr::select(.data$SYMBOL, .data$GENE_NAME, .data$CANCER_TYPE, .data$CLINICAL_SIGNIFICANCE,
                      .data$EVIDENCE_LEVEL, .data$RATING, .data$EVIDENCE_DIRECTION, .data$CITATION,
                      .data$THERAPEUTIC_CONTEXT, .data$EVIDENCE_TYPE, .data$DESCRIPTION, .data$BIOMARKER_MAPPING,
                      .data$CDS_CHANGE, .data$LOSS_OF_FUNCTION, .data$GENOMIC_CHANGE) %>%
        dplyr::distinct()
    }
  }

  return(clin_eitems_list)

}

#' Function that reports protein-coding geneset that overlaps BED file
#'
#' @param bed_file BED file name with selected transcripts from panel 0
#' @param pcgr_data object with PCGR annotation data
#' @export
custom_bed_genes <- function(bed_file, pcgr_data) {

  invisible(assertthat::assert_that(file.exists(bed_file),
                                    msg = paste0("BED file", bed_file, " does not exist")))
  bed_df <- as.data.frame(
    utils::read.table(file = bed_file, header = F, stringsAsFactors = F,
                           comment.char = "", quote = "", sep = "\t") %>%
    magrittr::set_colnames(c("chromosome", "segment_start", "segment_end", "onco_xref")) %>%
    dplyr::rowwise() %>%
    dplyr::mutate(symbol = unlist(strsplit(.data$onco_xref, "\\|"))[4]) %>%
    dplyr::select(.data$symbol) %>%
    dplyr::left_join(dplyr::select(pcgr_data$gene_xref$gencode, .data$symbol, .data$ENTREZ_ID, .data$GENENAME), by = "symbol") %>%
    dplyr::rename(genename = .data$GENENAME, entrezgene = .data$ENTREZ_ID) %>%
    dplyr::select(.data$symbol, .data$genename, .data$entrezgene) %>%
    dplyr::distinct()
  )
  return(bed_df)
}

