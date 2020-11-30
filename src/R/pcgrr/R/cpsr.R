## Function that generates predisposition_report - CPSR
#'
#' @param project_directory name of project directory
#' @param query_vcf2tsv name of gzipped TSV file (vcf2tsv)
#' with annotated query SNVs/InDels
#' @param custom_bed BED file with custom panels
#' @param pcgr_data List of data frames with PCGR data annotations
#' @param pcgr_version PCGR software version
#' @param cpsr_version CPSR software version
#' @param cpsr_config Object with CPSR configuration parameters
#' @param virtual_panel_id Identifier for virtual panel
#' @param sample_name sample identifier
#'

options(error = traceback)

generate_cpsr_report <- function(project_directory, query_vcf2tsv,
                                           custom_bed, pcgr_data,
                                           pcgr_version, cpsr_version,
                                           cpsr_config = NULL,
                                           virtual_panel_id = -1,
                                           sample_name = "SampleX") {

  invisible(assertthat::assert_that(
    !is.null(cpsr_config),
    msg = "Object 'cpsr_config' cannot be NULL"))
  cps_report <-
    pcgrr::init_report(
      cpsr_config, sample_name, class = NULL,
      pcgr_data = pcgr_data, pcgr_version = pcgr_version,
      cpsr_version = cpsr_version, type = "germline",
      virtual_panel_id = virtual_panel_id,
      custom_bed = custom_bed)
  if (is.null(cps_report$metadata$gene_panel)) {
    return(NULL)
  }

  class_1_5_display <-
    c("SYMBOL", "SOURCE", "CLINVAR_PHENOTYPE", "CONSEQUENCE",
      "PROTEIN_CHANGE", "GENOTYPE", "GENE_NAME", "PROTEIN_DOMAIN",
      "HGVSp", "HGVSc", "NCBI_REFSEQ", "CDS_CHANGE",
      "REFSEQ_MRNA", "PROB_GNOMAD_LOF_INTOLERANT",
      "PROB_GNOMAD_LOF_INTOLERANT_HOM", "PROB_GNOMAD_LOF_TOLERANT_NULL",
      "MUTATION_HOTSPOT", "RMSK_HIT", "PROTEIN_FEATURE", "PREDICTED_EFFECT",
      "GERP_DBNSFP",
      "LOSS_OF_FUNCTION", "DBSNP", "CLINVAR", "CLINVAR_CLASSIFICATION",
      "CLINVAR_REVIEW_STATUS_STARS", "CLINVAR_CONFLICTED",
      "CLINVAR_VARIANT_ORIGIN",
      "CPSR_CLASSIFICATION", "CPSR_PATHOGENICITY_SCORE",
      "CPSR_CLASSIFICATION_DOC", "CPSR_CLASSIFICATION_CODE", "ONCOGENE",
      "TUMOR_SUPPRESSOR", "GLOBAL_AF_GNOMAD",
      cps_report[["metadata"]][["config"]][["popgen"]][["vcftag_gnomad"]],
      "GENOMIC_CHANGE", "GENOME_VERSION")

  ## define tags/variables to display in data tables (secondary findings)
  secondary_findings_display <-
    c("SYMBOL", "CONSEQUENCE", "CLINVAR_CLASSIFICATION", "CLINVAR_PHENOTYPE",
      "PROTEIN_CHANGE", "GENOTYPE", "GENE_NAME", "PROTEIN_DOMAIN",
      "HGVSp", "HGVSc", "NCBI_REFSEQ", "CDS_CHANGE", "REFSEQ_MRNA",
      "MUTATION_HOTSPOT", "RMSK_HIT", "PROTEIN_FEATURE", "PREDICTED_EFFECT",
      "LOSS_OF_FUNCTION", "DBSNP", "CLINVAR", "CLINVAR_REVIEW_STATUS_STARS",
      "CLINVAR_CONFLICTED", "CLINVAR_VARIANT_ORIGIN",
      "GWAS_CITATION", "ONCOGENE", "TUMOR_SUPPRESSOR",
      "GLOBAL_AF_GNOMAD",
      cps_report[["metadata"]][["config"]][["popgen"]][["vcftag_gnomad"]],
      "GENOMIC_CHANGE", "GENOME_VERSION")

  predispose_gwas_display <-
    c("SYMBOL", "CONSEQUENCE", "GWAS_CITATION", "PROTEIN_CHANGE", "GENOTYPE",
      "LOSS_OF_FUNCTION", "PROTEIN_CHANGE", "GENE_NAME", "GWAS_PHENOTYPE",
      "PROTEIN_DOMAIN", "HGVSp", "HGVSc", "NCBI_REFSEQ",
      "CDS_CHANGE", "CODING_STATUS",
      "REFSEQ_MRNA", "PROTEIN_FEATURE",
      "PREDICTED_EFFECT", "DBSNP", "GLOBAL_AF_GNOMAD",
      "GENOMIC_CHANGE", "GENOME_VERSION")

  if (query_vcf2tsv != "None.gz") {
    if (!file.exists(query_vcf2tsv) | file.size(query_vcf2tsv) == 0) {
      rlogging::warning(
        paste0("File ",
               query_vcf2tsv,
               " does not exist or has zero size"))
    }else{
      if (!is.null(cpsr_config) & query_vcf2tsv != "None.gz") {

        ## read calls
        calls <-
          pcgrr::get_calls(
            query_vcf2tsv, pcgr_data, sample_name,
            cpsr_config,
            oncotree =
              cps_report[["metadata"]][["phenotype_ontology"]][["oncotree_query"]],
            cpsr = TRUE)
        if (nrow(calls) == 0) {
          rlogging::warning("There are zero calls in input file ",
                            "- no report will be produced")
          return(NULL)
        }
        if (cpsr_config[["ignore_noncoding"]] == T) {
          n_noncoding_vars <- calls %>%
            dplyr::filter(CODING_STATUS == "noncoding") %>%
            NROW()
          rlogging::message(
            "Excluding n = ",
            n_noncoding_vars,
            " variants for classification (option --ignore_noncoding)")
          calls <- calls %>%
            dplyr::filter(CODING_STATUS == "coding")
          if (NROW(calls) == 0) {
            rlogging::warning(
              "There are zero remaining protein-coding ",
              "calls in input file - no report will be produced")
            return(NULL)
          }

        }
        ## get overall call statistics
        call_stats <-
          pcgrr::variant_stats_report(calls, name = "v_stat")

        calls <- dplyr::mutate(calls, SYMBOL = as.character(SYMBOL))
        cpg_calls <-
          dplyr::inner_join(calls,
                            cps_report[["metadata"]][["gene_panel"]][["genes"]],
                            by = c("SYMBOL" = "symbol"))
        cpg_call_stats <-
          pcgrr::variant_stats_report(
            cpg_calls, name = "v_stat_cpg")

        rlogging::message(
          paste0("Number of coding variants in cancer predisposition genes: ",
                 cpg_call_stats[["v_stat_cpg"]][["n_coding"]]))
        rlogging::message(
          paste0(
            "Number of non-coding variants in cancer predisposition genes: ",
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
            dplyr::filter(!is.na(CLINVAR_MSID)) %>%
            dplyr::filter(CANCER_PHENOTYPE == 0) %>%
            NROW()

          rlogging::message(
            "ClinVar variants related to non-cancer conditions excluded",
            " from report: ", n_clinvar_noncancer)

          if (n_clinvar_noncancer > 0) {
            cpg_calls_exclude <- cpg_calls %>%
              dplyr::filter(!is.na(CLINVAR_MSID)) %>%
              dplyr::filter(CANCER_PHENOTYPE == 0) %>%
              dplyr::select(VAR_ID)

            cpg_calls <- cpg_calls %>%
              dplyr::anti_join(cpg_calls_exclude, by = "VAR_ID")

            rlogging::message(
              "Variants remaining after exclusion of non-cancer related",
              " ClinVar variants: ", NROW(cpg_calls))

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
          pcgrr::assign_cpsr_tier(cpg_calls,
                                  cps_report[["metadata"]][["config"]],
                                  class_1_5_display)
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
          pcgrr::remove_cols_from_df(cnames = c("CIVIC_ID", "CIVIC_ID_SEGMENT",
                                                "AMINO_ACID_END","AMINO_ACID_START",
                                                "EXON"))

        gene_hits <- paste(
          unique(
            cps_report[["content"]][["snv_indel"]][["variant_set"]][["tsv"]]$SYMBOL),
          collapse = ", ")
        rlogging::message("Variants were found in the following cancer ",
                          "predisposition genes: ", gene_hits)

        ## secondary findings
        if (cpsr_config[["secondary_findings"]] == TRUE) {
          if (nrow(secondary_calls) > 0) {
            rlogging::message("Assignment of other variants in genes ",
                              "recommended for reporting as secondary ",
                              "findings (ACMG SF v2.0)")
            cps_report[["content"]][["snv_indel"]][["disp"]][["sf"]] <-
              secondary_calls %>%
              dplyr::arrange(LOSS_OF_FUNCTION, CODING_STATUS) %>%
              dplyr::select(dplyr::one_of(secondary_findings_display))
            rlogging::message(
              "Number of pathogenic variants in the secondaryome - other ",
              "genes of clinical significance: ",
              cps_report[["content"]][["snv_indel"]][["v_stat_secondary"]][["n_coding"]])
          }
        }

        cps_report[["content"]][["snv_indel"]][["eval"]] <- TRUE

        if (cpsr_config[["gwas_findings"]] == TRUE) {
          rlogging::message("Assignment of other variants to hits ",
                            "from genome-wide association studies")

          ## Assign GWAS hits to cps_report object
          cps_report[["content"]][["snv_indel"]][["disp"]][["gwas"]] <-
            dplyr::filter(calls, !is.na(GWAS_HIT) & !is.na(GWAS_CITATION))
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
                dplyr::select(dplyr::one_of(predispose_gwas_display)) %>%
                dplyr::arrange(LOSS_OF_FUNCTION, CODING_STATUS)
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
          dplyr::filter(SOURCE == "ClinVar") %>%
          nrow()
        num_rows_other <- t1 %>%
          dplyr::filter(SOURCE == "Other") %>%
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
get_insilico_prediction_statistics <- function(cpg_calls) {

  insilico_pathogenicity_pred_algos <-
    c("SIFT_DBNSFP", "PROVEAN_DBNSFP",
      "META_LR_DBNSFP", "FATHMM_DBNSFP",
      "MUTATIONTASTER_DBNSFP", "DEOGEN2_DBNSFP",
      "PRIMATEAI_DBNSFP", "MUTATIONASSESSOR_DBNSFP",
      "FATHMM_MKL_DBNSFP", "M_CAP_DBNSFP",
      "CLINPRED_DBNSFP", "LIST_S2_DBNSFP",
      "BAYESDEL_ADDAF_DBNSFP",
      "SPLICE_SITE_ADA_DBNSFP", "SPLICE_SITE_RF_DBNSFP")
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
                  (cpg_calls[, algo] == "damaging" |
                     cpg_calls[, algo] == "possibly_damaging"),
                "N_INSILICO_DAMAGING"] <-
        cpg_calls[!is.na(cpg_calls[, algo]) &
                    (cpg_calls[, algo] == "damaging" |
                       cpg_calls[, algo] == "possibly_damaging"),
                  "N_INSILICO_DAMAGING"] + 1
      cpg_calls[!is.na(cpg_calls[, algo]) &
                  cpg_calls[, algo] == "tolerated",
                "N_INSILICO_TOLERATED"] <-
        cpg_calls[!is.na(cpg_calls[, algo]) &
                    cpg_calls[, algo] == "tolerated",
                  "N_INSILICO_TOLERATED"] + 1
      cpg_calls[!is.na(cpg_calls[, algo]) &
                  cpg_calls[, algo] == "affect_splicing",
                "N_INSILICO_SPLICING_AFFECTED"] <-
        cpg_calls[!is.na(cpg_calls[, algo]) &
                    cpg_calls[, algo] == "affect_splicing",
                  "N_INSILICO_SPLICING_AFFECTED"] + 1
      cpg_calls[!is.na(cpg_calls[, algo]) &
                  cpg_calls[, algo] == "splicing_neutral",
                "N_INSILICO_SPLICING_NEUTRAL"] <-
        cpg_calls[!is.na(cpg_calls[, algo]) &
                    cpg_calls[, algo] == "splicing_neutral",
                  "N_INSILICO_SPLICING_NEUTRAL"] + 1
    }
  }
  return(cpg_calls)
}

#' Function that assigns variant pathogenicity evidence based on ACMG guidelines
#'
#' @param cpg_calls sample calls with dbnsfp annotations
#' @param cpsr_config pcgr configuration object
#' @param pcgr_data pcgr data object
#'
#' @return cpg_calls
#'
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
        dplyr::if_else(stringr::str_detect(CANCER_PREDISPOSITION_MOI,
                                           "AD|AD/AR"),
                       "AD", as.character(NA), as.character(NA))) %>%
    dplyr::mutate(
      cpsr_gene_moi =
        dplyr::if_else(!stringr::str_detect(CANCER_PREDISPOSITION_MOI,
                                            "AD|Mosaic"),
                       "AR", cpsr_gene_moi, cpsr_gene_moi))

  predisposition_gene_info <-
    dplyr::select(pcgr_data[["predisposition"]][["genes"]], symbol,
                  mechanism_of_disease, path_truncation_rate,
                  benign_missense_rate) %>%
    dplyr::rename(SYMBOL = symbol,
                  MOD = mechanism_of_disease,
                  PATH_TRUNCATION_RATE = path_truncation_rate,
                  BENIGN_MISSENSE_RATE = benign_missense_rate)

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
  # MutationAssessor,M_CAP,MutPred,FATHMM,FATHMM-mkl,DBNSFP_LogReg,dbscSNV_RF,
  # dbscSNV_AdaBoost
  # Default scheme (from default TOML file):
  # 1) Damaging: Among all possible protein variant effect predictions, at
  #              least six algorithms must have made a call,
  #              with at least 5 predicted as damaging, and at most two
  #              predicted as tolerated (PP3)
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
        dplyr::if_else(N_INSILICO_CALLED >= dbnsfp_min_called &
                         N_INSILICO_DAMAGING >= dbnsfp_min_majority &
                         N_INSILICO_TOLERATED <= dbnsfp_max_minority &
                         N_INSILICO_SPLICING_NEUTRAL <= 1, TRUE,
                       FALSE, FALSE)) %>%
    dplyr::mutate(
      ACMG_BP4 = dplyr::if_else(N_INSILICO_CALLED >= dbnsfp_min_called &
                                  N_INSILICO_TOLERATED >= dbnsfp_min_majority &
                                  N_INSILICO_DAMAGING <= dbnsfp_max_minority &
                                  N_INSILICO_SPLICING_AFFECTED == 0, TRUE,
                                FALSE, FALSE)) %>%
    dplyr::mutate(ACMG_PP3 = dplyr::case_when(
      N_INSILICO_SPLICING_AFFECTED == 2 ~ TRUE, TRUE ~ as.logical(ACMG_PP3)))

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
                           gad_af <= pathogenic_range_af,
                         TRUE, FALSE, FALSE)) %>%
      dplyr::mutate(
        ACMG_PM2_2 = dplyr::if_else(is.na(!!rlang::sym(gad_AC_tag)),
                                    TRUE, FALSE, FALSE)) %>%
      dplyr::mutate(
        ACMG_BA1_AD = dplyr::if_else(ACMG_PM2_2 == FALSE &
                                       gad_af >= 0.005 &
                                       cpsr_gene_moi == "AD",
                                     TRUE, FALSE, FALSE)) %>%
      dplyr::mutate(
        ACMG_BS1_1_AD = dplyr::if_else(ACMG_BA1_AD == FALSE &
                                         ACMG_PM2_2 == FALSE &
                                         gad_af >= 0.001 &
                                         cpsr_gene_moi == "AD",
                                       TRUE, FALSE, FALSE)) %>%
      dplyr::mutate(
        ACMG_BS1_2_AD = dplyr::if_else(ACMG_BS1_1_AD == FALSE &
                                         ACMG_BA1_AD == FALSE &
                                         ACMG_PM2_2 == FALSE &
                                         gad_af > pathogenic_range_af &
                                         cpsr_gene_moi == "AD",
                                       TRUE, FALSE, FALSE)) %>%
      dplyr::mutate(
        ACMG_BA1_AR = dplyr::if_else(ACMG_PM2_2 == FALSE &
                                       gad_af >= 0.01 &
                                       (cpsr_gene_moi == "AR" |
                                          is.na(cpsr_gene_moi)),
                                     TRUE, FALSE, FALSE)) %>%
      dplyr::mutate(
        ACMG_BS1_1_AR = dplyr::if_else(ACMG_BA1_AR == FALSE &
                                         ACMG_PM2_2 == FALSE &
                                         gad_af >= 0.003 &
                                         (cpsr_gene_moi == "AR" |
                                            is.na(cpsr_gene_moi)),
                                       TRUE, FALSE, FALSE)) %>%
      dplyr::mutate(
        ACMG_BS1_2_AR = dplyr::if_else(ACMG_BA1_AR == FALSE &
                                         ACMG_BS1_1_AR == FALSE &
                                         ACMG_PM2_2 == FALSE &
                                         gad_af > pathogenic_range_af &
                                         (cpsr_gene_moi == "AR" |
                                            is.na(cpsr_gene_moi)),
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
        dplyr::if_else(NULL_VARIANT == T & LOSS_OF_FUNCTION == T &
                         MOD == "LoF" &
                         (ACMG_PM2_1 == TRUE | ACMG_PM2_2 == TRUE),
                       TRUE, FALSE, FALSE)) %>%
    dplyr::mutate(
      ACMG_PVS1_3 =
        dplyr::if_else(NULL_VARIANT == T & LOSS_OF_FUNCTION == T &
                         (is.na(MOD) | MOD != "LoF") &
                         (ACMG_PM2_1 == TRUE | ACMG_PM2_2 == TRUE),
                       TRUE, FALSE, FALSE)) %>%
    dplyr::mutate(
      ACMG_PVS1_2 =
        dplyr::if_else(NULL_VARIANT == T & LOSS_OF_FUNCTION == F &
                         MOD == "LoF" &
                         (ACMG_PM2_1 == TRUE | ACMG_PM2_2 == TRUE),
                       TRUE, FALSE, FALSE)) %>%
    dplyr::mutate(
      ACMG_PVS1_4 =
        dplyr::if_else(NULL_VARIANT == T & LOSS_OF_FUNCTION == F &
                         (is.na(MOD) | MOD != "LoF") &
                         (ACMG_PM2_1 == TRUE | ACMG_PM2_2 == TRUE),
                       TRUE, FALSE, FALSE)) %>%
    dplyr::mutate(
      ACMG_PVS1_5 =
        dplyr::if_else(CONSEQUENCE == "start_lost" & MOD == "LoF" &
                         (ACMG_PM2_1 == TRUE | ACMG_PM2_2 == TRUE),
                       TRUE, FALSE, FALSE)) %>%
    dplyr::mutate(
      ACMG_PVS1_6 =
        dplyr::if_else(CONSEQUENCE == "start_lost" & (is.na(MOD) |
                                                        MOD != "LoF") &
                         (ACMG_PM2_1 == TRUE | ACMG_PM2_2 == TRUE),
                       TRUE, FALSE, FALSE)) %>%
    dplyr::mutate(
      ACMG_PVS1_7 =
        dplyr::if_else(LOSS_OF_FUNCTION == T &
                         stringr::str_detect(CONSEQUENCE, "_donor|_acceptor") &
                         LAST_INTRON == F & MOD == "LoF" &
                         (ACMG_PM2_1 == TRUE | ACMG_PM2_2 == TRUE),
                       TRUE, FALSE, FALSE)) %>%
    dplyr::mutate(
      ACMG_PVS1_8 =
        dplyr::if_else(LOSS_OF_FUNCTION == T &
                         stringr::str_detect(CONSEQUENCE, "_donor|_acceptor") &
                         LAST_INTRON == T & MOD == "LoF" &
                         (ACMG_PM2_1 == TRUE | ACMG_PM2_2 == TRUE),
                       TRUE, FALSE, FALSE)) %>%
    dplyr::mutate(
      ACMG_PVS1_9 =
        dplyr::if_else(LOSS_OF_FUNCTION == T &
                         stringr::str_detect(CONSEQUENCE, "_donor|_acceptor") &
                         LAST_INTRON == F & (is.na(MOD) | MOD != "LoF") &
                         (ACMG_PM2_1 == TRUE | ACMG_PM2_2 == TRUE),
                       TRUE, FALSE, FALSE)) %>%
    dplyr::mutate(
      ACMG_PVS1_10 =
        dplyr::if_else(SPLICE_DONOR_RELEVANT == T & ACMG_PP3 == TRUE &
                         (ACMG_PM2_1 == TRUE | ACMG_PM2_2 == TRUE),
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
              CONSEQUENCE, "stop_lost|inframe_deletion|inframe_insertion") &
              is.na(RMSK_HIT) & cpsr_gene_moi == "AD",
            TRUE, FALSE, FALSE)) %>%
      dplyr::mutate(
        ACMG_PPC1 =
          dplyr::if_else(
            stringr::str_detect(
              CONSEQUENCE, "stop_lost|inframe_deletion|inframe_insertion") &
              is.na(RMSK_HIT) & (cpsr_gene_moi == "AR" | is.na(cpsr_gene_moi)),
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
          (is.na(BENIGN_MISSENSE_RATE) | BENIGN_MISSENSE_RATE <= 0.1) &
            (is.na(PATH_TRUNCATION_RATE) | PATH_TRUNCATION_RATE < 0.5) &
            stringr::str_detect(CONSEQUENCE, "^missense_variant"),
          TRUE, FALSE, FALSE))

  ## Assign a logical ACMG evidence indicator
  # ACMG_BP1 - Missense variant in a gene for which primarily truncating
  # variants (> 90%, as given in Maxwell et al.) are known to cause disease
  cpg_calls <- cpg_calls %>%
    dplyr::mutate(
      ACMG_BP1 =
        dplyr::if_else(PATH_TRUNCATION_RATE > 0.90 &
                         stringr::str_detect(CONSEQUENCE, "^missense_variant"),
                       TRUE, FALSE, FALSE))

  ## Assign a logical ACMG evidence indicator
  # ACMG_BP7 - Silent/intronic variant outside of the splice site consensus
  cpg_calls <- cpg_calls %>%
    dplyr::mutate(
      ACMG_BP7 =
        dplyr::if_else((
          (INTRON_POSITION < 0 & INTRON_POSITION < -3) |
            (INTRON_POSITION > 0 & INTRON_POSITION > 6) |
            (EXON_POSITION < 0 & EXON_POSITION < -2) |
            (EXON_POSITION > 0 & EXON_POSITION > 1)) &
            stringr::str_detect(
              CONSEQUENCE,
              "^(synonymous_variant|intron_variant|splice_region_variant)"),
          TRUE, FALSE, FALSE))

  ## Assign a logical ACMG evidence indicator
  # ACMG_BP3 - Variants in promoter or untranslated regions
  cpg_calls <- cpg_calls %>%
    dplyr::mutate(
      ACMG_BP3 =
        dplyr::if_else(
          stringr::str_detect(
            CONSEQUENCE,
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
        dplyr::filter(dplyr::select(cpg_calls, VAR_ID, codon_prefix, SYMBOL),
                      !is.na(codon_prefix)),
        pcgr_data[["clinvar"]][["cpg_loci"]][["hi"]][["P"]][["codon"]],
        by = c("codon_prefix" = "codon_prefix", "SYMBOL" = "symbol")) %>%
      dplyr::filter(clinvar_pathogenic_codon == T)
    cpg_calls <- cpg_calls %>%
      dplyr::left_join(dplyr::select(cpg_calls_pathogenic_codon,
                                     VAR_ID, clinvar_pathogenic_codon),
                       by = c("VAR_ID"))

    cpg_calls_benign_codon <-
      dplyr::left_join(
        dplyr::filter(dplyr::select(cpg_calls, VAR_ID, codon_prefix, SYMBOL),
                      !is.na(codon_prefix)),
        pcgr_data[["clinvar"]][["cpg_loci"]][["hi"]][["B"]][["codon"]],
        by = c("codon_prefix" = "codon_prefix", "SYMBOL" = "symbol")) %>%
      dplyr::filter(clinvar_benign_codon == T)

    cpg_calls <- cpg_calls %>%
      dplyr::left_join(
        dplyr::select(cpg_calls_benign_codon, VAR_ID, clinvar_benign_codon),
        by = c("VAR_ID")) %>%
      dplyr::mutate(ACMG_PM5 = dplyr::if_else(
        clinvar_pathogenic_codon == TRUE, TRUE, FALSE, FALSE)) %>%
      dplyr::mutate(ACMG_BMC1 = dplyr::if_else(
        clinvar_benign_codon == TRUE, TRUE, FALSE, FALSE))
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
          dplyr::select(cpg_calls, VAR_ID, HGVSp_short, SYMBOL),
          !is.na(HGVSp_short)),
        pcgr_data[["clinvar"]][["cpg_loci"]][["hi"]][["P"]][["peptide_change"]],
        by = c("HGVSp_short" = "hgvs_p", "SYMBOL" = "symbol")) %>%
      dplyr::filter(clinvar_pathogenic == T)

    cpg_calls <- cpg_calls %>%
      dplyr::left_join(
        dplyr::select(cpg_calls_pathogenic_hgvsp,
                      VAR_ID, clinvar_pathogenic), by = c("VAR_ID"))

    cpg_calls_benign_hgvsp <-
      dplyr::left_join(
        dplyr::filter(
          dplyr::select(cpg_calls, VAR_ID, HGVSp_short, SYMBOL),
          !is.na(HGVSp_short)),
        pcgr_data[["clinvar"]][["cpg_loci"]][["hi"]][["B"]][["peptide_change"]],
        by = c("HGVSp_short" = "hgvs_p", "SYMBOL" = "symbol")) %>%
      dplyr::filter(clinvar_benign == T)

    cpg_calls <- cpg_calls %>%
      dplyr::left_join(
        dplyr::select(cpg_calls_benign_hgvsp, VAR_ID, clinvar_benign),
                       by = c("VAR_ID")) %>%
      dplyr::mutate(ACMG_PS1 =
                      dplyr::if_else(clinvar_pathogenic == TRUE,
                                     TRUE, FALSE, FALSE)) %>%
      dplyr::mutate(ACMG_BSC1 =
                      dplyr::if_else(clinvar_benign == TRUE,
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
    dplyr::mutate(ACMG_PM5 =
                    dplyr::case_when(ACMG_PM5 == T & ACMG_PS1 == T ~ FALSE,
                                              TRUE ~ as.logical(ACMG_PM5))) %>%
    ## if previously found coinciding with benign variant (ACMG_BSC1),
    ##  set ACMG_BMC1 to false
    dplyr::mutate(ACMG_BMC1 =
                    dplyr::case_when(ACMG_BMC1 == T & ACMG_BSC1 == T ~ FALSE,
                                     TRUE ~ as.logical(ACMG_BMC1))) %>%
    dplyr::select(-c(clinvar_pathogenic_codon, clinvar_benign_codon,
                     clinvar_pathogenic,
                     clinvar_benign, cpsr_gene_moi, gad_af))

  ##Assign logical ACMG level
  # PM1 - missense variant in a somatic mutation hotspot as
  # determined by cancerhotspots.org (v2)
  cpg_calls <- cpg_calls %>%
    tidyr::separate(MUTATION_HOTSPOT,
                    c("hotspot_symbol", "hotspot_codon", "hotspot_pvalue"),
                    sep = "\\|", remove = F, extra = "drop") %>%
    dplyr::mutate(hotspot_codon =
                    dplyr::if_else(!is.na(hotspot_codon),
                                                 paste0("p.", hotspot_codon),
                                                 as.character(NA))) %>%
    dplyr::mutate(ACMG_PM1 =
                    dplyr::if_else(!is.na(hotspot_codon) &
                                     !is.na(hotspot_symbol) &
                                     !is.na(codon_prefix) &
                                     SYMBOL == hotspot_symbol &
                                     hotspot_codon == codon_prefix,
                                   TRUE, FALSE))

  return(cpg_calls)
}


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
          cpsr_score_benign = cpsr_score_benign +
            dplyr::if_else(pole == "B" & !!rlang::sym(cpsr_evidence_code) == T,
                           score, 0)) %>%
        dplyr::mutate(
          cpsr_score_pathogenic = cpsr_score_pathogenic +
            dplyr::if_else(pole == "P" & !!rlang::sym(cpsr_evidence_code) == T,
                           score, 0)) %>%
        dplyr::mutate(
          CPSR_CLASSIFICATION_DOC =
            paste0(CPSR_CLASSIFICATION_DOC,
                   dplyr::if_else(!!rlang::sym(cpsr_evidence_code) == T,
                                  paste0("- ", description), ""),
                   sep = "<br>")) %>%
        dplyr::mutate(
          CPSR_CLASSIFICATION_CODE =
            paste0(CPSR_CLASSIFICATION_CODE,
                   dplyr::if_else(!!rlang::sym(cpsr_evidence_code) == T,
                                  cpsr_evidence_code, ""), sep = "|"))
    }
    i <- i + 1
  }

  cpg_calls <- cpg_calls %>%
    dplyr::mutate(
      CPSR_CLASSIFICATION_CODE =
        stringr::str_replace_all(
          stringr::str_replace_all(
            CPSR_CLASSIFICATION_CODE,
            "(\\|{2,})", "|"), "(^\\|)|(\\|$)", "")) %>%
    dplyr::mutate(
      CPSR_CLASSIFICATION_DOC =
        stringr::str_replace_all(
          stringr::str_replace_all(
            CPSR_CLASSIFICATION_DOC,
            "(<br>){2,}", "<br>"), "(^(<br>))|((<br>)$)", "")) %>%

    ## Adjust scores in cases where critera are acting as a
    ## prerequisite for other criteria
    dplyr::mutate(
      cpsr_score_pathogenic =
        dplyr::if_else(
          stringr::str_detect(CPSR_CLASSIFICATION_CODE, "ACMG_PVS") &
            stringr::str_detect(CPSR_CLASSIFICATION_CODE, "ACMG_PM2_2"),
          cpsr_score_pathogenic - 1, cpsr_score_pathogenic)) %>%
    dplyr::mutate(
      cpsr_score_pathogenic =
        dplyr::if_else(
          stringr::str_detect(CPSR_CLASSIFICATION_CODE, "ACMG_PVS") &
            stringr::str_detect(CPSR_CLASSIFICATION_CODE, "ACMG_PM2_1"),
          cpsr_score_pathogenic - 0.5, cpsr_score_pathogenic)) %>%
    dplyr::mutate(
      cpsr_score_pathogenic =
        dplyr::if_else(
          stringr::str_detect(CPSR_CLASSIFICATION_CODE, "ACMG_PVS1_10") &
            stringr::str_detect(CPSR_CLASSIFICATION_CODE, "ACMG_PP3"),
          cpsr_score_pathogenic - 0.5, cpsr_score_pathogenic)) %>%

    ## Add scores accumulated with benign criteria and pathogenic criteria
    dplyr::mutate(CPSR_PATHOGENICITY_SCORE =
                    dplyr::if_else(cpsr_score_benign == 0,
                                   cpsr_score_pathogenic,
                                   cpsr_score_benign)) %>%
    dplyr::mutate(CPSR_PATHOGENICITY_SCORE =
                    dplyr::if_else(cpsr_score_benign < 0 &
                                     cpsr_score_pathogenic > 0,
                                   cpsr_score_benign + cpsr_score_pathogenic,
                                   CPSR_PATHOGENICITY_SCORE)) %>%
    dplyr::mutate(
      CPSR_CLASSIFICATION =
        dplyr::case_when(CPSR_PATHOGENICITY_SCORE <= -1.5 &
                           CPSR_PATHOGENICITY_SCORE >= -4.5 ~ "Likely_Benign",
                         CPSR_PATHOGENICITY_SCORE <= -5 ~ "Benign",
                         CPSR_PATHOGENICITY_SCORE <= 2.0 &
                           CPSR_PATHOGENICITY_SCORE >= -1.0 ~ "VUS",
                         CPSR_PATHOGENICITY_SCORE >= 5 ~ "Pathogenic",
                         CPSR_PATHOGENICITY_SCORE >= 2.5 &
                           CPSR_PATHOGENICITY_SCORE <= 4.5 ~ "Likely_Pathogenic",
                         TRUE ~ as.character("VUS"))) %>%
    dplyr::select(-c(cpsr_score_benign, cpsr_score_pathogenic))

  return(cpg_calls)

}

#' Function that assign variants to different tiers for
#' prioritization of germline variants
#'
#' @param cpg_calls data frame with variants in predisposition_genes
#'
assign_cpsr_tier <- function(cpg_calls, cpsr_config, display_tags) {

  predispose_tsv_tags <-
    c("GENOMIC_CHANGE", "VAR_ID", "GENOTYPE",
      "GENOME_VERSION", "VCF_SAMPLE_ID", "VARIANT_CLASS",
      "CODING_STATUS", "SYMBOL", "GENE_NAME", "CCDS",
      "ENTREZ_ID", "UNIPROT_ID", "ENSEMBL_GENE_ID", "ENSEMBL_TRANSCRIPT_ID",
      "REFSEQ_MRNA", "ONCOGENE", "TUMOR_SUPPRESSOR",
      #"MOD", "CANCER_PREDISPOSITION_MOI",
      "AMINO_ACID_START","AMINO_ACID_END",
      "CONSEQUENCE",
      "VEP_ALL_CSQ", "CIVIC_ID", "CIVIC_ID_SEGMENT", "PROTEIN_CHANGE",
      "PROTEIN_DOMAIN", "HGVSp", "HGVSc", "LAST_EXON",
      "EXON","CODING_STATUS","EXON_POSITION",
      "INTRON_POSITION","CDS_CHANGE","CANCER_PHENOTYPE",
      "MUTATION_HOTSPOT", "RMSK_HIT", "PROTEIN_FEATURE", "EFFECT_PREDICTIONS",
      "LOSS_OF_FUNCTION", "DBSNP", "CLINVAR_CLASSIFICATION", "CLINVAR_MSID",
      "CLINVAR_VARIANT_ORIGIN", "CLINVAR_CONFLICTED", "CLINVAR_PHENOTYPE",
      "CLINVAR_REVIEW_STATUS_STARS", "N_INSILICO_CALLED", "N_INSILICO_DAMAGING",
      "N_INSILICO_TOLERATED", "N_INSILICO_SPLICING_NEUTRAL",
      "N_INSILICO_SPLICING_AFFECTED",
      "GLOBAL_AF_GNOMAD", cpsr_config[["popgen"]][["vcftag_gnomad"]])

  evidence_codes <- pcgrr::cpsr_acmg[['evidence_codes']] %>%
    dplyr::filter(cpsr_evidence_code != "ACMG_BS2_1" &
                    cpsr_evidence_code != "ACMG_BS2_2" &
                    cpsr_evidence_code != "ACMG_BS2_3")
  predispose_tsv_tags <-
    c(predispose_tsv_tags, evidence_codes$cpsr_evidence_code,
      c("CPSR_CLASSIFICATION", "CPSR_PATHOGENICITY_SCORE",
        "CPSR_CLASSIFICATION_CODE", "CPSR_CLASSIFICATION_DOC", "SOURCE"))

  ## Add custom annotation tags to lists of tags to display
  if (!is.null(cpsr_config[["custom_tags"]])) {
    if (cpsr_config[["custom_tags"]][["custom_tags"]] != "") {
      tags <- stringr::str_split(cpsr_config[["custom_tags"]][["custom_tags"]],
                                 pattern = ",")[[1]]
      for (t in tags) {
        t <- stringr::str_trim(t)
        if (t %in% colnames(cpg_calls)) {
          predispose_tsv_tags <- c(predispose_tsv_tags, t)
        }else{
          rlogging::wwarning(
            paste0("Could NOT detect the following tag in query VCF: ", tag))
        }
      }
    }
  }


  rlogging::message("Generating tiered set of result variants for ",
                    "output in tab-separated values (TSV) file")

  snv_indel_report <- pcgrr::init_report(config = cpsr_config,
                                         class = "snv_indel",
                                         type = "predisposition")

  predispose_tags <- predispose_tsv_tags

  snv_indel_report[["variant_set"]][["class5"]] <- cpg_calls %>%
    dplyr::filter(!is.na(CLINVAR_CLASSIFICATION) &
                    CLINVAR_CLASSIFICATION == "Pathogenic") %>%
    dplyr::mutate(SOURCE = "ClinVar")

  snv_indel_report[["variant_set"]][["class4"]] <- cpg_calls %>%
    dplyr::filter(!is.na(CLINVAR_CLASSIFICATION) &
                    CLINVAR_CLASSIFICATION == "Likely_Pathogenic") %>%
    dplyr::mutate(SOURCE = "ClinVar")

  snv_indel_report[["variant_set"]][["class3"]] <- cpg_calls %>%
    dplyr::filter(!is.na(CLINVAR_CLASSIFICATION) &
                    CLINVAR_CLASSIFICATION == "VUS") %>%
    dplyr::mutate(SOURCE = "ClinVar")

  snv_indel_report[["variant_set"]][["class2"]] <- cpg_calls %>%
    dplyr::filter(!is.na(CLINVAR_CLASSIFICATION) &
                    CLINVAR_CLASSIFICATION == "Likely_Benign") %>%
    dplyr::mutate(SOURCE = "ClinVar")

  snv_indel_report[["variant_set"]][["class1"]] <- cpg_calls %>%
    dplyr::filter(!is.na(CLINVAR_CLASSIFICATION) &
                    CLINVAR_CLASSIFICATION == "Benign") %>%
    dplyr::mutate(SOURCE = "ClinVar")

  ## identify remaining calls not registered in clinvar
  all_clinvar_calls <- data.frame()
  for (c in c("class1", "class2", "class3", "class4", "class5")) {
    all_clinvar_calls <- all_clinvar_calls %>%
      dplyr::bind_rows(dplyr::select(snv_indel_report[["variant_set"]][[c]],
                                     VAR_ID))
  }
  cpg_calls <- cpg_calls %>%
    dplyr::anti_join(all_clinvar_calls, by = c("VAR_ID"))

  n_nonclinvar <- nrow(cpg_calls)

  cpg_calls <- cpg_calls %>%
    dplyr::filter(is.na(GLOBAL_AF_GNOMAD) |
                    GLOBAL_AF_GNOMAD < cpsr_config[["maf_upper_threshold"]])
  n_maf_filtered <- n_nonclinvar - nrow(cpg_calls)
  rlogging::message(
    paste0("Ignoring n = ", n_maf_filtered,
           " unclassified variants with a global MAF frequency above ",
           cpsr_config[["maf_upper_threshold"]]))

  n_after_maf_filtering <- nrow(cpg_calls)

  non_clinvar_calls <- list()
  non_clinvar_calls[["class5"]] <- cpg_calls %>%
    dplyr::filter(CPSR_CLASSIFICATION == "Pathogenic") %>%
    dplyr::mutate(SOURCE = "Other")

  non_clinvar_calls[["class4"]] <- cpg_calls %>%
    dplyr::filter(CPSR_CLASSIFICATION == "Likely_Pathogenic") %>%
    dplyr::mutate(SOURCE = "Other")

  non_clinvar_calls[["class3"]] <- cpg_calls %>%
    dplyr::filter(CPSR_CLASSIFICATION == "VUS") %>%
    dplyr::mutate(SOURCE = "Other")

  non_clinvar_calls[["class2"]] <- cpg_calls %>%
    dplyr::filter(CPSR_CLASSIFICATION == "Likely_Benign") %>%
    dplyr::mutate(SOURCE = "Other")

  non_clinvar_calls[["class1"]] <- cpg_calls %>%
    dplyr::filter(CPSR_CLASSIFICATION == "Benign") %>%
    dplyr::mutate(SOURCE = "Other")

  for (c in c("class1", "class2", "class3", "class4", "class5")) {
    snv_indel_report[["variant_set"]][[c]] <-
      dplyr::bind_rows(non_clinvar_calls[[c]],
                       snv_indel_report[["variant_set"]][[c]])

    ## If not 'classify_all' is turned on,
    ## remove CPSR classifications for existing ClinVar classifications
    if (cpsr_config[["classify_all"]] == F) {
      snv_indel_report[["variant_set"]][[c]] <-
        snv_indel_report[["variant_set"]][[c]] %>%
        dplyr::mutate(
          CPSR_CLASSIFICATION =
            dplyr::if_else(!is.na(CLINVAR_CLASSIFICATION),
                           "", as.character(CPSR_CLASSIFICATION))) %>%
        dplyr::mutate(
          CPSR_CLASSIFICATION_DOC =
            dplyr::if_else(!is.na(CLINVAR_CLASSIFICATION),
                           "", as.character(CPSR_CLASSIFICATION_DOC))) %>%
        dplyr::mutate(
          CPSR_CLASSIFICATION_CODE =
            dplyr::if_else(!is.na(CLINVAR_CLASSIFICATION),
                           "", as.character(CPSR_CLASSIFICATION_CODE))) %>%
        dplyr::mutate(
          CPSR_PATHOGENICITY_SCORE =
            dplyr::if_else(!is.na(CLINVAR_CLASSIFICATION),
                           as.numeric(NA),
                           as.numeric(CPSR_PATHOGENICITY_SCORE)))
    }
    snv_indel_report[["variant_set"]][[c]] <-
      snv_indel_report[["variant_set"]][[c]] %>%
      dplyr::arrange(SOURCE, desc(CANCER_PHENOTYPE),
                     desc(CPSR_PATHOGENICITY_SCORE))
      #dplyr::select(-CANCER_PHENOTYPE)

    snv_indel_report[["disp"]][[c]] <-
      dplyr::select(snv_indel_report[["variant_set"]][[c]],
                    dplyr::one_of(display_tags))
    snv_indel_report[["variant_set"]][[c]] <-
      dplyr::select(snv_indel_report[["variant_set"]][[c]],
                    dplyr::one_of(predispose_tags))

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
      dplyr::select(c("GENOMIC_CHANGE", "VAR_ID", "GENOTYPE", "SOURCE",
                      "GENOME_VERSION",
                      "VCF_SAMPLE_ID", "VARIANT_CLASS", "CODING_STATUS",
                      "SYMBOL", "GENE_NAME", "CCDS", "ENTREZ_ID", "UNIPROT_ID",
                      "ENSEMBL_GENE_ID", "ENSEMBL_TRANSCRIPT_ID", "REFSEQ_MRNA",
                      "ONCOGENE", "TUMOR_SUPPRESSOR", "CONSEQUENCE",
                      "VEP_ALL_CSQ",
                      "PROTEIN_CHANGE", "PROTEIN_DOMAIN", "DBSNP", "HGVSp",
                      "HGVSc", "LAST_EXON",
                      "EXON_POSITION","INTRON_POSITION","CDS_CHANGE",
                      "MUTATION_HOTSPOT", "RMSK_HIT",
                      "PROTEIN_FEATURE", "EFFECT_PREDICTIONS",
                      "LOSS_OF_FUNCTION", "DBSNP","CANCER_PHENOTYPE",
                      "CLINVAR_CLASSIFICATION", "CLINVAR_MSID",
                      "CLINVAR_VARIANT_ORIGIN",
                      "CLINVAR_CONFLICTED", "CLINVAR_PHENOTYPE",
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
      c("GLOBAL_AF_GNOMAD", cpsr_config[["popgen"]][["vcftag_gnomad"]]))
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
retrieve_secondary_calls <- function(calls, umls_map) {

  secondary_calls <- calls %>%
      dplyr::filter(!is.na(CANCER_PREDISPOSITION_SOURCE) &
                      CANCER_PREDISPOSITION_SOURCE == "ACMG_SF20" &
                      !is.na(CLINVAR_CLASSIFICATION)) %>%
      dplyr::filter(CLINVAR_CLASSIFICATION == "Pathogenic" |
                      CLINVAR_CLASSIFICATION == "Likely_Pathogenic")

  if (nrow(secondary_calls) > 0) {
    secondary_calls_per_trait <-
      tidyr::separate_rows(secondary_calls, CLINVAR_UMLS_CUI, sep = ",") %>%
      dplyr::select(VAR_ID, CLINVAR_UMLS_CUI) %>%
      dplyr::left_join(umls_map, by = c("CLINVAR_UMLS_CUI" = "cui")) %>%
      dplyr::rename(CLINVAR_PHENOTYPE = cui_name) %>%
      dplyr::group_by(VAR_ID) %>%
      dplyr::summarise(CLINVAR_PHENOTYPE =
                         paste(unique(CLINVAR_PHENOTYPE), collapse="; ")) %>%
      dplyr::distinct()

    secondary_calls <-
      dplyr::inner_join(secondary_calls,
                        dplyr::select(secondary_calls_per_trait,
                                      VAR_ID, CLINVAR_PHENOTYPE),
                        by = c("VAR_ID" = "VAR_ID"))
  }

  return(secondary_calls)

}

#' Function that retrieves variants in cancer predisposition genes linked
#' to cancer-related conditions according to ClinVar
#'
#' @param calls data frame with variants in predisposition_genes
#' @param oncotree data frame with hereditary cancer phenotypes from OncoTree
#' @param umls_map
#'
detect_cancer_traits_clinvar <- function(cpg_calls, oncotree, umls_map) {
  if (nrow(cpg_calls) > 0 & "CLINVAR_UMLS_CUI" %in% colnames(cpg_calls)
      & "VAR_ID" %in% colnames(cpg_calls)) {

    oncotree <- oncotree %>%
      dplyr::select(cui) %>%
      dplyr::mutate(CANCER_PHENOTYPE = 1) %>%
      dplyr::distinct()

    n_clinvar <- cpg_calls %>% dplyr::filter(!is.na(CLINVAR_UMLS_CUI)) %>%
      nrow()

    if (n_clinvar > 0) {
      cpg_calls_traits <- as.data.frame(
        tidyr::separate_rows(cpg_calls, CLINVAR_UMLS_CUI, sep = ",") %>%
        dplyr::select(VAR_ID, CLINVAR_UMLS_CUI) %>%
        dplyr::left_join(umls_map, by = c("CLINVAR_UMLS_CUI" = "cui")) %>%
        dplyr::distinct() %>%
        dplyr::filter(!is.na(cui_name)) %>%
        dplyr::left_join(oncotree, by = c("CLINVAR_UMLS_CUI" = "cui")) %>%
        dplyr::mutate(
          CANCER_PHENOTYPE = dplyr::if_else(is.na(CANCER_PHENOTYPE),
                                            as.integer(0),
                                            as.integer(CANCER_PHENOTYPE))) %>%
        dplyr::mutate(
          CANCER_PHENOTYPE =
            dplyr::if_else(stringr::str_detect(tolower(cui_name),
                                               pcgrr::cancer_phenotypes_regex),
                           as.integer(1),
                           as.integer(CANCER_PHENOTYPE))) %>%
        dplyr::group_by(VAR_ID) %>%
        dplyr::summarise(
          CLINVAR_PHENOTYPE = paste(unique(cui_name), collapse = "; "),
          CANCER_PHENOTYPE = max(CANCER_PHENOTYPE),
          .groups = "drop") %>%
          dplyr::mutate(
            CANCER_PHENOTYPE =
              dplyr::if_else(
                stringr::str_detect(
                  CLINVAR_PHENOTYPE,
                  "^(not specified; not provided|not specified|not provided)"),
                as.integer(1),
                as.integer(CANCER_PHENOTYPE)))

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
summary_donut_chart <- function(variants_tsv, plot_type = "ClinVar") {

  title <- "ClinVar variants"
  p <- NULL

  if (nrow(variants_tsv) > 0) {

    set_clinvar <- variants_tsv %>%
      dplyr::filter(!is.na(CLINVAR_CLASSIFICATION) &
                      !(CLINVAR_CLASSIFICATION == "NA"))
    set_other <- variants_tsv %>%
      dplyr::filter(nchar(CPSR_CLASSIFICATION) > 0 &
                      (is.na(CLINVAR_CLASSIFICATION) |
                         CLINVAR_CLASSIFICATION == "NA"))

    if ((plot_type == "ClinVar" & nrow(set_clinvar) > 0) |
        (plot_type != "ClinVar" & nrow(set_other) > 0)) {

      m <- data.frame()

      if (plot_type == "ClinVar") {
        if (nrow(set_clinvar) > 0) {
          t <- paste0("n = ", nrow(set_clinvar))
          title <- bquote("ClinVar variants, "~bold(.(t)))
          m <- as.data.frame(set_clinvar %>%
            dplyr::group_by(CLINVAR_CLASSIFICATION) %>%
            dplyr::summarise(n = dplyr::n()) %>%
            dplyr::rename(level = CLINVAR_CLASSIFICATION)) %>%
            dplyr::mutate(
              level =
                factor(
                  level, levels =
                    pcgrr::color_palette[["pathogenicity"]][["levels"]])) %>%
            dplyr::arrange(level) %>%
            dplyr::mutate(prop = as.numeric(n / sum(n))) %>%
            dplyr::mutate(lab.ypos = cumsum(prop) - 0.5 * prop) %>%
            dplyr::mutate(n = as.character(n))
        }
      }else{
        if (nrow(set_other) > 0) {
          t <- paste0("n = ", nrow(set_other))
          title <- bquote("Other variants, CPSR-classified, "~bold(.(t)))
          m <- as.data.frame(set_other %>%
            dplyr::group_by(CPSR_CLASSIFICATION) %>%
            dplyr::summarise(n = dplyr::n()) %>%
            dplyr::rename(level = CPSR_CLASSIFICATION)) %>%
            dplyr::mutate(
              level = factor(
                level,
                levels = pcgrr::color_palette[["pathogenicity"]][["levels"]])) %>%
            dplyr::arrange(level) %>%
            dplyr::mutate(prop = as.numeric(n / sum(n))) %>%
            dplyr::mutate(lab.ypos = cumsum(prop) - 0.5 * prop) %>%
            dplyr::mutate(n = as.character(n))

        }
      }


      p <- ggplot2::ggplot(m, ggplot2::aes(x = 2, y = prop, fill = level)) +
        ggplot2::geom_bar(stat = "identity", color = "white") +
        ggplot2::coord_polar(theta = "y", start = 0) +
        ggplot2::geom_text(ggplot2::aes(y = 1 - lab.ypos, label = n),
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


virtual_panel_display_html <- function(gene_df) {

  i <- 1
  gene_df <- gene_df %>%
    dplyr::arrange(desc(confidence_level), symbol)

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


#' Function that retrieves clinical evidence items (CIVIC, CBMDB) for
#'somatic cancer variants
#'
#' @param sample_calls data frame with germline variant callset from
#' query sample
#' @param colset vector with column names to display for each report element
#' @param eitems all clinical evidence items linking germline variants with
#' impact on therapeutic response/diagnosis/prognosis etc
#'
#' @return list

get_germline_biomarkers <- function(sample_calls,
                                    colset = NULL,
                                    eitems = NULL) {

  rlogging::message("Matching variant set with existing genomic",
                    " biomarkers from CIViC (germline)")

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
    dplyr::filter(!is.na(CPSR_CLASSIFICATION) |
                  !is.na(CLINVAR_CLASSIFICATION)) %>%
    dplyr::filter(stringr::str_detect(CPSR_CLASSIFICATION, "Pathogenic") |
                  stringr::str_detect(CLINVAR_CLASSIFICATION, "Pathogenic"))


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
        dplyr::filter(EVIDENCE_TYPE == stringr::str_to_title(type)) %>%
        dplyr::arrange(EVIDENCE_LEVEL, RATING) %>%
        dplyr::select(SYMBOL, GENE_NAME, CANCER_TYPE, CLINICAL_SIGNIFICANCE,
                      EVIDENCE_LEVEL, RATING, EVIDENCE_DIRECTION, CITATION,
                      THERAPEUTIC_CONTEXT, DESCRIPTION, BIOMARKER_MAPPING,
                      CDS_CHANGE, LOSS_OF_FUNCTION, GENOMIC_CHANGE) %>%
        dplyr::distinct()
    }
  }

  return(clin_eitems_list)

}

#' Function that reports protein-coding geneset that overlaps BED file
#'
#' @param bed_file BED file name with selected transcripts from panel 0
#' @param pcgr_data object with PCGR annotation data
#'
custom_bed_genes <- function(bed_file, pcgr_data) {

  invisible(assertthat::assert_that(file.exists(bed_file),
                                    msg = paste0("BED file", bed_file, " does not exist")))
  bed_df <- as.data.frame(
    read.table(file = bed_file, header = F, stringsAsFactors = F,
                           comment.char = "", quote = "", sep = "\t") %>%
    magrittr::set_colnames(c("chromosome", "segment_start", "segment_end", "onco_xref")) %>%
    dplyr::rowwise() %>%
    dplyr::mutate(symbol = unlist(strsplit(onco_xref, "\\|"))[4]) %>%
    dplyr::select(symbol) %>%
    dplyr::left_join(dplyr::select(pcgr_data$gene_xref$gencode, symbol, ENTREZ_ID, GENENAME), by = "symbol") %>%
    dplyr::rename(genename = GENENAME, entrezgene = ENTREZ_ID) %>%
    dplyr::select(symbol, genename, entrezgene) %>%
    dplyr::distinct()
  )
  return(bed_df)
}

