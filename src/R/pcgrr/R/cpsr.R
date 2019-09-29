
#' Function that generates predisposition_report - CPSR
#'
#' @param project_directory name of project directory
#' @param query_vcf2tsv name of gzipped TSV file (vcf2tsv) with annotated query SNVs/InDels
#' @param pcgr_data List of data frames with PCGR data annotations
#' @param cpsr_config Object with CPSR configuration parameters
#' @param virtual_panel_id Identifier for virtual panel
#' @param sample_name sample identifier
#'

generate_predisposition_report <- function(project_directory, query_vcf2tsv, pcgr_data, cpsr_config = NULL, virtual_panel_id = -1, diagnostic_grade_only = 0, sample_name = "SampleX"){

  cps_report <- pcgrr::init_pcg_report(cpsr_config, sample_name, class = NULL, pcgr_data = pcgr_data, type = "predisposition", virtual_panel_id = virtual_panel_id, diagnostic_grade_only = diagnostic_grade_only)

  fnames <- list()
  fnames[["tsv"]] <- paste0(project_directory, "/", sample_name, ".cpsr.snvs_indels.tiers.",pcgr_data[['assembly']][['grch_name']],".tsv")

  ## define tags/variables to display in data tables (class 1 - 5)
  class_1_5_display <- c("SYMBOL", "SOURCE", "CLINVAR_PHENOTYPE", "CONSEQUENCE", "PROTEIN_CHANGE", "GENOTYPE", "GENE_NAME", "PROTEIN_DOMAIN",
                          "HGVSp", "HGVSc", "CDS_CHANGE", "REFSEQ_MRNA","MUTATION_HOTSPOT", "RMSK_HIT","PROTEIN_FEATURE", "PREDICTED_EFFECT",
                          "LOSS_OF_FUNCTION", "DBSNP", "CLINVAR","CLINVAR_CLASSIFICATION","CLINVAR_REVIEW_STATUS_STARS","CLINVAR_CONFLICTED",
                          "CLINVAR_VARIANT_ORIGIN", "CPSR_CLASSIFICATION","CPSR_CLASSIFICATION_SCORE","CPSR_CLASSIFICATION_DOC",
                          "CPSR_CLASSIFICATION_CODE","ONCOGENE", "TUMOR_SUPPRESSOR",
                          "GLOBAL_AF_GNOMAD", cps_report[["metadata"]][["config"]][["popgen"]][["vcftag_gnomad"]],
                          "GENOMIC_CHANGE", "GENOME_VERSION")

  ## define tags/variables to display in data tables (secondary findings)
  secondary_findings_display <- c("SYMBOL", "CLIN_SIGNIFICANCE", "CLINVAR_PHENOTYPE", "CONSEQUENCE", "PROTEIN_CHANGE", "GENOTYPE", "GENE_NAME", "PROTEIN_DOMAIN",
                               "HGVSp", "HGVSc", "CDS_CHANGE", "REFSEQ_MRNA","MUTATION_HOTSPOT","RMSK_HIT","PROTEIN_FEATURE", "PREDICTED_EFFECT",
                               "LOSS_OF_FUNCTION", "DBSNP", "CLINVAR", "CLINVAR_REVIEW_STATUS_STARS","CLINVAR_CONFLICTED",
                               "CLINVAR_CLINICAL_SIGNIFICANCE", "CLINVAR_VARIANT_ORIGIN",
                               "GWAS_CITATION","ONCOGENE", "TUMOR_SUPPRESSOR",
                               "GLOBAL_AF_GNOMAD", cps_report[["metadata"]][["config"]][["popgen"]][["vcftag_gnomad"]],
                               "GENOMIC_CHANGE", "GENOME_VERSION")

  predispose_gwas_display <- c("SYMBOL","CONSEQUENCE", "GWAS_CITATION","PROTEIN_CHANGE","GENOTYPE","LOSS_OF_FUNCTION",
                                 "PROTEIN_CHANGE","GENE_NAME", "GWAS_PHENOTYPE","PROTEIN_DOMAIN","HGVSp", "HGVSc", "CDS_CHANGE","CODING_STATUS", "REFSEQ_MRNA",
                               "PROTEIN_FEATURE", "PREDICTED_EFFECT",
                                 "DBSNP","GLOBAL_AF_GNOMAD",
                               cps_report[["metadata"]][["config"]][["popgen"]][["vcftag_gnomad"]],
                                 "GENOMIC_CHANGE","GENOME_VERSION")


  phenotype_medgen_cancer <- pcgr_data[['phenotype_ontology']][['medgen_cancer']] %>%
    dplyr::filter(group == "Hereditary_Cancer_Syndrome_NOS" | group == "Hereditary_Cancer_Susceptibility_NOS") %>%
    dplyr::filter(!is.na(cui_name)) %>%
    dplyr::select(cui, cui_name) %>%
    dplyr::mutate(cancer_phenotype = 1) %>%
    dplyr::distinct()

  medgen <- pcgr_data[['phenotype_ontology']][['medgen_all']] %>%
    dplyr::left_join(phenotype_medgen_cancer, by = c("cui", "cui_name"))


  if (query_vcf2tsv != "None.gz"){
    if (!file.exists(query_vcf2tsv) | file.size(query_vcf2tsv) == 0){
      rlogging::warning(paste0("File ", query_vcf2tsv, " does not exist or has zero size"))
    }else{
      if (!is.null(cpsr_config) & query_vcf2tsv != "None.gz"){

        ## read calls
        calls <- pcgrr::get_calls(query_vcf2tsv, pcgr_data, sample_name, cpsr_config, medgen_ont = cps_report[['metadata']][['medgen_ontology']][['all']], cpsr = TRUE)
        calls <- dplyr::rename(calls, CLINVAR_CLINICAL_SIGNIFICANCE = CLINVAR_CLNSIG)
        call_stats <- pcgrr::variant_stats_report(calls, name = "variant_statistic")

        cpg_calls <- dplyr::inner_join(calls, cps_report[['metadata']][['gene_panel']][['genes']], by = c("SYMBOL" = "symbol"))
        cpg_call_stats <- pcgrr::variant_stats_report(cpg_calls, name = "variant_statistic_cpg")

        rlogging::message(paste0("Number of coding variants in cancer predisposition genes: ",cpg_call_stats[['variant_statistic_cpg']][['n_coding']]))
        rlogging::message(paste0("Number of non-coding variants in cancer predisposition genes: ",cpg_call_stats[['variant_statistic_cpg']][['n_noncoding']]))
        sf_calls <- pcgrr::retrieve_sf_calls(calls, medgen)
        sf_call_stats <- pcgrr::variant_stats_report(sf_calls, name = "variant_statistic_sf")

        #gene_hits <- paste(unique(cpg_calls$SYMBOL), collapse = ", ")
        #if(nrow(cpg_calls) > 0){
          #rlogging::message("Variants were found in the following cancer predisposition genes: ", gene_hits)
        if(nrow(cpg_calls) == 0){
          cps_report$variant_statistic_cpg <- cpg_call_stats$variant_statistic_cpg
          return(cps_report)
        }

        cpg_calls <- pcgrr::assign_pathogenicity_evidence(cpg_calls, cpsr_config, pcgr_data) %>%
          pcgrr::determine_pathogenicity_classification() %>%
          pcgrr::detect_cancer_traits_clinvar(medgen)
        snv_indel_report <- pcgrr::assign_cpsr_tier(cpg_calls, cps_report[['metadata']][['config']], class_1_5_display)
        snv_indel_report$variant_statistic <- call_stats$variant_statistic
        snv_indel_report$variant_statistic_cpg <- cpg_call_stats$variant_statistic_cpg
        snv_indel_report$variant_statistic_sf <- sf_call_stats$variant_statistic_sf
        snv_indel_report <- pcgrr::get_germline_biomarkers(snv_indel_report, pcgr_data)



        cps_report <- pcgrr::update_pcg_report(cps_report, report_data = snv_indel_report)


        gene_hits <- paste(unique(cps_report[["content"]][["snv_indel"]][["variant_set"]][['tsv']]$SYMBOL),collapse=", ")
        rlogging::message("Variants were found in the following cancer predisposition genes: ", gene_hits)

        if(cpsr_config[['secondary_findings']][['show_sf']] == TRUE){
          rlogging::message("Assignment of other variants in genes recommended for reporting as incidental findings (ACMG SF v2.0)")
          cps_report[['content']][["snv_indel"]][["variant_display"]][['sf']] <- sf_calls %>%
            dplyr::arrange(LOSS_OF_FUNCTION, CODING_STATUS) %>%
            dplyr::select(dplyr::one_of(secondary_findings_display))
          rlogging::message(paste0("Number of pathogenic variants in the incidentalome - other genes of clinical significance: ",cps_report[['content']][['snv_indel']][['variant_statistic_sf']][['n_coding']]))
        }

        cps_report[['content']][['snv_indel']][['eval']] <- TRUE

        if(cpsr_config[['gwas']][['gwas_hits']] == TRUE){
          rlogging::message("Assignment of other variants to hits from genome-wide association studies")
        }
        ## Assign GWAS hits to cps_report object
        cps_report[['content']][["snv_indel"]][["variant_display"]][["gwas"]] <-
          dplyr::filter(calls, !is.na(GWAS_HIT) & !is.na(GWAS_CITATION))
        if(nrow(cps_report[['content']][["snv_indel"]][["variant_display"]][["gwas"]]) > 0){
          if(nrow(cps_report[['content']][['snv_indel']][['variant_set']][['tsv']]) > 0){
            cps_report[['content']][["snv_indel"]][["variant_display"]][["gwas"]] <-
              dplyr::anti_join(cps_report[['content']][["snv_indel"]][["variant_display"]][["gwas"]],
                               cps_report[['content']][['snv_indel']][['variant_set']][['tsv']], by=c("GENOMIC_CHANGE"))
          }
          if(nrow(cps_report[['content']][["snv_indel"]][["variant_display"]][["gwas"]]) > 0){
            cps_report[['content']][["snv_indel"]][["variant_display"]][["gwas"]] <-
              cps_report[['content']][["snv_indel"]][["variant_display"]][["gwas"]] %>%
              dplyr::select(dplyr::one_of(predispose_gwas_display)) %>%
              dplyr::arrange(LOSS_OF_FUNCTION, CODING_STATUS)
          }
        }
      }
    }

    fname_key <- "tsv"
    if (!is.null(cps_report[["content"]][["snv_indel"]][["variant_set"]][[fname_key]])){
      if (nrow(cps_report[["content"]][["snv_indel"]][["variant_set"]][[fname_key]]) > 0){
        write.table(cps_report[["content"]][["snv_indel"]][["variant_set"]][[fname_key]], file = fnames[[fname_key]], sep = "\t", col.names = T, row.names = F, na = "NA", quote = F)
      }
    }
    cps_report[['metadata']][['medgen_ontology']] <- list()
  }
  return(cps_report)
}

#' Function that counts insilico predictions of variant effects (i.e. damaging/tolerated) from dbNSFP
#'
#' @param cpg_calls sample calls with dbnsfp annotations
#'
#' @return cpg_calls
#'
get_insilico_prediction_statistics <- function(cpg_calls){

  insilico_pathogenicity_pred_algos <- c('SIFT_DBNSFP','PROVEAN_DBNSFP','META_LR_DBNSFP','FATHMM_DBNSFP','MUTATIONTASTER_DBNSFP','DEOGEN2_DBNSFP','PRIMATEAI_DBNSFP',
                                         'MUTATIONASSESSOR_DBNSFP','FATHMM_MKL_DBNSFP','M_CAP_DBNSFP','SPLICE_SITE_ADA_DBNSFP','SPLICE_SITE_RF_DBNSFP')
  for(v in c('CALLED','DAMAGING','TOLERATED','SPLICING_NEUTRAL','SPLICING_AFFECTED')){
    cpg_calls[,paste0('N_INSILICO_',v)] <- 0
  }

  for(algo in insilico_pathogenicity_pred_algos){
    if(algo %in% colnames(cpg_calls)){
      cpg_calls[!is.na(cpg_calls[,algo]) & cpg_calls[,algo] != '.','N_INSILICO_CALLED'] <-  cpg_calls[!is.na(cpg_calls[,algo]) & cpg_calls[,algo] != '.','N_INSILICO_CALLED'] + 1
      cpg_calls[!is.na(cpg_calls[,algo]) & (cpg_calls[,algo] == 'damaging' | cpg_calls[,algo] == 'possibly_damaging'),'N_INSILICO_DAMAGING'] <-  cpg_calls[!is.na(cpg_calls[,algo]) & (cpg_calls[,algo] == 'damaging' | cpg_calls[,algo] == 'possibly_damaging'),'N_INSILICO_DAMAGING'] + 1
      cpg_calls[!is.na(cpg_calls[,algo]) & cpg_calls[,algo] == 'tolerated','N_INSILICO_TOLERATED'] <-  cpg_calls[!is.na(cpg_calls[,algo]) & cpg_calls[,algo] == 'tolerated','N_INSILICO_TOLERATED'] + 1
      cpg_calls[!is.na(cpg_calls[,algo]) & cpg_calls[,algo] == 'affect_splicing','N_INSILICO_SPLICING_AFFECTED'] <-  cpg_calls[!is.na(cpg_calls[,algo]) & cpg_calls[,algo] == 'affect_splicing','N_INSILICO_SPLICING_AFFECTED'] + 1
      cpg_calls[!is.na(cpg_calls[,algo]) & cpg_calls[,algo] == 'splicing_neutral','N_INSILICO_SPLICING_NEUTRAL'] <-  cpg_calls[!is.na(cpg_calls[,algo]) & cpg_calls[,algo] == 'splicing_neutral','N_INSILICO_SPLICING_NEUTRAL'] + 1
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
assign_pathogenicity_evidence <- function(cpg_calls, cpsr_config, pcgr_data){

  gad_population <- toupper(cpsr_config[['popgen']][['pop_gnomad']])
  gad_AN_tag <- paste0('NON_CANCER_AN_',gad_population)
  gad_AF_tag <- paste0('NON_CANCER_AF_',gad_population)
  gad_NHOMALT_tag <- paste0('NON_CANCER_NHOMALT_',gad_population)
  gad_AC_tag <- paste0('NON_CANCER_AC_',gad_population)

  acmg_ev_codes <- c('ACMG_BA1_AD',   ## Very high MAF (> 0.5% in gnomAD non-cancer pop subset) - min AN = 12,000, min AC = 12 - Dominant mechanism of disease
                     'ACMG_BS1_1_AD', ## High MAF (> 0.1% in gnomAD non-cancer pop subset) - min AN = 12,000, min AC = 12 - Dominant mechanism of disease
                     'ACMG_BS1_2_AD', ## Somewhat high AF (> 8 alleles in gnomAD non-cancer pop subset) - Dominant mechanism of disease
                     'ACMG_BA1_AR',   ## Very high MAF (> 1% in gnomAD non-cancer pop subset) - min AN = 12,000, min AC = 12 - Recessive mechanism of disease
                     'ACMG_BS1_1_AR', ## High MAF (> 0.3% in gnomAD non-cancer pop subset) - min AN = 12,000, min AC = 12 - Recessive mechanism of disease
                     'ACMG_BS1_2_AR', ## Somewhat high AF (> 8 alleles in gnomAD non-cancer pop subset) - Recessive mechanism of disease
                     #'ACMG_BS2_1',    ## 1 homozygote in gnomAD non-cancer pop subset - severe, early onset, highly penetrant
                     #'ACMG_BS2_2',    ## 2 homozygotes in gnomAD non-cancer pop subset - severe, early onset, highly penetrant
                     #'ACMG_BS2_3',    ## 2 homozygotes in gnomAD non-cancer pop subset - moderate, early onset, variably penetrant
                     'ACMG_PM2_1',    ## Allele count within pathogenic range (8 or fewer alleles in the population-specific non-cancer gnomAD subset)
                     'ACMG_PM2_2',    ## Alternate allele absent in the population-specific non-cancer gnomAD subset
                     'ACMG_PVS1_1',   ## Null variant - predicted as LoF by LOFTEE - within pathogenic range - LoF established for gene
                     'ACMG_PVS1_2',   ## Null variant - not predicted as LoF by LOFTEE - within pathogenic range - LoF established for gene
                     'ACMG_PVS1_3',   ## Null variant - predicted as LoF by LOFTEE - within pathogenic range - LoF not established for gene
                     'ACMG_PVS1_4',   ## Null variant - not predicted as LoF by LOFTEE -- within pathogenic range - LoF not established for gene
                     'ACMG_PVS1_5',   ## start lost - within pathogenic range - Lof established for gene
                     'ACMG_PVS1_6',   ## start lost - within pathogenic range - LoF not established for gene
                     'ACMG_PVS1_7',   ## donor/acceptor variant - predicted as LoF by LOFTEE - within pathogenic range
                                      ## - not last intron - LoF established for gene
                     'ACMG_PVS1_8',   ## donor/acceptor variant - last intron - within pathogenic range - LoF established for gene
                     'ACMG_PVS1_9',   ## donor/acceptor variant - not last intron - within pathogenic range - LoF not established for gene
                     'ACMG_PVS1_10',  ## donor variant at located at the +3, +4 or +5 position of the intron -  within the pathogenic range (i.e. <9 alleles in ExAC))
                     'ACMG_PS1',      ## Same amino acid change as a previously established pathogenic variant (ClinVar) regardless of nucleotide change
                     'ACMG_PP2',      ## Missense variant in a gene that has a relatively low rate of benign missense variation (<20%) and
                                      ## where missense variants are a common mechanism of disease (>50% of high-confidence pathogenic variants (ClinVar))
                     'ACMG_PM4',      ## Protein length changes due to inframe indels or nonstop variant in non-repetitive regions of genes
                                      ## that harbor variants with a dominant mode of inheritance.
                     'ACMG_PPC1',     ## Protein length changes due to inframe indels or nonstop variant in non-repetitive regions of genes
                                      ## that harbor variants with a recessive mode of inheritance.
                     'ACMG_PM5',      ## Novel missense change at an amino acid residue where a different missense change determined to be pathogenic
                                      ## has been seen before (ClinVar)
                     'ACMG_PP3',      ## Multiple lines of computational evidence support a deleterious effect on the gene or gene product
                                      ## (conservation, evolutionary, splicing impact, etc. - from dbNSFP
                     'ACMG_BP4',      ## Multiple lines of computational evidence support a benign effect on the gene or gene product
                                      ## (conservation, evolutionary, splicing impact, etc. - from dbNSFP
                     'ACMG_BMC1',     ## Peptide change is at the same location of a known benign change (ClinVar)
                     'ACMG_BSC1',     ## Peptide change is reported as benign (ClinVar),
                     'ACMG_BP1')      ## Missense variant in a gene for which primarily truncating variants are known to cause disease (ClinVar)


    path_columns <- c(acmg_ev_codes,"N_INSILICO_CALLED",
                    "N_INSILICO_DAMAGING", "N_INSILICO_TOLERATED", "N_INSILICO_SPLICING_NEUTRAL",
                    "N_INSILICO_SPLICING_AFFECTED", "codon_prefix", "clinvar_pathogenic_codon",
                    "clinvar_pathogenic", "clinvar_benign","clinvar_benign_codon","hotspot_symbol", "hotspot_codon","hotspot_pvalue",
                    "MOD","PATH_TRUNCATION_RATE","BENIGN_MISSENSE_RATE")
  cpg_calls <- cpg_calls[, !(colnames(cpg_calls) %in% path_columns)]


  cpg_calls <- cpg_calls %>%
    dplyr::mutate(cpsr_gene_moi = dplyr::if_else(stringr::str_detect(CANCER_PREDISPOSITION_MOI,"AD|AD/AR"),"AD",as.character(NA),as.character(NA))) %>%
    dplyr::mutate(cpsr_gene_moi = dplyr::if_else(!stringr::str_detect(CANCER_PREDISPOSITION_MOI,"AD|Mosaic"),"AR",cpsr_gene_moi, cpsr_gene_moi))

  predisposition_gene_info <- dplyr::select(pcgr_data[['predisposition']][['genes']], symbol, mechanism_of_disease, path_truncation_rate, benign_missense_rate) %>%
    dplyr::rename(SYMBOL = symbol,
                  MOD = mechanism_of_disease,
                  PATH_TRUNCATION_RATE = path_truncation_rate,
                  BENIGN_MISSENSE_RATE = benign_missense_rate)

  cpg_calls <- dplyr::left_join(cpg_calls, predisposition_gene_info, by=c("SYMBOL" = "SYMBOL"))

  ## Assign logical ACMG evidence indicators
  # ACMG_PP3 - Multiple lines of computational evidence support a deleterious effect on the gene or gene product (conservation, evolutionary, splicing impact, etc.)
  # ACMG_BP4 - Multiple lines (>1) of in silico evidence of none deleterious effect.
  #
  # Computational evidence for deleterious/benign effect is taken from invidual algorithm predictions in dbNSFP
  # SIFT,Provean,MutationTaster,MutationAssessor,M_CAP,MutPred,FATHMM,FATHMM-mkl,DBNSFP_LogReg,dbscSNV_RF,dbscSNV_AdaBoost
  # Default scheme (from default TOML file):
  # 1) Damaging: Among 8 possible protein variant effect predictions, at least seven algorithms must have made a call, with at least 6 predicted as damaging,
  #    and at most one predicted as tolerated (PP3)
  #    - at most 1 prediction for a splicing neutral effect
  # 2) Tolerated: Among 8 possible protein variant effect predictions, at least seven algorithms must have made a call, with at least 6 predicted as tolerated,
  #    and at most one predicted as damaging (BP4)
  #    - 0 predictions of splice site affected

  dbnsfp_min_majority <- 6
  dbnsfp_max_minority <- 2
  dbnsfp_min_called <- dbnsfp_min_majority + 1

  cpg_calls <- pcgrr::get_insilico_prediction_statistics(cpg_calls)

  cpg_calls <- cpg_calls %>% dplyr::mutate(ACMG_PP3 = dplyr::if_else(N_INSILICO_CALLED >= dbnsfp_min_called &
                                                                       N_INSILICO_DAMAGING >= dbnsfp_min_majority &
                                                                       N_INSILICO_TOLERATED <= dbnsfp_max_minority &
                                                                       N_INSILICO_SPLICING_NEUTRAL <= 1,TRUE,FALSE,FALSE))
  cpg_calls <- cpg_calls %>% dplyr::mutate(ACMG_BP4 = dplyr::if_else(N_INSILICO_CALLED >= dbnsfp_min_called &
                                                                       N_INSILICO_TOLERATED >= dbnsfp_min_majority &
                                                                       N_INSILICO_DAMAGING <= dbnsfp_max_minority &
                                                                       N_INSILICO_SPLICING_AFFECTED == 0,TRUE,FALSE,FALSE))

  cpg_calls <- cpg_calls %>% dplyr::mutate(ACMG_PP3 = dplyr::case_when(N_INSILICO_SPLICING_AFFECTED == 2 ~ TRUE, TRUE ~ as.logical(ACMG_PP3)))

  ## Assign logical ACMG evidence indicators based on population frequency data in non-cancer samples from gnomAD (Dominant vs. recessive modes of inheritance)
  # 'ACMG_BA1_AD'   - Very high MAF (> 0.5% in gnomAD non-cancer pop subset) - min AN = 12,000, min AC = 12 - Dominant mechanism of disease
  # 'ACMG_BS1_1_AD' - High MAF (> 0.1% in gnomAD non-cancer pop subset) - min AN = 12,000, min AC = 12 - Dominant mechanism of disease
  # 'ACMG_BS1_2_AD' - Somewhat high AF (> 8 alleles in gnomAD non-cancer pop subset) - Dominant mechanism of disease
  # 'ACMG_BA1_AR'   - Very high MAF (> 1% in gnomAD non-cancer pop subset) - min AN = 12,000, min AC = 12 - Recessive mechanism of disease
  # 'ACMG_BS1_1_AR' - High MAF (> 0.3% in gnomAD non-cancer pop subset) - min AN = 12,000, min AC = 12 - Recessive mechanism of disease
  # 'ACMG_BS1_2_AR' - Somewhat high AF (> 8 alleles in gnomAD non-cancer pop subset) - Recessive mechanism of disease
  # 'ACMG_PM2_1'    - Allele count within pathogenic range (8 or fewer alleles in the population-specific non-cancer gnomAD subset)
  # 'ACMG_PM2_2'    - Alternate allele absent in the population-specific non-cancer gnomAD subset
  if(gad_AN_tag %in% colnames(cpg_calls) & gad_AC_tag %in% colnames(cpg_calls) & gad_NHOMALT_tag %in% colnames(cpg_calls)){
    cpg_calls <- cpg_calls %>%
      dplyr::mutate(ACMG_PM2_1 = dplyr::if_else(!is.na(!!rlang::sym(gad_AC_tag)) & !!rlang::sym(gad_AC_tag) <= 8,TRUE,FALSE,FALSE))
    cpg_calls <- cpg_calls %>%
      dplyr::mutate(ACMG_PM2_2 = dplyr::if_else(is.na(!!rlang::sym(gad_AC_tag)),TRUE,FALSE,FALSE))

    cpg_calls <- cpg_calls %>%
      dplyr::mutate(gad_an_ac_sufficient = dplyr::if_else(!!rlang::sym(gad_AN_tag) >= 12000 & !!rlang::sym(gad_AC_tag) >= 12,TRUE,FALSE,FALSE))
    cpg_calls <- cpg_calls %>%
      dplyr::mutate(gad_af = dplyr::if_else(gad_an_ac_sufficient == TRUE, as.numeric(!!rlang::sym(gad_AC_tag)/!!rlang::sym(gad_AN_tag)),as.double(NA),as.double(NA)))

    cpg_calls <- cpg_calls %>%
      dplyr::mutate(ACMG_BA1_AD = dplyr::if_else(ACMG_PM2_2 == FALSE & gad_af >= 0.005 & cpsr_gene_moi == "AD",TRUE,FALSE,FALSE))
    cpg_calls <- cpg_calls %>%
      dplyr::mutate(ACMG_BS1_1_AD = dplyr::if_else(ACMG_BA1_AD == FALSE & ACMG_PM2_2 == FALSE & gad_af >= 0.001 & cpsr_gene_moi == "AD",TRUE,FALSE,FALSE))
    cpg_calls <- cpg_calls %>%
      dplyr::mutate(ACMG_BS1_2_AD = dplyr::if_else(ACMG_BS1_1_AD == FALSE & ACMG_BA1_AD == FALSE & ACMG_PM2_2 == FALSE & !!rlang::sym(gad_AC_tag) > 8 & cpsr_gene_moi == "AD",TRUE,FALSE,FALSE))
    cpg_calls <- cpg_calls %>%
      dplyr::mutate(ACMG_BA1_AR = dplyr::if_else(ACMG_PM2_2 == FALSE & gad_af >= 0.01 & (cpsr_gene_moi == "AR" | is.na(cpsr_gene_moi)),TRUE,FALSE,FALSE))
    cpg_calls <- cpg_calls %>%
      dplyr::mutate(ACMG_BS1_1_AR = dplyr::if_else(ACMG_BA1_AR == FALSE & ACMG_PM2_2 == FALSE & gad_af >= 0.003 & (cpsr_gene_moi == "AR" | is.na(cpsr_gene_moi)),TRUE,FALSE,FALSE))
    cpg_calls <- cpg_calls %>%
      dplyr::mutate(ACMG_BS1_2_AR = dplyr::if_else(ACMG_BA1_AR == FALSE & ACMG_BS1_1_AR == FALSE & ACMG_PM2_2 == FALSE & !!rlang::sym(gad_AC_tag) > 8 & (cpsr_gene_moi == "AR" | is.na(cpsr_gene_moi)),TRUE,FALSE,FALSE))
  }

  ## Assign logical ACMG evidence indicators on NULL variants in known predisposition genes (LoF established as mechanism of disease or not, presumed loss of mRNA/protein (LOFTEE) or not)
  # 'ACMG_PVS1_1' - Null variant (frameshift, nonsense) - predicted as LoF by LOFTEE - within pathogenic range - LoF established for gene
  # 'ACMG_PVS1_2' - Null variant (frameshift, nonsense) - not predicted as LoF by LOFTEE - within pathogenic range - LoF established for gene
  # 'ACMG_PVS1_3' - Null variant (frameshift, nonsense) - predicted as LoF by LOFTEE - within pathogenic range - LoF not established for gene
  # 'ACMG_PVS1_4' - Null variant (frameshift, nonsense) - not predicted as LoF by LOFTEE -- within pathogenic range - LoF not established for gene
  # 'ACMG_PVS1_5' - start lost - within pathogenic range - Lof established for gene
  # 'ACMG_PVS1_6' - start lost - within pathogenic range - LoF not established for gene
  cpg_calls <- cpg_calls %>%
    dplyr::mutate(ACMG_PVS1_1 = dplyr::if_else(NULL_VARIANT == T & LOSS_OF_FUNCTION == T & MOD == "LoF" & (ACMG_PM2_1 == TRUE | ACMG_PM2_2 == TRUE),TRUE,FALSE,FALSE))
  cpg_calls <- cpg_calls %>%
    dplyr::mutate(ACMG_PVS1_3 = dplyr::if_else(NULL_VARIANT == T & LOSS_OF_FUNCTION == T & (is.na(MOD) | MOD != "LoF") & (ACMG_PM2_1 == TRUE | ACMG_PM2_2 == TRUE),TRUE,FALSE,FALSE))
  cpg_calls <- cpg_calls %>%
    dplyr::mutate(ACMG_PVS1_2 = dplyr::if_else(NULL_VARIANT == T & LOSS_OF_FUNCTION == F & MOD == "LoF" & (ACMG_PM2_1 == TRUE | ACMG_PM2_2 == TRUE),TRUE,FALSE,FALSE))
  cpg_calls <- cpg_calls %>%
    dplyr::mutate(ACMG_PVS1_4 = dplyr::if_else(NULL_VARIANT == T & LOSS_OF_FUNCTION == F & (is.na(MOD) | MOD != "LoF") & (ACMG_PM2_1 == TRUE | ACMG_PM2_2 == TRUE),TRUE,FALSE,FALSE))
  cpg_calls <- cpg_calls %>%
    dplyr::mutate(ACMG_PVS1_5 = dplyr::if_else(CONSEQUENCE == "start_lost" & MOD == "LoF" & (ACMG_PM2_1 == TRUE | ACMG_PM2_2 == TRUE),TRUE,FALSE,FALSE))
  cpg_calls <- cpg_calls %>%
    dplyr::mutate(ACMG_PVS1_6 = dplyr::if_else(CONSEQUENCE == "start_lost" & (is.na(MOD) | MOD != "LoF") & (ACMG_PM2_1 == TRUE | ACMG_PM2_2 == TRUE),TRUE,FALSE,FALSE))
  cpg_calls <- cpg_calls %>%
    dplyr::mutate(ACMG_PVS1_7 = dplyr::if_else(LOSS_OF_FUNCTION == T & stringr::str_detect(CONSEQUENCE,"_donor|_acceptor") &
                                                    LAST_INTRON == F & MOD == "LoF" & (ACMG_PM2_1 == TRUE | ACMG_PM2_2 == TRUE),TRUE,FALSE,FALSE))
  cpg_calls <- cpg_calls %>%
    dplyr::mutate(ACMG_PVS1_8 = dplyr::if_else(LOSS_OF_FUNCTION == T & stringr::str_detect(CONSEQUENCE,"_donor|_acceptor") &
                                                   LAST_INTRON == T & MOD == "LoF" & (ACMG_PM2_1 == TRUE | ACMG_PM2_2 == TRUE),TRUE,FALSE,FALSE))
  cpg_calls <- cpg_calls %>%
    dplyr::mutate(ACMG_PVS1_9 = dplyr::if_else(LOSS_OF_FUNCTION == T & stringr::str_detect(CONSEQUENCE,"_donor|_acceptor") &
                                                   LAST_INTRON == F & (is.na(MOD) | MOD != "LoF") & (ACMG_PM2_1 == TRUE | ACMG_PM2_2 == TRUE),TRUE,FALSE,FALSE))
  cpg_calls <- cpg_calls %>%
    dplyr::mutate(ACMG_PVS1_10 = dplyr::if_else(SPLICE_DONOR_RELEVANT == T & ACMG_PP3 == TRUE & (ACMG_PM2_1 == TRUE | ACMG_PM2_2 == TRUE),TRUE,FALSE,FALSE))



  # ## Assign logical ACMG evidence indicators
  # # BA1 -  high population germline frequency (1000G/gnomAD)
  # if('GLOBAL_AF_1KG' %in% colnames(cpg_calls) & 'GLOBAL_AF_GNOMAD' %in% colnames(cpg_calls)){
  #   cpg_calls <- cpg_calls %>% dplyr::mutate(BA1 = dplyr::if_else(GLOBAL_AF_GNOMAD > 0.05 | GLOBAL_AF_1KG > 0.05,TRUE,FALSE,FALSE))
  #   cpg_calls <- cpg_calls %>% dplyr::mutate(BA1 = dplyr::if_else((SYMBOL == 'HFE' | SYMBOL == 'SERPINA1') & BA1 == TRUE & !is.na(GLOBAL_AF_1KG) & GLOBAL_AF_1KG < 0.25,FALSE,BA1))
  #   cpg_calls <- cpg_calls %>% dplyr::mutate(BA1 = dplyr::if_else((SYMBOL == 'HFE' | SYMBOL == 'SERPINA1') & BA1 == TRUE & !is.na(GLOBAL_AF_GNOMAD) & GLOBAL_AF_GNOMAD < 0.25,FALSE,BA1))
  # }

  ## Assign logical ACMG evidence indicator
  # PM4 - Protein length changes (in non-repetitive regions) due to inframe indels or nonstop variant of genes that harbor variants with a dominant mode of inheritance - PM4
  # PCC1 - Protein length changes (in non-repetitive regions) due to inframe indels or nonstop variant of genes that
  #        harbor variants with a recessive mode of inheritance (and unknown MOI) - PPC1
  if("RMSK_HIT" %in% colnames(cpg_calls)){
    cpg_calls <- cpg_calls %>%
      dplyr::mutate(ACMG_PM4 = dplyr::if_else(stringr::str_detect(CONSEQUENCE,'stop_lost|inframe_deletion|inframe_insertion') &
                                           is.na(RMSK_HIT) & cpsr_gene_moi == 'AD',TRUE,FALSE,FALSE))
    cpg_calls <- cpg_calls %>%
      dplyr::mutate(ACMG_PPC1 = dplyr::if_else(stringr::str_detect(CONSEQUENCE,'stop_lost|inframe_deletion|inframe_insertion') &
                                           is.na(RMSK_HIT) & (cpsr_gene_moi == 'AR' | is.na(cpsr_gene_moi)),TRUE,FALSE,FALSE))
  }

  ## Assign logical ACMG evidence indicator
  # ACMG_PP2 - Missense variant in a gene that has a relatively low rate of benign missense variation and where missense variants are a common mechanism of disease
  cpg_calls <- cpg_calls %>%
    dplyr::mutate(ACMG_PP2 = dplyr::if_else(BENIGN_MISSENSE_RATE <= 0.2 & PATH_TRUNCATION_RATE < 0.5 & stringr::str_detect(CONSEQUENCE,"^missense_variant"),TRUE,FALSE,FALSE))

  ## Assign a logical ACMG evidence indicator
  # ACMG_BP1 - Missense variant in a gene for which primarily truncating variants (> 90%, as given in Maxwell et al.) are known to cause disease
  cpg_calls <- cpg_calls %>%
    dplyr::mutate(ACMG_BP1 = dplyr::if_else(PATH_TRUNCATION_RATE > 0.90 & stringr::str_detect(CONSEQUENCE,'^missense_variant'),TRUE,FALSE,FALSE))


  ## Assign logical ACMG evidence indicators
  # ACMG_PS1 - coinciding with known pathogenic missense variants (yet with different nucleotide change)
  # ACMG_PM5 - occurs at the same codon as a known pathogenic missense variant
  # ACMG_BSC1 - coinciding with known benign missense variants
  # ACMG_BMC1 - occurs at the same codon as a known benign missense variant

  cpg_calls
  cpg_calls$codon_prefix <- NA
  if(nrow(cpg_calls[!is.na(cpg_calls$CONSEQUENCE) & cpg_calls$CONSEQUENCE == 'missense_variant',]) > 0){
    cpg_calls[cpg_calls$CONSEQUENCE == 'missense_variant',]$codon_prefix <-
      stringr::str_match(cpg_calls[cpg_calls$CONSEQUENCE == 'missense_variant',]$HGVSp_short,"p\\.[A-Z]{1}[0-9]{1,}")
  }

  if(nrow(cpg_calls[!is.na(cpg_calls$codon_prefix),]) > 0){
    cpg_calls_pathogenic_codon <-
      dplyr::left_join(dplyr::filter(dplyr::select(cpg_calls,VAR_ID,codon_prefix,SYMBOL), !is.na(codon_prefix)), pcgr_data[['clinvar']][['cpg_loci']][['high_confidence']][['pathogenic']][['codon']], by=c("codon_prefix" = "codon_prefix","SYMBOL" = "symbol"))
    cpg_calls <- cpg_calls %>%
      dplyr::left_join(dplyr::select(cpg_calls_pathogenic_codon, VAR_ID, clinvar_pathogenic_codon), by=c("VAR_ID"))

    cpg_calls_benign_codon <-
      dplyr::left_join(dplyr::filter(dplyr::select(cpg_calls,VAR_ID,codon_prefix,SYMBOL), !is.na(codon_prefix)), pcgr_data[['clinvar']][['cpg_loci']][['high_confidence']][['benign']][['codon']], by=c("codon_prefix" = "codon_prefix","SYMBOL" = "symbol"))
    cpg_calls <- cpg_calls %>%
      dplyr::left_join(dplyr::select(cpg_calls_benign_codon, VAR_ID, clinvar_benign_codon), by=c("VAR_ID"))

    cpg_calls <- cpg_calls %>%
      dplyr::mutate(ACMG_PM5 = dplyr::if_else(clinvar_pathogenic_codon == TRUE,TRUE,FALSE,FALSE))
    cpg_calls <- cpg_calls %>%
      dplyr::mutate(ACMG_BMC1 = dplyr::if_else(clinvar_benign_codon == TRUE,TRUE,FALSE,FALSE))
  }else{
    cpg_calls$ACMG_PM5 <- FALSE
    cpg_calls$ACMG_BMC1 <- FALSE
    cpg_calls$clinvar_pathogenic_codon <- NA
    cpg_calls$clinvar_benign_codon <- NA
  }


  if(nrow(cpg_calls[!is.na(cpg_calls$HGVSp_short),]) > 0){
    cpg_calls_pathogenic_hgvsp <-
      dplyr::left_join(dplyr::filter(dplyr::select(cpg_calls, VAR_ID, HGVSp_short, SYMBOL), !is.na(HGVSp_short)), pcgr_data[['clinvar']][['cpg_loci']][['high_confidence']][['pathogenic']][['peptide_change']], by=c("HGVSp_short" = "hgvs_p","SYMBOL" = "symbol"))
    cpg_calls <- cpg_calls %>%
      dplyr::left_join(dplyr::select(cpg_calls_pathogenic_hgvsp, VAR_ID, clinvar_pathogenic), by=c("VAR_ID"))

    cpg_calls_benign_hgvsp <-
      dplyr::left_join(dplyr::filter(dplyr::select(cpg_calls, VAR_ID, HGVSp_short, SYMBOL), !is.na(HGVSp_short)), pcgr_data[['clinvar']][['cpg_loci']][['high_confidence']][['benign']][['peptide_change']], by=c("HGVSp_short" = "hgvs_p","SYMBOL" = "symbol"))
    cpg_calls <- cpg_calls %>%
      dplyr::left_join(dplyr::select(cpg_calls_benign_hgvsp, VAR_ID, clinvar_benign), by=c("VAR_ID"))

    cpg_calls <- cpg_calls %>% dplyr::mutate(ACMG_PS1 = dplyr::if_else(clinvar_pathogenic == TRUE,TRUE,FALSE,FALSE))
    cpg_calls <- cpg_calls %>% dplyr::mutate(ACMG_BSC1 = dplyr::if_else(clinvar_benign == TRUE,TRUE,FALSE,FALSE))
  }else{
    cpg_calls$ACMG_PS1 <- FALSE
    cpg_calls$ACMG_BSC1 <- FALSE
    cpg_calls$clinvar_pathogenic <- NA
    cpg_calls$clinvar_benign <- NA
  }

  ## if previously found coinciding with pathogenic variant (ACMG_PS1), set ACMG_PM5 to false
  cpg_calls <- cpg_calls %>% dplyr::mutate(ACMG_PM5 = dplyr::case_when(ACMG_PM5 == T & ACMG_PS1 == T ~ FALSE, TRUE ~ as.logical(ACMG_PM5)))
  ## if previously found coinciding with benign variant (ACMG_BSC1), set ACMG_BMC1 to false
  cpg_calls <- cpg_calls %>% dplyr::mutate(ACMG_BMC1 = dplyr::case_when(ACMG_BMC1 == T & ACMG_BSC1 == T ~ FALSE, TRUE ~ as.logical(ACMG_BMC1)))
  cpg_calls <- cpg_calls %>% dplyr::select(-c(gad_an_ac_sufficient,clinvar_pathogenic_codon, clinvar_benign_codon, clinvar_pathogenic, clinvar_benign, cpsr_gene_moi, gad_af))

  ##Assign logical ACMG level
  # PM1 - missense variant in a somatic mutation hotspot as determined by cancerhotspots.org (v2)
  cpg_calls <- cpg_calls %>% tidyr::separate(MUTATION_HOTSPOT, c("hotspot_symbol","hotspot_codon","hotspot_pvalue"),sep="\\|",remove=F,extra="drop")
  if(nrow(cpg_calls[!is.na(cpg_calls$hotspot_codon),]) > 0){
    cpg_calls[!is.na(cpg_calls$hotspot_codon),]$hotspot_codon <- paste0('p.',cpg_calls[!is.na(cpg_calls$hotspot_codon),]$hotspot_codon)
  }

  cpg_calls <- cpg_calls %>%
    dplyr::mutate(ACMG_PM1 = dplyr::if_else(!is.na(hotspot_codon) & !is.na(hotspot_symbol) & !is.na(codon_prefix) &
                                              SYMBOL == hotspot_symbol & hotspot_codon == codon_prefix,TRUE,FALSE))

  return(cpg_calls)
}


determine_pathogenicity_classification <- function(cpg_calls){

  evidence_codes <- pcgrr::acmg_evidence_codes

  path_cols <- c('CPSR_CLASSIFICATION','CPSR_CLASSIFICATION_DOC','CPSR_CLASSIFICATION_CODE',
                'cpsr_score_pathogenic','cpsr_score_benign')
  cpg_calls <- cpg_calls[, !(colnames(cpg_calls) %in% path_cols)]

  cpg_calls$CPSR_CLASSIFICATION <- 'VUS'
  cpg_calls$CPSR_CLASSIFICATION_DOC <- ''
  cpg_calls$CPSR_CLASSIFICATION_CODE <- ''
  cpg_calls$cpsr_score_pathogenic <- 0
  cpg_calls$cpsr_score_benign <- 0

  i <- 1
  while(i <= nrow(evidence_codes)){
    category <- evidence_codes[i,]$category
    pole <- evidence_codes[i,]$pathogenicity_pole
    description <- evidence_codes[i,]$description
    cpsr_evidence_code <- evidence_codes[i,]$cpsr_evidence_code
    score <- evidence_codes[i,]$path_score
    if(cpsr_evidence_code %in% colnames(cpg_calls)){
      cpg_calls <- cpg_calls %>%
        dplyr::mutate(cpsr_score_benign = cpsr_score_benign + dplyr::if_else(pole == 'B' & !!rlang::sym(cpsr_evidence_code) == T,score,0)) %>%
        dplyr::mutate(cpsr_score_pathogenic = cpsr_score_pathogenic + dplyr::if_else(pole == 'P' & !!rlang::sym(cpsr_evidence_code) == T,score,0)) %>%
        dplyr::mutate(CPSR_CLASSIFICATION_DOC = paste0(CPSR_CLASSIFICATION_DOC,dplyr::if_else(!!rlang::sym(cpsr_evidence_code) == T,paste0("- ",description),""),sep="<br>")) %>%
        dplyr::mutate(CPSR_CLASSIFICATION_CODE = paste0(CPSR_CLASSIFICATION_CODE,dplyr::if_else(!!rlang::sym(cpsr_evidence_code) == T,cpsr_evidence_code,""),sep="|"))
    }
    i <- i + 1
  }

  cpg_calls <- cpg_calls %>%
    dplyr::mutate(CPSR_CLASSIFICATION_CODE = stringr::str_replace_all(stringr::str_replace_all(CPSR_CLASSIFICATION_CODE,"(\\|{2,})","|"),"(^\\|)|(\\|$)","")) %>%
    dplyr::mutate(CPSR_CLASSIFICATION_DOC = stringr::str_replace_all(stringr::str_replace_all(CPSR_CLASSIFICATION_DOC,"(<br>){2,}","<br>"),"(^(<br>))|((<br>)$)","")) %>%
    dplyr::mutate(cpsr_score_pathogenic = dplyr::if_else(stringr::str_detect(CPSR_CLASSIFICATION_CODE,"ACMG_PVS") & stringr::str_detect(CPSR_CLASSIFICATION_CODE,"ACMG_PM2_2"),cpsr_score_pathogenic - 1,cpsr_score_pathogenic)) %>%
    dplyr::mutate(cpsr_score_pathogenic = dplyr::if_else(stringr::str_detect(CPSR_CLASSIFICATION_CODE,"ACMG_PVS") & stringr::str_detect(CPSR_CLASSIFICATION_CODE,"ACMG_PM2_1"),cpsr_score_pathogenic - 0.5,cpsr_score_pathogenic)) %>%
    dplyr::mutate(CPSR_CLASSIFICATION = dplyr::if_else(cpsr_score_pathogenic >= 3.5 & cpsr_score_pathogenic < 5,"Likely_Pathogenic", CPSR_CLASSIFICATION)) %>%
    dplyr::mutate(CPSR_CLASSIFICATION = dplyr::if_else(cpsr_score_pathogenic >= 5,"Pathogenic", CPSR_CLASSIFICATION)) %>%
    dplyr::mutate(CPSR_CLASSIFICATION = dplyr::if_else(cpsr_score_benign <= -3 & cpsr_score_pathogenic <= 0.5,"Likely_Benign", CPSR_CLASSIFICATION)) %>%
    dplyr::mutate(CPSR_CLASSIFICATION = dplyr::if_else(cpsr_score_benign <= -5,"Benign", CPSR_CLASSIFICATION)) %>%
    dplyr::mutate(CPSR_CLASSIFICATION_SCORE = dplyr::if_else(cpsr_score_benign == 0, cpsr_score_pathogenic, cpsr_score_benign)) %>%
    dplyr::mutate(CPSR_CLASSIFICATION_SCORE = dplyr::if_else(cpsr_score_benign < 0 & cpsr_score_pathogenic > 0,cpsr_score_benign + cpsr_score_pathogenic,CPSR_CLASSIFICATION_SCORE)) %>%
    dplyr::select(-c(cpsr_score_benign,cpsr_score_pathogenic))

  return(cpg_calls)

}

#' Function that assign variants to different tiers for prioritisation of germline variants
#'
#' @param cpg_calls data frame with variants in predisposition_genes
#'
assign_cpsr_tier <- function(cpg_calls, cpsr_config, display_tags){

  predispose_tsv_tags <- c("GENOMIC_CHANGE","VAR_ID","GENOTYPE",
                           "GENOME_VERSION","VCF_SAMPLE_ID","VARIANT_CLASS","CODING_STATUS",
                           "SYMBOL","GENE_NAME","CCDS","ENTREZ_ID","UNIPROT_ID","ENSEMBL_GENE_ID","ENSEMBL_TRANSCRIPT_ID","REFSEQ_MRNA",
                           "ONCOGENE", "TUMOR_SUPPRESSOR","MOD","CONSEQUENCE","VEP_ALL_CSQ","CIVIC_ID","CIVIC_ID_2",
                           "PROTEIN_CHANGE","PROTEIN_DOMAIN", "HGVSp","HGVSc","LAST_EXON","CDS_CHANGE","MUTATION_HOTSPOT","RMSK_HIT",
                           "PROTEIN_FEATURE","EFFECT_PREDICTIONS", "LOSS_OF_FUNCTION", "DBSNP","CLINVAR_CLASSIFICATION",
                           "CLINVAR_MSID","CLINVAR_VARIANT_ORIGIN","CLINVAR_CONFLICTED", "CLINVAR_PHENOTYPE","CLINVAR_REVIEW_STATUS_STARS",
                           "N_INSILICO_CALLED","N_INSILICO_DAMAGING","N_INSILICO_TOLERATED",
                           "N_INSILICO_SPLICING_NEUTRAL","N_INSILICO_SPLICING_AFFECTED","GLOBAL_AF_GNOMAD",
                           cpsr_config[['popgen']][['vcftag_gnomad']])

  evidence_codes <- pcgrr::acmg_evidence_codes %>%
    dplyr::filter(cpsr_evidence_code != 'ACMG_BS2_1' & cpsr_evidence_code != 'ACMG_BS2_2' & cpsr_evidence_code != 'ACMG_BS2_3')
  predispose_tsv_tags <- c(predispose_tsv_tags, evidence_codes$cpsr_evidence_code,
                           c("CPSR_CLASSIFICATION","CPSR_CLASSIFICATION_SCORE","CPSR_CLASSIFICATION_CODE","CPSR_CLASSIFICATION_DOC","SOURCE"))

  rlogging::message("Generating tiered set of result variants for output in tab-separated values (TSV) file")

  snv_indel_report <- pcgrr::init_pcg_report(config = cpsr_config, class = "snv_indel", type = "predisposition")


  predispose_tags <- predispose_tsv_tags
  if(!is.null(cpsr_config[['custom_tags']])){
    if(cpsr_config[['custom_tags']][['custom_tags']] != ""){
      tags <- stringr::str_split(cpsr_config[['custom_tags']][['custom_tags']],pattern = ",")[[1]]
      for(t in tags){
        t <- stringr::str_trim(t)
        predispose_tags <- c(predispose_tags,t)
      }
    }
  }

  snv_indel_report[['variant_set']][['class5']] <- cpg_calls %>%
    dplyr::filter(!is.na(CLINVAR_CLASSIFICATION) & CLINVAR_CLASSIFICATION == "Pathogenic") %>%
    dplyr::mutate(SOURCE = "ClinVar")

  snv_indel_report[['variant_set']][['class4']] <- cpg_calls %>%
    dplyr::filter(!is.na(CLINVAR_CLASSIFICATION) & CLINVAR_CLASSIFICATION == "Likely_Pathogenic") %>%
    dplyr::mutate(SOURCE = "ClinVar")

  snv_indel_report[['variant_set']][['class3']] <- cpg_calls %>%
    dplyr::filter(!is.na(CLINVAR_CLASSIFICATION) & CLINVAR_CLASSIFICATION == "VUS") %>%
    dplyr::mutate(SOURCE = "ClinVar")

  snv_indel_report[['variant_set']][['class2']] <- cpg_calls %>%
    dplyr::filter(!is.na(CLINVAR_CLASSIFICATION) & CLINVAR_CLASSIFICATION == "Likely_Benign") %>%
    dplyr::mutate(SOURCE = "ClinVar")

  snv_indel_report[['variant_set']][['class1']] <- cpg_calls %>%
    dplyr::filter(!is.na(CLINVAR_CLASSIFICATION) & CLINVAR_CLASSIFICATION == "Benign") %>%
    dplyr::mutate(SOURCE = "ClinVar")

  ## identify remaining calls not registered in clinvar
  all_clinvar_calls <- data.frame()
  for(c in c('class1','class2','class3','class4','class5')){
    all_clinvar_calls <- dplyr::bind_rows(all_clinvar_calls, dplyr::select(snv_indel_report[['variant_set']][[c]], VAR_ID))
  }
  cpg_calls <- cpg_calls %>% dplyr::anti_join(all_clinvar_calls, by=c("VAR_ID"))

  n_nonclinvar = nrow(cpg_calls)

  cpg_calls <- cpg_calls %>%
    dplyr::filter(is.na(GLOBAL_AF_GNOMAD) | GLOBAL_AF_GNOMAD < cpsr_config[['maf_limits']][['maf_gnomad']])
  n_maf_filtered <- n_nonclinvar - nrow(cpg_calls)
  rlogging::message(paste0("Ignoring n = ",n_maf_filtered," unclassified variants with a global MAF frequency above ",cpsr_config[['maf_limits']][['maf_gnomad']]))

  n_after_maf_filtering <- nrow(cpg_calls)

  cpg_calls <- cpg_calls %>%
    dplyr::filter(CODING_STATUS == "coding")

  n_noncoding_filtered <- n_after_maf_filtering - nrow(cpg_calls)
  rlogging::message(paste0("Ignoring n = ",n_noncoding_filtered," unclassified variants with a noncoding variant consequence"))

  non_clinvar_calls <- list()
  non_clinvar_calls[['class5']] <- cpg_calls %>%
    dplyr::filter(CPSR_CLASSIFICATION == "Pathogenic") %>%
    dplyr::mutate(SOURCE = "Other")

  non_clinvar_calls[['class4']] <- cpg_calls %>%
    dplyr::filter(CPSR_CLASSIFICATION == "Likely_Pathogenic") %>%
    dplyr::mutate(SOURCE = "Other")

  non_clinvar_calls[['class3']] <- cpg_calls %>%
    dplyr::filter(CPSR_CLASSIFICATION == "VUS") %>%
    dplyr::mutate(SOURCE = "Other")

  non_clinvar_calls[['class2']] <- cpg_calls %>%
    dplyr::filter(CPSR_CLASSIFICATION == "Likely_Benign") %>%
    dplyr::mutate(SOURCE = "Other")

  non_clinvar_calls[['class1']] <- cpg_calls %>%
    dplyr::filter(CPSR_CLASSIFICATION == "Benign") %>%
    dplyr::mutate(SOURCE = "Other")

  for(c in c('class1','class2','class3','class4','class5')){
    snv_indel_report[['variant_set']][[c]] <- dplyr::bind_rows(non_clinvar_calls[[c]],snv_indel_report[['variant_set']][[c]])

    if(cpsr_config[['classification']][['clinvar_cpsr']] == F){
      snv_indel_report[['variant_set']][[c]] <- snv_indel_report[['variant_set']][[c]] %>%
        dplyr::mutate(CPSR_CLASSIFICATION = dplyr::if_else(!is.na(CLINVAR_CLASSIFICATION),"",as.character(CPSR_CLASSIFICATION))) %>%
        dplyr::mutate(CPSR_CLASSIFICATION_DOC = dplyr::if_else(!is.na(CLINVAR_CLASSIFICATION),"",as.character(CPSR_CLASSIFICATION_DOC))) %>%
        dplyr::mutate(CPSR_CLASSIFICATION_CODE = dplyr::if_else(!is.na(CLINVAR_CLASSIFICATION),"",as.character(CPSR_CLASSIFICATION_CODE))) %>%
        dplyr::mutate(CPSR_CLASSIFICATION_SCORE = dplyr::if_else(!is.na(CLINVAR_CLASSIFICATION),as.numeric(NA),as.numeric(CPSR_CLASSIFICATION_SCORE)))
    }
    snv_indel_report[['variant_set']][[c]] <- snv_indel_report[['variant_set']][[c]] %>%
      dplyr::arrange(SOURCE, CANCER_PHENOTYPE, desc(CPSR_CLASSIFICATION_SCORE)) %>%
      dplyr::select(-CANCER_PHENOTYPE)

    snv_indel_report[['variant_display']][[c]] <- dplyr::select(snv_indel_report[['variant_set']][[c]], dplyr::one_of(display_tags))
    snv_indel_report[['variant_set']][[c]] <- dplyr::select(snv_indel_report[['variant_set']][[c]], dplyr::one_of(predispose_tags))

    snv_indel_report[['variant_set']][[c]]$DBSNP <- unlist(lapply(stringr::str_match_all(snv_indel_report[['variant_set']][[c]]$DBSNP,">rs[0-9]{1,}<"),paste,collapse=","))
    snv_indel_report[['variant_set']][[c]]$DBSNP <- stringr::str_replace_all(snv_indel_report[['variant_set']][[c]]$DBSNP,">|<", "")
    snv_indel_report[['variant_set']][[c]]$GENE_NAME <- unlist(lapply(stringr::str_match_all(snv_indel_report[['variant_set']][[c]]$GENE_NAME,">.+<"),paste,collapse=","))
    snv_indel_report[['variant_set']][[c]]$GENE_NAME <- stringr::str_replace_all(snv_indel_report[['variant_set']][[c]]$GENE_NAME,">|<", "")
    snv_indel_report[['variant_set']][[c]]$PROTEIN_DOMAIN <- unlist(lapply(stringr::str_match_all(snv_indel_report[['variant_set']][[c]]$PROTEIN_DOMAIN,">.+<"),paste,collapse=","))
    snv_indel_report[['variant_set']][[c]]$PROTEIN_DOMAIN <- stringr::str_replace_all(snv_indel_report[['variant_set']][[c]]$PROTEIN_DOMAIN,">|<", "")
    snv_indel_report[['variant_set']][[c]]$CPSR_CLASSIFICATION_DOC <- stringr::str_replace_all(snv_indel_report[['variant_set']][[c]]$CPSR_CLASSIFICATION_DOC,"<br>-", ",")
    snv_indel_report[['variant_set']][[c]]$CPSR_CLASSIFICATION_DOC <- stringr::str_replace_all(snv_indel_report[['variant_set']][[c]]$CPSR_CLASSIFICATION_DOC,"^, ","")

    snv_indel_report[['variant_set']][[c]] <- snv_indel_report[['variant_set']][[c]] %>%
      dplyr::select(c("GENOMIC_CHANGE","VAR_ID","GENOTYPE","SOURCE","GENOME_VERSION","VCF_SAMPLE_ID","VARIANT_CLASS","CODING_STATUS",
                      "SYMBOL","GENE_NAME","CCDS","ENTREZ_ID","UNIPROT_ID","ENSEMBL_GENE_ID","ENSEMBL_TRANSCRIPT_ID","REFSEQ_MRNA",
                      "ONCOGENE", "TUMOR_SUPPRESSOR","MOD","CONSEQUENCE","VEP_ALL_CSQ",
                      "PROTEIN_CHANGE","PROTEIN_DOMAIN", "DBSNP","HGVSp","HGVSc","LAST_EXON","CDS_CHANGE","MUTATION_HOTSPOT","RMSK_HIT",
                      "PROTEIN_FEATURE","EFFECT_PREDICTIONS", "LOSS_OF_FUNCTION", "DBSNP","CLINVAR_CLASSIFICATION",
                      "CLINVAR_MSID","CLINVAR_VARIANT_ORIGIN","CLINVAR_CONFLICTED", "CLINVAR_PHENOTYPE","CLINVAR_REVIEW_STATUS_STARS"),
                    dplyr::everything())

    for(col in colnames(snv_indel_report[['variant_set']][[c]])){
      if(nrow(snv_indel_report[['variant_set']][[c]][!is.na(snv_indel_report[['variant_set']][[c]][,col]) & snv_indel_report[['variant_set']][[c]][,col] == "",]) > 0){
        snv_indel_report[['variant_set']][[c]][!is.na(snv_indel_report[['variant_set']][[c]][,col]) & snv_indel_report[['variant_set']][[c]][,col] == "",col] <- NA
      }
    }

    population_tags <- unique(c("GLOBAL_AF_GNOMAD", cpsr_config[["popgen"]][["vcftag_gnomad"]]))
    for (tag in population_tags){
      if(tag %in% colnames(snv_indel_report[["variant_display"]][[c]])){
        if (nrow(snv_indel_report[["variant_display"]][[c]][is.na(snv_indel_report[["variant_display"]][[c]][,tag]), ]) > 0){
          snv_indel_report[["variant_display"]][[c]][is.na(snv_indel_report[["variant_display"]][[c]][, tag]), tag] <- 0.00
        }
      }
    }
  }



  snv_indel_report[['variant_set']][['tsv']] <- dplyr::bind_rows(snv_indel_report[['variant_set']][['class5']],
                                                                 snv_indel_report[['variant_set']][['class4']],
                                                                 snv_indel_report[['variant_set']][['class3']],
                                                                 snv_indel_report[['variant_set']][['class2']],
                                                                 snv_indel_report[['variant_set']][['class1']])


  return(snv_indel_report)
}

#' Function that generate snv/indel + coding/noncoding stats for a given variant set
#'
#' @param calls data frame with variants in predisposition_genes
#' @param name type of variant group
#'
variant_stats_report <- function(calls, name = "variant_statistic"){

  call_stats <- list()
  call_stats[[name]] <- list()
  for(n in c('n','n_snv','n_indel','n_coding','n_noncoding')){
    call_stats[[name]][[n]] <- 0
  }

  call_stats[[name]][['n']] <- calls %>% nrow()
  call_stats[[name]][['n_snv']] <- calls %>% dplyr::filter(VARIANT_CLASS == 'SNV') %>% nrow()
  call_stats[[name]][['n_indel']] <- calls %>% dplyr::filter(VARIANT_CLASS != 'SNV') %>% nrow()
  call_stats[[name]][['n_coding']] <- calls %>% dplyr::filter(CODING_STATUS == 'coding') %>% nrow()
  call_stats[[name]][['n_noncoding']] <- calls %>% dplyr::filter(CODING_STATUS != 'coding') %>% nrow()

  return(call_stats)
}

#' Function that retrieves variants in genes recommended for incidental (secondary) findings
#'
#' @param calls data frame with variants in predisposition_genes
#' @param medgen object with medgen
#'
retrieve_sf_calls <- function(calls, medgen){

  calls$CLINVAR_PHENOTYPE <- NA
  sf_calls <- calls %>%
      dplyr::filter(!is.na(CANCER_PREDISPOSITION_SOURCE) & CANCER_PREDISPOSITION_SOURCE == 'ACMG_SF20' & !is.na(CLINVAR_CLINICAL_SIGNIFICANCE)) %>%
      dplyr::rename(CLIN_SIGNIFICANCE = CLINVAR_CLASSIFICATION) %>%
      dplyr::filter(CLIN_SIGNIFICANCE == "Pathogenic" | CLIN_SIGNIFICANCE == "Likely_Pathogenic")

  if(nrow(sf_calls) > 0){
    sf_calls_per_trait <- tidyr::separate_rows(sf_calls, CLINVAR_MEDGEN_CUI, sep = ",") %>%
      dplyr::select(VAR_ID, CLINVAR_MEDGEN_CUI) %>%
      dplyr::left_join(medgen, by = c("CLINVAR_MEDGEN_CUI" = "cui")) %>%
      dplyr::rename(CLINVAR_PHENOTYPE = cui_name) %>%
      dplyr::distinct()

    ## check that phenotype of ClinVar-registered variant matches the phenotype given in secondary findings
    sf_calls <- dplyr::inner_join(dplyr::select(sf_calls, -CLINVAR_PHENOTYPE), dplyr::select(sf_calls_per_trait, VAR_ID), by=c("VAR_ID" = "VAR_ID"))
    #sf_calls <- dplyr::inner_join(dplyr::select(sf_calls, -CLINVAR_PHENOTYPE), dplyr::select(sf_calls_per_trait, VAR_ID, CLINVAR_MEDGEN_CUI), by=c("CANCER_SYNDROME_CUI" = "CLINVAR_MEDGEN_CUI", "VAR_ID" = "VAR_ID"))
  }

  return(sf_calls)

}

#' Function that retrieves variants in cancer predisposition genes linked to cancer-related conditions according to ClinVar
#'
#' @param calls data frame with variants in predisposition_genes
#' @param medgen object with medgen
#'
detect_cancer_traits_clinvar <- function(cpg_calls, medgen){
  if(nrow(cpg_calls) > 0 & "CLINVAR_MEDGEN_CUI" %in% colnames(cpg_calls) & "VAR_ID" %in% colnames(cpg_calls)){

    #cpg_calls$CLINVAR_PHENOTYPE <- NA
    cpg_calls_per_trait <- tidyr::separate_rows(cpg_calls, CLINVAR_MEDGEN_CUI, sep = ",") %>%
      dplyr::select(VAR_ID, CLINVAR_MEDGEN_CUI) %>%
      dplyr::left_join(medgen, by = c("CLINVAR_MEDGEN_CUI" = "cui")) %>%
      dplyr::distinct()

    variants_with_cancer_assoc <- dplyr::filter(cpg_calls_per_trait, cancer_phenotype == 1) %>%
      dplyr::select(VAR_ID) %>%
      dplyr::mutate(CANCER_PHENOTYPE = 1) %>%
      dplyr::distinct()

    cpg_calls_with_phenotype <- cpg_calls_per_trait %>%
      dplyr::group_by(VAR_ID) %>%
      dplyr::summarise(CLINVAR_PHENOTYPE = paste(unique(cui_name), collapse = "; "))

    cpg_calls <- dplyr::left_join(cpg_calls, cpg_calls_with_phenotype, by = c("VAR_ID"))
    if(nrow(variants_with_cancer_assoc) > 0){
      cpg_calls <- dplyr::left_join(cpg_calls, variants_with_cancer_assoc, by = c("VAR_ID"))
    }else{
      cpg_calls$CANCER_PHENOTYPE <- NA
    }
  }
  return(cpg_calls)
}

#' Function that makes a piechart showing the number of variants at each significance level
#'
#' @param variants_tsv data frame with variants in predisposition_genes
#' @param plot_type ClinVar or Other
#'
summary_donut_chart <- function(variants_tsv, plot_type = 'ClinVar'){

  significance_colors <- c("#9E0142", "#D53E4F", "#000000", "#78C679", "#077009")
  significance_levels <- c('Pathogenic','Likely_Pathogenic','VUS','Likely_Benign','Benign')

  title = 'ClinVar variants'
  p <- NULL

  if(nrow(variants_tsv) > 0){

    set_clinvar <- variants_tsv %>% dplyr::filter(!is.na(CLINVAR_CLASSIFICATION))
    set_other <- variants_tsv %>% dplyr::filter(nchar(CPSR_CLASSIFICATION) > 0 & is.na(CLINVAR_CLASSIFICATION))

    if((plot_type == 'ClinVar' & nrow(set_clinvar) > 0) | (plot_type != "ClinVar" & nrow(set_other) > 0)){

      m <- data.frame()

      if(plot_type == 'ClinVar'){
        if(nrow(set_clinvar) > 0){
          t <- paste0('n = ',nrow(set_clinvar))
          title <- bquote('ClinVar variants, '~bold(.(t)))
          m <- as.data.frame(set_clinvar %>%
            dplyr::group_by(CLINVAR_CLASSIFICATION) %>%
            dplyr::summarise(n = dplyr::n()) %>%
            dplyr::rename(level = CLINVAR_CLASSIFICATION)) %>%
            dplyr::mutate(level = factor(level, levels = significance_levels)) %>%
            dplyr::arrange(level) %>%
            dplyr::mutate(prop = as.numeric(n/sum(n))) %>%
            dplyr::mutate(lab.ypos = cumsum(prop) - 0.5*prop) %>%
            dplyr::mutate(n = as.character(n))
        }
      }else{
        if(nrow(set_other) > 0){
          t <- paste0('n = ',nrow(set_other))
          title <- bquote('Other variants, CPSR-classified, '~bold(.(t)))
          m <- as.data.frame(set_other %>%
            dplyr::group_by(CPSR_CLASSIFICATION) %>%
            dplyr::summarise(n = dplyr::n()) %>%
            dplyr::rename(level = CPSR_CLASSIFICATION)) %>%
            dplyr::mutate(level = factor(level, levels = significance_levels)) %>%
            dplyr::arrange(level) %>%
            dplyr::mutate(prop = as.numeric(n/sum(n))) %>%
            dplyr::mutate(lab.ypos = cumsum(prop) - 0.5*prop) %>%
            dplyr::mutate(n = as.character(n))

        }
      }


      p <- ggplot2::ggplot(m, ggplot2::aes(x = 2, y = prop, fill = level)) +
        ggplot2::geom_bar(stat = "identity", color = "white") +
        ggplot2::coord_polar(theta = "y", start = 0) +
        ggplot2::geom_text(ggplot2::aes(y = 1-lab.ypos, label = n), color = "white", family = "Helvetica", size = 6)+
        ggplot2::scale_fill_manual(values = significance_colors, labels = significance_levels, drop = F) +
        ggplot2::theme_void() +
        ggplot2::xlim(0.5, 2.5) +
        #ggplot2::ylim(0,1)+
        ggplot2::ggtitle(title) +
        ggplot2::theme(plot.title = ggplot2::element_text(family = "Helvetica", size = 16, vjust = -1, hjust = 0.5),
                       legend.title = ggplot2::element_blank(),
                       plot.margin = ggplot2::margin(0, 0, 0, 0, "cm"),
                       #plot.margin = grid::unit(c(0.0,0.0,0.0,0.0), "mm"),
                       legend.text = ggplot2::element_text(family = "Helvetica", size = 14))
    }


  }
  return(p)

}


get_x_coords <- function(box_w, n_col, space = 0.2){
  i <- 1
  x_coords <- c(0)
  pos <- 0
  while(i <= (n_col - 1)){
    x <- pos + box_w + space
    x_coords <- c(x_coords, x)
    pos <- x
    i <- i + 1
  }
  return(x_coords)
}

#' Function that makes a tile chart for a gene list
#'
#' @param genes vector of gene symbols
#' @param confidence vector of confidence level (1-4)
#'

gene_selection_tiles <- function(genes = NULL, confidence = NULL, box_w = 2, box_h = 0.5, col_numbers = 8, space = 0.2){


  tile_numbers <- list()
  tile_numbers[['full']] <- list()
  tile_numbers[['full']][['n_col']] <- 0
  tile_numbers[['full']][['n_row']] <- 0
  tile_numbers[['custom']] <- list()
  tile_numbers[['custom']][['n_col']] <- length(genes)
  tile_numbers[['custom']][['n_row']] <- 1

  color_df <- data.frame('color' = c("#b8b8ba","#d9534f","#f0ad4e","#3fad46","#000000"), stringsAsFactors = F)
  color_df$confidence <- seq(0,4,1)

  confidence_df <- data.frame('color' = rep('black',length(genes)), stringsAsFactors = F)
  if(!is.null(confidence) & length(confidence) == length(genes)){
    confidence_df <- data.frame('confidence' = as.integer(confidence), stringsAsFactors = F) %>%
      dplyr::left_join(color_df, by = "confidence")
  }


  y_step <- box_h + space

  if(length(genes) >= col_numbers){
    tile_numbers[['full']][['n_row']] <- as.integer(length(genes) / col_numbers)
    tile_numbers[['full']][['n_col']] <- col_numbers
    if(length(genes) %% col_numbers > 0){
      tile_numbers[['custom']][['n_col']] <- length(genes) %% col_numbers
    }else{
      tile_numbers[['custom']][['n_row']] <- 0
    }
  }

  max_y <- ((tile_numbers[['full']][['n_row']] + tile_numbers[['custom']][['n_row']]) - 1) * y_step

  df_all <- data.frame()
  j <- 1
  y_coord <- max_y

  if(tile_numbers[['full']][['n_row']] > 0){
    m <- 1
    while(m <= tile_numbers[['full']][['n_row']]){
      k <- col_numbers * m
      if(m > 1){
        y_coord <- max_y - (y_step * (m - 1))
      }
      #xc = get_x_coords(box_w = 2, n_col = col_numbers)
      df <- data.frame(
        x = pcgrr::get_x_coords(box_w = box_w, n_col = col_numbers, space = space),
        y = rep(y_coord,col_numbers),
        h = rep(box_h, col_numbers),
        w = rep(box_w, col_numbers),
        info = genes[j:k], stringsAsFactors = F
      )
      df_all <- dplyr::bind_rows(df_all, df)
      m <- m + 1
      j <- j + col_numbers
    }
  }
  if(tile_numbers$custom$n_row > 0){
    df <- data.frame(
      x = pcgrr::get_x_coords(box_w = box_w, n_col = tile_numbers[['custom']][['n_col']], space = space),
      y = rep(y_coord - y_step,tile_numbers[['custom']][['n_col']]),
      h = rep(box_h, tile_numbers[['custom']][['n_col']]),
      w = rep(box_w, tile_numbers[['custom']][['n_col']]),
      info = genes[j:length(genes)], stringsAsFactors = F
    )
    df_all <- dplyr::bind_rows(df_all, df)
  }

  p <- ggplot2::ggplot(df_all, ggplot2::aes(x, y, height = h, width = w, label = info)) +
    ggplot2::geom_tile(colour = confidence_df$color, fill = confidence_df$color, size = 0.8) +
    ggplot2::geom_text(color = "white", fontface = "bold", size=6) +
    ggplot2::coord_fixed() +
    ggplot2::xlim(-(box_w / 2),(box_w * col_numbers) + 0.4) +
    ggplot2::theme_void() +
    ggplot2::theme(
      plot.margin = (grid::unit(c(0, 0, 0, 0), "mm"))
    )

  return(p)
}

#' Function that retrieves clinical evidence items (CIVIC, CBMDB) for somatic cancer variants
#'
#' @param snv_indel_report report structure
#' @param pcgr_data object with PCGR annotation data
#'
#' @return list

get_germline_biomarkers <- function(snv_indel_report, pcgr_data){

  rlogging::message('Matching variant set with existing genomic biomarkers from CIViC (germline)')

  clin_eitems_list <- list()
  for(type in c('diagnostic','prognostic','predictive','predisposing')){
    clin_eitems_list[[type]] <- data.frame()
  }

  germline_biomarkers <- pcgr_data[['biomarkers']][['civic']] %>%
    dplyr::filter((alteration_type == 'MUT' | alteration_type == 'MUT_LOF') & (variant_origin == "Germline Mutation" | variant_origin == 'Germline Polymorphism')) %>%
    dplyr::filter(mapping_category == 'exact' | mapping_category == 'codon' | mapping_category == 'gene') %>%
    dplyr::select(genesymbol, alteration_type, clinical_significance, evidence_id, evidence_direction, pubmed_html_link, rating,
                  therapeutic_context, cancer_type, evidence_description, evidence_level, evidence_type, mapping_category,
                  eitem_consequence, eitem_exon, eitem_codon) %>%
      dplyr::rename(citation = pubmed_html_link, description = evidence_description)


  clinical_evidence_items <- data.frame()
  bmarker_mapping_levels <- c('exact','gene')

  for(mapping in bmarker_mapping_levels){
    bmarkers <- germline_biomarkers %>% dplyr::filter(mapping_category == mapping)

    if(mapping == 'exact'){
      sample_calls_civic <- snv_indel_report[['variant_set']][['tsv']] %>% dplyr::filter(!is.na(CIVIC_ID))
      if(nrow(sample_calls_civic) > 0){
        eitems <- dplyr::select(sample_calls_civic,CIVIC_ID,VAR_ID, SYMBOL) %>%
          dplyr::rename(evidence_id = CIVIC_ID) %>%
          tidyr::separate_rows(evidence_id,sep=",") %>%
          dplyr::inner_join(bmarkers, by=c("evidence_id" = "evidence_id", "SYMBOL" = "genesymbol")) %>%
          dplyr::distinct() %>%
          dplyr::rename_all(toupper) %>%
          dplyr::select(-c(EITEM_CONSEQUENCE,MAPPING_CATEGORY,EITEM_CODON,EITEM_EXON))

        if(nrow(eitems) > 0){
          clinical_evidence_items <- dplyr::bind_rows(clinical_evidence_items, eitems) %>%
            dplyr::mutate(BIOMARKER_MAPPING = mapping)
        }
      }
    }
    else{
      sample_calls_civic <- snv_indel_report[['variant_set']][['tsv']] %>%
        dplyr::filter(!is.na(CIVIC_ID_2)) %>%
        dplyr::filter(!is.na(CPSR_CLASSIFICATION) | !is.na(CLINVAR_CLASSIFICATION)) %>%
        dplyr::filter(stringr::str_detect(CPSR_CLASSIFICATION,"Pathogenic") | stringr::str_detect(CLINVAR_CLASSIFICATION,"Pathogenic"))

      if(nrow(sample_calls_civic) > 0){
        eitems <- dplyr::select(sample_calls_civic,CIVIC_ID_2,VAR_ID, SYMBOL, LOSS_OF_FUNCTION) %>%
          dplyr::mutate(BIOMARKER_MAPPING = mapping) %>%
          tidyr::separate_rows(CIVIC_ID_2,sep=",") %>%
          dplyr::rename(evidence_id = CIVIC_ID_2) %>%
          dplyr::inner_join(bmarkers, by=c("evidence_id" = "evidence_id", "SYMBOL" = "genesymbol")) %>%
          dplyr::filter((alteration_type == 'MUT_LOF' & LOSS_OF_FUNCTION == T) | alteration_type == 'MUT') %>%
          dplyr::rename_all(toupper) %>%
          dplyr::select(-c(EITEM_CONSEQUENCE,MAPPING_CATEGORY,EITEM_CODON,EITEM_EXON,LOSS_OF_FUNCTION)) %>%
          dplyr::distinct()

        if(nrow(eitems) > 0){
          clinical_evidence_items <- dplyr::bind_rows(clinical_evidence_items, eitems)
        }
      }
    }
  }

  snv_indel_report[['variant_set']][['tsv']] <- snv_indel_report[['variant_set']][['tsv']] %>%
    dplyr::select(-c(CIVIC_ID, CIVIC_ID_2))

  if(nrow(clinical_evidence_items) > 0){
    clinical_evidence_items <- clinical_evidence_items %>%
      dplyr::mutate(BIOMARKER_MAPPING = dplyr::if_else(is.na(BIOMARKER_MAPPING),"gene",as.character(BIOMARKER_MAPPING)))


    supporting_variants <- as.data.frame(
      clinical_evidence_items %>%
        dplyr::select(VAR_ID, EVIDENCE_ID) %>%
        dplyr::left_join(dplyr::select(snv_indel_report[['variant_set']][['tsv']], VAR_ID, CDS_CHANGE, GENE_NAME, GENOMIC_CHANGE, SYMBOL), by = c("VAR_ID")) %>%
        dplyr::group_by(EVIDENCE_ID) %>%
        dplyr::summarise(SYMBOL = paste(unique(SYMBOL),collapse=", "),
                         GENE_NAME = paste(unique(GENE_NAME), collapse=", "),
                         CDS_CHANGE = paste(unique(CDS_CHANGE), collapse = ", "),
                         GENOMIC_CHANGE = paste(unique(GENOMIC_CHANGE), collapse=", ")
                         )
    )


    clinical_evidence_items <- clinical_evidence_items %>%
      dplyr::select(-VAR_ID) %>%
      dplyr::distinct() %>%
      dplyr::left_join(supporting_variants, by=c("EVIDENCE_ID","SYMBOL"))

    rlogging::message(paste0(nrow(clinical_evidence_items),' clinical evidence item(s) found .. mapping = ',mapping))

  }else{
    rlogging::message(paste0(nrow(clinical_evidence_items),' clinical evidence item(s) found .. mapping = ',mapping))
  }


  if(nrow(clinical_evidence_items) > 0){
    for(type in c('prognostic','diagnostic','predictive','predisposing')){
      clin_eitems_list[[type]] <- clinical_evidence_items %>%
        dplyr::filter(EVIDENCE_TYPE == stringr::str_to_title(type)) %>%
        dplyr::arrange(EVIDENCE_LEVEL) %>%
        dplyr::select(SYMBOL, GENE_NAME, CANCER_TYPE, CLINICAL_SIGNIFICANCE, EVIDENCE_LEVEL, RATING, dplyr::everything()) %>%
        dplyr::select(-c(ALTERATION_TYPE,EVIDENCE_ID, EVIDENCE_TYPE))
    }
  }

  snv_indel_report$clinical_evidence_item <- clin_eitems_list
  return(snv_indel_report)

}

