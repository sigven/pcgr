
#' Function that generates predisposition_report - CPSR
#'
#' @param project_directory name of project directory
#' @param query_vcf2tsv name of gzipped TSV file (vcf2tsv) with annotated query SNVs/InDels
#' @param pcgr_data List of data frames with PCGR data annotations
#' @param pcgr_config Object with CPSR configuration parameters
#' @param sample_name sample identifier
#' @param pcgr_version PCGR software version
#' @param genome_assembly human genome assembly version
#'

generate_predisposition_report <- function(project_directory, query_vcf2tsv, pcgr_data, pcgr_config = NULL, sample_name = "SampleX",
                            pcgr_version = "0.1.0", genome_assembly = "grch37"){

  pcg_report <- pcgrr::init_pcg_report(pcgr_config, sample_name, pcgr_version, genome_assembly, class = NULL, pcgr_data = pcgr_data, type = "predisposition")

  genome_seq <- BSgenome.Hsapiens.UCSC.hg38
  assembly <- "hg38"
  if(genome_assembly == "grch37"){
    genome_seq <- BSgenome.Hsapiens.UCSC.hg19
    assembly <- "hg19"
  }

  fnames <- list()
  fnames[["tsv"]] <- paste0(project_directory, "/", sample_name, ".cpsr.snvs_indels.tiers.",genome_assembly,".tsv")
  fnames[["json"]] <- paste0(project_directory, "/", sample_name, ".cpsr.",genome_assembly, ".json")

  ## define tags/variables to display in data tables (tier 1,2,3A,3B)
  predispose_tier1_2_display <- c("SYMBOL", "CONSEQUENCE", "PROTEIN_CHANGE", "CLINVAR_PHENOTYPE", "GENOTYPE", "GENE_NAME", "PROTEIN_DOMAIN",
                          "HGVSp", "HGVSc", "CDS_CHANGE", "MUTATION_HOTSPOT", "RMSK_HIT","PROTEIN_FEATURE", "PREDICTED_EFFECT", "LOSS_OF_FUNCTION", "DBSNP", "CLINVAR",
                          "CLINVAR_CLINICAL_SIGNIFICANCE", "CLINVAR_VARIANT_ORIGIN", "GWAS_CITATION","ONCOGENE", "ONCOSCORE", "TUMOR_SUPPRESSOR",
                          "GLOBAL_AF_GNOMAD", pcg_report[["pcgr_config"]][["popgen"]][["vcftag_gnomad"]], "GLOBAL_AF_1KG",
                          pcg_report[["pcgr_config"]][["popgen"]][["vcftag_tgp"]], "GENOMIC_CHANGE", "GENOME_VERSION")

  predispose_tier3A_display <- unique(c("SYMBOL", "CONSEQUENCE", "PROTEIN_CHANGE", "CLINVAR_PHENOTYPE", "PATHRANK", "GENOTYPE",
                                        "PATHDOC", "PATHSCORE", predispose_tier1_2_display))
  predispose_tier3B_display <- predispose_tier3A_display[!predispose_tier3A_display %in% c("CLINVAR", "CLINVAR_PHENOTYPE",
                                                                                           "CLINVAR_CLINICAL_SIGNIFICANCE",
                                                                                           "CLINVAR_VARIANT_ORIGIN")]

  predispose_tier3B_display <- unique(c("SYMBOL","CONSEQUENCE", "PROTEIN_CHANGE","PATHRANK","PATHSCORE","GENOTYPE","PATHDOC", predispose_tier3B_display))

  predispose_gwas_display <- c("SYMBOL","CONSEQUENCE", "GWAS_CITATION","PROTEIN_CHANGE","GENOTYPE","LOSS_OF_FUNCTION",
                                 "PROTEIN_CHANGE","GENE_NAME", "GWAS_PHENOTYPE","PROTEIN_DOMAIN","HGVSp", "HGVSc", "CDS_CHANGE", "PROTEIN_FEATURE", "PREDICTED_EFFECT",
                                 "DBSNP","CLINVAR","CLINVAR_PHENOTYPE","CLINVAR_CLINICAL_SIGNIFICANCE","GLOBAL_AF_GNOMAD",
                                  pcg_report[["pcgr_config"]][["popgen"]][["vcftag_gnomad"]],
                                 "GLOBAL_AF_1KG",pcg_report[["pcgr_config"]][["popgen"]][["vcftag_tgp"]],"GENOMIC_CHANGE","GENOME_VERSION")



  if (query_vcf2tsv != "None.gz"){
    if (!file.exists(query_vcf2tsv) | file.size(query_vcf2tsv) == 0){
      rlogging::warning(paste0("File ", query_vcf2tsv, " does not exist or has zero size"))
    }else{
      if (!is.null(pcgr_config) & query_vcf2tsv != "None.gz"){

        ##read calls
        calls <- pcgrr::get_calls(query_vcf2tsv, pcgr_data, pcgr_version, sample_name, pcgr_config, genome_seq, assembly)
        rlogging::message("Filtering against cancer predisposition genes:")
        sample_calls <- dplyr::semi_join(calls, pcg_report[["snv_indel"]][["predisposition_genes"]], by = c("SYMBOL" = "symbol"))
        sample_calls <- dplyr::rename(sample_calls, CLINVAR_CLINICAL_SIGNIFICANCE = CLINVAR_CLNSIG)
        rlogging::message(nrow(sample_calls), " variants remaining")
        sample_calls_coding <- dplyr::filter(sample_calls, CODING_STATUS == "coding")
        gene_hits_coding <- paste(unique(sample_calls_coding$SYMBOL), collapse = ", ")
        rlogging::message("Found coding variants in the following cancer predisposition genes: ", gene_hits_coding)

        pcg_report$snv_indel$eval <- TRUE


        rlogging::message("Limiting variants to hereditary cancer-predisposing syndromes/cancer conditions")
        phenotype_medgen_oncology <- pcgr_data$phenotype_medgen_oncology %>%
          dplyr::filter(group == "Hereditary_Cancer_Syndrome_NOS" | group == "Hereditary_Cancer_Susceptibility_NOS") %>%
          dplyr::filter(!is.na(cui_name)) %>%
          dplyr::select(cui, cui_name) %>%
          dplyr::mutate(cancer_phenotype = 1) %>%
          dplyr::distinct()

        medgen <- pcgr_data$medgen_map %>%
          dplyr::left_join(phenotype_medgen_oncology, by = c("cui", "cui_name"))

        sample_calls_per_trait <- tidyr::separate_rows(sample_calls, CLINVAR_MEDGEN_CUI, sep = ",") %>%
          dplyr::select(VAR_ID, CLINVAR_MEDGEN_CUI) %>%
          dplyr::left_join(medgen, by = c("CLINVAR_MEDGEN_CUI" = "cui")) %>%
          dplyr::distinct()

        variants_with_cancer_assoc <- dplyr::filter(sample_calls_per_trait, cancer_phenotype == 1) %>%
          dplyr::select(VAR_ID) %>%
          dplyr::distinct()

        sample_calls_with_phenotype <- sample_calls_per_trait %>%
          dplyr::group_by(VAR_ID) %>%
          dplyr::summarise(CLINVAR_PHENOTYPE = paste(unique(cui_name), collapse = "; "))

        sample_calls <- dplyr::left_join(sample_calls, sample_calls_with_phenotype, by = c("VAR_ID"))
        sample_calls_all <- list()
        sample_calls_all[["cancer_phenotype"]] <- dplyr::inner_join(sample_calls, variants_with_cancer_assoc, by = c("VAR_ID"))
        sample_calls_all[["noncancer_phenotype"]] <- dplyr::anti_join(sample_calls, dplyr::select(sample_calls_all[["cancer_phenotype"]], VAR_ID), by = c("VAR_ID"))

        if (nrow(sample_calls_all[["cancer_phenotype"]]) > 0 | nrow(sample_calls_all[["noncancer_phenotype"]]) > 0){

          rlogging::message("Assignment of variants to tier 1/tier 2/tier 3")

          pathogenic_clnsig_regex <- "^(pathogenic|pathogenic,_other|pathogenic,_risk_factor|pathogenic,_drug_response|pathogenic,_other,_risk_factor)$"
          likely_pathogenic_clnsig_regex <- "^(likely_pathogenic|pathogenic,likely_pathogenic|pathogenic,likely_pathogenic,_risk_factor|likely_pathogenic,_other|likely_pathogenic,_risk_factor|pathogenic,likely_pathogenic,_other)$"
          vus_clnsig_regex <- "^(uncertain_significance|risk_factor|uncertain_significance,_other)$"
          all_hits <- data.frame()

          ## Assign TIER 1 variants to pcg_report object
          for (ph in c("cancer_phenotype", "noncancer_phenotype")){
            pcg_report[["snv_indel"]][["variant_display"]][["tier1"]][[ph]] <-  dplyr::filter(sample_calls_all[[ph]], CLINVAR_CONFLICTED == 0 & stringr::str_detect(CLINVAR_VARIANT_ORIGIN,"germline") & !is.na(CLINVAR_CLINICAL_SIGNIFICANCE) & stringr::str_detect(CLINVAR_CLINICAL_SIGNIFICANCE, pathogenic_clnsig_regex))
            if(nrow(pcg_report[["snv_indel"]][["variant_display"]][["tier1"]][[ph]]) > 0){
              pcg_report[["snv_indel"]][["variant_display"]][["tier1"]][[ph]] <-
                pcg_report[["snv_indel"]][["variant_display"]][["tier1"]][[ph]] %>%
                dplyr::mutate(PATHRANK = NA, PATHSCORE = NA, PATHDOC = NA)
              all_hits <- dplyr::bind_rows(all_hits, dplyr::select(pcg_report[["snv_indel"]][["variant_display"]][["tier1"]][[ph]],GENOMIC_CHANGE))
            }
            rlogging::message(paste0("TIER 1: Pathogenic variants - ",ph,": n = ",nrow(pcg_report[["snv_indel"]][["variant_display"]][["tier1"]][[ph]])))

          }

          ## Assign TIER 2 variants (likely pathogenic) to pcg_report object
          for (ph in c("cancer_phenotype", "noncancer_phenotype")){
            pcg_report[["snv_indel"]][["variant_display"]][["tier2"]][[ph]] <-  dplyr::filter(sample_calls_all[[ph]], CLINVAR_CONFLICTED == 0 & stringr::str_detect(CLINVAR_VARIANT_ORIGIN,"germline") & !is.na(CLINVAR_CLINICAL_SIGNIFICANCE) & stringr::str_detect(CLINVAR_CLINICAL_SIGNIFICANCE, likely_pathogenic_clnsig_regex))
            if(nrow(pcg_report[["snv_indel"]][["variant_display"]][["tier2"]][[ph]]) > 0){
              pcg_report[["snv_indel"]][["variant_display"]][["tier2"]][[ph]] <-
                pcg_report[["snv_indel"]][["variant_display"]][["tier2"]][[ph]] %>%
                dplyr::mutate(PATHRANK = NA, PATHSCORE = NA, PATHDOC = NA)
              all_hits <- dplyr::bind_rows(all_hits, dplyr::select(pcg_report[["snv_indel"]][["variant_display"]][["tier2"]][[ph]],GENOMIC_CHANGE))
            }
            rlogging::message(paste0("TIER 2: Likely pathogenic variants - ",ph,": n = ",nrow(pcg_report[["snv_indel"]][["variant_display"]][["tier2"]][[ph]])))
          }

          ## Assign TIER 3A (VUS in ClinVar) variants to pcg_report object
          for (ph in c("cancer_phenotype", "noncancer_phenotype")){
            pcg_report[["snv_indel"]][["variant_display"]][["tier3A"]][[ph]] <-  dplyr::filter(sample_calls_all[[ph]], CLINVAR_CONFLICTED == 0 & stringr::str_detect(CLINVAR_VARIANT_ORIGIN, "germline") & !is.na(CLINVAR_CLINICAL_SIGNIFICANCE) & stringr::str_detect(CLINVAR_CLINICAL_SIGNIFICANCE, vus_clnsig_regex))
            if(nrow(pcg_report[["snv_indel"]][["variant_display"]][["tier3A"]][[ph]]) > 0){
              pcg_report[["snv_indel"]][["variant_display"]][["tier3A"]][[ph]] <-
                pcgrr::assign_pathogenicity_score(pcg_report[["snv_indel"]][["variant_display"]][["tier3A"]][[ph]], pcgr_config, pcgr_data)
              pcg_report[["snv_indel"]][["variant_display"]][["tier3A"]][[ph]] <-
                pcg_report[["snv_indel"]][["variant_display"]][["tier3A"]][[ph]] %>%
                dplyr::arrange(desc(PATHSCORE), LOSS_OF_FUNCTION, CODING_STATUS, desc(ONCOSCORE))
              all_hits <- dplyr::bind_rows(all_hits, dplyr::select(pcg_report[["snv_indel"]][["variant_display"]][["tier3A"]][[ph]],GENOMIC_CHANGE))
            }
            rlogging::message(paste0("TIER 3: Variants of uncertain significance - ",ph,": n = ",nrow(pcg_report[["snv_indel"]][["variant_display"]][["tier3A"]][[ph]])))
          }

          ## Assign TIER 3B (Non-classified variants) to pcg_report object
          pcg_report[["snv_indel"]][["variant_display"]][["tier3B"]] <-  dplyr::filter(sample_calls_all[["noncancer_phenotype"]], is.na(CLINVAR_CLINICAL_SIGNIFICANCE) & CODING_STATUS == "coding" & (is.na(GLOBAL_AF_GNOMAD) | GLOBAL_AF_GNOMAD < pcgr_config$maf_limits$maf_gnomad) & (is.na(GLOBAL_AF_1KG) | GLOBAL_AF_1KG < pcgr_config$maf_limits$maf_tgp))
          if(nrow(pcg_report[["snv_indel"]][["variant_display"]][["tier3B"]]) > 0){
            pcg_report[["snv_indel"]][["variant_display"]][["tier3B"]] <-
              pcgrr::assign_pathogenicity_score(pcg_report[["snv_indel"]][["variant_display"]][["tier3B"]], pcgr_config, pcgr_data)
            pcg_report[["snv_indel"]][["variant_display"]][["tier3B"]] <-
              pcg_report[["snv_indel"]][["variant_display"]][["tier3B"]] %>%
              dplyr::arrange(desc(PATHSCORE), LOSS_OF_FUNCTION, CODING_STATUS, desc(ONCOSCORE))
            all_hits <- dplyr::bind_rows(all_hits, dplyr::select(pcg_report[["snv_indel"]][["variant_display"]][["tier3B"]],GENOMIC_CHANGE))
          }
          rlogging::message(paste0("TIER 3: Other unclassified variants - ",ph,": n = ",nrow(pcg_report[["snv_indel"]][["variant_display"]][["tier3B"]])))

          if(pcgr_config[['gwas']][['gwas_hits']] == TRUE){
            rlogging::message("Assignment of other variants to hits from genome-wide association studies")
          }
          ## Assign GWAS hits to pcg_report object
          pcg_report[["snv_indel"]][["variant_display"]][["gwas"]] <-  dplyr::filter(calls, !is.na(GWAS_HIT) & !is.na(GWAS_CITATION))
          pcg_report[["snv_indel"]][["variant_display"]][["gwas"]] <- dplyr::anti_join(pcg_report[["snv_indel"]][["variant_display"]][["gwas"]], all_hits, by=c("GENOMIC_CHANGE"))
          if(nrow(pcg_report[["snv_indel"]][["variant_display"]][["gwas"]]) > 0){
            pcg_report[["snv_indel"]][["variant_display"]][["gwas"]] <-
              pcg_report[["snv_indel"]][["variant_display"]][["gwas"]] %>%
              dplyr::mutate(PATHRANK = NA, PATHSCORE = NA, PATHDOC = NA, CLINVAR_PHENOTYPE = NA) %>%
              dplyr::rename(CLINVAR_CLINICAL_SIGNIFICANCE = CLINVAR_CLNSIG) %>%
              dplyr::arrange(LOSS_OF_FUNCTION, CODING_STATUS, desc(ONCOSCORE))
          }
          if(pcgr_config[['gwas']][['gwas_hits']] == TRUE){
            rlogging::message(paste0("GWAS hits - cancer phenotypes: n = ",nrow(pcg_report[["snv_indel"]][["variant_display"]][["gwas"]])))
          }

          pcg_report <- pcgrr::generate_tier_tsv_cpsr(pcg_report, sample_name = sample_name)
          pcg_report <- pcgrr::summary_findings_cpsr(pcg_report)
          population_tags <- unique(c("GLOBAL_AF_1KG", pcg_report[["pcgr_config"]][["popgen"]][["vcftag_tgp"]],
                                      "GLOBAL_AF_GNOMAD", pcg_report[["pcgr_config"]][["popgen"]][["vcftag_gnomad"]]))
          for (class in c("tier1", "tier2", "tier3A", "tier3B","gwas")){
            if (class != "tier3B" & class != "gwas"){
              for (c in c("noncancer_phenotype", "cancer_phenotype")){
                if (nrow(pcg_report[["snv_indel"]][["variant_display"]][[class]][[c]]) > 0){
                  pcg_report[["snv_indel"]][["variant_display"]][[class]][[c]] <- pcg_report[["snv_indel"]][["variant_display"]][[class]][[c]] %>%
                    dplyr::mutate(CLINVAR = paste0("<a href=\"http://www.ncbi.nlm.nih.gov/clinvar/variation/", CLINVAR_MSID, "\" target=\"_blank\">", CLINVAR_MSID, "</a>")) %>%
                    dplyr::select(-CLINVAR_MSID)
                  for (tag in population_tags){
                    if (nrow(pcg_report[["snv_indel"]][["variant_display"]][[class]][[c]][is.na(pcg_report[["snv_indel"]][["variant_display"]][[class]][[c]][, tag]), ]) > 0){
                      pcg_report[["snv_indel"]][["variant_display"]][[class]][[c]][is.na(pcg_report[["snv_indel"]][["variant_display"]][[class]][[c]][, tag]), ][, tag] <- 0.00
                    }
                  }
                  if (class == "tier1" | class == "tier2"){
                    pcg_report[["snv_indel"]][["variant_display"]][[class]][[c]] <- dplyr::select(pcg_report[["snv_indel"]][["variant_display"]][[class]][[c]], dplyr::one_of(predispose_tier1_2_display))
                  }
                  else{
                    pcg_report[["snv_indel"]][["variant_display"]][[class]][[c]] <- dplyr::select(pcg_report[["snv_indel"]][["variant_display"]][[class]][[c]], dplyr::one_of(predispose_tier3A_display))
                  }
                }
              }
            }else{
              if (nrow(pcg_report[["snv_indel"]][["variant_display"]][[class]]) > 0){
                for (tag in population_tags){
                  if (nrow(pcg_report[["snv_indel"]][["variant_display"]][[class]][is.na(pcg_report[["snv_indel"]][["variant_display"]][[class]][, tag]), ]) > 0){
                    pcg_report[["snv_indel"]][["variant_display"]][[class]][is.na(pcg_report[["snv_indel"]][["variant_display"]][[class]][, tag]), ][, tag] <- 0.00
                  }
                }
                if(class == "tier3B"){
                  pcg_report[["snv_indel"]][["variant_display"]][[class]] <- dplyr::select(pcg_report[["snv_indel"]][["variant_display"]][[class]], dplyr::one_of(predispose_tier3B_display))
                }
                else{
                  pcg_report[["snv_indel"]][["variant_display"]][[class]] <- dplyr::select(pcg_report[["snv_indel"]][["variant_display"]][[class]], dplyr::one_of(predispose_gwas_display))

                }
              }
            }
          }
        }
      }
    }

    fname_key <- "tsv"
    if (!is.null(pcg_report[["snv_indel"]][["variant_set"]][[fname_key]])){
      if (nrow(pcg_report[["snv_indel"]][["variant_set"]][[fname_key]]) > 0){
        write.table(pcg_report[["snv_indel"]][["variant_set"]][[fname_key]], file = fnames[[fname_key]], sep = "\t", col.names = T, row.names = F, quote = F)
      }
    }

    pcgr_json <- jsonlite::toJSON(pcg_report, pretty = T, na = "string", null = "null")
    write(pcgr_json, fnames[["json"]])
    gzip_command <- paste0("gzip -f ", fnames[["json"]])
    system(gzip_command, intern = F)
    rmarkdown::render(system.file("templates", "report_predisposition.Rmd", package = "pcgrr"), output_format = rmarkdown::html_document(theme = pcg_report[["pcgr_config"]][["visual"]][["report_theme"]], toc = T, toc_depth = 3, toc_float = T, number_sections = F, includes = rmarkdown::includes(after_body = "disclaimer_predisposition.md")), output_file = paste0(sample_name, ".cpsr.", genome_assembly, ".html"), output_dir = project_directory, clean = T, intermediates_dir = project_directory, quiet = T)
  }
}

#' Function that generates summary findings for CPSR
#'
#' @param pcg_report report object
#'
summary_findings_cpsr <- function(pcg_report){

  if(nrow(pcg_report[["snv_indel"]][["variant_display"]][["tier1"]][["cancer_phenotype"]]) > 0){
    pcg_report[["summary"]][["tier1"]] <- as.data.frame(dplyr::select(pcg_report[["snv_indel"]][["variant_display"]][["tier1"]][["cancer_phenotype"]], CLINVAR_MSID, CLINVAR_PHENOTYPE, SYMBOL, HGVSp_short, CONSEQUENCE, ONCOSCORE) %>%
      dplyr::mutate(MUTATION = paste0("<a href=\"http://www.ncbi.nlm.nih.gov/clinvar/variation/",CLINVAR_MSID, "\" target=\"_blank\">", SYMBOL,":",CONSEQUENCE,":",HGVSp_short,"</a>")) %>%
      dplyr::mutate(MUTATION = stringr::str_replace(MUTATION,":NA$","")) %>%
      dplyr::select(-c(HGVSp_short,SYMBOL,CONSEQUENCE, CLINVAR_MSID)) %>%
      tidyr::separate_rows(CLINVAR_PHENOTYPE,sep="; ") %>%
      dplyr::filter(CLINVAR_PHENOTYPE != 'not provided' & CLINVAR_PHENOTYPE != 'not specified') %>%
      dplyr::group_by(CLINVAR_PHENOTYPE) %>%
      dplyr::summarise(VARIANTS = paste(unique(MUTATION), collapse=", "), n = n()) %>%
      dplyr::arrange(desc(n))) %>%
      head(2)
  }

  if(nrow(pcg_report[["snv_indel"]][["variant_display"]][["tier2"]][["cancer_phenotype"]]) > 0){
    pcg_report[["summary"]][["tier2"]] <- as.data.frame(dplyr::select(pcg_report[["snv_indel"]][["variant_display"]][["tier2"]][["cancer_phenotype"]], CLINVAR_MSID, CLINVAR_PHENOTYPE, SYMBOL, HGVSp_short, CONSEQUENCE, ONCOSCORE) %>%
      dplyr::mutate(MUTATION = paste0("<a href=\"http://www.ncbi.nlm.nih.gov/clinvar/variation/",CLINVAR_MSID, "\" target=\"_blank\">", SYMBOL,":",CONSEQUENCE,":",HGVSp_short,"</a>")) %>%
      dplyr::mutate(MUTATION = stringr::str_replace(MUTATION,":NA$","")) %>%
      dplyr::select(-c(HGVSp_short,SYMBOL,CONSEQUENCE,CLINVAR_MSID)) %>%
      tidyr::separate_rows(CLINVAR_PHENOTYPE,sep="; ") %>%
      dplyr::filter(CLINVAR_PHENOTYPE != 'not provided' & CLINVAR_PHENOTYPE != 'not specified') %>%
      dplyr::group_by(CLINVAR_PHENOTYPE) %>%
      dplyr::summarise(VARIANTS = paste(unique(MUTATION), collapse=", "), n = n()) %>%
      dplyr::arrange(desc(n))) %>%
      head(2)
  }

  if(nrow(pcg_report[["snv_indel"]][["variant_display"]][["tier3A"]][["cancer_phenotype"]]) > 0){
    pcg_report[["summary"]][["tier3A"]] <- dplyr::select(pcg_report[["snv_indel"]][["variant_display"]][["tier3A"]][["cancer_phenotype"]], CLINVAR_MSID, CLINVAR_PHENOTYPE, SYMBOL, PATHRANK, HGVSp_short, CONSEQUENCE, ONCOSCORE) %>%
      dplyr::mutate(MUTATION = paste0("<a href=\"http://www.ncbi.nlm.nih.gov/clinvar/variation/",CLINVAR_MSID, "\" target=\"_blank\">", SYMBOL,":",CONSEQUENCE,":",HGVSp_short,"</a>")) %>%
      dplyr::mutate(MUTATION = stringr::str_replace(MUTATION,":NA$","")) %>%
      dplyr::filter(PATHRANK == 'HIGH')

      if(nrow(pcg_report[["summary"]][["tier3A"]]) > 0){
        pcg_report[["summary"]][["tier3A"]] <- pcg_report[["summary"]][["tier3A"]] %>%
          dplyr::select(-c(HGVSp_short,SYMBOL,CONSEQUENCE, CLINVAR_MSID)) %>%
          tidyr::separate_rows(CLINVAR_PHENOTYPE,sep="; ") %>%
          dplyr::filter(CLINVAR_PHENOTYPE != 'not provided' & CLINVAR_PHENOTYPE != 'not specified')
      }

      if(nrow(pcg_report[["summary"]][["tier3A"]]) > 0){
        pcg_report[["summary"]][["tier3A"]] <- pcg_report[["summary"]][["tier3A"]] %>%
          dplyr::group_by(CLINVAR_PHENOTYPE) %>%
          dplyr::summarise(VARIANTS = paste(unique(MUTATION), collapse=", "), n = n()) %>%
          dplyr::arrange(desc(n))
      }
  }

  if(nrow(pcg_report[["snv_indel"]][["variant_display"]][["tier3B"]]) > 0){
    pcg_report[["summary"]][["tier3B"]] <- dplyr::select(pcg_report[["snv_indel"]][["variant_display"]][["tier3B"]],  SYMBOL, PATHRANK, HGVSp_short, CONSEQUENCE, ONCOSCORE) %>%
        dplyr::mutate(MUTATION =paste0(SYMBOL,":",CONSEQUENCE,":",HGVSp_short)) %>%
        dplyr::mutate(MUTATION = stringr::str_replace(MUTATION,":NA$","")) %>%
        dplyr::select(-c(HGVSp_short,SYMBOL,CONSEQUENCE)) %>%
        dplyr::filter(PATHRANK == 'HIGH')
  }

  return(pcg_report)
}


#' Function that counts insilico predictions of variant effects (i.e. damaging/tolerated) from dbNSFP
#'
#' @param sample_calls sample calls with dbnsfp annotations
#'
#' @return sample_calls
#'
get_insilico_prediction_statistics <- function(sample_calls){

  insilico_pathogenicity_pred_algos <- c('SIFT_DBNSFP','PROVEAN_DBNSFP','META_LR_DBNSFP','FATHMM_DBNSFP','MUTATIONTASTER_DBNSFP',
                                         'MUTATIONASSESSOR_DBNSFP','FATHMM_MKL_DBNSF','M_CAP_DBNSFP','SPLICE_SITE_ADA_DBNSFP','SPLICE_SITE_RF_DBNSFP')
  for(v in c('CALLED','DAMAGING','TOLERATED','SPLICING_NEUTRAL','SPLICING_AFFECTED')){
    sample_calls[,paste0('N_INSILICO_',v)] <- 0
  }

  for(algo in insilico_pathogenicity_pred_algos){
    if(algo %in% colnames(sample_calls)){
      sample_calls[!is.na(sample_calls[,algo]) & sample_calls[,algo] != '.','N_INSILICO_CALLED'] <-  sample_calls[!is.na(sample_calls[,algo]) & sample_calls[,algo] != '.','N_INSILICO_CALLED'] + 1
      sample_calls[!is.na(sample_calls[,algo]) & (sample_calls[,algo] == 'damaging' | sample_calls[,algo] == 'possibly_damaging'),'N_INSILICO_DAMAGING'] <-  sample_calls[!is.na(sample_calls[,algo]) & (sample_calls[,algo] == 'damaging' | sample_calls[,algo] == 'possibly_damaging'),'N_INSILICO_DAMAGING'] + 1
      sample_calls[!is.na(sample_calls[,algo]) & sample_calls[,algo] == 'tolerated','N_INSILICO_TOLERATED'] <-  sample_calls[!is.na(sample_calls[,algo]) & sample_calls[,algo] == 'tolerated','N_INSILICO_TOLERATED'] + 1
      sample_calls[!is.na(sample_calls[,algo]) & sample_calls[,algo] == 'affect_splicing','N_INSILICO_SPLICING_AFFECTED'] <-  sample_calls[!is.na(sample_calls[,algo]) & sample_calls[,algo] == 'affect_splicing','N_INSILICO_SPLICING_AFFECTED'] + 1
      sample_calls[!is.na(sample_calls[,algo]) & sample_calls[,algo] == 'splicing_neutral','N_INSILICO_SPLICING_NEUTRAL'] <-  sample_calls[!is.na(sample_calls[,algo]) & sample_calls[,algo] == 'splicing_neutral','N_INSILICO_SPLICING_NEUTRAL'] + 1
    }
  }
  return(sample_calls)
}

#' Function that assigns variant pathogenicity scores based on ACMG guidelines
#'
#' @param sample_calls sample calls with dbnsfp annotations
#' @param pcgr_config pcgr configuration object
#' @param pcgr_data pcgr data object
#'
#' @return sample_calls
#'
assign_pathogenicity_score <- function(sample_calls, pcgr_config, pcgr_data){

  path_columns <- c("PVS1","PSC1","PS1","PM1", "PM2", "BA1","PP3","PM4","PPC1", "BP4", "PP2","BMC1","BSC1","N_INSILICO_CALLED",
                    "N_INSILICO_DAMAGING", "N_INSILICO_TOLERATED", "N_INSILICO_SPLICING_NEUTRAL",
                    "N_INSILICO_SPLICING_AFFECTED", "codon_prefix", "clinvar_pathogenic_codon",
                    "clinvar_pathogenic", "clinvar_benign","clinvar_benign_codon","hotspot_symbol", "hotspot_codon","hotspot_pvalue",
                    "LOF_KNOWN_MOI","PATH_TRUNCATION_RATE","BENIGN_MISSENSE_RATE",
                    "PATHSCORE","PATHDOC","PATHRANK")
  sample_calls <- sample_calls[, !(colnames(sample_calls) %in% path_columns)]

  sample_calls$PVS1 <- FALSE #PVS1 null variant (nonsense, frameshift, canonical ±1 or 2 splice sites, initiation codon, single or multiexon deletion)
  #in a gene where LOF is a known mechanism of disease (Dominant mode of inheritance)
  sample_calls$PSC1 <- FALSE #null variant (nonsense, frameshift, canonical ±1 or 2 splice sites, initiation codon, single or multiexon deletion)
  #in a gene where LOF is a known mechanism of disease (Recessive mode of inheritance)
  sample_calls$PS1 <- FALSE #Same amino acid change as a previously established pathogenic variant regardless of nucleotide change
  sample_calls$PMC1 <- FALSE #null variant in a gene where LoF is not a known mechanism of disease
  sample_calls$PM2 <- FALSE #Absence or extremely low frequency (MAF < 0.0005) in 1000 Genomes Project, or gnomAD
  sample_calls$PP2 <- FALSE # Missense variant in a gene that has a relatively low rate of benign missense variation (<20%) and where missense variants are a common mechanism of disease (>50% of pathogenic variants)
  sample_calls$PM4 <- FALSE #Protein length changes due to inframe indels or nonstop variant of genes that harbor variants with a dominant mode of inheritance.
  sample_calls$PPC1 <- FALSE #Protein length changes due to inframe indels or nonstop variant of genes that harbor variants with a recessive mode of inheritance.
  sample_calls$PM5 <- FALSE #Novel missense change at an amino acid residue where a different missense change determined to be pathogenic has been seen before
  sample_calls$PP3 <- FALSE #Multiple lines (>=5) of computational evidence support a deleterious effect on the gene or gene product (conservation, evolutionary, splicing impact, etc.)
  sample_calls$BP4 <- FALSE #Multiple lines (>=5) of in silico evidence of none deleterious effect.
  sample_calls$BMC1 <- FALSE #Peptide change is at the same location of a known benign change
  sample_calls$BSC1 <- FALSE #Peptide change is known to be benign
  sample_calls$BA1 <- FALSE #Common population frequency (MAF > 0.05) in 1000 Genomes Project, or gnomAD
  sample_calls$BP1 <- FALSE #Missense variant in a gene for which primarily truncating variants are known to cause disease

  predisposition_gene_info <- dplyr::select(pcgr_data$predisposition_genes, symbol, lof_known_moi, path_truncation_rate, benign_missense_rate) %>%
    dplyr::rename(SYMBOL = symbol, LOF_KNOWN_MOI = lof_known_moi, PATH_TRUNCATION_RATE = path_truncation_rate, BENIGN_MISSENSE_RATE = benign_missense_rate)

  sample_calls <- dplyr::left_join(sample_calls, predisposition_gene_info, by=c("SYMBOL" = "SYMBOL"))

  ## Assign logical ACMG level
  # PP3 - Multiple lines of computational evidence support a deleterious effect on the gene or gene product (conservation, evolutionary, splicing impact, etc.)
  # BP4 - Multiple lines (>1) of in silico evidence of none deleterious effect.
  #
  # Computational evidence for deleterious/benign effect is taken from invidual algorithm predictions in dbNSFP
  # SIFT,Provean,MutationTaster,MutationAssessor,M_CAP,MutPred,FATHMM,FATHMM-mkl,DBNSFP_LogReg,dbscSNV_RF,dbscSNV_AdaBoost
  # Default scheme (from default TOML file):
  # 1) Damaging: Among 8 possible protein variant effect predictions, at least six algorithms must have made a call, with at least 5 predicted as damaging,
  #    and at most one predicted as tolerated (PP3)
  #    - at most 1 prediction for a splicing neutral effect
  # 2) Tolerated: Among 8 possible protein variant effect predictions, at least six algorithms must have made a call, with at least 5 predicted as tolerated,
  #    and at most one predicted as damaging (BP4)
  #    - 0 predictions of splice site affected

  dbnsfp_min_majority <- pcgr_config[['dbnsfp']][['min_majority']]
  dbnsfp_max_minority <- pcgr_config[['dbnsfp']][['max_minority']]
  dbnsfp_min_called <- dbnsfp_min_majority

  sample_calls <- pcgrr::get_insilico_prediction_statistics(sample_calls)
  if(nrow(sample_calls[sample_calls$N_INSILICO_CALLED >= dbnsfp_min_called &
                       sample_calls$N_INSILICO_DAMAGING >= dbnsfp_min_majority &
                       sample_calls$N_INSILICO_TOLERATED <= dbnsfp_max_minority &
                       sample_calls$N_INSILICO_SPLICING_NEUTRAL <= 1,]) > 0){
    sample_calls[sample_calls$N_INSILICO_CALLED >= dbnsfp_min_called &
                   sample_calls$N_INSILICO_DAMAGING >= dbnsfp_min_majority &
                   sample_calls$N_INSILICO_TOLERATED <= dbnsfp_max_minority &
                   sample_calls$N_INSILICO_SPLICING_NEUTRAL <= 1,]$PP3 <- TRUE
  }
  if(nrow(sample_calls[sample_calls$N_INSILICO_CALLED >= dbnsfp_min_called &
                       sample_calls$N_INSILICO_TOLERATED >= dbnsfp_min_majority &
                       sample_calls$N_INSILICO_DAMAGING <= dbnsfp_max_minority &
                       sample_calls$N_INSILICO_SPLICING_AFFECTED == 0,]) > 0){
    sample_calls[sample_calls$N_INSILICO_CALLED >= dbnsfp_min_called &
                   sample_calls$N_INSILICO_TOLERATED >= dbnsfp_min_majority &
                   sample_calls$N_INSILICO_DAMAGING <= dbnsfp_max_minority &
                   sample_calls$N_INSILICO_SPLICING_AFFECTED == 0,]$BP4 <- TRUE
  }
  if(nrow(sample_calls[sample_calls$N_INSILICO_SPLICING_AFFECTED == 2,]) > 0){
    sample_calls[sample_calls$N_INSILICO_SPLICING_AFFECTED == 2,]$PP3 <- TRUE
  }

  ## Assign logical ACMG level
  # PM2 -  absence/extremely low germline population frequency (1000G/gnomAD)
  if('GLOBAL_AF_1KG' %in% colnames(sample_calls) & 'GLOBAL_AF_GNOMAD' %in% colnames(sample_calls)){
    if(nrow(sample_calls[is.na(sample_calls$GLOBAL_AF_1KG) & is.na(sample_calls$GLOBAL_AF_GNOMAD),]) > 0){
      sample_calls[is.na(sample_calls$GLOBAL_AF_1KG) & is.na(sample_calls$GLOBAL_AF_GNOMAD),]$PM2 <- TRUE
    }
    if(nrow(sample_calls[!is.na(sample_calls$GLOBAL_AF_GNOMAD) & sample_calls$GLOBAL_AF_GNOMAD < 0.0005,]) > 0){
      sample_calls[!is.na(sample_calls$GLOBAL_AF_GNOMAD) & sample_calls$GLOBAL_AF_GNOMAD < 0.0005,]$PM2 <- TRUE
    }

  }

  ## Assign logical ACMG level
  # BA1 -  high population germline frequency (1000G/gnomAD)
  if('GLOBAL_AF_1KG' %in% colnames(sample_calls) & 'GLOBAL_AF_GNOMAD' %in% colnames(sample_calls)){
    if(nrow(sample_calls[!is.na(sample_calls$GLOBAL_AF_GNOMAD) & sample_calls$GLOBAL_AF_GNOMAD > 0.05,]) > 0){
      sample_calls[!is.na(sample_calls$GLOBAL_AF_GNOMAD) & sample_calls$GLOBAL_AF_GNOMAD > 0.05,]$BA1 <- TRUE
    }
    if(nrow(sample_calls[!is.na(sample_calls$GLOBAL_AF_1KG) & sample_calls$GLOBAL_AF_1KG > 0.05,]) > 0){
      sample_calls[!is.na(sample_calls$GLOBAL_AF_1KG) & sample_calls$GLOBAL_AF_1KG > 0.05,]$BA1 <- TRUE
    }
    if(nrow(sample_calls[(sample_calls$SYMBOL == 'HFE' | sample_calls$SYMBOL == 'SERPINA1') & sample_calls$BA1 == TRUE & !is.na(sample_calls$GLOBAL_AF_1KG) & sample_calls$GLOBAL_AF_1KG < 0.25,]) > 0){
      sample_calls[(sample_calls$SYMBOL == 'HFE' | sample_calls$SYMBOL == 'SERPINA1') & sample_calls$BA1 == TRUE & !is.na(sample_calls$GLOBAL_AF_1KG) & sample_calls$GLOBAL_AF_1KG < 0.25,]$BA1 <- FALSE
    }
    if(nrow(sample_calls[(sample_calls$SYMBOL == 'HFE' | sample_calls$SYMBOL == 'SERPINA1') & sample_calls$BA1 == TRUE & !is.na(sample_calls$GLOBAL_AF_GNOMAD) & sample_calls$GLOBAL_AF_GNOMAD < 0.25,]) > 0){
      sample_calls[(sample_calls$SYMBOL == 'HFE' | sample_calls$SYMBOL == 'SERPINA1') & sample_calls$BA1 == TRUE & !is.na(sample_calls$GLOBAL_AF_GNOMAD) & sample_calls$GLOBAL_AF_GNOMAD < 0.25,]$BA1 <- FALSE
    }

  }

  ## Assign logical ACMG level on loss-of-function variant in known predisposition gene
  # PVS1 - Truncations in susceptibility genes where LOF is a known mechanism of the disease and harbor variants with a dominant mode of inheritance.
  # PSC1 - Truncations in susceptibility genes where LOF is a known mechanism of the disease and harbor variants with a recessive mode of inheritance.
  if(nrow(sample_calls[!is.na(sample_calls$LOSS_OF_FUNCTION) & sample_calls$LOSS_OF_FUNCTION == TRUE & !is.na(sample_calls$LOF_KNOWN_MOI) & sample_calls$LOF_KNOWN_MOI == "LoF" & !is.na(sample_calls$CANCER_PREDISPOSITION_MOI),]) > 0){
    if(nrow(sample_calls[!is.na(sample_calls$LOSS_OF_FUNCTION) & sample_calls$LOSS_OF_FUNCTION == TRUE & !is.na(sample_calls$LOF_KNOWN_MOI) & sample_calls$LOF_KNOWN_MOI == "LoF" & !is.na(sample_calls$CANCER_PREDISPOSITION_MOI) & stringr::str_detect(sample_calls$CANCER_PREDISPOSITION_MOI,"AD|AD/AR"),]) > 0){
      sample_calls[!is.na(sample_calls$LOSS_OF_FUNCTION) & sample_calls$LOSS_OF_FUNCTION == TRUE & !is.na(sample_calls$LOF_KNOWN_MOI) & sample_calls$LOF_KNOWN_MOI == "LoF" & !is.na(sample_calls$CANCER_PREDISPOSITION_MOI) & stringr::str_detect(sample_calls$CANCER_PREDISPOSITION_MOI,"AD|AD/AR"),]$PVS1 <- TRUE
    }

    if(nrow(sample_calls[!is.na(sample_calls$LOSS_OF_FUNCTION) & sample_calls$LOSS_OF_FUNCTION == TRUE & !is.na(sample_calls$LOF_KNOWN_MOI) & sample_calls$LOF_KNOWN_MOI == "LoF" & !is.na(sample_calls$CANCER_PREDISPOSITION_MOI) & stringr::str_detect(sample_calls$CANCER_PREDISPOSITION_MOI,"^AR$"),]) > 0){
      sample_calls[!is.na(sample_calls$LOSS_OF_FUNCTION) & sample_calls$LOSS_OF_FUNCTION == TRUE & !is.na(sample_calls$LOF_KNOWN_MOI) & sample_calls$LOF_KNOWN_MOI == "LoF" & !is.na(sample_calls$CANCER_PREDISPOSITION_MOI) & stringr::str_detect(sample_calls$CANCER_PREDISPOSITION_MOI,"^AR$"),]$PSC1 <- TRUE
      ## skip genes that are associated with gain-of-function
    }
  }

  ## Assign a logical ACMG level
  # PM4 - Protein length changes due to inframe indels or nonstop variant of genes that harbor variants with a dominant mode of inheritance - PM4
  # PCC1 - Protein length changes due to inframe indels or nonstop variant of genes that harbor variants with a recessive mode of inheritance (and unknown MOI) - PPC1
  if(nrow(sample_calls[!is.na(sample_calls$CONSEQUENCE) & stringr::str_detect(sample_calls$CONSEQUENCE,'stop_lost|inframe_deletion|inframe_insertion'),]) > 0){
    if(nrow(sample_calls[!is.na(sample_calls$CONSEQUENCE) & stringr::str_detect(sample_calls$CONSEQUENCE,'stop_lost|inframe_deletion|inframe_insertion') & !is.na(sample_calls$CANCER_PREDISPOSITION_MOI) & stringr::str_detect(sample_calls$CANCER_PREDISPOSITION_MOI,"AD/AR|AD"),]) > 0){
      sample_calls[!is.na(sample_calls$CONSEQUENCE) & stringr::str_detect(sample_calls$CONSEQUENCE,'stop_lost|inframe_deletion|inframe_insertion') & !is.na(sample_calls$CANCER_PREDISPOSITION_MOI) & stringr::str_detect(sample_calls$CANCER_PREDISPOSITION_MOI,"AD/AR|AD"),]$PM4 <- TRUE
    }
    if(nrow(sample_calls[!is.na(sample_calls$CONSEQUENCE) & stringr::str_detect(sample_calls$CONSEQUENCE,'stop_lost|inframe_deletion|inframe_insertion') & !is.na(sample_calls$CANCER_PREDISPOSITION_MOI) & !stringr::str_detect(sample_calls$CANCER_PREDISPOSITION_MOI,"AD|Mosaic"),]) > 0){
      sample_calls[!is.na(sample_calls$CONSEQUENCE) & stringr::str_detect(sample_calls$CONSEQUENCE,'stop_lost|inframe_deletion|inframe_insertion') & !is.na(sample_calls$CANCER_PREDISPOSITION_MOI) & !stringr::str_detect(sample_calls$CANCER_PREDISPOSITION_MOI,"AD|Mosaic"),]$PPC1 <- TRUE
    }
  }

  ## Assign a logical ACMG level
  # PP2 - Missense variant in a gene that has a relatively low rate of benign missense variation and where missense variants are a common mechanism of disease
  if(nrow(sample_calls[!is.na(sample_calls$CONSEQUENCE) & !is.na(sample_calls$BENIGN_MISSENSE_RATE) & sample_calls$BENIGN_MISSENSE_RATE <= 0.2 & sample_calls$PATH_TRUNCATION_RATE < 0.5 & stringr::str_detect(sample_calls$CONSEQUENCE,'^missense_variant'),]) > 0){
    sample_calls[!is.na(sample_calls$CONSEQUENCE) & !is.na(sample_calls$BENIGN_MISSENSE_RATE) & sample_calls$BENIGN_MISSENSE_RATE <= 0.2 & sample_calls$PATH_TRUNCATION_RATE < 0.5 & stringr::str_detect(sample_calls$CONSEQUENCE,'^missense_variant'),]$PP2 <- TRUE
  }

  ## Assign a logical ACMG level
  # BP1 - Missense variant in a gene for which primarily truncating variants are known to cause disease
  if(nrow(sample_calls[!is.na(sample_calls$CONSEQUENCE) & !is.na(sample_calls$PATH_TRUNCATION_RATE) & sample_calls$PATH_TRUNCATION_RATE > 0.90 & stringr::str_detect(sample_calls$CONSEQUENCE,'^missense_variant'),]) > 0){
    sample_calls[!is.na(sample_calls$CONSEQUENCE) & !is.na(sample_calls$PATH_TRUNCATION_RATE) & sample_calls$PATH_TRUNCATION_RATE > 0.90 & stringr::str_detect(sample_calls$CONSEQUENCE,'^missense_variant'),]$BP1 <- TRUE
  }

  ## Assign a logical ACMG level
  # PMC1 - LoF variant in susceptibility gene, yet unknown or not LoF mechanism of disease
  if(nrow(sample_calls[!is.na(sample_calls$LOSS_OF_FUNCTION) & sample_calls$LOSS_OF_FUNCTION == TRUE & (is.na(sample_calls$LOF_KNOWN_MOI) | sample_calls$LOF_KNOWN_MOI != "LoF"),]) > 0){
    sample_calls[!is.na(sample_calls$LOSS_OF_FUNCTION) & sample_calls$LOSS_OF_FUNCTION == TRUE & (is.na(sample_calls$LOF_KNOWN_MOI) | sample_calls$LOF_KNOWN_MOI != "LoF"),]$PMC1 <- TRUE
  }

  ## Assign logical ACMG level
  # PS1 - coinciding with known pathogenic missense variants (yet with different nucleotide change)
  # PM5 - occurs at the same codon as a known pathogenic missense variant
  sample_calls$codon_prefix <- NA
  if(nrow(sample_calls[!is.na(sample_calls$CONSEQUENCE) & sample_calls$CONSEQUENCE == 'missense_variant',]) > 0){
    sample_calls[sample_calls$CONSEQUENCE == 'missense_variant',]$codon_prefix <- stringr::str_match(sample_calls[sample_calls$CONSEQUENCE == 'missense_variant',]$HGVSp_short,"p\\.[A-Z]{1}[0-9]{1,}")
  }

  if(nrow(sample_calls[!is.na(sample_calls$codon_prefix),]) > 0){
    sample_calls_pathogenic_codon <- dplyr::left_join(dplyr::filter(dplyr::select(sample_calls,VAR_ID,codon_prefix,SYMBOL), !is.na(codon_prefix)), pcgr_data$clinvar_cpg_loci[['pathogenic']][['codon']], by=c("codon_prefix" = "codon_prefix","SYMBOL" = "symbol"))
    sample_calls <- dplyr::left_join(sample_calls, dplyr::select(sample_calls_pathogenic_codon, VAR_ID, clinvar_pathogenic_codon), by=c("VAR_ID"))

    sample_calls_benign_codon <- dplyr::left_join(dplyr::filter(dplyr::select(sample_calls,VAR_ID,codon_prefix,SYMBOL), !is.na(codon_prefix)), pcgr_data$clinvar_cpg_loci[['benign']][['codon']], by=c("codon_prefix" = "codon_prefix","SYMBOL" = "symbol"))
    sample_calls <- dplyr::left_join(sample_calls, dplyr::select(sample_calls_benign_codon, VAR_ID, clinvar_benign_codon), by=c("VAR_ID"))
  }
  if(nrow(sample_calls[!is.na(sample_calls$HGVSp_short),]) > 0){
    sample_calls_pathogenic_hgvsp <- dplyr::left_join(dplyr::filter(dplyr::select(sample_calls, VAR_ID, HGVSp_short, SYMBOL), !is.na(HGVSp_short)), pcgr_data$clinvar_cpg_loci[['pathogenic']][['peptide_change']], by=c("HGVSp_short" = "hgvs_p","SYMBOL" = "symbol"))
    sample_calls <- dplyr::left_join(sample_calls, dplyr::select(sample_calls_pathogenic_hgvsp, VAR_ID, clinvar_pathogenic), by=c("VAR_ID"))

    sample_calls_benign_hgvsp <- dplyr::left_join(dplyr::filter(dplyr::select(sample_calls, VAR_ID, HGVSp_short, SYMBOL), !is.na(HGVSp_short)), pcgr_data$clinvar_cpg_loci[['benign']][['peptide_change']], by=c("HGVSp_short" = "hgvs_p","SYMBOL" = "symbol"))
    sample_calls <- dplyr::left_join(sample_calls, dplyr::select(sample_calls_benign_hgvsp, VAR_ID, clinvar_benign), by=c("VAR_ID"))

  }

  if(nrow(sample_calls[!is.na(sample_calls$clinvar_pathogenic) & sample_calls$clinvar_pathogenic == T,]) > 0){
    sample_calls[!is.na(sample_calls$clinvar_pathogenic) & sample_calls$clinvar_pathogenic == T,]$PS1 <- TRUE
  }
  if(nrow(sample_calls[!is.na(sample_calls$clinvar_pathogenic_codon) & sample_calls$clinvar_pathogenic_codon == T,]) > 0){
    sample_calls[!is.na(sample_calls$clinvar_pathogenic_codon) & sample_calls$clinvar_pathogenic_codon == T,]$PM5 <- TRUE
  }

  ## if previously found coinciding with pathogenic variant, set PM5 to false
  if(nrow(sample_calls[is.na(sample_calls$PS1) & is.na(sample_calls$PM5) & sample_calls$PM5 == T & sample_calls$PS1 == T,]) > 0){
    sample_calls[is.na(sample_calls$PS1) & is.na(sample_calls$PM5) & sample_calls$PM5 == T & sample_calls$PS1 == T,]$PM5 <- FALSE
  }

  if(nrow(sample_calls[!is.na(sample_calls$clinvar_benign) & sample_calls$clinvar_benign == T,]) > 0){
    sample_calls[!is.na(sample_calls$clinvar_benign) & sample_calls$clinvar_benign == T,]$BSC1 <- TRUE
  }
  if(nrow(sample_calls[!is.na(sample_calls$clinvar_benign_codon) & sample_calls$clinvar_benign_codon == T,]) > 0){
    sample_calls[!is.na(sample_calls$clinvar_benign_codon) & sample_calls$clinvar_benign_codon == T,]$BMC1 <- TRUE
  }

  ## if previously found coinciding with benign variant, set BMC1 to false
  if(nrow(sample_calls[is.na(sample_calls$BSC1) & is.na(sample_calls$BMC1) & sample_calls$BSC1 == T & sample_calls$BMC1 == T,]) > 0){
    sample_calls[is.na(sample_calls$BSC1) & is.na(sample_calls$BMC1) & sample_calls$BSC1 == T & sample_calls$BMC1 == T,]$BMC1 <- FALSE
  }

  ##Assign logical ACMG level
  # PM1 - missense variant in a somatic mutation hotspot as determined by cancerhotspots.org (v2)
  sample_calls <- sample_calls %>% tidyr::separate(MUTATION_HOTSPOT, c("hotspot_symbol","hotspot_codon","hotspot_pvalue"),sep="\\|",remove=F,extra="drop")
  if(nrow(sample_calls[!is.na(sample_calls$hotspot_codon),]) > 0){
    sample_calls[!is.na(sample_calls$hotspot_codon),]$hotspot_codon <- paste0('p.',sample_calls[!is.na(sample_calls$hotspot_codon),]$hotspot_codon)
  }

  sample_calls_somatic_hotspots <- sample_calls %>%
    dplyr::filter(!is.na(hotspot_codon) & !is.na(hotspot_symbol)) %>%
    dplyr::filter(!is.na(codon_prefix) & !is.na(SYMBOL)) %>%
    dplyr::filter(SYMBOL == hotspot_symbol & hotspot_codon == codon_prefix) %>%
    dplyr::select(VAR_ID) %>%
    dplyr::distinct() %>%
    dplyr::mutate(PM1 = TRUE)

  if(nrow(sample_calls_somatic_hotspots) > 0){
    sample_calls <- dplyr::left_join(sample_calls, sample_calls_somatic_hotspots, by=c("VAR_ID"))
  }else{
    sample_calls$PM1 <- FALSE
  }
  if(nrow(sample_calls[is.na(sample_calls$PM1),]) > 0){
    sample_calls[is.na(sample_calls$PM1),]$PM1 <- FALSE
  }

  sample_calls$PATHSCORE <- 0
  sample_calls$PATHDOC <- ""
  if(nrow(sample_calls[sample_calls$PVS1 == TRUE,]) > 0){
    sample_calls[sample_calls$PVS1 == TRUE,]$PATHSCORE <- sample_calls[sample_calls$PVS1 == TRUE,]$PATHSCORE + 8
    sample_calls[sample_calls$PVS1 == TRUE,]$PATHDOC <- paste(sample_calls[sample_calls$PVS1 == TRUE,]$PATHDOC,"- Loss-of-function variant in known susceptibility/syndrome gene (dominant MoI) - PVS1",sep="<br>")
  }
  if(nrow(sample_calls[sample_calls$PS1 == TRUE,]) > 0){
    sample_calls[sample_calls$PS1 == TRUE,]$PATHSCORE <- sample_calls[sample_calls$PS1 == TRUE,]$PATHSCORE + 7
    sample_calls[sample_calls$PS1 == TRUE,]$PATHDOC <- paste(sample_calls[sample_calls$PS1 == TRUE,]$PATHDOC,"- Same peptide change as a previously established pathogenic variant (ClinVar) - PS1",sep="<br>")
  }
  if(nrow(sample_calls[sample_calls$PSC1 == TRUE,]) > 0){
    sample_calls[sample_calls$PSC1 == TRUE,]$PATHSCORE <- sample_calls[sample_calls$PSC1 == TRUE,]$PATHSCORE + 4
    sample_calls[sample_calls$PSC1 == TRUE,]$PATHDOC <- paste(sample_calls[sample_calls$PSC1 == TRUE,]$PATHDOC,"- Loss-of-function variant in known susceptibility/syndrome gene (recessive MoI) - PSC1",sep="<br>")
  }
  if(nrow(sample_calls[sample_calls$PMC1 == TRUE,]) > 0){
    sample_calls[sample_calls$PMC1 == TRUE,]$PATHSCORE <- sample_calls[sample_calls$PMC1 == TRUE,]$PATHSCORE + 2
    sample_calls[sample_calls$PMC1 == TRUE,]$PATHDOC <- paste(sample_calls[sample_calls$PMC1 == TRUE,]$PATHDOC,"- Loss-of-function variant in known susceptibility/syndrome gene (unknown/LoF not a known mechanism of disease) - PMC1",sep="<br>")
  }

  if(nrow(sample_calls[sample_calls$PM4 == TRUE,]) > 0){
    sample_calls[sample_calls$PM4 == TRUE,]$PATHSCORE <- sample_calls[sample_calls$PM4 == TRUE,]$PATHSCORE + 2
    sample_calls[sample_calls$PM4 == TRUE,]$PATHDOC <- paste(sample_calls[sample_calls$PM4 == TRUE,]$PATHDOC,"- Protein length changes due to inframe indels/stoploss variants in known susceptibility genes (dominant MoI) - PM4",sep="<br>")
  }

  if(nrow(sample_calls[sample_calls$PM1 == TRUE,]) > 0){
    sample_calls[sample_calls$PM1 == TRUE,]$PATHSCORE <- sample_calls[sample_calls$PM1 == TRUE,]$PATHSCORE + 2
    sample_calls[sample_calls$PM1 == TRUE,]$PATHDOC <- paste(sample_calls[sample_calls$PM1 == TRUE,]$PATHDOC,"- Variant located in somatic mutation hotspot (cancerhotspots.org) - PM1",sep="<br>")
  }
  if(nrow(sample_calls[sample_calls$PM2 == TRUE,]) > 0){
    sample_calls[sample_calls$PM2 == TRUE,]$PATHSCORE <- sample_calls[sample_calls$PM2 == TRUE,]$PATHSCORE + 2
    sample_calls[sample_calls$PM2 == TRUE,]$PATHDOC <- paste(sample_calls[sample_calls$PM2 == TRUE,]$PATHDOC,"- Absent or extremely low allele frequency in the general population (gnomAD/1KG global MAF < 0.0005) - PM2",sep="<br>")
  }
  if(nrow(sample_calls[sample_calls$PM5 == TRUE,]) > 0){
    sample_calls[sample_calls$PM5 == TRUE,]$PATHSCORE <- sample_calls[sample_calls$PM5 == TRUE,]$PATHSCORE + 2
    sample_calls[sample_calls$PM5 == TRUE,]$PATHDOC <- paste(sample_calls[sample_calls$PM5 == TRUE,]$PATHDOC,"- Different peptide change of a pathogenic variant at the same reference peptide (ClinVar) - PM5",sep="<br>")
  }
  if(nrow(sample_calls[sample_calls$PP2 == TRUE,]) > 0){
    sample_calls[sample_calls$PP2 == TRUE,]$PATHSCORE <- sample_calls[sample_calls$PP2 == TRUE,]$PATHSCORE + 1
    sample_calls[sample_calls$PP2 == TRUE,]$PATHDOC <- paste(sample_calls[sample_calls$PP2 == TRUE,]$PATHDOC,"- Missense variant in susceptibility genes where pathogenic missense variants are relatively common - PP2",sep="<br>")
  }

  if(nrow(sample_calls[sample_calls$PPC1 == TRUE,]) > 0){
    sample_calls[sample_calls$PPC1 == TRUE,]$PATHSCORE <- sample_calls[sample_calls$PPC1 == TRUE,]$PATHSCORE + 1
    sample_calls[sample_calls$PPC1 == TRUE,]$PATHDOC <- paste(sample_calls[sample_calls$PPC1 == TRUE,]$PATHDOC,"- Protein length changes due to inframe indels/stoploss variants in known susceptibility genes (recessive MoI) - PPC1",sep="<br>")
  }
  if(nrow(sample_calls[sample_calls$PP3 == TRUE,]) > 0){
    sample_calls[sample_calls$PP3 == TRUE,]$PATHSCORE <- sample_calls[sample_calls$PP3 == TRUE,]$PATHSCORE + 1
    sample_calls[sample_calls$PP3 == TRUE,]$PATHDOC <- paste(sample_calls[sample_calls$PP3 == TRUE,]$PATHDOC,paste0("- Multiple lines (>=",pcgr_config[['dbnsfp']][['min_majority']],") of in silico evidence of deleterious effect - PP3"),sep="<br>")
  }

  if(nrow(sample_calls[sample_calls$BP1 == TRUE,]) > 0){
    sample_calls[sample_calls$BP1 == TRUE,]$PATHSCORE <- sample_calls[sample_calls$BP1 == TRUE,]$PATHSCORE - 1
    sample_calls[sample_calls$BP1 == TRUE,]$PATHDOC <- paste(sample_calls[sample_calls$BP1 == TRUE,]$PATHDOC,paste0("- Missense variant in a gene for which primarily (>90%) truncating variants are known to cause disease - BP1"),sep="<br>")
  }

  if(nrow(sample_calls[sample_calls$BP4 == TRUE,]) > 0){
    sample_calls[sample_calls$BP4 == TRUE,]$PATHSCORE <- sample_calls[sample_calls$BP4 == TRUE,]$PATHSCORE - 1
    sample_calls[sample_calls$BP4 == TRUE,]$PATHDOC <- paste(sample_calls[sample_calls$BP4 == TRUE,]$PATHDOC,paste0("- Multiple lines (>=",pcgr_config[['dbnsfp']][['min_majority']],") of in silico evidence of benign effect - BP4"),sep="<br>")
  }

  if(nrow(sample_calls[sample_calls$BMC1 == TRUE,]) > 0){
    sample_calls[sample_calls$BMC1 == TRUE,]$PATHSCORE <- sample_calls[sample_calls$BMC1 == TRUE,]$PATHSCORE - 2
    sample_calls[sample_calls$BMC1 == TRUE,]$PATHDOC <- paste(sample_calls[sample_calls$BMC1 == TRUE,]$PATHDOC,paste0("- Peptide change is at the same location of a known benign change (ClinVar) - BMC1"),sep="<br>")
  }

  if(nrow(sample_calls[sample_calls$BSC1 == TRUE,]) > 0){
    sample_calls[sample_calls$BSC1 == TRUE,]$PATHSCORE <- sample_calls[sample_calls$BSC1 == TRUE,]$PATHSCORE - 6
    sample_calls[sample_calls$BSC1 == TRUE,]$PATHDOC <- paste(sample_calls[sample_calls$BSC1 == TRUE,]$PATHDOC,paste0("- Peptide change is known to be benign (ClinVar) - BSC1"),sep="<br>")
  }

  if(nrow(sample_calls[sample_calls$BA1 == TRUE,]) > 0){
    sample_calls[sample_calls$BA1 == TRUE,]$PATHSCORE <- sample_calls[sample_calls$BA1 == TRUE,]$PATHSCORE - 8
    sample_calls[sample_calls$BA1 == TRUE,]$PATHDOC <- paste(sample_calls[sample_calls$BA1 == TRUE,]$PATHDOC,paste0("- High allele frequency in the general population (gnomAD/1KG global MAF > 0.05 (exceptions for HFE/SERPINA1, MAF > 0.25)) - BA1"),sep="<br>")
  }

  sample_calls$PATHDOC <- stringr::str_replace(sample_calls$PATHDOC,"^<br>","")
  sample_calls$PATHDOC <- paste0("CPSR pathogenicity score: ",sample_calls$PATHSCORE, "<br>",sample_calls$PATHDOC)
  sample_calls$PATHDOC <- stringr::str_replace(sample_calls$PATHDOC,"<br>$","")

  sample_calls$PATHRANK <- NA
  if(nrow(sample_calls[sample_calls$PATHSCORE > 8,]) > 0){
    sample_calls[sample_calls$PATHSCORE > 8,]$PATHRANK <- 'HIGH'
  }
  if(nrow(sample_calls[sample_calls$PATHSCORE > 4 & sample_calls$PATHSCORE <= 8,]) > 0){
    sample_calls[sample_calls$PATHSCORE > 4 & sample_calls$PATHSCORE <= 8,]$PATHRANK <- 'MODERATE'
  }
  if(nrow(sample_calls[sample_calls$PATHSCORE <= 4,]) > 0){
    sample_calls[sample_calls$PATHSCORE <= 4,]$PATHRANK <- 'LOW'
  }
  if(nrow(sample_calls[sample_calls$PATHSCORE < 0,]) > 0){
    sample_calls[sample_calls$PATHSCORE < 0,]$PATHRANK <- 'BENIGN'
  }
  return(sample_calls)

}



#' Function that generates tiered annotated variant datasets for CPSR
#'
#' @param pcg_report List with tiered variants
#'
#' @return tsv_variants data frame with tier-annotated list of variants for tab-separated output
#'
generate_tier_tsv_cpsr <- function(pcg_report, sample_name = "test"){

  predispose_tsv_tags <- c("VAR_ID","VCF_SAMPLE_ID","CODING_STATUS","SYMBOL","ENSEMBL_GENE_ID","ENSEMBL_TRANSCRIPT_ID","GENE_NAME",
                           "LOF_KNOWN_MOI","PATH_TRUNCATION_RATE","BENIGN_MISSENSE_RATE","ONCOGENE", "ONCOSCORE","TUMOR_SUPPRESSOR",
                           "GENOTYPE","CONSEQUENCE","PROTEIN_CHANGE","PROTEIN_DOMAIN", "HGVSp","HGVSc", "CDS_CHANGE","MUTATION_HOTSPOT","RMSK_HIT","PROTEIN_FEATURE","EFFECT_PREDICTIONS",
                           "LOSS_OF_FUNCTION", "DBSNP","CLINVAR_CLINICAL_SIGNIFICANCE", "CLINVAR_MSID","CLINVAR_VARIANT_ORIGIN",
                           "CLINVAR_CONFLICTED", "CLINVAR_PHENOTYPE","VEP_ALL_CONSEQUENCE",
                           "PATHSCORE","PATHRANK", "PATHDOC","PVS1","PSC1","PS1","PMC1","PM1","PM2","PM4","PM5","PP2","PPC1","PP3","BP3","BP4","BMC1","BSC1","BA1","BP1",
                           "N_INSILICO_CALLED","N_INSILICO_DAMAGING","N_INSILICO_TOLERATED","N_INSILICO_SPLICING_NEUTRAL","N_INSILICO_SPLICING_AFFECTED",
                           "GLOBAL_AF_GNOMAD", pcg_report[['pcgr_config']][['popgen']][['vcftag_gnomad']],
                           "GLOBAL_AF_1KG", pcg_report[['pcgr_config']][['popgen']][['vcftag_tgp']],"TIER","TIER_DESCRIPTION",
                           "GENOMIC_CHANGE", "GENOME_VERSION")

  rlogging::message("Generating tiered set of result variants for output in tab-separated values (TSV) file")

  tsv_variants <- data.frame()
  for(tier in c("tier1", "tier2", "tier3A", "tier3B","gwas")){
    if(tier != 'tier3B' & tier != "gwas"){
      predispose_tags <- predispose_tsv_tags
      for(ph in c('cancer_phenotype','noncancer_phenotype')){
        tierset <- data.frame()
        if(nrow(pcg_report[['snv_indel']][['variant_display']][[tier]][[ph]]) > 0){
          tierset <- pcg_report[['snv_indel']][['variant_display']][[tier]][[ph]]
          tierset$VCF_SAMPLE_ID <- sample_name
          if(tier == 'tier1'){
            tierset$TIER <- 'TIER 1'
            if(ph == 'cancer_phenotype'){
              tierset$TIER_DESCRIPTION <- 'Tier 1 - Pathogenic variant (ClinVar) associated with cancer phenotype'
            }else{
              tierset$TIER_DESCRIPTION <- 'Tier 1 - Pathogenic variant (ClinVar) associated with undefined/noncancer phenotype'
            }
          }
          if(tier == 'tier2'){
            tierset$TIER <- 'TIER 2'
            if(ph == 'cancer_phenotype'){
              tierset$TIER_DESCRIPTION <- 'Tier 2 - Likely pathogenic variant (ClinVar) associated with cancer phenotype'
            }else{
              tierset$TIER_DESCRIPTION <- 'Tier 2 - Likely pathogenic variant (ClinVar) associated with undefined/noncancer phenotype'
            }
          }
          if(tier == 'tier3A'){
            tierset$TIER <- 'TIER 3A'
            if(ph == 'cancer_phenotype'){
              tierset$TIER_DESCRIPTION <- 'Tier 3A - Variant of uncertain significance (VUS in ClinVar) associated with cancer phenotype'
            }else{
              tierset$TIER_DESCRIPTION <- 'Tier 3A - Variant of uncertain significance (VUS in ClinVar) associated with undefined/noncancer phenotype'
            }
          }
        }
        if(!is.null(pcg_report[['pcgr_config']][['custom_tags']])){
          if(pcg_report[['pcgr_config']][['custom_tags']][['custom_tags']] != ""){
            tags <- stringr::str_split(pcg_report[['pcgr_config']][['custom_tags']][['custom_tags']],pattern = ",")[[1]]
            for(t in tags){
              t <- stringr::str_trim(t)
              if(t %in% colnames(tierset)){
                predispose_tags <- c(predispose_tags,t)
              }
            }
          }
        }
        if(nrow(tierset) > 0){
          tsv_variants <- as.data.frame(dplyr::bind_rows(tsv_variants, dplyr::select(tierset, dplyr::one_of(predispose_tags))))
        }
      }
    }else{
      predispose_tags <- predispose_tsv_tags
      if(!is.null(pcg_report[['pcgr_config']][['custom_tags']])){
        if(pcg_report[['pcgr_config']][['custom_tags']][['custom_tags']] != ""){
          tags <- stringr::str_split(pcg_report[['pcgr_config']][['custom_tags']][['custom_tags']],pattern = ",")[[1]]
          for(t in tags){
            t <- stringr::str_trim(t)
            if(t %in% colnames(tierset)){
              predispose_tags <- c(predispose_tags,t)
            }
          }
        }
      }

      #cat(predispose_tags,'\n')
      if(nrow(pcg_report[['snv_indel']][['variant_display']][[tier]]) > 0){
        tierset <- pcg_report[['snv_indel']][['variant_display']][[tier]]
        tierset$VCF_SAMPLE_ID <- sample_name
        if(tier == 'tier3B'){
          tierset$TIER_DESCRIPTION <- 'Tier 3B - Unclassified variants'
          tierset$TIER <- 'TIER 3B'
        }
        if(tier == 'gwas'){
          tierset$TIER_DESCRIPTION <- 'GWAS hit'
          tierset$TIER <- 'GWAS'
        }
        if(nrow(tierset) > 0){
          tsv_variants <- as.data.frame(dplyr::bind_rows(tsv_variants, dplyr::select(tierset, dplyr::one_of(predispose_tags))))
        }
      }
    }
  }
  tsv_variants$DBSNP <- unlist(lapply(stringr::str_match_all(tsv_variants$DBSNP,">rs[0-9]{1,}<"),paste,collapse=","))
  tsv_variants$DBSNP <- stringr::str_replace_all(tsv_variants$DBSNP,">|<", "")
  tsv_variants$GENE_NAME <- unlist(lapply(stringr::str_match_all(tsv_variants$GENE_NAME,">.+<"),paste,collapse=","))
  tsv_variants$GENE_NAME <- stringr::str_replace_all(tsv_variants$GENE_NAME,">|<", "")
  tsv_variants$PROTEIN_DOMAIN <- unlist(lapply(stringr::str_match_all(tsv_variants$PROTEIN_DOMAIN,">.+<"),paste,collapse=","))
  tsv_variants$PROTEIN_DOMAIN <- stringr::str_replace_all(tsv_variants$PROTEIN_DOMAIN,">|<", "")
  tsv_variants$PATHDOC <- stringr::str_replace_all(tsv_variants$PATHDOC,"<br>-", ",")
  tsv_variants <- tsv_variants %>% dplyr::distinct()

  pcg_report[['snv_indel']][['variant_set']][['tsv']] <- tsv_variants
  for(t in c('TIER 1','TIER 2','TIER 3A','TIER 3B','GWAS')){
    if(t == 'TIER 1'){
      pcg_report[['snv_indel']][['variant_set']][['tier1']] <- dplyr::filter(tsv_variants, TIER == 'TIER 1')
    }
    if(t == 'TIER 2'){
      pcg_report[['snv_indel']][['variant_set']][['tier2']] <- dplyr::filter(tsv_variants, TIER == 'TIER 2')
    }
    if(t == 'TIER 3A'){
      pcg_report[['snv_indel']][['variant_set']][['tier3A']] <- dplyr::filter(tsv_variants, TIER == 'TIER 3A')
    }
    if(t == 'TIER 3B'){
      pcg_report[['snv_indel']][['variant_set']][['tier3B']] <- dplyr::filter(tsv_variants, TIER == 'TIER 3B')
    }
    if(t == 'GWAS'){
      pcg_report[['snv_indel']][['variant_set']][['gwas']] <- dplyr::filter(tsv_variants, TIER == 'GWAS')
    }
  }
  return(pcg_report)
}

