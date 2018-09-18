
#' Function that generates cancer genome report - Tier model pcgr.0
#'
#' @param project_directory name of project directory
#' @param query_vcf2tsv name of gzipped TSV file (vcf2tsv) with annotated query SNVs/InDels
#' @param pcgr_data List of data frames with PCGR data annotations
#' @param pcg_config Object with PCGR configuration parameters
#' @param sample_name sample identifier
#' @param pcgr_version PCGR software version
#' @param genome_assembly human genome assembly version
#'

generate_predisposition_report <- function(project_directory, query_vcf2tsv, pcgr_data, pcgr_config = NULL, sample_name = "SampleX",
                            pcgr_version = "0.1.0", genome_assembly = "grch37"){

  pcg_report <- pcgrr::init_pcg_report(pcgr_config, sample_name, pcgr_version, genome_assembly, class = NULL, pcgr_data = pcgr_data, type = "predispose")

  genome_seq <- BSgenome.Hsapiens.UCSC.hg38
  assembly <- "hg38"
  if(genome_assembly == "grch37"){
    genome_seq <- BSgenome.Hsapiens.UCSC.hg19
    assembly <- "hg19"
  }

  fnames <- list()
  fnames[["tsv"]] <- paste0(project_directory, "/", sample_name, ".cpsr.snvs_indels.tiers.tsv")
  fnames[["json"]] <- paste0(project_directory, "/", sample_name, ".cpsr.",genome_assembly, ".json")

  ## define tags/variables to display in data tables (tier 1,2,3A,3B)
  predispose_tier1_2_display <- c("SYMBOL", "CONSEQUENCE", "PROTEIN_CHANGE", "CLINVAR_PHENOTYPE", "GENOTYPE", "GENE_NAME", "PROTEIN_DOMAIN",
                          "HGVSp", "HGVSc", "CDS_CHANGE", "PROTEIN_FEATURE", "PREDICTED_EFFECT", "LOSS_OF_FUNCTION", "DBSNP", "CLINVAR",
                          "CLINVAR_CLINICAL_SIGNIFICANCE", "CLINVAR_VARIANT_ORIGIN", "ONCOGENE", "ONCOSCORE", "TUMOR_SUPPRESSOR",
                          "GLOBAL_AF_GNOMAD", pcg_report[["pcgr_config"]][["popgen"]][["vcftag_gnomad"]], "GLOBAL_AF_1KG",
                          pcg_report[["pcgr_config"]][["popgen"]][["vcftag_tgp"]], "GENOMIC_CHANGE", "GENOME_VERSION")

  predispose_tier3A_display <- unique(c("SYMBOL", "CONSEQUENCE", "PROTEIN_CHANGE", "CLINVAR_PHENOTYPE", "PATHRANK", "GENOTYPE",
                                        "PATHDOC", "PATHSCORE", predispose_tier1_2_display))
  predispose_tier3B_display <- predispose_tier3A_display[!predispose_tier3A_display %in% c("CLINVAR", "CLINVAR_PHENOTYPE",
                                                                                           "CLINVAR_CLINICAL_SIGNIFICANCE",
                                                                                           "CLINVAR_VARIANT_ORIGIN")]

  predispose_tier3B_display <- unique(c("SYMBOL","CONSEQUENCE", "PROTEIN_CHANGE","PATHRANK","PATHSCORE","GENOTYPE","PATHDOC", predispose_tier3B_display))


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

          ## Assign TIER 1 variants to pcg_report object
          for (ph in c("cancer_phenotype", "noncancer_phenotype")){
            pcg_report[["snv_indel"]][["variant_display"]][["tier1"]][[ph]] <-  dplyr::filter(sample_calls_all[[ph]], CLINVAR_CONFLICTED == 0 & stringr::str_detect(CLINVAR_VARIANT_ORIGIN,"germline") & !is.na(CLINVAR_CLINICAL_SIGNIFICANCE) & stringr::str_detect(CLINVAR_CLINICAL_SIGNIFICANCE, "^pathogenic$"))
            if(nrow(pcg_report[["snv_indel"]][["variant_display"]][["tier1"]][[ph]]) > 0){
              pcg_report[["snv_indel"]][["variant_display"]][["tier1"]][[ph]] <-
                pcg_report[["snv_indel"]][["variant_display"]][["tier1"]][[ph]] %>%
                dplyr::mutate(PATHRANK = NA, PATHSCORE = NA, PATHDOC = NA)
            }
          }

          ## Assign TIER 2 variants (likely pathogenic) to pcg_report object
          for (ph in c("cancer_phenotype", "noncancer_phenotype")){
            pcg_report[["snv_indel"]][["variant_display"]][["tier2"]][[ph]] <-  dplyr::filter(sample_calls_all[[ph]], CLINVAR_CONFLICTED == 0 & stringr::str_detect(CLINVAR_VARIANT_ORIGIN,"germline") & !is.na(CLINVAR_CLINICAL_SIGNIFICANCE) & stringr::str_detect(CLINVAR_CLINICAL_SIGNIFICANCE, "^(likely_pathogenic|pathogenic,likely_pathogenic)$"))
            if(nrow(pcg_report[["snv_indel"]][["variant_display"]][["tier2"]][[ph]]) > 0){
              pcg_report[["snv_indel"]][["variant_display"]][["tier2"]][[ph]] <-
                pcg_report[["snv_indel"]][["variant_display"]][["tier2"]][[ph]] %>%
                dplyr::mutate(PATHRANK = NA, PATHSCORE = NA, PATHDOC = NA)
            }
          }

          ## Assign TIER 3A (VUS in ClinVar) variants to pcg_report object
          for (ph in c("cancer_phenotype", "noncancer_phenotype")){
            pcg_report[["snv_indel"]][["variant_display"]][["tier3A"]][[ph]] <-  dplyr::filter(sample_calls_all[[ph]], CLINVAR_CONFLICTED == 0 & stringr::str_detect(CLINVAR_VARIANT_ORIGIN, "germline") & !is.na(CLINVAR_CLINICAL_SIGNIFICANCE) & stringr::str_detect(CLINVAR_CLINICAL_SIGNIFICANCE, "^(uncertain_significance|risk_factor)$"))
            if(nrow(pcg_report[["snv_indel"]][["variant_display"]][["tier3A"]][[ph]]) > 0){
              pcg_report[["snv_indel"]][["variant_display"]][["tier3A"]][[ph]] <-
                pcgrr::assign_pathogenicity_score(pcg_report[["snv_indel"]][["variant_display"]][["tier3A"]][[ph]], pcgr_config, pcgr_data)
              pcg_report[["snv_indel"]][["variant_display"]][["tier3A"]][[ph]] <-
                pcg_report[["snv_indel"]][["variant_display"]][["tier3A"]][[ph]] %>%
                dplyr::arrange(desc(PATHSCORE), LOSS_OF_FUNCTION, CODING_STATUS, desc(ONCOSCORE))
            }
          }

          ## Assign TIER 3B (Non-classified variants) to pcg_report object
          pcg_report[["snv_indel"]][["variant_display"]][["tier3B"]] <-  dplyr::filter(sample_calls_all[["noncancer_phenotype"]], is.na(CLINVAR_CLINICAL_SIGNIFICANCE) & CODING_STATUS == "coding" & (is.na(GLOBAL_AF_GNOMAD) | GLOBAL_AF_GNOMAD < pcgr_config$maf_limits$maf_gnomad) & (is.na(GLOBAL_AF_1KG) | GLOBAL_AF_1KG < pcgr_config$maf_limits$maf_tgp))
          if(nrow(pcg_report[["snv_indel"]][["variant_display"]][["tier3B"]]) > 0){
            pcg_report[["snv_indel"]][["variant_display"]][["tier3B"]] <-
              pcgrr::assign_pathogenicity_score(pcg_report[["snv_indel"]][["variant_display"]][["tier3B"]], pcgr_config, pcgr_data)
            pcg_report[["snv_indel"]][["variant_display"]][["tier3B"]] <-
              pcg_report[["snv_indel"]][["variant_display"]][["tier3B"]] %>%
              dplyr::arrange(desc(PATHSCORE), LOSS_OF_FUNCTION, CODING_STATUS, desc(ONCOSCORE))
          }


          pcg_report <- pcgrr::generate_tier_tsv_cpsr(pcg_report, sample_name = sample_name)
          pcg_report <- pcgrr::summary_findings_cpsr(pcg_report)
          population_tags <- unique(c("GLOBAL_AF_1KG", pcg_report[["pcgr_config"]][["popgen"]][["vcftag_tgp"]],
                                      "GLOBAL_AF_GNOMAD", pcg_report[["pcgr_config"]][["popgen"]][["vcftag_gnomad"]]))
          for (class in c("tier1", "tier2", "tier3A", "tier3B")){
            if (class != "tier3B"){
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
                pcg_report[["snv_indel"]][["variant_display"]][[class]] <- dplyr::select(pcg_report[["snv_indel"]][["variant_display"]][[class]], dplyr::one_of(predispose_tier3B_display))
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
    rmarkdown::render(system.file("templates", "report_predispose.Rmd", package = "pcgrr"), output_format = rmarkdown::html_document(theme = pcg_report[["pcgr_config"]][["visual"]][["report_theme"]], toc = T, toc_depth = 3, toc_float = T, number_sections = F, includes = rmarkdown::includes(after_body = "disclaimer.md")), output_file = paste0(sample_name, ".cpsr.", genome_assembly, ".html"), output_dir = project_directory, clean = T, intermediates_dir = project_directory, quiet = T)
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
