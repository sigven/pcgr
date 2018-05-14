
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

generate_predisposition_report <- function(project_directory, query_vcf2tsv, pcgr_data, pcgr_config = NULL, sample_name = 'SampleX',
                            pcgr_version = '0.1.0', genome_assembly = 'grch37'){

  pcg_report <- pcgrr::init_pcg_report(pcgr_config, sample_name, pcgr_version, genome_assembly, class = NULL, pcgr_data = pcgr_data, type = 'predispose')

  genome_seq = BSgenome.Hsapiens.UCSC.hg38
  assembly <- 'hg38'
  if(genome_assembly == 'grch37'){
    genome_seq = BSgenome.Hsapiens.UCSC.hg19
    assembly = 'hg19'
  }

  fnames <- list()
  fnames[['tsv']] <- paste0(project_directory, '/',sample_name,'.pcgr_predispose.snvs_indels.tiers.tsv')
  fnames[['json']] <- paste0(project_directory, '/',sample_name,'.pcgr_predispose.json')

  predispose_tier_tags_display <- c("SYMBOL",
                         "CONSEQUENCE",
                          "PROTEIN_CHANGE",
                          "GENE_NAME",
                          "PROTEIN_DOMAIN",
                          "HGVSp",
                          "HGVSc",
                          "NFE_AF_GNOMAD",
                          "GLOBAL_AF_GNOMAD",
                          "CDS_CHANGE",
                          "PROTEIN_FEATURE",
                          "PREDICTED_EFFECT",
                          "DBSNP",
                          "CLINVAR",
                          "CLINVAR_SIG",
                          "CLINVAR_VARIANT_ORIGIN",
                          "VEP_ALL_CONSEQUENCE",
                          "ONCOGENE",
                          "TUMOR_SUPPRESSOR",
                          "KEGG_PATHWAY",
                          "GENOMIC_CHANGE",
                          "GENOME_VERSION")


  if(query_vcf2tsv != 'None.gz'){
    if(!file.exists(query_vcf2tsv) | file.size(query_vcf2tsv) == 0){
      rlogging::warning(paste0("File ",query_vcf2tsv," does not exist or has zero size"))
    }
    else{
      if(!is.null(pcgr_config) & query_vcf2tsv != 'None.gz'){
        sample_calls <- pcgrr::get_calls(query_vcf2tsv, pcgr_data, pcgr_version, sample_name, pcgr_config, genome_seq, assembly)
        rlogging::message("Filtering against cancer predisposition genes:")
        sample_calls <- dplyr::semi_join(sample_calls,pcg_report[['snv_indel']][['predisposition_genes']],by=c("SYMBOL" = "symbol"))
        rlogging::message(nrow(sample_calls)," variants remaining")
        gene_hits <- paste(unique(sample_calls$SYMBOL), collapse=", ")
        rlogging::message("Found variants in the following cancer predisposition genes: ", gene_hits)
        pcg_report$snv_indel$eval <- TRUE


        rlogging::message("Limiting variants to hereditary cancer-predisposing syndromes/cancer conditions")
        phenotype_medgen_oncology <- read.table("../../../../../DB/var_annotation_tracks/phenotype_ontology/phenotype_ontology_oncology.tsv", header=T, stringsAsFactors = F, quote="", sep="\t",na.strings=c("","-",NA),comment.char="")
        phenotype_medgen_oncology <- phenotype_medgen_oncology %>% dplyr::select(cui,cui_name) %>% dplyr::distinct()

        sample_calls_per_trait <- tidyr::separate_rows(sample_calls, CLINVAR_MEDGEN_CUI, sep=",")
        sample_calls_per_cancer_trait <- dplyr::inner_join(sample_calls_per_trait, phenotype_medgen_oncology, by=c("CLINVAR_MEDGEN_CUI" = "cui")) %>% dplyr::distinct()
        # exclude 'Not provided' and 'Not specified'
        sample_calls_per_cancer_trait <- sample_calls_per_cancer_trait %>% dplyr::filter(CLINVAR_MEDGEN_CUI != 'CN169374' & CLINVAR_MEDGEN_CUI != 'CN517202')
        sample_calls_per_cancer_trait <- sample_calls_per_cancer_trait %>% dplyr::select(dplyr::one_of(predispose_tier_tags_display),cui_name)
        sample_calls_traits <- as.data.frame(sample_calls_per_cancer_trait %>% dplyr::group_by(predispose_tier_tags_display) %>% dplyr::summarise(PHENOTYPE = paste(unique(cui_name), collapse="; "))) %>% dplyr::distinct()

        # tmp <- dplyr::select(sample_calls, CLINVAR_MEDGEN_CUI, CLINVAR_TRAITS_ALL, VAR_ID, CLINVAR_SIG, CLINVAR_VARIANT_ORIGIN)
        # tmp_all_traits <- tidyr::separate_rows(tmp, CLINVAR_MEDGEN_CUI, sep=",")
        # all_cancer_traits <- dplyr::inner_join(tmp_all_traits, dplyr::select(phenotype_medgen_oncology,-c(parent_cui,top_cui)), by=c("CLINVAR_MEDGEN_CUI" = "cui")) %>% dplyr::distinct()
        # non_cancer_traits <- dplyr::anti_join(tmp_all_traits, dplyr::select(phenotype_medgen_oncology,-c(parent_cui,top_cui)), by=c("CLINVAR_MEDGEN_CUI" = "cui")) %>% dplyr::distinct()
        # ## exclude 'Not provided' and 'Not specified'
        # non_cancer_traits <- non_cancer_traits %>% dplyr::filter(CLINVAR_MEDGEN_CUI != 'CN169374' & CLINVAR_MEDGEN_CUI != 'CN517202')
        # non_cancer_traits <- dplyr::left_join(non_cancer_traits, medgen_map, by=c("CLINVAR_MEDGEN_CUI" = "cui"))
        #
        # all_cancer_traits$predisposition_stratification <- NA
        # all_cancer_traits[!is.na(all_cancer_traits$group) & all_cancer_traits$group == 'Hereditary_Cancer_Susceptibility_NOS',]$predisposition_stratification <- 'Hereditary cancer-predisposing syndromes/cancer susceptibility'
        # all_cancer_traits[!is.na(all_cancer_traits$group) & all_cancer_traits$group != 'Hereditary_Cancer_Susceptibility_NOS',]$predisposition_stratification <- 'Sporadic cancers'
        # all_cancer_traits <- dplyr::select(all_cancer_traits,-c(group,top_name,do_id,do_name)) %>% dplyr::distinct()


        # hereditary_cancer_traits <- dplyr::filter(all_cancer_traits, predisposition_stratification == 'Hereditary cancer-predisposing syndromes/cancer susceptibility') %>% dplyr::filter(!is.na(cui_name))
        # other_cancer_traits <- dplyr::filter(all_cancer_traits, predisposition_stratification != 'Hereditary cancer-predisposing syndromes/cancer susceptibility')
        # other_cancer_traits <- dplyr::anti_join(other_cancer_traits, dplyr::select(hereditary_cancer_traits, VAR_ID, CLINVAR_MEDGEN_CUI)) %>% dplyr::filter(!is.na(cui_name))
        #
        # cancer_traits <- rbind(hereditary_cancer_traits, other_cancer_traits)
        # balle <- as.data.frame(cancer_traits %>% dplyr::group_by(CLINVAR_MEDGEN_CUI, predisposition_stratification, VAR_ID, CLINVAR_SIG, CLINVAR_VARIANT_ORIGIN) %>% dplyr::summarise(cui_names = paste(cui_name, collapse="; ")))

        if(nrow(sample_calls_traits) > 0){

          rlogging::message("Assignment of variants to tier 1/tier 2/tier 3/unclassified")

          pcg_report[['snv_indel']][['variant_display']][['tier1']] <-  dplyr::filter(sample_calls_traits, stringr::str_detect(CLINVAR_VARIANT_ORIGIN,"germline|de_novo|inherited") & !is.na(CLINVAR_SIG) & stringr::str_detect("^pathogenic",CLINVAR_SIG)) %>% dplyr::select(dplyr::one_of(predispose_tier_tags_display)) %>% dplyr::rename(CLINVAR_SIGNIFICANCE = CLINVAR_SIG)


          pcg_report[['snv_indel']][['variant_display']][['tier2']] <-  dplyr::filter(sample_calls_traits, stringr::str_detect(CLINVAR_VARIANT_ORIGIN,"germline|de_novo|inherited") & !is.na(CLINVAR_SIG) & stringr::str_detect("^likely_pathogenic",CLINVAR_SIG)) %>% dplyr::select(dplyr::one_of(predispose_tier_tags_display)) %>% dplyr::rename(CLINVAR_SIGNIFICANCE = CLINVAR_SIG)
          pcg_report[['snv_indel']][['variant_display']][['tier3']] <-  dplyr::filter(sample_calls_traits, stringr::str_detect(CLINVAR_VARIANT_ORIGIN,"germline|de_novo|inherited") & !is.na(CLINVAR_SIG) & stringr::str_detect("^(uncertain_significance|risk_factor)",CLINVAR_SIG)) %>% dplyr::select(dplyr::one_of(predispose_tier_tags_display)) %>% dplyr::rename(CLINVAR_SIGNIFICANCE = CLINVAR_SIG)
          pcg_report[['snv_indel']][['variant_display']][['unclassified']] <-  dplyr::filter(sample_calls_traits, stringr::str_detect(CLINVAR_VARIANT_ORIGIN,"germline|de_novo|inherited")) %>% dplyr::filter(is.na(CLINVAR_SIG) & CODING_STATUS == 'coding' & (is.na(NFE_AF_GNOMAD) | NFE_AF_GNOMAD <= 0.001)) %>% dplyr::select(dplyr::one_of(predispose_tier_tags_display)) %>% dplyr::rename(CLINVAR_SIGNIFICANCE = CLINVAR_SIG)

          for(class in c('tier1','tier2','tier3','unclassified')){
            if(nrow(pcg_report[['snv_indel']][['variant_display']][[class]]) > 0){
              if(!any(!is.na(pcg_report[['snv_indel']][['variant_display']][[class]]$NFE_AF_GNOMAD))){
                pcg_report[['snv_indel']][['variant_display']][[class]]$NFE_AF_GNOMAD <- 0
              }
              if(!any(!is.na(pcg_report[['snv_indel']][['variant_display']][[class]]$GLOBAL_AF_GNOMAD))){
                pcg_report[['snv_indel']][['variant_display']][[class]]$GLOBAL_AF_GNOMAD <- 0
              }
            }
          }
        }
      }
    }
  }

  fname_key <- 'tsv'
  if(!is.null(pcg_report[['snv_indel']][['variant_set']][[fname_key]])){
    if(nrow(pcg_report[['snv_indel']][['variant_set']][[fname_key]]) > 0){
      write.table(pcg_report[['snv_indel']][['variant_set']][[fname_key]],file=fnames[[fname_key]], sep="\t",col.names = T,row.names = F,quote = F)
    }
  }

  pcgr_json <- jsonlite::toJSON(pcg_report, pretty=T,na='string',null='null')
  write(pcgr_json, fnames[['json']])
  rmarkdown::render(system.file("templates","report_predispose.Rmd", package="pcgrr"), output_format = rmarkdown::html_document(theme = pcg_report[['pcgr_config']][['visual']][['report_theme']], toc = T, toc_depth = 3, toc_float = T, number_sections = F), output_file = paste0(sample_name,'.pcgr_predispose.',genome_assembly,'.html'), output_dir = project_directory, clean = T, intermediates_dir = project_directory, quiet=T)


}
