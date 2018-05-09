
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
        if(nrow(sample_calls) > 0){

          rlogging::message("Assignment of variants to tier 1/tier 2/tier 3/unclassified")

          pcg_report[['snv_indel']][['variant_display']][['tier1']] <-  dplyr::filter(sample_calls, stringr::str_detect(CLINVAR_VARIANT_ORIGIN,"germline") & !is.na(CLINVAR_SIG) & stringr::str_detect("^pathogenic$",CLINVAR_SIG)) %>% dplyr::select(dplyr::one_of(predispose_tier_tags_display)) %>% dplyr::rename(CLINVAR_SIGNIFICANCE = CLINVAR_SIG)


          pcg_report[['snv_indel']][['variant_display']][['tier2']] <-  dplyr::filter(sample_calls, stringr::str_detect(CLINVAR_VARIANT_ORIGIN,"germline") & !is.na(CLINVAR_SIG) & stringr::str_detect("^likely_pathogenic$",CLINVAR_SIG)) %>% dplyr::select(dplyr::one_of(predispose_tier_tags_display)) %>% dplyr::rename(CLINVAR_SIGNIFICANCE = CLINVAR_SIG)
          pcg_report[['snv_indel']][['variant_display']][['tier3']] <-  dplyr::filter(sample_calls, stringr::str_detect(CLINVAR_VARIANT_ORIGIN,"germline") & !is.na(CLINVAR_SIG) & stringr::str_detect("^uncertain_significance$",CLINVAR_SIG)) %>% dplyr::select(dplyr::one_of(predispose_tier_tags_display)) %>% dplyr::rename(CLINVAR_SIGNIFICANCE = CLINVAR_SIG)
          pcg_report[['snv_indel']][['variant_display']][['unclassified']] <-  dplyr::filter(sample_calls, is.na(CLINVAR_SIG) & CODING_STATUS == 'coding' & (is.na(NFE_AF_GNOMAD) | NFE_AF_GNOMAD <= 0.001)) %>% dplyr::select(dplyr::one_of(predispose_tier_tags_display)) %>% dplyr::rename(CLINVAR_SIGNIFICANCE = CLINVAR_SIG)

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
