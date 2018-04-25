
#' Function that plots four value boxes with the most important findings in the cancer genome
#'
#' @param pcg_report pcg report with list elements
#' @return p
#'
#'

plot_value_boxes <- function(pcg_report){
  df <- data.frame(
    x = rep(seq(2, 15, 6.5), 2),
    y = c(rep(2,3), rep(6.5, 3)),
    h = rep(4, 6),
    w = rep(6, 6),
    info = c(pcg_report[['value_box']][['tmb_tertile']],
             pcg_report[['value_box']][['msi']],
             pcg_report[['value_box']][['scna']],
             pcg_report[['value_box']][['signatures']],
             pcg_report[['value_box']][['tier1']],
             pcg_report[['value_box']][['tier2']]),
    color = factor(1:6)
  )

  p <- ggplot2::ggplot(df, ggplot2::aes(x, y, height = h, width = w, label = info, fill = color)) +
    ggplot2::geom_tile() +
    ggplot2::geom_text(color = "white", fontface = "bold", size=7) +
    ggplot2::coord_fixed() +
    ggplot2::scale_fill_brewer(type = "qual",palette = "Dark2") +
    ggplot2::theme_void() +
    ggplot2::guides(fill = F)

  return(p)
}


#' Function that generates value box data for PCGR report
#'
#' @param pcg_report object with existing PCGR report data elements
#' @param pcgr_data object with PCGR annotation data
#' @param pcgr_version PCGR software version
#' @param sample_name sample identifier
#' @param pcgr_config Object with PCGR configuration parameters
#' @param genome_seq BSgenome object
#' @param genome_assembly human genome assembly version
#'
generate_report_data_value_box <- function(pcg_report, pcgr_data, pcgr_version, sample_name, pcgr_config, genome_seq, genome_assembly){

  pcg_report_value_box <- pcgrr::init_pcg_report(pcgr_config, sample_name, pcgr_version, genome_assembly, class = 'value_box')
  rlogging::message('------')
  rlogging::message("Assigning elements to PCGR value boxes")

  if(!pcg_report[['snv_indel']][['eval']]){
    return(pcg_report_value_box)
  }
  if(pcg_report[['m_signature']][['eval']]){
    if(!is.null(pcg_report[['m_signature']][['result']])){
      if(nrow(pcg_report[['m_signature']][['result']][['cancertypes_aetiologies']]) > 0){
        dominant_signatures <- as.data.frame(pcg_report[['m_signature']][['result']][['cancertypes_aetiologies']] %>% dplyr::filter(Keyword != 'Aging') %>% head(2) %>% dplyr::group_by(Keyword) %>% dplyr::summarise(sigs = paste(Signature_ID, collapse = ','))) %>% dplyr::mutate(v = paste0(Keyword," (",sigs,")"))
        dominant_aetiologies <- paste(unique(dominant_signatures$v),collapse="\n")
        pcg_report_value_box[['signatures']] <- paste0('Mutational signatures:\n',dominant_aetiologies)
      }
    }
  }

  if(pcg_report[['tmb']][['eval']]){
    if(!is.null(pcg_report[['tmb']][['variant_statistic']])){
      pcg_report_value_box[['tmb_tertile']] <- pcg_report[['tmb']][['variant_statistic']][['tmb_tertile']]
    }
  }
  if(pcg_report[['msi']][['eval']]){
    if(length(pcg_report[['msi']][['prediction']]) > 0){
      pcg_report_value_box[['msi']] <- pcg_report[['msi']][['prediction']][['msi_stats']][['vb']]
    }
  }
  if(!is.null(pcg_report[['cna']])){
    if(pcg_report[['cna']][['eval']]){
      if(nrow(pcg_report[['cna']][['variant_display']][['biomarker']]) > 0){
        pcg_report_value_box[['scna']] <- paste0('SCNAs:\n',paste(unique(head(pcg_report[['cna']][['variant_display']][['biomarker']]$SYMBOL,2),collapse=", ")))
      }
      else{
        pcg_report_value_box[['scna']] <- 'SCNAs:\nNone of strong\nclinical significance'
      }
    }
  }

  if(pcg_report[['snv_indel']][['eval']]){
    if(length(pcg_report[['snv_indel']][['variant_set']]) > 0){
      if(nrow(pcg_report[['snv_indel']][['variant_set']][['tier1']]) > 0){
        tier1_genes <- unique(unlist(pcg_report[['snv_indel']][['variant_set']][['tier1']]$SYMBOL))
        pcg_report_value_box['tier1'] <- paste(head(tier1_genes,2),collapse=", ")
        if(length(tier1_genes) > 2){
          pcg_report_value_box['tier1'] <- paste(paste(head(tier1_genes,2),collapse=", "),paste(tier1_genes[3:min(4,length(tier1_genes))],collapse=", "),sep="\n")
          if(length(tier1_genes) > 4){
            pcg_report_value_box['tier1'] <- paste0(pcg_report_value_box['tier1'],"++")
          }
        }
        pcg_report_value_box[['tier1']] <- paste0('Tier 1 variants:\n',pcg_report_value_box[['tier1']])
      }else{
        pcg_report_value_box[['tier1']] <- paste0('Tier 1 variants:0\n')
      }
      if(nrow(pcg_report[['snv_indel']][['variant_set']][['tier2']]) > 0){
        tier2_genes <- unique(unlist(pcg_report[['snv_indel']][['variant_set']][['tier2']]$SYMBOL))
        pcg_report_value_box['tier2'] <- paste(head(tier2_genes,2),collapse=", ")
        if(length(tier2_genes) > 2){
          pcg_report_value_box['tier2'] <- paste(paste(head(tier2_genes,2),collapse=", "),paste(tier2_genes[3:min(4,length(tier2_genes))],collapse=", "),sep="\n")
          if(length(tier2_genes) > 4){
            pcg_report_value_box['tier2'] <- paste0(pcg_report_value_box['tier2'],"++")
          }
        }
        pcg_report_value_box[['tier2']] <- paste0('Tier 2 variants:\n',pcg_report_value_box[['tier2']])
      }else{
        pcg_report_value_box[['tier2']] <- paste0('Tier 2 variants:0\n')
      }
    }
  }
  pcg_report_value_box$eval <- TRUE

  return(pcg_report_value_box)
}
