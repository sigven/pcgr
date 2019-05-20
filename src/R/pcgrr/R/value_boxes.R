
#' Function that plots four value boxes with the most important findings in the cancer genome
#'
#' @param pcg_report pcg report with list elements
#' @return p
#'
#'

plot_value_boxes <- function(pcg_report){
  df <- data.frame(
    x = rep(seq(0, 16, 8), 3),
    y = c(rep(1,3), rep(4.5, 3), rep(8,3)),
    h = rep(3, 9),
    w = rep(7, 9),
    info = c(pcg_report[['value_box']][['tmb_tertile']],
             pcg_report[['value_box']][['signatures']],
             pcg_report[['value_box']][['tumor_n']],
             pcg_report[['value_box']][['tier1']],
             pcg_report[['value_box']][['tier2']],
             pcg_report[['value_box']][['scna']],
             pcg_report[['value_box']][['tumor_purity']],
             pcg_report[['value_box']][['tumor_ploidy']],
             pcg_report[['value_box']][['msi']]

             ),
    color = factor(1:9)
  )
  df <- dplyr::filter(df, color != 3)

  custom_palette_set1 <- RColorBrewer::brewer.pal(9, "Set1")
  custom_palette_set1 <- custom_palette_set1[-6]

  p <- ggplot2::ggplot(df, ggplot2::aes(x, y, height = h, width = w, label = info, fill = color)) +
    ggplot2::geom_tile() +
    ggplot2::geom_text(color = "white", fontface = "bold", size=7) +
    ggplot2::coord_fixed() +
    ggplot2::scale_fill_manual(values = custom_palette_set1)+
    ggplot2::theme_void() +
    ggplot2::guides(fill = F)

  return(p)
}


#' Function that generates value box data for PCGR report
#'
#' @param pcg_report object with existing PCGR report data elements
#' @param pcgr_data object with PCGR annotation data
#' @param sample_name sample identifier
#' @param pcgr_config Object with PCGR configuration parameters
#'
generate_report_data_value_box <- function(pcg_report, pcgr_data, sample_name, pcgr_config){

  pcg_report_value_box <- pcgrr::init_pcg_report(pcgr_config, sample_name, class = 'value_box')
  rlogging::message('------')
  rlogging::message("Assigning elements to PCGR value boxes")

  if(!pcg_report[['content']][['snv_indel']][['eval']]){
    return(pcg_report_value_box)
  }
  if(pcg_report[['content']][['m_signature']][['eval']]){
    if(!is.null(pcg_report[['content']][['m_signature']][['result']])){
      if(nrow(pcg_report[['content']][['m_signature']][['result']][['cancertypes_aetiologies']]) > 0){
        dominant_signatures <- as.data.frame(pcg_report[['content']][['m_signature']][['result']][['cancertypes_aetiologies']] %>% dplyr::filter(Keyword != 'Aging') %>% head(2) %>% dplyr::group_by(Keyword) %>% dplyr::summarise(sigs = paste(Signature_ID, collapse = ','))) %>% dplyr::mutate(v = paste0(Keyword," (",sigs,")"))
        dominant_aetiologies <- paste(unique(dominant_signatures$v),collapse="\n")
        pcg_report_value_box[['signatures']] <- paste0('Mutational signatures:\n',dominant_aetiologies)
      }
    }
  }

  pcg_report_value_box[['tumor_n']] <- "NA"

  if(pcg_report[['content']][['tumor_purity']][['eval']]){
    if(!is.null(pcg_report[['content']][['tumor_purity']][['estimate']])){
      pcg_report_value_box[['tumor_purity']] <- paste0("Tumor purity:\n",pcg_report[['content']][['tumor_purity']][['estimate']])
    }
  }

  if(pcg_report[['content']][['tumor_ploidy']][['eval']]){
    if(!is.null(pcg_report[['content']][['tumor_ploidy']][['estimate']])){
      pcg_report_value_box[['tumor_ploidy']] <- paste0("Tumor ploidy:\n",pcg_report[['content']][['tumor_ploidy']][['estimate']])
    }
  }

  if(pcg_report[['content']][['tmb']][['eval']]){
    if(!is.null(pcg_report[['content']][['tmb']][['variant_statistic']])){
      pcg_report_value_box[['tmb_tertile']] <- pcg_report[['content']][['tmb']][['variant_statistic']][['tmb_tertile']]
    }
  }
  if(pcg_report[['content']][['msi']][['eval']]){
    if(length(pcg_report[['content']][['msi']][['prediction']]) > 0){
      pcg_report_value_box[['msi']] <- pcg_report[['content']][['msi']][['prediction']][['msi_stats']][['vb']]
    }
  }
  if(!is.null(pcg_report[['content']][['cna']])){
    if(pcg_report[['content']][['cna']][['eval']]){
      if(nrow(pcg_report[['content']][['cna']][['variant_display']][['biomarker']]) > 0){
        pcg_report_value_box[['scna']] <- paste0('SCNAs:\n',paste(unique(head(pcg_report[['content']][['cna']][['variant_display']][['biomarker']]$SYMBOL,2)),collapse=", "))
      }
      else{
        pcg_report_value_box[['scna']] <- 'SCNAs:\nNone of strong\nclinical significance'
      }
    }
  }

  if(pcg_report[['content']][['snv_indel']][['eval']]){
    if(length(pcg_report[['content']][['snv_indel']][['variant_set']]) > 0){
      if(nrow(pcg_report[['content']][['snv_indel']][['variant_set']][['tier1']]) > 0){
        tier1_genes <- unique(unlist(pcg_report[['content']][['snv_indel']][['variant_set']][['tier1']]$SYMBOL))
        pcg_report_value_box['tier1'] <- paste(head(tier1_genes,2),collapse=", ")
        if(length(tier1_genes) > 2){
          pcg_report_value_box['tier1'] <- paste(paste(head(tier1_genes,2),collapse=", "),paste(tier1_genes[3:min(4,length(tier1_genes))],collapse=", "),sep="\n")
          if(length(tier1_genes) > 4){
            pcg_report_value_box['tier1'] <- paste0(pcg_report_value_box['tier1'],"++")
          }
        }
        pcg_report_value_box[['tier1']] <- paste0('Tier 1 variants:\n',pcg_report_value_box[['tier1']])
      }else{
        pcg_report_value_box[['tier1']] <- paste0('Tier 1 variants:\nNone')
      }
      if(nrow(pcg_report[['content']][['snv_indel']][['variant_set']][['tier2']]) > 0){
        tier2_genes <- unique(unlist(pcg_report[['content']][['snv_indel']][['variant_set']][['tier2']]$SYMBOL))
        pcg_report_value_box['tier2'] <- paste(head(tier2_genes,2),collapse=", ")
        if(length(tier2_genes) > 2){
          pcg_report_value_box['tier2'] <- paste(paste(head(tier2_genes,2),collapse=", "),paste(tier2_genes[3:min(4,length(tier2_genes))],collapse=", "),sep="\n")
          if(length(tier2_genes) > 4){
            pcg_report_value_box['tier2'] <- paste0(pcg_report_value_box['tier2'],"++")
          }
        }
        pcg_report_value_box[['tier2']] <- paste0('Tier 2 variants:\n',pcg_report_value_box[['tier2']])
      }else{
        pcg_report_value_box[['tier2']] <- paste0('Tier 2 variants:\nNone')
      }
    }
  }
  pcg_report_value_box$eval <- TRUE

  return(pcg_report_value_box)
}

