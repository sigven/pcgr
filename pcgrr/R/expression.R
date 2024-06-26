#' Function that generates expression data for PCGR report
#'
#' @param ref_data PCGR reference data object
#' @param settings PCGR run/configuration settings
#'
#' @export
generate_report_data_expression <-
  function(ref_data = NULL,
           settings = NULL) {

  pcg_report_expression <-
    pcgrr::init_expression_content()

  pcg_report_expression[["eval"]] <- TRUE

  if(as.logical(settings$conf$expression$similarity_analysis) == TRUE){
    pcg_report_expression[["similarity_analysis"]] <-
      load_expression_similarity(settings = settings)
  }

  if(settings$molecular_data$fname_expression_outliers_tsv != "None" &
     file.exists(settings$molecular_data$fname_expression_outliers_tsv)){

    pcg_report_expression[["outliers"]] <-
     pcgrr::load_expression_outliers(settings = settings,
                                 ref_data = ref_data)
  }

  if(settings$molecular_data$fname_expression_tsv != "None" &
     file.exists(settings$molecular_data$fname_expression_tsv)){

    exp_data <-
      readr::read_tsv(
        settings$molecular_data$fname_expression_tsv,
        show_col_types = F, na = "."
      )

    pcg_report_expression[["expression"]] <- exp_data

    if("SYMBOL" %in% colnames(exp_data) == FALSE |
       "TPM" %in% colnames(exp_data) == FALSE |
       "BIOTYPE" %in% colnames(exp_data) == FALSE){
      pcgrr::log4r_warn(
        "Missing a required column in expression file: SYMBOL, TPM, BIOTYPE")
    }else{

      n_pc <- sum(exp_data$BIOTYPE == "protein_coding")

      if(n_pc > 0){
        pcgrr::log4r_info(
          "Estimating immune contexture of tumor sample from RNA-seq data")
        exp_protein_coding <- exp_data |>
          dplyr::filter(.data$BIOTYPE == "protein_coding") |>
          dplyr::group_by(.data$SYMBOL) |>
          dplyr::summarise(TPM = sum(.data$TPM, na.rm = TRUE)) |>
          dplyr::select(c("SYMBOL", "TPM")) |>
          dplyr::distinct()

        if(NROW(exp_protein_coding) > 0){
          rown <- exp_protein_coding$SYMBOL
          mat <- as.matrix(exp_protein_coding$TPM)
          rownames(mat) <- rown
          colnames(mat) <- "TPM"

          pcg_report_expression[["immune_contexture"]] <-
            suppressMessages(quantiseqr::run_quantiseq(
              expression_data = mat,
              is_tumordata = TRUE,
            ))

          if(is.data.frame(pcg_report_expression[["immune_contexture"]]) &
             "Sample" %in% colnames(pcg_report_expression[["immune_contexture"]])){
            pcg_report_expression[["immune_contexture"]] <-
              pcg_report_expression[["immune_contexture"]] |>
              dplyr::rename(sample_id = "Sample") |>
              dplyr::mutate(sample_id = settings$sample_id)
            rownames(pcg_report_expression[["immune_contexture"]]) <- NULL

            pcg_report_expression[["immune_contexture"]] <-
              pcg_report_expression[["immune_contexture"]] |>
              tidyr::pivot_longer(
                !.data$sample_id, names_to = "method_cell_type",
                values_to = "fraction") |>
              dplyr::mutate(fraction = round(
                .data$fraction, digits = 3)) |>
              dplyr::left_join(
                pcgrr::immune_celltypes, by = "method_cell_type") |>
              dplyr::distinct()

          }
        }
      }
    }

  }


  return(pcg_report_expression)
}


