#' Function that assigns one of six mutation types to a list of mutations
#'
#' @param var_df data frame with variants
#'
#' @return var_df
#'
#' @export
assign_mutation_type <- function(var_df) {

  invisible(
    assertthat::assert_that(
      is.data.frame(var_df),
      msg = "Argument 'var_df' must be a valid data.frame"))
  assertable::assert_colnames(var_df, c("VARIANT_CLASS", "REF", "ALT"),
                              only_colnames = F, quiet = T)
  var_df <- var_df |>
    dplyr::mutate(
      MUTATION_TYPE =
        dplyr::case_when(VARIANT_CLASS == "SNV" & REF == "G" &
                           ALT == "A" ~ "C>T:G>A",
                         VARIANT_CLASS == "SNV" & REF == "C" &
                           ALT == "T" ~ "C>T:G>A",
                         VARIANT_CLASS == "SNV" & REF == "G" &
                           ALT == "C" ~ "C>G:G>C",
                         VARIANT_CLASS == "SNV" & REF == "C" &
                           ALT == "G" ~ "C>G:G>C",
                         VARIANT_CLASS == "SNV" & REF == "C" &
                           ALT == "A" ~ "C>A:G>T",
                         VARIANT_CLASS == "SNV" & REF == "G" &
                           ALT == "T" ~ "C>A:G>T",
                         VARIANT_CLASS == "SNV" & REF == "A" &
                           ALT == "G" ~ "A>G:T>C",
                         VARIANT_CLASS == "SNV" & REF == "T" &
                           ALT == "C" ~ "A>G:T>C",
                         VARIANT_CLASS == "SNV" & REF == "A" &
                           ALT == "T" ~ "A>T:T>A",
                         VARIANT_CLASS == "SNV" & REF == "T" &
                           ALT == "A" ~ "A>T:T>A",
                         VARIANT_CLASS == "SNV" & REF == "A" &
                           ALT == "C" ~ "A>C:T>G",
                         VARIANT_CLASS == "SNV" & REF == "T" &
                           ALT == "G" ~ "A>C:T>G",
                         FALSE ~ as.character(NA)))
  return(var_df)
}

