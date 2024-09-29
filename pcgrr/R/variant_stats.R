#' Function that computes various variant statistics from a data frame
#' with variant records
#'
#' @param var_df data frame with variants
#' @param pct_other_limit numeric value specifying the percentage limit
#' for the 'Other' category
#'
#' @export
#'
get_variant_statistics <- function(var_df = NULL, pct_other_limit = 4){

  assertthat::assert_that(
    !is.null(var_df),
    is.data.frame(var_df),
    msg = "Argument 'var_df' must be a valid data.frame"
  )

  assertable::assert_colnames(
    var_df, c("VARIANT_CLASS", "CONSEQUENCE","CODING_STATUS"),
    only_colnames = F, quiet = T
  )

  consequence_stats <-
    var_df |>
    dplyr::mutate(CONSEQUENCE = stringr::str_replace_all(
      .data$CONSEQUENCE, "(, [0-9A-Za-z_]{1,}){1,}$",""
    )) |>
    dplyr::group_by(.data$CONSEQUENCE) |>
    dplyr::summarise(
      N = dplyr::n(),
      .groups = "drop"
    ) |>
    dplyr::arrange(dplyr::desc(.data$N))

  if(NROW(consequence_stats) > 5) {
    consequence_stats_top <- utils::head(consequence_stats, 4)
    consequence_stats_other <- consequence_stats |>
      dplyr::slice_tail(n = -4) |>
      dplyr::summarise(
        N = sum(.data$N),
        CONSEQUENCE = "other_consequences"
      )
    consequence_stats <- dplyr::bind_rows(
      consequence_stats_top, consequence_stats_other) |>
      dplyr::arrange(dplyr::desc(.data$N))
  }

  consequence_stats <- consequence_stats |>
    dplyr::mutate(Pct = .data$N / sum(.data$N) * 100)

  consequence_stats_coding <-
    var_df |>
    dplyr::filter(.data$CODING_STATUS == "coding")

  if(NROW(consequence_stats_coding) > 0) {
    consequence_stats_coding <-
      consequence_stats_coding |>
      dplyr::mutate(CONSEQUENCE = stringr::str_replace_all(
        .data$CONSEQUENCE, "(, [0-9A-Za-z_]{1,}){1,}$",""
      )) |>
      dplyr::group_by(.data$CONSEQUENCE) |>
      dplyr::summarise(
        N = dplyr::n(),
        .groups = "drop"
      ) |>
      dplyr::arrange(dplyr::desc(.data$N))

    if(NROW(consequence_stats_coding) > 5) {
      consequence_stats_coding_top <- utils::head(consequence_stats_coding, 4)
      consequence_stats_coding_other <- consequence_stats_coding |>
        dplyr::slice_tail(n = -4) |>
        dplyr::summarise(
          N = sum(.data$N),
          CONSEQUENCE = "other_consequences"
        )
      consequence_stats_coding <- dplyr::bind_rows(
        consequence_stats_coding_top,
        consequence_stats_coding_other) |>
        dplyr::arrange(dplyr::desc(.data$N))
    }

    consequence_stats_coding <-
      consequence_stats_coding |>
      dplyr::mutate(Pct = .data$N / sum(.data$N) * 100)
  }


  variant_class_stats <-
    var_df |>
    dplyr::group_by(.data$VARIANT_CLASS) |>
    dplyr::summarise(
      N = dplyr::n(),
      .groups = "drop"
    ) |>
    dplyr::mutate(Pct = .data$N / sum(.data$N) * 100) |>
    dplyr::arrange(dplyr::desc(.data$Pct))

  coding_stats <-
    var_df |>
    dplyr::group_by(.data$CODING_STATUS) |>
    dplyr::summarise(
      N = dplyr::n(),
      .groups = "drop"
    ) |>
    dplyr::mutate(Pct = .data$N / sum(.data$N) * 100) |>
    dplyr::arrange(dplyr::desc(.data$Pct))

  result <- list()
  result[['consequence']] <- consequence_stats
  result[['consequence_coding']] <- consequence_stats_coding
  result[['variant_class']] <- variant_class_stats
  result[['coding']] <- coding_stats

  return(result)
}
