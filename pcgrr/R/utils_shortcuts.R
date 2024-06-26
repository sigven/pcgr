#' Tidy eval helpers
#'
#' <https://cran.r-project.org/web/packages/dplyr/vignettes/programming.html>
#'
#' @name tidyeval
#' @keywords internal
#' @importFrom rlang .data :=
NULL

utils::globalVariables(c("."))

#' @noRd
dummy1 <- function() {
  # Solves R CMD check: Namespaces in Imports field not imported from
  DT::datatable
  htmltools::br
  plotly::plot_ly
  caret::createModel
  randomForest::importance
  shiny::bootstrapPage
  bslib::page_fillable
  crosstalk::filter_select
}
