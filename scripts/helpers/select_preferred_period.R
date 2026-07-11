#' Determine the preferred weather-accumulation window (period) from the results
#'
#' The analysis estimates the piecewise degree-day model across several growing-season
#' weather-accumulation windows. Rather than hard-coding which window is "preferred",
#' this function derives it from the estimated results: among the candidate windows that
#' produced valid degree-day knots (i.e. a correctly signed threshold pair survived the
#' selection in `003_..._knots.R` and therefore appears in `optimal_knots.rds`), it
#' returns the window with the highest in-sample R-squared.
#'
#' @param r2_period Integer vector of period codes for which an R-squared is available.
#' @param r2_value  Numeric vector of R-squared values, aligned element-wise with
#'   `r2_period`.
#' @param valid_periods Integer vector of period codes that produced valid knots
#'   (typically `unique(optimal_knots$target_periods[optimal_knots$crop == crop])`).
#' @param candidate_periods Integer vector of allowed windows to choose among
#'   (e.g. `105:110` for the five- to ten-month windows).
#'
#' @return The single period code (integer) with the highest R-squared among the
#'   candidate windows that have valid knots. Errors if the eligible set is empty.
#'
#' @examples
#' \dontrun{
#'   select_preferred_period(r2$period, r2$Estimate, valid, 105:110)
#' }
#'
#' @export
select_preferred_period <- function(r2_period, r2_value, valid_periods, candidate_periods) {
  eligible <- intersect(valid_periods, candidate_periods)
  keep <- r2_period %in% eligible
  if (!any(keep)) {
    stop("select_preferred_period: no candidate window has both valid knots and an ",
         "R-squared value (candidates: ", paste(candidate_periods, collapse = ", "),
         "; valid: ", paste(valid_periods, collapse = ", "), ").", call. = FALSE)
  }
  p <- r2_period[keep]
  v <- r2_value[keep]
  p[which.max(v)]
}
