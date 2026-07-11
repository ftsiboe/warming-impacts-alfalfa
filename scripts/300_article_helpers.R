# 300_article_helpers.R
# Shared paths, constants, and formatting helpers for the alfalfa warming-effects
# journal article (narrative). Sourced first by run_article.R. Run with the working
# directory = the repository ROOT (paths below are repo-root relative).
#
# Figures/exhibit CSVs are NOT built here: they are produced by the analysis pipeline
# (scripts 000-006 + 100_WarmingImpactsAlfalfa_exhibits.R) into output/exhibits/ and
# output/summary/. This article layer only consumes those outputs, computes the inline
# numbers (301_article_objects.R), and renders the document (302_render_article.R).

## --- Paths -----------------------------------------------------------------
DATA         <- "data"
OUTPUT       <- "output"
EXHIBITS     <- file.path(OUTPUT, "exhibits")
FIGDATA      <- file.path(EXHIBITS, "figure_data")
SUMMARY      <- file.path(OUTPUT, "summary")
NARRATIVE    <- "narrative"
OBJECTS_JSON <- file.path(NARRATIVE, "article_objects.json")

## --- Analysis constants ----------------------------------------------------
# Preferred specification used for the headline numbers.
#   PREFERRED_PERIOD is a FALLBACK default only. 301_article_objects.R recomputes it
#   data-drivenly via select_preferred_period() - the window (5-10 months) with the highest
#   in-sample R-squared among those that produced valid degree-day knots - and overrides the
#   value below, matching the identical selection in 100_..._exhibits.R. The default here is
#   used only for standalone sourcing of this file.
PREFERRED_CROP   <- "hay_alfalfa"
PREFERRED_PERIOD <- 107L
PREFERRED_BASE   <- "1991_2020"

# Warming scenarios (deg C) and their JSON key suffixes.
SCENARIOS   <- c(0.5, 1.0, 1.5, 2.0, 2.5, 3.0)
SCEN_KEYS   <- c("s05", "s10", "s15", "s20", "s25", "s30")

# Season-window robustness: analysis period code -> window length in months.
WINDOW_PERIODS <- c(105L, 106L, 107L, 108L, 109L, 110L)
WINDOW_MONTHS  <- c(5L,   6L,   7L,   8L,   9L,   10L)

# Baseline-climate robustness (climate_base NAME values in the summary objects).
BASELINE_KEYS <- c("1991_2020", "1981_2010", "1971_2000", "1961_1990")

# Piecewise degree-day segment keys, in order. The human-readable labels are built in
# 301 from the estimated thresholds (objs$knots) so they can never diverge from the
# endogenous cut points in optimal_knots.rds.
DD_SEGMENTS <- c("DD1", "DD2", "DD3")

## --- Formatting helpers ----------------------------------------------------
# Values in the summary objects are already expressed in percent (e.g. -6.72), so the
# percent formatters do NOT multiply by 100.
fmt_num     <- function(x, d = 2) formatC(x, format = "f", digits = d, big.mark = ",")
fmt_pct     <- function(x, d = 1) paste0(formatC(x, format = "f", digits = d), "%")
fmt_abs_pct <- function(x, d = 1) paste0(formatC(abs(x), format = "f", digits = d), "%")
fmt_signed  <- function(x, d = 1) paste0(ifelse(x >= 0, "+", "−"),
                                         formatC(abs(x), format = "f", digits = d), "%")
# Temperature threshold formatter (thresholds are whole degrees Celsius).
fmt_degC    <- function(x, d = 0) paste0(formatC(x, format = "f", digits = d), "°C")
# Scientific-notation formatter on a fixed 10^p scale (default 10^-3) for the small
# degree-day coefficients, so tiny values are legible (e.g. 0.0000 -> 0.025 x 10^-3).
# Returns an inline-math snippet, e.g. "$0.326 \\times 10^{-3}$".
fmt_sci     <- function(x, p = -3L, d = 3)
  paste0("$", formatC(x / 10^p, format = "f", digits = d), " \\times 10^{", p, "}$")
sig_stars   <- function(p) ifelse(p < 0.01, "***", ifelse(p < 0.05, "**",
                            ifelse(p < 0.10, "*", "")))
en_dash     <- "–"

# Guardrail: stop the build if any value destined for the prose is missing.
assert_present <- function(x, what) {
  if (length(x) == 0 || any(is.na(x)))
    stop("article_objects: missing/NA value for '", what,
         "' - upstream exhibit/summary is incomplete.", call. = FALSE)
  invisible(x)
}
