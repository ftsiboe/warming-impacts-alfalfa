# 301_article_objects.R
# Compute EVERY number quoted inline in the article from the analysis outputs and write
# them to narrative/article_objects.json. The master article and its section children
# read this file via  objs <- jsonlite::fromJSON("article_objects.json")  so that no
# result is hard-coded in the prose.
#
# Inputs (produced by scripts 006_...summary.R and 100_...exhibits.R):
#   output/summary/summary_impact_yield.rds   national + county yield-shock means
#   output/summary/summary_piecewise.rds      piecewise degree-day coefficients
#   data/spatial_representation.rds            county cattle/alfalfa panel (correlations)
#   output/exhibits/figure_data/table1.csv     summary statistics (Table 1)
#   output/exhibits/figure_data/spatial_Rep.csv state production/inventory means
#   output/exhibits/figure_data/regression_coefficients.csv  (Table 2)
#
# Sourced after 300_article_helpers.R. Run from the repository ROOT via run_article.R.
# Uses base R + data.table + jsonlite.

suppressPackageStartupMessages({library(data.table); library(jsonlite)})

# Allow standalone runs: load helper paths/constants/formatters if not already sourced
# (run_article.R sources 300 first, so this is a no-op in the normal flow).
if (!exists("OBJECTS_JSON")) source(file.path("scripts", "300_article_helpers.R"))
source(file.path("scripts", "helpers", "select_preferred_period.R"))

rd_csv <- function(f) utils::read.csv(file.path(FIGDATA, f), check.names = FALSE,
                                      stringsAsFactors = FALSE)
rd_rds <- function(f) as.data.table(readRDS(file.path(SUMMARY, f)))

## --- Consensus specification (no per-spec filtering) -----------------------
# 005/006 already reduce the 50 GW specs to a per-county CONSENSUS, so the summary
# tables carry no p/theta/longlat/DistName/kernel columns. keep_gw() is therefore an
# identity pass-through (the former optimal_gw single-spec pick was removed).
keep_gw <- function(dt) dt

## --- Data-driven preferred window + endogenous thresholds ------------------
# Preferred weather-accumulation window: the window (5-10 months) with the highest
# in-sample R-squared among those that produced valid degree-day knots in 003. This is
# computed from the results and OVERRIDES the default PREFERRED_PERIOD from 300, matching
# the identical selection in 100_..._exhibits.R so figures and text agree.
ok_all <- as.data.table(readRDS(file.path(OUTPUT, "optimal_knots.rds")))
valid_periods <- unique(ok_all[crop == PREFERRED_CROP, target_periods])
r2tab <- rd_csv("regression_coefficients.csv")
r2tab <- r2tab[r2tab$name == "r.squared", c("period", "Estimate")]
PREFERRED_PERIOD <- select_preferred_period(r2tab$period, r2tab$Estimate,
                                            valid_periods, WINDOW_PERIODS)

# Endogenous thresholds actually used by the fitted model: 005 and 100 select a single
# knot pair by minimum CV error over the pooled window set {104:112, 0}. Mirror that here
# so the article's stated thresholds equal the ones the model and figures actually use.
kn    <- ok_all[crop == PREFERRED_CROP & target_periods %in% c(104:112, 0)][order(cv_error)]
knots <- list(lower = assert_present(kn$Tmin[1], "knots$lower"),
              upper = assert_present(kn$Tmax[1], "knots$upper"))
# Degree-day segment labels, built from the thresholds (used by relation + both tables).
dd_labels <- c(
  DD1 = paste0("below ", fmt_degC(knots$lower)),
  DD2 = paste0(fmt_degC(knots$lower), en_dash, fmt_degC(knots$upper)),
  DD3 = paste0("above ", fmt_degC(knots$upper)))

# Cross-check: window that minimizes out-of-sample CV error among the valid-knot candidates
# (each window's selected-knot cv_error lives in optimal_knots.rds). Used to show the preferred
# window is robust to the selection metric.
cvtab   <- ok_all[crop == PREFERRED_CROP & target_periods %in% WINDOW_PERIODS]
cv_pref <- cvtab$target_periods[which.min(cvtab$cv_error)]

# objs$window_selection - transparency on how the preferred window was chosen.
window_selection <- list(
  criterion     = "highest in-sample R-squared among windows with valid degree-day knots",
  chosen_period = PREFERRED_PERIOD,
  chosen_months = WINDOW_MONTHS[match(PREFERRED_PERIOD, WINDOW_PERIODS)],
  chosen_r2     = round(r2tab$Estimate[r2tab$period == PREFERRED_PERIOD][1], 3),
  cv_chosen_months = WINDOW_MONTHS[match(cv_pref, WINDOW_PERIODS)],
  robust_to_cv  = identical(as.integer(cv_pref), as.integer(PREFERRED_PERIOD)),
  n_candidates  = length(WINDOW_PERIODS),
  n_valid       = length(intersect(valid_periods, WINDOW_PERIODS)))

# Note (not a hard guardrail): the preferred window is chosen by in-sample R-squared; the
# CV-optimal window may differ. The Methods text now reports BOTH windows (the metrics
# differ only marginally across the candidates), so a mismatch is expected - just log it.
if (!window_selection$robust_to_cv)
  message(sprintf(paste0("note: preferred window by R-squared (%d mo) differs from the ",
                         "minimum-CV window (%d mo); Methods reports both."),
                  window_selection$chosen_months, window_selection$cv_chosen_months))

## --- Yield-shock impacts (summary_impact_yield.rds) ------------------------
# Columns: crop,period,climate_base,region,state_code,county_code,fip,
#   warming_scenario,Estimate_0000,Estimate_mean,Estimate_sd,Estimate_n (consensus -
#   no spec cols). boot_summary names the point estimate Estimate_0000 (NOT bare
#   Estimate); alias it to Estimate so the rest of this script reads one name.
# Estimate is the point (boot "0000") percentage change in yield; national row has
# county_code == 0 (fip "00000"). warming_scenario in {0.5,1.0,1.5,2.0,2.5,3.0}.
iy <- keep_gw(rd_rds("summary_impact_yield.rds"))
iy[, Estimate := Estimate_0000]
iy[, warming_scenario := as.numeric(warming_scenario)]

# The NATIONAL aggregate is the row with county_code==0 AND state_code==0 AND
# region=="". Filtering on county_code==0 alone also matches the per-state and
# per-region aggregates (all have county_code==0); some of those are NaN, and row
# order is not guaranteed, so [1] could pick a NaN. Pin all three keys.
nat <- function(cr, per, base, scen)
  iy[crop == cr & period == per & climate_base == base &
       county_code == 0 & state_code == 0 & region == "" &
       warming_scenario == scen, Estimate][1]

# objs$impact - national mean decline by warming scenario (preferred spec)
impact <- list()
for (i in seq_along(SCENARIOS))
  impact[[SCEN_KEYS[i]]] <- assert_present(
    nat(PREFERRED_CROP, PREFERRED_PERIOD, PREFERRED_BASE, SCENARIOS[i]),
    paste0("impact$", SCEN_KEYS[i]))
impact$preferred_window_months <-
  WINDOW_MONTHS[match(PREFERRED_PERIOD, WINDOW_PERIODS)]

# objs$impact_state - per-state mean over counties, by scenario, preferred spec.
xw <- unique(rd_csv("spatial_Rep.csv")[, c("state_code", "State.Name", "STUSPS")])
st_map <- setNames(xw$STUSPS, xw$state_code)
state_roll <- function(scen) {
  d <- iy[crop == PREFERRED_CROP & period == PREFERRED_PERIOD &
            climate_base == PREFERRED_BASE & county_code != 0 &
            warming_scenario == scen,
          .(val = mean(Estimate, na.rm = TRUE)), by = state_code]
  d[, STUSPS := st_map[as.character(state_code)]]
  d[!is.na(STUSPS)]
}
impact_state <- list()
for (i in seq_along(SCENARIOS)) {
  d <- state_roll(SCENARIOS[i])
  for (j in seq_len(nrow(d)))
    impact_state[[d$STUSPS[j]]][[SCEN_KEYS[i]]] <- round(d$val[j], 3)
}

# objs$county_share - per-region share of counties with a positive +1C impact (and the
# national total), so the "X of Y counties positive" statements are data-driven.
cs <- iy[crop == PREFERRED_CROP & period == PREFERRED_PERIOD &
           climate_base == PREFERRED_BASE & county_code != 0 &
           warming_scenario == 1.0 & region != "",
         .(n = .N, n_pos = sum(Estimate > 0, na.rm = TRUE)), by = region]
cs[, pct_pos := round(100 * n_pos / n, 0)]
county_share <- list()
for (i in seq_len(nrow(cs)))
  county_share[[cs$region[i]]] <- list(n = cs$n[i], n_pos = cs$n_pos[i],
                                       pct_pos = cs$pct_pos[i])
county_share$national <- list(
  n = sum(cs$n), n_pos = sum(cs$n_pos),
  pct_pos = round(100 * sum(cs$n_pos) / sum(cs$n), 0))

## --- Guardrails: prose state claims must still match the data ---------------
# The section prose names specific states as gaining or losing under warming. Region
# membership is geographic (stable), but the SIGN of each state's impact is a result and
# could drift on a re-estimate. These checks stop the render with a clear message if any
# named state flips, so the prose is corrected deliberately rather than shipping a
# contradiction. (Ranking/superlative claims are guarded by sign only; see the report note
# on the Northeast "most significant increases" list.)
.imp <- function(st, scen) { v <- impact_state[[st]][[scen]]
  if (is.null(v)) NA_real_ else v }
assert_state_sign <- function(st, dir, scen = "s10") {
  v <- .imp(st, scen)
  if (is.na(v))
    stop(sprintf("guardrail: prose names %s at %s but it has no impact value.", st, scen),
         call. = FALSE)
  if ((dir == "gain") != (v > 0))
    stop(sprintf(paste0("guardrail: prose calls %s a %s at %s but its impact is %+.3f%%. ",
                        "Update the prose (or the claim) to match the data."),
                 st, dir, scen, v), call. = FALSE)
  invisible(TRUE)
}
for (s in c("CT","VT","RI","MA","ME","MN","WI","MI","WA","OR","MT"))
  assert_state_sign(s, "gain", "s10")
for (s in c("TX","LA","OK","MS","FL","OH","SD","IA","IN","IL",
            "ID","UT","CO","CA","MO","KS"))   # WY dropped: now +0.0% (rounding-level), no longer a loss
  assert_state_sign(s, "loss", "s10")
for (s in c("CT","ME","MA")) assert_state_sign(s, "gain", "s05")   # +0.5C gradient
for (s in c("AL","AZ","FL")) assert_state_sign(s, "loss", "s05")
for (s in c("AZ","FL"))      assert_state_sign(s, "loss", "s15")   # +1.5C gradient
for (s in c("FL","TX","MO","KS")) assert_state_sign(s, "loss", "s30")  # +3.0C gradient

# objs$baseline - +1 deg C national impact under alternative baseline climate periods.
baseline <- list()
for (b in BASELINE_KEYS)
  baseline[[b]] <- assert_present(nat(PREFERRED_CROP, PREFERRED_PERIOD, b, 1.0),
                                  paste0("baseline$", b))

# objs$window - +1 deg C national impact by accumulation window, plus non-alfalfa hay.
window <- list()
for (i in seq_along(WINDOW_PERIODS)) {
  key <- paste0("m", WINDOW_MONTHS[i])
  window[[key]] <- assert_present(
    nat(PREFERRED_CROP, WINDOW_PERIODS[i], PREFERRED_BASE, 1.0),
    paste0("window$", key))
}
window$hay_other <- nat("hay_other", PREFERRED_PERIOD, PREFERRED_BASE, 1.0)

## --- Piecewise degree-day response (summary_piecewise.rds) -----------------
# name in {ppt,ppt2,DD1,DD2,DD3,...}. boot_summary emits Estimate_0000 (point
# estimate) + Estimate_sd (bootstrap SE) and no longer stores p_value, so alias the
# coefficient to Estimate and reconstruct p_value from the bootstrap SE
# (two-sided normal), matching the regression_coefficients.csv built in 100.
# DD1/DD2/DD3 = below 14 / 14-29 / above 29 deg C segments.
pw <- keep_gw(rd_rds("summary_piecewise.rds"))
pw <- pw[crop == PREFERRED_CROP & period == PREFERRED_PERIOD]
pw[, Estimate := Estimate_0000]
pw[, p_value  := 2 * stats::pnorm(-abs(Estimate_0000 / Estimate_sd))]
pget <- function(nm, col) pw[name == nm, ..col][[1]][1]
relation <- list()
for (seg in DD_SEGMENTS) {
  relation[[seg]] <- list(
    label    = unname(dd_labels[seg]),
    estimate = assert_present(pget(seg, "Estimate"), paste0("relation$", seg)),
    p_value  = pget(seg, "p_value"),
    stars    = sig_stars(pget(seg, "p_value")))
}
# Precipitation turning point (inches): -b_ppt / (2*b_ppt2).
b_ppt  <- pget("ppt", "Estimate"); b_ppt2 <- pget("ppt2", "Estimate")
relation$ppt_turning_point <- round(-b_ppt / (2 * b_ppt2), 3)

## --- Regression table (Table 2: coefficients x six aggregation windows) -----
coef_raw <- rd_csv("regression_coefficients.csv")
# County count derived from the preferred national spec (unique county fips, excluding
# the county_code == 0 national aggregate) rather than hard-coded, so it tracks the panel.
N_COUNTIES <- length(unique(
  iy[crop == PREFERRED_CROP & period == PREFERRED_PERIOD &
       climate_base == PREFERRED_BASE & county_code != 0, fip]))
star <- function(p) ifelse(is.na(p) | p >= 0.10, "",
                    ifelse(p < 0.01, "***", ifelse(p < 0.05, "**", "*")))
cget <- function(nm, per, col) {
  v <- coef_raw[[col]][coef_raw$name == nm & coef_raw$period == per]
  if (length(v) == 0) NA_real_ else v[1]
}
# Coefficients and SEs are reported x1000 (matching the manuscript scale).
cell_coef <- function(nm, per) {
  e <- cget(nm, per, "Estimate"); se <- cget(nm, per, "StdError"); p <- cget(nm, per, "p_value")
  if (is.na(e)) return("")
  sprintf("%s%s (%s)", formatC(e * 1000, format = "f", digits = 3),
          star(p), formatC(se * 1000, format = "f", digits = 3))
}
cell_int  <- function(nm, per) { v <- cget(nm, per, "Estimate")
  if (is.na(v)) "" else formatC(v, format = "d", big.mark = ",") }
cell_num  <- function(nm, per, d = 3) { v <- cget(nm, per, "Estimate")
  if (is.na(v)) "" else formatC(v, format = "f", digits = d) }
cell_test <- function(nm, per) { v <- cget(nm, per, "Estimate"); p <- cget(nm, per, "p_value")
  if (is.na(v)) "" else paste0(formatC(v, format = "f", digits = 3), star(p)) }

coef_rows <- list(
  c(paste0("Degree-day (", dd_labels["DD1"], ")"), "DD1",   "coef"),
  c(paste0("Degree-day (", dd_labels["DD2"], ")"), "DD2",   "coef"),
  c(paste0("Degree-day (", dd_labels["DD3"], ")"), "DD3",   "coef"),
  c("Precipitation (inches)",              "ppt",          "coef"),
  c("Precipitation (inches) squared",      "ppt2",         "coef"),
  c("No. of observations",                 "n",            "int"),
  c("No. of counties",                     "",             "cty"),
  c("R-squared",                           "r.squared",    "num"),
  c("Test for no degree-day fixed effect", "test_temp",    "test"),
  c("Test for no precipitation effect",    "test_ppt",     "test"),
  c("Test for no weather effect",          "test_weather", "test"),
  c("Test for no county fixed effect",     "test_fe",      "test"))
mk <- function(nm, kind) vapply(WINDOW_PERIODS, function(per)
  switch(kind, coef = cell_coef(nm, per), int = cell_int(nm, per),
               num = cell_num(nm, per),  test = cell_test(nm, per),
               cty = formatC(N_COUNTIES, big.mark = ",")), character(1))
coef_tab <- as.data.frame(
  do.call(rbind, lapply(coef_rows, function(r) c(r[1], mk(r[2], r[3])))),
  stringsAsFactors = FALSE)
names(coef_tab) <- c(" ", paste0(WINDOW_MONTHS, " months"))
# Mark the data-driven preferred window column with the dagger (not a hard-coded column).
.pref_idx <- match(PREFERRED_PERIOD, WINDOW_PERIODS)
names(coef_tab)[.pref_idx + 1] <- paste0(WINDOW_MONTHS[.pref_idx], " months††")
coef <- list(tab = coef_tab)

## --- Summary statistics (Table 1) ------------------------------------------
t1 <- rd_csv("table1.csv")
sget <- function(v) t1$mean[t1$variable == v][1]
summ <- list(
  n_counties = N_COUNTIES,                   # derived from the panel (see above)
  rsq        = round(cget("r.squared", PREFERRED_PERIOD, "Estimate"), 3),
  area_mean  = round(sget("area"), 2),
  yield_mean = round(sget("yield"), 2),
  prod_mean  = round(sget("production"), 2),
  cattle_mean = round(sget("cattle"), 2),
  ppt_mean   = round(sget("ppt"), 2),
  dd1_mean   = round(sget("DD1"), 2),
  dd2_mean   = round(sget("DD2"), 2),
  dd3_mean   = round(sget("DD3"), 2))

# Table 1 as a labeled two-column frame (value with sd in parentheses), matching the
# manuscript layout — no raw R variable codes.
t1_lab <- c(area       = "Alfalfa production area (1,000 acres)",
            production = "Alfalfa production output (1,000 tons)",
            yield      = "Alfalfa yield (tons/acre)",
            cattle     = "Cattle inventory (1,000 head)",
            ppt        = "Alfalfa growing season precipitation (inches)",
            DD1        = paste0("Alfalfa growing season degree-day (", dd_labels["DD1"], ")"),
            DD2        = paste0("Alfalfa growing season degree-day (", dd_labels["DD2"], ")"),
            DD3        = paste0("Alfalfa growing season degree-day (", dd_labels["DD3"], ")"))
sd_get <- function(v) t1$sd[t1$variable == v][1]
n_obs  <- cget("n", PREFERRED_PERIOD, "Estimate")
summ$tab <- rbind(
  data.frame(
    Variables = unname(t1_lab),
    `Mean (Standard deviation)` = vapply(names(t1_lab), function(v)
      sprintf("%s (%s)",
              formatC(sget(v),  format = "f", digits = 2, big.mark = ","),
              formatC(sd_get(v), format = "f", digits = 2, big.mark = ",")),
      character(1)),
    check.names = FALSE, row.names = NULL),
  data.frame(
    Variables = c("Number of counties", "Number of observations"),
    `Mean (Standard deviation)` = c(formatC(summ$n_counties, big.mark = ","),
                                    formatC(n_obs, format = "d", big.mark = ",")),
    check.names = FALSE))

## --- Spatial production means by state and region (spatial_Rep.csv) --------
sr <- as.data.table(rd_csv("spatial_Rep.csv"))
# Region crosswalk taken from the impact object (state_code -> region).
reg_map <- unique(iy[county_code != 0 & region != "", .(state_code, region)])
sr <- merge(sr, reg_map, by = "state_code", all.x = TRUE)
spatial <- list(state = list(), region = list())
for (i in seq_len(nrow(sr))) {
  k <- sr$STUSPS[i]; if (is.na(k)) next
  spatial$state[[k]] <- list(area = round(sr$`area.mean`[i], 2),
                             production = round(sr$`production.mean`[i], 2),
                             yield = round(sr$`yield.mean`[i], 2))
}
for (rg in unique(stats::na.omit(sr$region))) {
  d <- sr[region == rg]
  spatial$region[[rg]] <- list(
    area = round(mean(d$`area.mean`, na.rm = TRUE), 2),
    production = round(mean(d$`production.mean`, na.rm = TRUE), 2),
    yield = round(mean(d$`yield.mean`, na.rm = TRUE), 2))
}

## --- Cattle-alfalfa correlations (data/spatial_representation.rds) ----------
# County-level correlations quoted in the intro/discussion (production ~0.53,
# yield ~0.32). Computed directly from the county panel, matching the calculation in
# 100_WarmingImpactsAlfalfa_exhibits.R:
#   cor(spa_rep[c("inventory","production","area","yield")], use = "pairwise.complete.obs")
sr_cty <- as.data.table(readRDS(file.path(DATA, "spatial_representation.rds")))
cmat   <- stats::cor(as.matrix(sr_cty[, .(inventory, production, area, yield)]),
                     use = "pairwise.complete.obs")
assoc <- list(
  cor_production = assert_present(round(cmat["inventory", "production"], 2), "assoc$cor_production"),
  cor_yield      = assert_present(round(cmat["inventory", "yield"], 2),      "assoc$cor_yield"))

## --- Availability + cattle warming impacts (summary_availability / summary_cattle) ---
# County-level % impacts (Estimate_0000 / cattleA_0000). These summaries carry NO national
# aggregate row, so the national figure is the mean over counties by scenario (matching the
# .national() reduction in 100_..._exhibits.R). Feeds Figures 7-10 and the availability text.
av <- keep_gw(rd_rds("summary_availability.rds"))
ca <- keep_gw(rd_rds("summary_cattle.rds"))
av[, warming_scenario := as.numeric(warming_scenario)]
ca[, warming_scenario := as.numeric(warming_scenario)]
nat_mean <- function(dt, val, scen)
  mean(dt[crop == PREFERRED_CROP & period == PREFERRED_PERIOD & climate_base == PREFERRED_BASE &
            !fip %in% c("0", "00000") & warming_scenario == scen, get(val)], na.rm = TRUE)
pct_neg <- function(dt, val)
  round(100 * mean(dt[crop == PREFERRED_CROP & period == PREFERRED_PERIOD &
                        climate_base == PREFERRED_BASE & !fip %in% c("0", "00000") &
                        warming_scenario == 1.0, get(val)] < 0, na.rm = TRUE), 0)

impact_avail <- list(); impact_cattle <- list()
for (i in seq_along(SCENARIOS)) {
  impact_avail[[SCEN_KEYS[i]]]  <- assert_present(nat_mean(av, "Estimate_0000", SCENARIOS[i]),
                                                  paste0("impact_avail$", SCEN_KEYS[i]))
  impact_cattle[[SCEN_KEYS[i]]] <- assert_present(nat_mean(ca, "cattleA_0000", SCENARIOS[i]),
                                                  paste0("impact_cattle$", SCEN_KEYS[i]))
}

# Cattle responsiveness to alfalfa availability + its own/neighbour production split
# (003 cross-sectional elasticities in availability_associations$associations; boot-invariant).
aa_assoc <- as.data.frame(readRDS(file.path(OUTPUT, "availability_associations.rds"))$associations)
aa_assoc <- aa_assoc[!aa_assoc$fip %in% c("0", "00000") & is.finite(aa_assoc$est), ]
el <- function(nm) stats::median(aa_assoc$est[aa_assoc$name == nm], na.rm = TRUE)
livestock <- list(
  responsiveness       = assert_present(round(el("avail00"), 2), "livestock$responsiveness"),
  own_production       = round(el("prod00"), 2),
  neighbour_production = round(el("prod00_LM"), 2),
  avail_pct_neg_s10    = pct_neg(av, "Estimate_0000"),
  cattle_pct_neg_s10   = pct_neg(ca, "cattleA_0000"))

## --- Assemble + write ------------------------------------------------------
objs <- list(
  meta = list(
    crop = PREFERRED_CROP, period = PREFERRED_PERIOD,
    window_months = WINDOW_MONTHS[match(PREFERRED_PERIOD, WINDOW_PERIODS)],
    baseline = PREFERRED_BASE, generated = as.character(Sys.time())),
  summary = summ, knots = knots, window_selection = window_selection,
  relation = relation, coef = coef, impact = impact,
  impact_state = impact_state, county_share = county_share, spatial = spatial,
  baseline = baseline, window = window, assoc = assoc,
  impact_avail = impact_avail, impact_cattle = impact_cattle, livestock = livestock)

jsonlite::write_json(objs, OBJECTS_JSON, pretty = TRUE, auto_unbox = TRUE, digits = NA)
message("301_article_objects: wrote ", OBJECTS_JSON,
        "  (+1C national = ", round(impact$s10, 3), "%,",
        " preferred window = ", objs$meta$window_months, " months)")
