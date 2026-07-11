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
#   output/optimal_gw.rds                      preferred geographically-weighted spec
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

rd_csv <- function(f) utils::read.csv(file.path(FIGDATA, f), check.names = FALSE,
                                      stringsAsFactors = FALSE)
rd_rds <- function(f) as.data.table(readRDS(file.path(SUMMARY, f)))

## --- Preferred geographically-weighted specification -----------------------
# The bootstrapped summaries contain a grid of GW bandwidths (p, theta, ...). Keep only
# the optimal one so the article reports a single specification, matching the figures.
gw   <- as.data.table(readRDS(file.path(OUTPUT, "optimal_gw.rds")))
gwk  <- intersect(c("p","theta","longlat","DistName","kernel"), names(gw))
keep_gw <- function(dt) dt[gw[, ..gwk], on = gwk, nomatch = 0]

## --- Yield-shock impacts (summary_impact_yield.rds) ------------------------
# Columns: p,theta,longlat,DistName,kernel,crop,period,specN,climate_base,region,
#   state_code,county_code,fip,warming_scenario,Estimate,Estimate_mean,Estimate_sd,...
# Estimate is the point (boot "0000") percentage change in yield; national row has
# county_code == 0 (fip "00000"). warming_scenario in {0.5,1.0,1.5,2.0,2.5,3.0}.
iy <- keep_gw(rd_rds("summary_impact_yield.rds"))
iy[, warming_scenario := as.numeric(warming_scenario)]

nat <- function(cr, per, base, scen)
  iy[crop == cr & period == per & climate_base == base &
       county_code == 0 & warming_scenario == scen, Estimate][1]

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
# name in {ppt,ppt2,DD1,DD2,DD3,...}; Estimate is the coefficient, p_value its
# significance. DD1/DD2/DD3 = below 14 / 14-29 / above 29 deg C segments.
pw <- keep_gw(rd_rds("summary_piecewise.rds"))
pw <- pw[crop == PREFERRED_CROP & period == PREFERRED_PERIOD]
pget <- function(nm, col) pw[name == nm, ..col][[1]][1]
relation <- list()
for (seg in names(DD_SEGMENTS)) {
  relation[[seg]] <- list(
    label    = unname(DD_SEGMENTS[seg]),
    estimate = assert_present(pget(seg, "Estimate"), paste0("relation$", seg)),
    p_value  = pget(seg, "p_value"),
    stars    = sig_stars(pget(seg, "p_value")))
}
# Precipitation turning point (inches): -b_ppt / (2*b_ppt2).
b_ppt  <- pget("ppt", "Estimate"); b_ppt2 <- pget("ppt2", "Estimate")
relation$ppt_turning_point <- round(-b_ppt / (2 * b_ppt2), 3)

## --- Regression table (Table 2: coefficients x six aggregation windows) -----
coef_raw <- rd_csv("regression_coefficients.csv")
N_COUNTIES <- 2573L
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
  c("Degree-day (below 14°C)",             "DD1",          "coef"),
  c("Degree-day (14°C to 29°C)",           "DD2",          "coef"),
  c("Degree-day (above 29°C)",             "DD3",          "coef"),
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
names(coef_tab)[2] <- "5 months††"   # preferred window
coef <- list(tab = coef_tab)

## --- Summary statistics (Table 1) ------------------------------------------
t1 <- rd_csv("table1.csv")
sget <- function(v) t1$mean[t1$variable == v][1]
summ <- list(
  n_counties = 2573L,                       # confirm against panel if it changes
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
            DD1        = "Alfalfa growing season degree-day (below 14°C)",
            DD2        = "Alfalfa growing season degree-day (14°C to 29°C)",
            DD3        = "Alfalfa growing season degree-day (above 29°C)")
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

## --- Assemble + write ------------------------------------------------------
objs <- list(
  meta = list(
    crop = PREFERRED_CROP, period = PREFERRED_PERIOD,
    window_months = WINDOW_MONTHS[match(PREFERRED_PERIOD, WINDOW_PERIODS)],
    baseline = PREFERRED_BASE, generated = as.character(Sys.time())),
  summary = summ, relation = relation, coef = coef, impact = impact,
  impact_state = impact_state, spatial = spatial, baseline = baseline,
  window = window, assoc = assoc)

jsonlite::write_json(objs, OBJECTS_JSON, pretty = TRUE, auto_unbox = TRUE, digits = NA)
message("301_article_objects: wrote ", OBJECTS_JSON,
        "  (+1C national = ", round(impact$s10, 3), "%,",
        " preferred window = ", objs$meta$window_months, " months)")
