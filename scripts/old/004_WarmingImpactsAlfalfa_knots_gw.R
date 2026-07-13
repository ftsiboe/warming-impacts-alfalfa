#-----------------------------------------------
# Preliminaries                              ####
# County-level degree-day KNOTS by geographically weighted fixed-effects
# regression. Run ONCE on the full sample (not bootstrapped): for the main
# analysis cell, search (Tmin,Tmax) pairs within +/- knot_band of the national
# knot, fit the demeaned FE model locally under all 50 gwkit specs, and per
# county pick the valid pair with the best local within-R^2. The selected knots
# are then consensused across specs into one Tmin(i)/Tmax(i) map.
#
# Output: output/optimal_knots_gw.rds  (fip, Tmin, Tmax, + diagnostics)
# Runs after 003 (national knot) and before 005 stage-2 (which trims the PRISM
# climate dday columns to the +/- knot_band band).
#-----------------------------------------------
rm(list=ls(all=TRUE));gc();library(magrittr);library(future.apply);library(tidyverse);library(data.table)
library(plm);library(terra);library(GWmodel);library(sp);library(sf)
study_environment <- readRDS("data/study_environment.rds")
invisible(lapply(list.files("scripts/helpers", pattern = "[.]R$", full.names = TRUE), source))
sysname  <- tolower(as.character(Sys.info()[["sysname"]]))
user     <- Sys.info()[["user"]]
udomain  <- Sys.info()[["udomain"]]

gwkit_src <- if (grepl("windows", sysname) && grepl("TSIB",toupper(udomain))) {
  file.path(dirname(dirname(getwd())), "packages/gwkit")
} else {
  file.path(dirname(getwd()), "gwkit")
}

devtools::install(gwkit_src, upgrade = FALSE, quick = TRUE, quiet = TRUE)

library(gwkit)

options(future.globals.maxSize = 8000 * 1024^2)

set.seed(08032024)

#-----------------------------------------------
# Configuration                              ####
knot_band         <- 5        # county knots searched within +/- this of national
main_crop         <- "hay_alfalfa"
candidate_periods <- 105:110  # eligible accumulation windows (5-10 months)

# Main analysis window: data-driven, not hard-coded. Among the candidate windows
# that produced valid knots in 003, pick the one with the highest in-sample R^2 -
# the same rule (select_preferred_period) used elsewhere in the pipeline.
.nk_all     <- as.data.frame(readRDS("output/optimal_knots.rds"))
.nk_all     <- .nk_all[.nk_all$crop %in% main_crop, ]
main_period <- select_preferred_period(
  r2_period         = .nk_all$target_periods,
  r2_value          = .nk_all$R,
  valid_periods     = unique(.nk_all$target_periods),
  candidate_periods = candidate_periods)
message("Preferred window (data-driven): period ", main_period)

# NOTE: the knot SEARCH absorbs the county fixed effect (demeaning) and keeps
# ppt/ppt2 controls, but omits the state-specific time trends for tractability
# (a 100-column local WLS x 40 pairs x 50 specs x 2,500 counties is infeasible).
# R^2-based knot selection is driven by the degree-day fit; the state trends are
# restored in the SLOPE stage (006 GWFE block).

#-----------------------------------------------
# National anchor knot (from 003)            ####
nk <- as.data.frame(readRDS("output/optimal_knots.rds"))
nk <- nk[nk$crop %in% main_crop & nk$target_periods %in% c(104:112, 0), ]
nk <- nk[order(nk$cv_error), ][1, ]
Tmin_star <- as.integer(nk$Tmin); Tmax_star <- as.integer(nk$Tmax)
message("National anchor knot: Tmin=", Tmin_star, "  Tmax=", Tmax_star)

#-----------------------------------------------
# County-year panel (full sample)            ####
data <- build_hay_weather_panel(
  crop = main_crop, target_periods = main_period,
  prism_weather_directory = "data/prism_weather")

dd_all   <- grep("^dday[0-9]+$", names(data), value = TRUE)
avail_dd <- as.integer(gsub("dday", "", dd_all))

panel <- doBy::summaryBy(
  as.formula(paste0("lny + ppt + ppt2 + Trend + Trend2 +", paste(dd_all, collapse = " + "),
                    " ~ fip + county_fips + state_code + commodity_year")),
  data = data[state_code %in% 20], FUN = mean, keep.names = TRUE, na.rm = TRUE)
setDT(panel)
#-----------------------------------------------
# National estimate                          ####
nat_panel <- as.data.frame(panel)
nat_panel$DD1 <- nat_panel[,"dday00"] - nat_panel[,paste0("dday",stringr::str_pad(Tmin_star,pad="0",2))]
nat_panel$DD2 <- nat_panel[,paste0("dday",stringr::str_pad(Tmin_star,pad="0",2))] - nat_panel[,paste0("dday",stringr::str_pad(Tmax_star,pad="0",2))]
nat_panel$DD3 <- nat_panel[,paste0("dday",stringr::str_pad(Tmax_star,pad="0",2))]
nat_panel$DD1 <- ifelse(nat_panel$DD1<0,0,nat_panel$DD1); nat_panel$DD2 <- ifelse(nat_panel$DD2<0,0,nat_panel$DD2); nat_panel$DD3 <- ifelse(nat_panel$DD3<0,0,nat_panel$DD3)

for(st in unique(nat_panel$state_code)){
  nat_panel[,paste0("trend1",st)] <- ifelse(nat_panel$state_code %in% st,nat_panel$Trend,0)
  nat_panel[,paste0("trend2",st)] <- ifelse(nat_panel$state_code %in% st,nat_panel$Trend2,0)
}
nat_panel$panel <- nat_panel$fip; nat_panel$time <- nat_panel$commodity_year
nat_panel <- pdata.frame(nat_panel, index = c("panel", "time"), drop.index = TRUE)
fit <- plm(as.formula(paste0("lny ~",paste0(c("ppt","ppt2","DD1","DD2","DD3",
                                              names(nat_panel)[grepl("trend",names(nat_panel))]),collapse = "+"))),
           data = nat_panel,model = "within")
coef_star <- coef(fit)
coef_star[["DD1"]]
#-----------------------------------------------
# Candidate (Tmin,Tmax) grid: national +/- band, spread >= 3, clamped ####
pad <- function(x) stringr::str_pad(x, 2, pad = "0")
pairs <- data.table::CJ(Tmin = (Tmin_star - knot_band):(Tmin_star + knot_band),
                        Tmax = (Tmax_star - knot_band):(Tmax_star + knot_band))
pairs <- pairs[(Tmax - Tmin) >= 3 & Tmin %in% avail_dd & Tmax %in% avail_dd]
message("County knot search: ", nrow(pairs), " candidate pairs within +/-", knot_band)

#-----------------------------------------------
# County polygons                            ####
Counties <- urbnmapr::get_urbn_map("counties", sf = TRUE)
Counties <- Counties[
  !Counties$state_abbv %in% c("AK", "HI", "PR", "GU", "VI", "MP", "AS"),
]

#-----------------------------------------------
# GW specifications (10 presets x 5 kernels)  ####
spec_gw <- as.data.frame(data.table::rbindlist(
  lapply(gw_distance_metric_names(), function(dm){
    data.frame(distance_metric = dm,
               kernel = c("gaussian","exponential","bisquare","boxcar","tricube"),
               stringsAsFactors = FALSE)
  }), fill = TRUE))
spec_gw$specN <- 1:nrow(spec_gw)
n_spec_total  <- nrow(spec_gw)          # full ensemble size (before array subsetting)
kernels_all   <- c("gaussian","exponential","bisquare","boxcar","tricube")

# The SLURM array is now BY DISTANCE METRIC (10 tasks): each task owns one
# distance metric = all 5 kernels. This lets estimate_gwr_kernels() build the
# obs->target distance matrix ONCE per (metric, pair) and reuse it across the 5
# kernels, and select the bandwidth ONCE per metric (reused across all pairs).
metrics <- gw_distance_metric_names()
if(!is.na(as.numeric(Sys.getenv("SLURM_ARRAY_TASK_ID")))){
  mtask   <- rep(as.numeric(Sys.getenv("SLURM_ARRAY_TASK_MIN")):as.numeric(Sys.getenv("SLURM_ARRAY_TASK_MAX")), length = length(metrics))
  metrics <- metrics[mtask %in% as.numeric(Sys.getenv("SLURM_ARRAY_TASK_ID"))]
}

# Ensure the per-spec scratch dir exists (guards against a stale study_environment.rds
# saved before `knots_gw` was added to the wd list in 001).
if(is.null(study_environment$wd$knots_gw))
  study_environment$wd$knots_gw <- file.path(dirname(study_environment$wd$knots), "knots_gw")
dir.create(study_environment$wd$knots_gw, showWarnings = FALSE, recursive = TRUE)

#-----------------------------------------------
# Per-metric county knot selection (distance + bandwidth reused across kernels) ####
# For a candidate (Tmin,Tmax): build DD1/DD2/DD3, fit all 5 kernels in ONE
# estimate_gwr_kernels() call (one distance build, shared across kernels), and
# return the valid per-county within-R^2 for each kernel. The bandwidth is
# selected per-kernel on the FIRST pair and reused for the rest (the bw is chosen
# on the FE-demeaned response intercept-only, so it is identical across pairs).
fit_pair_km <- function(tmin, tmax, distance_metric, bw){
  d <- data.table::copy(panel)
  d[, DD1 := pmax(get("dday00")                  - get(paste0("dday", pad(tmin))), 0)]
  d[, DD2 := pmax(get(paste0("dday", pad(tmin))) - get(paste0("dday", pad(tmax))), 0)]
  d[, DD3 := pmax(get(paste0("dday", pad(tmax))), 0)]
  out <- tryCatch(estimate_gwr_kernels(
    data = d, unit = "county_fips", geometry = Counties,
    formula = lny ~ DD1 + DD2 + DD3 + ppt + ppt2,
    panel = "county_fips", time = "commodity_year",
    distance_metric = distance_metric, kernel = kernels_all,
    adaptive = TRUE, bw = bw, bandwidth = "per_kernel", fit_stats = TRUE,
    terms = c("DD1","DD2","DD3")), error = function(e) NULL)
  if(is.null(out)) return(list(bw = bw, data = NULL))
  bw_new <- attr(out, "bandwidth")                # named per-kernel bandwidths
  w <- data.table::dcast(out[estimand == "mean"],
                         kernel + unit_id + r_squared ~ term, value.var = "estimate")
  # validity: beneficial-then-harmful sign pattern
  # w <- w[is.finite(DD1) & is.finite(DD2) & is.finite(DD3) & DD1 >= 0 & DD2 > 0 & DD3 < 0]
  dat <- if(nrow(w)) w[, .(kernel, county_fips = unit_id, r_squared,
                           Tmin = tmin, Tmax = tmax)] else NULL
  list(bw = bw_new, data = dat)
}

lapply(metrics, function(dm_name){
  # dm_name <- metrics[1]
  specN_map <- spec_gw[spec_gw$distance_metric == dm_name, c("kernel","specN")]
  out_files <- file.path(study_environment$wd$knots_gw,
                         paste0("optimal_gw_knots_", stringr::str_pad(specN_map$specN, 3, pad = "0"), ".rds"))
  if(all(file.exists(out_files))) return(NULL)

  bw_dm <- NULL                          # per-metric bandwidth cache (across pairs)
  acc   <- vector("list", nrow(pairs))
  for(pp in 1:nrow(pairs)){
    # pp <- 1
    r <- fit_pair_km(pairs$Tmin[pp], pairs$Tmax[pp], distance_metric = dm_name, bw = bw_dm)
    if(is.null(bw_dm) && !is.null(r$bw) && !is.null(names(r$bw))) bw_dm <- r$bw
    acc[[pp]] <- r$data
  }

  best <- data.table::rbindlist(acc, fill = TRUE)
  if(nrow(best) == 0) return(NULL)
  # max valid within-R^2 per (kernel, county)
  best <- best[order(kernel, county_fips, -r_squared)][, .SD[1], by = .(kernel, county_fips)]

  # one file per spec (= metric x kernel), matching the downstream naming
  for(kk in unique(best$kernel)){
    sN <- specN_map$specN[specN_map$kernel == kk]
    if(length(sN) != 1L) next
    bk <- best[kernel == kk][, specN := sN]
    saveRDS(bk[, .(county_fips, r_squared, Tmin, Tmax, specN)],
            file.path(study_environment$wd$knots_gw,
                      paste0("optimal_gw_knots_", stringr::str_pad(sN, 3, pad = "0"), ".rds")))
  }
})

#-----------------------------------------------
# Consensus across specs -> one knot per county ####
# Read the per-spec files back (they may be written across several array tasks);
# only build the final consensus once ALL specs are present.
spec_files <- list.files(study_environment$wd$knots_gw,
                         pattern = "^optimal_gw_knots_[0-9]+\\.rds$", full.names = TRUE)
if(length(spec_files) < n_spec_total){
  message("Spec files present: ", length(spec_files), "/", n_spec_total,
          " - skipping consensus (run 004 again after all array tasks finish).")
} else {

  per_spec  <- data.table::rbindlist(lapply(spec_files, readRDS), fill = TRUE)
  consensus <- per_spec[, .(Tmin = as.integer(round(stats::median(Tmin, na.rm = TRUE))),
                            Tmax = as.integer(round(stats::median(Tmax, na.rm = TRUE))),
                            n_spec_valid = .N,
                            r_squared = stats::median(r_squared, na.rm = TRUE)), by = county_fips]

  # counties with no valid county knot fall back to the national anchor
  all_fips <- data.table::data.table(county_fips = as.character(Counties$county_fips))
  consensus <- consensus[all_fips, on = "county_fips"]
  consensus[is.na(Tmin), `:=`(Tmin = Tmin_star, Tmax = Tmax_star, n_spec_valid = 0L)]
  consensus[, `:=`(Tmin_national = Tmin_star, Tmax_national = Tmax_star,
                   knot_band = knot_band, crop = main_crop, period = main_period)]

  # 006 consumes this keyed on `fip` (its own 5-char state+county code); our
  # county_fips holds the identical values, so save under the name 006 expects.
  data.table::setnames(consensus, "county_fips", "fip")
  saveRDS(as.data.frame(consensus), "output/optimal_knots_gw.rds")
  message("Wrote county GW knots for ", nrow(consensus), " counties; ",
          "fell back to national for ", sum(consensus$n_spec_valid == 0L), " counties.")
}
#-----------------------------------------------
