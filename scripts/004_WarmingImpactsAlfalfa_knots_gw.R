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
if(grepl("windows", sysname)){
  devtools::load_all(file.path(dirname(dirname(getwd())),"packages/gwkit"))
}else{
  devtools::load_all(file.path(dirname(getwd()),"packages/gwkit"))
}

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
# restored in the SLOPE stage (006 GWFE block). Set absorb_trends = TRUE to
# globally partial out the trends before the search if strict consistency is
# preferred.
absorb_trends <- FALSE

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

dd_all <- grep("^dday[0-9]+$", names(data), value = TRUE)
avail_dd <- as.integer(gsub("dday", "", dd_all))

panel <- doBy::summaryBy(
  as.formula(paste0("lny + ppt + ppt2 + ", paste(dd_all, collapse = " + "),
                    " ~ fip + state_code + commodity_year")),
  data = data, FUN = mean, keep.names = TRUE, na.rm = TRUE)
setDT(panel)

#-----------------------------------------------
# Candidate (Tmin,Tmax) grid: national +/- band, spread >= 3, clamped ####
pad <- function(x) stringr::str_pad(x, 2, pad = "0")
pairs <- data.table::CJ(Tmin = (Tmin_star - knot_band):(Tmin_star + knot_band),
                        Tmax = (Tmax_star - knot_band):(Tmax_star + knot_band))
pairs <- pairs[(Tmax - Tmin) >= 3 & Tmin %in% avail_dd & Tmax %in% avail_dd]
message("County knot search: ", nrow(pairs), " candidate pairs within +/-", knot_band)

#-----------------------------------------------
# County polygons                            ####
USMUR    <- rast(file.path(study_environment$gssurgo_archive,"MURASTER_30m.tif"))
Counties <- vect(file.path(study_environment$usaPolygons_archive,"USA_Counties.shp"))
Counties <- crop(project(Counties, crs(USMUR)), ext(USMUR)); rm(USMUR); gc()
Counties$fip <- as.character(paste0(stringr::str_pad(as.numeric(as.character(Counties$STATEFP)),2,pad="0"),
                                    stringr::str_pad(as.numeric(as.character(Counties$COUNTYFP)),3,pad="0")))

#-----------------------------------------------
# GW specifications (10 presets x 5 kernels)  ####
spec_gw <- as.data.frame(data.table::rbindlist(
  lapply(gw_distance_metric_names(), function(dm){
    data.frame(distance_metric = dm,
               kernel = c("gaussian","exponential","bisquare","boxcar","tricube"),
               stringsAsFactors = FALSE)
  }), fill = TRUE))
spec_gw$specN <- 1:nrow(spec_gw)

#-----------------------------------------------
# Per-spec county knot selection (parallel)   ####
nw <- max(1, as.numeric(Sys.getenv("SLURM_CPUS_PER_TASK", unset = "1")))
if(nw > 1) future::plan(future::multicore, workers = nw) else future::plan(future::sequential)

fit_pair <- function(tmin, tmax, spec, bw){
  d <- data.table::copy(panel)
  d[, DD1 := pmax(get("dday00")        - get(paste0("dday", pad(tmin))), 0)]
  d[, DD2 := pmax(get(paste0("dday", pad(tmin))) - get(paste0("dday", pad(tmax))), 0)]
  d[, DD3 := pmax(get(paste0("dday", pad(tmax))), 0)]
  out <- tryCatch(estimate_gwfe_coefficients_by_polygon(
    data = d, unit = "fip", polygons = Counties,
    formula = lny ~ DD1 + DD2 + DD3 + ppt + ppt2,
    panel = "fip", time = "commodity_year",
    distance_metric = spec$distance_metric, kernel = spec$kernel,
    adaptive = TRUE, bw = bw, terms = c("DD1","DD2","DD3")), error = function(e) NULL)
  if(is.null(out)) return(NULL)
  attr_bw <<- attr(out, "bandwidth")
  w <- data.table::dcast(out, unit_id + r_squared ~ term, value.var = "est")
  # validity: beneficial-then-harmful sign pattern
  w <- w[is.finite(DD1) & is.finite(DD2) & is.finite(DD3) &
           DD1 >= 0 & DD2 > 0 & DD3 < 0]
  if(nrow(w) == 0) return(NULL)
  w[, `:=`(Tmin = tmin, Tmax = tmax)]
  w[, .(fip = unit_id, r_squared, Tmin, Tmax)]
}

per_spec <- future.apply::future_lapply(1:nrow(spec_gw), function(s){
  spec <- spec_gw[s, ]
  bw_s <- NULL; attr_bw <<- NULL
  acc <- vector("list", nrow(pairs))
  for(pp in 1:nrow(pairs)){
    acc[[pp]] <- fit_pair(pairs$Tmin[pp], pairs$Tmax[pp], spec, bw_s)
    if(is.null(bw_s) && !is.null(attr_bw)) bw_s <- attr_bw   # reuse spec bandwidth
  }
  best <- data.table::rbindlist(acc, fill = TRUE)
  if(nrow(best) == 0) return(NULL)
  best <- best[order(fip, -r_squared)][, .SD[1], by = fip]   # max valid R^2 per county
  best[, specN := spec$specN]
  best
}, future.seed = TRUE)

per_spec <- data.table::rbindlist(per_spec, fill = TRUE)

#-----------------------------------------------
# Consensus across specs -> one knot per county ####
consensus <- per_spec[, .(Tmin = as.integer(round(stats::median(Tmin, na.rm = TRUE))),
                          Tmax = as.integer(round(stats::median(Tmax, na.rm = TRUE))),
                          n_spec_valid = .N,
                          r_squared = stats::median(r_squared, na.rm = TRUE)), by = fip]

# counties with no valid county knot fall back to the national anchor
all_fips <- data.table::data.table(fip = as.character(Counties$fip))
consensus <- consensus[all_fips, on = "fip"]
consensus[is.na(Tmin), `:=`(Tmin = Tmin_star, Tmax = Tmax_star, n_spec_valid = 0L)]
consensus[, `:=`(Tmin_national = Tmin_star, Tmax_national = Tmax_star,
                 knot_band = knot_band, crop = main_crop, period = main_period)]

saveRDS(as.data.frame(consensus), "output/optimal_knots_gw.rds")
message("Wrote county GW knots for ", nrow(consensus), " counties; ",
        "fell back to national for ", sum(consensus$n_spec_valid == 0L), " counties.")
#-----------------------------------------------
