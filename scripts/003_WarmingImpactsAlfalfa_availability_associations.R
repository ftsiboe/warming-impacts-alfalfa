#-----------------------------------------------
# Preliminaries                              ####
# Availability & Associations (geographically weighted, per county) - pulled OUT
# of the bootstrap. These are built from the FIXED production cross-section
# (spatial_representation: area, yield, cattle) and do NOT depend on the knots,
# the climate scenarios, or the bootstrap, so they are estimated ONCE here:
#   * avail00 = prod00 + m(prod00), the neighbour-weighted (self-excluded) mean of
#     surrounding production via gwkit::estimate_gwlag;
#   * associations cattle~avail00 and cattle~prod00+prod00_LM via gwkit::estimate_gwr,
#   both under all 50 gwkit specs (10 distance presets x 5 kernels) -> per-county
#   MEDIAN consensus (+ mean retained), plus a national global OLS row (fip 00000).
#
# The WARMING-scenario availability (avail05..avail30) still depends on the yield
# impacts and is produced by the bootstrap stage.
#
# Output: output/availability_associations.rds
#         $availability  (fip, prod00, prod00_LM, avail00 + *_specmean)
#         $associations  (fip, name, est, est_specmean, est_specsd, se, tv, pv, sign_agreement)
#-----------------------------------------------
rm(list=ls(all=TRUE)); gc()
library(magrittr); library(future.apply); library(tidyverse); library(data.table)
library(terra); library(sp); library(sf)
study_environment <- readRDS("data/study_environment.rds")
invisible(lapply(list.files("scripts/helpers", pattern = "[.]R$", full.names = TRUE), source))

# gwkit (sibling package): install if missing/older so parallel workers match.
sysname   <- tolower(as.character(Sys.info()[["sysname"]]))
gwkit_src <- if (grepl("windows", sysname)) file.path(dirname(dirname(getwd())), "packages/gwkit") else file.path(dirname(getwd()), "packages/gwkit")
.src_ver  <- tryCatch(as.package_version(read.dcf(file.path(gwkit_src, "DESCRIPTION"))[, "Version"]), error = function(e) NULL)
.inst_ver <- tryCatch(utils::packageVersion("gwkit"), error = function(e) NULL)
# probe the API so a stale install (same version, missing estimate_gwlag_kernels)
# is refreshed even when the DESCRIPTION version was not bumped.
.has_api  <- tryCatch("estimate_gwlag_kernels" %in% getNamespaceExports("gwkit"), error = function(e) FALSE)
if (is.null(.inst_ver) || (!is.null(.src_ver) && .inst_ver < .src_ver) || !.has_api) {
  message("Installing local gwkit from ", gwkit_src, " ...")
  devtools::install(gwkit_src, upgrade = FALSE, quick = TRUE, quiet = TRUE)
}
library(gwkit)
set.seed(08032024)

#-----------------------------------------------
# GW specifications: 10 distance presets x 5 kernels (metric-outer)  ####
metrics     <- gw_distance_metric_names()
kernels_all <- c("gaussian","exponential","bisquare","boxcar","tricube")

#-----------------------------------------------
# Fixed production cross-section + county centroids  ####
sp_data <- as.data.frame(readRDS("data/spatial_representation.rds"))
names(sp_data)[names(sp_data) %in% "inventory"] <- "cattle"
sp_data$fip    <- stringr::str_pad(as.character(sp_data$fip), 5, pad = "0")
sp_data$prod00 <- sp_data$area * sp_data$yield
sp_data$prod00 <- ifelse(sp_data$prod00 %in% c(NA, Inf, NaN, -Inf), 0, sp_data$prod00)

USMUR    <- rast(file.path(study_environment$gssurgo_archive, "MURASTER_30m.tif"))
Counties <- vect(file.path(study_environment$usaPolygons_archive, "USA_Counties.shp"))
Counties <- crop(project(Counties, crs(USMUR)), ext(USMUR)); rm(USMUR); gc()
Counties$fip <- as.character(paste0(stringr::str_pad(as.numeric(as.character(Counties$STATEFP)), 2, pad = "0"),
                                    stringr::str_pad(as.numeric(as.character(Counties$COUNTYFP)), 3, pad = "0")))
cxy <- geom(centroids(Counties))[, c("x","y")]
centroids_all <- data.frame(fip = Counties$fip, longitude = cxy[, "x"], latitude = cxy[, "y"], stringsAsFactors = FALSE)
sp_data <- dplyr::inner_join(sp_data, centroids_all, by = "fip")

#-----------------------------------------------
# One distance metric: availability + associations across its 5 kernels   ####
# estimate_gwlag_kernels builds the obs->target distance ONCE per metric and reuses
# it across the 5 kernels (availability). The associations (cattle ~ availability)
# use each kernel's OWN avail00, so they loop the kernels; within a spec b1 selects
# the adaptive bandwidth and b2 reuses it - same response/coords, so exact.
estimate_metric <- function(dm_name){
  tryCatch({
    src  <- sp_data[, c("fip","longitude","latitude","prod00")]
    lagk <- as.data.frame(gwkit::estimate_gwlag_kernels(
      data = src, unit = "fip", value_cols = "prod00",
      coords = c("longitude","latitude"), predict = centroids_all,
      distance_metric = dm_name, kernel = kernels_all,
      adaptive = FALSE, bw = NULL, bandwidth = "per_kernel",
      bw_response = "prod00", include_self = FALSE))

    lapply(kernels_all, function(km){
      tryCatch({
        lk <- lagk[lagk$kernel == km, c("fip","prod00_LM")]
        av <- dplyr::left_join(as.data.frame(centroids_all["fip"]), sp_data[, c("fip","prod00")], by = "fip")
        av <- dplyr::inner_join(av, lk, by = "fip")
        data.table::setDT(av)
        av[, prod00  := ifelse(is.na(prod00), 0, prod00)]
        av[, avail00 := prod00 + prod00_LM]
        av <- av[is.finite(avail00)]
        availability <- as.data.frame(av)[, c("fip","prod00","prod00_LM","avail00")]

        aj <- merge(data.table::as.data.table(availability),
                    data.table::as.data.table(sp_data[, c("fip","longitude","latitude","cattle")]), by = "fip")
        aj <- aj[is.finite(cattle) & is.finite(avail00) & is.finite(prod00) & is.finite(prod00_LM)]
        aj[, `:=`(lcattle = log(cattle + 1e-9), lavail00 = log(avail00 + 1e-9),
                  lprod00 = log(prod00 + 1e-9), lprod00_LM = log(prod00_LM + 1e-9))]
        aj <- as.data.frame(aj)

        b1  <- gwkit::estimate_gwr(aj, unit = "fip", formula = lcattle ~ lavail00,
          coords = c("longitude","latitude"), predict = centroids_all,
          distance_metric = dm_name, kernel = km, adaptive = TRUE, bw = NULL, terms = "lavail00")
        bwg <- attr(b1, "bandwidth")                       # reuse in b2 (same response/coords -> same bw)
        b2  <- gwkit::estimate_gwr(aj, unit = "fip", formula = lcattle ~ lprod00 + lprod00_LM,
          coords = c("longitude","latitude"), predict = centroids_all,
          distance_metric = dm_name, kernel = km, adaptive = TRUE, bw = bwg, terms = c("lprod00","lprod00_LM"))
        assoc <- data.table::rbindlist(list(b1, b2), fill = TRUE)
        assoc <- assoc[estimand %in% "mean"]
        assoc[, name := c(lavail00 = "avail00", lprod00 = "prod00", lprod00_LM = "prod00_LM")[term]]
        associations <- assoc[, .(fip = unit_id, name, est = estimate, se, tv, pv)]

        list(availability = availability, associations = as.data.frame(associations))
      }, error = function(e) NULL)
    })
  }, error = function(e) NULL)
}

# metric-outer, SEQUENTIAL (no parallel): 10 metrics x 5 kernels = 50 specs
spec_res <- unlist(lapply(metrics, estimate_metric), recursive = FALSE)
spec_res <- Filter(Negate(is.null), spec_res)
if (length(spec_res) == 0) stop("No GW spec produced results.")

#-----------------------------------------------
# Across-spec consensus (median primary + mean retained)   ####
av <- data.table::rbindlist(lapply(seq_along(spec_res), function(i){
  d <- data.table::as.data.table(spec_res[[i]]$availability); d[, .specN := i]; d }), fill = TRUE)
av_cols <- setdiff(names(av), c("fip",".specN"))
av_med  <- av[, lapply(.SD, stats::median, na.rm = TRUE), by = fip, .SDcols = av_cols]
av_mn   <- av[, lapply(.SD, mean,          na.rm = TRUE), by = fip, .SDcols = av_cols]
data.table::setnames(av_mn, av_cols, paste0(av_cols, "_specmean"))
availability <- as.data.frame(merge(av_med, av_mn, by = "fip"))
nat <- as.data.frame(as.list(colMeans(availability[, -1], na.rm = TRUE))); nat$fip <- "00000"
availability <- rbind(nat[, names(availability)], availability)

as_ <- data.table::rbindlist(lapply(seq_along(spec_res), function(i){
  d <- data.table::as.data.table(spec_res[[i]]$associations); d[, .specN := i]; d }), fill = TRUE)
associations <- as_[, .(est = stats::median(est, na.rm = TRUE),
                        est_specmean = mean(est, na.rm = TRUE),
                        est_specsd = stats::sd(est, na.rm = TRUE),
                        se = stats::median(se, na.rm = TRUE),
                        tv = stats::median(tv, na.rm = TRUE),
                        pv = stats::median(pv, na.rm = TRUE),
                        sign_agreement = mean(sign(est) == sign(stats::median(est, na.rm = TRUE)), na.rm = TRUE)),
                    by = .(fip, name)]

# national (unweighted) global fits on the consensus availability
ajn <- merge(data.table::as.data.table(availability[availability$fip != "00000", c("fip","avail00","prod00","prod00_LM")]),
             data.table::as.data.table(sp_data[, c("fip","cattle")]), by = "fip")
ajn <- ajn[is.finite(cattle) & is.finite(avail00) & is.finite(prod00) & is.finite(prod00_LM)]
ajn[, `:=`(lcattle = log(cattle + 1e-9), lavail00 = log(avail00 + 1e-9),
           lprod00 = log(prod00 + 1e-9), lprod00_LM = log(prod00_LM + 1e-9))]
g1 <- summary(lm(lcattle ~ lavail00, data = ajn))$coef["lavail00", ]
g2 <- summary(lm(lcattle ~ lprod00 + lprod00_LM, data = ajn))$coef[c("lprod00","lprod00_LM"), ]
nat_assoc <- data.table::data.table(
  fip = "00000", name = c("avail00","prod00","prod00_LM"),
  est = c(g1["Estimate"], g2[,"Estimate"]), est_specmean = NA_real_, est_specsd = NA_real_,
  se  = c(g1["Std. Error"], g2[,"Std. Error"]), tv = c(g1["t value"], g2[,"t value"]),
  pv  = c(g1["Pr(>|t|)"], g2[,"Pr(>|t|)"]), sign_agreement = 1)
associations <- rbind(nat_assoc, associations, fill = TRUE)

#-----------------------------------------------
# Save                                       ####
saveRDS(list(availability = availability, associations = as.data.frame(associations)),
        "output/availability_associations.rds")
message("Wrote output/availability_associations.rds: ", nrow(availability),
        " availability rows, ", nrow(associations), " association rows.")
#-----------------------------------------------
