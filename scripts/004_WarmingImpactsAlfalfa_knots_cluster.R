#-----------------------------------------------
# Preliminaries                              ####
# County CLIMATE/PRODUCTION partition for the two-stage knot estimator.
#
# Stage 0 (this script): group counties into a small number of agro-climatic
# clusters from features available for (nearly) every county - the weather
# climatology in `data/prism_weather` and the production cross-section underlying
# spatial_Rep.png (`data/spatial_representation.rds`). NO yield time series and NO
# estimated response enters the features (that would be circular and would break
# out-of-sample coverage).
#
# Downstream (004b, next script): pool each cluster's counties and estimate one
# (Tmin, Tmax, DD1, DD2, DD3) per cluster; attach to member counties; smooth the
# SLOPES (knots stay at the cluster level so 005 only needs a small threshold set).
#
# Output: output/knot_clusters.rds - county_fips -> cluster + features, plus the
#         standardized centers and scaling (to assign any county later).
#-----------------------------------------------
rm(list = ls(all = TRUE)); gc()
library(magrittr); library(tidyverse); library(data.table); library(stringr)
study_environment <- readRDS("data/study_environment.rds")
invisible(lapply(list.files("scripts/helpers", pattern = "[.]R$", full.names = TRUE), source))
set.seed(08032024)

#-----------------------------------------------
# Configuration                              ####
crop_code   <- 107L        # prism_weather commodity_code: 107 = alfalfa (332 = other hay)
clim_period <- 0L          # weather window for the climatology: 0 = full agronomic
                           # season (broadest thermal signal). Alternatives: 105:110.
dd_step     <- 5L          # thin the dday curve to every dd_step degrees, so the
                           # (collinear) thermal block does not swamp production
K_grid      <- 4:20        # candidate cluster counts to screen
K_final     <- NA_integer_ # set to a chosen K to lock the partition; NA = auto-pick
min_cluster <- 25L         # identification floor: min counties per cluster
n_start     <- 50L         # k-means restarts
eps         <- 1e-6        # log offset for skewed production features
use_production <- TRUE     # include area/production/cattle descriptors as features
heat_thresh <- 29L         # upper-knot region: a cluster should contain counties whose
                           # climatology reaches ~this dday threshold to identify Tmax
                           # downstream (diagnostic only here)

heat_col <- paste0("dday", stringr::str_pad(heat_thresh, 2, pad = "0"))

#-----------------------------------------------
# Weather climatology from prism_weather      ####
# prism_weather is one row per county_fips x period x commodity_year (already
# area-weighted in 002). Filter to crop + window per file (keeps memory down),
# then average the dday curve over years and keep interannual SDs.
read_one <- function(f){
  d <- tryCatch(readRDS(f), error = function(e) NULL)
  if (is.null(d)) return(NULL)
  data.table::setDT(d)
  if ("commodity_code" %in% names(d)) d <- d[commodity_code %in% crop_code]
  d[period %in% clim_period]
}

weatherlist <- list.files("data/prism_weather", full.names = TRUE)
weatherlist <- weatherlist[grepl(paste0(2002:2022,collapse = "|"),weatherlist)]
weather <- data.table::rbindlist(
  lapply(weatherlist, read_one), fill = TRUE)
if (nrow(weather) == 0) stop("No prism_weather rows after filtering to crop/period.")
stopifnot("precipitation" %in% names(weather), "county_fips" %in% names(weather))
weather[, county_fips := stringr::str_pad(as.character(county_fips), 5, pad = "0")]


# degree-day curve, thinned to every dd_step degrees
all_dd  <- grep("^dday[0-9]+$", names(weather), value = TRUE)
if (length(all_dd) == 0) stop("No dday columns found in prism_weather.")
dd_thr  <- as.integer(gsub("dday", "", all_dd))
dd_cols <- all_dd[dd_thr %in% seq(0, max(dd_thr, na.rm = TRUE), by = dd_step)]
base_dd <- if ("dday00" %in% all_dd) "dday00" else all_dd[which.min(dd_thr)]
keep_dd <- union(dd_cols, if (heat_col %in% all_dd) heat_col else character(0))

clim_mean <- weather[, lapply(.SD, mean, na.rm = TRUE),
                     by = county_fips, .SDcols = c(keep_dd, "precipitation")]
clim_sd <- weather[, .(dday_sd   = stats::sd(get(base_dd), na.rm = TRUE),
                       precip_sd = stats::sd(precipitation, na.rm = TRUE),
                       n_year    = .N), by = county_fips]
clim <- merge(clim_mean, clim_sd, by = "county_fips")
clim[, aridity := precipitation / (get(base_dd) + 1)]   # precip relative to heat

#-----------------------------------------------
# Production descriptors (spatial_Rep data)   ####
sp <- as.data.frame(readRDS("data/spatial_representation.rds"))
sp$county_fips <- stringr::str_pad(as.character(sp$fip), 5, pad = "0")
# log-transform skewed scale variables; yield deliberately excluded (circularity)
sp$l_area   <- log(sp$area       + eps)
sp$l_prod   <- log(sp$production + eps)
sp$l_cattle <- log(sp$inventory  + eps)
prod <- data.table::as.data.table(sp)[, .(county_fips, area, production, inventory,
                                          l_area, l_prod, l_cattle)]

#-----------------------------------------------
# Assemble the feature matrix                 ####
feat <- merge(clim, prod, by = "county_fips", all.x = TRUE)

climate_features <- c(dd_cols, "precipitation", "aridity", "dday_sd", "precip_sd")
prod_features    <- c("l_area", "l_prod", "l_cattle")
feature_cols     <- if (isTRUE(use_production)) c(climate_features, prod_features) else climate_features

# clustered counties must have the climate block (available ~everywhere); missing
# production descriptors are median-imputed so the county is still placed, with a
# flag recording that its production side was imputed.
feat[, has_production := stats::complete.cases(feat[, ..prod_features])]
for (v in prod_features) {
  m <- stats::median(feat[[v]], na.rm = TRUE)
  feat[[v]] <- ifelse(is.finite(feat[[v]]), feat[[v]], m)
}
feat <- feat[stats::complete.cases(feat[, ..climate_features])]

# drop any near-constant feature (e.g. a dday tail that is 0 everywhere) so scale() is safe
X0  <- as.matrix(feat[, ..feature_cols])
sds <- apply(X0, 2, stats::sd, na.rm = TRUE)
const <- names(sds)[!is.finite(sds) | sds < 1e-8]
if (length(const)) {
  message("Dropping near-constant feature(s): ", paste(const, collapse = ", "))
  feature_cols <- setdiff(feature_cols, const)
}

sc <- scale(as.matrix(feat[, ..feature_cols]))
Xs <- sc[, , drop = FALSE]
scale_center <- attr(sc, "scaled:center")
scale_scale  <- attr(sc, "scaled:scale")
message("Feature table: ", nrow(feat), " counties x ", length(feature_cols), " features (",
        sum(feat$has_production), " with observed production).")

#-----------------------------------------------
# K-selection diagnostics                     ####
# Per K: total within-cluster SS, smallest cluster size, clusters below the floor,
# and the share of clusters containing a county that reaches the upper-knot heat
# region (so Tmax is identifiable downstream). Silhouette when `cluster` is present.
has_heat <- if (heat_col %in% names(feat)) feat[[heat_col]] > 0 else rep(NA, nrow(feat))

diag <- data.table::rbindlist(lapply(K_grid, function(k){
  km <- stats::kmeans(Xs, centers = k, nstart = n_start, iter.max = 100)
  sizes   <- as.integer(table(km$cluster))
  heat_ok <- if (all(is.na(has_heat))) NA_real_ else
    mean(tapply(has_heat, km$cluster, any), na.rm = TRUE)
  sil <- NA_real_
  if (requireNamespace("cluster", quietly = TRUE) && nrow(Xs) <= 6000) {
    ss  <- cluster::silhouette(km$cluster, stats::dist(Xs))
    sil <- mean(ss[, "sil_width"])
  }
  data.table::data.table(K = k, tot_withinss = km$tot.withinss,
                         min_size = min(sizes), n_below_floor = sum(sizes < min_cluster),
                         frac_clusters_reach_heat = heat_ok, avg_silhouette = sil)
}), fill = TRUE)
print(diag)

#-----------------------------------------------
# Final partition                             ####
# Auto rule (override by setting K_final): among Ks whose every cluster clears the
# floor, take the best silhouette (or the largest K if silhouette unavailable).
if (is.na(K_final)) {
  ok <- diag[n_below_floor == 0]
  if (nrow(ok) == 0) ok <- diag
  K_final <- if (all(is.na(ok$avg_silhouette))) ok[which.max(K)]$K else ok[which.max(avg_silhouette)]$K
}
message("Chosen K = ", K_final)

km <- stats::kmeans(Xs, centers = K_final, nstart = n_start, iter.max = 100)
feat[, cluster := km$cluster]

#-----------------------------------------------
# Save partition + reusable scaling/centers   ####
out <- list(
  clusters     = as.data.frame(feat[, c("county_fips", "cluster", "has_production",
                                         feature_cols), with = FALSE]),
  centers      = km$centers,             # standardized space
  scale_center = scale_center,
  scale_scale  = scale_scale,
  feature_cols = feature_cols,
  K            = K_final,
  config       = list(crop_code = crop_code, clim_period = clim_period, dd_step = dd_step,
                      min_cluster = min_cluster, heat_thresh = heat_thresh),
  diagnostics  = as.data.frame(diag))
saveRDS(out, "output/knot_clusters.rds")

# cluster profile (means on the ORIGINAL scale) for a quick sanity read
profile <- feat[, c(list(n = .N), lapply(.SD, mean, na.rm = TRUE)),
                by = cluster, .SDcols = feature_cols][order(cluster)]
message("Wrote output/knot_clusters.rds: ", nrow(feat), " counties in ", K_final, " clusters.")
print(profile)

#-----------------------------------------------
# Stage 1: per-cluster knot + coefficient estimation ####
# Incorporate the national anchor + preferred window with the cluster map:
#   * national anchor (003) CENTERS the (Tmin,Tmax) search band (+/- knot_band);
#   * main_period is the data-driven accumulation window for the panel;
#   * each cluster is the POOLING unit - we run the 003 model (county FE + ppt +
#     state trends, 5-fold CV) over the band, pick the lowest-CV valid pair, and
#     keep its DD1-DD3 slopes. Knots/coefs attach to every county in the cluster.
# Pooling (not per-county GW) means we can afford the state trends here.
library(plm)

main_crop         <- "hay_alfalfa"
knot_band         <- 5L         # cluster knots searched within +/- this of national
candidate_periods <- 105:110    # eligible accumulation windows (5-10 months)

# Preferred window (data-driven), same rule as elsewhere in the pipeline
.nk_all     <- as.data.frame(readRDS("output/optimal_knots.rds"))
.nk_all     <- .nk_all[.nk_all$crop %in% main_crop, ]
main_period <- select_preferred_period(
  r2_period         = .nk_all$target_periods,
  r2_value          = .nk_all$R,
  valid_periods     = unique(.nk_all$target_periods),
  candidate_periods = candidate_periods)

# National anchor knot (lowest CV error across the eligible windows)
nk <- as.data.frame(readRDS("output/optimal_knots.rds"))
nk <- nk[nk$crop %in% main_crop & nk$target_periods %in% c(104:112, 0), ]
nk <- nk[order(nk$cv_error), ][1, ]
Tmin_star <- as.integer(nk$Tmin); Tmax_star <- as.integer(nk$Tmax)
message("National anchor knot: Tmin=", Tmin_star, "  Tmax=", Tmax_star, " ; window ", main_period)

# County-year panel at the preferred window, tagged with its cluster
pan <- build_hay_weather_panel(crop = main_crop, target_periods = main_period,
                               prism_weather_directory = "data/prism_weather")
data.table::setDT(pan)
pan[, county_fips := stringr::str_pad(as.character(county_fips), 5, pad = "0")]
cl_map <- data.table::as.data.table(readRDS("output/knot_clusters.rds")$clusters)[, .(county_fips, cluster)]
pan <- merge(pan, cl_map, by = "county_fips")

# Candidate (Tmin,Tmax) grid: national +/- band, spread >= 3, clamped to available dday
padc     <- function(x) stringr::str_pad(x, 2, pad = "0")
avail_dd <- as.integer(gsub("dday", "", grep("^dday[0-9]+$", names(pan), value = TRUE)))
pairs    <- data.table::CJ(Tmin = (Tmin_star - knot_band):(Tmin_star + knot_band),
                           Tmax = (Tmax_star - knot_band):(Tmax_star + knot_band))
pairs    <- pairs[(Tmax - Tmin) >= 3 & Tmin %in% avail_dd & Tmax %in% avail_dd]
message("Cluster knot search: ", nrow(pairs), " candidate pairs within +/-", knot_band)

# One (cluster x pair) pooled FE fit with 5-fold CV (mirrors 003)
fit_pair_cluster <- function(dc, tmin, tmax){
  tryCatch({
    d <- data.table::copy(dc)
    d[, DD1 := pmax(get("dday00")                  - get(paste0("dday", padc(tmin))), 0)]
    d[, DD2 := pmax(get(paste0("dday", padc(tmin))) - get(paste0("dday", padc(tmax))), 0)]
    d[, DD3 := pmax(get(paste0("dday", padc(tmax))), 0)]
    d <- doBy::summaryBy(lny + ppt + ppt2 + Trend + Trend2 + DD1 + DD2 + DD3 ~
                           fip + state_code + commodity_year,
                         data = d, FUN = mean, keep.names = TRUE, na.rm = TRUE)
    data.table::setDT(d)
    for(st in unique(d$state_code)){
      d[, (paste0("trend1", st)) := ifelse(state_code %in% st, Trend, 0)]
      d[, (paste0("trend2", st)) := ifelse(state_code %in% st, Trend2, 0)]
    }
    rhs  <- c("ppt","ppt2","DD1","DD2","DD3", grep("^trend", names(d), value = TRUE))
    form <- stats::as.formula(paste("lny ~", paste(rhs, collapse = " + ")))
    fit  <- plm::plm(form, data = d, index = c("fip","commodity_year"), model = "within")
    co   <- coef(fit)
    res  <- data.frame(Tmin = tmin, Tmax = tmax, R = summary(fit)$r.squared["rsq"],
                       DD1 = co[["DD1"]], DD2 = co[["DD2"]], DD3 = co[["DD3"]])

    # 5-fold CV error (out-of-sample SSE with the county FE added back), as in 003
    d[, fold := sample(1:5, .N, replace = TRUE)]
    pf <- plm::pdata.frame(as.data.frame(d), index = c("fip","commodity_year"), drop.index = TRUE)
    cv <- data.table::rbindlist(lapply(1:5, function(f){
      tryCatch({
        tr <- pf[!pf$fold %in% f, ]; te <- pf[pf$fold %in% f, ]
        ft <- plm::plm(form, data = tr, index = c("fip","commodity_year"), model = "within")
        b  <- coef(ft)
        X  <- model.matrix(form, te)[, names(b), drop = FALSE]
        ps <- as.numeric(X %*% b)
        fe <- plm::fixef(ft)
        a  <- unname(fe[match(as.character(plm::index(te)[[1]]), names(fe))])
        k  <- !is.na(a)
        data.frame(n = sum(k), e = sum((ps[k] + a[k] - te$lny[k])^2, na.rm = TRUE))
      }, error = function(e) NULL)
    }), fill = TRUE)
    cv <- cv[!cv$e %in% c(NaN, Inf, -Inf, NA), ]
    res$cv_error <- if(nrow(cv)) stats::weighted.mean(cv$e, cv$n) else NA_real_
    res
  }, error = function(e) NULL)
}

# Search each cluster; pick lowest CV, then max R, then widest spread (003's rule)
clusters <- sort(unique(pan$cluster))
cluster_knots <- data.table::rbindlist(lapply(clusters, function(cc){
  dc   <- pan[cluster %in% cc]
  grid <- data.table::rbindlist(lapply(seq_len(nrow(pairs)), function(i)
    fit_pair_cluster(dc, pairs$Tmin[i], pairs$Tmax[i])), fill = TRUE)
  # beneficial-then-harmful sign pattern (as in 003)
  grid <- grid[is.finite(DD1) & is.finite(DD2) & is.finite(DD3) &
                 DD1 >= 0 & DD2 > 0 & DD3 < 0]
  if(nrow(grid) == 0)
    return(data.table::data.table(cluster = cc, Tmin = NA_integer_, Tmax = NA_integer_,
                                  DD1 = NA_real_, DD2 = NA_real_, DD3 = NA_real_,
                                  R = NA_real_, cv_error = NA_real_,
                                  n_county = data.table::uniqueN(dc$fip)))
  grid <- grid[order(cv_error, -R, -(Tmax - Tmin))][1]
  data.table::data.table(cluster = cc, Tmin = as.integer(grid$Tmin), Tmax = as.integer(grid$Tmax),
                         DD1 = grid$DD1, DD2 = grid$DD2, DD3 = grid$DD3,
                         R = grid$R, cv_error = grid$cv_error,
                         n_county = data.table::uniqueN(dc$fip))
}), fill = TRUE)
message("Per-cluster knots:"); print(cluster_knots)

#-----------------------------------------------
# Descriptive cluster labels (from the feature profile) ####
# k-means numbers clusters arbitrarily, so derive names from each cluster's mean
# heat + moisture so they travel with the cluster id (verify against the map).
prof_lab <- feat[, .(heat = mean(get(base_dd), na.rm = TRUE),
                     aridity = mean(aridity, na.rm = TRUE)), by = cluster][order(cluster)]
lab <- stats::setNames(paste0("Cluster ", prof_lab$cluster), as.character(prof_lab$cluster))
if (K_final == 4L) {
  cool <- prof_lab$cluster[order(prof_lab$heat)][1]                    # coolest
  arid <- prof_lab[cluster != cool][order(aridity)][1]$cluster         # most arid of the rest
  hot2 <- prof_lab[!cluster %in% c(cool, arid)][order(-heat)]$cluster  # 2 hottest remaining
  lab[as.character(cool)]    <- "Northern Cool-Continental"
  lab[as.character(arid)]    <- "Arid Southwest Desert"
  lab[as.character(hot2[1])] <- "Southern Hot Long-Season"
  lab[as.character(hot2[2])] <- "Humid Temperate Transition"
}
cluster_labels <- data.table::data.table(cluster = prof_lab$cluster,
                                          cluster_label = unname(lab[as.character(prof_lab$cluster)]))
cluster_knots  <- merge(cluster_knots, cluster_labels, by = "cluster", all.x = TRUE)
message("Cluster labels: ", paste(cluster_labels$cluster, cluster_labels$cluster_label,
                                  sep = "=", collapse = "; "))

# Fall back to the national anchor for any cluster with no valid pooled knot
cluster_knots[is.na(Tmin), `:=`(Tmin = Tmin_star, Tmax = Tmax_star)]

# Attach cluster (Tmin,Tmax,DD1-3) to EVERY county in the partition; key on the
# 5-char `fip` that 005/006 expect, with the national anchor recorded for fallback.
county_knots <- merge(cl_map, cluster_knots, by = "cluster", all.x = TRUE)
county_knots[is.na(Tmin), `:=`(Tmin = Tmin_star, Tmax = Tmax_star)]
county_knots[, `:=`(Tmin_national = Tmin_star, Tmax_national = Tmax_star,
                    knot_band = knot_band, crop = main_crop, period = main_period)]
county_knots[, fip := county_fips]
saveRDS(as.data.frame(county_knots), "output/optimal_knots_cluster.rds")
saveRDS(as.data.frame(cluster_knots), "output/optimal_knots_cluster_byzone.rds")
message("Wrote output/optimal_knots_cluster.rds for ", nrow(county_knots), " counties across ",
        nrow(cluster_knots), " clusters (distinct knot pairs: ",
        nrow(unique(cluster_knots[, .(Tmin, Tmax)])), ").")

#-----------------------------------------------
# Map the partition                          ####
tryCatch({
  Counties <- urbnmapr::get_urbn_map("counties", sf = TRUE)
  Counties <- Counties[!Counties$state_abbv %in% c("AK","HI","PR","GU","VI","MP","AS"), ]
  Counties$county_fips <- stringr::str_pad(as.character(Counties$county_fips), 5, pad = "0")
  Counties <- merge(Counties,
                    as.data.frame(feat[, .(county_fips, cluster = factor(cluster))]),
                    by = "county_fips", all.x = TRUE)

  States <- urbnmapr::get_urbn_map("states", sf = TRUE)
  States <- States[!States$state_abbv %in% c("AK","HI","PR","GU","VI","MP","AS"), ]

  pal <- grDevices::hcl.colors(K_final, palette = "Dynamic")
  g <- ggplot2::ggplot() +
    ggplot2::geom_sf(data = Counties, ggplot2::aes(fill = cluster), color = NA) +
    ggplot2::geom_sf(data = States, fill = NA, color = "grey30", linewidth = 0.2) +
    ggplot2::scale_fill_manual(values = pal, na.value = "grey90", name = "Cluster") +
    ggplot2::labs(
      title    = paste0("Agro-climatic knot clusters (K = ", K_final, ")"),
      subtitle = paste0(nrow(feat), " counties; climate (dday curve, ppt/aridity) + production")) +
    ggplot2::theme_void(base_size = 11) +
    ggplot2::theme(legend.position = "right",
                   plot.title = ggplot2::element_text(face = "bold"))
  dir.create("output/exhibits", showWarnings = FALSE, recursive = TRUE)
  ggplot2::ggsave("output/exhibits/knot_clusters_map.png", g, width = 10, height = 6.5, dpi = 300)
  message("Wrote output/exhibits/knot_clusters_map.png")
}, error = function(e) message("Cluster map skipped: ", conditionMessage(e)))
#-----------------------------------------------
