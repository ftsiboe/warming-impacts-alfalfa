#-----------------------------------------------
# Preliminaries                              ####
# WARMING-scenario alfalfa availability + cattle shifts, computed OUTSIDE the
# bootstrap.
#
# The neighbourhood operator (distances, kernel weights, bandwidths) does not vary
# by bootstrap draw - only the production values do - so there is no reason to
# re-estimate it inside the boot loop. Here the lag is applied ONCE per
# (distance metric x cell) to EVERY boot's production vector at the same time:
# estimate_gwlag_kernels() builds the obs->target distance once per metric, reuses
# it across the 5 kernels, and lags all value columns in a single BLAS pass.
#
#   availXX(boot) = prodXX(boot) + m(prodXX(boot)),  prodXX = prod00 * (1 + ccXX/100)
#
# Inputs : 003 availability_associations.rds (baseline avail00 + associations)
#          007 consensus_cluster_*.rds       (per-boot impact_yield)
#          data/spatial_representation.rds   (fixed prod00 = area x yield, cattle)
# Outputs: output/availability_warming.rds          (boot x cell x county x scenario)
#          output/summary/summary_availability.rds
#          output/summary/summary_cattle.rds
#-----------------------------------------------
rm(list = ls(all = TRUE)); gc()
library(magrittr); library(tidyverse); library(data.table); library(terra)
study_environment <- readRDS("data/study_environment.rds")
invisible(lapply(list.files("scripts/helpers", pattern = "[.]R$", full.names = TRUE), source))

sysname   <- tolower(as.character(Sys.info()[["sysname"]]))
gwkit_src <- if (grepl("windows", sysname)) file.path(dirname(dirname(getwd())), "packages/gwkit") else file.path(dirname(getwd()), "packages/gwkit")
.src_ver  <- tryCatch(as.package_version(read.dcf(file.path(gwkit_src, "DESCRIPTION"))[, "Version"]), error = function(e) NULL)
.inst_ver <- tryCatch(utils::packageVersion("gwkit"), error = function(e) NULL)
.has_api  <- tryCatch("estimate_gwlag_kernels" %in% getNamespaceExports("gwkit"), error = function(e) FALSE)
if (is.null(.inst_ver) || (!is.null(.src_ver) && .inst_ver < .src_ver) || !.has_api) {
  message("Installing local gwkit from ", gwkit_src, " ...")
  devtools::install(gwkit_src, upgrade = FALSE, quick = TRUE, quiet = TRUE)
}
library(gwkit)
set.seed(08032024)

CELL_KEYS   <- c("crop","period","climate_base")
SFX         <- c("05","10","15","20","25","30")
metrics     <- gw_distance_metric_names()
kernels_all <- c("gaussian","exponential","bisquare","boxcar","tricube")

#-----------------------------------------------
# Fixed cross-section: prod00, cattle, centroids  ####
sp <- as.data.frame(readRDS("data/spatial_representation.rds"))
names(sp)[names(sp) %in% "inventory"] <- "cattle"
sp$fip    <- stringr::str_pad(as.character(sp$fip), 5, pad = "0")
sp$prod00 <- sp$area * sp$yield
sp$prod00 <- ifelse(sp$prod00 %in% c(NA, Inf, NaN, -Inf), 0, sp$prod00)

USMUR    <- rast(file.path(study_environment$gssurgo_archive, "MURASTER_30m.tif"))
Counties <- vect(file.path(study_environment$usaPolygons_archive, "USA_Counties.shp"))
Counties <- crop(project(Counties, crs(USMUR)), ext(USMUR)); rm(USMUR); gc()
Counties$fip <- as.character(paste0(stringr::str_pad(as.numeric(as.character(Counties$STATEFP)), 2, pad = "0"),
                                    stringr::str_pad(as.numeric(as.character(Counties$COUNTYFP)), 3, pad = "0")))
cxy <- geom(centroids(Counties))[, c("x","y")]
centroids_all <- data.frame(fip = Counties$fip, longitude = cxy[, "x"], latitude = cxy[, "y"], stringsAsFactors = FALSE)
sp <- dplyr::inner_join(sp[, c("fip","prod00","cattle")], centroids_all, by = "fip")

#-----------------------------------------------
# 003: baseline availability + associations   ####
aa         <- readRDS("output/availability_associations.rds")
base_avail <- data.table::as.data.table(aa$availability)[, .(fip, prod00_base = prod00,
                                                             prod00_LM, avail00)]
assoc      <- as.data.frame(aa$associations)[, c("fip","name","est")]
assoc$name <- paste0("b_", assoc$name)
assoc      <- data.table::as.data.table(tidyr::spread(assoc, name, est))

#-----------------------------------------------
# 007: per-boot county yield impacts          ####
bfiles <- list.files(study_environment$wd$boots, recursive = TRUE, full.names = TRUE,
                     pattern = "^consensus_cluster_.*\\.rds$")
if (length(bfiles) == 0) stop("No 007 boot files found under ", study_environment$wd$boots)
imp <- data.table::rbindlist(lapply(bfiles, function(f)
  tryCatch(data.table::as.data.table(readRDS(f)$impact_yield), error = function(e) NULL)), fill = TRUE)
imp <- imp[!fip %in% c("0","00000")]                       # county rows only
imp[, fip := stringr::str_pad(as.character(fip), 5, pad = "0")]

#-----------------------------------------------
# Per cell: lag EVERY boot x scenario in one pass   ####
cells <- unique(imp[, ..CELL_KEYS])
message("Cells: ", nrow(cells), " ; boots: ", data.table::uniqueN(imp$boot))

avail_long <- data.table::rbindlist(lapply(seq_len(nrow(cells)), function(ci){
  tryCatch({
    cl <- cells[ci]
    d  <- imp[crop %in% cl$crop & period %in% cl$period & climate_base %in% cl$climate_base]
    if (nrow(d) == 0) return(NULL)

    # long -> wide: one column per (boot x scenario), holding the impact %
    m <- data.table::melt(d[, c("boot","fip", paste0("cc", SFX)), with = FALSE],
                          id.vars = c("boot","fip"), variable.name = "cc", value.name = "impact")
    m[, colnm := paste0("p_", boot, "_", gsub("cc", "", as.character(cc)))]
    w  <- data.table::dcast(m, fip ~ colnm, value.var = "impact")
    vc <- setdiff(names(w), "fip")

    # attach the fixed cross-section, convert impacts -> shocked production
    w <- merge(data.table::as.data.table(sp), w, by = "fip", all.x = TRUE)
    for (v in vc) {
      z <- w[[v]]; z[!is.finite(z)] <- 0
      data.table::set(w, j = v, value = w$prod00 * (1 + z / 100))
    }

    # one call per distance metric: distance built once, reused across the 5 kernels,
    # and all (boot x scenario) columns lagged in a single BLAS pass
    lagk <- data.table::rbindlist(lapply(metrics, function(dm){
      tryCatch(data.table::as.data.table(gwkit::estimate_gwlag_kernels(
        data = as.data.frame(w), unit = "fip", value_cols = vc,
        coords = c("longitude","latitude"), predict = centroids_all,
        distance_metric = dm, kernel = kernels_all, adaptive = FALSE,
        bw = NULL, bandwidth = "per_kernel", bw_response = vc[1],
        include_self = FALSE)), error = function(e) NULL)
    }), fill = TRUE)
    if (nrow(lagk) == 0) return(NULL)

    # availXX = prodXX + m(prodXX); then the median across the 50 specs
    lagl <- data.table::melt(lagk, id.vars = c("fip","kernel"),
                             variable.name = "colnm", value.name = "lag")
    lagl[, colnm := gsub("_LM$", "", as.character(colnm))]
    prodl <- data.table::melt(w[, c("fip", vc), with = FALSE], id.vars = "fip",
                              variable.name = "colnm", value.name = "prod")
    prodl[, colnm := as.character(colnm)]
    # LEFT join onto the targets: extrapolation (non-production) counties get their
    # availability from the neighbour lag alone (own production = 0), matching 007.
    a <- merge(lagl, prodl, by = c("fip","colnm"), all.x = TRUE)
    a[!is.finite(prod), prod := 0]
    a[, avail := prod + lag]
    # keep the consensus lag too: the neighbour-production channel needs it
    a <- a[, .(avail          = stats::median(avail, na.rm = TRUE),
               avail_specmean = mean(avail, na.rm = TRUE),
               prodlag        = stats::median(lag, na.rm = TRUE),
               prod           = prod[1]), by = .(fip, colnm)]

    # unpack boot / scenario out of the column name
    a[, boot := sub("^p_(.*)_[0-9]+$", "\\1", colnm)]
    a[, warming_scenario := as.numeric(sub("^.*_([0-9]+)$", "\\1", colnm)) / 10]
    a[, colnm := NULL]
    data.frame(cl, as.data.frame(a), row.names = NULL)
  }, error = function(e) NULL)
}), fill = TRUE)
if (nrow(avail_long) == 0) stop("No warming availability produced.")

data.table::setDT(avail_long)
avail_long <- merge(avail_long, base_avail[, .(fip, avail00)], by = "fip", all.x = TRUE)
avail_long[, pct := ((avail / avail00) - 1) * 100]         # % change vs baseline availability
saveRDS(as.data.frame(avail_long), "output/availability_warming.rds")

#-----------------------------------------------
# Bootstrap summaries (point estimate + across-boot spread)  ####
boot_summary <- function(dt, value, keys, out_prefix = value){
  dt <- data.table::as.data.table(dt)
  dt[, .obs := ifelse(is.finite(get(value)), 1L, 0L)]
  p0 <- dt[boot %in% "0000", c(keys, value), with = FALSE]
  data.table::setnames(p0, value, paste0(out_prefix, "_0000"))
  bb <- dt[!boot %in% "0000"][
    , setNames(list(mean(get(value), na.rm = TRUE),
                    stats::sd(get(value), na.rm = TRUE),
                    sum(.obs, na.rm = TRUE)),
               paste0(out_prefix, c("_mean","_sd","_n"))), by = keys]
  merge(p0, bb, by = keys, all = TRUE)
}

dir.create("output/summary", showWarnings = FALSE, recursive = TRUE)
summary_availability <- as.data.frame(boot_summary(
  avail_long, "pct", c(CELL_KEYS, "fip", "warming_scenario"), out_prefix = "Estimate"))
saveRDS(summary_availability, "output/summary/summary_availability.rds")

#-----------------------------------------------
# Cattle shifts (003 associations x availability)  ####
# cattleA: via total availability; cattleB/C: via own and neighbouring production.
cat_in <- merge(avail_long, assoc, by = "fip")                  # avail_long already carries avail00
cat_in <- merge(cat_in, base_avail[, .(fip, prod00_base, prod00_LM)], by = "fip")
cat_in[, `:=`(
  cattleA = pct * b_avail00,                                              # total availability
  cattleB = ((prod    / prod00_base) - 1) * 100 * b_prod00,               # own production
  cattleC = ((prodlag / prod00_LM)   - 1) * 100 * b_prod00_LM)]           # neighbour production
summary_cattle <- Reduce(function(a,b) merge(a, b, by = c(CELL_KEYS,"fip","warming_scenario"), all = TRUE),
  list(boot_summary(cat_in, "cattleA", c(CELL_KEYS,"fip","warming_scenario")),
       boot_summary(cat_in, "cattleB", c(CELL_KEYS,"fip","warming_scenario")),
       boot_summary(cat_in, "cattleC", c(CELL_KEYS,"fip","warming_scenario"))))
saveRDS(as.data.frame(summary_cattle), "output/summary/summary_cattle.rds")

message("Wrote availability_warming.rds (", nrow(avail_long), " rows), ",
        "summary_availability.rds and summary_cattle.rds.")
#-----------------------------------------------
