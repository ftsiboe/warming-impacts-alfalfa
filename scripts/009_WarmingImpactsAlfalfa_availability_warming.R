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
# RESUMABLE: the expensive per-cell lag loop checkpoints each cell to its own file
# under output/checkpoints/availability_warming/. On restart, cells whose checkpoint
# already exists are skipped, so a broken run continues from where it stopped. The
# 007 boot read is cached too (output/checkpoints/imp.rds). To force a clean start,
# delete output/checkpoints/. Writes are atomic (temp file + rename) so a break
# mid-save never leaves a truncated checkpoint that would be mistaken for "done".
#
# OPTIONAL PARALLELISM: the cell loop runs sequentially by default. Turn on parallel
# execution over cells with either
#     options(aw.parallel = TRUE, aw.workers = 6)   # before sourcing, or
#     Sys.setenv(AW_PARALLEL = "1", AW_WORKERS = "6")
# It uses future.apply (multicore fork on Linux/Mac, multisession on Windows) and
# pins BLAS + data.table to one thread per worker so the workers don't oversubscribe
# the cores that gwkit's BLAS lag already uses. Checkpointing is unchanged: each
# worker skips finished cells and writes its own cell atomically, so resume and
# parallelism compose. `imp` is sliced to one small file per cell up front, and each
# worker reads only its own slice from disk - so the big table is never exported as a
# per-worker global (that is what tripped future.globals.maxSize) and RAM stays bounded
# to one cell per worker on both fork and multisession.
#
# Inputs : 003 availability_associations.rds (baseline avail00 + associations)
#          007 consensus_cluster_*.rds       (per-boot impact_yield)
#          data/spatial_representation.rds   (fixed prod00 = area x yield, cattle)
# Outputs: output/availability_warming.rds          (boot x cell x county x scenario)
#          output/summary/summary_availability.rds
#          output/summary/summary_cattle.rds
#-----------------------------------------------
rm(list = ls(all = TRUE)); gc()
library(tidyverse); library(data.table); library(terra)
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

options(
  aw.parallel = max(1L, round(parallel::detectCores()*0.5)) > 1, 
  aw.workers = max(1L, round(parallel::detectCores()*0.5)))

CELL_KEYS   <- c("crop","period","climate_base")
SFX         <- c("05","10","15","20","25","30")
metrics     <- gw_distance_metric_names()
kernels_all <- c("gaussian","exponential","bisquare","boxcar","tricube")

# checkpoint locations (created up front so both the imp cache and the cell loop can use them)
CKPT_ROOT <- "output/checkpoints"
CKPT_DIR  <- file.path(CKPT_ROOT, "availability_warming")   # per-cell OUTPUT checkpoints
IMP_DIR   <- file.path(CKPT_ROOT, "imp_by_cell")            # per-cell INPUT slices
dir.create(CKPT_DIR, showWarnings = FALSE, recursive = TRUE)
dir.create(IMP_DIR,  showWarnings = FALSE, recursive = TRUE)

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
# 007: per-boot county yield impacts (cached)  ####
# Cache the assembled `imp` once; restarts skip the (potentially slow) multi-file read.
imp_cache <- file.path(CKPT_ROOT, "imp.rds")
if (file.exists(imp_cache)) {
  imp <- readRDS(imp_cache)
  message("Loaded cached imp (", nrow(imp), " rows) from ", imp_cache)
} else {
  bfiles <- list.files(study_environment$wd$boots, recursive = TRUE, full.names = TRUE,
                       pattern = "^consensus_cluster_.*\\.rds$")
  if (length(bfiles) == 0) stop("No 007 boot files found under ", study_environment$wd$boots)
  imp <- data.table::rbindlist(lapply(bfiles, function(f)
    tryCatch(data.table::as.data.table(readRDS(f)$impact_yield), error = function(e) NULL)), fill = TRUE)
  imp <- imp[!fip %in% c("0","00000")]                       # county rows only
  imp[, fip := stringr::str_pad(as.character(fip), 5, pad = "0")]
  tmp <- paste0(imp_cache, ".tmp"); saveRDS(imp, tmp); file.rename(tmp, imp_cache)
}

#-----------------------------------------------
# Per cell: lag EVERY boot x scenario in one pass   ####
cells <- unique(imp[, ..CELL_KEYS])
message("Cells: ", nrow(cells), " ; boots: ", data.table::uniqueN(imp$boot))

# ---- stable, filename-safe tag + paths per cell ----
cell_tag       <- function(cl) gsub("[^A-Za-z0-9_.-]", "-", paste(cl$crop, cl$period, cl$climate_base, sep = "__"))
cell_file      <- function(cl) file.path(CKPT_DIR, paste0(cell_tag(cl), ".rds"))   # per-cell OUTPUT checkpoint
imp_slice_file <- function(cl) file.path(IMP_DIR,  paste0(cell_tag(cl), ".rds"))   # per-cell INPUT slice

# Split the big `imp` (can be >1 GiB) into one small file per cell, ONCE, then drop it
# from memory. Workers (and the sequential loop) read only their own cell's slice from
# disk - so `imp` is never exported as a per-worker global and RAM stays bounded to one
# cell at a time. Idempotent + atomic, so this step also resumes cleanly.
for (ci in seq_len(nrow(cells))) {
  cl <- cells[ci]; sf <- imp_slice_file(cl)
  if (file.exists(sf)) next
  d   <- imp[crop %in% cl$crop & period %in% cl$period & climate_base %in% cl$climate_base]
  tmp <- paste0(sf, ".tmp"); saveRDS(d, tmp); file.rename(tmp, sf)
}
rm(imp); gc()

# ---- one cell's work, factored out so the loop can checkpoint around it ----
# Returns a data.frame for the cell, or NULL when the cell legitimately has no data.
compute_cell <- function(cl){
  d <- data.table::as.data.table(readRDS(imp_slice_file(cl)))
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
}

# ---- process one cell index: skip if done, else compute + atomic checkpoint ----
# Safe to call from a parallel worker: each cell owns its own file, and the
# file.exists() guard + temp-then-rename write make it idempotent under a break.
process_cell <- function(ci){
  cl <- cells[ci]; f <- cell_file(cl)
  if (file.exists(f)) return(paste0("skip ", ci, "/", nrow(cells), " (already done)"))
  res <- tryCatch(compute_cell(cl), error = function(e) e)
  if (inherits(res, "error"))               # leave uncheckpointed -> retried on next run
    return(paste0("ERROR ", ci, "/", nrow(cells), ": ", conditionMessage(res)))
  tmp <- paste0(f, ".tmp"); saveRDS(res, tmp); file.rename(tmp, f)
  paste0("done ", ci, "/", nrow(cells))
}

# ---- optional parallelism toggle (default: sequential) ----
PARALLEL  <- isTRUE(getOption("aw.parallel",
                              tolower(Sys.getenv("AW_PARALLEL", "false")) %in% c("1","true","yes")))
N_WORKERS <- as.integer(getOption("aw.workers",
                                  Sys.getenv("AW_WORKERS", max(1L, parallel::detectCores() - 1L))))

todo <- which(!file.exists(vapply(seq_len(nrow(cells)), function(ci) cell_file(cells[ci]), character(1))))
message(length(todo), " of ", nrow(cells), " cells to compute",
        if (PARALLEL) paste0(" (parallel: ", N_WORKERS, " workers)") else " (sequential)")

if (PARALLEL && length(todo) > 1L &&
    requireNamespace("future", quietly = TRUE) &&
    requireNamespace("future.apply", quietly = TRUE)) {
  # pin math libraries to 1 thread per worker so workers don't oversubscribe the
  # cores gwkit's BLAS lag already uses
  if (requireNamespace("RhpcBLASctl", quietly = TRUE)) try(RhpcBLASctl::blas_set_num_threads(1), silent = TRUE)
  old_dt <- data.table::getDTthreads(); data.table::setDTthreads(1)
  old_max <- getOption("future.globals.maxSize")            # slices keep globals tiny, but leave headroom
  options(future.globals.maxSize = 512 * 1024^2)
  plan_kind <- if (grepl("windows", sysname)) "multisession" else "multicore"
  future::plan(plan_kind, workers = N_WORKERS)
  on.exit({ future::plan("sequential"); data.table::setDTthreads(old_dt)
            options(future.globals.maxSize = old_max) }, add = TRUE)

  invisible(future.apply::future_lapply(todo, function(ci){
    data.table::setDTthreads(1)             # re-assert inside each worker
    message(process_cell(ci)); NULL
  }, future.seed = TRUE))                    # reproducible RNG across workers
} else {
  if (PARALLEL) message("future/future.apply not available (or nothing to do) - running sequentially.")
  for (ci in todo) message(process_cell(ci))
}

# ---- reassemble from checkpoints ----
done <- list.files(CKPT_DIR, pattern = "\\.rds$", full.names = TRUE)
if (length(done) < nrow(cells))
  warning(nrow(cells) - length(done), " cell(s) still incomplete - re-run to finish.")

avail_long <- data.table::rbindlist(lapply(done, readRDS), fill = TRUE)  # NULL cells are dropped
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
