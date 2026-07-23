#-----------------------------------------------
# Preliminaries                              ####
rm(list=ls(all=TRUE));gc();library(tidyverse);library(data.table)
library(plm);library(car);library(lmtest);library(terra);library(sp);library(sf)
study_environment <- readRDS("data/study_environment.rds")
invisible(lapply(list.files("scripts/helpers", pattern = "[.]R$", full.names = TRUE), source))

# ==============================================================================
# Summarise the per-cell boot files written by 007 (yield-side quantities only).
# ------------------------------------------------------------------------------
# 007 writes ONE file per (boot x cell), named
#   consensus_cluster_<crop>_period<period>_<NAME>.rds
# carrying exposure / piecewise / relation / impact_yield (+ sp_data). No GW runs
# in 007 anymore, so there is no availability / cattle / bandwidth block here:
# warming-scenario availability and cattle shifts are summarised downstream by
# 009_..._availability_warming.R. The associations pass-through below is a
# boot-invariant copy of 003's output.
#
# The bootstrap summary keys are the analysis cell (crop, period, climate_base)
# plus each block's own dimensions (fip / warming_scenario / name / Temp). For
# every quantity we report the full-sample point estimate (boot "0000") together
# with the across-boot mean / sd / n over the 100 draws.
# ==============================================================================

CELL_KEYS <- c("crop","period","climate_base")

# read one block ($<name>) from every consensus file across all boot folders
read_block <- function(block){
  files <- list.files(study_environment$wd$boots, recursive = TRUE, full.names = TRUE,
                      pattern = "^consensus_.*\\.rds$")
  data.table::rbindlist(lapply(files, function(f){
    tryCatch(readRDS(f)[[block]], error = function(e) NULL)
  }), fill = TRUE)
}

# bootstrap summary: point estimate at boot "0000" + mean/sd/n over draws
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

#-----------------------------------------------
# Exposure                                   ####
exposure <- read_block("exposure")
summary_exposure <- as.data.frame(boot_summary(exposure, "exp", c(CELL_KEYS, "Temp")))
saveRDS(summary_exposure, "output/summary/summary_exposure.rds")

#-----------------------------------------------
# Piecewise linear function                  ####
piecewise <- read_block("piecewise")
summary_piecewise <- as.data.frame(boot_summary(piecewise, "Estimate", c(CELL_KEYS, "name")))
saveRDS(summary_piecewise, "output/summary/summary_piecewise.rds")

#-----------------------------------------------
# Nonlinear Relation                         ####
relation <- read_block("relation")
summary_relation <- as.data.frame(boot_summary(relation, "Piece", c(CELL_KEYS, "Temp")))
saveRDS(summary_relation, "output/summary/summary_relation.rds")

#-----------------------------------------------
# Yield impact                               ####
impact_yield <- read_block("impact_yield")
impact_yield <- impact_yield[, c("boot", CELL_KEYS, "region","state_code","county_code","fip",
                                 "cc05","cc10","cc15","cc20","cc25","cc30"), with = FALSE]
impact_yield <- impact_yield |>
  tidyr::gather(warming_scenario, Estimate, c("cc05","cc10","cc15","cc20","cc25","cc30"))
impact_yield$warming_scenario <- as.numeric(gsub("cc","",impact_yield$warming_scenario))/10
summary_impact_yield <- as.data.frame(boot_summary(
  impact_yield, "Estimate", c(CELL_KEYS, "region","state_code","county_code","fip","warming_scenario")))
saveRDS(summary_impact_yield, "output/summary/summary_impact_yield.rds")

#-----------------------------------------------
# Alfalfa availability + cattle shifts        ####
# Moved to 009_..._availability_warming.R: both need the GW neighbourhood lag,
# which is now applied once OUTSIDE the bootstrap. 007 no longer emits an
# `availability` block, so nothing is read here. 009 writes
# summary_availability.rds and summary_cattle.rds.

#-----------------------------------------------
# Associations  (boot-invariant -> from 003)  ####
# The cattle ~ availability associations no longer live in the per-boot payload;
# they are produced once by 003_..._availability_associations.R. They are
# boot-invariant, so we stamp the canonical cell keys and carry both a plain `est`
# column and the boot_summary-style `est_0000/_mean/_sd/_n` (0000 = mean = the
# estimate; sd = the across-spec dispersion) so downstream reads either schema.
assoc003 <- as.data.frame(readRDS("output/availability_associations.rds")$associations)
summary_associations <- data.frame(
  crop = "hay_alfalfa", period = 107L, climate_base = "1991_2020",
  fip = assoc003$fip, name = assoc003$name,
  est = assoc003$est, est_0000 = assoc003$est, est_mean = assoc003$est,
  est_sd = assoc003$est_specsd, est_n = NA_integer_,
  se = assoc003$se, tv = assoc003$tv, pv = assoc003$pv,
  sign_agreement = assoc003$sign_agreement,
  row.names = NULL, stringsAsFactors = FALSE)
saveRDS(summary_associations, "output/summary/summary_associations.rds")

#-----------------------------------------------
