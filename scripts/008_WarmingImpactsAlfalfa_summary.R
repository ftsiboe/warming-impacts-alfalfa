#-----------------------------------------------
# Preliminaries                              ####
rm(list=ls(all=TRUE));gc();library(magrittr);library(future.apply);library(tidyverse);library(data.table)
library(plm);library(car);library(lmtest);library(terra);library(GWmodel);library(sp);library(sf)
study_environment <- readRDS("data/study_environment.rds")
invisible(lapply(list.files("scripts/helpers", pattern = "[.]R$", full.names = TRUE), source))

# ==============================================================================
# Summarise the per-cell CONSENSUS boot files written by 007.
# ------------------------------------------------------------------------------
# 007 now writes ONE consensus file per (boot x cell), named
#   consensus_<crop>_period<period>_<NAME>.rds
# The 20/50 GW specifications have already been reduced to a per-county consensus
# inside 007, so there are NO p/theta/longlat/DistName/kernel/specN columns here.
# The bootstrap summary keys are just the analysis cell (crop, period,
# climate_base) plus each block's own dimensions (fip / warming_scenario / name /
# Temp). For every quantity we report the full-sample point estimate (boot
# "0000") together with the across-boot mean / sd / n over the 100 draws.
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
# Alfalfa availability                       ####
availability <- read_block("availability")
setDT(availability)
# percentage change in availability under each warming scenario
for(sfx in c("05","10","15","20","25","30"))
  availability[[paste0("avail",sfx)]] <- ((availability[[paste0("avail",sfx)]]/availability$avail00)-1)*100
availability <- availability[, c("boot", CELL_KEYS, "fip",
                                 "avail00","avail05","avail10","avail15","avail20","avail25","avail30"), with = FALSE]
availability <- as.data.frame(availability) |>
  tidyr::gather(warming_scenario, Estimate, c("avail00","avail05","avail10","avail15","avail20","avail25","avail30"))
availability$warming_scenario <- as.numeric(gsub("[^0-9]","",availability$warming_scenario))/10
summary_availability <- as.data.frame(boot_summary(
  availability, "Estimate", c(CELL_KEYS, "fip","warming_scenario")))
saveRDS(summary_availability, "output/summary/summary_availability.rds")

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
# Cattle shifts                              ####
# consensus association coefficients (spread) x consensus availability
assoc <- as.data.frame(readRDS("output/availability_associations.rds")$associations)[, c("fip","name","est")]
assoc$name <- paste0("b_", assoc$name)
assoc <- assoc |> tidyr::spread(name, est)
setDT(assoc)

cattle <- read_block("availability")
setDT(cattle)
cattle <- cattle[assoc, on = "fip", nomatch = 0]   # associations are boot-invariant -> join on fip only

for(sfx in c("05","10","15","20","25","30")){
  a <- paste0("avail",sfx); p <- paste0("prod",sfx); pl <- paste0("prod",sfx,"_LM")
  cattle[[paste0("cattleA",sfx)]] <- ((cattle[[a]]/cattle$avail00)-1)*100*cattle$b_avail00
  cattle[[paste0("cattleB",sfx)]] <- ((cattle[[p]]/cattle$prod00)-1)*100*cattle$b_prod00
  cattle[[paste0("cattleC",sfx)]] <- ((cattle[[p]]/cattle$prod00)-1)*100*cattle$b_prod00 +
                                     ((cattle[[pl]]/cattle$prod00_LM)-1)*100*cattle$b_prod00_LM
}
cattle <- as.data.frame(cattle)[, c("boot", CELL_KEYS, "fip", names(cattle)[grepl("cattle",names(cattle))])]
cattle <- cattle |> tidyr::gather(warming_scenario, Estimate, names(cattle)[grepl("cattle",names(cattle))])
cattle$cattle <- gsub("[0-9]","",cattle$warming_scenario)
cattle$warming_scenario <- as.numeric(gsub("[^0-9]","",cattle$warming_scenario))/10
cattle <- cattle |> tidyr::spread(cattle, Estimate)

summary_cattle <- Reduce(function(a,b) merge(a,b,by=c(CELL_KEYS,"fip","warming_scenario"),all=TRUE),
  list(boot_summary(cattle, "cattleA", c(CELL_KEYS,"fip","warming_scenario")),
       boot_summary(cattle, "cattleB", c(CELL_KEYS,"fip","warming_scenario")),
       boot_summary(cattle, "cattleC", c(CELL_KEYS,"fip","warming_scenario"))))
saveRDS(as.data.frame(summary_cattle), "output/summary/summary_cattle.rds")

#-----------------------------------------------
# Bandwidths (diagnostic)                     ####
bw <- read_block("bw")
saveRDS(as.data.frame(bw), "output/summary/summary_bw.rds")
#-----------------------------------------------
