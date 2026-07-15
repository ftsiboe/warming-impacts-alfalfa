#-----------------------------------------------
# Preliminaries                              ####
rm(list=ls(all=TRUE));gc();library(magrittr);library(future.apply);library(tidyverse);library(data.table)
library(plm);library(car);library(lmtest);library(terra);library(GWmodel);library(sp);library(sf)
study_environment <- readRDS("data/study_environment.rds")
invisible(lapply(list.files("scripts/helpers", pattern = "[.]R$", full.names = TRUE), source))
# devtools::document(file.path(dirname(dirname(getwd())),"packages/rAgroClimate"))

# gwkit (sibling package) is the canonical GW toolkit for this project.
sysname   <- tolower(as.character(Sys.info()[["sysname"]]))
gwkit_src <- if (grepl("windows", sysname)) {
  file.path(dirname(dirname(getwd())), "packages/gwkit")
} else {
  file.path(dirname(getwd()), "packages/gwkit")
}

# INSTALL the local gwkit so parallel WORKERS load the same version as the main
# session. devtools::load_all() injects gwkit into this process only; multisession
# workers fall back to the installed copy (multicore/fork DOES inherit load_all, so
# on the cluster this is belt-and-braces). Install when gwkit is missing, older than
# the source, OR predates the consensus refactor (dev version may be unchanged, so
# also probe the API). For the SLURM array, install gwkit ONCE before submitting so
# all tasks find it current and skip this (avoids a concurrent-install race).
.src_ver  <- tryCatch(as.package_version(read.dcf(file.path(gwkit_src, "DESCRIPTION"))[, "Version"]),
                      error = function(e) NULL)
.inst_ver <- tryCatch(utils::packageVersion("gwkit"), error = function(e) NULL)
.has_api  <- tryCatch(
  "gw_consensus_scalar" %in% getNamespaceExports("gwkit") &&
    "terms" %in% names(formals(getExportedValue("gwkit", "estimate_gwr"))),
  error = function(e) FALSE)
if (is.null(.inst_ver) || (!is.null(.src_ver) && .inst_ver < .src_ver) || !.has_api) {
  message("Installing local gwkit from ", gwkit_src, " ...")
  devtools::install(gwkit_src, upgrade = FALSE, quick = TRUE, quiet = TRUE)
}
library(gwkit)
set.seed(08032024)

# ==============================================================================
# STRUCTURE (CLUSTER variant: county knots+slopes from 005 optimal_knots_cluster)
# ------------------------------------------------------------------------------
# The SLURM array is BY BOOTSTRAP: one array task = one boot ("0000" + 100 draws).
# Within a task, for every analysis cell (crop x window x climate baseline):
#   * the spec-INVARIANT yield model is computed ONCE - the national piecewise
#     fit, plus per-county KNOTS + boot-refit SLOPES from the CLUSTER estimator
#     (005). County impacts use each county's own knot & cluster slope, built
#     from prism_climate at the union of thresholds 006 retains; aggregates are
#     acreage-weighted; then
#   * the GW AVAILABILITY/ASSOCIATIONS quantities (cattle ~ alfalfa availability)
#     are estimated under ALL 50 gwkit distance-metric presets x kernels (10 x 5)
#     in parallel across cores and reduced to a per-county CONSENSUS in-task.
# Only the consensus is written (one file per cell per boot).
#
# gwkit tools used (boot):
#   gw_distance_metric_names()               -> the 10 distance presets
#   estimate_gwlag()                         -> WARMING-scenario availability m(z_i)
#                                               (neighbour-weighted, self-excluded)
# Baseline avail00 + the cattle~availability associations (estimate_gwr) are
# produced once by 003_..._availability_associations.R and are NOT in this loop.
#
# The optimal_gw cross-validation pick is REMOVED (consensus replaces it).
# Bandwidths are precomputed once per (cell x spec) on the full sample and cached
# in output/gw_bandwidths.rds, then reused across all boots.
# ==============================================================================

#-----------------------------------------------
# boot list                                  #### (run once)
function(){
  boot_list <- readRDS("data/nass_hay_production.rds")
  boot_list <- list("0000"=unique(boot_list$commodity_year))
  for(bt in 1:100){
    boot_list[[stringr::str_pad(bt,pad="0",4)]] <- sample(boot_list$`0000`, length(boot_list$`0000`), replace = TRUE)
  }
  saveRDS(boot_list,file="output/boot_list.rds")
}

#-----------------------------------------------
# GW specifications: gwkit distance presets x kernels (10 x 5 = 50)  ####
spec_gw <- as.data.frame(data.table::rbindlist(
  lapply(gw_distance_metric_names(), function(dm){
    data.frame(distance_metric = dm,
               kernel = c("gaussian","exponential","bisquare","boxcar","tricube"),
               stringsAsFactors = FALSE)
  }), fill = TRUE))
spec_gw$specN <- 1:nrow(spec_gw)

#-----------------------------------------------
# Analysis cells (crop x window x climate baseline)         ####
CELLS <- data.frame(crop="hay_alfalfa",period=107,NAME=c("1991_2020","1981_2010","1971_2000","1961_1990","1951_2022"))
CELLS$climate_base[1] <- list(1991:2020) # NOAA [https://www.ncei.noaa.gov/products/land-based-station/us-climate-normals]
CELLS$climate_base[2] <- list(1981:2010) # NOAA [https://journals.ametsoc.org/view/journals/bams/93/11/bams-d-11-00197.1.xml]
CELLS$climate_base[3] <- list(1971:2000) # NOAA
CELLS$climate_base[4] <- list(1961:1990) # IPCC
CELLS$climate_base[5] <- list(1951:2022)

CELLS_crop <- data.frame(crop=c("hay_other","hay_all"),period=107,NAME="1991_2020")
CELLS_crop$climate_base[1] <- list(1991:2020)
CELLS_crop$climate_base[2] <- list(1991:2020)

CELLS_period <- data.frame(crop=c("hay_alfalfa"),period=c(0:12,101:106,108:112),NAME="1991_2020")
for(iii in 1:nrow(CELLS_period)){
  CELLS_period$climate_base[iii] <- list(1991:2020)
}

CELLS <- unique(rbind(CELLS,CELLS_crop,CELLS_period))
CELLS$cellN <- 1:nrow(CELLS)
rm(CELLS_crop,CELLS_period)

#-----------------------------------------------
# County degree-day yield slopes from the CLUSTER estimator            ####
# Per-county b_DD1/b_DD2/b_DD3 at each county's CLUSTER knot (from 005's
# optimal_knots_cluster). Knots are FIXED (from 005); only the SLOPES are refit
# per boot: for each cluster, pool its counties on the bootstrap sample and fit
# the piecewise FE panel lny ~ DD1+DD2+DD3+ppt+ppt2 + state trends at the cluster
# knot, then attach the cluster's slopes to its member counties. Returns per-county
# fip/Tmin/Tmax/b_DD1-3, or NULL (=> national PWM.COEF fallback) when no cluster
# knots exist for `crop`. No GW here - a handful of FE fits per boot.
.cluster_county_betas <- function(boot_panel, crop){
  cf <- "output/optimal_knots_cluster.rds"
  if(!file.exists(cf)) return(NULL)
  ck <- data.table::as.data.table(readRDS(cf))
  if("crop" %in% names(ck)) ck <- ck[crop %in% ..crop]
  ck[, county_fips := stringr::str_pad(as.character(fip), 5, pad = "0")]
  ck <- unique(ck[, .(county_fips, cluster, Tmin = as.integer(Tmin), Tmax = as.integer(Tmax))])
  ck <- ck[is.finite(Tmin) & is.finite(Tmax)]
  if(nrow(ck) == 0) return(NULL)

  pad    <- function(x) stringr::str_pad(x, 2, pad = "0")
  dd_all <- grep("^dday[0-9]+$", names(boot_panel), value = TRUE)
  d0 <- data.table::as.data.table(boot_panel)
  d0[, county_fips := stringr::str_pad(as.character(county_fips), 5, pad = "0")]
  d0 <- d0[, lapply(.SD, mean, na.rm = TRUE), by = c("state_code","county_fips","year"),
           .SDcols = c("lny","ppt","ppt2","Trend","Trend2", dd_all)]
  d0 <- merge(d0, ck[, .(county_fips, cluster, Tmin, Tmax)], by = "county_fips")
  if(nrow(d0) == 0) return(NULL)

  # refit each cluster's piecewise FE at its FIXED knot on the boot sample
  betas <- data.table::rbindlist(lapply(sort(unique(d0$cluster)), function(cc){
    tryCatch({
      dc <- data.table::copy(d0[cluster == cc])
      tmin <- dc$Tmin[1]; tmax <- dc$Tmax[1]
      dc[, DD1 := pmax(get("dday00")                  - get(paste0("dday", pad(tmin))), 0)]
      dc[, DD2 := pmax(get(paste0("dday", pad(tmin))) - get(paste0("dday", pad(tmax))), 0)]
      dc[, DD3 := pmax(get(paste0("dday", pad(tmax))), 0)]
      for(st in unique(dc$state_code)){
        dc[, (paste0("trend1", st)) := ifelse(state_code %in% st, Trend,  0)]
        dc[, (paste0("trend2", st)) := ifelse(state_code %in% st, Trend2, 0)]
      }
      rhs  <- c("ppt","ppt2","DD1","DD2","DD3", grep("^trend", names(dc), value = TRUE))
      form <- stats::as.formula(paste("lny ~", paste(rhs, collapse = " + ")))
      fit  <- plm::plm(form, data = as.data.frame(dc), index = c("county_fips","year"), model = "within")
      co   <- coef(fit)
      data.table::data.table(cluster = cc, b_DD1 = co[["DD1"]], b_DD2 = co[["DD2"]], b_DD3 = co[["DD3"]])
    }, error = function(e) NULL)
  }), fill = TRUE)
  if(nrow(betas) == 0) return(NULL)

  out <- merge(ck[, .(fip = county_fips, cluster, Tmin, Tmax)], betas, by = "cluster")
  out[, .(fip, Tmin, Tmax, b_DD1, b_DD2, b_DD3)]
}

#-----------------------------------------------
# BASE (spec-invariant) block for one boot x cell           ####
# Yield model + county production panel + projected county centroids. None of
# this depends on the GW specification, so it is computed once per cell.
estimate_cell_base <- function(boot, cell){
  tryCatch({
    period       <- cell$period
    crop         <- cell$crop
    climate_base <- cell$climate_base[[1]]

    data <- build_hay_weather_panel(
      crop = crop, target_periods = period,
      prism_weather_directory = "data/prism_weather")

    sampled_years <- readRDS("output/boot_list.rds")[[boot]]
    data <- as.data.frame(data.table::rbindlist(lapply(1:length(sampled_years), function(yr){
      tryCatch({ data <- data[data$commodity_year %in% sampled_years[yr],]; data$year <- yr; data },
               error = function(e){NULL}) }), fill = TRUE))

    # Exposure
    Expos <- data; Expos$id <- 1:nrow(Expos)
    Expos <- Expos[c("id",names(Expos)[grepl("exp",names(Expos))])]
    Expos <- Expos %>% tidyr::gather(Temp, exp,2:ncol(Expos))
    Expos$Temp <- as.numeric(gsub("[^0-9]","",gsub("exp","",Expos$Temp)))
    Expos <- doBy::summaryBy(exp~Temp,FUN=mean,keep.names=T,data=Expos,na.rm=T)

    # Piecewise linear
    knots <- readRDS("output/optimal_knots.rds")
    knots <- knots[(knots$target_periods %in% c(104:112,0) & knots$crop %in% crop),]
    knots <- knots[order(knots$cv_error),]
    Tmin <- knots$Tmin[1]; Tmax <- knots$Tmax[1]
    data$DD1 <- data[,"dday00"] - data[,paste0("dday",stringr::str_pad(Tmin,pad="0",2))]
    data$DD2 <- data[,paste0("dday",stringr::str_pad(Tmin,pad="0",2))] - data[,paste0("dday",stringr::str_pad(Tmax,pad="0",2))]
    data$DD3 <- data[,paste0("dday",stringr::str_pad(Tmax,pad="0",2))]
    data$DD1 <- ifelse(data$DD1<0,0,data$DD1); data$DD2 <- ifelse(data$DD2<0,0,data$DD2); data$DD3 <- ifelse(data$DD3<0,0,data$DD3)

    fit_data <- doBy::summaryBy(lny + ppt + ppt2 + ppt3 + Trend + Trend2 + DD1 + DD2 + DD3 + freeze~state_code+fip+commodity_year+year,
                                data=data,FUN=mean,keep.names = T,na.rm=T)
    for(st in unique(fit_data$state_code)){
      fit_data[,paste0("trend1",st)] <- ifelse(fit_data$state_code %in% st,fit_data$Trend,0)
      fit_data[,paste0("trend2",st)] <- ifelse(fit_data$state_code %in% st,fit_data$Trend2,0)
    }
    fit_data$panel <- fit_data$fip; fit_data$time <- fit_data$year
    fit_data<- pdata.frame(fit_data, index = c("panel", "time"), drop.index = TRUE)
    fit <- plm(as.formula(paste0("lny ~",paste0(c("ppt","ppt2","DD1","DD2","DD3",
                                                  names(fit_data)[grepl("trend",names(fit_data))]),collapse = "+"))),
               data = fit_data,model = "within")
    PWM.COEF <- coef(fit)
    G <- length(unique(fit_data$year))
    PWM.VCOV <- G/(G - 1) * vcovHC(fit, type = "HC1", cluster = "time")
    test1 <- linearHypothesis(fit, c("DD1=0","DD2=0","DD3=0"),vcov. = PWM.VCOV)
    test2 <- linearHypothesis(fit, c("ppt=0","ppt2=0"),vcov. = PWM.VCOV)
    test3 <- linearHypothesis(fit, c("DD1=0","DD2=0","DD3=0","ppt=0","ppt2=0"),vcov. = PWM.VCOV)
    test4 <- plmtest(fit$formula, data=fit_data, effect="individual", type="honda",vcov. = PWM.VCOV)
    PWM.TAB <- as.data.frame(coeftest(fit,G/(G - 1) * vcovHC(fit, type = "HC1", cluster = "time"))[,])
    PWM.TAB$name <- rownames(PWM.TAB); names(PWM.TAB) <- c("Estimate","StdError","t_value","p_value","name")
    fit <- summary(fit)
    PWM.TAB <- rbind(PWM.TAB, data.frame(
      name = c("n","n_farms","n_year","r.squared","test_temp","test_ppt","test_weather","test_fe"),
      Estimate=c(nrow(fit$model),length(unique(data$fip)),length(unique(data$year)),fit$r.squared[1],
                 test1$`Chisq`[2],test2$`Chisq`[2],test3$`Chisq`[2],test4$statistic),
      StdError=NA,t_value=NA,
      p_value=c(0,0,0,0,test1$`Pr(>Chisq)`[2],test2$`Pr(>Chisq)`[2],test3$`Pr(>Chisq)`[2],test4$p.value)))
    rm(test1,test2,test3,test4,fit,G,fit_data);gc()

    # Nonlinear relation
    Relation <- data.frame(Temp=1:45)
    Relation$I1<-as.numeric(Relation$Temp>=0); Relation$I2<-as.numeric(Relation$Temp>=Tmin); Relation$I3<-as.numeric(Relation$Temp>=Tmax)
    funcx <- paste0("~(1-",Relation$I2,")*(",Relation$Temp,"-0)*DD1+",
                    Relation$I2,"*(1-",Relation$I3,")*(DD1*",Tmin," + DD2*(",Relation$Temp,"-",Tmin,")) +",
                    Relation$I3,"*(DD1*",Tmin," + DD2*(",Tmax,"-",Tmin,") + DD3*(",Relation$Temp,"-",Tmax,"))")
    Form <- as.formula(funcx[1]); for(i in 2:length(funcx)){Form<-c(Form,as.formula(funcx[i]))}
    Relation <- compute_delta_method(func=Form,vcMat=PWM.VCOV,coefs=PWM.COEF)[c("Estimate","SE")]
    names(Relation) <- c("Piece","PieceSE"); Relation <- as.data.frame(Relation); Relation$Temp <- 1:45
    rm(Form);gc()

    # County polygons + centroids (spec-invariant; also needed by the GWFE block)
    USMUR    <- rast(file.path(study_environment$gssurgo_archive,"MURASTER_30m.tif"))
    Counties <- vect(file.path(study_environment$usaPolygons_archive,"USA_Counties.shp"))
    Counties <- crop(project(Counties, crs(USMUR)), ext(USMUR)); rm(USMUR);gc()
    Counties$fip <- as.character(paste0(stringr::str_pad(as.numeric(as.character(Counties$STATEFP)),2,pad="0"),
                                        stringr::str_pad(as.numeric(as.character(Counties$COUNTYFP)),3,pad="0")))
    cxy <- geom(centroids(Counties))[, c("x","y")]
    centroids_all <- data.frame(fip = Counties$fip, longitude = cxy[,"x"], latitude = cxy[,"y"], stringsAsFactors = FALSE)

    # County knots + boot-refit slopes from the CLUSTER estimator (005). Counties
    # with no cluster fall back to the national knot + national PWM.COEF.
    Tmin_nat <- as.integer(Tmin); Tmax_nat <- as.integer(Tmax)
    beta_county <- tryCatch(.cluster_county_betas(data, crop), error = function(e) NULL)

    # Warming-scenario climate exposure at EACH county's own (cluster) knot
    sim_climate <- data.table::rbindlist(lapply(list.files("data/prism_climate", full.names = TRUE),
      function(file){ tryCatch(readRDS(file), error = function(e) NULL) }), fill = TRUE); gc()
    per_sel <- period                       # scalar; avoid data.table column/var name clash
    sim_climate <- sim_climate[period %in% per_sel]
    sim_climate <- as.data.frame(data.table::rbindlist(lapply(sampled_years, function(nm){
      sim_climate[sim_climate$commodity_year %in% nm, ] }), fill = TRUE)); gc()
    sim_climate <- sim_climate[sim_climate$commodity_year %in% climate_base, ]
    data.table::setDT(sim_climate)
    sim_climate[, county_fips := stringr::str_pad(as.character(county_fips), 5, pad = "0")]
    sim_climate[, warming_scenario := as.numeric(as.character(warming_scenario))]

    # per-county knot + slope table (national fallback where a county has no cluster)
    pad <- function(x) stringr::str_pad(x, 2, pad = "0")
    kb  <- data.table::data.table(county_fips = unique(sim_climate$county_fips))
    if(!is.null(beta_county)){
      bc <- data.table::as.data.table(beta_county)
      bc[, county_fips := stringr::str_pad(as.character(fip), 5, pad = "0")]
      kb <- merge(kb, bc[, .(county_fips, Tmin, Tmax, b_DD1, b_DD2, b_DD3)], by = "county_fips", all.x = TRUE)
    } else {
      kb[, `:=`(Tmin = NA_integer_, Tmax = NA_integer_, b_DD1 = NA_real_, b_DD2 = NA_real_, b_DD3 = NA_real_)]
    }
    kb[is.na(Tmin) | is.na(Tmax), `:=`(Tmin = Tmin_nat, Tmax = Tmax_nat,
                                       b_DD1 = PWM.COEF["DD1"], b_DD2 = PWM.COEF["DD2"], b_DD3 = PWM.COEF["DD3"])]
    sim_climate <- merge(sim_climate, kb, by = "county_fips")

    # degree-day exposure at each county's knot (loop the few distinct knot pairs)
    prs <- unique(kb[, .(Tmin, Tmax)])
    sim_climate <- data.table::rbindlist(lapply(seq_len(nrow(prs)), function(pp){
      tmin <- prs$Tmin[pp]; tmax <- prs$Tmax[pp]
      d <- sim_climate[Tmin == tmin & Tmax == tmax]
      if(nrow(d) == 0) return(NULL)
      d[, DD1 := pmax(get("dday00")                  - get(paste0("dday", pad(tmin))), 0)]
      d[, DD2 := pmax(get(paste0("dday", pad(tmin))) - get(paste0("dday", pad(tmax))), 0)]
      d[, DD3 := pmax(get(paste0("dday", pad(tmax))), 0)]
      d
    }), fill = TRUE)

    # mean exposure by county x scenario, baseline (warming 0) delta, county impact
    sim_climate <- sim_climate[, .(DD1 = mean(DD1, na.rm=TRUE), DD2 = mean(DD2, na.rm=TRUE), DD3 = mean(DD3, na.rm=TRUE),
                                   b_DD1 = b_DD1[1], b_DD2 = b_DD2[1], b_DD3 = b_DD3[1]),
                               by = .(county_fips, warming_scenario)]
    base0 <- sim_climate[warming_scenario == 0, .(county_fips, DD1b = DD1, DD2b = DD2, DD3b = DD3)]
    sim_climate <- base0[sim_climate, on = "county_fips", nomatch = 0]
    sim_climate <- sim_climate[!warming_scenario %in% 0]
    sim_climate[, impact := (exp((DD1-DD1b)*b_DD1 + (DD2-DD2b)*b_DD2 + (DD3-DD3b)*b_DD3) - 1) * 100]

    # region + state/county codes
    sim_climate[, `:=`(state_code = as.numeric(substr(county_fips, 1, 2)),
                       county_code = as.numeric(substr(county_fips, 3, 5)))]
    us_states <- as.data.frame(urbnmapr::get_urbn_map(map = "states", sf = TRUE))
    us_states$state_code <- as.numeric(as.character(us_states$state_fips))
    sim_climate <- merge(sim_climate, data.table::as.data.table(us_states)[, .(state_code, state_abbv)],
                         by = "state_code", all.x = TRUE)
    sim_climate[, region := NA_character_]
    sim_climate[state_abbv %in% c("CT","DE","ME","MD","MA","NH","NJ","NY","PA","RI","VT"), region := "Northeast"]
    sim_climate[state_abbv %in% c("IL","IN","IA","KS","MI","MN","MO","NE","ND","OH","SD","WI"), region := "Midwest"]
    sim_climate[state_abbv %in% c("AL","AR","FL","GA","KY","LA","MS","NC","OK","SC","TN","TX","VA","WV"), region := "South"]
    sim_climate[state_abbv %in% c("AZ","CA","CO","HI","ID","MT","NV","NM","OR","UT","WA","WY"), region := "West"]

    # county alfalfa acreage (production-area weight for the aggregate rows)
    prod <- as.data.frame(readRDS("data/nass_hay_production.rds"))
    prod <- prod[prod$commodity_name %in% crop, ]
    prod <- as.data.frame(data.table::rbindlist(lapply(sampled_years, function(nm){ prod[prod$commodity_year %in% nm, ] }), fill = TRUE))
    prod <- prod[prod$commodity_year %in% climate_base, ]
    prod <- doBy::summaryBy(area ~ state_code + county_code, data = prod, FUN = mean, keep.names = TRUE, na.rm = TRUE)
    data.table::setDT(prod); data.table::setnames(prod, "area", "acre")
    sim_climate <- merge(sim_climate, prod, by = c("state_code","county_code"), all.x = TRUE)
    sim_climate[!is.finite(acre), acre := 0]

    # county rows + ACREAGE-WEIGHTED aggregates (state, region, national)
    mkfip <- function(s, c) as.character(paste0(stringr::str_pad(s, 2, pad="0"), stringr::str_pad(c, 3, pad="0")))
    ci  <- sim_climate[, .(region, state_code, county_code, warming_scenario, impact)]
    st  <- sim_climate[, .(impact = stats::weighted.mean(impact, acre, na.rm=TRUE)),
                       by = .(region, state_code, warming_scenario)][, county_code := 0]
    rg  <- sim_climate[, .(impact = stats::weighted.mean(impact, acre, na.rm=TRUE)),
                       by = .(region, warming_scenario)][, `:=`(state_code = 0, county_code = 0)]
    nat <- sim_climate[, .(impact = stats::weighted.mean(impact, acre, na.rm=TRUE)),
                       by = .(warming_scenario)][, `:=`(region = "", state_code = 0, county_code = 0)]
    impact_yield <- data.table::rbindlist(list(ci, st, rg, nat), fill = TRUE, use.names = TRUE)
    impact_yield[, fip := mkfip(state_code, county_code)]
    impact_yield <- as.data.frame(impact_yield)
    impact_yield$warming_scenario <- paste0("cc",stringr::str_pad(as.numeric(as.character(round(impact_yield$warming_scenario*10))),2,pad="0"))
    impact_yield <- impact_yield[c("region","state_code","county_code","fip","warming_scenario","impact")]
    impact_yield <- impact_yield %>% tidyr::spread(warming_scenario, impact)

    # County production panel (Counties + centroids already built above)
    sp_data <- as.data.frame(readRDS("data/spatial_representation.rds"))
    sp_data <- dplyr::inner_join(impact_yield,sp_data,by=c("state_code","county_code","fip"))
    names(sp_data)[names(sp_data) %in% "inventory"] <- "cattle"
    sp_data <- dplyr::inner_join(sp_data, centroids_all, by = "fip")
    sp_data$prod00 <- sp_data$area*sp_data$yield
    sp_data$prod00 <- ifelse(sp_data$prod00 %in% c(NA,Inf,NaN,-Inf),0,sp_data$prod00)
    for(sfx in c("05","10","15","20","25","30")){
      cc <- paste0("cc",sfx)
      sp_data[[paste0("prod",sfx)]] <- tryCatch(
        ifelse(sp_data$prod00 %in% c(NA,Inf,NaN,-Inf,0),0,sp_data$area*sp_data$yield*(1+(sp_data[[cc]]/100))),
        error = function(e) NA_real_)
    }

    list(exposure=as.data.frame(Expos), piecewise=as.data.frame(PWM.TAB),
         relation=as.data.frame(Relation), impact_yield=as.data.frame(impact_yield),
         sp_data=as.data.frame(sp_data), centroids=centroids_all)
  }, error = function(e){return(NULL)})
}

#-----------------------------------------------
# GW block for one spec: WARMING-scenario availability via gwkit   ####
# (baseline avail00 + the cattle~availability associations now live in 003)
PROD_SFX  <- c("00","05","10","15","20","25","30")
PROD_COLS <- paste0("prod", PROD_SFX)

estimate_cell_gw <- function(base, spec, bw_lag = NULL, bw_gwr = NULL){
  tryCatch({
    sp_data   <- base$sp_data
    centroids <- base$centroids

    # --- Availability: neighbour-weighted (self-excluded) production, gwkit ---
    # estimate_gwlag uses a FIXED-distance bandwidth (bw_lag); estimate_gwr below
    # uses an ADAPTIVE kNN bandwidth (bw_gwr). They are different units, so each
    # estimator selects/reuses its OWN bandwidth - they are never shared.
    src <- sp_data[, c("fip","longitude","latitude", PROD_COLS)]
    lag <- gwkit::estimate_gwlag(
      data = src, unit = "fip", value_cols = PROD_COLS,
      coords = c("longitude","latitude"), predict = centroids,
      distance_metric = spec$distance_metric, kernel = spec$kernel,
      adaptive = FALSE, bw = bw_lag, bw_response = "prod00", include_self = FALSE)
    if(is.null(bw_lag)) bw_lag <- attr(lag, "bandwidth")
    lag <- as.data.frame(lag)

    av <- dplyr::left_join(as.data.frame(centroids["fip"]),
                           sp_data[, c("fip", PROD_COLS)], by = "fip")
    av <- dplyr::inner_join(av, lag, by = "fip")
    setDT(av)
    for(sfx in PROD_SFX){
      p <- paste0("prod",sfx); lm_ <- paste0("prod",sfx,"_LM")
      av[[paste0("avail",sfx)]] <- ifelse(is.na(av[[p]]),0,av[[p]]) + av[[lm_]]
    }
    av <- av[is.finite(avail00)]
    keep <- c("fip", PROD_COLS, paste0(PROD_COLS,"_LM"), paste0("avail",PROD_SFX))
    availability <- as.data.frame(av)[, keep]
    # national mean row
    nat <- as.data.frame(as.list(colMeans(availability[,-1], na.rm = TRUE)))
    nat$fip <- "00000"; availability <- rbind(nat[, names(availability)], availability)

    # Associations (cattle ~ availability) are boot-INVARIANT and are produced once
    # by 003_..._availability_associations.R - not recomputed per boot here.
    list(specN = spec$specN, availability = availability, bw_lag = bw_lag, bw_gwr = NA_real_)
  }, error = function(e){return(NULL)})
}

#-----------------------------------------------
# Across-spec consensus (median primary + mean retained)    ####
reduce_consensus <- function(base, spec_list, cell, boot, Counties = NULL){
  spec_list <- Filter(Negate(is.null), spec_list)
  if(length(spec_list) == 0) return(NULL)
  tag <- data.frame(boot = boot, crop = cell$crop, period = cell$period,
                    climate_base = cell$NAME, cellN = cell$cellN,
                    n_spec = length(spec_list), stringsAsFactors = FALSE)

  exposure     <- data.frame(tag, base$exposure,     row.names = NULL)
  piecewise    <- data.frame(tag, base$piecewise,    row.names = NULL)
  relation     <- data.frame(tag, base$relation,     row.names = NULL)
  impact_yield <- data.frame(tag, base$impact_yield, row.names = NULL)

  # availability consensus
  av <- data.table::rbindlist(lapply(seq_along(spec_list), function(i){
    d <- data.table::as.data.table(spec_list[[i]]$availability); d[, .specN := i]; d }), fill = TRUE)
  av_cols <- setdiff(names(av), c("fip",".specN"))
  av_med  <- av[, lapply(.SD, stats::median, na.rm=TRUE), by=fip, .SDcols=av_cols]
  av_mn   <- av[, lapply(.SD, mean, na.rm=TRUE), by=fip, .SDcols=av_cols]
  data.table::setnames(av_mn, av_cols, paste0(av_cols,"_specmean"))
  availability <- data.frame(tag, as.data.frame(merge(av_med, av_mn, by="fip")), row.names = NULL)

  # Associations are produced once by 003 (boot-invariant); not recomputed here.

  list(exposure=exposure, piecewise=piecewise, relation=relation, impact_yield=impact_yield,
       availability=availability,
       bw = data.frame(tag, specN = vapply(spec_list, function(z) z$specN, numeric(1)),
                       bw_lag = vapply(spec_list, function(z) z$bw_lag, numeric(1)), row.names = NULL))
}

#-----------------------------------------------
# Bandwidth cache (run once, before the array)              ####
# Precomputes bw per (cell x spec) on the FULL sample (boot "0000") and caches
# them in output/gw_bandwidths.rds for reuse across all boots.
function(){
  nw <- max(1, as.numeric(Sys.getenv("SLURM_CPUS_PER_TASK", unset = "1")))
  if(nw > 1) future::plan(future::multicore, workers = nw) else future::plan(future::sequential)
  bw_cache <- list()
  for(ci in 1:nrow(CELLS)){
    base <- estimate_cell_base(boot = "0000", cell = CELLS[ci, ])
    if(is.null(base)){
      bw_cache[[paste0("cell",CELLS$cellN[ci])]] <-
        data.frame(specN = spec_gw$specN, bw_lag = NA_real_, bw_gwr = NA_real_)
      next
    }
    bws <- future.apply::future_lapply(1:nrow(spec_gw), function(s){
      r <- estimate_cell_gw(base, spec_gw[s, ], bw_lag = NULL, bw_gwr = NULL)
      if(is.null(r)) c(NA_real_, NA_real_) else c(r$bw_lag, r$bw_gwr)
    }, future.seed = TRUE)
    bw_cache[[paste0("cell",CELLS$cellN[ci])]] <- data.frame(
      specN  = spec_gw$specN,
      bw_lag = vapply(bws, function(z) as.numeric(z[[1]]), numeric(1)),
      bw_gwr = vapply(bws, function(z) as.numeric(z[[2]]), numeric(1)))
    message("bandwidths cached: cell ", ci, "/", nrow(CELLS))
  }
  saveRDS(bw_cache, file = "output/gw_bandwidths.rds")
}

#-----------------------------------------------
# Array = one BOOT per task                       ####
boots <- names(readRDS("output/boot_list.rds"))
for(bt in boots){ dir.create(file.path(study_environment$wd$boots, bt), showWarnings = FALSE, recursive = TRUE) }
if(!is.na(as.numeric(Sys.getenv("SLURM_ARRAY_TASK_ID")))){
  boots <- boots[as.numeric(Sys.getenv("SLURM_ARRAY_TASK_ID"))]
}

bw_cache <- tryCatch(readRDS("output/gw_bandwidths.rds"), error = function(e) NULL)
nw <- max(1, as.numeric(Sys.getenv("SLURM_CPUS_PER_TASK", unset = "1")))
if(nw > 1) future::plan(future::multicore, workers = nw) else future::plan(future::sequential)

#-----------------------------------------------
# Estimations (boot x cell; 50 specs in parallel; consensus)  ####
for(boot in boots){
  Counties_geo <- tryCatch({
    USMUR <- rast(file.path(study_environment$gssurgo_archive,"MURASTER_30m.tif"))
    Cty   <- vect(file.path(study_environment$usaPolygons_archive,"USA_Counties.shp"))
    Cty   <- crop(project(Cty, crs(USMUR)), ext(USMUR))
    Cty$fip <- as.character(paste0(stringr::str_pad(as.numeric(as.character(Cty$STATEFP)),2,pad="0"),
                                   stringr::str_pad(as.numeric(as.character(Cty$COUNTYFP)),3,pad="0")))
    rm(USMUR); Cty
  }, error = function(e) NULL)

  for(ci in 1:nrow(CELLS)){
    tryCatch({
      cell <- CELLS[ci, ]
      out_put_file <- file.path(study_environment$wd$boots, boot,
        paste0("consensus_cluster_", cell$crop, "_period", cell$period, "_", cell$NAME, ".rds"))
      if(file.exists(out_put_file)) next

      base <- estimate_cell_base(boot = boot, cell = cell)
      if(is.null(base)) next

      bw_df <- if(!is.null(bw_cache)) bw_cache[[paste0("cell", cell$cellN)]] else NULL
      spec_res <- future.apply::future_lapply(1:nrow(spec_gw), function(s){
        bl <- NULL; bg <- NULL
        if(!is.null(bw_df)){
          row <- bw_df[bw_df$specN == spec_gw$specN[s], ]
          if(nrow(row)){
            bl <- suppressWarnings(as.numeric(row$bw_lag[1])); if(!isTRUE(is.finite(bl))) bl <- NULL
            bg <- suppressWarnings(as.numeric(row$bw_gwr[1])); if(!isTRUE(is.finite(bg))) bg <- NULL
          }
        }
        estimate_cell_gw(base, spec_gw[s, ], bw_lag = bl, bw_gwr = bg)
      }, future.seed = TRUE)

      res <- reduce_consensus(base, spec_res, cell = cell, boot = boot, Counties = Counties_geo)
      if(!is.null(res)) saveRDS(res, file = out_put_file)
      rm(base, spec_res, res); gc()
    }, error = function(e){return(NULL)})
  }
  rm(Counties_geo); gc()
}
#-----------------------------------------------
