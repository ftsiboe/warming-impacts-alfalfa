#-----------------------------------------------
# Preliminaries                              ####
rm(list=ls(all=TRUE));gc();library(magrittr);library(future.apply);library(tidyverse);library(data.table)
library(plm);library(car);library(lmtest);library(terra);library(GWmodel);library(sp);library(sf)
study_environment <- readRDS("data/study_environment.rds")
invisible(lapply(list.files("scripts/helpers", pattern = "[.]R$", full.names = TRUE), source))
# devtools::document(file.path(dirname(dirname(getwd())),"packages/rAgroClimate"))

# gwkit (sibling package) is the canonical GW toolkit for this project.
devtools::load_all(file.path(dirname(dirname(getwd())),"packages/gwkit"))
set.seed(08032024)

# ==============================================================================
# STRUCTURE (rewritten 2026-07 to adopt gwkit)
# ------------------------------------------------------------------------------
# The SLURM array is BY BOOTSTRAP: one array task = one boot ("0000" + 100 draws).
# Within a task, for every analysis cell (crop x window x climate baseline):
#   * the spec-INVARIANT yield model (exposure, piecewise, relation, yield
#     impacts, county production panel) is computed ONCE; then
#   * the GW quantities are estimated under ALL 50 gwkit distance-metric
#     presets x kernels (10 x 5) in parallel across cores and reduced to a
#     per-county CONSENSUS (median primary + mean) in-task.
# Only the consensus is written (one file per cell per boot).
#
# gwkit tools used:
#   gw_distance_metric_names()               -> the 10 distance presets
#   estimate_gwlag_by_point()                -> availability m(z_i) (neighbour-
#                                               weighted, self-excluded, multi-col)
#   estimate_gwr_coefficients_by_point()     -> cattle~avail00 and
#                                               cattle~prod00+prod00_LM (GWR)
#   gw_optimal_scalar_by_polygon()           -> across-spec consensus diagnostics
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
# BASE (spec-invariant) block for one boot x cell           ####
# Yield model + county production panel + projected county centroids. None of
# this depends on the GW specification, so it is computed once per cell.
estimate_cell_base <- function(boot, cell){
  tryCatch({
    period       <- cell$period
    crop         <- cell$crop
    climate_base <- cell$climate_base[[1]]

    data <- build_hay_weather_panel(
      crop = "hay_alfalfa", target_periods = period,
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
    knots <- knots[(knots$target_periods %in% c(104:112,0) & knots$crop %in% "hay_alfalfa"),]
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
      Estimate=c(nrow(fit$model),length(unique(data$panelid)),length(unique(data$year)),fit$r.squared[1],
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

    # Yield impacts
    sim_climate <- data.table::rbindlist(lapply(list.files("data/prism_climate",full.names = T),
      function(file){ tryCatch(readRDS(file), error=function(e) NULL) }), fill = TRUE);gc()
    sim_climate <- sim_climate[period %in% period]
    sim_climate <- as.data.frame(data.table::rbindlist(lapply(sampled_years, function(nm){
      sim_climate[sim_climate$commodity_year %in% nm,] }), fill = TRUE));gc()
    sim_climate <- sim_climate[sim_climate$commodity_year %in% climate_base,]
    sim_climate$DD1 <- sim_climate[,"dday00"] - sim_climate[,paste0("dday",Tmin)]
    sim_climate$DD2 <- sim_climate[,paste0("dday",Tmin)] - sim_climate[,paste0("dday",Tmax)]
    sim_climate$DD3 <- sim_climate[,paste0("dday",Tmax)]
    setDT(sim_climate)
    sim_climate <- sim_climate[, .(DD1 = mean(DD1, na.rm = TRUE),DD2 = mean(DD2, na.rm = TRUE),DD3 = mean(DD3, na.rm = TRUE)),
                               by = c("county_fips","warming_scenario")]
    sim_climate <- tidyr::separate(sim_climate,"county_fips", into=c("state_code", "county_code"),sep=c(2))
    sim_climate$state_code <- as.numeric(as.character(sim_climate$state_code))
    sim_climate$county_code <- as.numeric(as.character(sim_climate$county_code))
    setDT(sim_climate)
    sim_climate <- sim_climate[warming_scenario == 0, .(state_code, county_code, DD1b = DD1, DD2b = DD2, DD3b = DD3)][
      sim_climate, on=.(state_code, county_code), nomatch=0]
    sim_climate <- sim_climate[!warming_scenario %in% 0]
    us_states <- as.data.frame(urbnmapr::get_urbn_map(map = "states", sf = TRUE))
    us_states$state_code <- as.numeric(as.character(us_states$state_fips))
    sim_climate <- dplyr::inner_join(as.data.frame(sim_climate),us_states[c("state_code","state_name","state_abbv" )],by="state_code")
    sim_climate$region <- ifelse(sim_climate$state_abbv %in% c("CT","DE","ME","MD","MA","NH","NJ","NY","PA","RI","VT"),"Northeast",NA)
    sim_climate$region <- ifelse(sim_climate$state_abbv %in% c("IL","IN","IA","KS","MI","MN","MO","NE","ND","OH","SD","WI"),"Midwest",sim_climate$region)
    sim_climate$region <- ifelse(sim_climate$state_abbv %in% c("AL","AR","FL","GA","KY","LA","MS","NC","OK","SC","TN","TX","VA","WV"),"South",sim_climate$region)
    sim_climate$region <- ifelse(sim_climate$state_abbv %in% c("AZ","CA","CO","HI","ID","MT","NV","NM","OR","UT","WA","WY"),"West",sim_climate$region)
    setDT(sim_climate)
    sim_county <- sim_climate[, .(DD1s = mean(DD1-DD1b, na.rm=TRUE),DD2s = mean(DD2-DD2b, na.rm=TRUE),DD3s = mean(DD3-DD3b, na.rm=TRUE)),
                              by = c("region","state_code","county_code","warming_scenario")]
    prod <- as.data.frame(readRDS("data/nass_hay_production.rds"))
    prod <- prod[prod$commodity_name %in% crop,]
    prod <- as.data.frame(data.table::rbindlist(lapply(sampled_years, function(nm){ prod[prod$commodity_year %in% nm,] }), fill = TRUE))
    prod <- prod[prod$commodity_year %in% climate_base,]
    prod <- dplyr::inner_join(doBy::summaryBy(area~state_code+county_code+commodity_year,data=prod,FU=sum,keep.names = T),
                              unique(as.data.frame(sim_climate)[c("state_code","county_code")]),by=c("state_code","county_code"))
    prod <- doBy::summaryBy(area~state_code+asd_cd+county_code+commodity_year,data=prod,FU=mean,keep.names = T)
    prod <- dplyr::inner_join(doBy::summaryBy(area~state_code+county_code+commodity_year,data=prod,FU=sum,keep.names = T),
                              doBy::summaryBy(area~state_code+county_code,data=prod,FU=sum),by=c("state_code","county_code"))
    prod$weight <- prod$area/prod$area.sum
    prod <- doBy::summaryBy(weight~state_code+county_code,data=prod,FU=mean,keep.names = T)
    setDT(prod)
    sim_climate <- sim_climate[prod, on=.(state_code, county_code), nomatch=0]
    sim_state <- copy(sim_climate)[, .(DD1s = weighted.mean(DD1-DD1b,weight,na.rm=TRUE),DD2s = weighted.mean(DD2-DD2b,weight,na.rm=TRUE),DD3s = weighted.mean(DD3-DD3b,weight,na.rm=TRUE)),
                                   by = c("region","state_code","warming_scenario")]; sim_state[, county_code := 0]
    sim_region <- copy(sim_climate)[, .(DD1s = weighted.mean(DD1-DD1b,weight,na.rm=TRUE),DD2s = weighted.mean(DD2-DD2b,weight,na.rm=TRUE),DD3s = weighted.mean(DD3-DD3b,weight,na.rm=TRUE)),
                                    by = c("region","warming_scenario")]; sim_region[, state_code := 0]; sim_region[, county_code := 0]
    sim_climate <- sim_climate[, .(DD1s = weighted.mean(DD1-DD1b,weight,na.rm=TRUE),DD2s = weighted.mean(DD2-DD2b,weight,na.rm=TRUE),DD3s = weighted.mean(DD3-DD3b,weight,na.rm=TRUE)),
                               by = c("warming_scenario")]; sim_climate[, state_code := 0]; sim_climate[, county_code := 0]; sim_climate[, region := ""]
    impact_yield <- rbind(sim_climate,sim_state,sim_county,sim_region)
    setDT(impact_yield)
    impact_yield[, impact := (exp(DD1s*PWM.COEF["DD1"] + DD2s*PWM.COEF["DD2"] + DD3s*PWM.COEF["DD3"]) - 1)*100]
    impact_yield[, fip := as.character(paste0(stringr::str_pad(as.numeric(as.character(state_code)),2,pad="0"),
                                              stringr::str_pad(as.numeric(as.character(county_code)),3,pad="0")))]
    impact_yield <- as.data.frame(impact_yield)
    impact_yield$warming_scenario <- paste0("cc",stringr::str_pad(as.numeric(as.character(round(impact_yield$warming_scenario*10))),2,pad="0"))
    impact_yield <- impact_yield[c("region","state_code","county_code","fip","warming_scenario","impact")]
    impact_yield <- impact_yield %>% tidyr::spread(warming_scenario, impact)

    # County production panel + projected centroids (spec-invariant)
    USMUR    <- rast(file.path(study_environment$gssurgo_archive,"MURASTER_30m.tif"))
    Counties <- vect(file.path(study_environment$usaPolygons_archive,"USA_Counties.shp"))
    Counties <- crop(project(Counties, crs(USMUR)), ext(USMUR)); rm(USMUR);gc()
    Counties$fip <- as.character(paste0(stringr::str_pad(as.numeric(as.character(Counties$STATEFP)),2,pad="0"),
                                        stringr::str_pad(as.numeric(as.character(Counties$COUNTYFP)),3,pad="0")))
    cxy <- geom(centroids(Counties))[, c("x","y")]
    centroids_all <- data.frame(fip = Counties$fip, longitude = cxy[,"x"], latitude = cxy[,"y"], stringsAsFactors = FALSE)

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
# GW block for one spec (availability + associations) via gwkit   ####
PROD_SFX  <- c("00","05","10","15","20","25","30")
PROD_COLS <- paste0("prod", PROD_SFX)

estimate_cell_gw <- function(base, spec, bw = NULL){
  tryCatch({
    sp_data   <- base$sp_data
    centroids <- base$centroids

    # --- Availability: neighbour-weighted (self-excluded) production, gwkit ---
    src <- sp_data[, c("fip","longitude","latitude", PROD_COLS)]
    lag <- gwkit::estimate_gwlag_by_point(
      data = src, unit = "fip", value_cols = PROD_COLS,
      coords = c("longitude","latitude"), predict_data = centroids,
      distance_metric = spec$distance_metric, kernel = spec$kernel,
      adaptive = FALSE, bw = bw, bw_response = "prod00", include_self = FALSE)
    if(is.null(bw)) bw <- attr(lag, "bandwidth")
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

    # --- Associations: gwkit local GWR coefficients (logged) ---
    # source frame for the local regressions: consensus availability + cattle
    aj <- merge(as.data.table(availability[, c("fip","avail00","prod00","prod00_LM")]),
                as.data.table(sp_data[, c("fip","longitude","latitude","cattle")]), by = "fip")
    aj <- aj[!fip %in% "00000" & is.finite(cattle) & is.finite(avail00) &
               is.finite(prod00) & is.finite(prod00_LM)]
    aj[, `:=`(lcattle    = log(cattle + 1e-9),
              lavail00   = log(avail00 + 1e-9),
              lprod00    = log(prod00 + 1e-9),
              lprod00_LM = log(prod00_LM + 1e-9))]
    aj <- as.data.frame(aj)

    b1 <- gwkit::estimate_gwr_coefficients_by_point(
      aj, unit = "fip", formula = lcattle ~ lavail00,
      coords = c("longitude","latitude"), predict_data = centroids,
      distance_metric = spec$distance_metric, kernel = spec$kernel,
      adaptive = TRUE, bw = bw, terms = "lavail00")
    b2 <- gwkit::estimate_gwr_coefficients_by_point(
      aj, unit = "fip", formula = lcattle ~ lprod00 + lprod00_LM,
      coords = c("longitude","latitude"), predict_data = centroids,
      distance_metric = spec$distance_metric, kernel = spec$kernel,
      adaptive = TRUE, bw = bw, terms = c("lprod00","lprod00_LM"))
    assoc <- data.table::rbindlist(list(b1, b2), fill = TRUE)
    assoc[, name := c(lavail00="avail00", lprod00="prod00", lprod00_LM="prod00_LM")[term]]
    associations <- assoc[, .(fip = unit_id, name, est, se, tv, pv)]

    # national (unweighted) global fits
    g1 <- summary(lm(lcattle ~ lavail00, data = aj))$coef["lavail00", ]
    g2 <- summary(lm(lcattle ~ lprod00 + lprod00_LM, data = aj))$coef[c("lprod00","lprod00_LM"), ]
    nat_assoc <- data.table::data.table(
      fip = "00000", name = c("avail00","prod00","prod00_LM"),
      est = c(g1["Estimate"], g2[,"Estimate"]), se = c(g1["Std. Error"], g2[,"Std. Error"]),
      tv  = c(g1["t value"], g2[,"t value"]),  pv = c(g1["Pr(>|t|)"], g2[,"Pr(>|t|)"]))
    associations <- rbind(nat_assoc, associations)

    list(specN = spec$specN, availability = availability,
         associations = as.data.frame(associations), bw_opt = bw)
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

  # associations consensus
  as_ <- data.table::rbindlist(lapply(seq_along(spec_list), function(i){
    d <- data.table::as.data.table(spec_list[[i]]$associations); d[, .specN := i]; d }), fill = TRUE)
  associations <- as_[, .(est = stats::median(est, na.rm=TRUE),
                          est_specmean = mean(est, na.rm=TRUE),
                          est_specsd = stats::sd(est, na.rm=TRUE),
                          se = stats::median(se, na.rm=TRUE),
                          tv = stats::median(tv, na.rm=TRUE),
                          pv = stats::median(pv, na.rm=TRUE),
                          sign_agreement = mean(sign(est) == sign(stats::median(est, na.rm=TRUE)), na.rm=TRUE)),
                      by = .(fip, name)]
  associations <- data.frame(tag, as.data.frame(associations), row.names = NULL)

  # gwkit consensus diagnostic on the headline coefficient (county rows)
  gwkit_avail00 <- NULL
  if(!is.null(Counties) && exists("gw_optimal_scalar_by_polygon")){
    tryCatch({
      st <- as_[name %in% "avail00" & !fip %in% c("0","00000") & is.finite(est)]
      if(nrow(st) > 0)
        gwkit_avail00 <- as.data.frame(gw_optimal_scalar_by_polygon(
          value_dt = st, unit_col = "fip", polygons = Counties,
          value_col = "est", agg_fun = stats::median, queen_smooth = FALSE))
    }, error = function(e){})
  }

  list(exposure=exposure, piecewise=piecewise, relation=relation, impact_yield=impact_yield,
       availability=availability, associations=associations,
       bw = data.frame(tag, specN = vapply(spec_list, function(z) z$specN, numeric(1)),
                       bw_opt = vapply(spec_list, function(z) z$bw_opt, numeric(1)), row.names = NULL),
       gwkit_avail00 = gwkit_avail00)
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
    if(is.null(base)){ bw_cache[[paste0("cell",CELLS$cellN[ci])]] <- setNames(rep(NA_real_,nrow(spec_gw)), spec_gw$specN); next }
    bws <- future.apply::future_lapply(1:nrow(spec_gw), function(s){
      r <- estimate_cell_gw(base, spec_gw[s, ], bw = NULL); if(is.null(r)) NA_real_ else r$bw_opt
    }, future.seed = TRUE)
    bw_cache[[paste0("cell",CELLS$cellN[ci])]] <- setNames(unlist(bws), spec_gw$specN)
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
        paste0("consensus_", cell$crop, "_period", cell$period, "_", cell$NAME, ".rds"))
      if(file.exists(out_put_file)) next

      base <- estimate_cell_base(boot = boot, cell = cell)
      if(is.null(base)) next

      bw_vec <- if(!is.null(bw_cache)) bw_cache[[paste0("cell", cell$cellN)]] else NULL
      spec_res <- future.apply::future_lapply(1:nrow(spec_gw), function(s){
        bw_s <- if(!is.null(bw_vec)) suppressWarnings(as.numeric(bw_vec[as.character(spec_gw$specN[s])])) else NULL
        if(!is.null(bw_s) && !is.finite(bw_s)) bw_s <- NULL
        estimate_cell_gw(base, spec_gw[s, ], bw = bw_s)
      }, future.seed = TRUE)

      res <- reduce_consensus(base, spec_res, cell = cell, boot = boot, Counties = Counties_geo)
      if(!is.null(res)) saveRDS(res, file = out_put_file)
      rm(base, spec_res, res); gc()
    }, error = function(e){return(NULL)})
  }
  rm(Counties_geo); gc()
}
#-----------------------------------------------
                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                     