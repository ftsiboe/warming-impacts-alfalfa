#-----------------------------------------------
# Preliminaries                              ####
rm(list=ls(all=TRUE));gc();library(magrittr);library(tidyverse);library(data.table)
library(plm);library(car);library(lmtest);library(terra);library(sp);library(sf)
study_environment <- readRDS("data/study_environment.rds")
invisible(lapply(list.files("scripts/helpers", pattern = "[.]R$", full.names = TRUE), source))
# devtools::document(file.path(dirname(dirname(getwd())),"packages/rAgroClimate"))

# NOTE: no geographically weighted (GW) estimation happens in this script anymore.
# It emits ONLY the per-boot yield-impact quantities. Warming-scenario alfalfa
# availability and the cattle associations are produced downstream by
# 009_..._availability_warming.R (which reuses 003's fixed neighbourhood operator),
# so gwkit / GWmodel and the 50-spec ensemble have been removed from the boot loop.
set.seed(08032024)

# ==============================================================================
# STRUCTURE (CLUSTER variant: county knots+slopes from 005 optimal_knots_cluster)
# ------------------------------------------------------------------------------
# The SLURM array is BY BOOTSTRAP: one array task = one boot ("0000" + 100 draws).
# Within a task, for every analysis cell (crop x window x climate baseline) the
# yield model is computed ONCE - the national piecewise fit, plus per-county
# KNOTS + boot-refit SLOPES from the CLUSTER estimator (005). County impacts use
# each county's own knot & cluster slope, built from prism_climate at the union
# of thresholds 006 retains; aggregates are acreage-weighted. One file per cell
# per boot is written, carrying exposure / piecewise / relation / impact_yield /
# sp_data.
#
# NO geographically weighted estimation happens here. Warming-scenario alfalfa
# availability = prodXX + m(prodXX) and the cattle~availability associations are
# produced downstream by 009_..._availability_warming.R, which applies 003's
# fixed (boot-invariant) neighbourhood operator to prod00*(1+impact/100) for
# every boot at once. Baseline avail00 + associations come from
# 003_..._availability_associations.R.
# ==============================================================================

#-----------------------------------------------
# boot list                                  #### (build once if missing)
# One "0000" full-sample entry + 100 year-block resamples. Deterministic given the
# seed above, so any task that regenerates it produces the identical list. For the
# SLURM array, generate it ONCE before submitting (the first task writes it, the
# rest read it) to avoid a concurrent write.
if(!file.exists("output/boot_list.rds")){
  dir.create("output", showWarnings = FALSE, recursive = TRUE)
  boot_list <- readRDS("data/nass_hay_production.rds")
  boot_list <- list("0000"=unique(boot_list$commodity_year))
  for(bt in 1:100){
    boot_list[[stringr::str_pad(bt,pad="0",4)]] <- sample(boot_list$`0000`, length(boot_list$`0000`), replace = TRUE)
  }
  saveRDS(boot_list,file="output/boot_list.rds")
}

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
# Assemble one boot x cell output (no GW)                    ####
# The yield quantities from estimate_cell_base are simply tagged and written.
# Warming-scenario availability + the cattle associations are produced downstream
# by 009_..._availability_warming.R.
assemble_output <- function(base, cell, boot){
  if(is.null(base)) return(NULL)
  tag <- data.frame(boot = boot, crop = cell$crop, period = cell$period,
                    climate_base = cell$NAME, cellN = cell$cellN,
                    stringsAsFactors = FALSE)
  list(
    exposure     = data.frame(tag, base$exposure,     row.names = NULL),
    piecewise    = data.frame(tag, base$piecewise,    row.names = NULL),
    relation     = data.frame(tag, base$relation,     row.names = NULL),
    impact_yield = data.frame(tag, base$impact_yield, row.names = NULL),
    sp_data      = data.frame(tag, base$sp_data,      row.names = NULL))
}

#-----------------------------------------------
# Array = one BOOT per task                       ####
boots <- names(readRDS("output/boot_list.rds"))
for(bt in boots){ dir.create(file.path(study_environment$wd$boots, bt), showWarnings = FALSE, recursive = TRUE) }
if(!is.na(as.numeric(Sys.getenv("SLURM_ARRAY_TASK_ID")))){
  boots <- boots[as.numeric(Sys.getenv("SLURM_ARRAY_TASK_ID"))]
}

#-----------------------------------------------
# Estimations (boot x cell; yield impacts only, no GW)       ####
for(boot in boots){
  for(ci in 1:nrow(CELLS)){
    tryCatch({
      cell <- CELLS[ci, ]
      out_put_file <- file.path(study_environment$wd$boots, boot,
        paste0("consensus_cluster_", cell$crop, "_period", cell$period, "_", cell$NAME, ".rds"))
      if(file.exists(out_put_file)) next

      base <- estimate_cell_base(boot = boot, cell = cell)
      if(is.null(base)) next

      res <- assemble_output(base, cell = cell, boot = boot)
      if(!is.null(res)) saveRDS(res, file = out_put_file)
      rm(base, res); gc()
    }, error = function(e){return(NULL)})
  }
}
#-----------------------------------------------
