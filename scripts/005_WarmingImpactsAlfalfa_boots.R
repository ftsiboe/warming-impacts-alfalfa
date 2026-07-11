#-----------------------------------------------
# Preliminaries                              ####
rm(list=ls(all=TRUE));gc();library(magrittr);library(future.apply);library(tidyverse);library(data.table)
library(plm);library(car);library(lmtest);library(terra);library(GWmodel);library(sp);library(sf)
study_environment <- readRDS("data/study_environment.rds")
invisible(lapply(list.files("scripts/helpers", pattern = "[.]R$", full.names = TRUE), source))
# devtools::document(file.path(dirname(dirname(getwd())),"packages/rAgroClimate"))
set.seed(08032024)
#-----------------------------------------------
# boot list                                  ####
function(){
  boot_list <- readRDS("data/nass_hay_production.rds")
  boot_list <- list("0000"=unique(boot_list$commodity_year))
  for(bt in 1:100){
    boot_list[[stringr::str_pad(bt,pad="0",4)]] <- sample(boot_list$`0000`, length(boot_list$`0000`), replace = TRUE)
  }
  saveRDS(boot_list,file="output/boot_list.rds")
}
#-----------------------------------------------
# optimal gw Specifications                  ####
function(){

  USMUR    <- rast(file.path(study_environment$gssurgo_archive,"MURASTER_30m.tif"))
  Counties <- vect(file.path(study_environment$usaPolygons_archive,"USA_Counties.shp"))
  Counties <- project(Counties, crs(USMUR))
  Counties <- crop(Counties, ext(USMUR))
  Counties$fip <-as.character(paste0(stringr::str_pad(as.numeric(as.character(Counties$STATEFP)), 2, pad = "0"),
                                     stringr::str_pad(as.numeric(as.character(Counties$COUNTYFP)), 3, pad = "0")))

  spec_gw <- rbind(data.frame(p=2.00, theta=0.0, longlat=F,DistName="Euclidean distance metric"),
                   data.frame(p=1.00, theta=0.5, longlat=F,DistName="Manhattan distance metric"),
                   data.frame(p=2.00, theta=0.0, longlat=T,DistName="Great Circle distance metric"),
                   data.frame(p=0.75, theta=0.8, longlat=F,DistName="Coordinate system is rotated by an angle 0.8 in radian"))

  spec_gw <- as.data.frame(
    data.table::rbindlist(lapply(c("gaussian","exponential","bisquare","boxcar","tricube"),
                                 function(kernel){return(data.frame(spec_gw,kernel=kernel))}), fill = TRUE))

  sp_data <- readRDS("data/spatial_representation.rds")[c("fip","yield","inventory")]
  sp_data <- sp_data[complete.cases(sp_data),]
  sp_data <- terra::merge(Counties, sp_data, by="fip")
  sp_data$fold <- sample(1:5, length(sp_data), replace = TRUE,prob=rep(0.2,5))
  # sp::plot(sp_data,"fold")

  res <- as.data.frame(
    data.table::rbindlist(
      lapply(
        1:nrow(spec_gw),
        function(task){

          tryCatch({
            # task <- 1

            res <- as.data.frame(
              data.table::rbindlist(
                lapply(
                  1:5,
                  function(f,task){
                    # f <- 1

                    resAll <- list()
                    data_train <- sp_data[!sp_data$fold %in% f,]
                    data_test  <- sp_data[sp_data$fold %in% f,]

                    tryCatch({
                      capture.output(bw_opt <- bw.gwr(yield~1, data=as(st_as_sf(data_train), "Spatial"),
                                                      approach="CV",kernel=spec_gw$kernel[task],adaptive=FALSE, p=spec_gw$p[task],
                                                      theta=spec_gw$theta[task],longlat=spec_gw$longlat[task],
                                                      dMat=gw.dist(dp.locat=geom(centroids(data_train))[, c("x", "y")],
                                                                   rp.locat=geom(centroids(data_train))[, c("x", "y")],
                                                                   focus=0, p=spec_gw$p[task],theta=spec_gw$theta[task], longlat=spec_gw$longlat[task])), file = nullfile())
                      res <- gwr.predict(yield~1,
                                         data=as(st_as_sf(data_train), "Spatial"),
                                         predictdata=as(st_as_sf(data_test), "Spatial"),
                                         bw=bw_opt, kernel=spec_gw$kernel[task],adaptive=FALSE,
                                         p=spec_gw$p[task],theta=spec_gw$theta[task], longlat=spec_gw$longlat[task])

                      resAll[[length(resAll)+1]] <- data.frame(outcome ="yield",f=f,n=nrow(data_test),e=sum((res$SDF$prediction - data_test$yield)^2,na.rm=T))
                    }, error = function(e){return(NULL)})
                    tryCatch({
                      capture.output(bw_opt <- bw.gwr(inventory~1, data=as(st_as_sf(data_train), "Spatial"),
                                                      approach="CV",kernel=spec_gw$kernel[task],adaptive=FALSE, p=spec_gw$p[task],
                                                      theta=spec_gw$theta[task],longlat=spec_gw$longlat[task],
                                                      dMat=gw.dist(dp.locat=geom(centroids(data_train))[, c("x", "y")],
                                                                   rp.locat=geom(centroids(data_train))[, c("x", "y")],
                                                                   focus=0, p=spec_gw$p[task],theta=spec_gw$theta[task], longlat=spec_gw$longlat[task])), file = nullfile())
                      res <- gwr.predict(inventory~1,
                                         data=as(st_as_sf(data_train), "Spatial"),
                                         predictdata=as(st_as_sf(data_test), "Spatial"),
                                         bw=bw_opt, kernel=spec_gw$kernel[task],adaptive=FALSE,
                                         p=spec_gw$p[task],theta=spec_gw$theta[task], longlat=spec_gw$longlat[task])

                      resAll[[length(resAll)+1]] <- data.frame(outcome ="inventory",f=f,n=nrow(data_test),e=sum((res$SDF$prediction - data_test$inventory)^2,na.rm=T))
                    }, error = function(e){return(NULL)})
                    resAll <- as.data.frame(data.table::rbindlist(resAll, fill = TRUE))

                    return(resAll)
                  },task=task), fill = TRUE))

            res <- res[!res$e %in% c(NaN,Inf,-Inf,NA),]
            res <- data.frame(spec_gw[task,],e=weighted.mean(x=res$e,w=res$n))
            return(res)
          }, error = function(e){return(NULL)})


        }), fill = TRUE))
  res <- res[!res$e %in% c(NaN,Inf,-Inf,NA),]
  res <- res[order(res$e),]
  saveRDS(res,file="output/optimal_gw.rds")
  rm(sp_data,spec_gw,Counties,USMUR,res)
}

#-----------------------------------------------
# Specifications                             ####
gc()
optimal_gw <- readRDS("output/optimal_gw.rds")

SPECS <- data.frame(crop="hay_alfalfa",period=107,NAME=c("1991_2020","1981_2010","1971_2000","1961_1990","1951_2022"))
SPECS$climate_base[1] <- list(1991:2020) # NOAA [https://www.ncei.noaa.gov/products/land-based-station/us-climate-normals]
SPECS$climate_base[2] <- list(1981:2010) # NOAA [https://journals.ametsoc.org/view/journals/bams/93/11/bams-d-11-00197.1.xml]
SPECS$climate_base[3] <- list(1971:2000) # NOAA
SPECS$climate_base[4] <- list(1961:1990) # IPCC
SPECS$climate_base[5] <- list(1951:2022)

SPECS_crop <- data.frame(crop=c("hay_other","hay_all"),period=107,NAME="1991_2020")
SPECS_crop$climate_base[1] <- list(1991:2020)
SPECS_crop$climate_base[2] <- list(1991:2020)

SPECS_period <- data.frame(crop=c("hay_alfalfa"),period=c(0:12,101:106,108:112),NAME="1991_2020")
for(iii in 1:nrow(SPECS_period)){
  SPECS_period$climate_base[iii] <- list(1991:2020)
}

SPECS <- unique(rbind(SPECS,SPECS_crop,SPECS_period))

SPECS <- rbind(data.frame(optimal_gw,SPECS[(SPECS$crop %in% "hay_alfalfa" & SPECS$period %in% 107 & SPECS$NAME %in% "1991_2020"),]),
               data.frame(optimal_gw[1,],SPECS[!(SPECS$crop %in% "hay_alfalfa" & SPECS$period %in% 107 & SPECS$NAME %in% "1991_2020"),]))

SPECS$specN <- 1:nrow(SPECS)

SPECS <- as.data.frame(
  data.table::rbindlist(
    lapply(names(readRDS("output/boot_list.rds")),
           function(boot){
             dir.create(file.path(study_environment$wd$boots,boot))
             return(data.frame(boot=boot,SPECS))}), fill = TRUE))

function(){
  SPECS <- SPECS[! file.path(
    study_environment$wd$boots,
    SPECS$boot,
    paste0("spec",stringr::str_pad(SPECS$specN,pad="0",2),"_",
           SPECS$crop,"_period",SPECS$period,"_",SPECS$NAME,".rds")
  ) %in% list.files(study_environment$wd$boots,full.names = T,recursive = T),]
}

if(!is.na(as.numeric(Sys.getenv("SLURM_ARRAY_TASK_ID")))){
  SPECS$TASK <- rep(as.numeric(Sys.getenv("SLURM_ARRAY_TASK_MIN")):as.numeric(Sys.getenv("SLURM_ARRAY_TASK_MAX")), length=nrow(SPECS))
  SPECS <- SPECS[SPECS$TASK %in% as.numeric(Sys.getenv("SLURM_ARRAY_TASK_ID")),]
}

rm(optimal_gw,SPECS_crop,SPECS_period)
#-----------------------------------------------
# Estimations                                ####

lapply(
  c(1:nrow(SPECS)),
  function(spec){
    tryCatch({
      # spec <- 1
      boot <- SPECS$boot[spec]
      period <- SPECS$period[spec]
      crop <- SPECS$crop[spec]
      climate_base <- SPECS$climate_base[spec][[1]]

      out_put_file <- file.path(
        study_environment$wd$boots,boot,
        paste0("spec",stringr::str_pad(SPECS$specN[spec],pad="0",2),"_",
               SPECS$crop[spec],"_period",SPECS$period[spec],"_",SPECS$NAME[spec],".rds"))

      if(!file.exists(out_put_file)){
        #-----------------------------------------------
        # data                                       ####
        print("data")

        data <- build_hay_weather_panel(
          crop = "hay_alfalfa",
          target_periods = period,
          prism_weather_directory = "data/prism_weather")

        sampled_years <- readRDS("output/boot_list.rds")[[boot]]
        data <- as.data.frame(
          data.table::rbindlist(
            lapply(
              1:length(sampled_years),
              function(yr){
                tryCatch({
                  data <- data[data$commodity_year %in% sampled_years[yr],]
                  data$year <- yr
                  return( data )
                }, error = function(e){return(NULL)})
              }), fill = TRUE))
        #-----------------------------------------------
        # Exposure                                   ####
        print("Exposure")
        Expos <- data
        Expos$id <- 1:nrow(Expos)
        Expos <- Expos[c("id",names(Expos)[grepl(paste0("exp"),names(Expos))])]
        Expos <- Expos %>% tidyr::gather(Temp, exp,2:ncol(Expos))
        Expos$Temp <- as.numeric(gsub("[^0-9]","",gsub("exp","",Expos$Temp)))
        Expos <- doBy::summaryBy(exp~Temp,FUN=mean,keep.names=T,data=Expos,na.rm=T)
        # plot(Expos$Temp,Expos$exp)
        #-----------------------------------------------
        # Piecewise linear                           ####
        print("Piecewise")
        knots <- readRDS("output/optimal_knots.rds")
        knots <- knots[(knots$target_periods %in% c(104:112,0) & knots$crop %in% "hay_alfalfa"),]
        knots <- knots[order(knots$cv_error),]
        Tmin <- knots$Tmin[1]
        Tmax <- knots$Tmax[1]
        data$DD1 <- data[,paste0("dday00")] - data[,paste0("dday",stringr::str_pad(Tmin,pad="0",2))]
        data$DD2 <- data[,paste0("dday",stringr::str_pad(Tmin,pad="0",2))] - data[,paste0("dday",stringr::str_pad(Tmax,pad="0",2))]
        data$DD3 <- data[,paste0("dday",stringr::str_pad(Tmax,pad="0",2))]
        data$DD1 <- ifelse(data$DD1<0,0,data$DD1)
        data$DD2 <- ifelse(data$DD2<0,0,data$DD2)
        data$DD3 <- ifelse(data$DD3<0,0,data$DD3)

        fit_data <- doBy::summaryBy(lny + ppt + ppt2 + ppt3 + Trend + Trend2 + DD1 + DD2 + DD3 + freeze~state_code+fip+commodity_year+year,
                                    data=data,FUN=mean,keep.names = T,na.rm=T)

        for(st in unique(fit_data$state_code)){
          fit_data[,paste0("trend1",st)] <- ifelse(fit_data$state_code %in% st,fit_data$Trend,0)
          fit_data[,paste0("trend2",st)] <- ifelse(fit_data$state_code %in% st,fit_data$Trend2,0)
        }

        fit_data$panel <- fit_data$fip
        fit_data$time  <- fit_data$year
        fit_data<- pdata.frame(fit_data, index = c("panel", "time"), drop.index = TRUE)

        fit <- plm(as.formula(
          paste0("lny ~",paste0(c("ppt","ppt2","DD1" ,"DD2" ,"DD3",
                                  names(fit_data)[grepl("trend",names(fit_data))]),collapse = "+"))),
          data = fit_data,model = "within")

        # summary(fit)

        PWM.COEF <- coef(fit)

        G <- length(unique(fit_data$year))
        PWM.VCOV <- G/(G - 1) * vcovHC(fit, type = "HC1", cluster = "time")
        test1 <- linearHypothesis(fit, c("DD1=0","DD2=0","DD3=0"),vcov. = PWM.VCOV)
        test2 <- linearHypothesis(fit, c("ppt=0","ppt2=0"),vcov. = PWM.VCOV)
        test3 <- linearHypothesis(fit, c("DD1=0","DD2=0","DD3=0","ppt=0","ppt2=0"),vcov. = PWM.VCOV)
        test4 <- plmtest(fit$formula, data=fit_data, effect="individual", type="honda",vcov. = PWM.VCOV)

        PWM.TAB <- as.data.frame(coeftest(fit,G/(G - 1) * vcovHC(fit, type = "HC1", cluster = "time"))[,])
        PWM.TAB$name <- rownames(PWM.TAB)
        names(PWM.TAB) <- c("Estimate","StdError","t_value","p_value","name")
        fit <- summary(fit)
        PWM.TAB <- rbind(PWM.TAB,
                         data.frame(
                           name = c("n","n_farms","n_year","r.squared","test_temp","test_ppt","test_weather","test_fe"),
                           Estimate=c(nrow(fit$model),length(unique(data$panelid)),length(unique(data$year)),fit$r.squared[1],
                                      test1$`Chisq`[2],test2$`Chisq`[2],test3$`Chisq`[2],test4$statistic),
                           StdError=NA,t_value=NA,
                           p_value=c(0,0,0,0,test1$`Pr(>Chisq)`[2],test2$`Pr(>Chisq)`[2],test3$`Pr(>Chisq)`[2],test4$p.value)))

        rm(test1,test2,test3,test4,fit,G,fit_data);gc()

        #-----------------------------------------------
        # Nonlinear Relation                         ####
        print("Nonlinear")
        Relation    <- data.frame(Temp=1:45)
        Relation$I1<-as.numeric(Relation$Temp>=0)
        Relation$I2<-as.numeric(Relation$Temp>=Tmin)
        Relation$I3<-as.numeric(Relation$Temp>=Tmax)

        funcx <- paste0("~(1-",Relation$I2,")*(",Relation$Temp,"-0)*DD1+",
                        Relation$I2,"*(1-",Relation$I3,")*(DD1*",Tmin," + DD2*(",Relation$Temp,"-",Tmin,")) +",
                        Relation$I3,"*(DD1*",Tmin," + DD2*(",Tmax,"-",Tmin,") + DD3*(",Relation$Temp,"-",Tmax,"))")
        Form <- as.formula(funcx[1])
        for(i in 2:length(funcx)){Form<-c(Form,as.formula(funcx[i]))}
        Relation <- compute_delta_method(func=Form,vcMat=PWM.VCOV,coefs=PWM.COEF)[c("Estimate","SE")]
        names(Relation) <- c("Piece","PieceSE")
        Relation <- as.data.frame(Relation)
        Relation$Temp <- 1:45
        rm(Form);gc()
        # plot(Relation$Temp,Relation$Piece)
        #-----------------------------------------------
        # Yield impacts                              ####
        print("Yield impacts ")

        sim_climate <- data.table::rbindlist(
          lapply(
            list.files("data/prism_climate",full.names = T),
            function(file){
              tryCatch({
                return(readRDS(file))
              }, error=function(e){NULL})
            }), fill = TRUE);gc()
        target_periods <- period
        sim_climate <- sim_climate[period %in% target_periods]

        sim_climate <- as.data.frame(
          data.table::rbindlist(
            lapply(
              sampled_years,
              function(nm){
                return(sim_climate[sim_climate$commodity_year %in% nm,])
              }), fill = TRUE));gc()

        sim_climate <- sim_climate[sim_climate$commodity_year %in% climate_base,]

        # xxx <- unique(as.data.table(sim_climate)[, c("state_code","county_code"),with=FALSE])

        sim_climate$DD1 <- sim_climate[,paste0("dday00")] - sim_climate[,paste0("dday",Tmin)]
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

        sim_climate$region <- ifelse(
          sim_climate$state_abbv %in% c("CT", "DE", "ME", "MD", "MA", "NH", "NJ", "NY", "PA", "RI", "VT"),"Northeast",
          NA)
        sim_climate$region <- ifelse(
          sim_climate$state_abbv %in% c("IL", "IN", "IA", "KS", "MI", "MN", "MO", "NE", "ND", "OH", "SD", "WI"),"Midwest",
          sim_climate$region)
        sim_climate$region <- ifelse(
          sim_climate$state_abbv %in% c("AL", "AR", "FL", "GA", "KY", "LA", "MS", "NC", "OK", "SC", "TN", "TX", "VA", "WV"),"South",
          sim_climate$region)
        sim_climate$region <- ifelse(
          sim_climate$state_abbv %in% c("AZ", "CA", "CO", "HI", "ID", "MT", "NV", "NM", "OR", "UT", "WA", "WY"),"West",
          sim_climate$region)

        setDT(sim_climate)

        sim_county <- sim_climate
        sim_county <- sim_county[, .(DD1s = mean(x=DD1-DD1b, na.rm = TRUE),
                                     DD2s = mean(x=DD2-DD2b, na.rm = TRUE),
                                     DD3s = mean(x=DD3-DD3b, na.rm = TRUE)),
                                 by = c("region","state_code","county_code","warming_scenario")]

        prod <- as.data.frame(readRDS("data/nass_hay_production.rds"))

        prod <- prod[prod$commodity_name %in% crop,]

        prod <- as.data.frame(
          data.table::rbindlist(
            lapply(
              sampled_years,
              function(nm){
                return( prod <- prod[prod$commodity_year %in% nm,])
              }), fill = TRUE))

        prod <- prod[prod$commodity_year %in% climate_base,]
        prod <- dplyr::inner_join(doBy::summaryBy(area~state_code+county_code+commodity_year,data=prod,FU=sum,keep.names = T),
                                  unique(as.data.frame(sim_climate)[c("state_code","county_code")]),by=c("state_code","county_code"))
        prod <- doBy::summaryBy(area~state_code+asd_cd+county_code+commodity_year,data=prod,FU=mean,keep.names = T)
        prod <- dplyr::inner_join(doBy::summaryBy(area~state_code+county_code+commodity_year,data=prod,FU=sum,keep.names = T),
                                  doBy::summaryBy(area~state_code+county_code,data=prod,FU=sum),
                                  by=c("state_code","county_code"))
        prod$weight <- prod$area/prod$area.sum
        prod <- doBy::summaryBy(weight~state_code+county_code,data=prod,FU=mean,keep.names = T)
        setDT(prod)

        sim_climate <- sim_climate[prod, on=.(state_code, county_code), nomatch=0]

        sim_state <- copy(sim_climate)[, .(DD1s = weighted.mean(x=DD1-DD1b,w=weight, na.rm = TRUE),
                                           DD2s = weighted.mean(x=DD2-DD2b,w=weight, na.rm = TRUE),
                                           DD3s = weighted.mean(x=DD3-DD3b,w=weight, na.rm = TRUE)),
                                       by = c("region","state_code","warming_scenario")]
        sim_state[, county_code := 0]

        sim_region <- copy(sim_climate)[, .(DD1s = weighted.mean(x=DD1-DD1b,w=weight, na.rm = TRUE),
                                            DD2s = weighted.mean(x=DD2-DD2b,w=weight, na.rm = TRUE),
                                            DD3s = weighted.mean(x=DD3-DD3b,w=weight, na.rm = TRUE)),
                                        by = c("region","warming_scenario")]
        sim_region[, state_code := 0]
        sim_region[, county_code := 0]

        sim_climate <- sim_climate[, .(DD1s = weighted.mean(x=DD1-DD1b,w=weight, na.rm = TRUE),
                                       DD2s = weighted.mean(x=DD2-DD2b,w=weight, na.rm = TRUE),
                                       DD3s = weighted.mean(x=DD3-DD3b,w=weight, na.rm = TRUE)),
                                   by = c("warming_scenario")]
        sim_climate[, state_code := 0]
        sim_climate[, county_code := 0]
        sim_climate[, region := ""]

        impact_yield <- rbind(sim_climate,sim_state,sim_county,sim_region)

        setDT(impact_yield)
        impact_yield[, impact := (exp(DD1s*PWM.COEF["DD1"] + DD2s*PWM.COEF["DD2"] + DD3s*PWM.COEF["DD3"]) - 1)*100]

        impact_yield[, fip := as.character(paste0(stringr::str_pad(as.numeric(as.character(state_code)), 2, pad = "0"),
                                                  stringr::str_pad(as.numeric(as.character(county_code)), 3, pad = "0")))]
        impact_yield <- as.data.frame(impact_yield)
        impact_yield$warming_scenario <- paste0("cc",stringr::str_pad(as.numeric(as.character(round(impact_yield$warming_scenario*10))), 2, pad = "0"))
        impact_yield <- impact_yield[c("region","state_code", "county_code","fip", "warming_scenario","impact")]
        impact_yield <- impact_yield %>% tidyr::spread(warming_scenario, impact)

        #rm(sim_climate,prod,sim_county);gc()

        #-----------------------------------------------
        # Availability                               ####
        print("Availability")
        USMUR    <- rast(file.path(study_environment$gssurgo_archive,"MURASTER_30m.tif"))
        Counties <- vect(file.path(study_environment$usaPolygons_archive,"USA_Counties.shp"))
        Counties <- project(Counties, crs(USMUR))
        Counties <- crop(Counties, ext(USMUR))
        rm(USMUR);gc()
        Counties$fip <- as.character(paste0(stringr::str_pad(as.numeric(as.character(Counties$STATEFP)), 2, pad = "0"),
                                            stringr::str_pad(as.numeric(as.character(Counties$COUNTYFP)), 3, pad = "0")))

        # Counties <- Counties[Counties$STATEFP %in% "20",]
        sp_data <- as.data.frame(readRDS("data/spatial_representation.rds"))
        sp_data <- dplyr::inner_join(impact_yield,sp_data,by=c("state_code", "county_code", "fip"))
        names(sp_data)[names(sp_data) %in% "inventory"] <- "cattle"

        sp_data <- terra::merge(Counties, sp_data, by="fip")
        sp_data$prod00 <- sp_data$area*sp_data$yield
        sp_data$prod00 <- ifelse(sp_data$prod00 %in% c(NA,Inf,NaN,-Inf),0,sp_data$prod00)
        tryCatch({sp_data$prod05 <- ifelse(sp_data$prod00 %in% c(NA,Inf,NaN,-Inf,0),0,sp_data$area*sp_data$yield*(1+(sp_data$cc05/100)))}, error=function(e){})
        tryCatch({sp_data$prod10 <- ifelse(sp_data$prod00 %in% c(NA,Inf,NaN,-Inf,0),0,sp_data$area*sp_data$yield*(1+(sp_data$cc10/100)))}, error=function(e){})
        tryCatch({sp_data$prod15 <- ifelse(sp_data$prod00 %in% c(NA,Inf,NaN,-Inf,0),0,sp_data$area*sp_data$yield*(1+(sp_data$cc15/100)))}, error=function(e){})
        tryCatch({sp_data$prod20 <- ifelse(sp_data$prod00 %in% c(NA,Inf,NaN,-Inf,0),0,sp_data$area*sp_data$yield*(1+(sp_data$cc20/100)))}, error=function(e){})
        tryCatch({sp_data$prod25 <- ifelse(sp_data$prod00 %in% c(NA,Inf,NaN,-Inf,0),0,sp_data$area*sp_data$yield*(1+(sp_data$cc25/100)))}, error=function(e){})
        tryCatch({sp_data$prod30 <- ifelse(sp_data$prod00 %in% c(NA,Inf,NaN,-Inf,0),0,sp_data$area*sp_data$yield*(1+(sp_data$cc30/100)))}, error=function(e){})

        # sp::plot(sp_data,"prod05")
        bw_opt <- 50778.4
        dmat <- GWmodel::gw.dist(dp.locat=geom(centroids(sp_data))[, c("x", "y")],
                                 rp.locat=geom(centroids(sp_data))[, c("x", "y")],
                                 focus=0, p=SPECS$p[spec],theta=SPECS$theta[spec], longlat=SPECS$longlat[spec])

        capture.output(bw_opt <- GWmodel::bw.gwr(prod00~1, data=as(st_as_sf(sp_data), "Spatial"),
                                                 approach="CV",kernel=SPECS$kernel[spec],adaptive=FALSE, p=SPECS$p[spec],
                                                 theta=SPECS$theta[spec],longlat=SPECS$longlat[spec],dMat=dmat), file = nullfile())

        rm(dmat);gc()

        county_i <- Counties
        county_j <- sp_data

        # sp::plot(county_j,"prod00")
        # sp::plot(county_i,"AWATER")

        dmat <- GWmodel::gw.dist(dp.locat=geom(centroids(county_j))[, c("x", "y")],
                                 rp.locat=geom(centroids(county_i))[, c("x", "y")],
                                 focus=0, p=SPECS$p[spec],theta=SPECS$theta[spec], longlat=SPECS$longlat[spec])

        wmat <- GWmodel::gw.weight(vdist=dmat,bw=bw_opt,kernel=SPECS$kernel[spec],adaptive=FALSE)

        row.names(wmat) <- sp_data$fip;colnames(wmat) <-  Counties$fip
        for (i in seq_len(nrow(wmat))) {
          for (j in seq_len(ncol(wmat))) {
            if (rownames(wmat)[i] == colnames(wmat)[j]) {
              wmat[i, j] <- 0
            }
          }
        }

        datasim <- as.data.frame(Counties)
        setDT(datasim)
        lm_fxn <- function(fip){
          tryCatch({
            # fip <- datasim[2,"fip"][[1]]

            gwr_fit <- as.data.frame(sp_data[c("prod00","prod05","prod10","prod15","prod20","prod25","prod30")],wij =wmat[,paste0(fip)])
            gwr_fit$fip <- fip
            setDT(gwr_fit)

            gwr_fit <- gwr_fit[
              , .(prod00_LM = weighted.mean(x=prod00,w=wij, na.rm= TRUE),
                  prod05_LM = weighted.mean(x=prod05,w=wij, na.rm= TRUE),
                  prod10_LM = weighted.mean(x=prod10,w=wij, na.rm= TRUE),
                  prod15_LM = weighted.mean(x=prod15,w=wij, na.rm= TRUE),
                  prod20_LM = weighted.mean(x=prod20,w=wij, na.rm= TRUE),
                  prod25_LM = weighted.mean(x=prod25,w=wij, na.rm= TRUE),
                  prod30_LM = weighted.mean(x=prod30,w=wij, na.rm= TRUE)),
              by = .(fip)]

            return(gwr_fit)
          }, error = function(e){return(NULL)})
        }

        availability <- rbindlist(datasim[, .(value = list(lm_fxn(fip=fip))), by = c("fip")]$value)
        availability <- dplyr::full_join(as.data.frame(sp_data)[c("fip",names(sp_data)[grepl("prod",names(sp_data))])],
                                         availability,by="fip")

        setDT(availability)
        availability[, avail00 := ifelse(prod00 %in% NA,0,prod00) + prod00_LM]
        availability[, avail05 := ifelse(prod05 %in% NA,0,prod05) + prod05_LM]
        availability[, avail10 := ifelse(prod10 %in% NA,0,prod10) + prod10_LM]
        availability[, avail15 := ifelse(prod15 %in% NA,0,prod15) + prod15_LM]
        availability[, avail20 := ifelse(prod20 %in% NA,0,prod20) + prod20_LM]
        availability[, avail25 := ifelse(prod25 %in% NA,0,prod25) + prod25_LM]
        availability[, avail30 := ifelse(prod30 %in% NA,0,prod30) + prod30_LM]

        availability <- availability[!avail00  %in% NA,
                                     .(fip,prod00,prod05,prod10,prod15,prod20,prod25,prod30,
                                       prod00_LM,prod05_LM,prod10_LM,prod15_LM,prod20_LM,prod25_LM,prod30_LM,
                                       avail00,avail05,avail10,avail15,avail20,avail25,avail30)]

        availability00 <- as.data.frame(availability)
        setDT(availability00)
        availability00[, fip := "00000"]

        availability00 <- availability00[
          , .(avail00 = mean(x=avail00, na.rm= TRUE),
              avail05 = mean(x=avail05, na.rm= TRUE),
              avail10 = mean(x=avail10, na.rm= TRUE),
              avail15 = mean(x=avail15, na.rm= TRUE),
              avail20 = mean(x=avail20, na.rm= TRUE),
              avail25 = mean(x=avail25, na.rm= TRUE),
              avail30 = mean(x=avail30, na.rm= TRUE),
              prod00 = mean(x=prod00, na.rm= TRUE),
              prod05 = mean(x=prod05, na.rm= TRUE),
              prod10 = mean(x=prod10, na.rm= TRUE),
              prod15 = mean(x=prod15, na.rm= TRUE),
              prod20 = mean(x=prod20, na.rm= TRUE),
              prod25 = mean(x=prod25, na.rm= TRUE),
              prod30 = mean(x=prod30, na.rm= TRUE),
              prod00_LM = mean(x=prod00_LM, na.rm= TRUE),
              prod05_LM = mean(x=prod05_LM, na.rm= TRUE),
              prod10_LM = mean(x=prod10_LM, na.rm= TRUE),
              prod15_LM = mean(x=prod15_LM, na.rm= TRUE),
              prod20_LM = mean(x=prod20_LM, na.rm= TRUE),
              prod25_LM = mean(x=prod25_LM, na.rm= TRUE),
              prod30_LM = mean(x=prod30_LM, na.rm= TRUE)),
          by = .(fip)]

        availability <- rbind(availability00,availability)

        rm(datasim,dmat,wmat,county_i,county_j);gc()

        #-----------------------------------------------
        # Associations                               ####
        print("Associations")
        sp_data$prod00_LM <- as.numeric(as.character(factor(sp_data$fip,levels = availability$fip, labels = availability$prod00_LM)))
        sp_data$prod00    <- as.numeric(as.character(factor(sp_data$fip,levels = availability$fip, labels = availability$prod00)))
        sp_data$avail00   <- as.numeric(as.character(factor(sp_data$fip,levels = availability$fip, labels = availability$avail00)))
        print("Associations")
        county_i <- Counties
        county_j <- sp_data[!sp_data$cattle %in% NA,]
        county_j <- county_j[!county_j$avail00 %in% NA,]

        # sp::plot(sp_data,"avail00")
        # sp::plot(sp_data,"cattle")
        print("Associations")
        dmat <- GWmodel::gw.dist(dp.locat=geom(centroids(county_j))[, c("x", "y")],
                                 rp.locat=geom(centroids(county_i))[, c("x", "y")],
                                 focus=0, p=SPECS$p[spec],theta=SPECS$theta[spec], longlat=SPECS$longlat[spec])

        wmat <- GWmodel::gw.weight(vdist=dmat,bw=bw_opt,kernel=SPECS$kernel[spec],adaptive=TRUE)
        print("Associations")
        datasim <- as.data.frame(Counties)
        setDT(datasim)
        lm_fxn <- function(fip,group_info){
          tryCatch({
            # fip <- datasim[1,"fip"]

            dataxx <- as.data.frame(county_j)
            dataxx$wij <- wmat[,c(1:length(county_i))[county_i$fip %in% fip]]

            dataxx <- dataxx[!dataxx$prod00_LM %in% c(NA,Inf,-Inf,NaN),]
            dataxx <- dataxx[!dataxx$cattle %in% c(NA,Inf,-Inf,NaN),]
            dataxx <- dataxx[!dataxx$avail00 %in% c(NA,Inf,-Inf,NaN),]
            dataxx <- dataxx[!dataxx$prod00 %in% c(NA,Inf,-Inf,NaN),]
            dataxx$cattle <- log(dataxx$cattle+0.000000001)
            dataxx$avail00 <- log(dataxx$avail00+0.000000001)
            dataxx$prod00 <- log(dataxx$prod00+0.000000001)
            dataxx$prod00_LM <- log(dataxx$prod00_LM+0.000000001)

            avail_fit <- as.data.frame(coef(summary(lm(cattle~avail00,data=dataxx,weights=wij))))
            avail_fit <- avail_fit["avail00",]
            names(avail_fit) <- c("est","se","tv","pv")
            avail_fit$name <- row.names(avail_fit)

            prod_fit <- as.data.frame(coef(summary(lm(cattle~prod00+prod00_LM,data=dataxx,weights=wij))))
            prod_fit <- prod_fit[c("prod00","prod00_LM"),]
            names(prod_fit) <- c("est","se","tv","pv")
            prod_fit$name <- row.names(prod_fit)

            gwr_fit <- rbind(avail_fit,prod_fit)
            setDT(gwr_fit)

            for (group_var in names(group_info)) {
              gwr_fit[, (group_var) := group_info[[group_var]]]
            }
            return(gwr_fit)
          }, error = function(e){return(NULL)})
        }
        print("Associations")
        associations <- rbindlist(datasim[, .(value = list(lm_fxn(fip=fip,group_info=mget(c("fip"))))), by = c("fip")]$value)
        print("Associations")
        dataxx <- as.data.frame(county_j)
        dataxx <- dataxx[!dataxx$prod00_LM %in% c(NA,Inf,-Inf,NaN),]
        dataxx <- dataxx[!dataxx$cattle %in% c(NA,Inf,-Inf,NaN),]
        dataxx <- dataxx[!dataxx$avail00 %in% c(NA,Inf,-Inf,NaN),]
        dataxx <- dataxx[!dataxx$prod00 %in% c(NA,Inf,-Inf,NaN),]
        dataxx$cattle <- log(dataxx$cattle+0.000000001)
        dataxx$avail00 <- log(dataxx$avail00+0.000000001)
        dataxx$prod00 <- log(dataxx$prod00+0.000000001)
        dataxx$prod00_LM <- log(dataxx$prod00_LM+0.000000001)

        print("Associations")
        avail_fit0 <- data.frame(coef(summary(lm(cattle~avail00,data=dataxx))))
        avail_fit0 <- avail_fit0["avail00",]
        names(avail_fit0) <- c("est","se","tv","pv")
        avail_fit0$name <- row.names(avail_fit0)
        print("Associations")
        prod_fit0 <- data.frame(coef(summary(lm(log(cattle)~prod00+prod00_LM,data=dataxx))))
        prod_fit0 <- prod_fit0[c("prod00","prod00_LM"),]
        names(prod_fit0) <- c("est","se","tv","pv")
        prod_fit0$name <- row.names(prod_fit0)
        acco_fit0 <- rbind(avail_fit0,prod_fit0)
        acco_fit0$fip <- "00000"
        setDT(acco_fit0)
        associations <- rbind(acco_fit0,associations)

        #-----------------------------------------------
        # Finalize                                   ####
        print("Finalize")
        res <- list(
          exposure=as.data.frame(Expos),
          piecewise = as.data.frame(PWM.TAB),
          relation = as.data.frame(Relation),
          impact_yield = as.data.frame(impact_yield),
          availability = as.data.frame(availability),
          associations = as.data.frame(associations)
        )

        for(xx in 1:length(res)){
          res[[xx]] <- data.frame(SPECS[spec,c("boot","p","theta","longlat","DistName","kernel","crop","period","specN" )],
                                  climate_base=SPECS$NAME[spec],res[[xx]])
        }

        res[["bw"]] <- data.frame(SPECS[spec,c("boot","p","theta","longlat","DistName","kernel","crop","period","specN" )],bw_opt=bw_opt)

        saveRDS(res,file=out_put_file)
        #-----------------------------------------------
      }

      return(spec)
      }, error = function(e){return(NULL)})
  })

#-----------------------------------------------
