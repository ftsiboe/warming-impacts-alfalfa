#-----------------------------------------------
# Preliminaries                              ####
rm(list=ls(all=TRUE));gc();library(magrittr);library(future.apply);library(tidyverse);library(data.table)
library(plm);library(car);library(lmtest);library(terra);library(GWmodel);library(sp);library(sf)
study_environment <- readRDS("data/study_environment.rds")
invisible(lapply(list.files("scripts/helpers", pattern = "[.]R$", full.names = TRUE), source))
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

if(!is.na(as.numeric(Sys.getenv("SLURM_ARRAY_TASK_ID")))){
  SPECS$TASK <- rep(as.numeric(Sys.getenv("SLURM_ARRAY_TASK_MIN")):as.numeric(Sys.getenv("SLURM_ARRAY_TASK_MAX")), length=nrow(SPECS))
  SPECS <- SPECS[SPECS$TASK %in% as.numeric(Sys.getenv("SLURM_ARRAY_TASK_ID")),]
}

rm(optimal_gw,SPECS_crop,SPECS_period)
#-----------------------------------------------
# BOOTS                                      ####

lapply(
  c(1:nrow(SPECS)),
  function(spec){
    tryCatch({
      # spec <- 1
      #-----------------------------------------------
      # Exposure                                   ####
      exposure <- data.table::rbindlist(
        lapply(
          list.files(study_environment$wd$boots,recursive = T,full.names = T,pattern = paste0("spec",stringr::str_pad(SPECS$specN[spec],pad="0",2))),
          function(file){
            tryCatch({
              return(readRDS(file)$exposure)
            }, error = function(e){return(NULL)})
          }), fill = TRUE)
      setDT(exposure)
      exposure[, obs := ifelse(exp %in% NA,0,1)]
      exposure0 <- exposure[boot %in% "0000", c("p","theta","longlat","DistName","kernel","crop","period","specN","climate_base","Temp","exp")]
      exposure <- exposure[!boot %in% "0000", ][
        , .(exp_mean = mean(exp, na.rm = TRUE),
            exp_sd = sd(exp, na.rm = TRUE),
            exp_n = sum(obs, na.rm = TRUE)), by = c("p","theta","longlat","DistName","kernel","crop","period","specN","climate_base","Temp")]
      exposure <- as.data.frame(exposure0[exposure, on=c("p","theta","longlat","DistName","kernel","crop","period","specN","climate_base","Temp"), nomatch=0])

      #-----------------------------------------------
      # Piecewise linear function                  ####
      piecewise <- data.table::rbindlist(
        lapply(
          list.files(study_environment$wd$boots,recursive = T,full.names = T,pattern = paste0("spec",stringr::str_pad(SPECS$specN[spec],pad="0",2))),
          function(file){
            tryCatch({
              return(readRDS(file)$piecewise)
            }, error = function(e){return(NULL)})
          }), fill = TRUE)

      setDT(piecewise)
      piecewise[, obs := ifelse(Estimate %in% NA,0,1)]
      piecewise0 <- piecewise[boot %in% "0000", c("p","theta","longlat","DistName","kernel","crop","period",
                                                  "specN","climate_base","name","Estimate","StdError","t_value","p_value")]
      piecewise <- piecewise[!boot %in% "0000", ][
        , .(Estimate_mean = mean(Estimate, na.rm = TRUE),
            Estimate_sd = sd(Estimate, na.rm = TRUE),
            Estimate_n = sum(obs, na.rm = TRUE)), by = c("p","theta","longlat","DistName","kernel","crop","period","specN","climate_base","name")]
      piecewise <- as.data.frame(piecewise0[piecewise,
                                            on=c("p","theta","longlat","DistName","kernel","crop","period","specN","climate_base","name"), nomatch=0])

      #-----------------------------------------------
      # Nonlinear Relation                         ####
      relation <- data.table::rbindlist(
        lapply(
          list.files(study_environment$wd$boots,recursive = T,full.names = T,pattern = paste0("spec",stringr::str_pad(SPECS$specN[spec],pad="0",2))),
          function(file){
            tryCatch({
              return(readRDS(file)$relation)
            }, error = function(e){return(NULL)})
          }), fill = TRUE)

      setDT(relation)
      relation[, obs := ifelse(Piece %in% NA,0,1)]
      relation0 <- relation[boot %in% "0000", c("p","theta","longlat","DistName","kernel","crop","period","specN","climate_base","Temp","Piece","PieceSE")]
      relation <- relation[!boot %in% "0000", ][
        , .(Piece_mean = mean(Piece, na.rm = TRUE),
            Piece_sd = sd(Piece, na.rm = TRUE),
            Piece_n = sum(obs, na.rm = TRUE)), by = c("p","theta","longlat","DistName","kernel","crop","period","specN","climate_base","Temp")]
      relation <- as.data.frame(relation0[relation, on=c("p","theta","longlat","DistName","kernel","crop","period","specN","climate_base","Temp"), nomatch=0])

      #-----------------------------------------------
      # yield impact                               ####
      impact_yield <- data.table::rbindlist(
        lapply(
          list.files(study_environment$wd$boots,recursive = T,full.names = T,pattern = paste0("spec",stringr::str_pad(SPECS$specN[spec],pad="0",2))),
          function(file){
            tryCatch({
              return(readRDS(file)$impact_yield)
            }, error = function(e){return(NULL)})
          }), fill = TRUE)

      impact_yield <- as.data.frame(impact_yield)
      impact_yield <- impact_yield[c("boot","p","theta","longlat","DistName","kernel","crop","period","specN","climate_base","region","state_code","county_code","fip",
                                     "cc05","cc10","cc15","cc20","cc25","cc30")]

      impact_yield <- impact_yield |> tidyr::gather(warming_scenario, Estimate, c("cc05","cc10","cc15","cc20","cc25","cc30"))
      impact_yield$warming_scenario <- as.numeric(gsub("cc","",impact_yield$warming_scenario))/10
      setDT(impact_yield)
      impact_yield[, obs := ifelse(Estimate %in% NA,0,1)]
      impact_yield0 <- impact_yield[boot %in% "0000", c("p","theta","longlat","DistName","kernel","crop","period","specN",
                                                        "climate_base","region","state_code","county_code","fip","warming_scenario","Estimate")]
      impact_yield <- impact_yield[!boot %in% "0000", ][
        , .(Estimate_mean = mean(Estimate, na.rm = TRUE),
            Estimate_sd = sd(Estimate, na.rm = TRUE),
            Estimate_n = sum(obs, na.rm = TRUE)), by = c("p","theta","longlat","DistName","kernel","crop","period","specN",
                                                         "climate_base","region","state_code","county_code","fip","warming_scenario")]
      impact_yield <- as.data.frame(impact_yield0[impact_yield, on=c("p","theta","longlat","DistName","kernel","crop","period",
                                                                     "specN","climate_base","region","state_code","county_code","fip","warming_scenario"), nomatch=0])

      #-----------------------------------------------
      # Alfalfa availability                       ####
      availability <- data.table::rbindlist(
        lapply(
          list.files(study_environment$wd$boots,recursive = T,full.names = T,pattern = paste0("spec",stringr::str_pad(SPECS$specN[spec],pad="0",2))),
          function(file){
            tryCatch({
              return(readRDS(file)$availability)
            }, error = function(e){return(NULL)})
          }), fill = TRUE)

      setDT(availability)

      availability[, avail05 := ((avail05/avail00)-1)*100]
      availability[, avail10 := ((avail10/avail00)-1)*100]
      availability[, avail15 := ((avail15/avail00)-1)*100]
      availability[, avail20 := ((avail20/avail00)-1)*100]
      availability[, avail25 := ((avail25/avail00)-1)*100]
      availability[, avail30 := ((avail30/avail00)-1)*100]

      availability <- as.data.frame(availability)
      availability <- availability[c("boot","p","theta","longlat","DistName","kernel","crop","period","specN","climate_base","fip",
                                   "avail00","avail05","avail10","avail15","avail20","avail25","avail30")]

      availability <- availability |>  tidyr::gather(warming_scenario, Estimate, c("avail00","avail05","avail10","avail15","avail20","avail25","avail30"))
      availability$warming_scenario <- as.numeric(gsub("[^0-9]","",availability$warming_scenario))/10
      setDT(availability)
      availability[, obs := ifelse(Estimate %in% NA,0,1)]
      availability0 <- availability[boot %in% "0000", c("p","theta","longlat","DistName","kernel","crop","period","specN","climate_base","fip","warming_scenario","Estimate")]
      availability <- availability[!boot %in% "0000", ][
        , .(Estimate_mean = mean(Estimate, na.rm = TRUE),
            Estimate_sd = sd(Estimate, na.rm = TRUE),
            Estimate_n = sum(obs, na.rm = TRUE)), by = c("p","theta","longlat","DistName","kernel","crop","period","specN","climate_base","fip","warming_scenario")]
      availability <- as.data.frame(availability0[availability, on=c("p","theta","longlat","DistName","kernel","crop","period","specN","climate_base","fip","warming_scenario"), nomatch=0])

      #-----------------------------------------------
      # Associations                               ####
      associations <- data.table::rbindlist(
        lapply(
          list.files(study_environment$wd$boots,recursive = T,full.names = T,pattern = paste0("spec",stringr::str_pad(SPECS$specN[spec],pad="0",2))),
          function(file){
            tryCatch({
              return(readRDS(file)$associations)
            }, error = function(e){return(NULL)})
          }), fill = TRUE)

      setDT(associations)

      associations[, obs := ifelse(est %in% NA,0,1)]
      associations0 <- associations[boot %in% "0000", c("p","theta","longlat","DistName","kernel","crop","period","specN","climate_base","fip","name","est","se", "tv","pv")]
      associations <- associations[!boot %in% "0000", ][
        , .(est_mean = mean(est, na.rm = TRUE),
            est_sd = sd(est, na.rm = TRUE),
            est_n = sum(obs, na.rm = TRUE)), by = c("p","theta","longlat","DistName","kernel","crop","period","specN","climate_base","fip","name")]
      associations <- as.data.frame(associations0[associations, on=c("p","theta","longlat","DistName","kernel","crop","period","specN","climate_base","fip","name"), nomatch=0])
      #-----------------------------------------------
      # Cattle shifts                              ####
      assoc <- data.table::rbindlist(
        lapply(
          list.files(study_environment$wd$boots,recursive = T,full.names = T,pattern = paste0("spec",stringr::str_pad(SPECS$specN[spec],pad="0",2))),
          function(file){
            tryCatch({
              return(readRDS(file)$associations)
            }, error = function(e){return(NULL)})
          }), fill = TRUE)

      assoc <- as.data.frame(assoc[,c("boot","p","theta","longlat","DistName","kernel","crop","period","specN","climate_base","fip","name","est")])
      assoc$name <- paste0("b_",assoc$name)
      assoc <- assoc |> tidyr::spread(name, est)
      setDT(assoc)

      cattle <- data.table::rbindlist(
        lapply(
          list.files(study_environment$wd$boots,recursive = T,full.names = T,pattern = paste0("spec",stringr::str_pad(SPECS$specN[spec],pad="0",2))),
          function(file){
            tryCatch({
              return(readRDS(file)$availability)
            }, error = function(e){return(NULL)})
          }), fill = TRUE)

      setDT(cattle)

      cattle <- cattle[assoc, on=c("boot","p","theta","longlat","DistName","kernel","crop","period","specN","climate_base","fip"), nomatch=0]

      cattle[, cattleA05 := ((avail05/avail00)-1)*100*b_avail00]
      cattle[, cattleA10 := ((avail10/avail00)-1)*100*b_avail00]
      cattle[, cattleA15 := ((avail15/avail00)-1)*100*b_avail00]
      cattle[, cattleA20 := ((avail20/avail00)-1)*100*b_avail00]
      cattle[, cattleA25 := ((avail25/avail00)-1)*100*b_avail00]
      cattle[, cattleA30 := ((avail30/avail00)-1)*100*b_avail00]

      cattle[, cattleB05 := ((prod05/prod00)-1)*100*b_prod00]
      cattle[, cattleB10 := ((prod10/prod00)-1)*100*b_prod00]
      cattle[, cattleB15 := ((prod15/prod00)-1)*100*b_prod00]
      cattle[, cattleB20 := ((prod20/prod00)-1)*100*b_prod00]
      cattle[, cattleB25 := ((prod25/prod00)-1)*100*b_prod00]
      cattle[, cattleB30 := ((prod30/prod00)-1)*100*b_prod00]

      cattle[, cattleC05 := ((prod05/prod00)-1)*100*b_prod00 + ((prod05_LM/prod00_LM)-1)*100*b_prod00_LM]
      cattle[, cattleC10 := ((prod10/prod00)-1)*100*b_prod00 + ((prod10_LM/prod00_LM)-1)*100*b_prod00_LM]
      cattle[, cattleC15 := ((prod15/prod00)-1)*100*b_prod00 + ((prod15_LM/prod00_LM)-1)*100*b_prod00_LM]
      cattle[, cattleC20 := ((prod20/prod00)-1)*100*b_prod00 + ((prod20_LM/prod00_LM)-1)*100*b_prod00_LM]
      cattle[, cattleC25 := ((prod25/prod00)-1)*100*b_prod00 + ((prod25_LM/prod00_LM)-1)*100*b_prod00_LM]
      cattle[, cattleC30 := ((prod30/prod00)-1)*100*b_prod00 + ((prod30_LM/prod00_LM)-1)*100*b_prod00_LM]

      cattle <- as.data.frame(cattle)
      cattle <- cattle[c("boot","p","theta","longlat","DistName","kernel","crop","period","specN","climate_base","fip",
                         names(cattle)[grepl("cattle",names(cattle))])]
      cattle <- cattle |>  tidyr::gather(warming_scenario, Estimate, names(cattle)[grepl("cattle",names(cattle))])
      cattle$cattle <- gsub("[0-9]","",cattle$warming_scenario)
      cattle$warming_scenario <- as.numeric(gsub("[^0-9]","",cattle$warming_scenario))/10
      cattle <- cattle |>  tidyr::spread(cattle, Estimate)

      setDT(cattle)
      cattle[, obs := ifelse(cattleA %in% NA,0,1)]
      cattle0 <- cattle[boot %in% "0000", c("p","theta","longlat","DistName","kernel","crop","period","specN","climate_base","fip","warming_scenario","cattleA","cattleB","cattleC")]
      cattle  <- cattle[!boot %in% "0000", ][
        , .(cattleA_mean = mean(cattleA, na.rm = TRUE),cattleA_sd = sd(cattleA, na.rm = TRUE),
            cattleB_mean = mean(cattleB, na.rm = TRUE),cattleB_sd = sd(cattleB, na.rm = TRUE),
            cattleC_mean = mean(cattleC, na.rm = TRUE),cattleC_sd = sd(cattleC, na.rm = TRUE),
            cattle_n = sum(obs, na.rm = TRUE)), by = c("p","theta","longlat","DistName","kernel","crop","period","specN","climate_base","fip","warming_scenario")]
      cattle <- as.data.frame(cattle0[cattle, on=c("p","theta","longlat","DistName","kernel","crop","period","specN","climate_base","fip","warming_scenario"), nomatch=0])

      #-----------------------------------------------
      # BW                                         ####
      bw <- data.table::rbindlist(
        lapply(
          list.files(study_environment$wd$boots,recursive = T,full.names = T,pattern = paste0("spec",stringr::str_pad(SPECS$specN[spec],pad="0",2))),
          function(file){
            tryCatch({
              return(readRDS(file)$bw)
            }, error = function(e){return(NULL)})
          }), fill = TRUE)

      setDT(bw)

      bw[, obs := ifelse(bw_opt %in% NA,0,1)]
      bw0 <- bw[boot %in% "0000", c("p","theta","longlat","DistName","kernel","crop","period","specN","bw_opt")]
      bw <- bw[!boot %in% "0000", ][
        , .(bw_opt_mean = mean(bw_opt, na.rm = TRUE),
            bw_opt_sd = sd(bw_opt, na.rm = TRUE),
            bw_opt_n = sum(obs, na.rm = TRUE)), by = c("p","theta","longlat","DistName","kernel","crop","period","specN")]
      bw <- as.data.frame(bw0[bw, on=c("p","theta","longlat","DistName","kernel","crop","period","specN"), nomatch=0])
      #-----------------------------------------------
      # Finalize                                   ####

      res <- list(
        exposure=as.data.frame(exposure),
        piecewise = as.data.frame(piecewise),
        relation = as.data.frame(relation),
        impact_yield = as.data.frame(impact_yield),
        availability = as.data.frame(availability),
        associations = as.data.frame(associations),
        cattle = as.data.frame(cattle),
        bw= as.data.frame(bw)
      )

      saveRDS(res,file=file.path(study_environment$wd$summary,
                                 paste0("spec",stringr::str_pad(SPECS$specN[spec],pad="0",2),"_",
                                        SPECS$crop[spec],"_period",SPECS$period[spec],"_",SPECS$NAME[spec],".rds")))
      #-----------------------------------------------
      return(spec)
    }, error = function(e){return(NULL)})
  })

function(){

  lapply(
    names(readRDS(list.files(study_environment$wd$summary,full.names = T)[1])),
    function(xx){
      tryCatch({

        res <- data.table::rbindlist(
          lapply(
            list.files(study_environment$wd$summary,full.names = T),
            function(file){
              tryCatch({
                return(readRDS(file)[[xx]])
              }, error = function(e){return(NULL)})
            }), fill = TRUE)

        saveRDS(res,file=paste0("output/summary/summary_",xx,".rds"))

        return(xx)
      }, error = function(e){return(NULL)})
    })

}



#-----------------------------------------------
