# Hard reset of workspace
rm(list = ls(all = TRUE));gc()
library(data.table)
#----------------------------------------------------
# Initialize environment                          ####
invisible(lapply(list.files("scripts/helpers", pattern = "[.]R$", full.names = TRUE), source))
study_environment <- setup_environment(
  year_beg = 1918, year_end = 2024, seed = 1980632,
  project_name="WarmingImpactsAlfalfa",
  local_directories = list(
    file.path("data"),
    file.path("output","summary"),
    file.path("output","exhibits","figure_data"),
    file.path("output","exhibits"),
    file.path("output","bootstraps"),
    file.path("output","releases")
  ),
  fastscratch_directories = NULL)


sysname <- tolower(as.character(Sys.info()[["sysname"]]))
user    <- Sys.info()[["user"]]

# Fast-scratch working directories

if(grepl("windows", sysname)){
  fastscratch_directory <- "fastscratch"
}else{
  fastscratch_directory <- file.path("/fastscratch", user,"warming-impacts-alfalfa")
}

study_environment$wd <- list(
  prism_weather = file.path(fastscratch_directory,"prism_weather"),
  prism_climate = file.path(fastscratch_directory,"prism_climate"),
  knots         = file.path(fastscratch_directory,"knots"),
  boots         = file.path(fastscratch_directory,"boots"),
  summary       = file.path(fastscratch_directory,"summary")
)

invisible(lapply(study_environment$wd, dir.create, recursive = TRUE, showWarnings = FALSE))

temporary_dir <- tempdir()


if(grepl("TSIB",toupper(user))){
  if(grepl("windows", sysname)){
    nass_database            <- file.path("C:/Users",user,"Dropbox (Personal)/database/usa/usda_nass")
    prism_archive            <- file.path("C:/Users",user,"Dropbox (Personal)/database/spatialData/archive/prism")
    gssurgo_archive          <- file.path("C:/Users",user,"Dropbox (Personal)/database/spatialData/archive/usda_nrcs_gssurgo/gSSURGO_CONUS_2026/FY2026_gSSURGO_mukey_grid")
    usaPolygons_archive      <- file.path("C:/Users",user,"Dropbox (Personal)/database/spatialData/archive/usaPolygons")
    county_weights_file_path <- file.path("C:/Users",user,"Dropbox (Personal)/database/spatialData/output/countySpatialWeights.rds")
  }else{
    nass_database            <- file.path("/homes",user,"database/usa/usda_nass")
    prism_archive            <- file.path("/homes",user,"database/spatialData/archive/prism")
    gssurgo_archive          <- file.path("/homes",user,"database/spatialData/archive/usda_nrcs_gssurgo/gSSURGO_CONUS_2026/FY2026_gSSURGO_mukey_grid")
    usaPolygons_archive      <- file.path("/homes",user,"database/spatialData/archive/usaPolygons")
    county_weights_file_path <- file.path("/homes",user,"database/spatialData/output/countySpatialWeights.rds")
  }
}

study_environment$nass_database  <- nass_database
study_environment$prism_archive <- prism_archive
study_environment$gssurgo_archive <- gssurgo_archive
study_environment$usaPolygons_archive <- usaPolygons_archive
study_environment$county_weights_file_path <- county_weights_file_path
saveRDS(study_environment,file ="data/study_environment.rds")

Keep.List<-c("Keep.List",ls())
#----------------------------------------------------
# Season Length                                   ####
rm(list= ls()[!(ls() %in% c(Keep.List))]);gc()
piggyback::pb_download(
  file = "nassSeasonLengthState.rds",
  dest = temporary_dir,
  repo = "ftsiboe/USFarmSafetyNetLab",
  tag  = "nass_extracts",
  overwrite = TRUE)
data <- readRDS(file.path(temporary_dir,"nassSeasonLengthState.rds"))
data <- data[grepl("HAY",commodity_name)]
data[, data_source :=NULL]
saveRDS(data,file="data/nass_hay_season_length.rds")
#----------------------------------------------------
# Hay Production                                  ####
rm(list= ls()[!(ls() %in% c(Keep.List))])
piggyback::pb_download(
  file = "nassSurveyHayProduction.rds",
  dest = temporary_dir,
  repo = "ftsiboe/USFarmSafetyNetLab",
  tag  = "nass_extracts",
  overwrite = TRUE)
data <- readRDS(file.path(temporary_dir,"nassSurveyHayProduction.rds"))
data[, data_source :=NULL]
saveRDS(data,file="data/nass_hay_production.rds")
#----------------------------------------------------
# ANIMAL TOTALS                                   ####
rm(list= ls()[!(ls() %in% c(Keep.List))])
piggyback::pb_download(
  file = "nassSurveyAnimalInventory.rds",
  dest = temporary_dir,
  repo = "ftsiboe/USFarmSafetyNetLab",
  tag  = "nass_extracts",
  overwrite = TRUE)
data <- readRDS(file.path(temporary_dir,"nassSurveyAnimalInventory.rds"))
data[, data_source :=NULL]
saveRDS(data,file="data/nass_animal_inventory.rds")
#----------------------------------------------------
# Census                                          ####
rm(list= ls()[!(ls() %in% c(Keep.List))])

dir_fastscratch  <- file.path(nass_database,"quick_stats_database","fastscratch")

if(!dir.exists(dir_fastscratch)) {
  dir.create(dir_fastscratch, recursive = TRUE)
}

downloaded_nass_large_datasets(
  large_datasets = c(
    paste0("census", c(2022, 2017, 2012, 2007, 2002))),
  dir_dest = dir_fastscratch)

data <- data.table::rbindlist(
  lapply(
    c("qs.census2022.txt.gz","qs.census2017.txt.gz","qs.census2007.txt.gz","qs.census2002.txt.gz"),
    function(xx){
      tryCatch({
        # xx <- "qs.census2022.txt.gz"
        data <- data.table::fread(list.files(dir_fastscratch,pattern = xx,full.names = T))
        data <- as.data.frame(data[data$SHORT_DESC %in% c("HAY, ALFALFA - PRODUCTION, MEASURED IN TONS",
                                                          "HAY, ALFALFA - ACRES HARVESTED",
                                                          "CATTLE, INCL CALVES - INVENTORY"),])
        data <- data[data$AGG_LEVEL_DESC %in% c("COUNTY"),]
        data <- data[data$DOMAIN_DESC %in% c("TOTAL"),]
        names(data) <- tolower(names(data))
        data$value <- as.numeric(gsub(",","",as.character(data$value)))
        data <- data[! data$value %in% c(0,NA,Inf,-Inf,NaN),]
        data <- data[!data$county_name %in% c("OTHER (COMBINED) COUNTIES"),]
        data <- doBy::summaryBy(value~year+state_fips_code+state_alpha+asd_code+county_code+short_desc,FUN=mean,na.rm=T,data=data,keep.names = T)
        data <- doBy::summaryBy(value~year+state_fips_code+state_alpha+county_code+short_desc,FUN=sum,na.rm=T,data=data)
        names(data) <- c("census_year","state_code","state_ab","county_code","short_desc","value")
        gc()
        return(data)
      }, error = function(e){return(NULL)})
    }),fill = T)
saveRDS(data,file="data/nass_census.rds")
#----------------------------------------------------
# Spatial Rep                                     ####
rm(list= ls()[!(ls() %in% c(Keep.List))])

Cattle<- readRDS("data/nass_animal_inventory.rds")
Cattle  <- Cattle[commodity_year %in% 2002:2022]
Cattle[,inventory := cattle]
Cattle<-doBy::summaryBy(inventory~state_code+county_code,data=Cattle,FUN=mean,keep.names = T,na.rm=T)
Cattle$state_code <- as.numeric(as.character(Cattle$state_code))
Cattle$county_code <- as.numeric(as.character(Cattle$county_code))

Alfalfa <- readRDS("data/nass_hay_production.rds")
Alfalfa <- Alfalfa[commodity_name %in% "hay_alfalfa"]
Alfalfa<-Alfalfa[commodity_year %in% 2002:2022]
Alfalfa<-doBy::summaryBy(production+area+yield~state_code+county_code,data=Alfalfa,FUN=mean,keep.names = T,na.rm=T)

survey<-dplyr::full_join(Alfalfa,Cattle,by=c("state_code","county_code"))
survey <- as.data.frame(survey)[c("county_code","state_code","area","production","yield","inventory")]
survey$area <- survey$area/1000
survey$production <- survey$production/1000
survey$inventory <- survey$inventory/1000
survey <- survey |>  tidyr::gather(Type, survey, c("area","production","yield","inventory"))
survey <- survey[!survey$survey %in% c(NaN,NA,-Inf,Inf),]

census <- readRDS("data/nass_census.rds")
census <- census[census$census_year %in% 2002:2022,]
census$short_desc <- as.character(factor(census$short_desc,
                                         levels = c("CATTLE, INCL CALVES - INVENTORY","HAY, ALFALFA - ACRES HARVESTED","HAY, ALFALFA - PRODUCTION, MEASURED IN TONS"),
                                         labels = c("inventory","area","production")))
census <- census |> tidyr::spread(short_desc, value)
census$yield <- census$production/census$area
census$area <- census$area/1000
census$production <- census$production/1000
census$inventory <- census$inventory/1000
census<-doBy::summaryBy(production+area+yield+inventory~state_code+county_code,data=census,FUN=mean,keep.names = T,na.rm=T)
census <- census |>  tidyr::gather(Type, census, c("area","production","yield","inventory"))
census <- census[!census$census %in% c(NaN,NA,-Inf,Inf),]

Plot.data <- dplyr::full_join(census,survey,by=c("state_code","county_code","Type"))

Plot.data$fip<-as.character(paste0(stringr::str_pad(as.numeric(as.character(Plot.data$state_code)), 2, pad = "0"),
                                   stringr::str_pad(as.numeric(as.character(Plot.data$county_code)), 3, pad = "0")))
Plot.data$Value <- ifelse(Plot.data$census %in% c(NaN,NA,-Inf,Inf,0),Plot.data$survey ,Plot.data$census )
Plot.data <- Plot.data[c("state_code","county_code","Type","fip","Value")]
Plot.data <- Plot.data |> tidyr::spread(Type, Value)
for(i in 1:10){
  Plot.data$area <- ifelse(Plot.data$area %in% c(NaN,NA,-Inf,Inf,0),Plot.data$production/Plot.data$yield,Plot.data$area)
  Plot.data$yield <- ifelse(Plot.data$yield %in% c(NaN,NA,-Inf,Inf,0),Plot.data$production/Plot.data$area,Plot.data$yield)
  Plot.data$production <- ifelse(Plot.data$production %in% c(NaN,NA,-Inf,Inf,0),Plot.data$yield*Plot.data$area,Plot.data$production)
}
Plot.data$area       <- ifelse(Plot.data$yield %in% c(NaN,NA,-Inf,Inf,0) | Plot.data$production %in% c(NaN,NA,-Inf,Inf,0) | Plot.data$area %in% c(NaN,NA,-Inf,Inf,0), NA,Plot.data$area)
Plot.data$yield      <- ifelse(Plot.data$yield %in% c(NaN,NA,-Inf,Inf,0) | Plot.data$production %in% c(NaN,NA,-Inf,Inf,0) | Plot.data$area %in% c(NaN,NA,-Inf,Inf,0), NA,Plot.data$yield)
Plot.data$production <- ifelse(Plot.data$yield %in% c(NaN,NA,-Inf,Inf,0) | Plot.data$production %in% c(NaN,NA,-Inf,Inf,0) | Plot.data$area %in% c(NaN,NA,-Inf,Inf,0), NA,Plot.data$production)
saveRDS(Plot.data,file="data/spatial_representation.rds")
#----------------------------------------------------

