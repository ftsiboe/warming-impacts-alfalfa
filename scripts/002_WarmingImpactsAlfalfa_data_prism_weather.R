#SBATCH --job-name=prism_county
#SBATCH --time=0-10:00:00
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=1
#SBATCH --mem=15G
#SBATCH --array=1-205
#SBATCH --mail-user=ftsiboe@ksu.edu
#SBATCH --mail-type=ALL
#SBATCH --output=/dev/null

rm(list=ls(all=TRUE));library(future.apply);library(data.table)

study_environment <- readRDS("data/study_environment.rds")

devtools::document(file.path(dirname(dirname(getwd())),"packages/rAgroClimate"))

prism_list <- data.frame(file=list.files(study_environment$prism_archive, recursive = TRUE, full.names = TRUE))
prism_list$date <- as.Date(gsub("prism_daily_all_4km_","",basename(prism_list$file)),format = "%Y%m%d")
prism_list$year <- as.numeric(format(prism_list$date ,format = "%Y"))
year_list <- unique(prism_list$year)
year_list <- max(year_list):min(year_list)

# plan(multisession)

# County Aggregation
countySpatialWeights <- readRDS(study_environment$county_weights_file_path)
countySpatialWeights <- countySpatialWeights[commodity_code %in% c(107,332)]
countySpatialWeights <- countySpatialWeights[!gridNumber_prism %in% c(NA,Inf,-Inf,NaN)]
countySpatialWeights <- countySpatialWeights[
  , .(weight = sum(weight,na.rm=TRUE)),
  by=c("state_code","county_code","gridNumber_prism")
];gc()

countySpatialWeights[
  , county_fips := as.character(paste0(stringr::str_pad(as.numeric(as.character(state_code)), 2, pad = "0"),
                               stringr::str_pad(as.numeric(as.character(county_code)), 3, pad = "0")))]

function(){

  state_fips <- unique(as.data.table(tigris::fips_codes)[, .(state_code, state)])
  state_fips[, state_code := as.integer(state_code)]

  ARRAY <- unique(countySpatialWeights$state_code)
  ARRAY <- ARRAY[ARRAY %in% unique(readRDS("data/nass_hay_production.rds")$state_code)]
  ARRAY <- data.table::rbindlist(
    lapply(
      year_list,
      function(year){
        dir.create(file.path(study_environment$wd$prism_weather,year))
        ARRAY <- data.frame(
          data.table::rbindlist(
            lapply(
              ARRAY,
              function(st){
                # year <- 2021 ; state <- 1
                ARRAY <- data.frame(year=year,state_code=st,state_ab=state_fips[match(st, state_code), state])
                ARRAY$name <- paste0("prism_",ARRAY$year,"_",ARRAY$state_ab,".rds")
                ARRAY$done <- ARRAY$name %in% list.files(path=file.path(study_environment$wd$prism_weather,year))
                return(ARRAY)
              }),fill = T))
        return(ARRAY)
      }),fill = T)

  ARRAY <- ARRAY[ARRAY$done %in% F,] # 1136/2911
  ARRAY <- ARRAY[order(-ARRAY$year),]
  # ARRAY <- ARRAY[state_ab %in% "KS"]

  ARRAY <- ARRAY[! name %in% basename(list.files(path=study_environment$wd$prism_weather, recursive = TRUE))]

  saveRDS(ARRAY,file="output/weather_ARRAY.rds")

}


if(Sys.getenv("SLURM_JOB_NAME") %in% "prism_wther1") {

  ARRAY <- readRDS("output/weather_ARRAY.rds")
  if(!is.na(as.numeric(Sys.getenv("SLURM_ARRAY_TASK_ID")))){
    ARRAY$TASK <- rep(as.numeric(Sys.getenv("SLURM_ARRAY_TASK_MIN")):as.numeric(Sys.getenv("SLURM_ARRAY_TASK_MAX")), length=nrow(ARRAY))
    ARRAY <- ARRAY[ARRAY$TASK %in% as.numeric(Sys.getenv("SLURM_ARRAY_TASK_ID")),]
  }

  lapply(
    1:nrow(ARRAY),
    function(task){
      # ARRAY <- ARRAY[ARRAY$state_ab %in% "KS",]; task <- 1
      year <- ARRAY$year[task]
      state <- ARRAY$state_code[task]
      name <- ARRAY$name[task]
      fiplist <- unique(countySpatialWeights[state_code %in% state][["county_fips"]])

      if(! name %in% list.files(path=file.path(study_environment$wd$prism_weather,year))){
        tryCatch({

          res <- data.table::rbindlist(
            lapply(
              fiplist,
              function(year,state,fip){
                cat(crayon::black("*** Cleaning",state,fip,year,Sys.time()),fill=T)
                tryCatch({
                  # fip <- fiplist[1]
                  spatial_cd <- countySpatialWeights[county_fips %in% fip]

                  if(nrow(spatial_cd)>0){
                    # unique(Final$crop)

                    # Determine when to aggregate
                    SATES <- as.data.frame(terra::vect(file.path(study_environment$usaPolygons_archive,"USA_States.shp")))

                    # June 1 to May 31 for all other States.
                    temporal_list <- as.Date(as.Date(paste0(year,"-06-01")):as.Date(paste0(year+1,"-05-31")),origin="1970-01-01")

                    # April 1 to March 31 for Arizona and California
                    state_cd <- as.numeric(SATES[toupper(SATES$NAME) %in% toupper(c("Arizona", "California")),"STATEFP"])
                    if(spatial_cd$state_code[1] %in% state_cd){
                      temporal_list <- as.Date(as.Date(paste0(year,"-04-01")):as.Date(paste0(year+1,"-03-31")),origin="1970-01-01")
                    }

                    # May 1 to April 30
                    S2 <- toupper(c("Kansas", "Kentucky", "Louisiana", "Mississippi", "Missouri", "Nevada", "New Mexico",
                                    "North Carolina", "Oklahoma", "Pennsylvania", "South Carolina", "Tennessee", "Texas", "Utah", "Virginia"))
                    # S2 <- S2[!S2 %in% toupper(SATES$NAME)]
                    state_cd <- as.numeric(SATES[toupper(SATES$NAME) %in% toupper(S2),"STATEFP"])
                    if(spatial_cd$state_code[1] %in% state_cd){
                      temporal_list <- as.Date(as.Date(paste0(year,"-05-01")):as.Date(paste0(year+1,"-04-30")),origin="1970-01-01")
                    }

                    # Determine where to aggregate
                    relevant_prism <- unique(prism_list[prism_list$date %in% temporal_list,])
                    relevant_prism <- data.table::rbindlist(
                      lapply(
                        1:nrow(relevant_prism),
                        function(nm){
                          tryCatch({
                            # nm <- 1
                            data <- readRDS(paste0(relevant_prism$file[nm]))
                            data <- data[gridNumber_prism %in% spatial_cd$gridNumber_prism]
                            return(data)
                          }, error=function(e){NULL})
                        }),fill = T)

                    # Determine where to aggregate for each crop
                    relevant_prism <- dplyr::inner_join(relevant_prism,spatial_cd,by=c("gridNumber_prism"))

                    relevant_prism <- rbind(
                      data.frame(relevant_prism,period=0),
                      data.frame(
                        data.table::rbindlist(
                          lapply(
                            unique(as.integer(as.factor(as.numeric(format(relevant_prism$date,"%Y%m"))))),
                            function(mth){
                              data <- relevant_prism
                              data$period <- as.integer(as.factor(as.numeric(format(relevant_prism$date,"%Y%m"))))
                              data <- data[data$period<=mth,]
                              data$period <- 100+mth
                              return(data)
                            }),fill = T)),
                      data.frame(relevant_prism,period=as.integer(as.factor(as.numeric(format(relevant_prism$date,"%Y%m"))))))

                    rm(spatial_cd,temporal_list);gc()

                    # --- Aggregate weather into exposure / degree-day metrics
                    relevant_prism <- aggregate_weather_variables(
                      data                       = relevant_prism,
                      date_col                   = "date",
                      precipitation_col          = "ppt",
                      temperature_minimum_col    = "tmin",
                      temperature_maximum_col    = "tmax",
                      vapor_pressure_minimum_col = "vpdmin",
                      vapor_pressure_maximum_col = "vpdmax",
                      weight_variable_col        = "weight",
                      identifiers                = intersect(c(
                        "state_code","county_code","county_fips","commodity_code",
                        "period","dataset_type", "network_type", "resolution"
                      ),
                      names(relevant_prism)
                      ),
                      warming_scenario = 0
                    )

                    relevant_prism <- relevant_prism[, c(names(relevant_prism)[!(grepl("ddayN",names(relevant_prism)) | grepl("expN",names(relevant_prism)))]), with = FALSE]

                  }

                  return(relevant_prism)
                }, error = function(e){return(NULL)})
              },year=year,state=state),fill = T)

          res[, commodity_year := year ]

          saveRDS(res, file = file.path(study_environment$wd$prism_weather,year, name))

          rm(res);gc()

        }, error=function(e){})
      }
      return(task)})
}

if(Sys.getenv("SLURM_JOB_NAME") %in% "prism_wther2") {


  dir.create(file.path("data/prism_weather"))

  lapply(
    list.files(study_environment$wd$prism_weather),
    function(year){
      tryCatch({

        df <- data.table::rbindlist(
          lapply(
            list.files(file.path(study_environment$wd$prism_weather,year),recursive = T,full.names = T),
            function(file){
              tryCatch({
                # file <- list.files(study_environment$wd$prism_weather,recursive = T,full.names = T)[1]
                data <- readRDS(file)
                data[,warming_scenario := NULL]
                return(data)
              }, error=function(e){NULL})
            }), fill = TRUE)

        saveRDS(df,file=file.path("data/prism_weather",paste0("prism_weather_",year,".rds")))

        invisible()

      }, error=function(e){NULL})
    })

}



