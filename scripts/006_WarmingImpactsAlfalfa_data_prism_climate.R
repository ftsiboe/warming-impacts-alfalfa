#SBATCH --job-name=prism_county
#SBATCH --time=0-10:00:00
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=1
#SBATCH --mem=15G
#SBATCH --array=1-205
#SBATCH --mail-user=ftsiboe@ksu.edu
#SBATCH --mail-type=ALL
#SBATCH --output=/dev/null

rm(list=ls(all=TRUE));library(future.apply);library(data.table);library(terra)

study_environment <- readRDS("data/study_environment.rds")

# Rerun check
# Stage 2 (prism_climate2) trims data/prism_climate to the degree-day columns
# implied by the knots. If the knots have changed since those files were last
# written (national 004, or the county-level cluster knots from 005,
# optimal_knots_cluster), the trim is stale and stage 2 must be re-run.
# This block only REPORTS status - it changes nothing.
# NOTE: climate dday columns are written zero-padded to two digits
# (aggregate_weather_variables -> "dday00".."dday45"); the +/-band search can
# introduce single-digit thresholds, so any column lookup must pad likewise.
# The check below compares thresholds as INTEGERS, so it is padding-agnostic.
local({
  tryCatch({
    req <- 0L
    if (file.exists("output/optimal_knots.rds")) {
      nk <- as.data.frame(readRDS("output/optimal_knots.rds"))
      req <- c(req, as.integer(nk$Tmin), as.integer(nk$Tmax))
    }
    if (file.exists("output/optimal_knots_cluster.rds")) {
      ck <- as.data.frame(readRDS("output/optimal_knots_cluster.rds"))
      req <- c(req, as.integer(ck$Tmin), as.integer(ck$Tmax))
    }
    req <- sort(unique(req[is.finite(req) & req >= 0]))

    files <- list.files("data/prism_climate", pattern = "[.]rds$", full.names = TRUE)
    if (length(req) == 0) {
      message(">>> 005 rerun check: no knot files found - run 003 (and 004) first.")
    } else if (length(files) == 0) {
      message(">>> RERUN 006: data/prism_climate is empty - run stage 1 (prism_climate1) then stage 2 (prism_climate2).")
    } else {
      have <- names(readRDS(files[1]))
      have_thr <- suppressWarnings(as.integer(gsub("dday", "", grep("^dday[0-9]+$", have, value = TRUE))))
      missing <- setdiff(req, have_thr)
      if (length(missing) > 0) {
        message(">>> RERUN 006 stage 2 (prism_climate2): data/prism_climate is missing dday thresholds: ",
                paste0("dday", missing, collapse = ", "),
                "  (knots changed since it was last written).")
      } else {
        message("005 rerun check: data/prism_climate already carries every knot threshold (",
                paste0("dday", req, collapse = ", "), ") - no rerun needed.")
      }
    }
  }, error = function(e) message("005 rerun check skipped: ", conditionMessage(e)))
})

# Load rAgroClimate (sibling package) READ-ONLY. Do NOT devtools::document() here:
# with a large SLURM array every task would regenerate the shared package man/ +
# NAMESPACE concurrently (race/corruption) and it is slow. Document once, offline.
devtools::load_all(file.path(dirname(dirname(getwd())),"packages/rAgroClimate"), quiet = TRUE)

prism_list <- data.frame(file=list.files(study_environment$prism_archive, recursive = TRUE, full.names = TRUE))
prism_list$date <- as.Date(gsub("prism_daily_all_4km_","",basename(prism_list$file)),format = "%Y%m%d")
prism_list$year <- as.numeric(format(prism_list$date ,format = "%Y"))
year_list <- unique(prism_list$year)
year_list <- max(year_list):min(year_list)

# plan(multisession)
Counties <- vect(file.path(study_environment$usaPolygons_archive,"USA_Counties.shp"))
Counties$county_fips <-as.character(paste0(stringr::str_pad(as.numeric(as.character(Counties$STATEFP)), 2, pad = "0"),
                                           stringr::str_pad(as.numeric(as.character(Counties$COUNTYFP)), 3, pad = "0")))
Counties$state_code <- as.numeric(as.character(Counties$STATEFP))

gridNumber_prism  = get_raster_grid(source = "prism",  force = FALSE)$prismRaster

function(){

  state_fips <- unique(as.data.table(tigris::fips_codes)[, .(state_code, state)])
  state_fips[, state_code := as.integer(state_code)]

  ARRAY <- unique(state_fips$state_code)
  ARRAY <- data.table::rbindlist(
    lapply(
      year_list,
      function(year){
        dir.create(file.path(study_environment$wd$prism_climate,year))
        ARRAY <- data.frame(
          data.table::rbindlist(
            lapply(
              ARRAY,
              function(st){
                # year <- 2021 ; state <- 1
                ARRAY <- data.frame(year=year,state_code=st,state_ab=state_fips[match(st, state_code), state])
                ARRAY$name <- paste0("prism_",ARRAY$year,"_",ARRAY$state_ab,".rds")
                ARRAY$done <- ARRAY$name %in% list.files(path=file.path(study_environment$wd$prism_climate,year))
                return(ARRAY)
              }),fill = T))
        return(ARRAY)
      }),fill = T)

  ARRAY <- ARRAY[ARRAY$done %in% F,] # 1136/2911
  ARRAY <- ARRAY[order(-ARRAY$year),]
  # ARRAY <- ARRAY[state_ab %in% "KS"]
  saveRDS(ARRAY,file="output/climate_ARRAY.rds")

}

if(Sys.getenv("SLURM_JOB_NAME") %in% "prism_climate1") {

  ARRAY <- readRDS("output/climate_ARRAY.rds")
  if(!is.na(as.numeric(Sys.getenv("SLURM_ARRAY_TASK_ID")))){
    ARRAY$TASK <- rep(as.numeric(Sys.getenv("SLURM_ARRAY_TASK_MIN")):as.numeric(Sys.getenv("SLURM_ARRAY_TASK_MAX")), length=nrow(ARRAY))
    ARRAY <- ARRAY[ARRAY$TASK %in% as.numeric(Sys.getenv("SLURM_ARRAY_TASK_ID")),]
  }

  # TASK <- 5
  # plan(list(
  #   tweak(multisession, workers = 5),
  #   tweak(multisession, workers = 1)
  # ))

  plan(multisession)

  lapply(
    1:nrow(ARRAY),
    function(task){
      # ARRAY <- ARRAY[ARRAY$state_ab %in% "KS",]; task <- 1
      year <- ARRAY$year[task]
      state <- ARRAY$state_code[task]
      name <- ARRAY$name[task]
      fiplist <- unique(Counties$county_fips[Counties$state_code %in% state])

      if(! name %in% list.files(path=file.path(study_environment$wd$prism_climate,year))){
        tryCatch({

          res <- data.table::rbindlist(
            lapply(
              fiplist,
              function(year,state,fip){
                cat(crayon::black("*** Cleaning",state,fip,year,Sys.time()),fill=T)
                tryCatch({
                  # fip <- fiplist[1]
                  spatial_cd <- terra::project(Counties[Counties$county_fips %in% fip,], terra::crs(gridNumber_prism))
                  spatial_cd <- terra::crop(terra::mask(gridNumber_prism, spatial_cd), spatial_cd)

                  if(nrow(spatial_cd)>0){
                    # unique(Final$crop)

                    # Determine when to aggregate
                    SATES <- as.data.frame(terra::vect(file.path(study_environment$usaPolygons_archive,"USA_States.shp")))

                    # June 1 to May 31 for all other States.
                    temporal_list <- as.Date(as.Date(paste0(year,"-06-01")):as.Date(paste0(year+1,"-05-31")),origin="1970-01-01")

                    # April 1 to March 31 for Arizona and California
                    state_cd <- as.numeric(SATES[toupper(SATES$NAME) %in% toupper(c("Arizona", "California")),"STATEFP"])
                    if(state %in% state_cd){
                      temporal_list <- as.Date(as.Date(paste0(year,"-04-01")):as.Date(paste0(year+1,"-03-31")),origin="1970-01-01")
                    }

                    # May 1 to April 30
                    S2 <- toupper(c("Kansas", "Kentucky", "Louisiana", "Mississippi", "Missouri", "Nevada", "New Mexico",
                                    "North Carolina", "Oklahoma", "Pennsylvania", "South Carolina", "Tennessee", "Texas", "Utah", "Virginia"))
                    # S2 <- S2[!S2 %in% toupper(SATES$NAME)]
                    state_cd <- as.numeric(SATES[toupper(SATES$NAME) %in% toupper(S2),"STATEFP"])
                    if(state %in% state_cd){
                      temporal_list <- as.Date(as.Date(paste0(year,"-05-01")):as.Date(paste0(year+1,"-04-30")),origin="1970-01-01")
                    }

                    # Determine where to aggregate
                    relevant_prism <- unique(prism_list[prism_list$date %in% temporal_list,])
                    relevant_prism <- data.table::rbindlist(
                      lapply(
                        c(1:nrow(relevant_prism)),
                        function(nm){
                          tryCatch({
                            # nm <- 1
                            data <- readRDS(paste0(relevant_prism$file[nm]))
                            data <- data[gridNumber_prism %in% c(spatial_cd[])]
                            return(data)
                          }, error=function(e){NULL})
                        }),fill = T)

                    # Determine where to aggregate for each crop

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

                    relevant_prism$weight <- 1
                    relevant_prism$county_fips <- fip
                    relevant_prism$state_code <- state

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
                        "state_code","county_code","county_fips",
                        "period","dataset_type", "network_type", "resolution"
                      ),
                      names(relevant_prism)
                      ),
                      warming_scenario = seq(0, 5, 0.5)
                    )

                    relevant_prism <- relevant_prism[, c(names(relevant_prism)[!(grepl("ddayN",names(relevant_prism)) | grepl("expN",names(relevant_prism)))]), with = FALSE]

                  }

                  return(relevant_prism)
                }, error = function(e){return(NULL)})
              },year=year,state=state),fill = T)

          res[, commodity_year := year ]

          saveRDS(res, file = file.path(study_environment$wd$prism_climate,year, name))

          rm(res);gc()

        }, error=function(e){})
      }
      return(task)})
}

if(Sys.getenv("SLURM_JOB_NAME") %in% "prism_climate2") {

  optimal_knots <- readRDS("output/optimal_knots.rds")
  # Retain every dday threshold implied by the national knots (004) AND the
  # county-level cluster knots (005, optimal_knots_cluster), so 007 can build
  # county-specific DD. Cluster knots search national +/- knot_band, so this is
  # the union over the band.
  knot_thr <- unique(c(optimal_knots$Tmin, optimal_knots$Tmax))
  if (file.exists("output/optimal_knots_cluster.rds")) {
    clk <- as.data.frame(readRDS("output/optimal_knots_cluster.rds"))
    knot_thr <- unique(c(knot_thr, clk$Tmin, clk$Tmax))
  }
  knot_thr <- sort(knot_thr[is.finite(knot_thr) & knot_thr >= 0])

  dir.create(file.path("data/prism_climate"))

  lapply(
    list.files(study_environment$wd$prism_climate),
    function(year){
      tryCatch({

        df <- data.table::rbindlist(
          lapply(
            list.files(file.path(study_environment$wd$prism_climate,year),recursive = T,full.names = T),
            function(file){
              tryCatch({
                # file <- list.files(study_environment$wd$prism_climate,recursive = T,full.names = T)[1]
                return(readRDS(file))
              }, error=function(e){NULL})
            }), fill = TRUE)

        df[,warming_scenario := as.numeric(as.character(warming_scenario))]

        df <- df[
          ,c("commodity_year","state_code", "county_fips", "period","warming_scenario","dday00",
             paste0("dday", stringr::str_pad(knot_thr, 2, pad = "0"))),
          with = FALSE
        ]

        saveRDS(df,file=file.path("data/prism_climate",paste0("prism_climate_",year,".rds")))

        invisible()

      }, error=function(e){NULL})
    })

}
