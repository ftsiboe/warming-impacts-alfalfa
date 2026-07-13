#' Build county-level hay yield and weather panel data
#'
#' Constructs a county-level panel dataset by combining PRISM weather data with
#' NASS hay production data for a selected crop. The function filters weather
#' records to the baseline warming scenario (`warming_scenario == 0`), aggregates
#' monthly or seasonal weather measures to the requested target period, merges
#' the results with county-level hay yield data, and creates additional derived
#' variables commonly used in regression analysis.
#'
#' @param crop Character scalar. Commodity name to extract from the hay
#'   production data. Default is `"hay_alfalfa"`.
#' @param target_periods Integer. Defines the weather periods to aggregate.
#'   If `target_periods` is in `1:12`, weather is aggregated from period 1
#'   up to and including `target_periods`. Otherwise, weather is filtered to
#'   the exact value(s) supplied in `target_periods`.
#' @param prism_weather_directory directory to prism weather data files
#'
#' @details
#' The function performs the following steps:
#' \itemize{
#'   \item Reads PRISM weather data from `"data/prism_weather.rds"` and keeps
#'   only records where `warming_scenario == 0`.
#'   \item Selects weather periods based on `target_periods`.
#'   \item Replaces negative values with zero for precipitation, freeze,
#'   growing degree day variables (`dday*`), and extreme weather variables
#'   (`exp*`).
#'   \item Aggregates selected weather variables to the
#'   `commodity_year`-`state_code`-`county_code` level.
#'   \item Reads hay production data from `"data/nass_hay_production.rds"` and
#'   keeps the requested crop.
#'   \item Merges the weather and yield data.
#'   \item Creates a county identifier `fip` using state, ASD, and county codes.
#'   \item Restricts the sample to counties with at least two yearly
#'   observations.
#'   \item Creates derived variables:
#'   \describe{
#'     \item{lny}{Log yield}
#'     \item{ppt}{Precipitation converted from mm to inches}
#'     \item{ppt2}{Squared precipitation}
#'     \item{ppt3}{Cubed precipitation}
#'     \item{Trend}{Linear time trend starting at 1}
#'     \item{Trend2}{Squared time trend}
#'   }
#' }
#'
#' @return A `data.frame` or `data.table`-like object containing:
#' \itemize{
#'   \item identifiers: `commodity_year`, `state_code`, `asd_code`,
#'   `county_code`, `fip`
#'   \item outcome: `yield`, `lny`
#'   \item aggregated weather variables
#'   \item derived controls such as `ppt`, `ppt2`, `ppt3`, `Trend`, `Trend2`
#' }
#'
#' @examples
#' \dontrun{
#'   df <- build_hay_weather_panel()
#'   df <- build_hay_weather_panel(crop = "hay_other", target_periods = 6)
#'   df <- build_hay_weather_panel(crop = "hay_alfalfa", target_periods = 107)
#' }
#'
#' @export
build_hay_weather_panel <- function(crop = "hay_alfalfa", target_periods = 107, prism_weather_directory){

  weather <- data.table::rbindlist(
    lapply(
      list.files(prism_weather_directory,full.names = T),
      function(file){
        tryCatch({
          return(readRDS(file))
        }, error=function(e){NULL})
      }), fill = TRUE);gc()

  if (!target_periods %in% 1:12) {
    weather <- weather[period %in% target_periods]
  } else {
    weather <- weather[period %in% 1:12]
    weather <- weather[period <= target_periods, ]
  }

  weather_vars <- c(
    "precipitation",
    "freeze",
    names(weather)[grepl("dday", names(weather))],
    names(weather)[grepl("exp", names(weather))]
  )

  for (xx in weather_vars) {
    weather[, (xx) := ifelse(get(xx) < 0, 0, get(xx))]
  }

  weather <- weather[
    ,
    lapply(.SD, function(x) sum(x, na.rm = TRUE)),
    by = c("commodity_year", "state_code", "county_code"),
    .SDcols = weather_vars
  ]

  data <- readRDS("data/nass_hay_production.rds")
  data <- data[
    commodity_name %in% crop,
    c("commodity_year", "state_code", "asd_code", "county_code", "yield"),
    with = FALSE
  ]

  data <- dplyr::inner_join(
    data,
    weather,
    by = c("commodity_year", "state_code", "county_code")
  )

  data[, fip := paste0(
    stringr::str_pad(state_code, pad = "0", 2),
    stringr::str_pad(asd_code,   pad = "0", 2),
    stringr::str_pad(county_code, pad = "0", 3)
  )]

  data[, county_fips := paste0(
    stringr::str_pad(state_code, pad = "0", 2),
    stringr::str_pad(county_code, pad = "0", 3)
  )]
  
  panel <- doBy::summaryBy(commodity_year ~ fip, data = data, FUN = length)
  panel <- panel[panel$commodity_year.length >= 2, ]
  data <- data[fip %in% panel$fip]

  data$lny <- log(data$yield)
  data$ppt <- data$precipitation / 25.4
  data$ppt2 <- data$ppt^2
  data$ppt3 <- data$ppt^3

  data$Trend <- data$commodity_year - min(data$commodity_year) + 1
  data$Trend2 <- data$Trend * data$Trend

  return(data)
}
