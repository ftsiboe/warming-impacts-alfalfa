# 100_WarmingImpactsAlfalfa_exhibits_v00_recovery.R
# ------------------------------------------------------------------------------
# Recovers the v00 alfalfa-availability and cattle-association figures that the
# current article no longer displays, using the CURRENT results in output/summary/.
# Extracted from 100_WarmingImpactsAlfalfa_exhibits.R (Preliminaries + the
# Associations, Alfalfa Availability, and cattle sections). Two adjustments were
# made: the one stale results-folder path was repointed to the current object,
# and every figure is written under the v00_recovery subfolder.
# NOTE: numbers reflect the CURRENT run (2,573 counties, data-driven preferred
# window, 10C lower threshold), NOT v00's original values.
# Run from the repository ROOT, after 006_WarmingImpactsAlfalfa_summary.R.
# ------------------------------------------------------------------------------

#-------------------------------
# Preliminaries              ####
rm(list=ls(all=TRUE))
library(ggplot2);library(terra);library(ggridges);library(gridExtra);library(gtable);library(data.table)
if(Sys.info()['sysname'] =="Windows"){library(gganimate);library(magick)}
if(grepl("windows", sysname)){
  devtools::load_all(file.path(dirname(dirname(getwd())),"packages/gwkit"))
}else{
  devtools::load_all(file.path(dirname(getwd()),"packages/gwkit"))
}
study_environment <- readRDS("data/study_environment.rds")
invisible(lapply(list.files("scripts/helpers", pattern = "[.]R$", full.names = TRUE), source))
myTheme <-   ers_theme() +
  theme(plot.title= element_text(size=10.5),
        axis.title= element_text(size=9,color="black"),
        axis.text = element_text(size=10,color="black"),
        axis.title.y= element_text(size=9,color="black"),
        legend.title=element_blank(),
        legend.text=element_text(size=9),
        plot.caption = element_text(size=8),
        strip.text = element_text(size = 10),
        strip.background = element_rect(fill = "white", colour = "black", size = 1))

USMUR  <- rast(file.path(study_environment$gssurgo_archive,"MURASTER_30m.tif"))
States <- vect(file.path(study_environment$usaPolygons_archive,"USA_States.shp"))
States <- terra::project(States, terra::crs(USMUR))
States <- terra::crop(States, terra::ext(USMUR))

Counties <- vect(file.path(study_environment$usaPolygons_archive,"USA_Counties.shp"))
Counties <- terra::project(Counties, terra::crs(USMUR))
Counties <- terra::crop(Counties, terra::ext(USMUR))
Counties$fip<-as.character(paste0(stringr::str_pad(as.numeric(as.character(Counties$STATEFP)), 2, pad = "0"),
                                  stringr::str_pad(as.numeric(as.character(Counties$COUNTYFP)), 3, pad = "0")))

stnames <- as.data.frame(vect(file.path(study_environment$usaPolygons_archive,"USA_States.shp")))
stnames$state_code <- as.numeric(as.character(stnames$STATEFP))
stnames$State.Name <- stnames$NAME
stnames$State.Abbreviation <- stnames$STUSPS

#-------------------------------
# Preferred window (data-driven) ####
# Highest in-sample R-squared among the candidate windows (5-10 months) that produced
# valid degree-day knots in 003 (i.e. appear in optimal_knots.rds for hay_alfalfa).
# Every "preferred spec" exhibit below filters to `preferred_period` rather than a
# hard-coded window, so the figures track this determination automatically.
.knots_all      <- as.data.frame(readRDS("output/optimal_knots.rds"))
.valid_periods  <- unique(.knots_all$target_periods[.knots_all$crop %in% "hay_alfalfa"])
.r2             <- as.data.frame(readRDS("output/summary/summary_piecewise.rds"))
.r2             <- dplyr::inner_join(.r2, as.data.frame(readRDS("output/optimal_gw.rds"))[1, ])
.r2             <- .r2[.r2$crop %in% "hay_alfalfa" & .r2$climate_base %in% "1991_2020" &
                         .r2$name %in% "r.squared", ]
preferred_period <- select_preferred_period(.r2$period, .r2$Estimate, .valid_periods, 105:110)
preferred_months <- preferred_period - 100
preferred_lab    <- paste0(preferred_months, " months**")
rm(.knots_all, .valid_periods, .r2)

Keep.List<-c("Keep.List",ls())

dir.create("output/exhibits/v00_recovery", recursive = TRUE, showWarnings = FALSE)

# Associations               ####
rm(list= ls()[!(ls() %in% c(Keep.List))])
# 20-spec consensus per county for each association coefficient, via gwkit, in place
# of picking optimal_gw[1,]. (Loaded here too so this section can run on its own.)

assoc <- as.data.frame(readRDS("output/summary/summary_associations.rds"))
assoc <- assoc[!assoc$fip %in% c("0","00000") &
                 assoc$crop %in% "hay_alfalfa" &
                 assoc$period %in% preferred_period &
                 assoc$climate_base %in% "1991_2020", ]
assoc <- assoc[is.finite(assoc$est), ]

# one consensus surface per coefficient (avail00 / prod00 / prod00_LM), then stack
cons <- do.call(rbind, lapply(c("avail00", "prod00", "prod00_LM"), function(nm) {
  cn <- as.data.frame(gw_optimal_scalar_by_polygon(
    value_dt = assoc[assoc$name %in% nm, ], unit_col = "fip",
    polygons = Counties, value_col = "est", agg_fun = stats::median,
    queen_smooth = FALSE))
  cn$name <- nm
  cn
}))

sf_object <- sf::st_as_sf(terra::merge(Counties, cons, by = "fip"))   # one-to-many by name
sf_object <- sf_object[is.finite(sf_object$consensus), ]

# Per-coefficient equal-count quantile classes. The three facets (avail00, prod00,
# prod00_LM) differ by orders of magnitude, so a single absolute scale collapses the
# within-facet variation; binning WITHIN each name makes every facet span all colours.
nbin  <- 6L
qlabs <- c("Lowest","Low","Mid-low","Mid-high","High","Highest")
sf_object$Value <- factor(
  qlabs[ave(sf_object$consensus, sf_object$name, FUN = function(x) dplyr::ntile(x, nbin))],
  levels = qlabs)
sf_object <- sf_object[!is.na(sf_object$Value), ]

# Descriptive panel labels. Availability decomposes as own-county production plus the
# neighbouring (spatially-lagged) production: avail00 = prod00 + prod00_LM.
sf_object$panel <- factor(sf_object$name,
  levels = c("avail00", "prod00", "prod00_LM"),
  labels = c("(a) Total alfalfa availability",
             "(b) Own-county production",
             "(c) Neighbouring production"))

Fig06 <- ggplot() +
  geom_sf(data = sf::st_as_sf(States), colour = "black", fill = "darkred",size = 0.2) +
  geom_sf(data = sf_object,aes(fill = Value), colour = NA,size = 0.2) +
  geom_sf(data = sf::st_as_sf(States), colour = "black", fill = NA,size = 0.2) +
  scale_fill_manual(drop=FALSE, values=colorRampPalette(c("yellow","darkgreen"))(nlevels(sf_object$Value)), na.value="#EEEEEE",
                    name="Within-panel\nquantile\n(low -> high)") +
  labs(title= "", x = "", y = "", caption = "") +
  guides(fill = guide_legend(ncol = 2, override.aes = list(size = 1))) +
  facet_wrap(~panel, ncol = 2) +                       # 2x2: three maps + empty 4th quadrant
  ers_theme() +
  theme_bw() +
  theme(panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        axis.ticks = element_blank(),
        plot.title = element_text(size = 10),
        legend.position = c(0.75, 0.25),               # legend in the empty bottom-right quadrant
        legend.justification = c(0.5, 0.5),
        legend.background = element_blank(),
        legend.key.size = unit(0.35, "cm"),
        legend.text = element_text(size = 9),
        legend.title = element_text(size = 10),
        axis.title = element_blank(),
        axis.text = element_blank(),
        strip.text = element_text(size = 10),
        strip.background = element_blank()) + coord_sf()

ggsave("output/exhibits/v00_recovery/associations.png", Fig06, dpi = 600, width = 9, height = 6)

#-------------------------------
# Alfalfa Availability       ####
rm(list= ls()[!(ls() %in% c(Keep.List))])
# 20-spec consensus of alfalfa availability (baseline, preferred window) via gwkit,
# in place of picking optimal_gw[1,]. The all-specs facet figure is dropped.
avail <- as.data.frame(readRDS("output/summary/summary_availability.rds"))
avail <- avail[avail$crop %in% "hay_alfalfa" & avail$period %in% preferred_period &
                 avail$climate_base %in% "1991_2020" & avail$warming_scenario %in% 0.0, ]
avail <- avail[!avail$fip %in% c("0","00000"), ]
avail$est <- avail$Estimate
avail <- avail[is.finite(avail$est), ]

avail_cons <- as.data.frame(gw_optimal_scalar_by_polygon(
  value_dt = avail, unit_col = "fip", polygons = Counties,
  value_col = "est", agg_fun = stats::median, queen_smooth = FALSE))

sf_object <- sf::st_as_sf(terra::merge(Counties, avail_cons, by = "fip"))
sf_object <- sf_object[is.finite(sf_object$consensus), ]

# equal-count quantile bins (6 classes, ~equal county counts) for max spatial variation
brks <- unique(stats::quantile(sf_object$consensus, probs = seq(0, 1, length.out = 7), na.rm = TRUE))
sf_object$Value <- cut(sf_object$consensus, breaks = brks, include.lowest = TRUE, dig.lab = 4)
sf_object <- sf_object[!is.na(sf_object$Value), ]


Fig_avail <- ggplot() +
  #geom_sf(data = sf::st_as_sf(States), colour = "black", fill = "darkred",size = 0.2) +
  geom_sf(data = sf_object,aes(fill = Value), colour = NA,size = 0.2) +
  geom_sf(data = sf::st_as_sf(States), colour = "black", fill = NA,size = 0.2) +
  scale_fill_manual(drop=FALSE, values=c(colorRampPalette(c("yellow","darkgreen"))(length(unique(as.character(sf_object$Value))))), na.value="#EEEEEE",
                    name="1,000 tons") +
  labs(title= "(a) Alfalfa availability (20-spec consensus)", x = "", y = "",fill ="", fill='',caption = "") +
  guides(fill = guide_legend(ncol=2,override.aes = list(size=1))) +
  ers_theme() +
  theme_bw() +
  theme(panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        axis.ticks.x = element_blank(),
        axis.ticks.y = element_blank(),
        plot.title=element_text(size=10),
        #legend.position="bottom",
        legend.background = element_blank(),
        legend.position=c(0.17,0.18),
        legend.key.size = unit(0.2,"cm"),
        legend.text = element_text(size=7),
        legend.title = element_text(size=7),
        axis.title.y = element_text(size=8),
        axis.title.x = element_text(size=8),
        axis.text.x  = element_blank(), #
        axis.text.y  = element_blank(),
        strip.text = element_text(size = 8),
        strip.background = element_blank())+coord_sf()

# ggsave("output/exhibits/v00_recovery/alfalfa_availability_main.png", Fig_avail, dpi = 600,width = 6.5, height = 5)


#-------------------------------
# cattle                     ####
# "Responsiveness of cattle inventory to alfalfa availability" is the GWR coefficient
# b_avail00, which lives in summary_associations (est where name == "avail00"), NOT in
# summary_cattle (that object holds the cattle *shift* cattleA/B/C, and has no corr_et).
corr_cat <- as.data.frame(readRDS("output/summary/summary_associations.rds"))
corr_cat <- corr_cat[corr_cat$name %in% "avail00" &
                       corr_cat$crop %in% "hay_alfalfa" &
                       corr_cat$period %in% preferred_period &
                       corr_cat$climate_base %in% "1991_2020", ]
corr_cat <- corr_cat[!corr_cat$fip %in% c("0","00000"),]
corr_cat <- corr_cat[is.finite(corr_cat$est), ]

# 20-spec consensus of the responsiveness coefficient via gwkit.
corr_cons <- as.data.frame(gw_optimal_scalar_by_polygon(
  value_dt = corr_cat, unit_col = "fip", polygons = Counties,
  value_col = "est", agg_fun = stats::median, queen_smooth = FALSE))

sf_object <- sf::st_as_sf(terra::merge(Counties, corr_cons, by = "fip"))
sf_object <- sf_object[is.finite(sf_object$consensus), ]

# equal-count quantile bins (6 classes, ~equal county counts) for max spatial variation
brks <- unique(stats::quantile(sf_object$consensus, probs = seq(0, 1, length.out = 7), na.rm = TRUE))
sf_object$Value <- cut(sf_object$consensus, breaks = brks, include.lowest = TRUE, dig.lab = 4)
sf_object <- sf_object[!is.na(sf_object$Value), ]

Fig_corr <- ggplot() +
  #geom_sf(data = sf::st_as_sf(States), colour = "black", fill = "darkred",size = 0.2) +
  geom_sf(data = sf_object,aes(fill = Value), colour = NA,size = 0.2) +
  geom_sf(data = sf::st_as_sf(States), colour = "black", fill = NA,size = 0.2) +
  scale_fill_manual(drop=FALSE, values=c(colorRampPalette(c("#FFEBCD","#800000"))(length(unique(as.character(sf_object$Value))))), na.value="#EEEEEE",
                    name="Percentage") +
  labs(title= "(b) Responsiveness of cattle inventory to alfalfa availability (20-spec consensus)", x = "", y = "",fill ="", fill='',caption = "") +
  guides(fill = guide_legend(ncol=2,override.aes = list(size=1))) +
  ers_theme() +
  theme_bw() +
  theme(panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        axis.ticks.x = element_blank(),
        axis.ticks.y = element_blank(),
        plot.title=element_text(size=10),
        #legend.position="bottom",
        legend.background = element_blank(),
        legend.position=c(0.17,0.18),
        legend.key.size = unit(0.2,"cm"),
        legend.text = element_text(size=7),
        legend.title = element_text(size=7),
        axis.title.y = element_text(size=8),
        axis.title.x = element_text(size=8),
        axis.text.x  = element_blank(), #
        axis.text.y  = element_blank(),
        strip.text = element_text(size = 8),
        strip.background = element_blank())+coord_sf()

# ggsave("output/exhibits/v00_recovery/corr_catle_avail_main.png", Fig_corr, dpi = 600,width = 6.5, height = 5)

marg <- c(-0.01,0.05,-0.01,0.05)

Fig <- cowplot::plot_grid(
  Fig_avail + theme(plot.margin=unit(marg, "cm")) ,
  Fig_corr + theme(plot.margin=unit(marg, "cm")) ,
  ncol=1, align="v",rel_heights=c(1,1),
  greedy=F)

ggsave("output/exhibits/v00_recovery/availability_cattle.png", Fig, dpi = 600,width = 6.5, height = 9)

#-------------------------------
# Predicted impacts (v00 Figures 5 & 6)   ####
# Yield + availability + cattle % impact by warming scenario, per county. Availability
# and cattle are GW outcomes (20 specs) -> per-county 20-spec consensus per scenario via
# gwkit; yield is spec-invariant, so its consensus equals its value. This recovers the
# v00 "Predicted county-level impacts" (Fig 6) and "Predicted impacts" (Fig 5), which the
# original code read from stale Results/ paths (summary_impact_avail / _catle) and which
# expected an Estimate column on cattle (the current object stores cattleA/B/C; the
# availability channel cattleA is used here).
rm(list = ls()[!(ls() %in% c(Keep.List))])
if (requireNamespace("gwkit", quietly = TRUE)) {
  library(gwkit)
} else if (grepl("windows", sysname)) {
  devtools::load_all(file.path(dirname(dirname(getwd())),"packages/gwkit"))
} else {
  devtools::load_all(file.path(dirname(getwd()),"packages/gwkit"))
}

.scen <- c(0.5, 1.0, 1.5, 2.0, 2.5, 3.0)

.impact_consensus <- function(path, value_col, outcome) {
  d <- as.data.frame(readRDS(path))
  d$warming_scenario <- suppressWarnings(as.numeric(as.character(d$warming_scenario)))
  d <- d[d$crop %in% "hay_alfalfa" & d$period %in% preferred_period &
           d$climate_base %in% "1991_2020" & d$warming_scenario %in% .scen &
           !d$fip %in% c("0", "00000"), ]
  d$est <- suppressWarnings(as.numeric(d[[value_col]]))
  d <- d[is.finite(d$est), ]
  do.call(rbind, lapply(.scen, function(sc) {
    ds <- d[d$warming_scenario %in% sc, ]
    if (nrow(ds) == 0) return(NULL)
    cn <- as.data.frame(gw_optimal_scalar_by_polygon(
      ds, unit_col = "fip", polygons = Counties, value_col = "est",
      agg_fun = stats::median, queen_smooth = FALSE))
    data.frame(fip = cn$fip, warming_scenario = sc, est = cn$consensus,
               outcome = outcome, stringsAsFactors = FALSE)
  }))
}

impacts <- rbind(
  .impact_consensus("output/summary/summary_impact_yield.rds", "Estimate", "(a) Alfalfa yield"),
  .impact_consensus("output/summary/summary_availability.rds", "Estimate", "(b) Alfalfa availability"),
  .impact_consensus("output/summary/summary_cattle.rds",       "cattleA",  "(c) Cattle inventory"))
impacts$SimCat <- factor(impacts$warming_scenario, levels = .scen,
                         labels = paste0("+", format(.scen, nsmall = 1), " °C"))
write.csv(impacts, "output/exhibits/figure_data/predicted_impacts_consensus.csv", row.names = FALSE)

# --- Figure 6: predicted county-level impacts (spatial, faceted by scenario) ----------
# Continuous diverging scale centred at 0 (impacts are +/-) for maximum spatial variation.
lim <- stats::quantile(abs(impacts$est), 0.98, na.rm = TRUE)
panel <- function(outcome) {
  sfo <- sf::st_as_sf(terra::merge(Counties, impacts[impacts$outcome %in% outcome, ], by = "fip"))
  sfo <- sfo[is.finite(sfo$est), ]
  ggplot() +
    geom_sf(data = sf::st_as_sf(States), colour = "black", fill = "grey85", size = 0.1) +
    geom_sf(data = sfo, aes(fill = est), colour = NA, size = 0.2) +
    geom_sf(data = sf::st_as_sf(States), colour = "black", fill = NA, size = 0.1) +
    scale_fill_gradient2(low = "#B2182B", mid = "#F7F7F7", high = "#2166AC", midpoint = 0,
                         limits = c(-lim, lim), oob = scales::squish, name = "% impact") +
    labs(title = outcome, x = "", y = "", caption = "") +
    facet_wrap(vars(SimCat), nrow = 1) +
    ers_theme() + theme_bw() +
    theme(panel.grid = element_blank(), axis.ticks = element_blank(),
          plot.title = element_text(size = 9), legend.position = "right",
          legend.key.width = unit(0.3, "cm"), legend.text = element_text(size = 6),
          axis.text = element_blank(), strip.text = element_text(size = 7),
          strip.background = element_blank()) + coord_sf()
}
Fig6 <- cowplot::plot_grid(panel("(a) Alfalfa yield"),
                           panel("(b) Alfalfa availability"),
                           panel("(c) Cattle inventory"),
                           ncol = 1, align = "v")
ggsave("output/exhibits/v00_recovery/predicted_county_impacts.png", Fig6, dpi = 600, width = 9, height = 6)

# --- Figure 5: predicted impacts (v00 three-column dot plot) ---------------------------
# National impact by kernel, distance metric, hay type, season window, climate baseline,
# and warming scenario, faceted by outcome (yield / availability / cattle). Availability
# and cattle come from the current objects (fip == "0" national row; cattle uses the
# availability-channel cattleA). Yield is spec-invariant, so its kernel/distance blocks are
# flat; availability and cattle vary across GW specs.
optimal_gw <- as.data.frame(readRDS("output/optimal_gw.rds"))
.common <- c("p","theta","longlat","DistName","kernel","crop","period","climate_base",
             "warming_scenario","est","se","outcome")

yield <- as.data.frame(readRDS("output/summary/summary_impact_yield.rds"))
yield$est <- yield$Estimate; yield$se <- yield$Estimate_sd
yield <- yield[yield$county_code %in% 0 & yield$state_code %in% 0 & yield$region %in% "", ]
yield$outcome <- "(a) Alfalfa yield"
yield <- yield[, .common]

# summary_availability / summary_cattle are per-county (GW) with no national aggregate row,
# so build the national series by averaging the per-county impact within each spec/scenario
# (se = across-county sd). Yield keeps its stored national bootstrap se.
.national <- function(path, val, outcome_label) {
  d <- data.table::as.data.table(readRDS(path))
  d <- d[!fip %in% c("0", "00000")]
  d[, est := suppressWarnings(as.numeric(get(val)))]
  d[, warming_scenario := suppressWarnings(as.numeric(as.character(warming_scenario)))]
  d <- d[is.finite(est)]
  keys <- intersect(c("p","theta","longlat","DistName","kernel","crop","period",
                      "climate_base","warming_scenario"), names(d))
  ag <- d[, .(est = mean(est, na.rm = TRUE), se = stats::sd(est, na.rm = TRUE)), by = keys]
  ag[, outcome := outcome_label]
  as.data.frame(ag)[, .common]
}
avail  <- .national("output/summary/summary_availability.rds", "Estimate", "(b) Alfalfa availability")
cattle <- .national("output/summary/summary_cattle.rds",       "cattleA",  "(c) Cattle inventory")

data <- rbind(yield, avail, cattle)
data$warming_scenario <- suppressWarnings(as.numeric(as.character(data$warming_scenario)))

data_main <- data[(data$crop %in% "hay_alfalfa" & data$period %in% preferred_period & data$climate_base %in% "1991_2020"),]
data_main <- dplyr::inner_join(data_main,optimal_gw[1,c("p", "theta", "longlat", "DistName", "kernel")])

data_crop <- data[(data$warming_scenario %in% 1.0 & data$crop %in% c("hay_other","hay_alfalfa") & data$period %in% preferred_period & data$climate_base %in% "1991_2020"),]
data_crop <- dplyr::inner_join(data_crop,optimal_gw[1,c("p", "theta", "longlat", "DistName", "kernel")])

data_year <- data[(data$warming_scenario %in% 1.0 &data$crop %in% "hay_alfalfa" & data$period %in% preferred_period & data$climate_base %in% c("1991_2020","1981_2010","1971_2000","1961_1990")),]
data_year <- dplyr::inner_join(data_year,optimal_gw[1,c("p", "theta", "longlat", "DistName", "kernel")])

data_wind <- data[(data$warming_scenario %in% 1.0 &data$crop %in% "hay_alfalfa" & data$period %in% c(105:110) & data$climate_base %in% "1991_2020"),]
data_wind <- dplyr::inner_join(data_wind,optimal_gw[1,c("p", "theta", "longlat", "DistName", "kernel")])

data_dist <- data[(data$warming_scenario %in% 1.0 & data$crop %in% "hay_alfalfa" & data$period %in% preferred_period &
                     data$climate_base %in% "1991_2020" & data$kernel %in% optimal_gw[1,"kernel"]),]
data_kenl <- data[(data$warming_scenario %in% 1.0 & data$crop %in% "hay_alfalfa" & data$period %in% preferred_period &
                     data$climate_base %in% "1991_2020" & data$DistName %in% optimal_gw[1,"DistName"]),]

data_main$type <- "Impact by warming scenario"
data_year$type <- "1+°C warming by climate baseline"
data_wind$type <- "1+°C warming by season window"
data_crop$type <- "1+°C warming by hay type"
data_kenl$type <- "1+°C warming by kernel"
data_dist$type <- "1+°C warming by distance metric"

data_wind$period <- ifelse(data_wind$period %in% 0,113,data_wind$period)
data_wind <- data_wind[order(-data_wind$period),]

data_wind$x_name <- paste0(data_wind$period-100," months")
data_wind$x_name <- ifelse(data_wind$x_name %in% "13 months","All months",data_wind$x_name)
data_year$x_name <- data_year$climate_base
data_crop$x_name <- as.character(factor(data_crop$crop,levels = c("hay_other","hay_all","hay_alfalfa"),
                           labels = c("Non-alfalfa","Combined","Alfalfa")))
data_main$x_name <- paste0(format(as.numeric(as.character(data_main$warming_scenario)), nsmall = 1) ,"+°C")
data_kenl$x_name <- as.character(factor(data_kenl$kernel,levels = c("exponential","gaussian","boxcar","bisquare","tricube"),
                                        labels = c("Exponential","Gaussian","Boxcar","Bisquare","Tricube")))
data_dist$x_name <- as.character(factor(data_dist$DistName,levels = c("Manhattan distance metric","Euclidean distance metric",
                                                                      "Coordinate system is rotated by an angle 0.8 in radian"),
                                        labels = c("Manhattan","Euclidean","Coordinate system")))

# x positions taken from the yield column, then shared across all three facets.
data_kenly <- data_kenl[data_kenl$outcome %in% "(a) Alfalfa yield",]
data_disty <- data_dist[data_dist$outcome %in% "(a) Alfalfa yield",]
data_cropy <- data_crop[data_crop$outcome %in% "(a) Alfalfa yield",]
data_windy <- data_wind[data_wind$outcome %in% "(a) Alfalfa yield",]
data_yeary <- data_year[data_year$outcome %in% "(a) Alfalfa yield",]
data_mainy <- data_main[data_main$outcome %in% "(a) Alfalfa yield",]

data_kenly$x <- 1:nrow(data_kenly) + 2
data_disty$x <- max(data_kenly$x) + c(1:nrow(data_disty)) + 2
data_cropy$x <- max(data_disty$x) + c(1:nrow(data_cropy)) + 2
data_windy$x <- max(data_cropy$x) + c(1:nrow(data_windy)) + 2
data_yeary$x <- max(data_windy$x) + c(1:nrow(data_yeary)) + 2
data_mainy$x <- max(data_yeary$x) + c(1:nrow(data_mainy)) + 2

datay <- rbind(data_cropy,data_mainy,data_yeary,data_windy,data_kenly,data_disty)
datay <- unique(datay[c("p","theta","longlat","DistName","kernel","crop","period","climate_base","warming_scenario","x_name","x")])
data <- rbind(data_crop,data_main,data_year,data_wind,data_kenl,data_dist)
data <- dplyr::inner_join(data,datay)

data_text <- doBy::summaryBy(x~type,data=data[data$outcome %in% "(a) Alfalfa yield",],FUN=max,keep.names = T,na.rm=T)
data_text$outcome <- "(a) Alfalfa yield"
.ylab <- min(data$est[data$outcome %in% "(a) Alfalfa yield"], na.rm = TRUE)

Fig04 <- ggplot(data,aes(x=x,y=est,group=1)) +
  geom_hline(yintercept = 0,size = 0.2,color = "black") +
  geom_errorbar(aes(ymin = est - se*1.96, ymax = est + se*1.96),color="purple") +
  geom_point(color="purple") +
  geom_text(data=data_text,aes(x=x+1.5,y=.ylab,label=type, hjust = 0), size = 2.5, col = "black",
            stat = "identity",check_overlap = TRUE) +
  facet_wrap(~outcome, nrow = 1, scales = "free") +
  scale_x_continuous(breaks=datay$x,labels=datay$x_name) +
  labs(title="", x="", y ="\nImpact (%)",caption = "") +
  myTheme +
  theme(axis.text.x = element_text(size=7,color="black"),
        axis.text.y = element_text(size=7,color="black"),
        strip.text = element_text(size = 9),
        legend.position="none")+ coord_flip()

write.csv(data, "output/exhibits/figure_data/predicted_impacts_mean.csv", row.names = FALSE)
ggsave("output/exhibits/v00_recovery/predicted_impacts_mean.png", Fig04, dpi = 600, width = 10, height = 7.5)
