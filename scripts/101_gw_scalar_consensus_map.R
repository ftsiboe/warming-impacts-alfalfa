# 101_gw_scalar_consensus_map.R
# ------------------------------------------------------------------------------
# One demonstration map built with the new gwkit scalar-consensus function.
#
# The article currently reduces its 20 geographically-weighted specifications
# (4 distance metrics x 5 kernels) to results by simply picking the first,
# optimal_gw[1,]. Here we instead form a per-county CONSENSUS across all 20
# specs with gwkit::gw_consensus_scalar(), for one continuous GW
# outcome: the responsiveness of cattle inventory to alfalfa availability
# (the "avail00" coefficient in summary_associations).
#
# Output: output/exhibits/v00_recovery/consensus_avail00_responsiveness.png
#         output/exhibits/figure_data/consensus_avail00_responsiveness.csv
# Run from the repository ROOT, after 006_WarmingImpactsAlfalfa_summary.R.
# ------------------------------------------------------------------------------

rm(list = ls(all = TRUE))
library(ggplot2); library(terra); library(sf); library(data.table)

study_environment <- readRDS("data/study_environment.rds")
# project helpers (ers_theme, select_preferred_period, ...)
invisible(lapply(list.files("scripts/helpers", pattern = "[.]R$", full.names = TRUE), source))

# --- load gwkit (sibling package under ../packages/gwkit) ---------------------
if(grepl("windows", sysname)){
  devtools::load_all(file.path(dirname(dirname(getwd())),"packages/gwkit"))
}else{
  devtools::load_all(file.path(dirname(getwd()),"packages/gwkit"))
}

dir.create("output/exhibits/v00_recovery", recursive = TRUE, showWarnings = FALSE)
dir.create("output/exhibits/figure_data",  recursive = TRUE, showWarnings = FALSE)

# --- data-driven preferred window (same rule as the article) -----------------
# The multi-spec (20-way) associations are computed only for the preferred window;
# other windows carry a single spec, so the consensus must be taken on this one.
.knots_all     <- as.data.frame(readRDS("output/optimal_knots.rds"))
.valid_periods <- unique(.knots_all$target_periods[.knots_all$crop %in% "hay_alfalfa"])
.r2            <- as.data.frame(readRDS("output/summary/summary_piecewise.rds"))
.r2            <- dplyr::inner_join(.r2, as.data.frame(readRDS("output/optimal_gw.rds"))[1, ])
.r2            <- .r2[.r2$crop %in% "hay_alfalfa" & .r2$climate_base %in% "1991_2020" &
                        .r2$name %in% "r.squared", ]
preferred_period <- select_preferred_period(.r2$period, .r2$Estimate, .valid_periods, 105:110)

# --- county geometry ----------------------------------------------------------
USMUR    <- rast(file.path(study_environment$gssurgo_archive, "MURASTER_30m.tif"))
States   <- vect(file.path(study_environment$usaPolygons_archive, "USA_States.shp"))
States   <- terra::crop(terra::project(States, terra::crs(USMUR)), terra::ext(USMUR))
Counties <- vect(file.path(study_environment$usaPolygons_archive, "USA_Counties.shp"))
Counties <- terra::crop(terra::project(Counties, terra::crs(USMUR)), terra::ext(USMUR))
Counties$fip <- as.character(paste0(
  stringr::str_pad(as.numeric(as.character(Counties$STATEFP)),  2, pad = "0"),
  stringr::str_pad(as.numeric(as.character(Counties$COUNTYFP)), 3, pad = "0")))

# --- stacked GW outcome across the 20 specs -----------------------------------
# One row per county per (kernel x distance_metric) spec; value = est of avail00.
assoc <- as.data.frame(readRDS("output/summary/summary_associations.rds"))
assoc <- assoc[assoc$name %in% "avail00" &
                 assoc$crop %in% "hay_alfalfa" &
                 assoc$period %in% preferred_period &
                 assoc$climate_base %in% "1991_2020", ]
assoc <- assoc[!assoc$fip %in% c("0", "00000"), ]
assoc <- assoc[is.finite(assoc$est), ]

stopifnot(nrow(assoc) > 0)
message("specs per county (summary): ",
        paste(range(table(assoc$fip)), collapse = "-"),
        "  (expected up to 20)")

# --- consensus across specs, per county (median; equal weights) ---------------
consensus <- gw_consensus_scalar(
  value_dt     = assoc,
  unit_col     = "fip",
  geometry     = Counties,
  value_col    = "est",
  agg_fun      = stats::median,   # user-suppliable reducer
  probs        = c(0.05, 0.95),
  queen_smooth = FALSE            # pure across-spec consensus for this map
)
consensus <- as.data.frame(consensus)
write.csv(consensus,
          "output/exhibits/figure_data/consensus_avail00_responsiveness.csv",
          row.names = FALSE)

# --- map the consensus --------------------------------------------------------
merged   <- terra::merge(Counties, consensus, by = "fip")
sf_object <- sf::st_as_sf(merged)
sf_object <- sf_object[is.finite(sf_object$consensus), ]

lim <- stats::quantile(abs(sf_object$consensus), 0.98, na.rm = TRUE)  # symmetric clip

Fig <- ggplot() +
  geom_sf(data = sf::st_as_sf(States), colour = "black", fill = "grey95", size = 0.2) +
  geom_sf(data = sf_object, aes(fill = consensus), colour = NA, size = 0.2) +
  geom_sf(data = sf::st_as_sf(States), colour = "black", fill = NA, size = 0.2) +
  scale_fill_gradient2(low = "#B2182B", mid = "#F7F7F7", high = "#2166AC",
                       midpoint = 0, limits = c(-lim, lim), oob = scales::squish,
                       name = "Consensus\ncoefficient") +
  labs(title = "Responsiveness of cattle inventory to alfalfa availability",
       subtitle = paste0("Per-county median consensus across all 20 GW specs ",
                         "(4 distance metrics x 5 kernels)"),
       x = "", y = "", caption = "gwkit::gw_consensus_scalar()") +
  ers_theme() + theme_bw() +
  theme(panel.grid = element_blank(),
        axis.text = element_blank(), axis.ticks = element_blank(),
        legend.position = "right",
        plot.title = element_text(size = 11),
        plot.subtitle = element_text(size = 9)) +
  coord_sf()

ggsave("output/exhibits/v00_recovery/consensus_avail00_responsiveness.png",
       Fig, dpi = 600, width = 7.5, height = 5)

message("Wrote consensus map for ", nrow(sf_object), " counties; ",
        "median sign-agreement across specs = ",
        round(stats::median(consensus$sign_agreement, na.rm = TRUE), 3))
