#-------------------------------
# Preliminaries              ####
rm(list=ls(all=TRUE))
library(ggplot2);library(terra);library(ggridges);library(gridExtra);library(gtable);library(data.table)
if(Sys.info()['sysname'] =="Windows"){library(gganimate);library(magick)}
study_environment <- readRDS("data/study_environment.rds")
invisible(lapply(list.files("scripts/helpers", pattern = "[.]R$", full.names = TRUE), source))
sysname <- tolower(as.character(Sys.info()[["sysname"]]))
if(grepl("windows", sysname)){
  devtools::load_all(file.path(dirname(dirname(getwd())),"packages/gwkit"))
}else{
  devtools::load_all(file.path(dirname(getwd()),"packages/gwkit"))
}
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
.r2 <- .r2
.r2             <- .r2[.r2$crop %in% "hay_alfalfa" & .r2$climate_base %in% "1991_2020" &
                         .r2$name %in% "r.squared", ]
preferred_period <- select_preferred_period(.r2$period, .r2$Estimate, .valid_periods, 105:110)
preferred_months <- preferred_period - 100
preferred_lab    <- paste0(preferred_months, " months**")
rm(.knots_all, .valid_periods, .r2)

Keep.List<-c("Keep.List",ls())
#-------------------------------
# Spatial Rep                ####
rm(list= ls()[!(ls() %in% c(Keep.List))])

Plot.data <- readRDS("data/spatial_representation.rds")
spa_rep <- Plot.data
spa_rep <- dplyr::inner_join(spa_rep,stnames[c("state_code","State.Name","STUSPS" )],by="state_code")

spa_rep$region <- ifelse(
  spa_rep$STUSPS %in% c("CT", "DE", "ME", "MD", "MA", "NH", "NJ", "NY", "PA", "RI", "VT"),"Northeast",
  NA)
spa_rep$region <- ifelse(
  spa_rep$STUSPS %in% c("IL", "IN", "IA", "KS", "MI", "MN", "MO", "NE", "ND", "OH", "SD", "WI"),"Midwest",
  spa_rep$region)
spa_rep$region <- ifelse(
  spa_rep$STUSPS %in% c("AL", "AR", "FL", "GA", "KY", "LA", "MS", "NC", "OK", "SC", "TN", "TX", "VA", "WV"),"South",
  spa_rep$region)
spa_rep$region <- ifelse(
  spa_rep$STUSPS %in% c("AZ", "CA", "CO", "HI", "ID", "MT", "NV", "NM", "OR", "UT", "WA", "WY"),"West",
  spa_rep$region)


rnk_reg <- doBy::summaryBy(area+inventory+production+yield~region,FUN=mean,data=spa_rep,na.rm=T,keep.names = T)
rnk_reg <- rnk_reg[complete.cases(rnk_reg),]
tops <- list()
for(xx in c("inventory","production","area","yield")){
  rnk_reg <- rnk_reg[order(-rnk_reg[,xx]),]
  tops[[xx]] <- paste0(paste0(rnk_reg$region," (",round(rnk_reg[,xx],2),")"),collapse = ", ")
}
tops

cor(spa_rep[c("inventory","production","area","yield")],use = "pairwise.complete.obs")

datax <- spa_rep[spa_rep$region %in% "South",]

cor(datax[c("inventory","production","area","yield")],use = "pairwise.complete.obs")
state_data <- doBy::summaryBy(area+inventory+production+yield~State.Name,FUN=mean,data=datax,na.rm=T,keep.names = T)
tops <- list()
for(xx in c("inventory","yield","area","production")){
  state_data <- state_data[order(-state_data[,xx]),]
  tops[[xx]]   <- paste0(paste0(state_data$State.Name[1:5]," (",round(state_data[,xx][1:5],2),")"),collapse = ", ")
}
tops
state_data <- state_data[order(state_data[,xx]),]
paste0(paste0(state_data$State.Name[1:5]," (",
              round(state_data[,"area"][1:5],2),", ",
              round(state_data[,"production"][1:5],2),", ",
              round(state_data[,"yield"][1:5],2),")"),collapse = ", ")

state_data <- doBy::summaryBy(area+inventory+production+yield~state_code,FUN=c(mean,sd),data=Plot.data,na.rm=T)
state_data <- dplyr::inner_join(state_data,stnames[c("state_code","State.Name","STUSPS" )],by="state_code")
write.csv(state_data,"output/exhibits/figure_data/spatial_Rep.csv")

Plot.data <- Plot.data |>  tidyr::gather(Type, Value, c("area","production","yield","inventory"))

USMUR  <- rast(file.path(study_environment$gssurgo_archive,"MURASTER_30m.tif"))
States <- vect(file.path(study_environment$usaPolygons_archive,"USA_States.shp"))
States <- terra::project(States, terra::crs(USMUR))
States <- terra::crop(States, terra::ext(USMUR))

Counties <- vect(file.path(study_environment$usaPolygons_archive,"USA_Counties.shp"))
Counties <- terra::project(Counties, terra::crs(USMUR))
Counties <- terra::crop(Counties, terra::ext(USMUR))
Counties$fip<-as.character(paste0(stringr::str_pad(as.numeric(as.character(Counties$STATEFP)), 2, pad = "0"),
                                  stringr::str_pad(as.numeric(as.character(Counties$COUNTYFP)), 3, pad = "0")))

plotlist <- list("area"=list("(a) Alfalfa Acres Harvested",c("cornsilk","saddlebrown"),"1,000 acres"),
                 "production"=list("(b) Alfalfa Production Output",c("yellow","darkgreen"),"1,000 tons"),
                 "yield"=list("(c) Alfalfa Yield",c("lavender","midnightblue"),"1,000 tons/acre"),
                 "inventory"=list("(d) Cattle inventory (including calves) as of first of january",c("#D8BFD8","#4B0082"),"1,000 head"))

fig_fxn <- function(varp){
  # varp <- "inventory"
  merged_spat_vector <- terra::merge(Counties, Plot.data[Plot.data$Type %in% varp,], by="fip")
  sf_object <- sf::st_as_sf(merged_spat_vector)
  sf_object$Value <- cut(sf_object$Value,c(0,max(sf_object$Value,na.rm=T),quantile(sf_object$Value[! sf_object$Value %in% c(NA,Inf,-Inf,NaN)], probs = seq(0.1,1,0.2))))
  sf_object <- sf_object[!is.na(sf_object$Value), ]
  ggplot() +
    geom_sf(data = sf_object,aes(fill = Value), colour = NA,size = 0.2) +
    geom_sf(data = sf::st_as_sf(States), colour = "black", fill = NA,size = 0.2) +
    scale_fill_manual(drop=FALSE, values=colorRampPalette(plotlist[[varp]][[2]])(length(unique(as.character(sf_object$Value)))), na.value="#EEEEEE", name=plotlist[[varp]][[3]]) +
    labs(title= plotlist[[varp]][[1]], x = "", y = "",fill ="", fill='',caption = "") +
    guides(fill = guide_legend(nrow=6,override.aes = list(size=1))) +
    ers_theme() +
    theme_bw() +
    theme(panel.grid.major = element_blank(),
          panel.grid.minor = element_blank(),
          axis.ticks.x = element_blank(),
          axis.ticks.y = element_blank(),
          plot.title=element_text(size=10),
          #legend.position="right",
          legend.background = element_blank(),
          legend.position=c(0.17,0.18),
          legend.key.size = unit(0.2,"cm"),
          legend.text=element_text(size=7),
          legend.title=element_text(size=7),
          axis.title.y = element_blank(),
          axis.text.x  = element_blank(), #
          axis.text.y  = element_blank())+coord_sf()
}

marg <- c(0.05,0.05,-0.01,0.05)

Fig1 <- cowplot::plot_grid(
  fig_fxn("area") + theme(plot.margin=unit(marg, "cm")) ,
  fig_fxn("production") + theme(plot.margin=unit(marg, "cm")) ,
  fig_fxn("yield") + theme(plot.margin=unit(marg, "cm")) ,
  fig_fxn("inventory")+ theme(plot.margin=unit(marg, "cm")) ,
  ncol=2, align="v",rel_heights=c(1,1,1,1),
  greedy=F)

ggsave("output/exhibits/spatial_Rep.png", Fig1, dpi = 600,width = 10, height = 7)

#-------------------------------
# Regression variables       ####
rm(list= ls()[!(ls() %in% c(Keep.List))])
knots <- readRDS("output/optimal_knots.rds")
knots <- knots[(knots$target_periods %in% c(104:112,0) & knots$crop %in% "hay_alfalfa"),]
knots <- knots[order(knots$cv_error),]
Tmin <- knots$Tmin[1]
Tmax <- knots$Tmax[1]

data <- build_hay_weather_panel(
  crop = "hay_alfalfa", target_periods = preferred_period,   # data-driven preferred window (Table 1 & Figure 1)
  prism_weather_directory = "data/prism_weather")
data <- as.data.frame(data)
data$DD1 <- data[,paste0("dday00")] - data[,paste0("dday",stringr::str_pad(Tmin,pad="0",2))]
data$DD2 <- data[,paste0("dday",stringr::str_pad(Tmin,pad="0",2))] - data[,paste0("dday",stringr::str_pad(Tmax,pad="0",2))]
data$DD3 <- data[,paste0("dday",stringr::str_pad(Tmax,pad="0",2))]
data$DD1 <- ifelse(data$DD1<0,0,data$DD1)
data$DD2 <- ifelse(data$DD2<0,0,data$DD2)
data$DD3 <- ifelse(data$DD3<0,0,data$DD3)

# Degree Days 3
DD1 <- ggplot(data=data, aes(x=commodity_year, y=DD1,group = commodity_year)) +
  geom_boxplot(color="blue",outlier.size=0.5) +
  labs(title=paste0("(b) Degree-day (below ",Tmin,"°C)"),
       x="\nYear",
       y ="inches\n",
       caption = "") +
  ers_theme() +
  theme_bw() +
  theme(panel.grid.major = element_blank(),panel.grid.minor = element_blank(), axis.ticks.x = element_blank()) +
  theme(legend.position="none")+
  theme(legend.text= element_text(size=11),
        legend.title= element_blank(),
        axis.title.y= element_text(size=11),
        axis.title.x= element_text(size = 8, colour="black"),
        axis.text.x = element_text(size = 6,colour="black"), #
        axis.text.y = element_text(size = 8, colour="black"),
        plot.caption = element_text(size=11,hjust = 0 ,vjust = 0, face = "italic"),
        strip.text = element_text(size = 12),
        strip.background = element_rect(fill = "white", colour = "black", size = 1))
ggsave("output/exhibits/DD1.png", DD1, dpi = 600,width = 9, height = 6)

# Degree Days 2
DD2 <- ggplot(data=data, aes(x=commodity_year, y=DD2,group = commodity_year)) +
  geom_boxplot(color="blue",outlier.size=0.5) +
  labs(title=paste0("(c) Degree-day (",Tmin,"°C to ",Tmax,"°C)"),
       x="\nYear",
       y ="inches\n",
       caption = "") +
  ers_theme() +
  theme_bw() +
  theme(panel.grid.major = element_blank(),panel.grid.minor = element_blank(), axis.ticks.x = element_blank()) +
  theme(legend.position="none")+
  theme(legend.text= element_text(size=11),
        legend.title= element_blank(),
        axis.title.y= element_text(size=11),
        axis.title.x= element_text(size = 8, colour="black"),
        axis.text.x = element_text(size = 6,colour="black"), #
        axis.text.y = element_text(size = 8, colour="black"),
        plot.caption = element_text(size=11,hjust = 0 ,vjust = 0, face = "italic"),
        strip.text = element_text(size = 12),
        strip.background = element_rect(fill = "white", colour = "black", size = 1))
ggsave("output/exhibits/DD2.png", DD2, dpi = 600,width = 9, height = 6)

# Degree Days 3
DD3 <- ggplot(data=data, aes(x=commodity_year, y=DD3,group = commodity_year)) +
  geom_boxplot(color="blue",outlier.size=0.5) +
  labs(title=paste0("(d) Degree-day (above ",Tmax,"°C)"),
       x="\nYear",
       y ="inches\n",
       caption = "") +
  ers_theme() +
  theme_bw() +
  theme(panel.grid.major = element_blank(),panel.grid.minor = element_blank(), axis.ticks.x = element_blank()) +
  theme(legend.position="none")+
  theme(legend.text= element_text(size=11),
        legend.title= element_blank(),
        axis.title.y= element_text(size=11),
        axis.title.x= element_text(size = 8, colour="black"),
        axis.text.x = element_text(size = 6,colour="black"), #
        axis.text.y = element_text(size = 8, colour="black"),
        plot.caption = element_text(size=11,hjust = 0 ,vjust = 0, face = "italic"),
        strip.text = element_text(size = 12),
        strip.background = element_rect(fill = "white", colour = "black", size = 1))
ggsave("output/exhibits/DD3.png", DD3, dpi = 600,width = 9, height = 6)

# Precipitation
ppt <- ggplot(data=data, aes(x=commodity_year, y=ppt,group = commodity_year)) +
  geom_boxplot(color="blue",outlier.size=0.5) +
  labs(title="(e) Precipitation (inches)",
       x="\nYear",
       y ="inches\n",
       caption = "") +
  ers_theme() +
  theme_bw() +
  theme(panel.grid.major = element_blank(),panel.grid.minor = element_blank(), axis.ticks.x = element_blank()) +
  theme(legend.position="none")+
  theme(legend.text= element_text(size=11),
        legend.title= element_blank(),
        axis.title.y= element_text(size=11),
        axis.title.x= element_text(size = 8, colour="black"),
        axis.text.x = element_text(size = 6,colour="black"), #
        axis.text.y = element_text(size = 8, colour="black"),
        plot.caption = element_text(size=11,hjust = 0 ,vjust = 0, face = "italic"),
        strip.text = element_text(size = 12),
        strip.background = element_rect(fill = "white", colour = "black", size = 1))
ggsave("output/exhibits/ppt.png", ppt, dpi = 600,width = 9, height = 6)


# Yield
yield <- ggplot(data=data[data$yield<15,], aes(x=commodity_year, y=yield,group = commodity_year)) +
  geom_boxplot(color="blue",outlier.size=0.5) +
  labs(title="(a) Yield (tons/acre)",
       x="\nYear", y ="inches\n",
       caption = "") +
  ers_theme() +
  theme_bw() +
  theme(panel.grid.major = element_blank(),panel.grid.minor = element_blank(), axis.ticks.x = element_blank()) +
  theme(legend.position="none")+
  theme(legend.text= element_text(size=11),
        legend.title= element_blank(),
        axis.title.y= element_text(size=11),
        axis.title.x= element_text(size = 8, colour="black"),
        axis.text.x = element_text(size = 6,colour="black"), #
        axis.text.y = element_text(size = 8, colour="black"),
        plot.caption = element_text(size=11,hjust = 0 ,vjust = 0, face = "italic"),
        strip.text = element_text(size = 12),
        strip.background = element_rect(fill = "white", colour = "black", size = 1))
ggsave("output/exhibits/yield.png", yield, dpi = 600,width = 9, height = 6)

fig <- cowplot::plot_grid(DD1 , DD2 , DD3,ppt,ncol=2, align="v",greedy=F)
fig <- cowplot::plot_grid(yield,fig,ncol=1, align="h",rel_heights=c(0.5,1),greedy=F)
ggsave("output/exhibits/regression_variables.png", fig, dpi = 600,width = 8, height = 11)


nass_county <- as.data.frame(readRDS("data/nass_hay_production.rds"))

table1 <- dplyr::full_join(
  nass_county[nass_county$commodity_name %in% "hay_alfalfa",
              c("commodity_year","state_code","asd_code","county_code","area","production")],
  data)

Cattle <- as.data.frame(readRDS("data/nass_animal_inventory.rds"))
table1 <- dplyr::full_join(table1,Cattle[c("commodity_year","state_code","county_code","cattle")])
table1 <- table1[!table1$yield %in% NA,]

table1$area <- table1$area/1000
table1$cattle <- table1$cattle/1000
table1$production <- table1$production/1000
table1$x<-1
table1 <- doBy::summaryBy(area+cattle+production+yield+DD1+DD2+DD3+ppt~x,FUN=c(mean,sd),data=table1,na.rm=T)
table1 <- table1 |> tidyr::gather(variable, value, 2:ncol(table1))
table1 <- tidyr::separate(table1,"variable",into = c("variable","stat"),sep="[.]")
table1 <- table1 |> tidyr::spread(stat, value)
table1
write.csv(table1,"output/exhibits/figure_data/table1.csv")


#-------------------------------
# Regression coefficients    ####
rm(list= ls()[!(ls() %in% c(Keep.List))])
res <- as.data.frame(readRDS("output/summary/summary_piecewise.rds"))
res <- res
res <- res[res$crop %in% "hay_alfalfa",]
res <- res[res$climate_base %in% "1991_2020",]
res <- res[res$period %in% c(0,104:112),c("name","period","Estimate","StdError","t_value","p_value")]
res <- res[!grepl("trend",res$name),]
write.csv(res,"output/exhibits/figure_data/regression_coefficients.csv")

#-------------------------------
# Cluster knot outcomes (map + panel ring)   ####
# Agro-climatic cluster partition (004_..._knots_cluster.R) shown as a labeled
# county map framed by per-cluster GROSS-MEAN outcome panels. Each panel has one
# bar per cluster plus a full-sample reference bar. 5 columns x 4 rows, map centered.
rm(list= ls()[!(ls() %in% c(Keep.List))])
library(ggplot2); library(patchwork)

byz <- as.data.frame(readRDS("output/optimal_knots_cluster_byzone.rds"))
byz <- byz[order(byz$cluster), ]
if (is.null(byz$cluster_label)) byz$cluster_label <- paste("Cluster", byz$cluster)
clw <- as.data.frame(readRDS("output/knot_clusters.rds")$clusters)[, c("county_fips","cluster")]
ckc <- as.data.frame(readRDS("output/optimal_knots_cluster.rds"))
main_period <- as.integer(unique(ckc$period))[1]
Tmin_star   <- as.integer(unique(ckc$Tmin_national))[1]
Tmax_star   <- as.integer(unique(ckc$Tmax_national))[1]

# national slope / R / cv reference (from 003)
nk <- as.data.frame(readRDS("output/optimal_knots.rds"))
nk <- nk[nk$crop %in% "hay_alfalfa" & nk$target_periods %in% c(104:112, 0), ]
nk <- nk[order(nk$cv_error), ][1, ]

# per-county production & cattle joined to cluster
sp <- as.data.frame(readRDS("data/spatial_representation.rds"))
sp$county_fips <- stringr::str_pad(as.character(sp$fip), 5, pad = "0")
sp <- merge(sp[, c("county_fips","yield","inventory","area")], clw, by = "county_fips")

# degree-day EXPOSURE means: each cluster at its own knot, full sample at national knot
pan <- build_hay_weather_panel(crop = "hay_alfalfa", target_periods = main_period,
                               prism_weather_directory = "data/prism_weather")
data.table::setDT(pan); pan[, county_fips := stringr::str_pad(as.character(county_fips), 5, pad = "0")]
pan <- merge(pan, data.table::as.data.table(clw), by = "county_fips")
padc     <- function(x) stringr::str_pad(x, 2, pad = "0")
dd_means <- function(dt, tmin, tmax){
  dt <- data.table::copy(dt)
  dt[, D1 := pmax(get("dday00") - get(paste0("dday", padc(tmin))), 0)]
  dt[, D2 := pmax(get(paste0("dday", padc(tmin))) - get(paste0("dday", padc(tmax))), 0)]
  dt[, D3 := pmax(get(paste0("dday", padc(tmax))), 0)]
  c(DD1 = mean(dt$D1, na.rm = TRUE), DD2 = mean(dt$D2, na.rm = TRUE), DD3 = mean(dt$D3, na.rm = TRUE))
}
ddm      <- t(sapply(byz$cluster, function(cc){
  kt <- byz[byz$cluster == cc, ]; dd_means(pan[cluster == cc], kt$Tmin, kt$Tmax) }))
ddm_full <- dd_means(pan, Tmin_star, Tmax_star)
ym <- tapply(sp$yield, sp$cluster, mean, na.rm = TRUE)
cm <- tapply(sp$inventory, sp$cluster, mean, na.rm = TRUE)
am <- tapply(sp$area, sp$cluster, sum, na.rm = TRUE)   # total alfalfa acres (1,000 acres) per cluster

# one row per series (clusters + full sample); slopes scaled x1000 for legible labels
co <- data.frame(
  series   = byz$cluster_label,
  Tmin = byz$Tmin, Tmax = byz$Tmax,
  DD1_mean = ddm[, "DD1"], DD2_mean = ddm[, "DD2"], DD3_mean = ddm[, "DD3"],
  DD1_slope = byz$DD1 * 1000, DD2_slope = byz$DD2 * 1000, DD3_slope = byz$DD3 * 1000,
  yield = as.numeric(ym[as.character(byz$cluster)]),
  cattle = as.numeric(cm[as.character(byz$cluster)]),
  acres = as.numeric(am[as.character(byz$cluster)]),
  n_county = byz$n_county, cv = byz$cv_error, R = byz$R, stringsAsFactors = FALSE)
co <- rbind(co, data.frame(
  series = "Full sample", Tmin = Tmin_star, Tmax = Tmax_star,
  DD1_mean = ddm_full["DD1"], DD2_mean = ddm_full["DD2"], DD3_mean = ddm_full["DD3"],
  DD1_slope = nk$DD1 * 1000, DD2_slope = nk$DD2 * 1000, DD3_slope = nk$DD3 * 1000,
  yield = mean(sp$yield, na.rm = TRUE), cattle = mean(sp$inventory, na.rm = TRUE),
  acres = sum(sp$area, na.rm = TRUE),
  n_county = nrow(clw), cv = nk$cv_error, R = nk$R, stringsAsFactors = FALSE))
co$series <- factor(co$series, levels = c(byz$cluster_label, "Full sample"))
write.csv(co, "output/exhibits/figure_data/cluster_knot_outcomes.csv", row.names = FALSE)

# consistent palette: clusters + a neutral grey full sample (shared by map and bars)
base_cols <- c("#2a78d6","#eb6834","#008300","#4a3aa7","#1baf7a","#eda100","#e34948","#c05780")
pal     <- setNames(c(base_cols[seq_len(nrow(byz))], "#8a8a86"), levels(co$series))
map_pal <- setNames(base_cols[seq_len(nrow(byz))], as.character(byz$cluster))

mk <- function(col, title, dp){
  d <- co; d$val <- d[[col]]; d$vj <- ifelse(d$val >= 0, -0.35, 1.25)
  ggplot(d, aes(stats::reorder(series, -val), val, fill = series)) +   # bars sorted within panel
    geom_col(width = 0.78) +
    geom_hline(yintercept = 0, linewidth = 0.2, colour = "grey55") +
    geom_text(aes(label = formatC(val, format = "f", digits = dp), vjust = vj),
              size = 2.3, colour = "grey25") +
    scale_fill_manual(values = pal, guide = "none") +
    scale_y_continuous(expand = expansion(mult = c(0.16, 0.20))) +
    labs(title = title, x = NULL, y = NULL) +
    theme_bw(base_size = 8) +
    theme(axis.text.x = element_blank(), axis.ticks.x = element_blank(),
          panel.grid = element_blank(), plot.title = element_text(size = 7.6),
          plot.margin = margin(2, 2, 2, 2))
}
p_tmin  <- mk("Tmin",     "Lower knot: Tmin (°C)", 0)
p_tmax  <- mk("Tmax",     "Upper knot: Tmax (°C)", 0)
p_dd1v  <- mk("DD1_mean", "Degree-day (below Tmin): mean", 0)
p_dd2v  <- mk("DD2_mean", "Degree-day (Tmin to Tmax): mean", 0)
p_dd3v  <- mk("DD3_mean", "Degree-day (above Tmax): mean", 0)
p_dd1s  <- mk("DD1_slope","Degree-day (below Tmin): slope (x10^3)", 2)
p_dd2s  <- mk("DD2_slope","Degree-day (Tmin to Tmax): slope (x10^3)", 2)
p_dd3s  <- mk("DD3_slope","Degree-day (above Tmax): slope (x10^3)", 2)
p_yield <- mk("yield",    "Yield (tons/acre)", 1)
p_catt  <- mk("cattle",   "Cattle inventory (1,000 head)", 0)
p_acres <- mk("acres",    "Alfalfa acres in cluster (1,000 acres)", 0)
p_cv    <- mk("cv",       "Cross-validation error", 2)
p_r2    <- mk("R",        "Within-county fit (R^2)", 2)

# legend as its own panel (map now carries colour only, no region text)
leg_df <- data.frame(y = rev(seq_len(nlevels(co$series))), lab = levels(co$series))
leg_df$wrapped <- stringr::str_wrap(leg_df$lab, 14)
gleg <- ggplot(leg_df, aes(0, y, colour = lab)) +
  geom_point(size = 3.2, shape = 15) +
  geom_text(aes(x = 0.1, label = wrapped), hjust = 0, size = 2.2, colour = "grey20", lineheight = 0.85) +
  scale_colour_manual(values = pal, guide = "none") +
  scale_x_continuous(limits = c(-0.08, 1.4)) +
  scale_y_continuous(expand = expansion(mult = c(0.3, 0.3))) +
  labs(title = "Series") +
  theme_void(base_size = 8) + theme(plot.title = element_text(size = 7.6, hjust = 0))

# centered county cluster map (county fill, STATE boundaries overlaid)
cnt_sf <- sf::st_as_sf(Counties)
cnt_sf$county_fips <- stringr::str_pad(as.character(cnt_sf$fip), 5, pad = "0")
cnt_sf <- merge(cnt_sf, clw, by = "county_fips", all.x = TRUE)
st_sf  <- sf::st_as_sf(States)
gmap <- ggplot() +
  geom_sf(data = cnt_sf, aes(fill = factor(cluster)), colour = NA) +
  geom_sf(data = st_sf, fill = NA, colour = "grey30", linewidth = 0.18) +
  scale_fill_manual(values = map_pal, na.value = "grey85", guide = "none") +
  theme_void()

# 5 columns x 4 rows; map spans the centre (rows 2-3, cols 2-4)
design <- "ABCDE
FOOOH
GOOOI
JKLMN"
fig <- gleg + p_tmin + p_tmax + p_dd1v + p_dd1s + p_dd2v + p_dd2s + p_dd3v + p_dd3s +
  p_yield + p_catt + p_acres + p_cv + p_r2 + gmap +
  plot_layout(design = design)
ggsave("output/exhibits/cluster_knot_outcomes.png", fig, width = 13, height = 8, dpi = 300)

#-------------------------------
# Nonlinear Relation         ####
rm(list= ls()[!(ls() %in% c(Keep.List))])
relation <- as.data.frame(readRDS("output/summary/summary_relation.rds"))
relation <- relation
relation <- relation[relation$crop %in% "hay_alfalfa",]
relation <- relation[relation$climate_base %in% "1991_2020",]
relation <- relation[relation$period %in% c(105:110),c("crop","period","Temp","Piece","PieceSE")]

exposure <- as.data.frame(readRDS("output/summary/summary_exposure.rds"))
exposure <- exposure
exposure <- exposure[exposure$crop %in% "hay_alfalfa",]
exposure <- exposure[exposure$climate_base %in% "1991_2020",]
exposure <- exposure[exposure$period %in% c(105:110),c("crop","period","Temp","exp","exp_sd")]

Result <- dplyr::inner_join(exposure,relation,by=c("Temp","crop","period"))

Result$x <- Result$Temp
#Result <- Result[Result$Temp <= 35,]
#plot(Result$Temp,Result$Piece)

Xlab <- unique(Result[c("x","Temp")])
Xlab <- Xlab[Xlab$Temp %in% seq(1,45,5),]

colors <- c("Step Function" = "#FF00FF", "95% confidence interval for **" = "thistle")
colors[preferred_lab] <- "purple"


# Label each window; the preferred window carries the "**" emphasis, the rest are plain.
Result$months <- ifelse(Result$period %in% 0, "all months",
                 ifelse(Result$period %in% preferred_period, preferred_lab,
                        paste0(Result$period - 100, " months")))
.oth <- sort(unique(Result$period[!Result$period %in% c(0, preferred_period)]))
Result$months <- factor(Result$months,
                        levels = c(preferred_lab, paste0(.oth - 100, " months"), "all months"))
Result <- Result[!Result$months %in% NA,]
write.csv(Result,"output/exhibits/figure_data/Nonlinear_Relation.csv")

fig.a <- ggplot(data=Result,aes(x=x, y=Piece,group=period)) +
  geom_hline(yintercept = 0,size = 0.2,color = "black") +
  geom_ribbon(data=Result[Result$period %in% preferred_period,], aes(x=x, ymin = (Piece-1.96*PieceSE) ,ymax = (Piece+1.96*PieceSE),
                               fill="95% confidence interval for **"),color="white") +
  geom_line(data=Result[Result$period %in% preferred_period,],aes(x=x,y=Piece,color=preferred_lab),size = 0.8,linetype="solid")  +
  geom_line(data=Result[!Result$period %in% preferred_period,],aes(x=x,y=Piece,color=months),size = 0.8,linetype="aa",alpha=0.5)  +
  labs(title="",x = "", y = "log-yield (Mt/ha)\n",fill ="",color ="growing season weather\naggregation window", caption = "") +
  scale_color_manual(values = c(colorRampPalette(c("yellow","darkgreen"))(length(unique(Result[!Result$period %in% preferred_period,"months"]))),"purple")) +
  scale_fill_manual(values = colors) +
  scale_x_continuous(breaks=Xlab$x, labels=Xlab$Temp) +
  scale_y_continuous(breaks=seq(-0.20,0.20,0.005), labels=sprintf(seq(-0.20,0.20,0.005),fmt="%#.3f")) +
  ers_theme() +
  theme_bw() +
  theme(panel.grid.major = element_blank(),panel.grid.minor = element_blank(), axis.ticks.x = element_blank()) +
  theme(legend.position=c(0.35,0.40),legend.key.size = unit(0.43,"cm"))+
  theme(legend.text=element_text(size=11),
        legend.title=element_text(size=11),
        axis.title.y=element_text(size=11),
        axis.title.x=element_blank(),
        #axis.text.x = element_blank(), #
        axis.text.y = element_text(size = 8, colour="black"),
        plot.caption = element_text(size=11,hjust = 0 ,vjust = 0, face = "italic"),
        strip.text = element_text(size = 12),
        strip.background = element_rect(fill = "white", colour = "black", size = 1))

fig.b <- ggplot(Result[Result$period %in% preferred_period,]) +
  geom_bar(aes(x=x, y=exp),stat="identity", color = "black", fill = "mediumpurple") +
  #geom_errorbar(aes(x=x,ymin = exp - exp_sd*1.96, ymax = exp + exp_sd*1.96),color = "red") +
  labs(title="",x = "\nTemperature (°C)", y = "",fill ="",color ="", caption = "") +
  scale_x_continuous(breaks=Xlab$x, labels=Xlab$Temp) +
  scale_y_continuous(breaks=0, labels="",
                     sec.axis = sec_axis(~., name = "Exposure (Days)\n",breaks=seq(0,40,2))) +
  ers_theme() +
  theme_bw() +
  theme(panel.grid.major = element_blank(),panel.grid.minor = element_blank()) +
  theme(legend.position="none")+
  theme(legend.text=element_text(size=11),
        legend.title=element_blank(),
        axis.title.y=element_text(size=11),
        axis.title.x=element_text(size=11),
        axis.text.x = element_text(size = 8),
        axis.text.y = element_text(size = 8),
        plot.caption = element_text(size=11,hjust = 0 ,vjust = 0, face = "italic"),
        strip.text = element_text(size = 12),
        strip.background = element_rect(fill = "white", colour = "black", size = 1))

fig <- cowplot::plot_grid(fig.a + theme(plot.margin=unit(c(0.05,0.5,-1,0.5), "cm")),
                          fig.b + theme(plot.margin=unit(c(-1,0.5,1,0.5), "cm")),
                          ncol=1, align="v",rel_heights=c(1, 0.60),
                          greedy=F)

ggsave("output/exhibits/Nonlinear_Relation.png", fig, dpi = 600,width = 6.7, height = 7.5)
#-------------------------------
# Impacts mean               ####
rm(list= ls()[!(ls() %in% c(Keep.List))])
data <- as.data.frame(readRDS("output/summary/summary_impact_yield.rds"))
data$est <- data$Estimate
data$se  <- data$Estimate_sd
data <- data[data$county_code %in% 0,]
data <- data[data$state_code %in% 0,]
data <- data[data$region %in% "",]

data_main <- data[(data$crop %in% "hay_alfalfa" & data$period %in% preferred_period & data$climate_base %in% "1991_2020"),]
data_main <- data_main

data_crop <- data[(data$warming_scenario %in% 1.0 & data$crop %in% c("hay_other","hay_alfalfa") & data$period %in% preferred_period & data$climate_base %in% "1991_2020"),]
data_crop <- data_crop

data_year <- data[(data$warming_scenario %in% 1.0 &data$crop %in% "hay_alfalfa" & data$period %in% preferred_period & data$climate_base %in% c("1991_2020","1981_2010","1971_2000","1961_1990")),]
data_year <- data_year

data_wind <- data[(data$warming_scenario %in% 1.0 &data$crop %in% "hay_alfalfa" & data$period %in% c(105:110) & data$climate_base %in% "1991_2020"),]
data_wind <- data_wind

data_main$type <- "Impact by warming scenario"
data_year$type <- "1+°C warming by climate baseline"
data_wind$type <- "1+°C warming by season window"
data_crop$type <- "1+°C warming by hay type"

data_main <- data_main[order(data_main$warming_scenario),]
data_crop <- data_crop[order(data_crop$est),]
data_year <- data_year[order(data_year$est),]
data_main <- data_main[order(data_main$est),]

data_wind$period <- ifelse(data_wind$period %in% 0,113,data_wind$period)
data_wind <- data_wind[order(-data_wind$period),]

data_wind$x_name <- paste0(data_wind$period-100," months")
data_wind$x_name <- ifelse(data_wind$x_name %in% "13 months","All months",data_wind$x_name)
data_year$x_name <- data_year$climate_base
data_crop$x_name <- factor(data_crop$crop,levels = c("hay_other","hay_all","hay_alfalfa"),
                           labels = c("Non-alfalfa","Combined","Alfalfa"))
data_main$x_name <- paste0(format(as.numeric(as.character(data_main$warming_scenario)), nsmall = 1) ,"+°C")

data_cropy <- data_crop
data_windy <- data_wind
data_yeary <- data_year
data_mainy <- data_main

data_cropy$x <- c(1:nrow(data_cropy)) + 2
data_windy$x <- max(data_cropy$x) + c(1:nrow(data_windy)) + 2
data_yeary$x <- max(data_windy$x) + c(1:nrow(data_yeary)) + 2
data_mainy$x <- max(data_yeary$x) + c(1:nrow(data_mainy)) + 2

datay <- rbind(data_cropy,data_mainy,data_yeary,data_windy)
datay <- unique(datay[c("p","theta","longlat","DistName","kernel","crop","period","climate_base","warming_scenario","x_name","x")])
data <- rbind(data_crop,data_main,data_year,data_wind)
data <- dplyr::inner_join(data,datay)

data_text <- doBy::summaryBy(x~type,data=data,FUN=max,keep.names = T,na.rm=T)

Fig04 <- ggplot(data,aes(x=x,y=est,group=1)) +
  geom_hline(yintercept = 0,size = 0.2,color = "blue") +
  geom_errorbar(aes(ymin = est - se*1.96, ymax = est + se*1.96),color="purple") +
  geom_point(color="purple") +
  geom_text(data=data_text,aes(x=x+1.5,y=-13,label=type, hjust = 0), size = 2.5, col = "black",
            stat = "identity",check_overlap = TRUE) +
  scale_y_continuous(breaks=seq(-60,0,2.5),labels=format(round(seq(-60,0,2.5),2), nsmall = 2)) +
  scale_x_continuous(breaks=data$x,labels=data$x_name) +
  labs(title="", x="", y ="\nYield impact (%)",caption = "") +
  myTheme +
  theme(axis.text.x = element_text(size=6,color="black"),
        axis.text.y = element_text(size=7,color="black"),
        strip.text = element_text(size = 9),
        legend.position="none")+ coord_flip()

write.csv(data,"output/exhibits/figure_data/impacts_mean.csv")
ggsave("output/exhibits/impacts_mean.png", Fig04, dpi = 600,width = 6.7, height = 7.5)
#ggsave("output/exhibits/impacts_mean.png", Fig04, dpi = 600,width = 7, height = 6)

#-------------------------------
# Impacts Spatial            ####
rm(list= ls()[!(ls() %in% c(Keep.List))])
data <- as.data.frame(readRDS("output/summary/summary_impact_yield.rds"))
data$est <- data$Estimate
data$se  <- data$Estimate_sd
data <- data[!data$county_code %in% 0,]
data <- data[(data$crop %in% "hay_alfalfa" & data$period %in% preferred_period & data$climate_base %in% "1991_2020"),]
data <- data[c("p","theta","longlat","DistName","kernel","crop","period","climate_base","warming_scenario","fip","est","se")]

data$SimCat<-factor(data$warming_scenario,levels=unique(data$warming_scenario),
                    labels = paste0("+",format(unique(data$warming_scenario), nsmall = 1)," °C"))

merged_spat_vector <- terra::merge(Counties, data, by="fip")
sf_object <- sf::st_as_sf(merged_spat_vector)

cutlist <- unique(c(min(sf_object$est,na.rm=T),0,max(sf_object$est,na.rm=T),
                    quantile(unique(sf_object$est[! sf_object$est %in% c(NA,Inf,-Inf,NaN) | sf_object$est <0.0001]), probs = seq(0.1,1,0.10))))

sf_object$Value <- cut(sf_object$est,cutlist)
table(sf_object$Value)
sf_object <- sf_object[!is.na(sf_object$Value), ]

Fig <- ggplot() +
  geom_sf(data = sf::st_as_sf(States), colour = "black", fill = "darkred",size = 0.1) +
  geom_sf(data = sf_object,aes(fill = Value), colour = NA,size = 0.2) +
  geom_sf(data = sf::st_as_sf(States), colour = "black", fill = NA,size = 0.1) +
  scale_fill_manual("Percentage",drop=FALSE, values=c(colorRampPalette(c("red","yellow"))(length(unique(as.character(sf_object$Value)))-2),
                                         "green","darkgreen"), na.value="#EEEEEE",
                    name="") +
  labs(title= "", x = "", y = "",fill ="Percentage",caption = "") +
  guides(fill = guide_legend(nrow=3)) +
  facet_wrap(vars(SimCat),ncol=2) +
  ers_theme() +
  theme_bw() +
  theme(panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        axis.ticks.x = element_blank(),
        axis.ticks.y = element_blank(),
        plot.title=element_text(size=8),
        legend.position="bottom",
        legend.background = element_blank(),
        legend.key.size = unit(0.2,"cm"),
        legend.text=element_text(size=4),
        legend.title=element_text(size=6),
        axis.title.y = element_blank(),
        axis.title.x = element_blank(),
        axis.text.x = element_blank(), #
        axis.text.y = element_blank(),
        strip.text = element_text(size=9),
        strip.background = element_blank())+coord_sf()

ggsave(paste0("output/exhibits/spatial_impact_yield.png"), Fig, dpi = 600,width = 6.7, height = 7.5)

#-------------------------------
# Regional estimates         ####
rm(list= ls()[!(ls() %in% c(Keep.List))]);gc()

yield <- as.data.frame(readRDS("output/summary/summary_impact_yield.rds"))
yield$impact_yield <- yield$Estimate
yield$fip<-as.character(paste0(stringr::str_pad(as.numeric(as.character(yield$state_code)), 2, pad = "0"),
                               stringr::str_pad(as.numeric(as.character(yield$county_code)), 3, pad = "0")))
yield <- yield[!yield$state_code %in% 0,]
yield <- yield[yield$warming_scenario %in% 1.0,]
yield <- yield[(yield$crop %in% "hay_alfalfa" & yield$period %in% preferred_period & yield$climate_base %in% "1991_2020"),]
yield <- yield
yield <- yield[c("state_code","county_code","fip","impact_yield")]
data <- yield

avail <- as.data.frame(readRDS("output/summary/summary_availability.rds"))
avail$level_avail <- avail$Estimate
avail <- avail[!avail$fip %in% "0000",]
avail <- avail[avail$warming_scenario %in% c(0.0),]
avail <- avail[(avail$crop %in% "hay_alfalfa" & avail$period %in% preferred_period & avail$climate_base %in% "1991_2020"),]
avail <- avail
avail <- avail[c("fip","level_avail")]
data <- dplyr::inner_join(data,avail,by="fip")

res <- as.data.frame(readRDS("output/summary/summary_associations.rds"))
res <- res[!res$fip %in% "00000",]
res <- res[(res$crop %in% "hay_alfalfa" & res$period %in% preferred_period & res$climate_base %in% "1991_2020"),
           c("p", "theta", "longlat", "DistName", "kernel","fip","name","est")]
res <- res |> tidyr::spread(name, est)
res <- res
res <- res[c("fip","avail00","prod00","prod00_LM")]
names(res) <- c("fip","b_avail00","b_prod00","b_prod00_LM")
data <- dplyr::inner_join(data,res,by="fip")

avail <- as.data.frame(readRDS("output/summary/summary_availability.rds"))
avail$impact_avail <- avail$Estimate
avail <- avail[!avail$fip %in% "0000",]
avail <- avail[avail$warming_scenario %in% 1.0,]
avail <- avail[(avail$crop %in% "hay_alfalfa" & avail$period %in% preferred_period & avail$climate_base %in% "1991_2020"),]
avail <- avail
avail <- avail[c("fip","impact_avail")]
data <- dplyr::inner_join(data,avail,by="fip")

catle <- as.data.frame(readRDS("output/summary/summary_cattle.rds"))
catle$impact_cattle <- catle$Estimate
catle <- catle[!catle$fip %in% "0000",]
catle <- catle[catle$warming_scenario %in% 1.0,]
catle <- catle[(catle$crop %in% "hay_alfalfa" & catle$period %in% preferred_period & catle$climate_base %in% "1991_2020"),]
catle <- catle
catle <- catle[c("fip","cattleA","cattleB","cattleC")]
data <- dplyr::inner_join(data,catle,by="fip")

data <- dplyr::inner_join(data,stnames[c("state_code","State.Name","State.Abbreviation" )],by="state_code")

data$region <- ifelse(
  data$State.Abbreviation %in% c("CT", "DE", "ME", "MD", "MA", "NH", "NJ", "NY", "PA", "RI", "VT"),"Northeast",
  NA)
data$region <- ifelse(
  data$State.Abbreviation %in% c("IL", "IN", "IA", "KS", "MI", "MN", "MO", "NE", "ND", "OH", "SD", "WI"),"Midwest",
  data$region)
data$region <- ifelse(
  data$State.Abbreviation %in% c("AL", "AR", "FL", "GA", "KY", "LA", "MS", "NC", "OK", "SC", "TN", "TX", "VA", "WV"),"South",
  data$region)
data$region <- ifelse(
  data$State.Abbreviation %in% c("AZ", "CA", "CO", "HI", "ID", "MT", "NV", "NM", "OR", "UT", "WA", "WY"),"West",
  data$region)
saveRDS(data,"output/regional_estimates.rds")

#-------------------------------
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
  cn <- as.data.frame(gw_consensus_scalar(
    value_dt = assoc[assoc$name %in% nm, ], unit_col = "fip",
    geometry = Counties, value_col = "est", agg_fun = stats::median,
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

ggsave("output/exhibits/associations.png", Fig06, dpi = 600, width = 9, height = 6)

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

avail_cons <- as.data.frame(gw_consensus_scalar(
  value_dt = avail, unit_col = "fip", geometry = Counties,
  value_col = "est", agg_fun = stats::median, queen_smooth = FALSE))

sf_object <- sf::st_as_sf(terra::merge(Counties, avail_cons, by = "fip"))
sf_object <- sf_object[is.finite(sf_object$consensus), ]

# equal-count quantile bins (6 classes, ~equal county counts) for max spatial variation
brks <- unique(stats::quantile(sf_object$consensus, probs = seq(0, 1, length.out = 7), na.rm = TRUE))
sf_object$Value <- cut(sf_object$consensus, breaks = brks, include.lowest = TRUE, dig.lab = 4)
sf_object <- sf_object[!is.na(sf_object$Value), ]

Fig_avail <- ggplot() +
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
        legend.background = element_blank(),
        legend.position=c(0.17,0.18),
        legend.key.size = unit(0.2,"cm"),
        legend.text = element_text(size=7),
        legend.title = element_text(size=7),
        axis.title.y = element_text(size=8),
        axis.title.x = element_text(size=8),
        axis.text.x  = element_blank(),
        axis.text.y  = element_blank(),
        strip.text = element_text(size = 8),
        strip.background = element_blank())+coord_sf()

#-------------------------------
# cattle                     ####
# "Responsiveness of cattle inventory to alfalfa availability" is the GWR coefficient
# b_avail00, which lives in summary_associations (est where name == "avail00").
corr_cat <- as.data.frame(readRDS("output/summary/summary_associations.rds"))
corr_cat <- corr_cat[corr_cat$name %in% "avail00" &
                       corr_cat$crop %in% "hay_alfalfa" &
                       corr_cat$period %in% preferred_period &
                       corr_cat$climate_base %in% "1991_2020", ]
corr_cat <- corr_cat[!corr_cat$fip %in% c("0","00000"),]
corr_cat <- corr_cat[is.finite(corr_cat$est), ]

# 20-spec consensus of the responsiveness coefficient via gwkit.
corr_cons <- as.data.frame(gw_consensus_scalar(
  value_dt = corr_cat, unit_col = "fip", geometry = Counties,
  value_col = "est", agg_fun = stats::median, queen_smooth = FALSE))

# --- (folded from 101) continuous diverging map + consensus CSV + sign-agreement diagnostic ---
# Same per-county consensus as Fig_corr below, rendered on a continuous diverging scale
# (symmetric 98% clip) and written with its spread / sign-agreement columns.
write.csv(corr_cons, "output/exhibits/figure_data/consensus_avail00_responsiveness.csv",
          row.names = FALSE)
sf_div  <- sf::st_as_sf(terra::merge(Counties, corr_cons, by = "fip"))
sf_div  <- sf_div[is.finite(sf_div$consensus), ]
lim_div <- stats::quantile(abs(sf_div$consensus), 0.98, na.rm = TRUE)   # symmetric clip
Fig_div <- ggplot() +
  geom_sf(data = sf::st_as_sf(States), colour = "black", fill = "grey95", size = 0.2) +
  geom_sf(data = sf_div, aes(fill = consensus), colour = NA, size = 0.2) +
  geom_sf(data = sf::st_as_sf(States), colour = "black", fill = NA, size = 0.2) +
  scale_fill_gradient2(low = "#B2182B", mid = "#F7F7F7", high = "#2166AC",
                       midpoint = 0, limits = c(-lim_div, lim_div), oob = scales::squish,
                       name = "Consensus\ncoefficient") +
  labs(title = "Responsiveness of cattle inventory to alfalfa availability",
       subtitle = "Per-county median consensus across all GW specs",
       x = "", y = "", caption = "gwkit::gw_consensus_scalar()") +
  ers_theme() + theme_bw() +
  theme(panel.grid = element_blank(), axis.text = element_blank(), axis.ticks = element_blank(),
        legend.position = "right", plot.title = element_text(size = 11),
        plot.subtitle = element_text(size = 9)) + coord_sf()
ggsave("output/exhibits/consensus_avail00_responsiveness.png", Fig_div, dpi = 600,
       width = 7.5, height = 5)
message("consensus_avail00: ", nrow(sf_div), " counties; median sign-agreement = ",
        round(stats::median(corr_cons$sign_agreement, na.rm = TRUE), 3))

sf_object <- sf::st_as_sf(terra::merge(Counties, corr_cons, by = "fip"))
sf_object <- sf_object[is.finite(sf_object$consensus), ]

brks <- unique(stats::quantile(sf_object$consensus, probs = seq(0, 1, length.out = 7), na.rm = TRUE))
sf_object$Value <- cut(sf_object$consensus, breaks = brks, include.lowest = TRUE, dig.lab = 4)
sf_object <- sf_object[!is.na(sf_object$Value), ]

Fig_corr <- ggplot() +
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
        legend.background = element_blank(),
        legend.position=c(0.17,0.18),
        legend.key.size = unit(0.2,"cm"),
        legend.text = element_text(size=7),
        legend.title = element_text(size=7),
        axis.title.y = element_text(size=8),
        axis.title.x = element_text(size=8),
        axis.text.x  = element_blank(),
        axis.text.y  = element_blank(),
        strip.text = element_text(size = 8),
        strip.background = element_blank())+coord_sf()

marg <- c(-0.01,0.05,-0.01,0.05)

Fig <- cowplot::plot_grid(
  Fig_avail + theme(plot.margin=unit(marg, "cm")) ,
  Fig_corr + theme(plot.margin=unit(marg, "cm")) ,
  ncol=1, align="v",rel_heights=c(1,1),
  greedy=F)

ggsave("output/exhibits/availability_cattle.png", Fig, dpi = 600,width = 6.5, height = 9)

#-------------------------------
# Predicted county-level impacts (spatial; not in article yet)   ####
# Yield + availability + cattle % impact by warming scenario, per county. Availability
# and cattle are GW outcomes -> per-county 20-spec consensus per scenario via gwkit;
# yield is spec-invariant, so its consensus equals its value.
rm(list = ls()[!(ls() %in% c(Keep.List))])

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
    cn <- as.data.frame(gw_consensus_scalar(
      ds, unit_col = "fip", geometry = Counties, value_col = "est",
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
ggsave("output/exhibits/predicted_county_impacts.png", Fig6, dpi = 600, width = 9, height = 6)

#-------------------------------
# Predicted impacts (mean; dot plot; not in article yet)   ####
# National impact by hay type, season window, climate baseline, and warming scenario,
# faceted by outcome (yield / availability / cattle). The v00 by-kernel and
# by-distance-metric panels are dropped: 006 now reduces the 50 GW specs to one
# per-county consensus, so there is no per-spec dimension to display.
rm(list = ls()[!(ls() %in% c(Keep.List))])

.common <- c("crop","period","climate_base","warming_scenario","est","se","outcome")

yield <- data.table::as.data.table(readRDS("output/summary/summary_impact_yield.rds"))
yield$est <- yield$Estimate; yield$se <- yield$Estimate_sd
yield <- yield[county_code %in% 0 & state_code %in% 0 & region %in% "", ]
yield$outcome <- "(a) Alfalfa yield"
yield <- as.data.frame(yield)[, .common]

# summary_availability / summary_cattle are per-county with no national aggregate row, so
# build the national series by averaging the per-county impact within each cell/scenario
# (se = across-county sd). Yield keeps its stored national bootstrap se.
.national <- function(path, val, outcome_label) {
  d <- data.table::as.data.table(readRDS(path))
  d <- d[!fip %in% c("0", "00000")]
  d[, est := suppressWarnings(as.numeric(get(val)))]
  d[, warming_scenario := suppressWarnings(as.numeric(as.character(warming_scenario)))]
  d <- d[is.finite(est)]
  keys <- intersect(c("crop","period","climate_base","warming_scenario"), names(d))
  ag <- d[, .(est = mean(est, na.rm = TRUE), se = stats::sd(est, na.rm = TRUE)), by = keys]
  ag[, outcome := outcome_label]
  as.data.frame(ag)[, .common]
}
avail  <- .national("output/summary/summary_availability.rds", "Estimate", "(b) Alfalfa availability")
cattle <- .national("output/summary/summary_cattle.rds",       "cattleA",  "(c) Cattle inventory")

data <- rbind(yield, avail, cattle)
data$warming_scenario <- suppressWarnings(as.numeric(as.character(data$warming_scenario)))

data_main <- data[(data$crop %in% "hay_alfalfa" & data$period %in% preferred_period & data$climate_base %in% "1991_2020"),]
data_crop <- data[(data$warming_scenario %in% 1.0 & data$crop %in% c("hay_other","hay_alfalfa") & data$period %in% preferred_period & data$climate_base %in% "1991_2020"),]
data_year <- data[(data$warming_scenario %in% 1.0 & data$crop %in% "hay_alfalfa" & data$period %in% preferred_period & data$climate_base %in% c("1991_2020","1981_2010","1971_2000","1961_1990")),]
data_wind <- data[(data$warming_scenario %in% 1.0 & data$crop %in% "hay_alfalfa" & data$period %in% c(105:110) & data$climate_base %in% "1991_2020"),]

data_main$type <- "Impact by warming scenario"
data_year$type <- "1+°C warming by climate baseline"
data_wind$type <- "1+°C warming by season window"
data_crop$type <- "1+°C warming by hay type"

data_wind$period <- ifelse(data_wind$period %in% 0,113,data_wind$period)
data_wind <- data_wind[order(-data_wind$period),]

data_wind$x_name <- paste0(data_wind$period-100," months")
data_wind$x_name <- ifelse(data_wind$x_name %in% "13 months","All months",data_wind$x_name)
data_year$x_name <- data_year$climate_base
data_crop$x_name <- as.character(factor(data_crop$crop,levels = c("hay_other","hay_all","hay_alfalfa"),
                           labels = c("Non-alfalfa","Combined","Alfalfa")))
data_main$x_name <- paste0(format(as.numeric(as.character(data_main$warming_scenario)), nsmall = 1) ,"+°C")

# x positions taken from the yield column, then shared across all three facets.
data_cropy <- data_crop[data_crop$outcome %in% "(a) Alfalfa yield",]
data_windy <- data_wind[data_wind$outcome %in% "(a) Alfalfa yield",]
data_yeary <- data_year[data_year$outcome %in% "(a) Alfalfa yield",]
data_mainy <- data_main[data_main$outcome %in% "(a) Alfalfa yield",]

data_cropy$x <- 1:nrow(data_cropy) + 2
data_windy$x <- max(data_cropy$x) + c(1:nrow(data_windy)) + 2
data_yeary$x <- max(data_windy$x) + c(1:nrow(data_yeary)) + 2
data_mainy$x <- max(data_yeary$x) + c(1:nrow(data_mainy)) + 2

datay <- rbind(data_cropy,data_mainy,data_yeary,data_windy)
datay <- unique(datay[c("crop","period","climate_base","warming_scenario","x_name","x")])
data <- rbind(data_crop,data_main,data_year,data_wind)
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
ggsave("output/exhibits/predicted_impacts_mean.png", Fig04, dpi = 600, width = 10, height = 7.5)
