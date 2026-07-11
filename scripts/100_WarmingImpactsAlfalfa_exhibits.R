#-------------------------------
# Preliminaries              ####
rm(list=ls(all=TRUE))
library(ggplot2);library(terra);library(ggridges);library(gridExtra);library(gtable);library(data.table)
if(Sys.info()['sysname'] =="Windows"){library(gganimate);library(magick)}
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
  crop = "hay_alfalfa", target_periods = 107,   # 7-month preferred window (Table 1 & Figure 1)
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
optimal_gw <- as.data.frame(readRDS("output/optimal_gw.rds"))
res <- as.data.frame(readRDS("output/summary/summary_piecewise.rds"))
res <- dplyr::inner_join(res,optimal_gw[1,])
res <- res[res$crop %in% "hay_alfalfa",]
res <- res[res$climate_base %in% "1991_2020",]
res <- res[res$period %in% c(0,104:112),c("name","period","Estimate","StdError","t_value","p_value")]
res <- res[!grepl("trend",res$name),]
write.csv(res,"output/exhibits/figure_data/regression_coefficients.csv")

#-------------------------------
# Nonlinear Relation         ####
rm(list= ls()[!(ls() %in% c(Keep.List))])
optimal_gw <- as.data.frame(readRDS("output/optimal_gw.rds"))
relation <- as.data.frame(readRDS("output/summary/summary_relation.rds"))
relation <- dplyr::inner_join(relation,optimal_gw[1,])
relation <- relation[relation$crop %in% "hay_alfalfa",]
relation <- relation[relation$climate_base %in% "1991_2020",]
relation <- relation[relation$period %in% c(105:110),c("crop","period","Temp","Piece","PieceSE")]

exposure <- as.data.frame(readRDS("output/summary/summary_exposure.rds"))
exposure <- dplyr::inner_join(exposure,optimal_gw[1,])
exposure <- exposure[exposure$crop %in% "hay_alfalfa",]
exposure <- exposure[exposure$climate_base %in% "1991_2020",]
exposure <- exposure[exposure$period %in% c(105:110),c("crop","period","Temp","exp","exp_sd")]

Result <- dplyr::inner_join(exposure,relation,by=c("Temp","crop","period"))

Result$x <- Result$Temp
#Result <- Result[Result$Temp <= 35,]
#plot(Result$Temp,Result$Piece)

Xlab <- unique(Result[c("x","Temp")])
Xlab <- Xlab[Xlab$Temp %in% seq(1,45,5),]

colors <- c("Step Function" = "#FF00FF", "5 months**" = "purple",
            "95% confidence interval for **" = "thistle")


Result$months <- as.numeric(as.character(factor(Result$period,levels = c(105,c(107,106,108:112,0)),labels = 1:9)))
Result$months <- factor(Result$months,levels = 1:9,labels = c("5 months**",paste0(c(6:8:12)," months"),"all months"))
Result <- Result[!Result$months %in% NA,]
write.csv(Result,"output/exhibits/figure_data/Nonlinear_Relation.csv")

fig.a <- ggplot(data=Result,aes(x=x, y=Piece,group=period)) +
  geom_hline(yintercept = 0,size = 0.2,color = "black") +
  geom_ribbon(data=Result[Result$period %in% 105,], aes(x=x, ymin = (Piece-1.96*PieceSE) ,ymax = (Piece+1.96*PieceSE),
                               fill="95% confidence interval for **"),color="white") +
  geom_line(data=Result[Result$period %in% 105,],aes(x=x,y=Piece,color="5 months**"),size = 0.8,linetype="solid")  +
  geom_line(data=Result[!Result$period %in% 105,],aes(x=x,y=Piece,color=months),size = 0.8,linetype="aa",alpha=0.5)  +
  labs(title="",x = "", y = "log-yield (Mt/ha)\n",fill ="",color ="growing season weather\naggregation window", caption = "") +
  scale_color_manual(values = c(colorRampPalette(c("yellow","darkgreen"))(length(unique(Result[!Result$period %in% 105,"months"]))),"purple")) +
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

fig.b <- ggplot(Result[Result$period %in% 105,]) +
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
optimal_gw <- as.data.frame(readRDS("output/optimal_gw.rds"))
data <- as.data.frame(readRDS("output/summary/summary_impact_yield.rds"))
data$est <- data$Estimate
data$se  <- data$Estimate_sd
data <- data[data$county_code %in% 0,]
data <- data[data$state_code %in% 0,]
data <- data[data$region %in% "",]

data_main <- data[(data$crop %in% "hay_alfalfa" & data$period %in% 105 & data$climate_base %in% "1991_2020"),]
data_main <- dplyr::inner_join(data_main,optimal_gw[1,c("p", "theta", "longlat", "DistName", "kernel")])

data_crop <- data[(data$warming_scenario %in% 1.0 & data$crop %in% c("hay_other","hay_alfalfa") & data$period %in% 105 & data$climate_base %in% "1991_2020"),]
data_crop <- dplyr::inner_join(data_crop,optimal_gw[1,c("p", "theta", "longlat", "DistName", "kernel")])

data_year <- data[(data$warming_scenario %in% 1.0 &data$crop %in% "hay_alfalfa" & data$period %in% 105 & data$climate_base %in% c("1991_2020","1981_2010","1971_2000","1961_1990")),]
data_year <- dplyr::inner_join(data_year,optimal_gw[1,c("p", "theta", "longlat", "DistName", "kernel")])

data_wind <- data[(data$warming_scenario %in% 1.0 &data$crop %in% "hay_alfalfa" & data$period %in% c(105:110) & data$climate_base %in% "1991_2020"),]
data_wind <- dplyr::inner_join(data_wind,optimal_gw[1,c("p", "theta", "longlat", "DistName", "kernel")])

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
optimal_gw <- as.data.frame(readRDS("output/optimal_gw.rds"))
data <- as.data.frame(readRDS("output/summary/summary_impact_yield.rds"))
data$est <- data$Estimate
data$se  <- data$Estimate_sd
data <- data[!data$county_code %in% 0,]
data <- data[(data$crop %in% "hay_alfalfa" & data$period %in% 105 & data$climate_base %in% "1991_2020"),]
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
optimal_gw <- as.data.frame(readRDS("output/optimal_gw.rds"))

yield <- as.data.frame(readRDS("output/summary/summary_impact_yield.rds"))
yield$impact_yield <- yield$Estimate
yield$fip<-as.character(paste0(stringr::str_pad(as.numeric(as.character(yield$state_code)), 2, pad = "0"),
                               stringr::str_pad(as.numeric(as.character(yield$county_code)), 3, pad = "0")))
yield <- yield[!yield$state_code %in% 0,]
yield <- yield[yield$warming_scenario %in% 1.0,]
yield <- yield[(yield$crop %in% "hay_alfalfa" & yield$period %in% 105 & yield$climate_base %in% "1991_2020"),]
yield <- dplyr::inner_join(yield,optimal_gw[1,c("p", "theta", "longlat", "DistName", "kernel")])
yield <- yield[c("state_code","county_code","fip","impact_yield")]
data <- yield

avail <- as.data.frame(readRDS("output/summary/summary_availability.rds"))
avail$level_avail <- avail$Estimate
avail <- avail[!avail$fip %in% "0000",]
avail <- avail[avail$warming_scenario %in% c(0.0),]
avail <- avail[(avail$crop %in% "hay_alfalfa" & avail$period %in% 105 & avail$climate_base %in% "1991_2020"),]
avail <- dplyr::inner_join(avail,optimal_gw[1,c("p", "theta", "longlat", "DistName", "kernel")])
avail <- avail[c("fip","level_avail")]
data <- dplyr::inner_join(data,avail,by="fip")

res <- as.data.frame(readRDS("output/summary/summary_associations.rds"))
res <- res[!res$fip %in% "00000",]
res <- res[(res$crop %in% "hay_alfalfa" & res$period %in% 105 & res$climate_base %in% "1991_2020"),
           c("p", "theta", "longlat", "DistName", "kernel","fip","name","est")]
res <- res |> tidyr::spread(name, est)
res <- dplyr::inner_join(res,optimal_gw[1,c("p", "theta", "longlat", "DistName", "kernel")])
res <- res[c("fip","avail00","prod00","prod00_LM")]
names(res) <- c("fip","b_avail00","b_prod00","b_prod00_LM")
data <- dplyr::inner_join(data,res,by="fip")

avail <- as.data.frame(readRDS("output/summary/summary_availability.rds"))
avail$impact_avail <- avail$Estimate
avail <- avail[!avail$fip %in% "0000",]
avail <- avail[avail$warming_scenario %in% 1.0,]
avail <- avail[(avail$crop %in% "hay_alfalfa" & avail$period %in% 105 & avail$climate_base %in% "1991_2020"),]
avail <- dplyr::inner_join(avail,optimal_gw[1,c("p", "theta", "longlat", "DistName", "kernel")])
avail <- avail[c("fip","impact_avail")]
data <- dplyr::inner_join(data,avail,by="fip")

catle <- as.data.frame(readRDS("output/summary/summary_cattle.rds"))
catle$impact_cattle <- catle$Estimate
catle <- catle[!catle$fip %in% "0000",]
catle <- catle[catle$warming_scenario %in% 1.0,]
catle <- catle[(catle$crop %in% "hay_alfalfa" & catle$period %in% 105 & catle$climate_base %in% "1991_2020"),]
catle <- dplyr::inner_join(catle,optimal_gw[1,c("p", "theta", "longlat", "DistName", "kernel")])
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
#-------------------------------
#-------------------------------
#-------------------------------
#-------------------------------
#-------------------------------
# Associations               ####
rm(list= ls()[!(ls() %in% c(Keep.List))])
optimal_gw <- as.data.frame(readRDS("output/optimal_gw.rds"))
data <- as.data.frame(readRDS("output/summary/summary_associations.rds"))
data <- data[!data$fip %in% "00000",]
data <- data[(data$crop %in% "hay_alfalfa" & data$period %in% 105 & data$climate_base %in% "1991_2020"),]
data <- dplyr::inner_join(data,optimal_gw[1,c("p", "theta", "longlat", "DistName", "kernel")])
#data <- data[data$name %in% c("prod00","prod00_LM"),]



merged_spat_vector <- terra::merge(Counties, data, by="fip")
sf_object <- sf::st_as_sf(merged_spat_vector)

cutlist <- unique(c(min(sf_object$est,na.rm=T),max(sf_object$est,na.rm=T),
                    quantile(round(unique(sf_object$est[! sf_object$est %in% c(NA,Inf,-Inf,NaN) | sf_object$est <0.0001]),3), probs = seq(0.056,0.32,0.001))))

avail00 <- sf_object$est[sf_object$name %in% "avail00"]
prod00 <- sf_object$est[sf_object$name %in% "prod00"]
prod00_LM <- sf_object$est[sf_object$name %in% "prod00_LM"]

cutlist <- unique(c(min(sf_object$est,na.rm=T),max(sf_object$est,na.rm=T),
                           quantile(unique(avail00), probs = seq(0.3,1,0.3)),
                           quantile(unique(prod00), probs = seq(0.3,1,0.3)),
                           quantile(unique(prod00_LM), probs = seq(0.3,1,0.3))))

cutlist <- c(0,unique(round(prod00,4)),unique(round(avail00,3)),unique(round(prod00_LM,3)))
sf_object$Value <- cut(sf_object$est,cutlist)
table(sf_object$Value)

sf_object <- sf_object[!is.na(sf_object$Value), ]


Fig06 <- ggplot() +
  geom_sf(data = sf::st_as_sf(States), colour = "black", fill = "darkred",size = 0.2) +
  geom_sf(data = sf_object,aes(fill = Value), colour = NA,size = 0.2) +
  geom_sf(data = sf::st_as_sf(States), colour = "black", fill = NA,size = 0.2) +
  scale_fill_manual(drop=FALSE, values=c(colorRampPalette(c("yellow","darkgreen"))(1+length(unique(as.character(sf_object$Value))))), na.value="#EEEEEE",
                    name="1,000 tons") +
  labs(title= "", x = "", y = "",fill ="", fill='',caption = "") +
  scale_y_continuous(sec.axis = sec_axis(~ . , name = "Kernel\n", breaks = NULL, labels = NULL)) +
  scale_x_continuous(sec.axis = sec_axis(~ . , name = "Distance metric\n",breaks = NULL, labels = NULL)) +
  guides(fill = guide_legend(nrow=2,override.aes = list(size=1))) +
  facet_wrap(~name,ncol=1) +
  ers_theme() +
  theme_bw() +
  theme(panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        axis.ticks.x = element_blank(),
        axis.ticks.y = element_blank(),
        plot.title=element_text(size=10),
        legend.position="bottom",
        legend.background = element_blank(),
        #legend.position=c(0.17,0.18),
        legend.key.size = unit(0.2,"cm"),
        legend.text = element_text(size=7),
        legend.title = element_text(size=7),
        axis.title.y = element_text(size=8),
        axis.title.x = element_text(size=8),
        axis.text.x  = element_blank(), #
        axis.text.y  = element_blank(),
        strip.text = element_text(size = 8),
        strip.background = element_blank())+coord_sf()

ggsave("output/exhibits/associations.png", Fig06, dpi = 600,width = 6.5, height = 8)

#-------------------------------
# Alfalfa Availability       ####
rm(list= ls()[!(ls() %in% c(Keep.List))])
optimal_gw <- as.data.frame(readRDS("output/optimal_gw.rds"))
avail <- as.data.frame(readRDS("output/summary/summary_availability.rds"))
avail$est <- avail$Estimate
avail <- avail[!avail$fip %in% "0",]
avail <- avail[avail$warming_scenario %in% c(0.0),]
avail <- avail[c("p","theta","longlat","DistName","kernel","crop","period","climate_base","fip","est")]
avail$outcome <- "(a) Alfalfa availability"

data <- avail

data <- dplyr::full_join(data,data.frame(optimal_gw[1,c("p", "theta", "longlat", "DistName", "kernel")],mainsest=T))

merged_spat_vector <- terra::merge(Counties, data, by="fip")
sf_object <- sf::st_as_sf(merged_spat_vector)

cutlist <- unique(c(min(sf_object$est,na.rm=T),max(sf_object$est,na.rm=T),
                    quantile(unique(sf_object$est[! sf_object$est %in% c(NA,Inf,-Inf,NaN) | sf_object$est <0.0001]), probs = seq(0.1,1,0.2))))

sf_object$Value <- cut(sf_object$est,cutlist)
table(sf_object$Value)

sf_object <- sf_object[!is.na(sf_object$Value), ]

sf_object$DistName <- as.numeric(as.character(factor(sf_object$DistName,
                                                     levels = c("Euclidean distance metric","Manhattan distance metric",
                                                                "Coordinate system is rotated by an angle 0.8 in radian"),
                                                     1:3)))

sf_object$DistName <- factor(sf_object$DistName,
                             levels = 1:3,
                             c("Euclidean","Manhattan",
                               "Coordinate system\nrotated by 0.8 radian"))


sf_object$kernel <- as.numeric(as.character(factor(sf_object$kernel,
                                                   levels = c("boxcar","bisquare","tricube","gaussian","exponential"),
                                                   1:5)))

sf_object$kernel <- factor(sf_object$kernel,
                           levels = 1:5,c("Boxcar","Bisquare","Tricube","Gaussian","Exponential"))

unique(as.data.frame(sf_object)[c("kernel")])


Fig06 <- ggplot() +
  geom_sf(data = sf::st_as_sf(States), colour = "black", fill = "darkred",size = 0.2) +
  geom_sf(data = sf_object,aes(fill = Value), colour = NA,size = 0.2) +
  geom_sf(data = sf::st_as_sf(States), colour = "black", fill = NA,size = 0.2) +
  scale_fill_manual(drop=FALSE, values=c(colorRampPalette(c("yellow","darkgreen"))(length(unique(as.character(sf_object$Value))))), na.value="#EEEEEE",
                    name="1,000 tons") +
  labs(title= "", x = "", y = "",fill ="", fill='',caption = "") +
  scale_y_continuous(sec.axis = sec_axis(~ . , name = "Kernel\n", breaks = NULL, labels = NULL)) +
  scale_x_continuous(sec.axis = sec_axis(~ . , name = "Distance metric\n",breaks = NULL, labels = NULL)) +
  guides(fill = guide_legend(nrow=2,override.aes = list(size=1))) +
  facet_grid(kernel~DistName) +
  ers_theme() +
  theme_bw() +
  theme(panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        axis.ticks.x = element_blank(),
        axis.ticks.y = element_blank(),
        plot.title=element_text(size=10),
        legend.position="bottom",
        legend.background = element_blank(),
        #legend.position=c(0.17,0.18),
        legend.key.size = unit(0.2,"cm"),
        legend.text = element_text(size=7),
        legend.title = element_text(size=7),
        axis.title.y = element_text(size=8),
        axis.title.x = element_text(size=8),
        axis.text.x  = element_blank(), #
        axis.text.y  = element_blank(),
        strip.text = element_text(size = 8),
        strip.background = element_blank())+coord_sf()

ggsave("output/exhibits/alfalfa_availability_all.png", Fig06, dpi = 600,width = 6.5, height = 8)


Fig_avail <- ggplot() +
  #geom_sf(data = sf::st_as_sf(States), colour = "black", fill = "darkred",size = 0.2) +
  geom_sf(data = sf_object[(sf_object$mainsest %in% TRUE),],aes(fill = Value), colour = NA,size = 0.2) +
  geom_sf(data = sf::st_as_sf(States), colour = "black", fill = NA,size = 0.2) +
  scale_fill_manual(drop=FALSE, values=c(colorRampPalette(c("yellow","darkgreen"))(length(unique(as.character(sf_object$Value))))), na.value="#EEEEEE",
                    name="1,000 tons") +
  labs(title= "(a) Alfalfa availability", x = "", y = "",fill ="", fill='',caption = "") +
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

ggsave("output/exhibits/alfalfa_availability_main.png", Fig_avail, dpi = 600,width = 6.5, height = 5)


#-------------------------------
# cattle                     ####
corr_cat <- as.data.frame(readRDS("output/summary/summary_cattle.rds"))
corr_cat$est <- corr_cat$corr_et
corr_cat <- corr_cat[!corr_cat$fip %in% "0",]
corr_cat <- corr_cat[c("p","theta","longlat","DistName","kernel","crop","period","climate_base","fip","est")]
corr_cat$outcome <- "(b) corr_cat"
optimal_gw <- as.data.frame(readRDS("Results/optimal_gw.rds"))
corr_cat <- dplyr::inner_join(corr_cat,data.frame(optimal_gw[1,c("p", "theta", "longlat", "DistName", "kernel")],mainsest=T))

merged_spat_vector <- terra::merge(Counties, corr_cat, by="fip")
sf_object <- sf::st_as_sf(merged_spat_vector)

cutlist <- unique(c(min(sf_object$est,na.rm=T),max(sf_object$est,na.rm=T),
                    quantile(unique(sf_object$est[! sf_object$est %in% c(NA,Inf,-Inf,NaN) | sf_object$est <0.0001]), probs = seq(0.1,1,0.2))))

sf_object$Value <- cut(sf_object$est,cutlist)
table(sf_object$Value)

sf_object <- sf_object[!is.na(sf_object$Value), ]

sf_object$DistName <- as.numeric(as.character(factor(sf_object$DistName,
                                                     levels = c("Euclidean distance metric","Manhattan distance metric",
                                                                "Coordinate system is rotated by an angle 0.8 in radian"),
                                                     1:3)))

sf_object$DistName <- factor(sf_object$DistName,
                             levels = 1:3,
                             c("Euclidean","Manhattan",
                               "Coordinate system\nrotated by 0.8 radian"))


sf_object$kernel <- as.numeric(as.character(factor(sf_object$kernel,
                                                   levels = c("boxcar","bisquare","tricube","gaussian","exponential"),
                                                   1:5)))

sf_object$kernel <- factor(sf_object$kernel,
                           levels = 1:5,c("Boxcar","Bisquare","Tricube","Gaussian","Exponential"))

Fig_corr <- ggplot() +
  #geom_sf(data = sf::st_as_sf(States), colour = "black", fill = "darkred",size = 0.2) +
  geom_sf(data = sf_object[(sf_object$mainsest %in% TRUE),],aes(fill = Value), colour = NA,size = 0.2) +
  geom_sf(data = sf::st_as_sf(States), colour = "black", fill = NA,size = 0.2) +
  scale_fill_manual(drop=FALSE, values=c(colorRampPalette(c("#FFEBCD","#800000"))(length(unique(as.character(sf_object$Value))))), na.value="#EEEEEE",
                    name="Percentage") +
  labs(title= "(b) Responsiveness of cattle inventory to alfalfa availability", x = "", y = "",fill ="", fill='',caption = "") +
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

ggsave("output/exhibits/corr_catle_avail_main.png", Fig_corr, dpi = 600,width = 6.5, height = 5)

marg <- c(-0.01,0.05,-0.01,0.05)

Fig <- cowplot::plot_grid(
  Fig_avail + theme(plot.margin=unit(marg, "cm")) ,
  Fig_corr + theme(plot.margin=unit(marg, "cm")) ,
  ncol=1, align="v",rel_heights=c(1,1),
  greedy=F)

ggsave("output/exhibits/availability_cattle.png", Fig, dpi = 600,width = 6.5, height = 9)

#-------------------------------
# Regional estimates         ####
rm(list= ls()[!(ls() %in% c(Keep.List))]);gc()
optimal_gw <- as.data.frame(readRDS("Results/optimal_gw.rds"))

yield <- as.data.frame(readRDS("Results/summary_impact_yield.rds"))
yield$impact_yield <- yield$Estimate
yield$fip<-as.character(paste0(stringr::str_pad(as.numeric(as.character(yield$state_code)), 2, pad = "0"),
                               stringr::str_pad(as.numeric(as.character(yield$county_code)), 3, pad = "0")))
yield <- yield[!yield$state_code %in% 0,]
yield <- yield[yield$warming_scenario %in% 1.0,]
yield <- yield[(yield$crop %in% "hay_alfalfa" & yield$period %in% 105 & yield$climate_base %in% "1991_2020"),]
yield <- dplyr::inner_join(yield,optimal_gw[1,c("p", "theta", "longlat", "DistName", "kernel")])
yield <- yield[c("state_code","county_code","fip","impact_yield")]

avail <- as.data.frame(readRDS("Results/summary_impact_avail.rds"))
avail$level_avail <- avail$Estimate
avail <- avail[!avail$fip %in% "0",]
avail <- avail[avail$warming_scenario %in% c(0.0),]
avail <- avail[(avail$crop %in% "hay_alfalfa" & avail$period %in% 105 & avail$climate_base %in% "1991_2020"),]
avail <- dplyr::inner_join(avail,optimal_gw[1,c("p", "theta", "longlat", "DistName", "kernel")])
avail <- avail[c("fip","level_avail")]
data <- dplyr::inner_join(data,avail,by="fip")

res <- as.data.frame(readRDS("Results/summary_corr_catle_avail.rds"))
res$corr_cat <- res$corr_et
res <- res[!res$fip %in% "0",]
res <- res[(res$crop %in% "hay_alfalfa" & res$period %in% 105 & res$climate_base %in% "1991_2020"),]
res <- dplyr::inner_join(res,optimal_gw[1,c("p", "theta", "longlat", "DistName", "kernel")])
res <- res[c("fip","corr_cat")]
data <- dplyr::inner_join(data,res,by="fip")

avail <- as.data.frame(readRDS("Results/summary_impact_avail.rds"))
avail$impact_avail <- avail$Estimate
avail <- avail[!avail$fip %in% "0",]
avail <- avail[avail$warming_scenario %in% 1.0,]
avail <- avail[(avail$crop %in% "hay_alfalfa" & avail$period %in% 105 & avail$climate_base %in% "1991_2020"),]
avail <- dplyr::inner_join(avail,optimal_gw[1,c("p", "theta", "longlat", "DistName", "kernel")])
avail <- avail[c("fip","impact_avail")]
data <- dplyr::inner_join(yield,avail,by="fip")

catle <- as.data.frame(readRDS("Results/summary_impact_catle.rds"))
catle$impact_cattle <- catle$Estimate
catle <- catle[!catle$fip %in% "0",]
catle <- catle[catle$warming_scenario %in% 1.0,]
catle <- catle[(catle$crop %in% "hay_alfalfa" & catle$period %in% 105 & catle$climate_base %in% "1991_2020"),]
catle <- dplyr::inner_join(catle,optimal_gw[1,c("p", "theta", "longlat", "DistName", "kernel")])
catle <- catle[c("fip","impact_cattle")]
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

write.csv(data,"output/exhibits/figure_data/regional_estimates.csv")
#-------------------------------
# Impacts Spatial            ####
rm(list= ls()[!(ls() %in% c(Keep.List))])
optimal_gw <- as.data.frame(readRDS("Results/optimal_gw.rds"))
yield <- as.data.frame(readRDS("Results/summary_impact_yield.rds"))
yield$est <- yield$Estimate
yield$se  <- yield$Estimate_sd
yield <- yield[!yield$county_code %in% 0,]
yield <- yield[c("p","theta","longlat","DistName","kernel","crop","period","climate_base","warming_scenario","fip","est","se")]
yield$outcome <- "(a) Alfalfa yield"

avail <- as.data.frame(readRDS("Results/summary_impact_avail.rds"))
avail$est <- avail$Estimate
avail$se <- avail$Estimate_sd
avail <- avail[!avail$fip %in% "0",]
avail <- avail[!avail$warming_scenario %in% c(0.0),]
avail <- avail[c("p","theta","longlat","DistName","kernel","crop","period","climate_base","warming_scenario","fip","est","se")]
avail$outcome <- "(b) Alfalfa availability"

cattle <- as.data.frame(readRDS("Results/summary_impact_catle.rds"))
cattle$est <- cattle$Estimate
cattle$se <- cattle$Estimate_sd
cattle <- cattle[!cattle$fip %in% "0",]
cattle <- cattle[!cattle$warming_scenario %in% c(0.0),]
cattle <- cattle[c("p","theta","longlat","DistName","kernel","crop","period","climate_base","warming_scenario","fip","est","se")]
cattle$outcome <- "(c) Cattle inventory"

data <- rbind(avail,yield,cattle)

data <- data[(data$crop %in% "hay_alfalfa" & data$period %in% 105 & data$climate_base %in% "1991_2020"),]

data$SimCat<-factor(data$warming_scenario,levels=unique(data$warming_scenario),
                    labels = paste0("+",format(unique(data$warming_scenario), nsmall = 1)," °C"))

USMUR  <- terra::rast(paste0(Dr.USMUR,"USAMU_90m.tif"))
States <- terra::vect(paste0(Dr.POLYG,"USA_States.shp"))
States <- terra::project(States, terra::crs(USMUR))
States <- terra::crop(States, terra::ext(USMUR))

Counties <- terra::vect(paste0(Dr.POLYG,"USA_Counties.shp"))
Counties <- terra::project(Counties, terra::crs(USMUR))
Counties <- terra::crop(Counties, terra::ext(USMUR))
Counties$fip<-as.character(paste0(stringr::str_pad(as.numeric(as.character(Counties$STATEFP)), 2, pad = "0"),
                                  stringr::str_pad(as.numeric(as.character(Counties$COUNTYFP)), 3, pad = "0")))

# state_data <- dplyr::inner_join(state_data,stnames[c("state_code","State.Name","State.Abbreviation" )],by="state_code")
# write.csv(state_data,"output/exhibits/figure_data/impacts_yield_spatial.csv")


fig_fxn <- function(outcome){
  #outcome <- "(a) Alfalfa yield"

  merged_spat_vector <- terra::merge(Counties, data[data$outcome %in% outcome,], by="fip")
  sf_object <- sf::st_as_sf(merged_spat_vector)

  cutlist <- unique(c(min(sf_object$est,na.rm=T),0,max(sf_object$est,na.rm=T),
                      quantile(unique(sf_object$est[! sf_object$est %in% c(NA,Inf,-Inf,NaN) | sf_object$est <0.0001]), probs = seq(0.1,1,0.10))))

  sf_object$Value <- cut(sf_object$est,cutlist)
  table(sf_object$Value)

  sf_object <- sf_object[!is.na(sf_object$Value), ]
  table(sf_object$Value,sf_object$SimCat)
  table(sf_object$SimCat)
  Fig <- ggplot() +
    geom_sf(data = sf::st_as_sf(States), colour = "black", fill = "darkred",size = 0.1) +
    geom_sf(data = sf_object,aes(fill = Value), colour = NA,size = 0.2) +
    geom_sf(data = sf::st_as_sf(States), colour = "black", fill = NA,size = 0.1) +
    scale_fill_manual(drop=FALSE, values=c(colorRampPalette(c("red","yellow"))(length(unique(as.character(sf_object$Value)))-2),
                                           "green","darkgreen"), na.value="#EEEEEE",
                      name="") +
    labs(title= paste0(outcome," impact in percentage"), x = "", y = "",fill ="", fill='',caption = "") +
    guides(fill = guide_legend(nrow=1)) +
    facet_wrap(vars(SimCat),nrow=1) +
    ers_theme() +
    theme_bw() +
    theme(panel.grid.major = element_blank(),
          panel.grid.minor = element_blank(),
          axis.ticks.x = element_blank(),
          axis.ticks.y = element_blank(),
          plot.title=element_text(size=8),
          legend.position="top",
          legend.background = element_blank(),
          #legend.position=c(0.17,0.18),
          legend.key.size = unit(0.2,"cm"),
          legend.text=element_text(size=4),
          legend.title=element_blank(),
          axis.title.y = element_blank(),
          axis.title.x = element_blank(),
          axis.text.x = element_blank(), #
          axis.text.y = element_blank(),
          strip.text = element_blank(),
          strip.background = element_blank())+coord_sf()

  ggsave(paste0("output/exhibits/spatial_impart_",gsub(" ","_",gsub("[(]a[)] ","",gsub("[(]b[)] ","",gsub("[(]c[)] ","",outcome)))),".png"),
         Fig, dpi = 600,width = 7, height = 2)
  gc()
  return(Fig)
}

Figa <- fig_fxn("(a) Alfalfa yield")
Figb <- fig_fxn("(b) Alfalfa availability")
Figc <- fig_fxn("(c) Cattle inventory")

marg <- c(-0.01,0.05,-0.01,0.05)

Fig <- cowplot::plot_grid(
  Figa + theme(plot.margin=unit(marg, "cm")) ,
  Figb + theme(plot.margin=unit(marg, "cm")) ,
  Figc + theme(plot.margin=unit(marg, "cm")),
  ncol=1, align="v",rel_heights=c(1,1,1),
  greedy=F)

stripT <- gtable_filter(ggplot_gtable(ggplot_build(
  Figa + theme(strip.background = element_blank(),strip.text = element_text(size=7)))), "strip-t")

Figx <- grid.arrange(stripT,Fig,heights=c(0.05,1),ncol = 1)

ggsave(paste0("output/exhibits/spatial_impart.png"), Figx, dpi = 600,width = 6.5, height = 4.3)

#-------------------------------
# Impacts mean               ####
rm(list= ls()[!(ls() %in% c(Keep.List))])
optimal_gw <- as.data.frame(readRDS("Results/optimal_gw.rds"))
yield <- as.data.frame(readRDS("Results/summary_impact_yield.rds"))
yield$est <- yield$Estimate
yield$se  <- yield$Estimate_sd
yield <- yield[yield$county_code %in% 0,]
yield <- yield[c("p","theta","longlat","DistName","kernel","crop","period","climate_base","warming_scenario","est","se")]
yield$outcome <- "(a) Alfalfa yield"

# avail <- as.data.frame(readRDS("Results/summary_impact_avail.rds"))
# avail$est <- avail$Estimate
# avail$se <- avail$Estimate_sd
# avail <- avail[avail$fip %in% "0",]
# avail <- avail[!avail$warming_scenario %in% c(0.0),]
# avail <- avail[c("p","theta","longlat","DistName","kernel","crop","period","climate_base","warming_scenario","est","se")]
# avail$outcome <- "(b) Alfalfa availability"
#
# cattle <- as.data.frame(readRDS("Results/summary_impact_catle.rds"))
# cattle$est <- cattle$Estimate
# cattle$se <- cattle$Estimate_sd
# cattle <- cattle[cattle$fip %in% "0",]
# cattle <- cattle[!cattle$warming_scenario %in% c(0.0),]
# cattle <- cattle[c("p","theta","longlat","DistName","kernel","crop","period","climate_base","warming_scenario","est","se")]
# cattle$outcome <- "(c) Cattle inventory"

# data <- rbind(avail,yield,cattle)
data <- yield
data_main <- data[(data$crop %in% "hay_alfalfa" & data$period %in% 105 & data$climate_base %in% "1991_2020"),]
data_main <- dplyr::inner_join(data_main,optimal_gw[1,c("p", "theta", "longlat", "DistName", "kernel")])

data_crop <- data[(data$warming_scenario %in% 1.0 & data$crop %in% c("hay_other","hay_alfalfa") & data$period %in% 105 & data$climate_base %in% "1991_2020"),]
data_crop <- dplyr::inner_join(data_crop,optimal_gw[1,c("p", "theta", "longlat", "DistName", "kernel")])

data_year <- data[(data$warming_scenario %in% 1.0 &data$crop %in% "hay_alfalfa" & data$period %in% 105 & data$climate_base %in% c("1991_2020","1981_2010","1971_2000","1961_1990")),]
data_year <- dplyr::inner_join(data_year,optimal_gw[1,c("p", "theta", "longlat", "DistName", "kernel")])

data_wind <- data[(data$warming_scenario %in% 1.0 &data$crop %in% "hay_alfalfa" & data$period %in% c(105:110) & data$climate_base %in% "1991_2020"),]
data_wind <- dplyr::inner_join(data_wind,optimal_gw[1,c("p", "theta", "longlat", "DistName", "kernel")])

data_dist <- data[(data$warming_scenario %in% 1.0 &data$crop %in% "hay_alfalfa" & data$period %in% 105 & data$climate_base %in% c("1991_2020") &
                     data$kernel %in% optimal_gw[1,"kernel"]),]

data_kenl <- data[(data$warming_scenario %in% 1.0 &data$crop %in% "hay_alfalfa" & data$period %in% 105 & data$climate_base %in% c("1991_2020") &
                     data$DistName %in% optimal_gw[1,"DistName"]),]

data_main$type <- "Impact by warming scenario"
data_year$type <- "1+°C warming by climate baseline"
data_wind$type <- "1+°C warming by season window"
data_crop$type <- "1+°C warming by hay type"
data_kenl$type <- "1+°C warming by kernel"
data_dist$type <- "1+°C warming by distance metric"

data_main <- data_main[order(data_main$warming_scenario),]
data_kenl <- data_kenl[order(data_kenl$est),]
data_dist <- data_dist[order(data_dist$est),]
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

data_kenl$x_name <- as.character(factor(data_kenl$kernel,levels = c("exponential","gaussian","boxcar","bisquare","tricube"),
                                        labels = c("Exponential","Gaussian","Boxcar","Bisquare","Tricube")))
data_dist$x_name <- as.character(factor(data_dist$DistName,levels = c("Manhattan distance metric","Euclidean distance metric",
                                                                      "Coordinate system is rotated by an angle 0.8 in radian"),
                                        labels = c("Manhattan","Euclidean","Coordinate system")))

unique(data_dist$DistName)

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

data_text <- doBy::summaryBy(x~type,data=data,FUN=max,keep.names = T,na.rm=T)
data_text$outcome <- "(a) Alfalfa yield"

Fig04 <- ggplot(data,aes(x=x,y=est,group=1)) +
  geom_hline(yintercept = 0,size = 0.2,color = "black") +
  geom_errorbar(aes(ymin = est - se*1.96, ymax = est + se*1.96),color="purple") +
  geom_point(color="purple") +
  geom_text(data=data_text,aes(x=x+1.5,y=-13,label=type, hjust = 0), size = 2.5, col = "black",
            stat = "identity",check_overlap = TRUE) +
  scale_y_continuous(breaks=seq(-60,0,2.5),labels=format(round(seq(-60,0,2.5),2), nsmall = 2)) +
  #facet_wrap(~outcome,nrow=1,scales = "free_x") +
  scale_x_continuous(breaks=data$x,labels=data$x_name) +
  labs(title="", x="", y ="\nImpact (%)",caption = "") +
  myTheme +
  theme(axis.text.x = element_text(size=7,color="black"),
        axis.text.y = element_text(size=7,color="black"),
        strip.text = element_text(size = 9),
        legend.position="none")+ coord_flip()

write.csv(data,"output/exhibits/figure_data/impacts_mean.csv")
ggsave("output/exhibits/impacts_mean.png", Fig04, dpi = 600,width = 4, height = 6)
#ggsave("output/exhibits/impacts_mean.png", Fig04, dpi = 600,width = 7, height = 6)

#-------------------------------
