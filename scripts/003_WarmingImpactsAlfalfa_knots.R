
rm(list=ls(all=TRUE));library(future.apply);library(data.table);library(gmm);library(tidyverse);gc()

invisible(lapply(list.files("scripts/helpers", pattern = "[.]R$", full.names = TRUE), source))

study_environment <- readRDS("data/study_environment.rds")

#
options(future.globals.maxSize = 8000 * 1024^2)

ARRAY <- data.frame(crop=c("hay_alfalfa","hay_other","hay_all"))
ARRAY <- as.data.frame(
  data.table::rbindlist(
    lapply(c(0:12,101:112),function(period){return(data.frame(period=period,ARRAY))}), fill = TRUE))

#ARRAY <- ARRAY[!paste0("optimal_knots_",ARRAY$crop,"_period",ARRAY$period,".rds") %in% basename(list.files(study_environment$wd$knots)),]

if(!is.na(as.numeric(Sys.getenv("SLURM_ARRAY_TASK_ID")))){
  ARRAY$TASK <- rep(as.numeric(Sys.getenv("SLURM_ARRAY_TASK_MIN")):as.numeric(Sys.getenv("SLURM_ARRAY_TASK_MAX")), length=nrow(ARRAY))
  ARRAY <- ARRAY[ARRAY$TASK %in% as.numeric(Sys.getenv("SLURM_ARRAY_TASK_ID")),]
}

lapply(
  1:nrow(ARRAY),
  function(task){
    tryCatch({

      # task <- 58

      target_periods <- ARRAY$period[task]
      crop   <- ARRAY$crop[task]

      data <- build_hay_weather_panel(
        crop = "hay_alfalfa",
        target_periods = target_periods,
        prism_weather_directory = "data/prism_weather")

      varlist <- c("ppt","ppt2","ppt3","Trend","Trend2","DD1","DD2","DD3")

      Pair <- data.frame(c1=NA,c2=NA)
      for(i in 5:35){for(j in (i+1):40){Pair<-rbind(Pair,data.frame(c1=i,c2=j))}}
      Pair <- Pair[complete.cases(Pair),]
      Pair <- Pair[Pair$c2-Pair$c1>=3,]
      Pair <- Pair[Pair$c1 %in% 10:20,]
      Pair <- Pair[Pair$c2 %in% 29:35,]

      Optknots <- data.table::rbindlist(
        lapply(
          unique(Pair$c1),
          function(c1){
            # c1 <- 11
            Pairset <- Pair[Pair$c1 %in% c1,]

            Optknots <- data.table::rbindlist(
              lapply(
                c(1:nrow(Pairset)), #[120:130]
                function(pair){
                  tryCatch({
                    # pair<-1;tmin <- 10;tmax<-29
                    print(paste0((round(pair/nrow(Pairset),2))*100,"%"))
                    datax <- data

                    tmin <- Pairset$c1[pair]
                    tmax <- Pairset$c2[pair]

                    datax[,DD1 := get(paste0("dday00")) - get(paste0("dday",stringr::str_pad(tmin,pad="0",2)))]
                    datax[,DD2 := get(paste0("dday",stringr::str_pad(tmin,pad="0",2))) - get(paste0("dday",stringr::str_pad(tmax,pad="0",2)))]
                    datax[,DD3 := get(paste0("dday",stringr::str_pad(tmax,pad="0",2)))]

                    datax <- doBy::summaryBy(lny + ppt + ppt2 + ppt3 + Trend + Trend2 + DD1 + DD2 + DD3~fip+state_code+commodity_year,
                                             data=datax,FUN=mean,keep.names = T,na.rm=T)

                    for(st in unique(datax$state_code)){
                      datax[,(paste0("trend1",st)) := ifelse(state_code %in% st,Trend,0)]
                      datax[,(paste0("trend2",st)) := ifelse(state_code %in% st,Trend2,0)]
                    }

                    fit <- plm::plm(as.formula(
                      paste0("lny ~",paste0(c("ppt","ppt2","DD1" ,"DD2" ,"DD3",
                                              names(datax)[grepl("trend",names(datax))]),collapse = "+"))) ,
                      data=datax, index=c("fip", "commodity_year"), model="within")

                    res <- data.frame(Tmin=tmin,Tmax=tmax,R = summary(fit)$r.squared["rsq"],t(coef(fit)[c("DD1","DD2","DD3")]))

                    datax[,fold := sample(1:5, nrow(datax), replace = TRUE,prob=rep(0.2,5))]

                    datax <-  plm::pdata.frame(datax, index = c("fip", "commodity_year"), drop.index = TRUE)
                    res_cv <- as.data.frame(
                      data.table::rbindlist(
                        lapply(
                          1:5, # [120:130]
                          function(f){
                            tryCatch({
                              # f <- 1

                              data_train <- datax[!datax$fold %in% f, ]
                              data_test  <- datax[datax$fold %in% f, ]

                              rhs_vars <- c(
                                "ppt", "ppt2", "DD1", "DD2", "DD3",
                                names(data_train)[grepl("trend", names(data_train))]
                              )

                              form <- as.formula(
                                paste("lny ~", paste(rhs_vars, collapse = " + "))
                              )

                              fit_f <- plm::plm(
                                formula = form,
                                data    = data_train,
                                index   = c("fip", "commodity_year"),
                                model   = "within"
                              )

                              beta <- coef(fit_f)

                              X0 <- model.matrix(form, data_test)
                              X_test <- X0[, names(beta), drop = FALSE]
                              pred_slope <- as.numeric(X_test %*% beta)

                              fe <- plm::fixef(fit_f)
                              fip_test <- as.character(plm::index(data_test)[[1]])
                              alpha_test <- unname(fe[match(fip_test, names(fe))])

                              keep <- !is.na(alpha_test)

                              out <- data.frame(
                                f = f,
                                n = sum(keep),
                                e = sum((pred_slope[keep] + alpha_test[keep] - data_test$lny[keep])^2, na.rm = TRUE)
                              )

                              out

                              return(out)
                            }, error = function(e){return(NULL)})
                          }), fill = TRUE))

                    res_cv <- res_cv[!res_cv$e %in% c(NaN,Inf,-Inf,NA),]
                    res$cv_error <- weighted.mean(x=res_cv$e,w=res_cv$n)

                    function(){
                      res <- data.frame(Temp=1:45,DD1=(coef(fit)["DD1"]),DD2=(coef(fit)["DD2"]),DD3=(coef(fit)["DD3"]))
                      res$I1<-as.numeric(res$Temp>=0)
                      res$I2<-as.numeric(res$Temp>=tmin)
                      res$I3<-as.numeric(res$Temp>=tmax)
                      res$response<-(1-res$I2)*(res$DD1*(res$Temp-0)) +
                        (res$I2)*(1-res$I3)*(res$DD1*(tmin) + res$DD2*(res$Temp-tmin)) +
                        (res$I3)*(res$DD1*(tmin) + res$DD2*(tmax-tmin) + res$DD3*(res$Temp-tmax))
                      plot(res$Temp,res$response)
                    }

                    return(res)
                  }, error = function(e){return(NULL)})
                }), fill = TRUE)

            return(Optknots)
          }), fill = TRUE)

      # OptknotsP <- Optknots[complete.cases(Optknots),]
      # # OptknotsP <- OptknotsP[OptknotsP$DD1 >=0,]
      # # OptknotsP <- OptknotsP[( OptknotsP$DD1 <= OptknotsP$DD2),]
      # OptknotsP <- OptknotsP[round(OptknotsP$R,6) == max(round(OptknotsP$R,6)),]
      # OptknotsP <- OptknotsP[(OptknotsP$Tmax-OptknotsP$Tmin) == max(OptknotsP$Tmax-OptknotsP$Tmin),]

      Optknots[,crop := crop]
      Optknots[,target_periods := target_periods]

      saveRDS(Optknots,file=file.path(study_environment$wd$knots,paste0("optimal_knots_",ARRAY$crop[task],"_period",ARRAY$period[task],".rds")))

      NULL
    }, error = function(e){return(NULL)})

  })

function(){

  OptknotsP <- as.data.frame(
    data.table::rbindlist(
    lapply(
      list.files(study_environment$wd$knots,full.names = T), # [120:130]
      function(file){
        tryCatch({
          return(readRDS(file))
        }, error = function(e){return(NULL)})
      }), fill = TRUE))

  OptknotsP <- OptknotsP[complete.cases(OptknotsP),]
  OptknotsP <- OptknotsP[OptknotsP$DD1 >=0,]
  OptknotsP <- OptknotsP[OptknotsP$DD2 >0,]
  OptknotsP <- OptknotsP[OptknotsP$DD3 <0,]
  # OptknotsP <- OptknotsP[( OptknotsP$DD1 <= OptknotsP$DD2),]

  OptknotsP <- OptknotsP |> group_by(crop,target_periods) |> mutate(cv_error_min = max(cv_error,na.rm=T)) |> as.data.frame(.)
  OptknotsP <- OptknotsP[round(OptknotsP$cv_error,6) == round(OptknotsP$cv_error_min,6),]

  OptknotsP <- OptknotsP |> group_by(crop,target_periods) |> mutate(R_max = max(R,na.rm=T)) |> as.data.frame(.)
  OptknotsP <- OptknotsP[round(OptknotsP$R_max,6) == round(OptknotsP$R,6),]

  OptknotsP <- OptknotsP |> group_by(crop,target_periods) |> mutate(tdiff = max(Tmax-Tmin,na.rm=T)) |> as.data.frame(.)

  OptknotsP <- OptknotsP[(OptknotsP$Tmax-OptknotsP$Tmin) == OptknotsP$tdiff,]

  # OptknotsP <- OptknotsP |> group_by(crop) |> mutate(R_max = max(R,na.rm=T)) |> as.data.frame(.)
  #
  # OptknotsP <- OptknotsP[round(OptknotsP$R,6) == round(OptknotsP$R_max,6),]
  #
  # OptknotsP <- OptknotsP |> group_by(crop) |> mutate(tdiff = max(Tmax-Tmin,na.rm=T)) |> as.data.frame(.)
  #
  # OptknotsP <- OptknotsP[(OptknotsP$Tmax-OptknotsP$Tmin) == OptknotsP$tdiff,]

  saveRDS(OptknotsP,file=paste0("output/optimal_knots.rds"))

}
# unlink(list.files(getwd(),pattern =".out",full.names = T))
