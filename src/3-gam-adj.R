rm(list=ls())

source(here::here("0-config.R"))

d<-readRDS(paste0(dropboxDir, "Data/Cleaned/Audrie/pregnancy_stress_data.RDS"))

#Set list of adjustment variables
#Make vectors of adjustment variable names
Wvars<-c("sex","birthord", "momage","momheight","momedu", 
  "hfiacat", "Nlt18","Ncomp", "watmin", "walls", "floor", "roof", "HHwealth",
  "tr", "life_viol_any_t3", "viol_any_preg")

Wvars[!(Wvars %in% colnames(d))]



#Add in time varying covariates:
Wvars2 <- c(Wvars, c("month_blood_t0", "month_ut2", "ageday_ut2")) 
Wvars3 <- c(Wvars, c("month_blood_t0", "ageday_t3_vital", "ageday_t3_oragene", "ageday_t3_salimetrics", "month_vt3", "month_ot3", "month_lt3"))
Wvars23 <- c(Wvars, c("month_blood_t0", "month_ut2", "ageday_ut2", "ageday_t3_vital", "ageday_t3_oragene", "ageday_t3_salimetrics", "month_vt3", "month_ot3", "month_lt3"))

pick_covariates <- function(j){
  # j is outcome as string
  # choose correct adjustment set based on outcome
  if(grepl("t2", j)){Wset = Wvars2}
  else if(grepl("t3", j)){Wset = Wvars3}
  else{Wset = Wvars23}
  return(Wset)
}


#Loop over exposure-outcome pairs

##Hypothesis 1
#Maternal cortisol levels during pregnancy are associated with child stress biomarkers at Year 1 and Year 2. 
Xvars <- c("ln_preg_cort")            
Yvars <- c("t2_f2_8ip",  "t2_f2_23d", "t2_f2_VI", "t2_f2_12i", "t2_iso_pca", 
           "t3_cort_slope", "t3_cort_z01", "t3_cort_z03",
           "t3_saa_slope", "t3_saa_z01", "t3_saa_z02",
           "t3_map", "t3_hr_mean",
           "t3_gcr_mean", "t3_gcr_cpg12") 

#Fit models
H1_adj_models <- NULL
for(i in Xvars){
  for(j in Yvars){
    print(i)
    print(j)
    Wset<-pick_covariates(j)
    res_adj <- fit_RE_gam(d=d, X=i, Y=j,  W=Wset)
    res <- data.frame(X=i, Y=j, fit=I(list(res_adj$fit)), dat=I(list(res_adj$dat)))
    H1_adj_models <- bind_rows(H1_adj_models, res)
  }
}



#Get primary contrasts
H1_adj_res <- NULL
for(i in 1:nrow(H1_adj_models)){
  res <- data.frame(X=H1_adj_models$X[i], Y=H1_adj_models$Y[i])
  if(grepl("_def", H1_adj_models$X[i])){
    preds <- predict_gam_diff(fit=H1_adj_models$fit[i][[1]], d=H1_adj_models$dat[i][[1]], quantile_diff=c(0.25,0.75), Xvar=res$X, Yvar=res$Y, binary=T)
  }else{
    preds <- predict_gam_diff(fit=H1_adj_models$fit[i][[1]], d=H1_adj_models$dat[i][[1]], quantile_diff=c(0.25,0.75), Xvar=res$X, Yvar=res$Y)
  }
  H1_adj_res <-  bind_rows(H1_adj_res , preds$res)
}

#Make list of plots
H1_adj_plot_list <- NULL
H1_adj_plot_data <- NULL
for(i in 1:nrow(H1_adj_models)){
  res <- data.frame(X=H1_adj_models$X[i], Y=H1_adj_models$Y[i])
  simul_plot <- gam_simul_CI(H1_adj_models$fit[i][[1]], H1_adj_models$dat[i][[1]], xlab=res$X, ylab=res$Y, title="")
  H1_adj_plot_list[[i]] <-  simul_plot$p
  H1_adj_plot_data <-  rbind(H1_adj_plot_data, data.frame(Xvar=res$X, Yvar=res$Y, adj=0, simul_plot$pred %>% subset(., select = c(Y,X,id,fit,se.fit,uprP, lwrP,uprS,lwrS))))
}


#Save models
#saveRDS(H1_adj_models, paste0(dropboxDir,"results/pregnancy-stress-models/models/H1_adj_models.RDS"))

#Save results
saveRDS(H1_adj_res, here("results/adjusted/H1_adj_res.RDS"))


#Save plots
#saveRDS(H1_adj_plot_list, paste0(dropboxDir,"results/pregnancy-stress-models/figure-objects/H1_adj_splines.RDS"))

#Save plot data
saveRDS(H1_adj_plot_data, paste0(dropboxDir,"results/pregnancy-stress-models/figure-data/H1_adj_spline_data.RDS"))


## Hypothesis 2
# Maternal inflammation during pregnancy is associated with child stress biomarkers at Year 1 and Year 2. 
# change in telomere length

Xvars <- c("logCRP", "logAGP", "sumscore_t0_mom_Z")            
Yvars <- c("t2_f2_8ip",  "t2_f2_23d", "t2_f2_VI", "t2_f2_12i", "t2_iso_pca", 
           "t3_cort_slope",  "t3_cort_z01", "t3_cort_z03",
           "t3_saa_slope",  "t3_saa_z01", "t3_saa_z02",
           "t3_map", "t3_hr_mean",
           "t3_gcr_mean", "t3_gcr_cpg12") 

#Fit models
H2_adj_models <- NULL
for(i in Xvars){
  for(j in Yvars){
    print(i)
    print(j)
    Wset<-pick_covariates(j)
    res_adj <- fit_RE_gam(d=d, X=i, Y=j,  W=Wset)
    res <- data.frame(X=i, Y=j, fit=I(list(res_adj$fit)), dat=I(list(res_adj$dat)))
    H2_adj_models <- bind_rows(H2_adj_models, res)
  }
}

#Get primary contrasts
H2_adj_res <- NULL
for(i in 1:nrow(H2_adj_models)){
  res <- data.frame(X=H2_adj_models$X[i], Y=H2_adj_models$Y[i])
  preds <- predict_gam_diff(fit=H2_adj_models$fit[i][[1]], d=H2_adj_models$dat[i][[1]], quantile_diff=c(0.25,0.75), Xvar=res$X, Yvar=res$Y)
  H2_adj_res <-  bind_rows(H2_adj_res , preds$res)
}

#Make list of plots
H2_adj_plot_list <- NULL
H2_adj_plot_data <- NULL
for(i in 1:nrow(H2_adj_models)){
  res <- data.frame(X=H2_adj_models$X[i], Y=H2_adj_models$Y[i])
  simul_plot <- gam_simul_CI(H2_adj_models$fit[i][[1]], H2_adj_models$dat[i][[1]], xlab=res$X, ylab=res$Y, title="")
  H2_adj_plot_list[[i]] <-  simul_plot$p
  H2_adj_plot_data <-  rbind(H2_adj_plot_data, data.frame(Xvar=res$X, Yvar=res$Y, adj=0, simul_plot$pred%>% subset(., select = c(Y,X,id,fit,se.fit,uprP, lwrP,uprS,lwrS))))
}


#Save models
#saveRDS(H2_adj_models, paste0(dropboxDir,"results/pregnancy-stress-models/models/adj_H2_adj_models.RDS"))

#Save results
saveRDS(H2_adj_res, here("results/adjusted/H2_adj_res.RDS"))


#Save plots
#saveRDS(H2_adj_plot_list, paste0(dropboxDir,"results/pregnancy-stress-models/figure-objects/H2_adj_adj_splines.RDS"))

#Save plot data
saveRDS(H2_adj_plot_data, paste0(dropboxDir,"results/pregnancy-stress-models/figure-data/H2_adj_adj_spline_data.RDS"))



##Hypothesis 3
#Maternal nutrition during pregnancy is associated with child stress biomarkers at Year 1 and Year 2. 

Xvars <- c("vitD_nmol_per_L", "logFERR_inf", "logSTFR_inf", "logRBP_inf")            
Yvars <- c("t2_f2_8ip",  "t2_f2_23d", "t2_f2_VI", "t2_f2_12i", "t2_iso_pca", 
           "t3_cort_slope",  "t3_cort_z01", "t3_cort_z03",
           "t3_saa_slope",  "t3_saa_z01", "t3_saa_z02",
           "t3_map", "t3_hr_mean",
           "t3_gcr_mean", "t3_gcr_cpg12") 

#Fit models
H3_adj_models <- NULL
for(i in Xvars){
  for(j in Yvars){
    print(i)
    print(j)
    Wset<-pick_covariates(j)
    res_adj <- fit_RE_gam(d=d, X=i, Y=j,  W=Wset)
    res <- data.frame(X=i, Y=j, fit=I(list(res_adj$fit)), dat=I(list(res_adj$dat)))
    H3_adj_models <- bind_rows(H3_adj_models, res)
  }
}

#Get primary contrasts
H3_adj_res <- NULL
for(i in 1:nrow(H3_adj_models)){
  res <- data.frame(X=H3_adj_models$X[i], Y=H3_adj_models$Y[i])
  preds <- predict_gam_diff(fit=H3_adj_models$fit[i][[1]], d=H3_adj_models$dat[i][[1]], quantile_diff=c(0.25,0.75), Xvar=res$X, Yvar=res$Y)
  H3_adj_res <-  bind_rows(H3_adj_res , preds$res)
}

#Make list of plots
H3_adj_plot_list <- NULL
H3_adj_plot_data <- NULL
for(i in 1:nrow(H3_adj_models)){
  res <- data.frame(X=H3_adj_models$X[i], Y=H3_adj_models$Y[i])
  simul_plot <- gam_simul_CI(H3_adj_models$fit[i][[1]], H3_adj_models$dat[i][[1]], xlab=res$X, ylab=res$Y, title="")
  H3_adj_plot_list[[i]] <-  simul_plot$p
  H3_adj_plot_data <-  rbind(H3_adj_plot_data, data.frame(Xvar=res$X, Yvar=res$Y, adj=0, simul_plot$pred%>% subset(., select = c(Y,X,id,fit,se.fit,uprP, lwrP,uprS,lwrS))))
}


#Save models
#saveRDS(H3_models, paste0(dropboxDir,"results/pregnancy-stress-models/models/adj_H3_models.RDS"))

#Save results
saveRDS(H3_adj_res, here("results/adjusted/H3_adj_res.RDS"))


#Save plots
#saveRDS(H3_adj_plot_list, paste0(dropboxDir,"results/pregnancy-stress-models/figure-objects/H3_adj_splines.RDS"))

#Save plot data
saveRDS(H3_adj_plot_data, paste0(dropboxDir,"results/pregnancy-stress-models/figure-data/H3_adj_spline_data.RDS"))


##Hypothesis 4
#Maternal estriol during pregnancy is associated with child stress biomarkers at Year 1 and Year 2. 

Xvars <- c("ln_preg_estri")            
Yvars <- c("t2_f2_8ip",  "t2_f2_23d", "t2_f2_VI", "t2_f2_12i", "t2_iso_pca", 
           "t3_cort_slope",  "t3_cort_z01", "t3_cort_z03",
           "t3_saa_slope",  "t3_saa_z01", "t3_saa_z02",
           "t3_map", "t3_hr_mean",
           "t3_gcr_mean", "t3_gcr_cpg12")


#Fit models
H4_adj_models <- NULL
for(i in Xvars){
  for(j in Yvars){
    print(i)
    print(j)
    Wset<-pick_covariates(j)
    res_adj <- fit_RE_gam(d=d, X=i, Y=j,  W=Wset)
    res <- data.frame(X=i, Y=j, fit=I(list(res_adj$fit)), dat=I(list(res_adj$dat)))
    H4_adj_models <- bind_rows(H4_adj_models, res)
  }
}

#Get primary contrasts
H4_adj_res <- NULL
for(i in 1:nrow(H4_adj_models)){
  res <- data.frame(X=H4_adj_models$X[i], Y=H4_adj_models$Y[i])
  preds <- predict_gam_diff(fit=H4_adj_models$fit[i][[1]], d=H4_adj_models$dat[i][[1]], quantile_diff=c(0.25,0.75), Xvar=res$X, Yvar=res$Y)
  H4_adj_res <-  bind_rows(H4_adj_res , preds$res)
}

#Make list of plots
H4_adj_plot_list <- NULL
H4_adj_plot_data <- NULL
for(i in 1:nrow(H4_adj_models)){
  res <- data.frame(X=H4_adj_models$X[i], Y=H4_adj_models$Y[i])
  simul_plot <- gam_simul_CI(H4_adj_models$fit[i][[1]], H4_adj_models$dat[i][[1]], xlab=res$X, ylab=res$Y, title="")
  H4_adj_plot_list[[i]] <-  simul_plot$p
  H4_adj_plot_data <-  rbind(H4_adj_plot_data, data.frame(Xvar=res$X, Yvar=res$Y, adj=0, simul_plot$pred%>% subset(., select = c(Y,X,id,fit,se.fit,uprP, lwrP,uprS,lwrS))))
}


#Save models
#saveRDS(H4_adj_models, paste0(dropboxDir,"results/pregnancy-stress-models/models/H4_adj_models.RDS"))

#Save results
saveRDS(H4_adj_res, here("results/adjusted/H4_adj_res.RDS"))


#Save plots
#saveRDS(H4_adj_plot_list, paste0(dropboxDir,"results/pregnancy-stress-models/figure-objects/H4_adj_splines.RDS"))

#Save plot data
saveRDS(H4_adj_plot_data, paste0(dropboxDir,"results/pregnancy-stress-models/figure-data/H4_adj_spline_data.RDS"))

