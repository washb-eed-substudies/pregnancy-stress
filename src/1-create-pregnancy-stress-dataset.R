


rm(list=ls())

source(here::here("0-config.R"))


d<-read.csv(paste0(dropboxDir, "Data/Cleaned/Audrie/washb-bd-pregnancy-serum-micronutrient-immun-cortisol-covariates-child-stress.csv"))

#Clean stress outcomes

#---------------------------------------------------------------------------------------------
# transform outcome distributions
#---------------------------------------------------------------------------------------------
d <- d %>% 
  mutate(
    t3_saa_z01_raw=t3_saa_z01, 
    t3_saa_z02_raw=t3_saa_z02, 
    t3_cort_z01_raw=t3_cort_z01, 
    t3_cort_z03_raw=t3_cort_z03, 
    t2_f2_8ip_raw=t2_f2_8ip, 
    t2_f2_23d_raw=t2_f2_23d, 
    t2_f2_VI_raw=t2_f2_VI,
    t2_f2_12i_raw=t2_f2_12i, 
    t3_gcr_mean_raw=t3_gcr_mean, 
    t3_gcr_cpg12_raw=t3_gcr_cpg12,
    t3_saa_z01=log(t3_saa_z01), 
    t3_saa_z02=log(t3_saa_z02), 
    t3_cort_z01=log(t3_cort_z01), 
    t3_cort_z03=log(t3_cort_z03), 
    t2_f2_8ip=log(t2_f2_8ip), 
    t2_f2_23d=log(t2_f2_23d), 
    t2_f2_VI=log(t2_f2_VI),
    t2_f2_12i=log(t2_f2_12i), 
    t3_gcr_mean2=logit(t3_gcr_mean/100), 
    t3_gcr_cpg12=logit(t3_gcr_cpg12/100))
#Clean inf values
d <- do.call(data.frame,lapply(d, function(x) replace(x, is.infinite(x),NA)))

#---------------------------------------------------------------------------------------------
# calc combined iso variable
#---------------------------------------------------------------------------------------------
# Combined exposures: To determine overall oxidative stress, we will combine our four measures of urinary F2-
# isoprostanes: iPF(2a)-III, 2,3-dinor-iPF(2a)-III, iPF(2a)-VI, and 8,12-iso-iPF(2a)-VI that are
# consistent with prior oxidative stress operationalization into a single score.29 We will use
# the first principal component of a principal components analysis of the four measures of
# urinary F2-isoprostanes as the oxidative stress score if all measures are correlated with each other (P-value < 0.2), 
# otherwise we will analyze the urinary F2-isoprostanes separately.
#Get correlation of isoprostanes
iso <- d %>% select(c("t2_f2_8ip", "t2_f2_23d", "t2_f2_VI", "t2_f2_12i"))       
cor(iso, use="pairwise.complete.obs")
cor.test(iso[,1], iso[,2])$p.value < 0.2
cor.test(iso[,1], iso[,3])$p.value < 0.2
cor.test(iso[,1], iso[,4])$p.value < 0.2
cor.test(iso[,2], iso[,3])$p.value < 0.2
cor.test(iso[,2], iso[,4])$p.value < 0.2
cor.test(iso[,3], iso[,4])$p.value < 0.2
#isoprostanes are significantly correlated, so calculate 1st principal component
df <-  d %>% select(c("childid","t2_f2_8ip", "t2_f2_23d", "t2_f2_VI", "t2_f2_12i")) %>% as.data.frame()
dim(df)
df <- df[complete.cases(df),]
dim(df)
# #Select assets and seperate out ID
id<-subset(df, select=c("childid")) #drop subjectid
df<-df[,which(!(colnames(df) %in% c("childid")))]
##Computing the principal component using eigenvalue decomposition ##
princ.return <- princomp(df) 
## To get the first principal component in a variable ##
load <- loadings(princ.return)[,1]   
pr.cp <- as.matrix(df) %*% load  ## Matrix multiplication of the input data with the loading for the 1st PC gives us the 1st PC in matrix form. 
df$t2_iso_pca <- as.numeric(pr.cp) ## Gives us the 1st PC in numeric form in pr.
#examine combined score
hist(df$t2_f2_8ip)
hist(df$t2_f2_23d)
hist(df$t2_f2_VI)
hist(df$t2_f2_12i)
hist(df$t2_iso_pca)
#check direction between individual isoprostanes and the combined score
ggplot(df, aes(x=t2_f2_8ip, y=t2_iso_pca)) + geom_point() + geom_smooth()
ggplot(df, aes(x=t2_f2_23d, y=t2_iso_pca)) + geom_point() + geom_smooth()
ggplot(df, aes(x=t2_f2_VI, y=t2_iso_pca)) + geom_point() + geom_smooth()
ggplot(df, aes(x=t2_f2_12i, y=t2_iso_pca)) + geom_point() + geom_smooth()
#merge combined score back into main dataset
df.pca <- data.frame(childid=id, t2_iso_pca=df$t2_iso_pca)
d <- left_join(d, df.pca, by="childid")
hist(d$t2_iso_pca)

## add household wealth
d_hhwealth <- read.csv(paste0(dropboxDir, "Data/Cleaned/Audrie/hhwealth.csv"))
d <- left_join(d, d_hhwealth, by="dataid")
## convert time of day of pre-stressor measurement of cortisol and sAA into continuous variable
time_day <- d$t3_col_time_z01
time_split <- str_split(time_day, ":")
cont_time <- function(list_hr_min){
  # takes in list of time
  # first number is hour of the day
  # second number in list is minute of the hour
  num_time <- as.numeric(unlist(list_hr_min))
  num_time[1]+num_time[2]/60
}
d$t3_col_time_z01_cont <- sapply(time_split, cont_time)


################################################

##### Pregnancy biomarker data cleaning

######### Add combined ratios of cytokines ##########

#drop Z-score, sd, and ratio measures
d <- d[,!(grepl("^(z_)",colnames(d)) | grepl("^(sd_)",colnames(d)))]

x=c("il1_mom_t0", "il6_mom_t0", "tnfa_mom_t0")[1]
summary(as.vector(scale(d[,x], center = FALSE, scale = apply(as.matrix(d[,x]), 2, sd, na.rm = TRUE))))
x=c("il1_mom_t0", "il6_mom_t0", "tnfa_mom_t0")[2]
summary(as.vector(scale(d[,x], center = FALSE, scale = apply(as.matrix(d[,x]), 2, sd, na.rm = TRUE))))
x=c("il1_mom_t0", "il6_mom_t0", "tnfa_mom_t0")[3]
summary(as.vector(scale(d[,x], center = FALSE, scale = apply(as.matrix(d[,x]), 2, sd, na.rm = TRUE))))


#function to create composite score
create_score <- function(d, numerator_vars=c("il1_mom_t0", "il6_mom_t0", "tnfa_mom_t0"), denominator_vars="il10_mom_t0", varname="mom_t0_ratio_pro_il10"){
  for(i in numerator_vars){
    if(i==numerator_vars[1]){
      x = as.vector(scale(d[,i], center = FALSE, scale = apply(as.matrix(d[,i]), 2, sd, na.rm = TRUE)))
    }else{
      x = x + as.vector(scale(d[,i], center = FALSE, scale = apply(as.matrix(d[,i]), 2, sd, na.rm = TRUE)))
    }
  }
  summary(x)
  
  
  for(i in denominator_vars){
    if(i==denominator_vars[1]){
      y = as.vector(scale(d[,i], center = FALSE, scale = apply(as.matrix(d[,i]), 2, sd, na.rm = TRUE)))
    }else{
      y = y + as.vector(scale(d[,i], center = FALSE, scale = apply(as.matrix(d[,i]), 2, sd, na.rm = TRUE)))
    }
  }
  summary(y)
  
  score=log(x/y)
  summary(score)
  d$score <- score
  colnames(d)[ncol(d)] <- varname
  return(d)
}

# *Pro-inflammatory cytokines / IL-10
# *(IL-1 + IL-6 + TNF-a) / IL-10
d <- create_score(d, numerator_vars=c("il1_mom_t0", "il6_mom_t0", "tnfa_mom_t0"), denominator_vars="il10_mom_t0", varname="mom_t0_ratio_pro_il10")
summary(d$mom_t0_ratio_pro_il10)
ggplot(d, aes(x=mom_t0_ratio_pro_il10)) + geom_density()

# *Th1 / IL-10
# *(IL-12 + IFN) / IL-10
# gen mom_t0_ratio_th1_il10 = (il12_mom_t0 + ifng_mom_t0) / il10_mom_t0
d <- create_score(d, numerator_vars=c("il12_mom_t0", "ifng_mom_t0"), denominator_vars="il10_mom_t0", varname="mom_t0_ratio_th1_il10")
summary(d$mom_t0_ratio_th1_il10)
ggplot(d, aes(x=mom_t0_ratio_th1_il10)) + geom_density()

# *Th2 / IL-10 
# *(IL-4 + IL-5 + IL-13) / IL-10
# gen mom_t0_ratio_th2_il10 = (il4_mom_t0 + il5_mom_t0 + il13_mom_t0) / il10_mom_t0
d <- create_score(d, numerator_vars=c("il4_mom_t0", "il5_mom_t0", "il13_mom_t0"), denominator_vars="il10_mom_t0", varname="mom_t0_ratio_th2_il10")
summary(d$mom_t0_ratio_th2_il10)
ggplot(d, aes(x=mom_t0_ratio_th2_il10)) + geom_density()


# *Th17 / IL-10
# *(IL-17A + IL-21) / IL-10
# gen mom_t0_ratio_th17_il10 = (il17_mom_t0 + il21_mom_t0) / il10_mom_t0
d <- create_score(d, numerator_vars=c("il17_mom_t0", "il21_mom_t0"), denominator_vars="il10_mom_t0", varname="mom_t0_ratio_th17_il10")
summary(d$mom_t0_ratio_th17_il10)
ggplot(d, aes(x=mom_t0_ratio_th17_il10)) + geom_density()


# *Th1 / Th2
# *(IL-12 + IFN) / (IL-4 + IL-5 + IL-13)
# gen mom_t0_ratio_th1_th2 = (il12_mom_t0 + ifng_mom_t0) / (il4_mom_t0 + il5_mom_t0 + il13_mom_t0)
d <- create_score(d, numerator_vars=c("il12_mom_t0", "ifng_mom_t0"), denominator_vars=c("il4_mom_t0", "il5_mom_t0", "il13_mom_t0"), varname="mom_t0_ratio_th1_th2")
summary(d$mom_t0_ratio_th1_th2)
ggplot(d, aes(x=mom_t0_ratio_th1_th2)) + geom_density()


# *Th1 / Th17
# *(IL-12+IFN) / (IL-17A + IL-21)
# gen mom_t0_ratio_th1_th17 = (il12_mom_t0 + ifng_mom_t0) / (il17_mom_t0 + il21_mom_t0)
d <- create_score(d, numerator_vars=c("il12_mom_t0", "ifng_mom_t0"), denominator_vars=c("il17_mom_t0", "il21_mom_t0"), varname="mom_t0_ratio_th1_th17")
summary(d$mom_t0_ratio_th1_th17)
ggplot(d, aes(x=mom_t0_ratio_th1_th17)) + geom_density()



########### Add deficiency cutoff exposures #############
# 1 if deficient, 0 if not deficient
d$vit_A_def <- ifelse(d$RBP_inf_preg < 0.83, 1, 0)
d$vit_D_def <- ifelse(d$vitD_nmol_per_L < 30, 1, 0)
d$iron_def <- ifelse(d$FERR_inf_preg < 12 | d$STFR_inf_preg > 8.3, 1, 0)


############## Merge in maternal sum scores ####################
d_sum <- read.csv(paste0(dropboxDir,"Data/Cleaned/Audrie/maternal inflammation sum score.csv"))%>%
  select(-X)


d <- left_join(d, d_sum, by="dataid")



############# Check covariate missingness ###################
exp <- c("vitD_nmol_per_L", "logFERR_inf", "logSTFR_inf", "logRBP_inf",
         "ln_preg_cort", "logCRP", "logAGP","ifng_mom_t0", "sumscore_t0_mom_Z", "ln_preg_estri")

out <- c("t2_f2_8ip",  "t2_f2_23d", "t2_f2_VI", "t2_f2_12i", "t2_iso_pca", "t3_cort_slope", "t3_residual_cort", "t3_cort_z01", "t3_cort_z03",
         "t3_saa_slope", "t3_residual_saa", "t3_saa_z01", "t3_saa_z02",
         "t3_map", "t3_hr_mean",
         "t3_gcr_mean", "t3_gcr_cpg12") 

Wvars<-c("sex","birthord", "momage","momheight","momedu", 
         "hfiacat", "Nlt18","Ncomp", "watmin", "walls", "floor", "roof", "HHwealth",
         "tr", "life_viol_any_t3", "viol_any_preg", "ageday_ut2", "ageday_t3_vital", "ageday_t3_oragene", "ageday_t3_salimetrics",  
         "month_blood_t0", "month_ut2", "month_vt3", "month_ot3", "month_lt3" ) %>% unique()

#Changed times from _ht2/3 to _at2/3 -> whats diff?
W <- d %>% select(all_of(Wvars))  


miss <- data.frame(name = names(W), missing = colSums(is.na(W))/nrow(W), row.names = c(1:ncol(W)))
for (i in 1:nrow(miss)) {
  miss$class[i] <- class(W[,which(colnames(W) == miss[i, 1])])
}
miss

for (w in Wvars){
  print(w)
  if(is.numeric(W[,w])){
    print(summary(W[,w]))
  }
  else{print(table(W[,w]))}
}

# roof has low variability
mean(W$roof, na.rm=T)
sd(W$roof, na.rm=T)
# remove roof from covariates

# add missingness category to IPV covariates
d$life_viol_any_t3<-as.factor(d$life_viol_any_t3)
d$life_viol_any_t3<-addNA(d$life_viol_any_t3)
levels(d$life_viol_any_t3)[length(levels(d$life_viol_any_t3))]<-"Missing"

d$viol_any_preg<-as.factor(d$viol_any_preg)
d$viol_any_preg<-addNA(d$viol_any_preg)
levels(d$viol_any_preg)[length(levels(d$viol_any_preg))]<-"Missing"

for(outcome in out){
  d_sub <- subset(d, !is.na(d[,outcome]))
  W_sub <- d_sub %>% select(all_of(Wvars))  
  
  miss_sub <- data.frame(name = names(W_sub), missing = colSums(is.na(W_sub)), row.names = c(1:ncol(W_sub)))
  for (i in 1:nrow(miss_sub)) {
    miss_sub$class[i] <- class(W_sub[,which(colnames(W_sub) == miss_sub[i, 1])])
  }
  print(outcome)
  print(miss_sub)
}



saveRDS(d, paste0(dropboxDir,"Data/Cleaned/Audrie/pregnancy_stress_data.RDS"))


