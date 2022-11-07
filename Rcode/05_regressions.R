#source("Rcode/00_pkg_functions.R")
source("Rcode/01_data_prep.R")#
source("Rcode/02_FD.R")#
source("Rcode/03_PD.R")# 
source("Rcode/04_combine.R") #

df = na.omit(data_raptor)# deletes pixels with SR <2 because no FD nor PD calculations are possible for these communities

df[, 38:47] = scale(df[, 38:47])

#correlations among environmental variables
pairs(df[, 38:47])
Hmisc::rcorr(as.matrix(df[, 38:47]),type="pearson")

## using partitioning 
library(remotePARTS)
library(dplyr)
library(ggplot2)

#adding quadratic terms

df$elevation2<-df$elevation_sd^2
df$hfp2<-df$hfp_mean ^2
df$npp2<-df$npp_mean_2000_2015^2
df$paleotemp2<-df$paleo_temp^2

# partitioning data in 4 subsets

pm <- sample_partitions(npix = nrow(df), npart = 4) #dividing in 4 parts ~3000 pixeles per part
dim(pm)#3175 pixeles per partition

# fitting GLS with potentially different nuggets for partition (the more flexible option)

## FULL MODELS 

names(df)

rich_full = fitGLS_partition(formula = sr ~ elevation_sd+elevation2+hfp_mean+hfp2+npp_mean_2000_2015+npp2+paleo_temp+paleotemp2,data = df,
                             partmat = pm,covar_FUN = "covar_exp", nugget = NA, 
                             coord.names=c("lon_wgs84", "lat_wgs84"), parallel = TRUE, ncores= 60)

sesfdis_morph_full = fitGLS_partition(formula = fdis_ses_morph ~ elevation_sd+elevation2+hfp_mean+hfp2+npp_mean_2000_2015+npp2+paleo_temp+paleotemp2,data = df,
                             partmat = pm,covar_FUN = "covar_exp", nugget = NA, 
                             coord.names=c("lon_wgs84", "lat_wgs84"), parallel = TRUE, ncores= 60)

sesfdis_niche_full = fitGLS_partition(formula = fdis_ses_niche ~ elevation_sd+elevation2+hfp_mean+hfp2+npp_mean_2000_2015+npp2+paleo_temp+paleotemp2,data = df,
                                      partmat = pm,covar_FUN = "covar_exp", nugget = NA, 
                                      coord.names=c("lon_wgs84", "lat_wgs84"), parallel = TRUE, ncores= 60)

sesfdis_dispe_full = fitGLS_partition(formula = fdis_ses_dispe ~ elevation_sd+elevation2+hfp_mean+hfp2+npp_mean_2000_2015+npp2+paleo_temp+paleotemp2,data = df,
                                      partmat = pm,covar_FUN = "covar_exp", nugget = NA, 
                                      coord.names=c("lon_wgs84", "lat_wgs84"), parallel = TRUE, ncores= 60)

sesfdis_diet_full = fitGLS_partition(formula = fdis_ses_diet ~ elevation_sd+elevation2+hfp_mean+hfp2+npp_mean_2000_2015+npp2+paleo_temp+paleotemp2,data = df,
                                      partmat = pm,covar_FUN = "covar_exp", nugget = NA, 
                                      coord.names=c("lon_wgs84", "lat_wgs84"), parallel = TRUE, ncores= 60)
sesfdis_forag_full = fitGLS_partition(formula = fdis_ses_forag ~ elevation_sd+elevation2+hfp_mean+hfp2+npp_mean_2000_2015+npp2+paleo_temp+paleotemp2,data = df,
                                      partmat = pm,covar_FUN = "covar_exp", nugget = NA, 
                                      coord.names=c("lon_wgs84", "lat_wgs84"), parallel = TRUE, ncores= 60)

sespd_full = fitGLS_partition(formula = pd.uroot.z ~ elevation_sd+elevation2+hfp_mean+hfp2+npp_mean_2000_2015+npp2+paleo_temp+paleotemp2,data = df,
                                      partmat = pm,covar_FUN = "covar_exp", nugget = NA, 
                                      coord.names=c("lon_wgs84", "lat_wgs84"), parallel = TRUE, ncores= 60)

sesmpd_full = fitGLS_partition(formula = mpd.z ~ elevation_sd+elevation2+hfp_mean+hfp2+npp_mean_2000_2015+npp2+paleo_temp+paleotemp2,data = df,
                              partmat = pm,covar_FUN = "covar_exp", nugget = NA, 
                              coord.names=c("lon_wgs84", "lat_wgs84"), parallel = TRUE, ncores= 60)

sesmntd_full = fitGLS_partition(formula = mntd.z ~ elevation_sd+elevation2+hfp_mean+hfp2+npp_mean_2000_2015+npp2+paleo_temp+paleotemp2,data = df,
                               partmat = pm,covar_FUN = "covar_exp", nugget = NA, 
                               coord.names=c("lon_wgs84", "lat_wgs84"), parallel = TRUE, ncores= 60)


## overall estimates 
rich_full$overall
sesfdis_morph_full$overall
sesfdis_niche_full$overall
sesfdis_dispe_full$overall
sesfdis_diet_full$overall
sesfdis_forag_full$overall
sespd_full$overall
sesmpd_full$overall
sesmntd_full$overall


# FINAL MODELS (ELIMINATING QUADRATIC TERMS WHEN NON-SIGNIFICANT, TO SIMPLIFY MODELS) # ---------------------------------------------------------

names(df)

rich_full2 = fitGLS_partition(formula = sr ~ elevation_sd+elevation2+hfp_mean+hfp2+npp_mean_2000_2015+paleo_temp,data = df,
                             partmat = pm,covar_FUN = "covar_exp", nugget = NA, 
                             coord.names=c("lon_wgs84", "lat_wgs84"), parallel = TRUE, ncores= 60)

sesfdis_morph_full2 = fitGLS_partition(formula = fdis_ses_morph ~ elevation_sd+elevation2+hfp_mean+npp_mean_2000_2015+npp2+paleo_temp,data = df,
                                      partmat = pm,covar_FUN = "covar_exp", nugget = NA, 
                                      coord.names=c("lon_wgs84", "lat_wgs84"), parallel = TRUE, ncores= 60)

sesfdis_niche_full2 = fitGLS_partition(formula = fdis_ses_niche ~ elevation_sd+elevation2+hfp_mean+hfp2+npp_mean_2000_2015+npp2+paleo_temp,data = df,
                                      partmat = pm,covar_FUN = "covar_exp", nugget = NA, 
                                      coord.names=c("lon_wgs84", "lat_wgs84"), parallel = TRUE, ncores= 60)

sesfdis_dispe_full2 = fitGLS_partition(formula = fdis_ses_dispe ~ elevation_sd+hfp_mean+hfp2+npp_mean_2000_2015+npp2+paleo_temp,data = df,
                                      partmat = pm,covar_FUN = "covar_exp", nugget = NA, 
                                      coord.names=c("lon_wgs84", "lat_wgs84"), parallel = TRUE, ncores= 60)

sesfdis_diet_full2 = fitGLS_partition(formula = fdis_ses_diet ~ elevation_sd+elevation2+hfp_mean+npp_mean_2000_2015+npp2+paleo_temp+paleotemp2,data = df,
                                     partmat = pm,covar_FUN = "covar_exp", nugget = NA, 
                                     coord.names=c("lon_wgs84", "lat_wgs84"), parallel = TRUE, ncores= 60)

sesfdis_forag_full2 = fitGLS_partition(formula = fdis_ses_forag ~ elevation_sd+elevation2+hfp_mean+npp_mean_2000_2015+paleo_temp+paleotemp2,data = df,
                                      partmat = pm,covar_FUN = "covar_exp", nugget = NA, 
                                      coord.names=c("lon_wgs84", "lat_wgs84"), parallel = TRUE, ncores= 60)

sespd_full2 = fitGLS_partition(formula = pd.uroot.z ~ elevation_sd+hfp_mean+npp_mean_2000_2015+npp2+paleo_temp,data = df,
                              partmat = pm,covar_FUN = "covar_exp", nugget = NA, 
                              coord.names=c("lon_wgs84", "lat_wgs84"), parallel = TRUE, ncores= 60)

sesmpd_full2 = fitGLS_partition(formula = mpd.z ~ elevation_sd+hfp_mean+hfp2+npp_mean_2000_2015+paleo_temp+paleotemp2,data = df,
                               partmat = pm,covar_FUN = "covar_exp", nugget = NA, 
                               coord.names=c("lon_wgs84", "lat_wgs84"), parallel = TRUE, ncores= 60)

sesmntd_full2 = fitGLS_partition(formula = mntd.z ~ elevation_sd+hfp_mean+npp_mean_2000_2015+npp2+paleo_temp+paleotemp2,data = df,
                                partmat = pm,covar_FUN = "covar_exp", nugget = NA, 
                                coord.names=c("lon_wgs84", "lat_wgs84"), parallel = TRUE, ncores= 60)


## overall estimates 
rich_full2$overall
sesfdis_morph_full2$overall
sesfdis_niche_full2$overall
sesfdis_dispe_full2$overall
sesfdis_diet_full2$overall
sesfdis_forag_full2$overall
sespd_full2$overall
sesmpd_full2$overall
sesmntd_full2$overall