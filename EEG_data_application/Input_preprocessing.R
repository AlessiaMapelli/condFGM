# Required packages
suppressPackageStartupMessages(library(yaml))
suppressPackageStartupMessages(library(abind))

############################################
#  INPUT PREPROCESSING: BIND and REFORMAT  #
############################################
load("EEG_data_application/alco_filtered_array.Rdata")
load("EEG_data_application/contrl_filtered_array.Rdata")
dim(alco.filtered.array)
dim(contrl.filtered.array)

discrete_fun_obs <- abind(contrl.filtered.array, alco.filtered.array, along = 1)
covariates_df <- data.frame(group=c(rep("control", dim(contrl.filtered.array)[1]),rep("alcoholic",dim(alco.filtered.array)[1] )))
covariates_df$group <- factor(covariates_df$group, levels= c("control", "alcoholic"))
save(discrete_fun_obs,covariates_df, file="EEG_data_application/EEG_application_input_data.RData")




###### TEST RESULTS #########
covariates_df <- data.frame(group=c(rep("control", dim(contrl.filtered.array)[1]),rep("alcoholic",dim(alco.filtered.array)[1] )))
covariates_df$group <- factor(covariates_df$group, levels= c("alcoholic", "control"))
save(discrete_fun_obs,covariates_df, file="EEG_data_application/EEG_application_input_data_test.RData")
