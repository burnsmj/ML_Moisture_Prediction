#Cleaner + Predictor

#####
# Required Packages
#####
library(pls)
library(dplyr)
library(ggplot2)
library(reshape2)
library(Hmisc)

#####
# Loading Training Set
#####
nir_training_spectra_full<-read_xlsx("Data/Spectra_Data/spectra_learning_sets.xlsx")
nir_training_spectra <- nir_training_spectra_full[,-c(3:23,145)]
head(nir_training_spectra)
#####
# Cleaning Prediction Set
#####
data_w_outlier<-read_csv("Data/Spectra_Data/WiDiv_Spectra_w_Time.csv")

outlier_rm<-function(data_w_outlier){
  inliers<-na.omit(data_w_outlier)
  for(q in 1:ncol(inliers)){
    if(is.double(inliers[,q])){
      data_mean<-mean(inliers[,q])
      data_sd<-sd(inliers[,q])
      for(n in 1:length(inliers[,1])){
        if(inliers[n,q]>(data_mean-(3*data_sd)) & inliers[n,q]<(data_mean+(3*data_sd))){
          inliers[n,q]<-inliers[n,q]
        }else{inliers[n,q]<-NA}
      }
    }
  }
  na.omit(inliers)
}

data_clean<-outlier_rm(data_w_outlier)

#####
# Removing Duplicates
#####
dim(data_clean)
dataset_sorted<-data_clean[order(data_clean$SampleID),]
for(i in 2:nrow(dataset_sorted)){
  if(is.na(dataset_sorted$SampleID[i-1])){}else{
    if(dataset_sorted$SampleID[i]==dataset_sorted$SampleID[i-1]){
      if(dataset_sorted$Scan_Time[i]>dataset_sorted$Scan_Time[i-1]){
        dataset_sorted$SampleID[i-1]<-NA
      }else{
        if(dataset_sorted$Scan_Time[i]==dataset_sorted$Scan_Time[i-1]){
          dataset_sorted$SampleID[i-1]<-NA
        }else{dataset_sorted$SampleID[i]<-NA}
      }
    }}
}
dataset_sorted<-na.omit(dataset_sorted)
dim(dataset_sorted)

#####
# PLS Learning
#####
spectra_readings<-as.matrix(nir_training_spectra[,c(3:143)]) #having the spectra in a matrix will allow a cleaner code when calling them as predictors.
moisture_spectra_model<-plsr(nir_training_spectra$Moisture_Avg~spectra_readings, ncomp = 141, validation = "CV")

#####
# Choosing Components
#####
R2_Y<-R2(moisture_spectra_model, estimate = "train")
R2_Y_max<-max(R2_Y$val)
R2_Y_mincomp<-R2_Y$val<0.8*R2_Y_max
R2_Y_compnum<-sum(R2_Y_mincomp)+1

RMSEP_mod<-RMSEP(moisture_spectra_model)
RMSEP_val<-RMSEP_mod$val
RMSEP_val_CV<-RMSEP_val[c(T,F)]
RMSEP_min<-which.min(RMSEP_val_CV)-1 # -1 corrects for intercept being considered in the RMSEP values list

plot(RMSEP_mod, legendpos = "topright", 
     main = "RMSE vs Number of Components Used", 
     ylab = "RMSE",
     xlab = "Number of Components",
     sub = paste("Minimum RMSE found with", RMSEP_min, "components", sep = " "),
     col.sub = "blue")
abline(v = RMSEP_min, col = "blue")
abline(v = RMSEP_min+10, col = 'red')

plot(R2_Y, 
     main = "R^2 vs Number of Components Used", 
     ylab = "R^2",
     xlab = "Number of Components",
     sub = paste(round(0.8*R2_Y_max*100,1), "% of variation found with ", R2_Y_compnum, " components", sep = ""),
     col.sub = "blue")
abline(v = R2_Y_compnum, col = "blue")
abline(v = RMSEP_min+10, col = 'red')

ncomp_onesigma<-selectNcomp(moisture_spectra_model, 
                            method = "onesigma", 
                            plot = T, 
                            xlim = c(0,15), 
                            sub = paste("Variance Explained: ", 
                                        round(R2_Y$val[selectNcomp(moisture_spectra_model,
                                                                   method ="onesigma")+1],3)*100,"%", sep =""),
                            col.sub = "blue")

ncomp_permut<-selectNcomp(moisture_spectra_model, 
                          method = "randomization", 
                          plot = T, 
                          xlim = c(0,15),
                          sub = paste("Variance Explained: ", 
                                      round(R2_Y$val[selectNcomp(moisture_spectra_model,
                                                                 method ="randomization")+1],3)*100,"%", sep =""),
                          col.sub = "blue")
# +1's corrects for intercept being included in the R^2 dataset.

###########################################################################################################################
##### Pause here and look at the plots made to justify the ncomp value found below
###########################################################################################################################

#####
# PLS Predictions
#####
nir_pred_spectra_mat<-as.matrix(dataset_sorted[,c(3:143)])
moisture_predictions_pls<-drop(predict(moisture_spectra_model, ncomp = 15, newdata = nir_pred_spectra_mat))

#####
# Printing Results
#####
moisture_predictions_w_spectra<-cbind(dataset_sorted,moisture_predictions_pls)
moisture_predictions_w_spectra_reordered<-moisture_predictions_w_spectra[,c(1,144,2:143)]
write.csv(moisture_predictions_w_spectra_reordered, "/Users/michael/Desktop/Grad_School/Research/Datasets/Machine Learning/Spectra_Data/Predictions/WiDiv_Spectra_Predictions.csv", row.names = F)


#####
# Picking 40 Individuals
#####
moisture_predictions_w_spectra_reordered_sorted<-moisture_predictions_w_spectra_reordered[order(moisture_predictions_w_spectra_reordered$moisture_predictions_pls),]
pick_every<-nrow(moisture_predictions_w_spectra_reordered_sorted)/40
samples_to_use<-as.matrix(moisture_predictions_w_spectra_reordered_sorted[c(T,rep(F,pick_every)),c("SampleID","moisture_predictions_pls")])
write.csv(samples_to_use, "/Users/michael/Desktop/Grad_School/Research/Datasets/Machine Learning/Spectra_Data/Predictions/WiDiv_Samples_to_Cook.csv", row.names = F)


#####
# Extras
#####
hist(moisture_predictions_pls)


samples_to_use$moisture_predictions_pls <- as.numeric(samples_to_use$moisture_predictions_pls)
samples_to_use[,moisture_predictions_pls]
is.matrix(samples_to_use)
as.numeric(samples_to_use[,2])
ggplot()+
  geom_histogram(data = moisture_predictions_w_spectra, aes(x = moisture_predictions_pls))+
  geom_rug(aes(x = as.numeric(samples_to_use[,2]), color = "red"))+
  xlab("Moisture Uptake")+
  labs(title = "Predicted Moisture Uptake of WiDiv Samples")



