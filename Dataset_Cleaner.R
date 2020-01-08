#Dataset Cleaner

data_w_outlier<-read.csv("/Users/michael/Desktop/Grad_School/Research/Datasets/Machine Learning/Spectra_Data/WiDiv_Spectra_For_Predictions.csv", check.names = F)

outlier_rm_3sd<-function(data_w_outlier){
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

data_clean<-outlier_rm_3sd(data_w_outlier)

write.csv(data_clean,"/Users/michael/Desktop/Grad_School/Research/Datasets/Machine Learning/Spectra_Data/WiDiv_Dataset_Inliers.csv", row.names = F)
