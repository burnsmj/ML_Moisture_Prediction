#Duplicate Remover

###STILL NEEDS WORK!

#####
# Load Data
#####
dataset<-read.csv("/Users/michael/Desktop/Grad_School/Research/Datasets/Machine Learning/Spectra_Data/WiDiv_Spectra_w_Time.csv", check.names = F, stringsAsFactors = F)
dim(dataset)
dataset$Scan_Time<-as.numeric(dataset$Scan_Time)

#####
# Order based on SampleID
#####
dataset_sorted<-dataset[order(dataset$SampleID),]
i<-4684
for(i in 2:nrow(dataset_sorted)){
  print(i)
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

#####
# Removing duplicates
#####
dataset_sorted<-na.omit(dataset_sorted)
dim(dataset_sorted)

#####
# Exporting as .csv file
#####
write.csv(dataset_sorted, "/Users/michael/Desktop/Grad_School/Research/Datasets/Machine Learning/Spectra_Data/WiDiv_wout_Dups.csv", row.names = F)
