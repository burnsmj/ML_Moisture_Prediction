### Non-Parallel Learning Curve on MSI
### Michael Burns
### 12/4/2020

# NOTE: It may be possible to write an if statement before the return function to detect if an iteration ran. If not, a simple NA could be inserted, and later removed before graphing.

#############
# Libraries #
#############
library(caret)
library(tidyverse)
library(readxl)

################
# Loading Data #
################
###############################
# Loading Macro Training Data #
###############################
macro_training_data <- read_csv("../Data/Moisture_Uptake_Master_Dataset_Cook_Macro_Spectra.csv") %>%
  select(c(1,2,21:26))

##################################
# Loading Spectral Training Data #
##################################
spectra_training_data <- read_csv("../Data/Moisture_Uptake_Master_Dataset_Cook_Macro_Spectra.csv") %>%
  select(-c(3:25))

Vis <- spectra_training_data
#Create a matrix with lower and upper limits equal to the min and max wavelengths in your dataset
#Select every 10th (you could also do every 5th, 3rd, etc.)wavelength to remove collinearity issues for outlier detection
minwv <- 950 #What is the minimum wavelength you scanned
maxwv <- 1650#What is the maxium wavelength you scanned
Rand <- as.matrix(seq(minwv,maxwv,10))
#Combine the sample name with the subset of columns selected
MData<-cbind(Vis[,1],Vis[,colnames(Vis) %in% Rand])
colnames(MData)[1]<-"Sample"
#Calculate the mahalanobis distance on all samples using subset of wavelengths
Spec_Start <- 4#The column where the spectra data starts in your dataset
m_dist <- as.matrix(mahalanobis(MData[,Spec_Start:ncol(MData)], colMeans(MData[,Spec_Start:ncol(MData)]), cov(MData[,Spec_Start:ncol(MData)]),na.rm=TRUE))
df <- as.data.frame(round(m_dist, 1))
colnames(df)[1]<-"MH"
#Combine the Mahalanobis distance with the original data
Combined <- cbind(df,Vis)
#Calculate  threshold for sample to be considered an outlier (here I use 3 times the mean)
Threshold <- 3*mean(m_dist)
#Only reatin samples that are less than outlier threshold
Final <- Combined[Combined[,1]<Threshold,]
#Delete column that has Mahalanobis distance
clean_spectra_training_data <- Final[,2:ncol(Final)]

norm_training_data <- clean_spectra_training_data %>%
  pivot_longer(cols = -c(Sample_ID, Genotype, Moisture_Uptake), names_to = "Waveband", values_to = "Absorbance") %>%
  group_by(Sample_ID) %>%
  mutate(sum_lxl_abs = sum(abs(Absorbance)),
         Norm_Abs = Absorbance / sum_lxl_abs) %>%
  select(-c(Absorbance, sum_lxl_abs)) %>%
  ungroup() %>%
  pivot_wider(id_cols = c(Sample_ID, Genotype, Moisture_Uptake), 
              values_from = Norm_Abs, 
              names_from = Waveband)

########################
# Macro Learning Curve #
########################
lc_reps_macro <- data.frame(NULL)
for(i in 1:100){
  set.seed(i)
  lc <- learning_curve_dat(dat = as.data.frame(macro_training_data) %>%
                             select_if(is.numeric),
                           outcome = "Moisture_Uptake",
                           test_prop = 0.1,
                           proportion = 5:20/20,
                           method = "lm", 
                           metric = "Rsquared", 
                           tuneGrid = expand.grid(intercept = 0.2),
                           trControl = trainControl(method = "cv",
                                                    index = groupKFold(macro_training_data$Genotype, k = 10), 
                                                    savePredictions = T, 
                                                    #predictionBounds = c(T,T),
                                                    allowParallel = T
                           )
  ) 
  lc_tibble <- lc %>%
    as_tibble() %>%
    filter(Data != "Resampling")
    
  lc_reps_macro <- rbind(lc_reps_macro, lc_tibble)
}

###################
# Writing Results #
###################
lc_reps_macro %>%
  as_tibble() %>%
  unnest(everything()) %>%
  filter(Data != "Resampling") %>%
  write_csv("../Results/LC_MSI_Macro_Reps_NP.csv")

##########################
# Spectra Learning Curve #
##########################
lc_reps_spectra <- data.frame(NULL)
for(i in 1:100){
  lc <- learning_curve_dat(dat = as.data.frame(norm_training_data) %>%
                             select_if(is.numeric),
                           outcome = "Moisture_Uptake",
                           test_prop = 0.1,
                           proportion = 5:20/20,
                           method = "svmLinear", 
                           metric = "Rsquared", 
                           tuneGrid = expand.grid(C = 71.407),
                           trControl = trainControl(method = "cv",
                                                    index = groupKFold(norm_training_data$Genotype, k = 10), 
                                                    savePredictions = T, 
                                                    #predictionBounds = c(T,T),
                                                    allowParallel = T
                           )
  ) 
  lc_tibble <- lc %>%
    as_tibble() %>%
    filter(Data != "Resampling")
  
  lc_reps_spectra <- rbind(lc_reps_spectra, lc_tibble)
}

###################
# Writing Results #
###################
lc_reps_spectra %>%
  as_tibble() %>%
  unnest(everything()) %>%
  filter(Data != "Resampling") %>%
  write_csv("../Results/LC_MSI_Spectra_Reps_NP.csv")
