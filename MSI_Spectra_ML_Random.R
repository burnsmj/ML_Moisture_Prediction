### MSI Spectra Machine Learning -- Random
### Michael Burns
### 2020-07-08

# Perhaps the most succint way to do this will be to put a parallel foreach loop inside a for loop that will iter through a list of the normalized
# spectra that we have.  This loop can then be inside a larger loop that will repeat it a number of times.

#############
# Libraries #
#############
library(caret)
library(pls)
library(parallel)
library(doParallel)
library(foreach)
library(ggplot2)
library(broom)
library(tidyverse)
library(readxl)
library(dplyr)
library(kernlab)
library(randomForest)


print("Random Partitioning")

################
# Loading Data #
################
dataset <- read_csv("../Data/Moisture_Uptake_Master_Dataset_Cook_Macro_Spectra.csv") %>%
  select(-c(3:4,6:25))

##############################
# Removing Spectral Outliers #
##############################
#This code was developed and given to me by a PepsiCo statistician.
#Read in your raw NIR matrix (unless you decided to pre-process the spectra already)
Vis<-dataset

#Create a matrix with lower and upper limits equal to the min and max wavelengths in your dataset
#Select every 10th (you could also do every 5th, 3rd, etc.)wavelength to remove collinearity issues for outlier detection
minwv<- 950 #What is the minimum wavelength you scanned
maxwv<- 1650#What is the maxium wavelength you scanned
Rand<-as.matrix(seq(minwv,maxwv,10))

#Combine the sample name with the subset of columns selected
MData<-cbind(Vis[,1],Vis[,colnames(Vis) %in% Rand])
colnames(MData)[1]<-"Sample"

#Calculate the mahalanobis distance on all samples using subset of wavelengths
Spec_Start<-4#The column where the spectra data starts in your dataset
m_dist <- as.matrix(mahalanobis(MData[,Spec_Start:ncol(MData)], colMeans(MData[,Spec_Start:ncol(MData)]), cov(MData[,Spec_Start:ncol(MData)]),na.rm=TRUE))
df <- as.data.frame(round(m_dist, 1))
colnames(df)[1]<-"MH"

#Combine the Mahalanobis distance with the original data
Combined<-cbind(df,Vis)

#Calculate  threshold for sample to be considered an outlier (here I use 3 times the mean)
Threshold<-3*mean(m_dist)

#Only reatin samples that are less than outlier threshold
Final<-Combined[Combined[,1]<Threshold,]

#Delete column that has Mahalanobis distance
dataset<-Final[,2:ncol(Final)]

###########################
# Splitting Data Randomly #
###########################
for(n in 1:20) {
  ######################
  # Splitting Randomly #
  ######################
  training_index <- createDataPartition(dataset$Moisture_Uptake,
                                        times = 1,
                                        p = 0.8,
                                        list = FALSE)
  training_set <- dataset[training_index,]
  validation_set <- dataset[-training_index,]

  #######################
  # Normalizing Spectra #
  #######################
  ### Normalization by maximum absorbance ###
  max_train <- training_set %>%
    pivot_longer(cols = -c(Sample_ID, Genotype, Env, Moisture_Uptake), names_to = "Waveband", values_to = "Absorbance") %>%
    group_by(Sample_ID) %>%
    mutate(Max_Abs = max(Absorbance),
           Norm_Abs = Absorbance / Max_Abs) %>%
    select(-c(Absorbance, Max_Abs)) %>%
    ungroup() %>%
    pivot_wider(id_cols = c(Sample_ID, Genotype, Env, Moisture_Uptake), values_from = Norm_Abs, names_from = Waveband)
  
  max_valid <- validation_set %>%
    pivot_longer(cols = -c(Sample_ID, Genotype, Env, Moisture_Uptake), names_to = "Waveband", values_to = "Absorbance") %>%
    group_by(Sample_ID) %>%
    mutate(Max_Abs = max(Absorbance),
           Norm_Abs = Absorbance / Max_Abs) %>%
    select(-c(Absorbance, Max_Abs)) %>%
    ungroup() %>%
    pivot_wider(id_cols = c(Sample_ID, Genotype, Env, Moisture_Uptake), values_from = Norm_Abs, names_from = Waveband)
  
  ### Normalization by sum of abs value ###
  abs_train <- training_set %>%
    pivot_longer(cols = -c(Sample_ID, Genotype, Env, Moisture_Uptake), names_to = "Waveband", values_to = "Absorbance") %>%
    group_by(Sample_ID) %>%
    mutate(sum_lxl_abs = sum(abs(Absorbance)),
           Norm_Abs = Absorbance / sum_lxl_abs) %>%
    select(-c(Absorbance, sum_lxl_abs)) %>%
    ungroup() %>%
    pivot_wider(id_cols = c(Sample_ID, Genotype, Env, Moisture_Uptake), 
                values_from = Norm_Abs, 
                names_from = Waveband)
  
  abs_valid <- validation_set %>%
    pivot_longer(cols = -c(Sample_ID, Genotype, Env, Moisture_Uptake), names_to = "Waveband", values_to = "Absorbance") %>%
    group_by(Sample_ID) %>%
    mutate(sum_lxl_abs = sum(abs(Absorbance)),
           Norm_Abs = Absorbance / sum_lxl_abs) %>%
    select(-c(Absorbance, sum_lxl_abs)) %>%
    ungroup() %>%
    pivot_wider(id_cols = c(Sample_ID, Genotype, Env, Moisture_Uptake), 
                values_from = Norm_Abs, 
                names_from = Waveband)
  
  ### Normalization by sum of sq ###
  sq_train <- training_set %>%
    pivot_longer(cols = -c(Sample_ID, Genotype, Env, Moisture_Uptake), names_to = "Waveband", values_to = "Absorbance") %>%
    group_by(Sample_ID) %>%
    mutate(sum_sq_abs = sum((Absorbance)^2),
           Norm_Abs = Absorbance / sum_sq_abs) %>%
    select(-c(Absorbance, sum_sq_abs)) %>%
    ungroup() %>%
    pivot_wider(id_cols = c(Sample_ID, Genotype, Env, Moisture_Uptake), 
                values_from = Norm_Abs, 
                names_from = Waveband)
  
  sq_valid <- validation_set %>%
    pivot_longer(cols = -c(Sample_ID, Genotype, Env, Moisture_Uptake), names_to = "Waveband", values_to = "Absorbance") %>%
    group_by(Sample_ID) %>%
    mutate(sum_sq_abs = sum((Absorbance)^2),
           Norm_Abs = Absorbance / sum_sq_abs) %>%
    select(-c(Absorbance, sum_sq_abs)) %>%
    ungroup() %>%
    pivot_wider(id_cols = c(Sample_ID, Genotype, Env, Moisture_Uptake), 
                values_from = Norm_Abs, 
                names_from = Waveband)
  
  ### Standardization ###
  std_train <- training_set %>%
    pivot_longer(cols = -c(Sample_ID, Genotype, Env, Genotype, Env, Moisture_Uptake), names_to = "Waveband", values_to = "Absorbance") %>%
    group_by(Sample_ID) %>%
    mutate(st_abs = (Absorbance - mean(Absorbance)) / sd(Absorbance)) %>%
    select(-c(Absorbance)) %>%
    ungroup() %>%
    pivot_wider(id_cols = c(Sample_ID, Genotype, Env, Moisture_Uptake), 
                values_from = st_abs, 
                names_from = Waveband)
  
  std_valid <- validation_set %>%
    pivot_longer(cols = -c(Sample_ID, Genotype, Env, Moisture_Uptake), names_to = "Waveband", values_to = "Absorbance") %>%
    group_by(Sample_ID) %>%
    mutate(st_abs = (Absorbance - mean(Absorbance)) / sd(Absorbance)) %>%
    select(-c(Absorbance)) %>%
    ungroup() %>%
    pivot_wider(id_cols = c(Sample_ID, Genotype, Env, Moisture_Uptake), 
                values_from = st_abs, 
                names_from = Waveband)

  ###########################
  # List of Normalized Data #
  ###########################

  datasets_list <- list(base_training = training_set %>% as_tibble(),
                        base_validation = validation_set %>% as_tibble(),
                        max_training = max_train,
                        max_validation = max_valid,
                        sq_training = sq_train,
                        sq_validation = sq_valid,
                        abs_training = abs_train,
                        abs_validation = abs_valid,
                        std_training = std_train,
                        std_validation = std_valid)

  for(d in c(1,3,5,7,9)) { # Use [[d]] for training and [[d+1]] for validation.
    #############################
    # Specifying Dataset to Use #
    ############################# 
    train_data <- datasets_list[[d]]
    valid_data <- datasets_list[[d+1]]
    if(d == 3) { # removes all near zero variance columns (very important for svms).  Only for one column when max normalization is used.  Used a for loop because I got lazy and didn't want to spend so much time figuring out how to not subset when there wasnt a near zero variance column.
      train_data <- train_data[,-nearZeroVar(train_data)]
      valid_data <- valid_data[,-nearZeroVar(valid_data)]
    }

    num_cores_limit <- 47
    Spectra_Random_Results <- tibble()
    if (detectCores() > num_cores_limit) {
      num_cores <- num_cores_limit
    } else {
      num_cores <- detectCores()
      }
    registerDoParallel(cores = num_cores)
    
    spectra_random_results <- foreach(i = 1:141, .combine = rbind) %dopar% {
      ##########
      # Models #
      ##########
      ptm <- proc.time()
      set.seed(n)
      spectra_lm <- train(Moisture_Uptake ~ .,
                          data = train_data %>%
                            select_if(is.numeric) %>%
                            select(-Env),
                          method = "lm",
                          metric = "Rsquared",
                          tuneGrid = expand.grid(intercept = i/141),
                          trControl = trainControl(method = "cv",
                                                   index = groupKFold(train_data$Genotype, k = 10), 
                                                   savePredictions = T, 
                                                   predictionBounds = c(T,T),
                                                   allowParallel = T
                                                   )
                          )
      set.seed(n)
      spectra_pcr <- train(Moisture_Uptake ~ .,
                           data = train_data %>%
                             select_if(is.numeric) %>%
                             select(-Env),
                           method = "pcr",
                           metric = "Rsquared",
                           tuneGrid = expand.grid(ncomp = i),
                           trControl = trainControl(method = "cv",
                                                    index = groupKFold(train_data$Genotype, k = 10), 
                                                    savePredictions = T, 
                                                    predictionBounds = c(T,T),
                                                    allowParallel = T
                                                    )
                           )
      set.seed(n)
      spectra_pls <- train(Moisture_Uptake ~ .,
                           data = train_data %>%
                             select_if(is.numeric) %>%
                             select(-Env),
                           method = "pls",
                           metric = "Rsquared",
                           tuneGrid = expand.grid(ncomp = i),
                           trControl = trainControl(method = "cv",
                                                    index = groupKFold(train_data$Genotype, k = 10), 
                                                    savePredictions = T, 
                                                    predictionBounds = c(T,T),
                                                    allowParallel = T
                                                    )
                            )
      set.seed(n)
      spectra_rf <- train(Moisture_Uptake ~ .,
                          data = train_data %>%
                            select_if(is.numeric) %>%
                            select(-Env),
                          method = "rf",
                          metric = "Rsquared",
                          tuneGrid = expand.grid(mtry = i),
                          trControl = trainControl(method = "cv",
                                                   index = groupKFold(train_data$Genotype, k = 10), 
                                                   savePredictions = T, 
                                                   predictionBounds = c(T,T),
                                                   allowParallel = T
                                                   )
                           )
      set.seed(n)
      spectra_svml <- train(Moisture_Uptake ~ .,
                            data = train_data %>%
                              select_if(is.numeric) %>%
                              select(-Env),
                            method = "svmLinear",
                            metric = "Rsquared",
                            tuneGrid = expand.grid(C = ((i^2) * 0.007)),
                            trControl = trainControl(method = "cv",
                                                     index = groupKFold(train_data$Genotype, k = 10), 
                                                     savePredictions = T, 
                                                     predictionBounds = c(T,T),
                                                     allowParallel = T
                                                     )
                            )
      set.seed(n)
      spectra_svmp <- train(Moisture_Uptake ~ .,
                            data = train_data %>%
                              select_if(is.numeric) %>%
                              select(-Env),
                            method = "svmPoly",
                            metric = "Rsquared",
                            tuneGrid = expand.grid(C = ((i^2) * 0.007), degree = seq(1,3,1), scale = c(0.001,0.01,0.1)),
                            trControl = trainControl(method = "cv",
                                                     index = groupKFold(train_data$Genotype, k = 10), 
                                                     savePredictions = T, 
                                                     predictionBounds = c(T,T),
                                                     allowParallel = T
                                                     )
                            )
      set.seed(n)
      spectra_svmr <- train(Moisture_Uptake ~ .,
                            data = train_data %>%
                              select_if(is.numeric) %>%
                              select(-Env),
                            method = "svmRadial",
                            metric = "Rsquared",
                            tuneGrid = expand.grid(C = ((i^2) * 0.007), sigma = seq(0.01,0.1,0.01)),
                            trControl = trainControl(method = "cv",
                                                     index = groupKFold(train_data$Genotype, k = 10), 
                                                     savePredictions = T, 
                                                     predictionBounds = c(T,T),
                                                     allowParallel = T
                                                     )
                            )
      ###############
      # Predictions #
      ###############
      lm_pred <- predict(spectra_lm, valid_data)
      pls_pred <- predict(spectra_pls, valid_data)
      pcr_pred <- predict(spectra_pcr, valid_data)
      rf_pred <- predict(spectra_rf, valid_data)
      svml_pred <- predict(spectra_svml, valid_data)
      svmp_pred <- predict(spectra_svmp, valid_data)
      svmr_pred <- predict(spectra_svmr, valid_data)
      
      proc.time() - ptm
      
      ###################
      # List of Results #
      ###################
      return(list(seed = n,
                  dataset = names(datasets_list[d]),
                  lm_spectra_cor = cor(lm_pred, valid_data$Moisture_Uptake, method = "spearman"),
                  pls_spectra_cor = cor(pls_pred, valid_data$Moisture_Uptake, method = "spearman"),
                  pcr_spectra_cor = cor(pcr_pred, valid_data$Moisture_Uptake, method = "spearman"),
                  rf_spectra_cor = cor(rf_pred, valid_data$Moisture_Uptake, method = "spearman"),
                  svml_spectra_cor = cor(svml_pred, valid_data$Moisture_Uptake, method = "spearman"),
                  svmp_spectra_cor = cor(svmp_pred, valid_data$Moisture_Uptake, method = "spearman"),
                  svmr_spectra_cor = cor(svmr_pred, valid_data$Moisture_Uptake, method = "spearman"),
                  lm_spectra_rmse = RMSE(lm_pred, valid_data$Moisture_Uptake),
                  pls_spectra_rmse = RMSE(pls_pred, valid_data$Moisture_Uptake),
                  pcr_spectra_rmse = RMSE(pcr_pred, valid_data$Moisture_Uptake),
                  rf_spectra_rmse = RMSE(rf_pred, valid_data$Moisture_Uptake),
                  svml_spectra_rmse = RMSE(svml_pred, valid_data$Moisture_Uptake),
                  svmp_spectra_rmse = RMSE(svmp_pred, valid_data$Moisture_Uptake),
                  svmr_spectra_rmse = RMSE(svmr_pred, valid_data$Moisture_Uptake),
                  lm_spectra_int = spectra_lm$bestTune[[1]],
                  pls_spectra_ncomp = spectra_pls$bestTune[[1]],
                  pcr_spectra_ncomp = spectra_pcr$bestTune[[1]],
                  rf_spectra_mtry = spectra_rf$bestTune[[1]],
                  svml_spectra_cost = spectra_svml$bestTune[[1]],
                  svmp_spectra_cost = spectra_svmp$bestTune[[3]],
                  svmp_spectra_scale = spectra_svmp$bestTune[[2]],
                  svmp_spectra_degree = spectra_svmp$bestTune[[1]],
                  svmr_spectra_cost = spectra_svmr$bestTune[[2]],
                  svmr_spectra_sigma = spectra_svmr$bestTune[[1]]))
    }
    stopImplicitCluster()
    
    spectra_random_results_tibble <- spectra_random_results %>%
      as_tibble() %>%
      unnest(cols = everything())

    Spectra_Random_Results <- bind_rows(Spectra_Random_Results, spectra_random_results_tibble)
  }
}

Spectra_Random_Results %>% 
  pivot_longer(cols = -c(seed, 
                         dataset,             # this section will "tidy" the data to make it more manageable
                         contains("ncomp"), 
                         contains("deg"), 
                         contains("scale"), 
                         contains("sigma"),
                         contains("_int"),
                         contains("_cost"),
                         contains("mtry")
  ),
  names_to = "Mod_Data_Met", 
  values_to = "Performance") %>%
  mutate(Model = case_when(str_detect(string = Mod_Data_Met, pattern = "lm") ~ "LM", # this section creates new columns for the different models
                           str_detect(string = Mod_Data_Met, pattern = "pcr") ~ "PCR",
                           str_detect(string = Mod_Data_Met, pattern = "rf") ~ "RF",
                           str_detect(string = Mod_Data_Met, pattern = "svml") ~ "SVML",
                           str_detect(string = Mod_Data_Met, pattern = "svmp") ~ "SVMP",
                           str_detect(string = Mod_Data_Met, pattern = "svmr") ~ "SVMR",
                           str_detect(string = Mod_Data_Met, pattern = "pls") ~ "PLS"),
         Metric = case_when(str_detect(string = Mod_Data_Met, pattern = "rmse") ~ "RMSE", # this section creates new columns for the different metrics
                            str_detect(string = Mod_Data_Met, pattern = "cor") ~ "Spearman_R"),
         Data = case_when(str_detect(string = Mod_Data_Met, pattern = "spectra") ~ "Spectra")) %>%
  select(-Mod_Data_Met) %>% # getting rid of the messy colum with mulitple pieces of information since they each have their own column  now
  mutate(Intercept = case_when(Model == "LM" ~ as.double(lm_spectra_int)),
         Mtry = case_when(Model == "RF" ~ as.double(rf_spectra_mtry)),
         C = case_when(Model == "SVML" ~ as.double(svml_spectra_cost),
                       Model == "SVMP" ~ as.double(svmp_spectra_cost),
                       Model == "SVMR" ~ as.double(svmr_spectra_cost),),
         Degree = case_when(Model == "SVMP" ~ as.double(svmp_spectra_degree),),
         Scale = case_when(Model == "SVMP" ~ as.double(svmp_spectra_scale),),
         Sigma = case_when(Model == "SVMR" ~ as.double(svmr_spectra_sigma)),
         Ncomp = case_when(Model == "PLS" ~ as.double(pls_spectra_ncomp) | Model == "PCR" ~ as.double(pcr_spectra_ncomp))) %>%
  select(seed, dataset, Model, Metric, Data, Intercept, Mtry, Ncomp, C, Degree, Scale, Sigma, Performance) %>% # getting rid of the extra columns we dont need
  group_by(Model, Metric, Data, Intercept, Mtry, Ncomp, C) %>%
  summarise(Mean_Degree = mean(Degree),
            Mean_Scale = mean(Scale),
            Mean_Sigma = mean(Sigma),
            Mean_Performance = mean(Performance),
            Std_Dev_Performance = sd(Performance),
            Range_Performance = max(Performance) - min(Performance),
            N = n(),
            Splitting = "Random") %>%
  arrange(desc(Mean_Performance)) %>% # arranging data from highest performance to lowest
  write_csv("../Results/spectra_randomsplit_results_clean.csv")

Spectra_Random_Results %>%
  write_csv("../Results/spectra_randomsplit_results_raw.csv")


