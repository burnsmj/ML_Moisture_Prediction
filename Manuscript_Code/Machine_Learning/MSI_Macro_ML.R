### MSI Macro Trait Machine Learning
### Michael Burns 
### 2020-07-07

# Doing machine learning on both macro and spectral traits at one time is a heafty project that can likely be extremely simplified by splitting them into
# their own tasks on MSI.

#############
# Libraries #
#############
library(caret)
library(pls)
library(spls)
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

################
# Loading Data #
################

dataset <- read_csv("../Data/Moisture_Uptake_Master_Dataset_Cook_Macro_Spectra.csv") %>%
  select(Sample_ID, Genotype, Protein_As_is, Fat_As_is, Fiber_As_is, Ash_As_is, Starch_As_is, Moisture_Uptake)

#################################
#################################
## Predictions with Boundaries ##
#################################
print("Predictions with Boundaries")
#################################

##########################################
# Splitting Data Randomly with Random CV #
##########################################

num_cores_limit <- 5
iterations <- 5

if (detectCores() > num_cores_limit) {
  num_cores <- num_cores_limit
} else {
  num_cores <- detectCores()
}

registerDoParallel(cores = num_cores)

Macro_Random_Results <- tibble()
print("Entering Random - Random Loop")
for(n in 1:100) {
  random_macro_results <- foreach(i = 1:iterations, .combine = rbind) %dopar% {
    ##################
    # Splitting Data #
    ##################
    training_index <- createDataPartition(dataset$Moisture_Uptake,
                                          times = 1,
                                          p = 0.8,
                                          list = FALSE)
    training_set <- dataset[training_index,]
    validation_set <- dataset[-training_index,]
    
    #############
    # Modelling #
    #############
    set.seed(n)
    macro_lm <- train(Moisture_Uptake ~ ., 
                      data = training_set %>%
                        select_if(is.numeric), 
                      method = "lm", 
                      metric = "Rsquared", 
                      tuneGrid = expand.grid(intercept = 1/i),
                      trControl = trainControl(method = "cv",
                                               #index = groupKFold(training_set$Genotype, k = 10), 
                                               savePredictions = T, 
                                               predictionBounds = c(T,T),
                                               allowParallel = T
                                               )
                      )   
    
    set.seed(n)
    macro_cart <- train(Moisture_Uptake ~ ., 
                        data = training_set %>%
                          select_if(is.numeric), 
                        method = "rpart", 
                        metric = "Rsquared", 
                        tuneGrid = expand.grid(cp = i*0.02),
                        trControl = trainControl(method = "cv",
                                                 #index = groupKFold(training_set$Genotype, k = 10), 
                                                 savePredictions = T, 
                                                 predictionBounds = c(T,T),
                                                 allowParallel = T
                                                 )
                        )
    
    set.seed(n)
    macro_knn <- train(Moisture_Uptake ~ ., 
                       data = training_set %>%
                         select_if(is.numeric), 
                       method = "knn", 
                       metric = "Rsquared", 
                       tuneGrid = expand.grid(k = i),
                       trControl = trainControl(method = "cv",
                                                #index = groupKFold(training_set$Genotype, k = 10), 
                                                savePredictions = T, 
                                                predictionBounds = c(T,T),
                                                allowParallel = T
                                                )
                        )
    
    set.seed(n)
    macro_svml <- train(Moisture_Uptake ~ ., 
                        data = training_set %>%
                          select_if(is.numeric), 
                        method = "svmLinear", 
                        metric = "Rsquared",
                        tuneGrid = expand.grid(C = i),
                        trControl = trainControl(method = "cv",
                                                 #index = groupKFold(training_set$Genotype, k = 10), 
                                                 savePredictions = T, 
                                                 predictionBounds = c(T,T),
                                                 allowParallel = T
                                                 )
                        )
    
    set.seed(n)
    macro_svmr <- train(Moisture_Uptake ~ ., 
                        data = training_set %>%
                          select_if(is.numeric), 
                        method = "svmRadial", 
                        metric = "Rsquared",
                        tuneGrid = expand.grid(C = i, sigma = c(0.001, 0.01, 0.1)),
                        trControl = trainControl(method = "cv",
                                                 #index = groupKFold(training_set$Genotype, k = 10), 
                                                 savePredictions = T, 
                                                 predictionBounds = c(T,T),
                                                 allowParallel = T
                                                 )
                        )
    
    set.seed(n)
    macro_svmp <- train(Moisture_Uptake ~ ., 
                        data = training_set %>%
                          select_if(is.numeric), 
                        method = "svmPoly", 
                        metric = "Rsquared", 
                        tuneGrid = expand.grid(C = i, scale = c(0.001, 0.01, 0.1), degree = c(1, 2, 3)),
                        trControl = trainControl(method = "cv",
                                                 #index = groupKFold(training_set$Genotype, k = 10), 
                                                 savePredictions = T, 
                                                 predictionBounds = c(T,T),
                                                 allowParallel = T
                                                 )
                        )
    
    set.seed(n)
    macro_rf <- train(Moisture_Uptake ~ ., 
                      data = training_set %>%
                        select_if(is.numeric), 
                      method = "rf", 
                      metric = "Rsquared", 
                      tuneGrid = expand.grid(mtry = i),
                      trControl = trainControl(method = "cv",
                                               #index = groupKFold(training_set$Genotype, k = 10), 
                                               savePredictions = T, 
                                               predictionBounds = c(T,T),
                                               allowParallel = T
                                               )
                      )
    
    ###############
    # Predictions #
    ###############
    lm_pred <- predict(macro_lm, validation_set)
    cart_pred <- predict(macro_cart, validation_set)
    knn_pred <- predict(macro_knn, validation_set)
    rf_pred <- predict(macro_rf, validation_set)
    svml_pred <- predict(macro_svml, validation_set)
    svmp_pred <- predict(macro_svmp, validation_set)
    svmr_pred <- predict(macro_svmr, validation_set)
    
    #####################
    # Returning Results #
    #####################
    return(list(lm_macro_random_cor = cor(lm_pred, validation_set$Moisture_Uptake, method = "spearman"),
                cart_macro_random_cor = cor(cart_pred, validation_set$Moisture_Uptake, method = "spearman"),
                knn_macro_random_cor = cor(knn_pred, validation_set$Moisture_Uptake, method = "spearman"),
                rf_macro_random_cor = cor(rf_pred, validation_set$Moisture_Uptake, method = "spearman"),
                svml_macro_random_cor = cor(svml_pred, validation_set$Moisture_Uptake, method = "spearman"),
                svmp_macro_random_cor = cor(svmp_pred, validation_set$Moisture_Uptake, method = "spearman"),
                svmr_macro_random_cor = cor(svmr_pred, validation_set$Moisture_Uptake, method = "spearman"),
                lm_macro_random_rmse = RMSE(lm_pred, validation_set$Moisture_Uptake),
                cart_macro_random_rmse = RMSE(cart_pred, validation_set$Moisture_Uptake),
                knn_macro_random_rmse = RMSE(knn_pred, validation_set$Moisture_Uptake),
                rf_macro_random_rmse = RMSE(rf_pred, validation_set$Moisture_Uptake),
                svml_macro_random_rmse = RMSE(svml_pred, validation_set$Moisture_Uptake),
                svmp_macro_random_rmse = RMSE(svmp_pred, validation_set$Moisture_Uptake),
                svmr_macro_random_rmse = RMSE(svmr_pred, validation_set$Moisture_Uptake),
                lm_macro_random_int = macro_lm$bestTune[[1]],
                cart_macro_random_cp = macro_cart$bestTune[[1]],
                knn_macro_random_k = macro_knn$bestTune[[1]],
                rf_macro_random_mtry = macro_rf$bestTune[[1]],
                svml_macro_random_cost = macro_svml$bestTune[[1]],
                svmp_macro_random_cost = macro_svmp$bestTune[[3]],
                svmp_macro_random_scale = macro_svmp$bestTune[[2]],
                svmp_macro_random_degree = macro_svmp$bestTune[[1]],
                svmr_macro_random_cost = macro_svmr$bestTune[[2]],
                svmr_macro_random_sigma = macro_svmr$bestTune[[1]],
                repetition = n
            ))
  }
  
  random_macro_results_tibble <- random_macro_results %>%
    as_tibble() %>%
    unnest(cols = everything())
  
  Macro_Random_Results <- bind_rows(Macro_Random_Results, random_macro_results_tibble)

}

########################################
# Stop Parallelization and Write Files #
########################################
stopImplicitCluster()

Macro_Random_Results %>% 
  pivot_longer(cols = -c(repetition,              # this section will "tidy" the data to make it more manageable
                         contains("ncomp"), 
                         contains("_cp"), 
                         contains("deg"), 
                         contains("scale"), 
                         contains("sigma"),
                         contains("_int"),
                         contains("_k"),
                         contains("_cost"),
                         contains("mtry")
                         ),
               names_to = "Mod_Data_Met", 
               values_to = "Performance") %>%
  mutate(Model = case_when(str_detect(string = Mod_Data_Met, pattern = "lm") ~ "LM", # this section creates new columns for the different models
                           str_detect(string = Mod_Data_Met, pattern = "knn") ~ "KNN",
                           str_detect(string = Mod_Data_Met, pattern = "rf") ~ "RF",
                           str_detect(string = Mod_Data_Met, pattern = "svml") ~ "SVML",
                           str_detect(string = Mod_Data_Met, pattern = "svmp") ~ "SVMP",
                           str_detect(string = Mod_Data_Met, pattern = "svmr") ~ "SVMR",
                           str_detect(string = Mod_Data_Met, pattern = "cart") ~ "CART"),
         Metric = case_when(str_detect(string = Mod_Data_Met, pattern = "rmse") ~ "RMSE", # this section creates new columns for the different metrics
                            str_detect(string = Mod_Data_Met, pattern = "cor") ~ "Spearman_R"),
         Data = case_when(str_detect(string = Mod_Data_Met, pattern = "macro") ~ "Macro")) %>%
  select(-Mod_Data_Met) %>% # getting rid of the messy colum with mulitple pieces of information since they each have their own column  now
  mutate(Intercept = case_when(Model == "LM" ~ as.double(lm_macro_random_int)),
         Mtry = case_when(Model == "RF" ~ as.double(rf_macro_random_mtry)),
         C = case_when(Model == "SVML" ~ as.double(svml_macro_random_cost),
                       Model == "SVMP" ~ as.double(svmp_macro_random_cost),
                       Model == "SVMR" ~ as.double(svmr_macro_random_cost),),
         Degree = case_when(Model == "SVMP" ~ as.double(svmp_macro_random_degree),),
         Scale = case_when(Model == "SVMP" ~ as.double(svmp_macro_random_scale),),
         Sigma = case_when(Model == "SVMR" ~ as.double(svmr_macro_random_sigma)),
         K = case_when(Model == "KNN" ~ as.double(knn_macro_random_k)),
         CP = case_when(Model == "CART" ~ as.double(cart_macro_random_cp))) %>%
  select(repetition, Model, Metric, Data, Intercept, Mtry, CP, K, C, Degree, Scale, Sigma, Performance) %>% # getting rid of the extra columns we dont need
  group_by(Model, Metric, Data, Intercept, Mtry, CP, K, C) %>%
  summarise(Mean_Degree = mean(Degree),
            Mean_Scale = mean(Scale),
            Mean_Sigma = mean(Sigma),
            Mean_Performance = mean(Performance),
            Std_Dev_Performance = sd(Performance),
            Range_Performance = max(Performance) - min(Performance),
            N = n(),
            Splitting = "Random") %>%
  arrange(desc(Mean_Performance)) %>% # arranging data from highest performance to lowest
  write_csv("../Results/macro_randomsplit_randomcv_clean.csv")
  
Macro_Random_Results %>%
  write_csv("../Results/macro_randomsplit_randomcv_raw.csv")

print('Randomly Split Random CV Models Finished')

#####################################################
# Splitting Data By Genotype with Genotype based CV #
#####################################################

num_cores_limit <- 5
iterations <- 5

if (detectCores() > num_cores_limit) {
  num_cores <- num_cores_limit
} else {
  num_cores <- detectCores()
}

registerDoParallel(cores = num_cores)

Macro_Genosplit_Genocv_Results <- tibble()
print("Entering Genotype - Genotype Loop")
for(n in 1:100) {
  genosplit_genocv_macro_results <- foreach(i = 1:iterations, .combine = rbind) %dopar% {
    ##################
    # Splitting Data #
    ##################
    training_genotypes <- dataset %>% 
      select(Genotype) %>% 
      unique() %>%
      sample_frac(size = 0.8) %>%
      pull()
    
    training_set <- dataset %>%
      filter(Genotype %in% training_genotypes)
    
    validation_set <- dataset %>%
      filter(!Genotype %in% training_genotypes)
    
    #############
    # Modelling #
    #############
    set.seed(n)
    macro_lm <- train(Moisture_Uptake ~ ., 
                      data = training_set %>%
                        select_if(is.numeric), 
                      method = "lm", 
                      metric = "Rsquared", 
                      tuneGrid = expand.grid(intercept = 1/i),
                      trControl = trainControl(method = "cv",
                                               index = groupKFold(training_set$Genotype, k = 10), 
                                               savePredictions = T, 
                                               predictionBounds = c(T,T),
                                               allowParallel = T
                      )
    )   
    
    set.seed(n)
    macro_cart <- train(Moisture_Uptake ~ ., 
                        data = training_set %>%
                          select_if(is.numeric), 
                        method = "rpart", 
                        metric = "Rsquared", 
                        tuneGrid = expand.grid(cp = i*0.02),
                        trControl = trainControl(method = "cv",
                                                 index = groupKFold(training_set$Genotype, k = 10), 
                                                 savePredictions = T, 
                                                 predictionBounds = c(T,T),
                                                 allowParallel = T
                        )
    )
    
    set.seed(n)
    macro_knn <- train(Moisture_Uptake ~ ., 
                       data = training_set %>%
                         select_if(is.numeric), 
                       method = "knn", 
                       metric = "Rsquared", 
                       tuneGrid = expand.grid(k = i),
                       trControl = trainControl(method = "cv",
                                                index = groupKFold(training_set$Genotype, k = 10), 
                                                savePredictions = T, 
                                                predictionBounds = c(T,T),
                                                allowParallel = T
                       )
    )
    
    set.seed(n)
    macro_svml <- train(Moisture_Uptake ~ ., 
                        data = training_set %>%
                          select_if(is.numeric), 
                        method = "svmLinear", 
                        metric = "Rsquared",
                        tuneGrid = expand.grid(C = i),
                        trControl = trainControl(method = "cv",
                                                 index = groupKFold(training_set$Genotype, k = 10), 
                                                 savePredictions = T, 
                                                 predictionBounds = c(T,T),
                                                 allowParallel = T
                        )
    )
    
    set.seed(n)
    macro_svmr <- train(Moisture_Uptake ~ ., 
                        data = training_set %>%
                          select_if(is.numeric), 
                        method = "svmRadial", 
                        metric = "Rsquared",
                        tuneGrid = expand.grid(C = i, sigma = c(0.001, 0.01, 0.1)),
                        trControl = trainControl(method = "cv",
                                                 index = groupKFold(training_set$Genotype, k = 10), 
                                                 savePredictions = T, 
                                                 predictionBounds = c(T,T),
                                                 allowParallel = T
                        )
    )
    
    set.seed(n)
    macro_svmp <- train(Moisture_Uptake ~ ., 
                        data = training_set %>%
                          select_if(is.numeric), 
                        method = "svmPoly", 
                        metric = "Rsquared", 
                        tuneGrid = expand.grid(C = i, scale = c(0.001, 0.01, 0.1), degree = c(1, 2, 3)),
                        trControl = trainControl(method = "cv",
                                                 index = groupKFold(training_set$Genotype, k = 10), 
                                                 savePredictions = T, 
                                                 predictionBounds = c(T,T),
                                                 allowParallel = T
                        )
    )
    
    set.seed(n)
    macro_rf <- train(Moisture_Uptake ~ ., 
                      data = training_set %>%
                        select_if(is.numeric), 
                      method = "rf", 
                      metric = "Rsquared", 
                      tuneGrid = expand.grid(mtry = i),
                      trControl = trainControl(method = "cv",
                                               index = groupKFold(training_set$Genotype, k = 10), 
                                               savePredictions = T, 
                                               predictionBounds = c(T,T),
                                               allowParallel = T
                      )
    )
    
    ###############
    # Predictions #
    ###############
    lm_pred <- predict(macro_lm, validation_set)
    cart_pred <- predict(macro_cart, validation_set)
    knn_pred <- predict(macro_knn, validation_set)
    rf_pred <- predict(macro_rf, validation_set)
    svml_pred <- predict(macro_svml, validation_set)
    svmp_pred <- predict(macro_svmp, validation_set)
    svmr_pred <- predict(macro_svmr, validation_set)
    
    #####################
    # Returning Results #
    #####################
    return(list(lm_macro_random_cor = cor(lm_pred, validation_set$Moisture_Uptake, method = "spearman"),
                cart_macro_random_cor = cor(cart_pred, validation_set$Moisture_Uptake, method = "spearman"),
                knn_macro_random_cor = cor(knn_pred, validation_set$Moisture_Uptake, method = "spearman"),
                rf_macro_random_cor = cor(rf_pred, validation_set$Moisture_Uptake, method = "spearman"),
                svml_macro_random_cor = cor(svml_pred, validation_set$Moisture_Uptake, method = "spearman"),
                svmp_macro_random_cor = cor(svmp_pred, validation_set$Moisture_Uptake, method = "spearman"),
                svmr_macro_random_cor = cor(svmr_pred, validation_set$Moisture_Uptake, method = "spearman"),
                lm_macro_random_rmse = RMSE(lm_pred, validation_set$Moisture_Uptake),
                cart_macro_random_rmse = RMSE(cart_pred, validation_set$Moisture_Uptake),
                knn_macro_random_rmse = RMSE(knn_pred, validation_set$Moisture_Uptake),
                rf_macro_random_rmse = RMSE(rf_pred, validation_set$Moisture_Uptake),
                svml_macro_random_rmse = RMSE(svml_pred, validation_set$Moisture_Uptake),
                svmp_macro_random_rmse = RMSE(svmp_pred, validation_set$Moisture_Uptake),
                svmr_macro_random_rmse = RMSE(svmr_pred, validation_set$Moisture_Uptake),
                lm_macro_random_int = macro_lm$bestTune[[1]],
                cart_macro_random_cp = macro_cart$bestTune[[1]],
                knn_macro_random_k = macro_knn$bestTune[[1]],
                rf_macro_random_mtry = macro_rf$bestTune[[1]],
                svml_macro_random_cost = macro_svml$bestTune[[1]],
                svmp_macro_random_cost = macro_svmp$bestTune[[3]],
                svmp_macro_random_scale = macro_svmp$bestTune[[2]],
                svmp_macro_random_degree = macro_svmp$bestTune[[1]],
                svmr_macro_random_cost = macro_svmr$bestTune[[2]],
                svmr_macro_random_sigma = macro_svmr$bestTune[[1]],
                repetition = n
    ))
  }
  
  genosplit_genocv_macro_results_tibble <- genosplit_genocv_macro_results %>%
    as_tibble() %>%
    unnest(cols = everything())
  
  Macro_Genosplit_Genocv_Results <- bind_rows(Macro_Genosplit_Genocv_Results, genosplit_genocv_macro_results_tibble)
  
}

########################################
# Stop Parallelization and Write Files #
########################################
stopImplicitCluster()

Macro_Genosplit_Genocv_Results %>% 
  pivot_longer(cols = -c(repetition,              # this section will "tidy" the data to make it more manageable
                         contains("ncomp"), 
                         contains("_cp"), 
                         contains("deg"), 
                         contains("scale"), 
                         contains("sigma"),
                         contains("_int"),
                         contains("_k"),
                         contains("_cost"),
                         contains("mtry")
  ),
  names_to = "Mod_Data_Met", 
  values_to = "Performance") %>%
  mutate(Model = case_when(str_detect(string = Mod_Data_Met, pattern = "lm") ~ "LM", # this section creates new columns for the different models
                           str_detect(string = Mod_Data_Met, pattern = "knn") ~ "KNN",
                           str_detect(string = Mod_Data_Met, pattern = "rf") ~ "RF",
                           str_detect(string = Mod_Data_Met, pattern = "svml") ~ "SVML",
                           str_detect(string = Mod_Data_Met, pattern = "svmp") ~ "SVMP",
                           str_detect(string = Mod_Data_Met, pattern = "svmr") ~ "SVMR",
                           str_detect(string = Mod_Data_Met, pattern = "cart") ~ "CART"),
         Metric = case_when(str_detect(string = Mod_Data_Met, pattern = "rmse") ~ "RMSE", # this section creates new columns for the different metrics
                            str_detect(string = Mod_Data_Met, pattern = "cor") ~ "Spearman_R"),
         Data = case_when(str_detect(string = Mod_Data_Met, pattern = "macro") ~ "Macro")) %>%
  select(-Mod_Data_Met) %>% # getting rid of the messy colum with mulitple pieces of information since they each have their own column  now
  mutate(Intercept = case_when(Model == "LM" ~ as.double(lm_macro_random_int)),
         Mtry = case_when(Model == "RF" ~ as.double(rf_macro_random_mtry)),
         C = case_when(Model == "SVML" ~ as.double(svml_macro_random_cost),
                       Model == "SVMP" ~ as.double(svmp_macro_random_cost),
                       Model == "SVMR" ~ as.double(svmr_macro_random_cost),),
         Degree = case_when(Model == "SVMP" ~ as.double(svmp_macro_random_degree),),
         Scale = case_when(Model == "SVMP" ~ as.double(svmp_macro_random_scale),),
         Sigma = case_when(Model == "SVMR" ~ as.double(svmr_macro_random_sigma)),
         K = case_when(Model == "KNN" ~ as.double(knn_macro_random_k)),
         CP = case_when(Model == "CART" ~ as.double(cart_macro_random_cp))) %>%
  select(repetition, Model, Metric, Data, Intercept, Mtry, CP, K, C, Degree, Scale, Sigma, Performance) %>% # getting rid of the extra columns we dont need
  group_by(Model, Metric, Data, Intercept, Mtry, CP, K, C) %>%
  summarise(Mean_Degree = mean(Degree),
            Mean_Scale = mean(Scale),
            Mean_Sigma = mean(Sigma),
            Mean_Performance = mean(Performance),
            Std_Dev_Performance = sd(Performance),
            Range_Performance = max(Performance) - min(Performance),
            N = n(),
            Splitting = "Genotype") %>%
  arrange(desc(Mean_Performance)) %>% # arranging data from highest performance to lowest
  write_csv("../Results/macro_genosplit_genocv_clean.csv")

Macro_Genosplit_Genocv_Results %>%
  write_csv("../Results/macro_genosplit_genocv_raw.csv")

print('Genotype Split Genotype CV Models Finished')

#############################################
# Splitting Data By Genotype with Random CV #
#############################################

num_cores_limit <- 5
iterations <- 5

if (detectCores() > num_cores_limit) {
  num_cores <- num_cores_limit
} else {
  num_cores <- detectCores()
}

registerDoParallel(cores = num_cores)

Macro_Genosplit_Randomcv_Results <- tibble()
print("Entering Genotype - Random Loop")
for(n in 1:100) {
  genosplit_randomcv_macro_results <- foreach(i = 1:iterations, .combine = rbind) %dopar% {
    ##################
    # Splitting Data #
    ##################
    training_genotypes <- dataset %>% 
      select(Genotype) %>% 
      unique() %>%
      sample_frac(size = 0.8) %>%
      pull()
    
    training_set <- dataset %>%
      filter(Genotype %in% training_genotypes)
    
    validation_set <- dataset %>%
      filter(!Genotype %in% training_genotypes)
    
    #############
    # Modelling #
    #############
    set.seed(n)
    macro_lm <- train(Moisture_Uptake ~ ., 
                      data = training_set %>%
                        select_if(is.numeric), 
                      method = "lm", 
                      metric = "Rsquared", 
                      tuneGrid = expand.grid(intercept = 1/i),
                      trControl = trainControl(method = "cv",
                                               #index = groupKFold(training_set$Genotype, k = 10), 
                                               savePredictions = T, 
                                               predictionBounds = c(T,T),
                                               allowParallel = T
                      )
    )   
    
    set.seed(n)
    macro_cart <- train(Moisture_Uptake ~ ., 
                        data = training_set %>%
                          select_if(is.numeric), 
                        method = "rpart", 
                        metric = "Rsquared", 
                        tuneGrid = expand.grid(cp = i*0.02),
                        trControl = trainControl(method = "cv",
                                                 #index = groupKFold(training_set$Genotype, k = 10), 
                                                 savePredictions = T, 
                                                 predictionBounds = c(T,T),
                                                 allowParallel = T
                        )
    )
    
    set.seed(n)
    macro_knn <- train(Moisture_Uptake ~ ., 
                       data = training_set %>%
                         select_if(is.numeric), 
                       method = "knn", 
                       metric = "Rsquared", 
                       tuneGrid = expand.grid(k = i),
                       trControl = trainControl(method = "cv",
                                                #index = groupKFold(training_set$Genotype, k = 10), 
                                                savePredictions = T, 
                                                predictionBounds = c(T,T),
                                                allowParallel = T
                       )
    )
    
    set.seed(n)
    macro_rf <- train(Moisture_Uptake ~ ., 
                      data = training_set %>%
                        select_if(is.numeric), 
                      method = "rf", 
                      metric = "Rsquared", 
                      tuneGrid = expand.grid(mtry = i),
                      trControl = trainControl(method = "cv",
                                               #index = groupKFold(training_set$Genotype, k = 10), 
                                               savePredictions = T, 
                                               predictionBounds = c(T,T),
                                               allowParallel = T
                      )
    )
    
    set.seed(n)
    macro_svml <- train(Moisture_Uptake ~ ., 
                        data = training_set %>%
                          select_if(is.numeric), 
                        method = "svmLinear", 
                        metric = "Rsquared",
                        tuneGrid = expand.grid(C = i),
                        trControl = trainControl(method = "cv",
                                                 #index = groupKFold(training_set$Genotype, k = 10), 
                                                 savePredictions = T, 
                                                 predictionBounds = c(T,T),
                                                 allowParallel = T
                        )
    )
    
    set.seed(n)
    macro_svmr <- train(Moisture_Uptake ~ ., 
                        data = training_set %>%
                          select_if(is.numeric), 
                        method = "svmRadial", 
                        metric = "Rsquared",
                        tuneGrid = expand.grid(C = i, sigma = c(0.001, 0.01, 0.1)),
                        trControl = trainControl(method = "cv",
                                                 #index = groupKFold(training_set$Genotype, k = 10), 
                                                 savePredictions = T, 
                                                 predictionBounds = c(T,T),
                                                 allowParallel = T
                        )
    )
    
    set.seed(n)
    macro_svmp <- train(Moisture_Uptake ~ ., 
                        data = training_set %>%
                          select_if(is.numeric), 
                        method = "svmPoly", 
                        metric = "Rsquared", 
                        tuneGrid = expand.grid(C = i, scale = c(0.001, 0.01, 0.1), degree = c(1, 2, 3)),
                        trControl = trainControl(method = "cv",
                                                 #index = groupKFold(training_set$Genotype, k = 10), 
                                                 savePredictions = T, 
                                                 predictionBounds = c(T,T),
                                                 allowParallel = T
                        )
    )
    
    ###############
    # Predictions #
    ###############
    lm_pred <- predict(macro_lm, validation_set)
    cart_pred <- predict(macro_cart, validation_set)
    knn_pred <- predict(macro_knn, validation_set)
    rf_pred <- predict(macro_rf, validation_set)
    svml_pred <- predict(macro_svml, validation_set)
    svmp_pred <- predict(macro_svmp, validation_set)
    svmr_pred <- predict(macro_svmr, validation_set)
    
    #####################
    # Returning Results #
    #####################
    return(list(lm_macro_random_cor = cor(lm_pred, validation_set$Moisture_Uptake, method = "spearman"),
                cart_macro_random_cor = cor(cart_pred, validation_set$Moisture_Uptake, method = "spearman"),
                knn_macro_random_cor = cor(knn_pred, validation_set$Moisture_Uptake, method = "spearman"),
                rf_macro_random_cor = cor(rf_pred, validation_set$Moisture_Uptake, method = "spearman"),
                svml_macro_random_cor = cor(svml_pred, validation_set$Moisture_Uptake, method = "spearman"),
                svmp_macro_random_cor = cor(svmp_pred, validation_set$Moisture_Uptake, method = "spearman"),
                svmr_macro_random_cor = cor(svmr_pred, validation_set$Moisture_Uptake, method = "spearman"),
                lm_macro_random_rmse = RMSE(lm_pred, validation_set$Moisture_Uptake),
                cart_macro_random_rmse = RMSE(cart_pred, validation_set$Moisture_Uptake),
                knn_macro_random_rmse = RMSE(knn_pred, validation_set$Moisture_Uptake),
                rf_macro_random_rmse = RMSE(rf_pred, validation_set$Moisture_Uptake),
                svml_macro_random_rmse = RMSE(svml_pred, validation_set$Moisture_Uptake),
                svmp_macro_random_rmse = RMSE(svmp_pred, validation_set$Moisture_Uptake),
                svmr_macro_random_rmse = RMSE(svmr_pred, validation_set$Moisture_Uptake),
                lm_macro_random_int = macro_lm$bestTune[[1]],
                cart_macro_random_cp = macro_cart$bestTune[[1]],
                knn_macro_random_k = macro_knn$bestTune[[1]],
                rf_macro_random_mtry = macro_rf$bestTune[[1]],
                svml_macro_random_cost = macro_svml$bestTune[[1]],
                svmp_macro_random_cost = macro_svmp$bestTune[[3]],
                svmp_macro_random_scale = macro_svmp$bestTune[[2]],
                svmp_macro_random_degree = macro_svmp$bestTune[[1]],
                svmr_macro_random_cost = macro_svmr$bestTune[[2]],
                svmr_macro_random_sigma = macro_svmr$bestTune[[1]],
                repetition = n
    ))
  }
  
  genosplit_randomcv_macro_results_tibble <- genosplit_randomcv_macro_results %>%
    as_tibble() %>%
    unnest(cols = everything())
  
  Macro_Genosplit_Randomcv_Results <- bind_rows(Macro_Genosplit_Randomcv_Results, genosplit_randomcv_macro_results_tibble)
  
}

########################################
# Stop Parallelization and Write Files #
########################################
stopImplicitCluster()

Macro_Genosplit_Randomcv_Results %>% 
  pivot_longer(cols = -c(repetition,              # this section will "tidy" the data to make it more manageable
                         contains("ncomp"), 
                         contains("_cp"), 
                         contains("_int"),
                         contains("_k"),
                         contains("mtry")
  ),
  names_to = "Mod_Data_Met", 
  values_to = "Performance") %>%
  mutate(Model = case_when(str_detect(string = Mod_Data_Met, pattern = "lm") ~ "LM", # this section creates new columns for the different models
                           str_detect(string = Mod_Data_Met, pattern = "knn") ~ "KNN",
                           str_detect(string = Mod_Data_Met, pattern = "rf") ~ "RF",
                           str_detect(string = Mod_Data_Met, pattern = "cart") ~ "CART"),
         Metric = case_when(str_detect(string = Mod_Data_Met, pattern = "rmse") ~ "RMSE", # this section creates new columns for the different metrics
                            str_detect(string = Mod_Data_Met, pattern = "cor") ~ "Spearman_R"),
         Data = case_when(str_detect(string = Mod_Data_Met, pattern = "macro") ~ "Macro")) %>%
  select(-Mod_Data_Met) %>% # getting rid of the messy colum with mulitple pieces of information since they each have their own column  now
  mutate(Intercept = case_when(Model == "LM" ~ as.double(lm_macro_random_int)),
         Mtry = case_when(Model == "RF" ~ as.double(rf_macro_random_mtry)),
         K = case_when(Model == "KNN" ~ as.double(knn_macro_random_k)),
         CP = case_when(Model == "CART" ~ as.double(cart_macro_random_cp))) %>%
  select(repetition, Model, Metric, Data, Intercept, Mtry, CP, K, Performance) %>% # getting rid of the extra columns we dont need
  group_by(Model, Metric, Data, Intercept, Mtry, CP, K) %>%
  summarise(Mean_Performance = mean(Performance),
            Std_Dev_Performance = sd(Performance),
            Range_Performance = max(Performance) - min(Performance),
            N = n(),
            Splitting = "Environment") %>%
  arrange(desc(Mean_Performance)) %>% # arranging data from highest performance to lowest
  write_csv("../Results/macro_genosplit_randomcv_clean.csv")

Macro_Genosplit_Randomcv_Results %>%
  write_csv("../Results/macro_genosplit_randomcv_raw.csv")

print('Genotype Split Random CV Models Finished')


####################################
####################################
## Predictions without Boundaries ##
####################################
####################################
print("Predictions without Boundaries")
##########################################
# Splitting Data Randomly with Random CV #
##########################################

num_cores_limit <- 5
iterations <- 5

if (detectCores() > num_cores_limit) {
  num_cores <- num_cores_limit
} else {
  num_cores <- detectCores()
}

registerDoParallel(cores = num_cores)

Macro_Random_Results <- tibble()
print("Entering Random - Random Loop")
for(n in 1:100) {
  random_macro_results <- foreach(i = 1:iterations, .combine = rbind) %dopar% {
    ##################
    # Splitting Data #
    ##################
    training_index <- createDataPartition(dataset$Moisture_Uptake,
                                          times = 1,
                                          p = 0.8,
                                          list = FALSE)
    training_set <- dataset[training_index,]
    validation_set <- dataset[-training_index,]
    
    #############
    # Modelling #
    #############
    set.seed(n)
    macro_lm <- train(Moisture_Uptake ~ ., 
                      data = training_set %>%
                        select_if(is.numeric), 
                      method = "lm", 
                      metric = "Rsquared", 
                      tuneGrid = expand.grid(intercept = 1/i),
                      trControl = trainControl(method = "cv",
                                               #index = groupKFold(training_set$Genotype, k = 10), 
                                               savePredictions = T, 
                                               predictionBounds = c(T,T),
                                               allowParallel = T
                      )
    )   
    
    set.seed(n)
    macro_cart <- train(Moisture_Uptake ~ ., 
                        data = training_set %>%
                          select_if(is.numeric), 
                        method = "rpart", 
                        metric = "Rsquared", 
                        tuneGrid = expand.grid(cp = i*0.02),
                        trControl = trainControl(method = "cv",
                                                 #index = groupKFold(training_set$Genotype, k = 10), 
                                                 savePredictions = T, 
                                                 predictionBounds = c(T,T),
                                                 allowParallel = T
                        )
    )
    
    set.seed(n)
    macro_knn <- train(Moisture_Uptake ~ ., 
                       data = training_set %>%
                         select_if(is.numeric), 
                       method = "knn", 
                       metric = "Rsquared", 
                       tuneGrid = expand.grid(k = i),
                       trControl = trainControl(method = "cv",
                                                #index = groupKFold(training_set$Genotype, k = 10), 
                                                savePredictions = T, 
                                                predictionBounds = c(T,T),
                                                allowParallel = T
                       )
    )
    
    set.seed(n)
    macro_svml <- train(Moisture_Uptake ~ ., 
                        data = training_set %>%
                          select_if(is.numeric), 
                        method = "svmLinear", 
                        metric = "Rsquared",
                        tuneGrid = expand.grid(C = i),
                        trControl = trainControl(method = "cv",
                                                 #index = groupKFold(training_set$Genotype, k = 10), 
                                                 savePredictions = T, 
                                                 predictionBounds = c(T,T),
                                                 allowParallel = T
                        )
    )
    
    set.seed(n)
    macro_svmr <- train(Moisture_Uptake ~ ., 
                        data = training_set %>%
                          select_if(is.numeric), 
                        method = "svmRadial", 
                        metric = "Rsquared",
                        tuneGrid = expand.grid(C = i, sigma = c(0.001, 0.01, 0.1)),
                        trControl = trainControl(method = "cv",
                                                 #index = groupKFold(training_set$Genotype, k = 10), 
                                                 savePredictions = T, 
                                                 predictionBounds = c(T,T),
                                                 allowParallel = T
                        )
    )
    
    set.seed(n)
    macro_svmp <- train(Moisture_Uptake ~ ., 
                        data = training_set %>%
                          select_if(is.numeric), 
                        method = "svmPoly", 
                        metric = "Rsquared", 
                        tuneGrid = expand.grid(C = i, scale = c(0.001, 0.01, 0.1), degree = c(1, 2, 3)),
                        trControl = trainControl(method = "cv",
                                                 #index = groupKFold(training_set$Genotype, k = 10), 
                                                 savePredictions = T, 
                                                 predictionBounds = c(T,T),
                                                 allowParallel = T
                        )
    )
    
    set.seed(n)
    macro_rf <- train(Moisture_Uptake ~ ., 
                      data = training_set %>%
                        select_if(is.numeric), 
                      method = "rf", 
                      metric = "Rsquared", 
                      tuneGrid = expand.grid(mtry = i),
                      trControl = trainControl(method = "cv",
                                               #index = groupKFold(training_set$Genotype, k = 10), 
                                               savePredictions = T, 
                                               predictionBounds = c(T,T),
                                               allowParallel = T
                      )
    )
    
    ###############
    # Predictions #
    ###############
    lm_pred <- predict(macro_lm, validation_set)
    cart_pred <- predict(macro_cart, validation_set)
    knn_pred <- predict(macro_knn, validation_set)
    rf_pred <- predict(macro_rf, validation_set)
    svml_pred <- predict(macro_svml, validation_set)
    svmp_pred <- predict(macro_svmp, validation_set)
    svmr_pred <- predict(macro_svmr, validation_set)
    
    #####################
    # Returning Results #
    #####################
    return(list(lm_macro_random_cor = cor(lm_pred, validation_set$Moisture_Uptake, method = "spearman"),
                cart_macro_random_cor = cor(cart_pred, validation_set$Moisture_Uptake, method = "spearman"),
                knn_macro_random_cor = cor(knn_pred, validation_set$Moisture_Uptake, method = "spearman"),
                rf_macro_random_cor = cor(rf_pred, validation_set$Moisture_Uptake, method = "spearman"),
                svml_macro_random_cor = cor(svml_pred, validation_set$Moisture_Uptake, method = "spearman"),
                svmp_macro_random_cor = cor(svmp_pred, validation_set$Moisture_Uptake, method = "spearman"),
                svmr_macro_random_cor = cor(svmr_pred, validation_set$Moisture_Uptake, method = "spearman"),
                lm_macro_random_rmse = RMSE(lm_pred, validation_set$Moisture_Uptake),
                cart_macro_random_rmse = RMSE(cart_pred, validation_set$Moisture_Uptake),
                knn_macro_random_rmse = RMSE(knn_pred, validation_set$Moisture_Uptake),
                rf_macro_random_rmse = RMSE(rf_pred, validation_set$Moisture_Uptake),
                svml_macro_random_rmse = RMSE(svml_pred, validation_set$Moisture_Uptake),
                svmp_macro_random_rmse = RMSE(svmp_pred, validation_set$Moisture_Uptake),
                svmr_macro_random_rmse = RMSE(svmr_pred, validation_set$Moisture_Uptake),
                lm_macro_random_int = macro_lm$bestTune[[1]],
                cart_macro_random_cp = macro_cart$bestTune[[1]],
                knn_macro_random_k = macro_knn$bestTune[[1]],
                rf_macro_random_mtry = macro_rf$bestTune[[1]],
                svml_macro_random_cost = macro_svml$bestTune[[1]],
                svmp_macro_random_cost = macro_svmp$bestTune[[3]],
                svmp_macro_random_scale = macro_svmp$bestTune[[2]],
                svmp_macro_random_degree = macro_svmp$bestTune[[1]],
                svmr_macro_random_cost = macro_svmr$bestTune[[2]],
                svmr_macro_random_sigma = macro_svmr$bestTune[[1]],
                repetition = n
    ))
  }
  
  random_macro_results_tibble <- random_macro_results %>%
    as_tibble() %>%
    unnest(cols = everything())
  
  Macro_Random_Results <- bind_rows(Macro_Random_Results, random_macro_results_tibble)
  
}

########################################
# Stop Parallelization and Write Files #
########################################
stopImplicitCluster()

Macro_Random_Results %>% 
  pivot_longer(cols = -c(repetition,              # this section will "tidy" the data to make it more manageable
                         contains("ncomp"), 
                         contains("_cp"), 
                         contains("deg"), 
                         contains("scale"), 
                         contains("sigma"),
                         contains("_int"),
                         contains("_k"),
                         contains("_cost"),
                         contains("mtry")
  ),
  names_to = "Mod_Data_Met", 
  values_to = "Performance") %>%
  mutate(Model = case_when(str_detect(string = Mod_Data_Met, pattern = "lm") ~ "LM", # this section creates new columns for the different models
                           str_detect(string = Mod_Data_Met, pattern = "knn") ~ "KNN",
                           str_detect(string = Mod_Data_Met, pattern = "rf") ~ "RF",
                           str_detect(string = Mod_Data_Met, pattern = "svml") ~ "SVML",
                           str_detect(string = Mod_Data_Met, pattern = "svmp") ~ "SVMP",
                           str_detect(string = Mod_Data_Met, pattern = "svmr") ~ "SVMR",
                           str_detect(string = Mod_Data_Met, pattern = "cart") ~ "CART"),
         Metric = case_when(str_detect(string = Mod_Data_Met, pattern = "rmse") ~ "RMSE", # this section creates new columns for the different metrics
                            str_detect(string = Mod_Data_Met, pattern = "cor") ~ "Spearman_R"),
         Data = case_when(str_detect(string = Mod_Data_Met, pattern = "macro") ~ "Macro")) %>%
  select(-Mod_Data_Met) %>% # getting rid of the messy colum with mulitple pieces of information since they each have their own column  now
  mutate(Intercept = case_when(Model == "LM" ~ as.double(lm_macro_random_int)),
         Mtry = case_when(Model == "RF" ~ as.double(rf_macro_random_mtry)),
         C = case_when(Model == "SVML" ~ as.double(svml_macro_random_cost),
                       Model == "SVMP" ~ as.double(svmp_macro_random_cost),
                       Model == "SVMR" ~ as.double(svmr_macro_random_cost),),
         Degree = case_when(Model == "SVMP" ~ as.double(svmp_macro_random_degree),),
         Scale = case_when(Model == "SVMP" ~ as.double(svmp_macro_random_scale),),
         Sigma = case_when(Model == "SVMR" ~ as.double(svmr_macro_random_sigma)),
         K = case_when(Model == "KNN" ~ as.double(knn_macro_random_k)),
         CP = case_when(Model == "CART" ~ as.double(cart_macro_random_cp))) %>%
  select(repetition, Model, Metric, Data, Intercept, Mtry, CP, K, C, Degree, Scale, Sigma, Performance) %>% # getting rid of the extra columns we dont need
  group_by(Model, Metric, Data, Intercept, Mtry, CP, K, C) %>%
  summarise(Mean_Degree = mean(Degree),
            Mean_Scale = mean(Scale),
            Mean_Sigma = mean(Sigma),
            Mean_Performance = mean(Performance),
            Std_Dev_Performance = sd(Performance),
            Range_Performance = max(Performance) - min(Performance),
            N = n(),
            Splitting = "Random") %>%
  arrange(desc(Mean_Performance)) %>% # arranging data from highest performance to lowest
  write_csv("../Results/macro_randomsplit_randomcv_nobound_clean.csv")

Macro_Random_Results %>%
  write_csv("../Results/macro_randomsplit_randomcv_nobound_raw.csv")

print('No Bound Randomly Split Random CV Models Finished')

#####################################################
# Splitting Data By Genotype with Genotype based CV #
#####################################################

num_cores_limit <- 5
iterations <- 5

if (detectCores() > num_cores_limit) {
  num_cores <- num_cores_limit
} else {
  num_cores <- detectCores()
}

registerDoParallel(cores = num_cores)

Macro_Genosplit_Genocv_Results <- tibble()
print("Entering Genotype - Genotype Loop")
for(n in 1:100) {
  genosplit_genocv_macro_results <- foreach(i = 1:iterations, .combine = rbind) %dopar% {
    ##################
    # Splitting Data #
    ##################
    training_genotypes <- dataset %>% 
      select(Genotype) %>% 
      unique() %>%
      sample_frac(size = 0.8) %>%
      pull()
    
    training_set <- dataset %>%
      filter(Genotype %in% training_genotypes)
    
    validation_set <- dataset %>%
      filter(!Genotype %in% training_genotypes)
    
    #############
    # Modelling #
    #############
    set.seed(n)
    macro_lm <- train(Moisture_Uptake ~ ., 
                      data = training_set %>%
                        select_if(is.numeric), 
                      method = "lm", 
                      metric = "Rsquared", 
                      tuneGrid = expand.grid(intercept = 1/i),
                      trControl = trainControl(method = "cv",
                                               index = groupKFold(training_set$Genotype, k = 10), 
                                               savePredictions = T, 
                                               predictionBounds = c(T,T),
                                               allowParallel = T
                      )
    )   
    
    set.seed(n)
    macro_cart <- train(Moisture_Uptake ~ ., 
                        data = training_set %>%
                          select_if(is.numeric), 
                        method = "rpart", 
                        metric = "Rsquared", 
                        tuneGrid = expand.grid(cp = i*0.02),
                        trControl = trainControl(method = "cv",
                                                 index = groupKFold(training_set$Genotype, k = 10), 
                                                 savePredictions = T, 
                                                 predictionBounds = c(T,T),
                                                 allowParallel = T
                        )
    )
    
    set.seed(n)
    macro_knn <- train(Moisture_Uptake ~ ., 
                       data = training_set %>%
                         select_if(is.numeric), 
                       method = "knn", 
                       metric = "Rsquared", 
                       tuneGrid = expand.grid(k = i),
                       trControl = trainControl(method = "cv",
                                                index = groupKFold(training_set$Genotype, k = 10), 
                                                savePredictions = T, 
                                                predictionBounds = c(T,T),
                                                allowParallel = T
                       )
    )
    
    set.seed(n)
    macro_svml <- train(Moisture_Uptake ~ ., 
                        data = training_set %>%
                          select_if(is.numeric), 
                        method = "svmLinear", 
                        metric = "Rsquared",
                        tuneGrid = expand.grid(C = i),
                        trControl = trainControl(method = "cv",
                                                 index = groupKFold(training_set$Genotype, k = 10), 
                                                 savePredictions = T, 
                                                 predictionBounds = c(T,T),
                                                 allowParallel = T
                        )
    )
    
    set.seed(n)
    macro_svmr <- train(Moisture_Uptake ~ ., 
                        data = training_set %>%
                          select_if(is.numeric), 
                        method = "svmRadial", 
                        metric = "Rsquared",
                        tuneGrid = expand.grid(C = i, sigma = c(0.001, 0.01, 0.1)),
                        trControl = trainControl(method = "cv",
                                                 index = groupKFold(training_set$Genotype, k = 10), 
                                                 savePredictions = T, 
                                                 predictionBounds = c(T,T),
                                                 allowParallel = T
                        )
    )
    
    set.seed(n)
    macro_svmp <- train(Moisture_Uptake ~ ., 
                        data = training_set %>%
                          select_if(is.numeric), 
                        method = "svmPoly", 
                        metric = "Rsquared", 
                        tuneGrid = expand.grid(C = i, scale = c(0.001, 0.01, 0.1), degree = c(1, 2, 3)),
                        trControl = trainControl(method = "cv",
                                                 index = groupKFold(training_set$Genotype, k = 10), 
                                                 savePredictions = T, 
                                                 predictionBounds = c(T,T),
                                                 allowParallel = T
                        )
    )
    
    set.seed(n)
    macro_rf <- train(Moisture_Uptake ~ ., 
                      data = training_set %>%
                        select_if(is.numeric), 
                      method = "rf", 
                      metric = "Rsquared", 
                      tuneGrid = expand.grid(mtry = i),
                      trControl = trainControl(method = "cv",
                                               index = groupKFold(training_set$Genotype, k = 10), 
                                               savePredictions = T, 
                                               predictionBounds = c(T,T),
                                               allowParallel = T
                      )
    )
    
    ###############
    # Predictions #
    ###############
    lm_pred <- predict(macro_lm, validation_set)
    cart_pred <- predict(macro_cart, validation_set)
    knn_pred <- predict(macro_knn, validation_set)
    rf_pred <- predict(macro_rf, validation_set)
    svml_pred <- predict(macro_svml, validation_set)
    svmp_pred <- predict(macro_svmp, validation_set)
    svmr_pred <- predict(macro_svmr, validation_set)
    
    #####################
    # Returning Results #
    #####################
    return(list(lm_macro_random_cor = cor(lm_pred, validation_set$Moisture_Uptake, method = "spearman"),
                cart_macro_random_cor = cor(cart_pred, validation_set$Moisture_Uptake, method = "spearman"),
                knn_macro_random_cor = cor(knn_pred, validation_set$Moisture_Uptake, method = "spearman"),
                rf_macro_random_cor = cor(rf_pred, validation_set$Moisture_Uptake, method = "spearman"),
                svml_macro_random_cor = cor(svml_pred, validation_set$Moisture_Uptake, method = "spearman"),
                svmp_macro_random_cor = cor(svmp_pred, validation_set$Moisture_Uptake, method = "spearman"),
                svmr_macro_random_cor = cor(svmr_pred, validation_set$Moisture_Uptake, method = "spearman"),
                lm_macro_random_rmse = RMSE(lm_pred, validation_set$Moisture_Uptake),
                cart_macro_random_rmse = RMSE(cart_pred, validation_set$Moisture_Uptake),
                knn_macro_random_rmse = RMSE(knn_pred, validation_set$Moisture_Uptake),
                rf_macro_random_rmse = RMSE(rf_pred, validation_set$Moisture_Uptake),
                svml_macro_random_rmse = RMSE(svml_pred, validation_set$Moisture_Uptake),
                svmp_macro_random_rmse = RMSE(svmp_pred, validation_set$Moisture_Uptake),
                svmr_macro_random_rmse = RMSE(svmr_pred, validation_set$Moisture_Uptake),
                lm_macro_random_int = macro_lm$bestTune[[1]],
                cart_macro_random_cp = macro_cart$bestTune[[1]],
                knn_macro_random_k = macro_knn$bestTune[[1]],
                rf_macro_random_mtry = macro_rf$bestTune[[1]],
                svml_macro_random_cost = macro_svml$bestTune[[1]],
                svmp_macro_random_cost = macro_svmp$bestTune[[3]],
                svmp_macro_random_scale = macro_svmp$bestTune[[2]],
                svmp_macro_random_degree = macro_svmp$bestTune[[1]],
                svmr_macro_random_cost = macro_svmr$bestTune[[2]],
                svmr_macro_random_sigma = macro_svmr$bestTune[[1]],
                repetition = n
    ))
  }
  
  genosplit_genocv_macro_results_tibble <- genosplit_genocv_macro_results %>%
    as_tibble() %>%
    unnest(cols = everything())
  
  Macro_Genosplit_Genocv_Results <- bind_rows(Macro_Genosplit_Genocv_Results, genosplit_genocv_macro_results_tibble)
  
}

########################################
# Stop Parallelization and Write Files #
########################################
stopImplicitCluster()

Macro_Genosplit_Genocv_Results %>% 
  pivot_longer(cols = -c(repetition,              # this section will "tidy" the data to make it more manageable
                         contains("ncomp"), 
                         contains("_cp"), 
                         contains("deg"), 
                         contains("scale"), 
                         contains("sigma"),
                         contains("_int"),
                         contains("_k"),
                         contains("_cost"),
                         contains("mtry")
  ),
  names_to = "Mod_Data_Met", 
  values_to = "Performance") %>%
  mutate(Model = case_when(str_detect(string = Mod_Data_Met, pattern = "lm") ~ "LM", # this section creates new columns for the different models
                           str_detect(string = Mod_Data_Met, pattern = "knn") ~ "KNN",
                           str_detect(string = Mod_Data_Met, pattern = "rf") ~ "RF",
                           str_detect(string = Mod_Data_Met, pattern = "svml") ~ "SVML",
                           str_detect(string = Mod_Data_Met, pattern = "svmp") ~ "SVMP",
                           str_detect(string = Mod_Data_Met, pattern = "svmr") ~ "SVMR",
                           str_detect(string = Mod_Data_Met, pattern = "cart") ~ "CART"),
         Metric = case_when(str_detect(string = Mod_Data_Met, pattern = "rmse") ~ "RMSE", # this section creates new columns for the different metrics
                            str_detect(string = Mod_Data_Met, pattern = "cor") ~ "Spearman_R"),
         Data = case_when(str_detect(string = Mod_Data_Met, pattern = "macro") ~ "Macro")) %>%
  select(-Mod_Data_Met) %>% # getting rid of the messy colum with mulitple pieces of information since they each have their own column  now
  mutate(Intercept = case_when(Model == "LM" ~ as.double(lm_macro_random_int)),
         Mtry = case_when(Model == "RF" ~ as.double(rf_macro_random_mtry)),
         C = case_when(Model == "SVML" ~ as.double(svml_macro_random_cost),
                       Model == "SVMP" ~ as.double(svmp_macro_random_cost),
                       Model == "SVMR" ~ as.double(svmr_macro_random_cost),),
         Degree = case_when(Model == "SVMP" ~ as.double(svmp_macro_random_degree),),
         Scale = case_when(Model == "SVMP" ~ as.double(svmp_macro_random_scale),),
         Sigma = case_when(Model == "SVMR" ~ as.double(svmr_macro_random_sigma)),
         K = case_when(Model == "KNN" ~ as.double(knn_macro_random_k)),
         CP = case_when(Model == "CART" ~ as.double(cart_macro_random_cp))) %>%
  select(repetition, Model, Metric, Data, Intercept, Mtry, CP, K, C, Degree, Scale, Sigma, Performance) %>% # getting rid of the extra columns we dont need
  group_by(Model, Metric, Data, Intercept, Mtry, CP, K, C) %>%
  summarise(Mean_Degree = mean(Degree),
            Mean_Scale = mean(Scale),
            Mean_Sigma = mean(Sigma),
            Mean_Performance = mean(Performance),
            Std_Dev_Performance = sd(Performance),
            Range_Performance = max(Performance) - min(Performance),
            N = n(),
            Splitting = "Genotype") %>%
  arrange(desc(Mean_Performance)) %>% # arranging data from highest performance to lowest
  write_csv("../Results/macro_genosplit_genocv_nobound_clean.csv")

Macro_Genosplit_Genocv_Results %>%
  write_csv("../Results/macro_genosplit_genocv_nobound_raw.csv")

print('No Bound Genotype Split Genotype CV Models Finished')

#############################################
# Splitting Data By Genotype with Random CV #
#############################################

num_cores_limit <- 5
iterations <- 5

if (detectCores() > num_cores_limit) {
  num_cores <- num_cores_limit
} else {
  num_cores <- detectCores()
}

registerDoParallel(cores = num_cores)

Macro_Genosplit_Randomcv_Results <- tibble()
print("Entering Genotype - Random Loop")
for(n in 1:100) {
  genosplit_randomcv_macro_results <- foreach(i = 1:iterations, .combine = rbind) %dopar% {
    ##################
    # Splitting Data #
    ##################
    training_genotypes <- dataset %>% 
      select(Genotype) %>% 
      unique() %>%
      sample_frac(size = 0.8) %>%
      pull()
    
    training_set <- dataset %>%
      filter(Genotype %in% training_genotypes)
    
    validation_set <- dataset %>%
      filter(!Genotype %in% training_genotypes)
    
    #############
    # Modelling #
    #############
    set.seed(n)
    macro_lm <- train(Moisture_Uptake ~ ., 
                      data = training_set %>%
                        select_if(is.numeric), 
                      method = "lm", 
                      metric = "Rsquared", 
                      tuneGrid = expand.grid(intercept = 1/i),
                      trControl = trainControl(method = "cv",
                                               #index = groupKFold(training_set$Genotype, k = 10), 
                                               savePredictions = T, 
                                               predictionBounds = c(T,T),
                                               allowParallel = T
                      )
    )   
    
    set.seed(n)
    macro_cart <- train(Moisture_Uptake ~ ., 
                        data = training_set %>%
                          select_if(is.numeric), 
                        method = "rpart", 
                        metric = "Rsquared", 
                        tuneGrid = expand.grid(cp = i*0.02),
                        trControl = trainControl(method = "cv",
                                                 #index = groupKFold(training_set$Genotype, k = 10), 
                                                 savePredictions = T, 
                                                 predictionBounds = c(T,T),
                                                 allowParallel = T
                        )
    )
    
    set.seed(n)
    macro_knn <- train(Moisture_Uptake ~ ., 
                       data = training_set %>%
                         select_if(is.numeric), 
                       method = "knn", 
                       metric = "Rsquared", 
                       tuneGrid = expand.grid(k = i),
                       trControl = trainControl(method = "cv",
                                                #index = groupKFold(training_set$Genotype, k = 10), 
                                                savePredictions = T, 
                                                predictionBounds = c(T,T),
                                                allowParallel = T
                       )
    )
    
    set.seed(n)
    macro_rf <- train(Moisture_Uptake ~ ., 
                      data = training_set %>%
                        select_if(is.numeric), 
                      method = "rf", 
                      metric = "Rsquared", 
                      tuneGrid = expand.grid(mtry = i),
                      trControl = trainControl(method = "cv",
                                               #index = groupKFold(training_set$Genotype, k = 10), 
                                               savePredictions = T, 
                                               predictionBounds = c(T,T),
                                               allowParallel = T
                      )
    )
    
    set.seed(n)
    macro_svml <- train(Moisture_Uptake ~ ., 
                        data = training_set %>%
                          select_if(is.numeric), 
                        method = "svmLinear", 
                        metric = "Rsquared",
                        tuneGrid = expand.grid(C = i),
                        trControl = trainControl(method = "cv",
                                                 #index = groupKFold(training_set$Genotype, k = 10), 
                                                 savePredictions = T, 
                                                 predictionBounds = c(T,T),
                                                 allowParallel = T
                        )
    )
    
    set.seed(n)
    macro_svmr <- train(Moisture_Uptake ~ ., 
                        data = training_set %>%
                          select_if(is.numeric), 
                        method = "svmRadial", 
                        metric = "Rsquared",
                        tuneGrid = expand.grid(C = i, sigma = c(0.001, 0.01, 0.1)),
                        trControl = trainControl(method = "cv",
                                                 #index = groupKFold(training_set$Genotype, k = 10), 
                                                 savePredictions = T, 
                                                 predictionBounds = c(T,T),
                                                 allowParallel = T
                        )
    )
    
    set.seed(n)
    macro_svmp <- train(Moisture_Uptake ~ ., 
                        data = training_set %>%
                          select_if(is.numeric), 
                        method = "svmPoly", 
                        metric = "Rsquared", 
                        tuneGrid = expand.grid(C = i, scale = c(0.001, 0.01, 0.1), degree = c(1, 2, 3)),
                        trControl = trainControl(method = "cv",
                                                 #index = groupKFold(training_set$Genotype, k = 10), 
                                                 savePredictions = T, 
                                                 predictionBounds = c(T,T),
                                                 allowParallel = T
                        )
    )
    
    ###############
    # Predictions #
    ###############
    lm_pred <- predict(macro_lm, validation_set)
    cart_pred <- predict(macro_cart, validation_set)
    knn_pred <- predict(macro_knn, validation_set)
    rf_pred <- predict(macro_rf, validation_set)
    svml_pred <- predict(macro_svml, validation_set)
    svmp_pred <- predict(macro_svmp, validation_set)
    svmr_pred <- predict(macro_svmr, validation_set)
    
    #####################
    # Returning Results #
    #####################
    return(list(lm_macro_random_cor = cor(lm_pred, validation_set$Moisture_Uptake, method = "spearman"),
                cart_macro_random_cor = cor(cart_pred, validation_set$Moisture_Uptake, method = "spearman"),
                knn_macro_random_cor = cor(knn_pred, validation_set$Moisture_Uptake, method = "spearman"),
                rf_macro_random_cor = cor(rf_pred, validation_set$Moisture_Uptake, method = "spearman"),
                svml_macro_random_cor = cor(svml_pred, validation_set$Moisture_Uptake, method = "spearman"),
                svmp_macro_random_cor = cor(svmp_pred, validation_set$Moisture_Uptake, method = "spearman"),
                svmr_macro_random_cor = cor(svmr_pred, validation_set$Moisture_Uptake, method = "spearman"),
                lm_macro_random_rmse = RMSE(lm_pred, validation_set$Moisture_Uptake),
                cart_macro_random_rmse = RMSE(cart_pred, validation_set$Moisture_Uptake),
                knn_macro_random_rmse = RMSE(knn_pred, validation_set$Moisture_Uptake),
                rf_macro_random_rmse = RMSE(rf_pred, validation_set$Moisture_Uptake),
                svml_macro_random_rmse = RMSE(svml_pred, validation_set$Moisture_Uptake),
                svmp_macro_random_rmse = RMSE(svmp_pred, validation_set$Moisture_Uptake),
                svmr_macro_random_rmse = RMSE(svmr_pred, validation_set$Moisture_Uptake),
                lm_macro_random_int = macro_lm$bestTune[[1]],
                cart_macro_random_cp = macro_cart$bestTune[[1]],
                knn_macro_random_k = macro_knn$bestTune[[1]],
                rf_macro_random_mtry = macro_rf$bestTune[[1]],
                svml_macro_random_cost = macro_svml$bestTune[[1]],
                svmp_macro_random_cost = macro_svmp$bestTune[[3]],
                svmp_macro_random_scale = macro_svmp$bestTune[[2]],
                svmp_macro_random_degree = macro_svmp$bestTune[[1]],
                svmr_macro_random_cost = macro_svmr$bestTune[[2]],
                svmr_macro_random_sigma = macro_svmr$bestTune[[1]],
                repetition = n
    ))
  }
  
  genosplit_randomcv_macro_results_tibble <- genosplit_randomcv_macro_results %>%
    as_tibble() %>%
    unnest(cols = everything())
  
  Macro_Genosplit_Randomcv_Results <- bind_rows(Macro_Genosplit_Randomcv_Results, genosplit_randomcv_macro_results_tibble)
  
}

########################################
# Stop Parallelization and Write Files #
########################################
stopImplicitCluster()

Macro_Genosplit_Randomcv_Results %>% 
  pivot_longer(cols = -c(repetition,              # this section will "tidy" the data to make it more manageable
                         contains("ncomp"), 
                         contains("_cp"), 
                         contains("_int"),
                         contains("_k"),
                         contains("mtry")
  ),
  names_to = "Mod_Data_Met", 
  values_to = "Performance") %>%
  mutate(Model = case_when(str_detect(string = Mod_Data_Met, pattern = "lm") ~ "LM", # this section creates new columns for the different models
                           str_detect(string = Mod_Data_Met, pattern = "knn") ~ "KNN",
                           str_detect(string = Mod_Data_Met, pattern = "rf") ~ "RF",
                           str_detect(string = Mod_Data_Met, pattern = "cart") ~ "CART"),
         Metric = case_when(str_detect(string = Mod_Data_Met, pattern = "rmse") ~ "RMSE", # this section creates new columns for the different metrics
                            str_detect(string = Mod_Data_Met, pattern = "cor") ~ "Spearman_R"),
         Data = case_when(str_detect(string = Mod_Data_Met, pattern = "macro") ~ "Macro")) %>%
  select(-Mod_Data_Met) %>% # getting rid of the messy colum with mulitple pieces of information since they each have their own column  now
  mutate(Intercept = case_when(Model == "LM" ~ as.double(lm_macro_random_int)),
         Mtry = case_when(Model == "RF" ~ as.double(rf_macro_random_mtry)),
         K = case_when(Model == "KNN" ~ as.double(knn_macro_random_k)),
         CP = case_when(Model == "CART" ~ as.double(cart_macro_random_cp))) %>%
  select(repetition, Model, Metric, Data, Intercept, Mtry, CP, K, Performance) %>% # getting rid of the extra columns we dont need
  group_by(Model, Metric, Data, Intercept, Mtry, CP, K) %>%
  summarise(Mean_Performance = mean(Performance),
            Std_Dev_Performance = sd(Performance),
            Range_Performance = max(Performance) - min(Performance),
            N = n(),
            Splitting = "Environment") %>%
  arrange(desc(Mean_Performance)) %>% # arranging data from highest performance to lowest
  write_csv("../Results/macro_genosplit_randomcv_nobound_clean.csv")

Macro_Genosplit_Randomcv_Results %>%
  write_csv("../Results/macro_genosplit_randomcv_nobound_raw.csv")

print('No Bound Genotype Split Random CV Models Finished')


