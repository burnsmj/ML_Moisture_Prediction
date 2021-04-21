### Moisture Uptake Prediction Model Setup
### Michael J Burns
### 4/21/21

# Print a welcome statement to introduce the user to the software
print("Welcome to the Moisture Uptake Prediction Model")

# Load required packages
print("Loading required packages (tidyverse, caret)")
require(tidyverse) 
require(caret)
require(stringr)
print("Packages successfully loaded")

# Load training data and normalize it with an absolute value norm
print("Loading and normalizing training data")
training_norm <- suppressMessages(read_csv("Moisture_Uptake_Master_Dataset_Cook_Macro_Spectra.csv")) %>%
  select(-c(3:25)) %>%
  pivot_longer(cols = -c(Sample_ID, Genotype, Moisture_Uptake), names_to = "Waveband", values_to = "Absorbance") %>%
  group_by(Sample_ID) %>%
  mutate(sum_lxl_abs = sum(abs(Absorbance)),
         Norm_Abs = Absorbance / sum_lxl_abs) %>%
  select(-c(Absorbance, sum_lxl_abs)) %>%
  ungroup() %>%
  pivot_wider(id_cols = c(Sample_ID, Genotype, Moisture_Uptake), 
              values_from = Norm_Abs, 
              names_from = Waveband)
print("Training data has been loaded and normalized, let's double check the dimensions")
if(dim(training_norm)[1] == 316 && dim(training_norm)[2] == 144){ # hard coded to ensure the same dataset is used.  will proceed anyways, but with a warning.
  print("Training dataset has correct dimensions")
} else{print("Warning: something is wrong with the training data")}

#Train the svmLinear model
print("Now training the prediction model - this will take a few moments")
r_svml_spectra <- train(Moisture_Uptake ~ ., 
                        data = training_norm %>%
                          select_if(is.numeric), 
                        method = "svmLinear", 
                        metric = "Rsquared", 
                        tuneGrid = expand.grid(C = 71.407),
                        trControl = trainControl(method = "cv",
                                                 index = groupKFold(training_norm$Genotype, k = 10), 
                                                 savePredictions = T,
                                                 allowParallel = T
                                                 )
                        )
print("Model has been trained")

# Predict on validation dataset to confirm it still predicts the same
print("Checking that the model performs the same")
validation_norm <- read_csv("ML_Master_Validation_Dataset.csv") %>%
  filter(SampleID != "YC16:1029") %>%
  select(-c(3:8)) %>%
  pivot_longer(cols = -c(SampleID, Genotype), names_to = "Waveband", values_to = "Absorbance") %>%
  group_by(SampleID) %>%
  mutate(sum_lxl_abs = sum(abs(Absorbance)),
         Norm_Abs = Absorbance / sum_lxl_abs) %>%
  select(-c(Absorbance, sum_lxl_abs)) %>%
  ungroup() %>%
  pivot_wider(id_cols = c(SampleID, Genotype), 
              values_from = Norm_Abs, 
              names_from = Waveband)

val_pred <-  predict(r_svml_spectra, validation_norm)

val_pred_table <- tibble(Predicted_Moisture_Uptake = val_pred)

suppressMessages(orig_val_pred <- suppressMessages(read_csv("Validation_Set_Predictions.csv")))

if(sum(val_pred_table == orig_val_pred) == nrow(val_pred_table)){
  print("Validation predictions match original predictions")
} else{print("Warning: validation predictions no longer match original predictions")}
