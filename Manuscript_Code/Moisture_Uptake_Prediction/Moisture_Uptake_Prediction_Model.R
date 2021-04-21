### Moisutre Uptake Prediction Model Function
### Michael Burns
### 4/21/21

# Function to predict moisture uptake
Moisture_Uptake_Prediction <- function(){
  # Source the setup file if needed
  if(!exists("r_svml_spectra")){
    source("Prediction_Model_Setup.R")
  }
  
  # Get input file from user
  input <- file.path(readline(prompt = "Enter Input CSV File Path: "))
  output <- file.path(readline(prompt = "Enter Output CSV File Path: "))
  
  # Remove special characters (\") from file paths
  if(str_detect(string = input, pattern = "\"")){
     new_input <- str_remove_all(string = input, pattern = "\"")
  }
  if(str_detect(string = output, pattern = "\"")){
    new_output <- str_remove_all(string = output, pattern = "\"")
  }
  
  # Read and normalize input file
  print("Reading in and normalizing input data")
  input_data <- suppressMessages(read_csv(new_input)) 
  normalized_data <- input_data %>%
    pivot_longer(cols = -c(Sample_ID), names_to = "Waveband", values_to = "Absorbance") %>%
    group_by(Sample_ID) %>%
    mutate(sum_lxl_abs = sum(abs(Absorbance)),
           Norm_Abs = Absorbance / sum_lxl_abs) %>%
    select(-c(Absorbance, sum_lxl_abs)) %>%
    ungroup() %>%
    pivot_wider(id_cols = c(Sample_ID), 
                values_from = Norm_Abs, 
                names_from = Waveband)
  print("Data has been normalized")
  
  # Make predictions
  print("Predicting on new data")
  predictions <- predict(r_svml_spectra, normalized_data)
  print("Predictions complete")
  
  # Generate table to write out
  print("Generating output file")
  input_data %>%
    mutate(Predicted_Moisture_Uptake = predictions) %>%
    write_csv(new_output)
  print("Output file has be written to the desired directory")
}
