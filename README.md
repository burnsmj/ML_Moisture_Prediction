# Nixtamalization Moisture Content Prediction in a Maize Inbred Diversity Panel
Michael Burns  
University of Minnesota - Applied Plant Sciences

## Purpose:
Food-grade maize (specifically maize used for producing tortillas and tortilla chips) is a small usage and production segment of
commercial maize in the United States. Due to this, food-grade maize breeding programs do not typically have the resources to 
develop tools for breeding platforms specific to food-grade purposes. One trait of interest is nixtamalization moisture
content, which is the amount of water absorbed by kernels during the nixtamalization process. Nixtamalization involves cooking 
maize kernels in an alkaline solution at high temperatures to remove the pericarp and soften the grain for grinding. 
The moisture content of the resulting masa affects its taste, texture, oil content, and more making it an important trait for 
food-grade maize. Multiple methods have been proposed for assessing nixtamalization moisture content, however, they tend to focus
on the manufacturing process and require anywhere from 100g to 1,000 lbs of seed. This is not a viable option in early maize breeding pipelines
where large numbers of small-quantity samples are present. To address this, we utilized Near-Infrared (NIR) spectroscopy combined 
with machine learning to develop a high-throughput method for predicting moisture content during nixtamalization. By training a 
machine learning model on NIR spectra and nixtamalization moisture content data we were able to train an accurate (Spearman R = 0.85) model on only 30g 
of grain per sample. This model was then used to analyze trait variation and genetic factors. The study's approach offers a rapid and scalable solution for 
evaluating nixtamalization moisture content, providing valuable insights for improving food-grade maize germplasm and informing 
breeding strategies for masa-based products.

## Publication:
To learn more about the utility and biological insights of this work, please see:

Burns et al. (2021). Predicting moisture content during maize nixtamalization using machine learning with NIR spectroscopy.
Theoretical and Applied Genetics. 134(11):3743-3757.
doi: 10.1007/s00122-021-03926-8.

## Repository Files:
- Manuscript_Code
  - Data Cleaning
    - Scripts for cleaning spectroscopic and cook test data.
  - Data_Analysis
    - The Rmarkdown that contains the analyses performed in the manuscript.
  - GWAS_Code
    - The code used to run GWAS on the WiDiv panel on the Minnesota Supercomputing Institute's High Performance Computing Clusters.
    - The general pipeline can be found in the readme within this directory.
  - Machine_Learning
    - Scripts developed to train many machine learning models on the Minnesota Supercomputing Institute's High Performance Computing Clusters.
  - Moisture_Uptake_Prediction
    - A proof of concept shiny application to that can be used to predict nixtamalization moisture content in inbred maize samples.
    - Directions for use can be found in the repository.
    - A more complete and applicable version of this application is coming in a future manuscript.
- Data
  - Raw_Data
    - A directory of the raw data collected during the experiment.  
  - Cleaned_Data
    - A directory of the cleaned data created during the analysis.
  - Supplemental_Table_Data
    - A directory of the data used to create the supplemental tables in the manuscript.
  - Note: HapMap data is not included in this directory due to its size. Please see O’Connor et al. 2020 for the hapmap data used.
 
## Pipeline:
The following scripts were utilized in this manuscript and were run locally unless otherwise specified.

  - Data Preparation:
    - Digital_Creation_of_Mastersheet.Rmd
    - digital_cleaning_n5000.Rmd
    - ML_Master_Validation_Dataset.Rmd
    - Digitally_Cleaning_Training_Set.Rmd
      - Used to create manuscript ready table.
  - Machine Learning:
    - MSI_Macro_ML.R
      - Run on HPC
    - MSI_Spectra_ML_NoBound_RandomSplit_RandomCV.R
      - Run on HPC
    - MSI_Spectra_ML_NoBound_GenoSplit_GenoCV.R
      - Run on HPC
    - Non-Parallel_LC.R
  - GWAS:
    - Numerical_GAPIT.R
      - Run on HPC
    - Splitting_BLUPs_by_Env.R
      - Run on HPC
    - GWAS_QC.sh
      - Run on HPC
    - GWAS_Env_Splitter.pl
      - Run on HPC
    - p_value_GAPIT.R
      - Run on HPC
    - Moisture_Uptake_FarmCPU.R
      - Run on HPC
  - Data Analysis:
    - Moisture_Uptake_Prediction_Paper_Figures.Rmd (Aspects of this script were created and used throughout the pipeline to inform decisions, but as it is set up, it can run after the pipeline entirely)
    
