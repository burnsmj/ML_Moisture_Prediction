### Michael Burns
### 2020-02-21
### Clustering of wavelengths

dataset <- read_xlsx("Data/Spectra_Data/spectra_learning_sets.xlsx")
view(dataset)

dataset %>% 
  select(-c(Genotype, Env, Rep, Block, ERB_Group, Amylopectin_average, Ankom_Crude_Fiber_average, Ash, 
            Crude.fat, Crude.fiber, Fructose, Glucose, N_combustion_average, N_Kjeltec_average, 
            Sucrose_average, Total_Sugars_average, Moisture, Protein_As_is, Fat_As_is, Fiber_As_is, 
            Ash_As_is, Starch_As_is, DML_Percent)) %>%
  pivot_longer(cols = -c(SampleID, Moisture_Avg), names_to = 'Wavelength', values_to = 'Absorbance') %>%
  ggplot(aes(x = Absorbance, y = Moisture_Avg, color = Wavelength, alpha = 0.1)) +
  geom_point() +
  theme(legend.position = "none")
  
dataset %>% 
  select(-c(Genotype, Env, Rep, Block, ERB_Group, Amylopectin_average, Ankom_Crude_Fiber_average, Ash, 
            Crude.fat, Crude.fiber, Fructose, Glucose, N_combustion_average, N_Kjeltec_average, 
            Sucrose_average, Total_Sugars_average, Moisture, Protein_As_is, Fat_As_is, Fiber_As_is, 
            Ash_As_is, Starch_As_is, DML_Percent)) %>%
  pivot_longer(cols = -c(SampleID, Moisture_Avg), names_to = 'Wavelength', values_to = 'Absorbance') %>%
  ggplot(aes(x = Wavelength, y = Moisture_Avg, alpha = Absorbance)) +
  geom_point() +
  theme(legend.position = "none")

# The two graphs made here don't seem to show a pattern between the moisture average, wavelength, 
# and absorbance.



