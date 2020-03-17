### Michael Burns
### 2020-02-19
### Width of Spectra Learning Set

dataset <- read_xlsx("Data/Spectra_Data/spectra_learning_sets.xlsx")
view(dataset)

dataset %>%
  ggplot(aes(x = Moisture_Avg)) +
  geom_histogram(binwidth = 0.01, fill = "darkblue", color = "gray")+
  xlab("Moisture Uptake")+
  labs(title = "Moisture Uptake Spread")

dataset %>%
  mutate(moisture_z = (Moisture_Avg - mean(Moisture_Avg)) / sd(Moisture_Avg)) %>%
  ggplot(aes(x = moisture_z)) +
  geom_histogram(binwidth = 0.4, fill = "darkblue", color = "gray")+
  xlab("Moisture Uptake")+
  labs(title = "Moisture Uptake Spread")


dataset %>%
  summarise(min = min(Moisture_Avg),
            max = max(Moisture_Avg))

dataset %>% 
  select(-c(Genotype, Env, Rep, Block, ERB_Group, Amylopectin_average, Ankom_Crude_Fiber_average, Ash, 
            Crude.fat, Crude.fiber, Fructose, Glucose, N_combustion_average, N_Kjeltec_average, 
            Sucrose_average, Total_Sugars_average, Moisture, Protein_As_is, Fat_As_is, Fiber_As_is, 
            Ash_As_is, Starch_As_is, DML_Percent)) %>%
  pivot_longer(cols = -c(SampleID, Moisture_Avg), names_to = 'Wavelength', values_to = 'Absorbance') %>%
  ggplot(aes(x = Absorbance)) +
  geom_histogram()
  #facet_grid(~Wavelength, cols = 12)

