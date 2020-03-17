dataset <- read_csv("Data/Raw_Data/Cooking_Param_Clean.csv")

library(ggridges)

dataset %>%
  mutate(Cook_Time = Cook_Time * 24 * 60,
         Steep_Time = Steep_Time * 24) %>%
  filter(Moisture_Avg < (mean(Moisture_Avg) + 3 * sd(Moisture_Avg)) & 
           Moisture_Avg > (mean(Moisture_Avg) - 3 * sd(Moisture_Avg))) %>%
  filter(Cook_Time < 35 & Cook_Time > 15) %>%
  filter(Steep_Time < 18 & Steep_Time > 14) %>%
  mutate(Cook_Time = (Cook_Time - mean(Cook_Time)) / sd(Cook_Time),
         Steep_Time = (Steep_Time - mean(Steep_Time)) / sd(Steep_Time),
         PH = (PH - mean(PH)) / sd(PH),
         Moisture_Avg = (Moisture_Avg - mean(Moisture_Avg)) / sd(Moisture_Avg)) %>%
  select(c(Genotype, SampleID, Cook_Time, Steep_Time, PH, Moisture_Avg)) %>%
  pivot_longer(cols = -c(Genotype, SampleID), names_to = "Variable", values_to = "Value") %>%
  ggplot(aes(x = Value, y = Variable))+
  geom_density_ridges()


