##################
# Load Libraries #
##################
library(tidyverse)
library(lme4)
library(magrittr)

#########################
# Set Working Directory #
#########################
#setwd("/home/hirschc1/burns756/Machine_Learning/Data/") # Only for use in MSI

################
# Loading Data #
################
spectra_training_data <- read_csv("Data/Moisture_Uptake_Master_Dataset_Cook_Macro_Spectra.csv") %>%
  select(c(1,2,26)) %>%
  mutate(Moisture_Uptake = Moisture_Uptake * 100)

predictions <- read_csv("Data/SNP_Data//N5000_Prediction_Dataset.csv") %>% # Loading N5000 dataset and keeping the prediction, genotype, env, rep, and block.
  select(2:5, 152) %>%
  mutate(SVML_Prediction = SVML_Prediction * 100) %>% # Putting moisture uptake into percentage rather than proportion
  mutate(Genotype = factor(Genotype),
         Block = factor(Block),
         Rep = factor(Rep),
         Env = factor(Env))

myG_Taxa <- read.table("Data/SNP_Data/GD_Taxa_List.txt", head = T)

############################
# Filtering by Environment #
############################
env1 <- predictions %>%
  filter(Env == 1) %>%
  filter(SVML_Prediction >= min(spectra_training_data$Moisture_Uptake) & SVML_Prediction <= max(spectra_training_data$Moisture_Uptake))

env2 <- predictions %>%
  filter(Env == 2) %>%
  filter(SVML_Prediction >= min(spectra_training_data$Moisture_Uptake) & SVML_Prediction <= max(spectra_training_data$Moisture_Uptake))

env3 <- predictions %>%
  filter(Env == 3) %>%
  filter(SVML_Prediction >= min(spectra_training_data$Moisture_Uptake) & SVML_Prediction <= max(spectra_training_data$Moisture_Uptake))

env4 <- predictions %>%
  filter(Env == 4) %>%
  filter(SVML_Prediction >= min(spectra_training_data$Moisture_Uptake) & SVML_Prediction <= max(spectra_training_data$Moisture_Uptake))

env5 <- predictions %>%
  filter(Env == 5) %>%
  filter(SVML_Prediction >= min(spectra_training_data$Moisture_Uptake) & SVML_Prediction <= max(spectra_training_data$Moisture_Uptake))

dim(env1)
dim(env2)
dim(env3)
dim(env4)
dim(env5)

####################
# Extracting BLUPs #
####################
lmer1 <- lmer(data = env1, SVML_Prediction ~ (1 | Genotype) + (1 | Rep/Block), REML = T)
blup1 <- ranef(lmer1)$Genotype
summary(lmer1)
lmer2 <- lmer(data = env2, SVML_Prediction ~ (1 | Genotype) + (1 | Rep/Block), REML = T)
blup2 <- ranef(lmer2)$Genotype
summary(lmer2)
lmer3 <- lmer(data = env3, SVML_Prediction ~ (1 | Genotype) + (1 | Rep/Block), REML = T)
blup3 <- ranef(lmer3)$Genotype
summary(lmer3)
lmer4 <- lmer(data = env4, SVML_Prediction ~ (1 | Genotype) + (1 | Rep/Block), REML = T)
blup4 <- ranef(lmer4)$Genotype
summary(lmer4)
lmer5 <- lmer(data = env5, SVML_Prediction ~ (1 | Genotype) + (1 | Rep/Block), REML = T)
blup5 <- ranef(lmer5)$Genotype
summary(lmer5)

dim(blup1)
dim(blup2)
dim(blup3)
dim(blup4)
dim(blup5)

##########################
# Creating BLUP Datasets #
##########################
Env_1_BLUPs <- tibble(taxa = rownames(blup1),
                      BLUPs_1 = blup1[[1]])
Env_2_BLUPs <- tibble(taxa = rownames(blup2),
                      BLUPs_2 = blup2[[1]])
Env_3_BLUPs <- tibble(taxa = rownames(blup3),
                      BLUPs_3 = blup3[[1]])
Env_4_BLUPs <- tibble(taxa = rownames(blup4),
                      BLUPs_4 = blup4[[1]])
Env_5_BLUPs <- tibble(taxa = rownames(blup5),
                      BLUPs_5 = blup5[[1]])

dim(Env_1_BLUPs)
dim(Env_2_BLUPs)
dim(Env_3_BLUPs)
dim(Env_4_BLUPs)
dim(Env_5_BLUPs)

#######################################
# List of Genotypes with Genomic Data #
#######################################
genotypes <- myG_Taxa$taxa

#########################################
# Matching BLUP and Genomic Information #
#########################################
Env_1_Y <- Env_1_BLUPs %>%
  filter(taxa %in% genotypes)

Env_2_Y <- Env_2_BLUPs %>%
  filter(taxa %in% genotypes)

Env_3_Y <- Env_3_BLUPs %>%
  filter(taxa %in% genotypes)

Env_4_Y <- Env_4_BLUPs %>%
  filter(taxa %in% genotypes)

Env_5_Y <- Env_5_BLUPs %>%
  filter(taxa %in% genotypes)

dim(Env_1_Y)
dim(Env_2_Y)
dim(Env_3_Y)
dim(Env_4_Y)
dim(Env_5_Y)

genos_1 <- tibble(taxa = genotypes) %>%
  filter(taxa %in% Env_1_Y$taxa)
genos_2 <- tibble(taxa = genotypes) %>%
  filter(taxa %in% Env_2_Y$taxa)
genos_3 <- tibble(taxa = genotypes) %>%
  filter(taxa %in% Env_3_Y$taxa)
genos_4 <- tibble(taxa = genotypes) %>%
  filter(taxa %in% Env_4_Y$taxa)
genos_5 <- tibble(taxa = genotypes) %>%
  filter(taxa %in% Env_5_Y$taxa)

dim(genos_1)
dim(genos_2)
dim(genos_3)
dim(genos_4)
dim(genos_5)

Env_1_Y_Match <- genos_1 %>%
  full_join(Env_1_Y)
Env_2_Y_Match <- genos_2 %>%
  full_join(Env_2_Y)
Env_3_Y_Match <- genos_3 %>%
  full_join(Env_3_Y)
Env_4_Y_Match <- genos_4 %>%
  full_join(Env_4_Y)
Env_5_Y_Match <- genos_5 %>%
  full_join(Env_5_Y)

dim(Env_1_Y_Match)
sum(Env_1_Y_Match$taxa == genos_1$taxa)
dim(Env_2_Y_Match)
sum(Env_2_Y_Match$taxa == genos_2$taxa)
dim(Env_3_Y_Match)
sum(Env_3_Y_Match$taxa == genos_3$taxa)
dim(Env_4_Y_Match)
sum(Env_4_Y_Match$taxa == genos_4$taxa)
dim(Env_5_Y_Match)
sum(Env_5_Y_Match$taxa == genos_5$taxa)

### Phenotypic BLUP Data ###
write_csv(Env_1_Y_Match,"Data/SNP_Data/MyY_Env_1.csv")
write_csv(Env_2_Y_Match,"Data/SNP_Data/MyY_Env_2.csv")
write_csv(Env_3_Y_Match,"Data/SNP_Data/MyY_Env_3.csv")
write_csv(Env_4_Y_Match,"Data/SNP_Data/MyY_Env_4.csv")
write_csv(Env_5_Y_Match,"Data/SNP_Data/MyY_Env_5.csv")


##########
### QC ###
##########
Env_1_Y_Match %>%
  summarise(min_blups = min(BLUPs_1),
            max_blups = max(BLUPs_1),
            range = max_blups-min_blups)
Env_2_Y_Match %>%
  summarise(min_blups = min(BLUPs_2),
            max_blups = max(BLUPs_2),
            range = max_blups-min_blups)
Env_3_Y_Match %>%
  summarise(min_blups = min(BLUPs_3),
            max_blups = max(BLUPs_3),
            range = max_blups-min_blups)
Env_4_Y_Match %>%
  summarise(min_blups = min(BLUPs_4),
            max_blups = max(BLUPs_4),
            range = max_blups-min_blups)
Env_5_Y_Match %>%
  summarise(min_blups = min(BLUPs_5),
            max_blups = max(BLUPs_5),
            range = max_blups-min_blups)

read_csv("/Users/michael/Downloads/Moisture_Prediction_Tables - Table_S1 (2).csv", skip = 1) %>%
  select(1,8) %>%
  rename(taxa = Genotype) %>%
  full_join(Env_1_Y_Match) %>%
  filter(!is.na(BLUPs_1)) %>%
  filter(!is.na(Env1_Moisture_BLUP)) %>%
  ggplot(aes(x = BLUPs_1, y = Env1_Moisture_BLUP))+
  geom_point()

read_csv("/Users/michael/Downloads/Moisture_Prediction_Tables - Table_S1 (2).csv", skip = 1) %>%
  select(1,8) %>%
  ggplot(aes(x = Env1_Moisture_BLUP))+
  geom_histogram()+
  xlim(c(-5,5))

read_csv("/Users/michael/Downloads/Moisture_Prediction_Tables - Table_S1 (2).csv", skip = 1) %>%
  select(1,8) %>%
  rename(taxa = Genotype) %>%
  full_join(Env_1_Y_Match) %>%
  filter(!is.na(BLUPs_1)) %>%
  filter(!is.na(Env1_Moisture_BLUP)) %>%
  ggplot(aes(x = BLUPs_1))+
  geom_histogram()+
  xlim(c(-5,5))
