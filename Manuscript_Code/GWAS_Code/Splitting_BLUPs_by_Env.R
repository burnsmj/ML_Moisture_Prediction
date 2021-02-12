##################
# Load Libraries #
##################
library(tidyverse)
library(lme4)

################
# Loading Data #
################
predictions <- read_csv("Data/N5000_Master_With_Predictions.csv") %>%
  select(2:5, 152) 
myG <- read.table("Data/SNP_Data/christine_common_final.hmp.txt", head = F) #normally this file is christine_common_final.hmp.txt

############################
# Filtering by Environment #
############################
env1 <- predictions %>%
  filter(Env == 1)

env2 <- predictions %>%
  filter(Env == 2)

env3 <- predictions %>%
  filter(Env == 3)

env4 <- predictions %>%
  filter(Env == 4)

env5 <- predictions %>%
  filter(Env == 5)

dim(env1)
dim(env2)
dim(env3)
dim(env4)
dim(env5)

####################
# Extracting BLUPs #
####################
lmer1 <- lmer(data = env1, SVML_Prediction ~ (1 | Genotype) + (1 | Rep) + (1 | Rep/Block) + (1 | Genotype:Env), REML = T)
blup1 <- ranef(lmer1)$Genotype
dim(blup1)

lmer2 <- lmer(data = env2, SVML_Prediction ~ (1 | Genotype) + (1 | Rep) + (1 | Rep/Block) + (1 | Genotype:Env), REML = T)
blup2 <- ranef(lmer2)$Genotype
dim(blup2)

lmer3 <- lmer(data = env3, SVML_Prediction ~ (1 | Genotype) + (1 | Rep) + (1 | Rep/Block) + (1 | Genotype:Env), REML = T)
blup3 <- ranef(lmer3)$Genotype
dim(blup3)

lmer4 <- lmer(data = env4, SVML_Prediction ~ (1 | Genotype) + (1 | Rep) + (1 | Rep/Block) + (1 | Genotype:Env), REML = T)
blup4 <- ranef(lmer4)$Genotype
dim(blup4)

lmer5 <- lmer(data = env5, SVML_Prediction ~ (1 | Genotype) + (1 | Rep) + (1 | Rep/Block) + (1 | Genotype:Env), REML = T)
blup5 <- ranef(lmer5)$Genotype
dim(blup5)

##########################
# Creating BLUP Datasets #
##########################
Env_1_BLUPs <- tibble(taxa = rownames(blup1),
                      BLUPs = blup1[[1]])
Env_2_BLUPs <- tibble(taxa = rownames(blup2),
                      BLUPs = blup2[[1]])
Env_3_BLUPs <- tibble(taxa = rownames(blup3),
                      BLUPs = blup3[[1]])
Env_4_BLUPs <- tibble(taxa = rownames(blup4),
                      BLUPs = blup4[[1]])
Env_5_BLUPs <- tibble(taxa = rownames(blup5),
                      BLUPs = blup5[[1]])

dim(Env_1_BLUPs)
dim(Env_2_BLUPs)
dim(Env_3_BLUPs)
dim(Env_4_BLUPs)
dim(Env_5_BLUPs)

#######################################
# List of Genotypes with Genomic Data #
#######################################
genotypes <- unlist(list(myG[1,c(12:ncol(myG))]))

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

g1_index <- which(genotypes %in% genos_1$geno) + 11
g2_index <- which(genotypes %in% genos_2$geno) + 11
g3_index <- which(genotypes %in% genos_3$geno) + 11
g4_index <- which(genotypes %in% genos_4$geno) + 11
g5_index <- which(genotypes %in% genos_5$geno) + 11

length(g1_index)
length(g2_index)
length(g3_index)
length(g4_index)
length(g5_index)

myG_Env_1 <- myG[,c(1:11, g1_index)]
myG_Env_2 <- myG[,c(1:11, g2_index)]
myG_Env_3 <- myG[,c(1:11, g3_index)]
myG_Env_4 <- myG[,c(1:11, g4_index)]
myG_Env_5 <- myG[,c(1:11, g5_index)]

dim(myG_Env_1) - 11
dim(myG_Env_2) - 11
dim(myG_Env_3) - 11
dim(myG_Env_4) - 11
dim(myG_Env_5) - 11

head(myG_Env_1)

#########################
# Write Out Subset Data #
#########################
### Genomic Data ###
write_delim(myG_Env_1,"Data/SNP_Data/MyG_Env_1.hmp.txt", delim = "\t", col_names = F)
write_delim(myG_Env_2,"Data/SNP_Data/MyG_Env_2.hmp.txt", delim = "\t", col_names = F)
write_delim(myG_Env_3,"Data/SNP_Data/MyG_Env_3.hmp.txt", delim = "\t", col_names = F)
write_delim(myG_Env_4,"Data/SNP_Data/MyG_Env_4.hmp.txt", delim = "\t", col_names = F)
write_delim(myG_Env_5,"Data/SNP_Data/MyG_Env_5.hmp.txt", delim = "\t", col_names = F)

### Phenotypic Data ###
write_csv(Env_1_Y_Match,"Data/SNP_Data/MyY_Env_1.csv")
write_csv(Env_2_Y_Match,"Data/SNP_Data/MyY_Env_2.csv")
write_csv(Env_3_Y_Match,"Data/SNP_Data/MyY_Env_3.csv")
write_csv(Env_4_Y_Match,"Data/SNP_Data/MyY_Env_4.csv")
write_csv(Env_5_Y_Match,"Data/SNP_Data/MyY_Env_5.csv")
