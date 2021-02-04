### Numerical GAPIT - Moisture Uptake
### Michael Burns
### 11/12/2020

#####################
# Loading Libraries #
#####################
library("multtest")
library("snpStats")
library("gplots")
library("LDheatmap")
library("genetics")
library("ape")
library("EMMREML")
library("compiler")
library("scatterplot3d")
library("bigmemory")
library("biganalytics")
library("tidyverse")

##############################
# Sourcing GAPIT and FarmCPU #
##############################
source("http://www.zzlab.net/GAPIT/emma.txt")
source("http://www.zzlab.net/GAPIT/gapit_functions.txt")
source("http://zzlab.net/FarmCPU/FarmCPU_functions.txt")

########################
# Setting memory limit #
########################
memory.limit()
memory.limit(size = 35000) # Had to increase device memory allocation to R in order to get this loaded.

###############################
# Getting Shell Script Inputs #
###############################
cli_arg <- commandArgs(trailingOnly = TRUE)
iter <- as.numeric(cli_arg[1])

################
# Loading Data #
################
#myY <- read.csv("my_Y.csv")
myG <- read.table(paste("../Data/MyG_Env_", iter, ".hmp.txt", sep = ""), sep = "\t", header = F) #normally this file is christine_common_final.hmp.txt

#######################################
# Matching Genotypes Between Datasets #
#######################################
#print("Matching Genotypes")

#myY_genos <- myY %>%
#  select(taxa) %>%
#  unique() %>%
#  pull()

#myG_Genos <- myG[1,c(12:ncol(myG))] %>%
#  pivot_longer(cols = everything(), names_to = "Column", values_to = "Genotype") %>%
#  filter(Genotype %in% myY_genos) %>%
#  select(Genotype) %>%
#  unique() %>%
#  pull()

#myG_columns <- myG[1,c(12:ncol(myG))] %>%
#  pivot_longer(cols = everything(), names_to = "Column", values_to = "Genotype") %>%
#  filter(Genotype %in% myY_genos) %>%
#  select(Column) %>%
#  unique() %>%
#  pull()

#myY <- myY %>%
#  filter(taxa %in% myG_Genos) %>%
#  as.data.frame()
#myG <- myG %>%
#  select(c(1:11, myG_columns)) %>%
#  as.data.frame()

############################################################
# Running Gapit for Conversion of HapMap to Numeric Format #
############################################################
print("Starting GAPIT")

myGAPIT <- GAPIT(
  G = myG,
  output.numerical=T
)

print("Finished GAPIT")

myGD <- myGAPIT$GD
myGM <- myGAPIT$GM
#myPCA <- myGAPIT$PCA

print("Saved Files")

#write.csv(myG, "myG.csv")
write.csv(myGD, paste("myGD_", iter, ".csv", sep = ""))
write.csv(myGM, paste("myGM_", iter, ".csv", sep = ""))

print("El Fin")