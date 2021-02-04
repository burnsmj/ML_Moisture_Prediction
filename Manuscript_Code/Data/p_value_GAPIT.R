### GAPIT P Value - Moisture Uptake
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
memory.limit(size = 35000)

###############################
# Getting Shell Script Inputs #
###############################
cli_arg <- commandArgs(trailingOnly = TRUE)
iter <- as.numeric(cli_arg[1])

################
# Loading Data #
################
#n5000_lmer_data <- read.csv("../Data/N5000_Master_With_Predictions.csv")
#snp_data <- read.table("../Data/christine_common_final.hmp.txt", head = F) # Had to increase device memory allocation to R in order to get this loaded.
myY <- read.csv(cat("../Data/MyY_Env_", iter, ".csv", sep = ""))
myGD <- read.big.matrix(cat("myGD_", iter, ".csv", sep = ""),type = "char", head=T)
myGM <- read.csv(cat("myGM_", iter, ".csv", sep = ""), header=T)

#################################
# Determining P-Value Threshold #
#################################
FarmCPU.P.Threshold(
  Y=myY, #only two columns allowed, the first column is taxa name and the second is phenotype value
  GD=myGD,
  GM=myGM,
  trait=cat("moisture_uptake_", iter, sep = ""), #name of the trait, only used for the output file name
  theRep=100 #number of permutation times
)
