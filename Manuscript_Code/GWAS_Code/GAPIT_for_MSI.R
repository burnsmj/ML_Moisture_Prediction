### GAPIT - Moisture Uptake
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
#cli_arg <- commandArgs(trailingOnly = TRUE)
#iter <- as.numeric(cli_arg[1])
#iter = 1

################
# Loading Data #
################
#myY <- read.csv(paste("../Data/MyY_Env_", iter, ".csv", sep = ""), header = T)
myG <- read.table("../Data/widiv_446g_christine_SNPs.hmp.txt", sep = "\t", header = F) #normally this file is christine_common_final.hmp.txt
#myGD <- read.csv(paste("myGD_", iter, ".csv", sep = ""), header=T)
#myGM <- read.csv(paste("myGM_", iter, ".csv", sep = ""), header=T)

#############################
# Running GAPIT for FarmCPU #
#############################
myGAPIT2 <- GAPIT(
  #Y = myY,
  G = myG,
  #GD = myGD,
  #GM = myGM,
  PCA.total = 5,
  #method.bin = "optimum",
  #model = "FarmCPU"
) 

myPCA <- myGAPIT2$PCA

write.csv(myPCA, "myPCA.widiv.csv", row.names = F)
