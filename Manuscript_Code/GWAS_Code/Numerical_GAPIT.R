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

#########################
# Set Working Directory #
#########################
setwd("/home/hirschc1/burns756/Machine_Learning/Data/")

################
# Loading Data #
################
myG <- read.table("christine_common_final.hmp.txt", sep = "\t", header = F) #normally this file is christine_common_final.hmp.txt

############################################################
# Running Gapit for Conversion of HapMap to Numeric Format #
############################################################
print("Starting GAPIT")

myGAPIT <- GAPIT(
  G = myG,
  output.numerical=T,
  PCA.total = 5
)

print("Finished GAPIT")

print("El Fin")