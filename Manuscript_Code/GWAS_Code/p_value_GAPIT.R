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
source("http://zzlab.net/GAPIT/emma.txt")
source("http://zzlab.net/GAPIT/gapit_functions.txt")
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
#iter = 1

#########################
# Set Working Directory #
#########################
setwd("/home/hirschc1/burns756/Machine_Learning/Data/")

################
# Loading Data #
################
myY <- read.csv(paste0("MyY_Env_", iter, ".csv"), header = T)
myGD <- read.big.matrix(paste0("GD_Env_", iter, ".txt"), sep = "\t",type = "char", head=T)
myGM <- read.delim("GAPIT.Genotype.map.txt", sep = "\t", header=T)

#################################
# Determining P-Value Threshold #
#################################
FarmCPU.P.Threshold(
  Y=myY, #only two columns allowed, the first column is taxa name and the second is phenotype value
  GD=myGD,
  GM=myGM,
  trait=paste0("moisture_uptake_", iter), #name of the trait, only used for the output file name
  theRep=30 #number of permutation times
)
