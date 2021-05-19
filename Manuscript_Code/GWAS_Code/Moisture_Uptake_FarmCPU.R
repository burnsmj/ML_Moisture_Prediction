### FarmCPU - Moisture Uptake
### Michael Burns
### 12/21/2020

#############
# Libraries #
#############
library("tidyverse")
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

###############################
# Getting Shell Script Inputs #
###############################
cli_arg <- commandArgs(trailingOnly = TRUE)
iter <- as.numeric(cli_arg[1])
#iter = 1

#########################
# Set Working Directory #
#########################
setwd(paste0("/home/hirschc1/burns756/Machine_Learning/Data/Environment_", iter))

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

################
# Loading Data #
################
myY <- read.csv(paste0("../MyY_Env_", iter, ".csv"), header = T)
myGD <- read.big.matrix(paste0("../GD_Env_", iter, ".txt"), sep = "\t", type = "char", head=T)
myGM <- read.delim("../GAPIT.Genotype.map.txt", sep = "\t", header=T)
myPCA <- read.csv("../GAPIT.PCA.csv", header=T)
myPValue <- read.delim(paste0("../FarmCPU.p.threshold.optimize.moisture_uptake_", iter, ".txt"), header = F, sep = "\t")

pval <- print(sort(myPValue$V1)[ceiling(30*0.05)]) #quantile(myPValue$V1, 0.05)

#########################
# Filtering PCA Dataset #
#########################
newPCA <- myPCA %>%
  filter(taxa %in% myY$taxa) %>%
  as_data_frame()

######
# QC #
######
length(myY$taxa)
length(newPCA$taxa)
sum(myY$taxa == newPCA$taxa) 
paste("P value: ", pval)
write(newPCA$taxa, paste("PCA_taxa_", iter, ".csv"))
write(myY$taxa, paste("MyY_Taxa_", iter, ".csv"))

#################
# FarmCPU Model #
#################
myFarmCPU <- FarmCPU(
  Y=myY, #phenotype
  GD=myGD, #Genotype matrix
  GM=myGM, #Genotypic map
  CV=newPCA[,-1], #Covariate variables (First 5 PCAs from GAPIT); taxa should not be included.
  threshold.output=1, #P value smaller than threshold.output will be output in GWAS table
  p.threshold=pval,
  MAF.calculate=TRUE, #Calculate minor allele frequency (MAF) or not, if set to TRUE, the SNPs with a lower MAF (<maf.threshold) will be deleted
  method.bin="optimum",
  maf.threshold=0.05, #When MAF.calculate=TRUE, the SNPs with a lower MAF (<maf.threshold) will be deleted
  maxLoop=50, #Maximum number of iterations allowed
  memo = paste("Env_", iter, sep = "") #Add extension to file name for parallel runs
)

