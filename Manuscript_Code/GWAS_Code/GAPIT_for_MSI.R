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

################
# Loading Data #
################
myY <- read.csv("my_Y.csv")

#######################
# Mixed Effects Model #
#######################
#trait_lmer <- n5000_lmer_data %$%
#  lmer(SVML_Prediction ~ (1 | Genotype) + (1 | Env/Rep) + (1 | Rep/Block) + (1 | Genotype:Env), REML = T) # (1|Env) is left out of the equation since the nesting with rep variable accounts for env.  If you add env itself to the equation, ranova wont work since there are technically two environments. After comparing with and without env models, the addition of env only accounts for 0.000012 units of stdev.

####################
# Extracting BLUPs #
####################
#summary(trait_lmer, correlation=FALSE)
#random_effects <- ranef(trait_lmer) # Write out BLUPs for Genotypes

#myY <- data.frame(taxa = rownames(random_effects$Genotype),
#                  Moisture_Uptake_BLUP = random_effects$Genotype[[1]])
#tail(myY)
####################
# Reading SNP Data #
####################
myGD <- read.big.matrix("Data/SNP_Data/myGD.csv",type = "char", head=T)
myGM <- read.csv("Data/SNP_Data/myGM.csv", header = T)

#############################
# Running GAPIT for FarmCPU #
#############################
myGAPIT2 <- GAPIT(
  Y = myY,
  GD = myGD,
  GM = myGM,
  PCA.total = 5,
  method.bin = "optimum",
  model = "FarmCPU"
) 

