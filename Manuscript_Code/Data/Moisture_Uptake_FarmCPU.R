### FarmCPU - Moisture Uptake
### Michael Burns
### 12/21/2020

#############
# Libraries #
#############
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
myY <- read.csv(cat("../Data/MyY_Env_", iter, ".csv", sep = ""))
myGD <- read.big.matrix(cat("myGD_", iter, ".csv", sep = ""),type = "char", head=T)
myGM <- read.csv(cat("myGM_", iter, ".csv", sep = ""), header=T)
myPCA <- read.csv(cat("GAPIT.PCA_", iter, ".csv", sep = ""), header=T)
myPValue <- read.delim(cat("FarmCPU.p.threshold.optimize.moisture_uptake_", iter, ".txt", sep = ""), header = F, sep = "\t")

pval <- quantile(myPValue$V1, 0.05)

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
#BLUPs <- random_effects$Genotype

########################
# Setting Up Variables #
########################
#myY <- tibble(taxa = rownames(BLUPs),
#              Moisture_Uptake_BLUPs = BLUPs[[1]])

#filtered_myY <- myY %>%
#  filter(taxa %in% myPCA$taxa)

#matched_myY <- tibble(taxa = myPCA$taxa) %>%
#  full_join(filtered_myY) %>%
#  as.data.frame() %>%
#  write_csv("my_Y.csv")

######
# QC #
######
length(myY$taxa)
length(myPCA$taxa)
sum(myY$taxa == myPCA$taxa) #should all be 408

#tibble(filtered_y = filtered_myY$taxa, 
#       pca_taxa = myPCA$taxa,
#       y_taxa = matched_myY$taxa) %>%
#  write_csv("FarmCPU_Taxa_Order_QC.csv")

#################
# FarmCPU Model #
#################
myFarmCPU <- FarmCPU(
  Y=myY, #phenotype
  GD=myGD, #Genotype matrix
  GM=myGM, #Genotypic map
  CV=myPCA[,-1], #Covariate variables (First 5 PCAs from GAPIT); taxa should not be included.
  threshold.output=0.01, #P value smaller than threshold.output will be output in GWAS table
  p.threshold=pval,
  MAF.calculate=TRUE, #Calculate minor allele frequency (MAF) or not, if set to TRUE, the SNPs with a lower MAF (<maf.threshold) will be deleted
  method.bin="optimum",
  maf.threshold=0.05, #When MAF.calculate=TRUE, the SNPs with a lower MAF (<maf.threshold) will be deleted
  maxLoop=50 #Maximum number of iterations allowed
)
