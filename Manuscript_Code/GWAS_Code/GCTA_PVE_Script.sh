#! /bin/bash

# Remember when using GCTA that the phenotype file must have at least three columns (Family ID, Individual ID, and Phenotype), and no header line. For our analysis, the family ID column is all -9 (missing)
# For this analysis, GCTA from cnsgenomics.com is used to determine the percent variance explained by significant SNPs and SNPs as a whole in each environment.
# Version 1.93.2 was used for this analysis and run locally

# All SNPs
for i in {1..5}; do
	./gcta64 --grm Data/gcta_grm_out_env_${i} --pheno Data/GCTA_Pheno_Env_${i}.txt --reml --out Data/GCTA_PVE_Env_${i}
done

# Significant SNPs
#./gcta64 --grm Data/sig_snps_env_1_grm --pheno Data/GCTA_Pheno_Env_1.txt --reml --out Data/Sig_SNPs_PVE_Env_1
#./gcta64 --grm Data/sig_snps_env_2_grm --pheno Data/GCTA_Pheno_Env_2.txt --reml --out Data/Sig_SNPs_PVE_Env_2
#./gcta64 --grm Data/sig_snps_env_4_grm --pheno Data/GCTA_Pheno_Env_4.txt --reml --out Data/Sig_SNPs_PVE_Env_4
#./gcta64 --grm Data/sig_snps_env_5_grm --pheno Data/GCTA_Pheno_Env_5.txt --reml --out Data/Sig_SNPs_PVE_Env_5
