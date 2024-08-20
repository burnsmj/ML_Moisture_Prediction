#! /bin/bash

#!/bin/bash
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=8
#SBATCH --cpus-per-task=1
#SBATCH --mem=10gb
#SBATCH --time=1:00:00
#SBATCH --mail-type=ALL
#SBATCH --mail-user=burns756@umn.edu
#SBATCH -o /home/hirschc1/burns756/Machine_Learning/Results/%j.out
#SBATCH -e /home/hirschc1/burns756/Machine_Learning/Results/%j.err

# Load module
module load plink/1.07

# Run Plink
for i in {1..5}; do
	echo Working on Environment ${i}
	plink --noweb --file Machine_Learning/Data/christine_final_env_${i}_binary_plk --make-bed --out Machine_Learning/Data/christine_final_env_${i}_binary_plk_bed
done

# For running significant snps:
#plink --noweb --file Machine_Learning/Data/sig_snps_env_1_binary_plk --make-bed --out Machine_Learning/Data/sig_snps_env_1_binary_plk_bed
#plink --noweb --file Machine_Learning/Data/sig_snps_env_2_binary_plk --make-bed --out Machine_Learning/Data/sig_snps_env_2_binary_plk_bed
#plink --noweb --file Machine_Learning/Data/sig_snps_env_4_binary_plk --make-bed --out Machine_Learning/Data/sig_snps_env_4_binary_plk_bed
#plink --noweb --file Machine_Learning/Data/sig_snps_env_5_binary_plk --make-bed --out Machine_Learning/Data/sig_snps_env_5_binary_plk_bed