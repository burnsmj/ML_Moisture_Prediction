#! /bin/bash
# All SNPs
for i in {1..5}; do 
./gcta64 --bfile Data/christine_common_final_binary_env_${i} --make-grm --thread-num 10 --out Data/gcta_grm_out_env_$i 
done

# Significant SNPs
# For running significant SNPs
#./gcta64 --bfile Data/sig_snps_env_1_binary_plk_bed --make-grm --thread-num 10 --out Data/sig_snps_env_1_grm
#./gcta64 --bfile Data/sig_snps_env_2_binary_plk_bed --make-grm --thread-num 10 --out Data/sig_snps_env_2_grm
#./gcta64 --bfile Data/sig_snps_env_4_binary_plk_bed --make-grm --thread-num 10 --out Data/sig_snps_env_4_grm
#./gcta64 --bfile Data/sig_snps_env_5_binary_plk_bed --make-grm --thread-num 10 --out Data/sig_snps_env_5_grm
