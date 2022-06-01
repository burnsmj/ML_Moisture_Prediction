#! /bin/bash

for i in {1..5}; do
	echo Working on Environment ${i}
	perl ../../Perl_Scripts/Sig_SNP_Splitter.pl -s Env_${i}_Sig_SNPs.txt -d christine_final_env_${i}_binary.hmp.txt -o sig_snps_filter_env_${i}.hmp.txt
	echo Environment ${i} is length:
	wc -l sig_snps_filter_env_${i}.hmp.txt
done