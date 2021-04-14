#!/bin/bash

# Extract the first column from the trait dataset (taxa) and save it as a file to be used in the gwas splitter later.
for i in {1..5}; do 
	cut -f1 -d',' MyY_Env_${i}.csv > MyY_Taxa_Env_${i}.txt;
done

# Extract the first column (taxa) from each of the numerical datasets and trait datasets and look for common lines. 
for i in {1..5}; do 
	cut -f1 GD_Env_${i}.txt > GD_Taxa_Env_${i}.txt; 
	cut -f1 -d',' MyY_Env_${i}.csv > Y_Taxa_Env_${i}.txt;
	comm -12 Y_Taxa_Env_${i}.txt GD_Taxa_Env_${i}.txt > Common_Taxa_Env_${i}.txt;
	comm -13 Y_Taxa_Env_${i}.txt GD_Taxa_Env_${i}.txt > GD_Only_Taxa_Env_${i}.txt;
	comm -23 Y_Taxa_Env_${i}.txt GD_Taxa_Env_${i}.txt > Trait_Only_Taxa_Env_${i}.txt;
done

# Check the length of the datasets created above (all three should be equal)
for i in {1..5}; do 
	echo Environment  $i
	wc -l GD_Taxa_Env_$i.txt; 
	wc -l Y_Taxa_Env_$i.txt;
	wc -l Common_Taxa_Env_$i.txt;
done

# Check dimensions of numerical datasets. Head line quickly checks for the number of snps since each snp is formatted as SNN_NNNNNNNN and only contains one _. This does not count the taxa in the final count.
for i in {1..5}; do
	echo Environment $i
	echo Columns:
	head -n1 GD_Env_${i}.txt | grep -o "_" | wc -l
	echo Rows:
	wc -l GD_Env_${i}.txt
done

# Check across environments to see how many genotypes are considered to be in common
for i in {1..5}; do
for n in {1..5}; do
	echo Number of Common Lines Between Environments $i and $n
	comm -12 GD_Taxa_Env_$i.txt GD_Taxa_Env_$n.txt | wc -l
done
done

# Check across environments to see how many lines are in common.  In theory it should be close (if not equal to) the size of the smaller dataset
for i in {1..5}; do
for n in {1..5}; do
	echo Number of Common Lines Between Environments $i and $n
	comm -12 GD_Taxa_Env_$i.txt GD_Taxa_Env_$n.txt | wc -l
	comm -12 GD_Env_${i}.txt GD_Env_${n}.txt | wc -l
done
done

# Since it appears something is misaligned between the different environments in the numerical datasets, repeat the loop above on the hapmap dataset.


