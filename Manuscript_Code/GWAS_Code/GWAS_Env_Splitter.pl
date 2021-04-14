#!/bin/perl
### GWAS Environment Splitter
### Michael Burns
### 2/22/21

# The purpose of this code is to read in two files (environment taxa and numericalized hapmap dataset) and extract from the hapmap data the genotypes found in the environment taxa list

#################
# Starter Stuff #
#################
use strict;
use warnings;
use Getopt::Std;

my $i; # Initialize iteration variable
for ($i=1; $i<=5; ++$i){ # For i in 1 through 5 (adding 1 to i each iteration) do the following...
	print "Currently Working on Environment $i\n"; # Tell me which environment is in progress 

	######################
	# Opening File Paths #
	######################
	open(my $Env_in_fh, '<', "MyY_Taxa_Env_$i.txt") or die("Env Dataset Couldn't be Found\n"); # Open path to environment taxa
	open(my $GD_in_fh, '<', "GAPIT.Genotype.Numerical.txt") or die("GD Dataset Couldn't be Found\n"); # Open path to genotypic dataset -- Change this later!
	open(my $Env_out_fh, '>', "GD_Env_$i.txt") or die("Output Dataset Couldn't be Found\n"); # Open path to output file

	###################
	# Populating Hash #
	###################
	my %env_hash; # Initialize hash
	while(my $line = <$GD_in_fh>){ # While there are lines in the genotype dataset...
		chomp $line; # Remove tailing \n  
		my @fields = split('\t', $line); # Split first line into separate elements by tabs
		my $taxa = $fields[0]; # Extract the first column (taxa) from the genotypic dataset
		my $num_cols = scalar @fields-1; # This should count the number of columns present, and subtract one. This is important since the indexing in perl starts at 0.
		my $geno = join "\t", @fields[1..$num_cols]; # Rejoin all of the genotypic data columns -- This will not throw an error if the wrong number of fields are selected, so be careful!
		$env_hash{$taxa} = $geno; # Add taxa key and genotypic values to hash -- all genotypic data needs to be in one element
	}

	########################################################
	# Matching Keys and Values in Environment Taxa Dataset #
	########################################################
	while(my $key = <$Env_in_fh>){ # While there are lines in environment taxa dataset...
		chomp $key; # Remove tailing \n 
		print $Env_out_fh $key, "\t", $env_hash{$key}, "\n"; # Print matches in order of Key    Numerical Data \n
	}

	####################
	# Close File Paths #
	####################
	close $Env_in_fh;
	close $GD_in_fh;
	close $Env_out_fh;
}