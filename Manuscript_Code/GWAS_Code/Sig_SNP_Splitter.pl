#! /bin/perl

# Perl based SNP filter to remove all significant SNPs 

####################
# Starting Options #
####################
use strict;
use warnings;
use Getopt::Std;
use Data::Dumper;

####################
# Usage Statements #
####################
my $usage = "\n$0 -s <LIST OF SIGNIFICANT SNPs> -d <HAPMAP DATASET> -o <OUTPUT FILE> -h <help>\n\n"; # statement to print if something goes wrong
our ($opt_s, $opt_d, $opt_o, $opt_h); # defining options
getopts("s:d:o:h") or die "\nGetopts Problem: $usage\n"; # requiring our options

if((!(defined $opt_s)) or (!(defined $opt_d)) or (!(defined $opt_o)) or (defined $opt_h)) { # checking to make sure all needed options are defined
	die "\n Definition Problem: $usage\n"; # stop and print usage statement if something goes wrong
}

###################
# Open File Paths #
###################
open(my $list_fh, '<', $opt_s) or die("Significant SNPs File Not Found"); # open significant snps list file path
open(my $hapmap_fh, '<', $opt_d) or die("HapMap Dataset Not Found"); # open hapmap dataset file path
open(my $output_fh, '>', $opt_o) or die("Output File Path Could Not Be Created"); # open output file path

my %hash;
<$list_fh>;
while (my $line = <$list_fh>){
	chomp $line;
	#print $line, "\n";
	$hash{$line} = '';
}

#print Dumper(%hash);

my $header = <$hapmap_fh>;
print $output_fh "$header";

while(my $line = <$hapmap_fh>){
	chomp $line;
	my @array = split('\t', $line);
	my $id = $array[0];
	if(defined $hash{$id}){
		print $output_fh "$line\n";
	}
}

close $list_fh;
close $hapmap_fh;
close $output_fh;
exit;