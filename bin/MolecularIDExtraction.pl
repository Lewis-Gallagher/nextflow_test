#!/usr/bin/perl -w
# extract ref and alt reads from a VCF file
# write out the VCF file with columns appended
# will also work with any VCF-like (e.g. merge) file
# James Campbell, 26th Jul 2012
# version 1.0

use strict;
use Getopt::Long;
use List::Util qw(sum);

$| = 1;

my $vcf="";
my $bed = "";
my $out = undef;
my $args = scalar(@ARGV);
my $help;
my %snps=();
my @keys=();
my $result = GetOptions (
  "infile=s" => \$vcf,
  "bed=s" => \$bed,
  "out=s" => \$out,
  "help" => \$help,
  );



# print usage message if requested or no args supplied
if(defined($help) || !$args) {
  &usage;
  exit(0);
}

$out = $vcf . '.table' unless defined($out);

open BED, "< $bed" or die "Can't read bed input file $bed: $!\n";
while(<BED>){

	next if (/^@/);
	chomp;
	my @data=split(/\t/);
	next if($data[4] !~/^rs/); ## working on Molecular ID SNPs not hotspots
	$data[0]=~m/chr(\d+)/;
	my $key=$1*10000000000+$data[1];
	$snps{$key}{name}=$data[4];
	push(@keys, $key);
	
}

close BED;
open VCF, "< $vcf" or die "Can't read input file $vcf: $!\n";


## Vardict output is not in VCF format! 
## columns infor: 0: sample name, 2: chr; 3:start; 4:end; 5: refallele; 6:Altallele; 7:depth; 8: altdepth; 13: genotype; 14: allele frequency;
## https://github.com/AstraZeneca-NGS/VarDict/wiki/VarDict-Raw-Output

while(<VCF>){
	next if (/^#/);
	## reset all parameters
	my $depth=0;
	my @infors=();
	my $genotype='./.';
	my @adp=();
	my $vaf=0;
	my $gcode=9; ## 9 for failed genotyping; 0=wildtype; 1=het; 2=hom
	
	chomp;
	my @data=split(/\t/);
	
	## filter out anything not in bed file!!!
	$data[2]=~m/chr(\d+)/;
	my $key=$1*10000000000+$data[3];
	next if(!defined($snps{$key}{name}));
	$key=$1*10000000000+$data[4];
	next if(!defined($snps{$key}{name}));
	
	$depth=$data[7];
	$vaf=$data[14];
	next if(defined($snps{$key}{VAF}) && $snps{$key}{VAF} > $vaf);
	
	$snps{$key}{VAF}=$vaf;
	if($vaf < 0.05){
		$gcode=0;
		$genotype='0/0';
	
	}elsif($vaf > 0.90){
		$gcode=2;
		$genotype='1/1';
	
	} else{
		$gcode=1;
		$genotype='0/1';
	
	}
	
	if($depth < 30){
		$gcode=9;
		$genotype='./.';
	
	}
	$snps{$key}{output}="$snps{$key}{name}\t$data[2]:$data[3]\t$data[5]\t$data[6]\t$depth\t$vaf\t$genotype\t$gcode\n";
	
}

close VCF;

open OUT, "> $out" or die "Can't write output file $out: $!\n";

foreach my $key (sort {$a <=> $b} @keys){
	
	print OUT $snps{$key}{output};
	
	
}

close OUT;


sub usage() {
        my $usage =<<END;
#------------------------#
#  get_allele_freqs.pl   #
#------------------------#

James Campbell (james.campbell\@icr.ac.uk); Lina Yuan (lina.yuan\@icr.ac.uk)

Usage:
perl get_allele_freqs.pl [options]

Options
--help          Display this message and quit
--vcf           Path to where the vcf file is located           [required]
--columns	Defaults to (zero-based) all samples in the VCF. List of comma-separated integers
--out           Path to output file                             [default: ./inFile.classes]

END

        print $usage;
}


