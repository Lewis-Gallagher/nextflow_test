#!/usr/bin/perl

##########################################################################
#                                                                        #
# Copyright 2017-02, Lina Yuan (Lina.Yuan@icr.ac.uk)        #
#                                                                        #
# This file is part of MDx.                                           #
#                                                                        #
# MDx is free software: you can redistribute it and/or modify         #
# it under the terms of the GNU General Public License as published by   #
# the Free Software Foundation, either version 3 of the License, or      #
# (at your option) any later version.                                    #
#                                                                        #
# MDx is distributed in the hope that it will be useful,              #
# but WITHOUT ANY WARRANTY; without even the implied warranty of         #
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the          #
# GNU General Public License for more details.                           #
#                                                                        #
# You should have received a copy of the GNU General Public License      #
# along with MDx.  If not, see <http://www.gnu.org/licenses/>.        #
##########################################################################

#use warnings;
use strict;
use Getopt::Long;
use List::MoreUtils qw(uniq);


my $infile='';
my $report='';
my $overall='';
my $samp='';
my $insuffix='.base.coverage';
my $reportsuffix='.gene-report.tsv';
my $oversuffix='.overall_coverage.txt';
my $dir='';

my $tc=0.5; ### default tumour content
my $mindepth=40;### default min coverage 40x
my $pcoverage=0.85; ### default min percentage of coverage = 0.85
my $hcut=0.95; ### default cut off of 'High Confidence'
my $mcut=0.85; ### default cut off of 'Medium Confidence', below 0.85 defined as 'Low Confidence'
my $depth=0;

my @genes=();
my %ROI=();
my %fexon=();
my %mindepth=();
my %meandepth=();
my %percentcoverage=();
my $panelcoverage=0;

my $total=0;
my $fail=0;


my $args = scalar(@ARGV);
my $help=0;
my $result = GetOptions (
  "dir=s" => \$dir,
  "samp=s" => \$samp,
  "insuffix=s" => \$insuffix,
  "reportsuffix=s" => \$reportsuffix,
  "oversuffix=s" => \$oversuffix,
  "infile=s" => \$infile,
  "tcontent=s" => \$tc,
  "help" => \$help,
  );

# print usage message if requested or no args supplied
if($help || !$args) {
  &usage;
  exit(0);
}

if($infile eq ''){

	$infile=$dir.'/'.$samp.$insuffix;

}
#print STDERR $infile, "\n";

$report=$dir.'/'.$samp.$reportsuffix;
$overall=$dir.'/'.$samp.$oversuffix;
$depth=int($mindepth/$tc);

my $wholeexon='';

## new modification to handle qc regions within gene!
my %genes=(); ## to store all exons related to a gene
my %mID=();
my %genelevelstats=(); ## to store all gene level stats
my %exonlevelstats=(); ## to store all exon level stats



open IN, "<$infile";


while(<IN>){
	
	next if(/^chrom/);
	chomp;
	my @data=split(/\t/);
	my $molecularID=0;

	if($data[2] =~/mut/ || $data[2] =~/\w+_\d/){
    	next if($data[2] =~/_cn_/); ## remove CNV regions
    	$molecularID=1 if($data[2] =~/_qc_/);
		$total++;
		$fail++ if($data[3] < $depth);

	## dealing with overlapped regions

		if($data[2]=~m/\|/){
      #print STDERR $data[2], "\n";
			my @minfor=split('\|', $data[2]);
			my $keep='';
			foreach my $m(@minfor){
	
				if($m =~/mut/ || $m =~/\w+_\d/){
					$keep =$m if($keep eq '');
					if($keep ne ''){

						my @d=split('_', $m);
						$keep=$keep.'-'.$d[$#d];

					}
				}
      		}
			$data[2]=$keep;
    	}
		my @infor=split('_', $data[2]);
	# gene level
		if(!defined($genes{$infor[0]})){
			$mID{$infor[0]}=1 if($molecularID);
			$genelevelstats{$infor[0]}{'fgn'}=0;
			$genelevelstats{$infor[0]}{'mindepth'}=100000;
			$genelevelstats{$infor[0]}{'meandepth'}=0;
			$genes{$infor[0]} = 'none';
		}
	# exon level
		if(!defined($exonlevelstats{$data[2]}{'fen'})){
			
			$exonlevelstats{$data[2]}{'fen'}=0;
			$exonlevelstats{$data[2]}{'exon'}=$infor[$#infor];
			
			if($genes{$infor[0]} eq 'none'){
				
				$genes{$infor[0]}=$data[2];
			} else {
			
				$genes{$infor[0]}=$genes{$infor[0]}.','.$data[2];
			}
		}
		$exonlevelstats{$data[2]}{'en'}++;
		$exonlevelstats{$data[2]}{'fen'}++ if($data[3] < $depth);
		$genelevelstats{$infor[0]}{'gn'}++;
		$genelevelstats{$infor[0]}{'fgn'}++ if($data[3] < $depth);
		$genelevelstats{$infor[0]}{'mindepth'}=$data[3] if($genelevelstats{$infor[0]}{'mindepth'}>$data[3]);
		$genelevelstats{$infor[0]}{'meandepth'}=$genelevelstats{$infor[0]}{'meandepth'}+$data[3];
		
	}
}

close IN;


## calculation steps


foreach my $gene (keys %genes){
	my @gexons=split(/,/, $genes{$gene});
	my @fexons=();
	my @exons=();
	## at exon level
	foreach my $gexon (@gexons){
	
		my $pfen=1-$exonlevelstats{$gexon}{'fen'}/$exonlevelstats{$gexon}{'en'};
		push(@exons, $exonlevelstats{$gexon}{'exon'});
		push(@fexons, $exonlevelstats{$gexon}{'exon'}) if($pfen < $pcoverage); ## check if to push in exon name or just the number!
			
	}
	$genelevelstats{$gene}{'meandepth'}=sprintf("%d", $genelevelstats{$gene}{'meandepth'}/$genelevelstats{$gene}{'gn'});
	$genelevelstats{$gene}{'percentcoverage'}=sprintf("%.2f",1-$genelevelstats{$gene}{'fgn'}/$genelevelstats{$gene}{'gn'});
	$genelevelstats{$gene}{'fexon'}='None' if($#fexons == -1);
	$genelevelstats{$gene}{'fexon'}=get_ROI('Exon', join('-', @fexons))if($#fexons >0);
	$genelevelstats{$gene}{'fexon'}='Exon '.$fexons[0] if($#fexons == 0);
	$genelevelstats{$gene}{'ROI'}=get_ROI('Exon', join('-', @exons))if($#exons >0);
	$genelevelstats{$gene}{'ROI'}='Exon '.$exons[0] if($#exons == 0);

}


open OUT, ">$overall";
$panelcoverage= 1- $fail/$total;
print OUT sprintf("%.2f", $panelcoverage*100), "\n";
close OUT;

open OUT, ">$report";
print OUT "Gene\tGene_ROI\tMin\tMean\tx$depth\tExon below $pcoverage covered\tGene_report\n";

## output gene coverage first
foreach my $gene (sort keys %genes){
	next if(defined($mID{$gene}));
	my $greport='Low Confidence';

	$greport='Medium Confidence' if($genelevelstats{$gene}{'percentcoverage'} > $mcut);
	$greport='High Confidence' if($genelevelstats{$gene}{'percentcoverage'} > $hcut);

	print OUT "$gene\t$genelevelstats{$gene}{'ROI'}\t$genelevelstats{$gene}{'mindepth'}\t$genelevelstats{$gene}{'meandepth'}\t$genelevelstats{$gene}{'percentcoverage'}\t$genelevelstats{$gene}{'fexon'}\t$greport\n";

}
## output mID markers

foreach my $gene (sort keys %mID){

	my $greport='Low Confidence';

	$greport='Medium Confidence' if($genelevelstats{$gene}{'percentcoverage'} > $mcut);
	$greport='High Confidence' if($genelevelstats{$gene}{'percentcoverage'} > $hcut);
	$genelevelstats{$gene}{'fexon'} =~s/Exon/Marker/;
	$genelevelstats{$gene}{'ROI'} =~s/Exon/Marker/;
	print OUT "$gene\t$genelevelstats{$gene}{'ROI'}\t$genelevelstats{$gene}{'mindepth'}\t$genelevelstats{$gene}{'meandepth'}\t$genelevelstats{$gene}{'percentcoverage'}\t$genelevelstats{$gene}{'fexon'}\t$greport\n";

}


close OUT;


sub get_ROI{
	my ($label, $ROI_string)=@_;
	#$ROI_string =~ s/[a-zA-Z]//g;
	my @data=split('-', $ROI_string);
	my @data_uniq= uniq @data;
	if($#data >0){

		$label=$label.'s';
	}
	my @data_ordered=sort { $a <=> $b || $a cmp $b } @data_uniq;
	my $ROI=$label.' '.$data_ordered[0];
	my $k1=$data_ordered[0]+1;
	my $k2=$data_ordered[0]+1;

	foreach my $i (1..$#data_ordered){
		if($data_ordered[$i] == $k2){
			$k2 = $data_ordered[$i]+1;
			next;
		} else {
			if($k1 == $k2){
				$ROI=$ROI.','.$data_ordered[$i];
			} else {
				$ROI=$ROI.'-'.($k2-1).','.$data_ordered[$i];
			}
			$k1=$data_ordered[$i]+1;
			$k2=$data_ordered[$i]+1;
		}
	}
	$ROI=$ROI.'-'.($k2 -1) if($k1 != $k2);
	return($ROI);
}



sub usage() {
        my $usage =<<END;

VCFToTable.pl

Usage: perl VCFToTable_A06.pl [options]

Options

 To fill later

END
        print $usage;
}



exit;
