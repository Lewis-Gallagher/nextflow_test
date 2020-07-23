#!/usr/bin/perl -w ###################################################
#                                                                    #
# qc.pl                                                              #
# script used to qc results from aligned bam files                   #
# 07-11-15 wei.yuan@icr.ac.uk - v1.1 add in amplicon analysis        #
# 08-04-18 Lina.yuan@icr.ac.uk - v1.3 modify to suit new pipeline        #
######################################################################



use strict;
use Getopt::Long;

# variables needed:
my $output_script = '';				# path to write the bwa job script
my $ref_fasta='';					# path to the fasta file
my $analysis_dir = '';				# path to analysis folder
my $samp_id = '';					# sample ID
my $java_tmp = '';					# temp folder for java
my $picard = '';					# path to picard
my $getexoncoverage = '';					# path to get_exon_covearge.pl
my $interval='';					# path to interval file (capture bed file with genome header)
my $seq_type=''; 					# sequencing type: exome (capture), wgs, amplicon, rna et al
my $tcontent=0.5;
my $mid_file='';

### temp fix hotspots file

my $hotspots='/scratch/DMP/DUDMP/TRANSGEN/transgen-mdx/ngs/7.resources/hg19/bedfiles/hotspots.v1.0.bed.intervals';

my $args = scalar(@ARGV);
my $help;
my $result = GetOptions (
  "output_script=s" => \$output_script,
  "java_tmp=s" => \$java_tmp,
  "picard=s" => \$picard,
  "getexoncoverage=s" => \$getexoncoverage,
  "ref_fasta=s" => \$ref_fasta,
  "analysis_dir=s" => \$analysis_dir,
  "samp_id=s" => \$samp_id,
  "interval_file=s" => \$interval,
  "MID_file=s" => \$mid_file,
  "seq_type=s" => \$seq_type,
  "tcontent=s" => \$tcontent,
  "help" => \$help,
  );

# print usage message if requested or no args supplied
if(defined($help) || !$args) {
  &usage;
  exit(0);
}


my $bed_interval = $interval;
$bed_interval=~s/\.CoverageCalculator//;

$seq_type = lc $seq_type;

my $commands='';
if($seq_type eq 'capture' || $seq_type eq 'exome'){
#if($seq_type eq 'exome' || $seq_type eq "capture"){

	$commands =<<END_OF_OUTPUT_SCRIPT_pt1;
#!/bin/bash

java -Xmx10g -XX:ParallelGCThreads=1 -jar $picard CalculateHsMetrics MAX_RECORDS_IN_RAM=500000 \\
REFERENCE_SEQUENCE=$ref_fasta \\
INPUT=$samp_id.sort.bam \\
OUTPUT=$samp_id.qc.metrics \\
BAIT_INTERVALS=$interval \\
TARGET_INTERVALS=$interval \\
BAIT_SET_NAME=$samp_id \\
METRIC_ACCUMULATION_LEVEL=ALL_READS \\
TMP_DIR=$java_tmp \\
PER_TARGET_COVERAGE=$samp_id.qc.coverage \\
PER_BASE_COVERAGE=$samp_id.base.coverage \\
VALIDATION_STRINGENCY=LENIENT

java -Xmx10g -XX:ParallelGCThreads=1 -jar $picard CalculateHsMetrics MAX_RECORDS_IN_RAM=500000 \\
REFERENCE_SEQUENCE=$ref_fasta \\
INPUT=$samp_id.sort.bam \\
OUTPUT=$samp_id.qc.metrics.targetbed \\
BAIT_INTERVALS=$bed_interval \\
TARGET_INTERVALS=$bed_interval \\
BAIT_SET_NAME=$samp_id \\
METRIC_ACCUMULATION_LEVEL=ALL_READS \\
TMP_DIR=$java_tmp \\
PER_TARGET_COVERAGE=$samp_id.qc.coverage.targetbed \\
PER_BASE_COVERAGE=$samp_id.base.coverage.targetbed \\
VALIDATION_STRINGENCY=LENIENT


perl $getexoncoverage --dir . --samp $samp_id --tcontent $tcontent --insuffix .base.coverage.targetbed

END_OF_OUTPUT_SCRIPT_pt1

} elsif ($seq_type eq 'amplicon' || $seq_type eq 'umi' || $seq_type eq 'cfdna'){
  $commands =<<END_OF_OUTPUT_SCRIPT_pt2;
#!/bin/bash

java -Xmx10g -XX:ParallelGCThreads=1 -jar $picard CollectTargetedPcrMetrics MAX_RECORDS_IN_RAM=500000 \\
REFERENCE_SEQUENCE=$ref_fasta \\
INPUT=$samp_id.sort.bam \\
OUTPUT=$samp_id.qc.metrics \\
AMPLICON_INTERVALS=$interval \\
TARGET_INTERVALS=$interval \\
CUSTOM_AMPLICON_SET_NAME=$samp_id \\
METRIC_ACCUMULATION_LEVEL=ALL_READS \\
TMP_DIR=$java_tmp \\
PER_TARGET_COVERAGE=$samp_id.qc.coverage \\
PER_BASE_COVERAGE=$samp_id.base.coverage \\
VALIDATION_STRINGENCY=LENIENT

perl $getexoncoverage --dir . --samp $samp_id --tcontent $tcontent

END_OF_OUTPUT_SCRIPT_pt2

} else {
	$commands =<<END_OF_OUTPUT_SCRIPT_pt1;

#!/bin/bash

java -Xmx10g -XX:ParallelGCThreads=1 -jar $picard CollectWgsMetrics MAX_RECORDS_IN_RAM=500000 \\
REFERENCE_SEQUENCE=$ref_fasta \\
INPUT=$samp_id.sort.bam \\
OUTPUT=$samp_id.qc.metrics \\
TMP_DIR=$java_tmp \\
VALIDATION_STRINGENCY=LENIENT

perl $getexoncoverage --dir . --samp $samp_id --tcontent $tcontent

END_OF_OUTPUT_SCRIPT_pt1


}

my $mid_cmd='';

## fix structure variant error from VarDict for cfDNA samples -- 2019.12.10
# add in -U 


$mid_cmd =<<END_OF_OUTPUT_SCRIPT_pt2;

VarDict -p -U -G $ref_fasta  -N $samp_id -K \\
    -b $samp_id.sort.bam \\
    -z -c 1 -S 2 -E 3 \\
    $mid_file  > $samp_id.mid.vardict

perl MolecularIDExtraction.pl --bed $mid_file --infile $samp_id.mid.vardict

END_OF_OUTPUT_SCRIPT_pt2

## extra QC for all sequencing types

my $cmd2='';

$cmd2=<<END_OF_OUTPUT_SCRIPT_pt3;

java -Xmx10g -XX:ParallelGCThreads=1 -jar $picard CalculateHsMetrics MAX_RECORDS_IN_RAM=500000 \\
REFERENCE_SEQUENCE=$ref_fasta \\
INPUT=$samp_id.sort.bam \\
OUTPUT=$samp_id.qc.metrics.hotspots \\
BAIT_INTERVALS=$hotspots \\
TARGET_INTERVALS=$hotspots \\
BAIT_SET_NAME=$samp_id \\
METRIC_ACCUMULATION_LEVEL=ALL_READS \\
TMP_DIR=$java_tmp \\
PER_TARGET_COVERAGE=$samp_id.qc.coverage.hotspots \\
PER_BASE_COVERAGE=$samp_id.base.coverage.hotspots \\
VALIDATION_STRINGENCY=LENIENT

samtools stats $samp_id.sort.bam > $samp_id.qc.samtools.stats
plot-bamstats -p $samp_id/ $samp_id.qc.samtools.stats


END_OF_OUTPUT_SCRIPT_pt3



open OUT, "> $output_script" or die "can't write script to run bwa in $output_script: $!\n\n";
print OUT $commands;
print OUT $cmd2;
print OUT $mid_cmd if(defined($mid_file) && -e $mid_file);
close OUT;

chmod 0755, $output_script or warn "Could not set permissions for file $output_script: $!\n\n";

sub usage() {
        my $usage =<<END;

qc.pl

Usage: perl qc.pl [options]

Options

 --java_tmp						Path to Java tmp folder
 --picard						Path to picard
 --interval_file			Path to interval file (capture bed file with genome header)
 --analysis_dir			Path to analysis directory
 --ref_fasta			Path to the reference fasta file
 --samp_id			Sample ID
 --help				Display this message and quit
 --output_script		Path to write script to

Updated: 11 Nov 2015 by Wei.Yuan\@icr.ac.uk
Updated: 18 Apr 2018 by Lina.Yuan\@icr.ac.uk
Updated: 13 March 2019 by Lina.Yuan\@icr.ac.uk

END
        print $usage;
}
