#!/usr/bin/perl -w ###################################################
#                                                                    #
# bwa.pl                                                             #
# script used to configure and run bwa alignments                    #
# recieves options needed to configure a bwa                         #
# write a shell script in the analysis/proj/samp dir                 #
#
#Wei.Yuan@icr.ac.uk Sep 2015                                         #
#   Update Lina.Yuan@icr.ac.uk Apr 2018                                                                 #
# change run_analysis logic to solve issue of mix pool with UMI and none-UMI samples -- 2019.12.13 Lina
######################################################################


use strict;
use Getopt::Long;

# variables needed:
my $output_script = '';				# path to write the bwa job script
my $num_cpu = 4; 					# default number of cpus
my $ref_fasta = '';					# path to reference fasta
my $analysis_dir = '';				# path to analysis folder
my $fastq_dir = '';					# path to fastq folder
my $samp_id = '';					# sample ID
my $seq_plat = 'Illumina';				# platform of sequencing, default:illumina
my $barcode_seq = 'NNNNNNN'; 			# barcode seq for RGPU=CTTGTA
my $java_tmp = '';					# temp folder for java
my $rm_dups='true';					# remove PCR duplicates
my $java_mem=32;					# default memory for java.
my $java_cpu=2;						# default CPU numbers for java.
my $picard='';						# path to picard
my $umi='none';						# $umi taking information from individual UMI information: none, IDT, Q33, Avenio
my $seq_type='capture';						# $umi taking information from Seq_type: capture,amplicon, wgs, lcwgs, umi or cfDNA
my $run_type ='none'; ## 0 for non-umi run, 1 for umi run;
my $args = scalar(@ARGV);
my $help;
my $umiqc='/scratch/DMP/DUDMP/TRANSGEN/lchen/ngs/9.Scripts/1.pipeline_modules/1.aligner/UMIQC.pl';
my $filetransfer='';
my $expected_umi_len=9; ## length of UMI: IDT=9; Q33=12; AVENIO=6+6=12
my $avenio_rstr='6M+T 6M+T'; ## Aveno read structure  --read-structure: B=Barcode, T=Template, M=Molecular Identifier) 
my $Q33_rstr='+T 12M12S+T'; ## Q33 read structure --read-structure: B=Barcode, T=Template, M=Molecular Identifier) 
my $seed =50; # default  bwa seed is 50bp

my $result = GetOptions (
  "output_script=s" => \$output_script,
  "java_tmp=s" => \$java_tmp,
  "picard=s" => \$picard,
  "num_cpu=i" => \$num_cpu,
  "ref_fasta=s" => \$ref_fasta,
  "analysis_dir=s" => \$analysis_dir,
  "fastq_dir=s" => \$fastq_dir,
  "samp_id=s" => \$samp_id,
  "seq_plat=s" => \$seq_plat,
  "seed=i" => \$seed,
  "barcode_seq=s" => \$barcode_seq,
  "rm_dups=s" => \$rm_dups,
  "run_type=s" => \$run_type,
  "seq_type=s" => \$seq_type,
  "umi=s" => \$umi,
  "help" => \$help,
  );

# print usage message if requested or no args supplied
if(defined($help) || !$args) {
  &usage;
  exit(0);
}

my $transfer_dir=$analysis_dir;
$transfer_dir=~s/3\.analysis/001\.reports/ ;



my $fq1_cmd = "\'<zcat $fastq_dir/$samp_id\*_R1_\*.gz\'";
my $fq2_cmd = "\'<zcat $fastq_dir/$samp_id\*_R2_\*.gz\'";
my $umi_seq = "$fastq_dir/$samp_id\*_R2_\*.gz";
my $rstr=$Q33_rstr; ## default is Q33;
$rstr=$avenio_rstr if($run_type eq 'AVENIO');

#$fq2_cmd = "\'<zcat $fastq_dir/$samp_id\*_R3_\*.gz\'" if($run_type ne 'none' && $run_type ne 'AVENIO');
$fq2_cmd = "\'<zcat $fastq_dir/$samp_id\*_R3_\*.gz\'" if($run_type eq 'IDT');


## increase CPUs for cfDNA test ...

#if($java_cpu >= $num_cpu){
#	$java_cpu =1;

#}

$java_cpu =$num_cpu-1;

#my $java_mem_loc = $num_cpu * 10; ## 16g per CPU on zygon and Cyberman, 12g for dalek, 12.5g for davros

#if($java_mem > $java_mem_loc){
#	$java_mem =$java_mem_loc;

#}
$java_mem = $num_cpu * 12;



my $RG=" \'\@RG\tID:$samp_id\tPL:$seq_plat\tPU:$barcode_seq\tLB:$samp_id\tSM:$samp_id\'";


#my $RG_stampy="ID:$samp_id,PL:$seq_plat,PU:$barcode_seq,LB:$samp_id,SM:$samp_id";


my $command1 =<<END_OF_OUTPUT_SCRIPT_pt1;
#!/bin/bash

bwa mem -M -c 1 -t $num_cpu -k $seed -R $RG \\
$ref_fasta $fq1_cmd $fq2_cmd \\
> $analysis_dir/$samp_id.fq.aln.sam


java -Xmx${java_mem}g -XX:ParallelGCThreads=$java_cpu -jar $picard SortSam MAX_RECORDS_IN_RAM=4000000 \\
SORT_ORDER=coordinate \\
INPUT=$analysis_dir/$samp_id.fq.aln.sam \\
OUTPUT=$analysis_dir/$samp_id.fq.aln.sorted.bam \\
CREATE_INDEX=true \\
TMP_DIR=$java_tmp \\
VALIDATION_STRINGENCY=LENIENT

java -Xmx${java_mem}g -XX:ParallelGCThreads=$java_cpu -jar $picard MarkDuplicates MAX_RECORDS_IN_RAM=4000000 \\
REMOVE_DUPLICATES=$rm_dups \\
INPUT=$analysis_dir/$samp_id.fq.aln.sorted.bam \\
OUTPUT=$analysis_dir/$samp_id.fq.aln.sorted.dedup.bam \\
METRICS_FILE=$analysis_dir/$samp_id.mark_dups.metrics \\
ASSUME_SORTED=true \\
CREATE_INDEX=true \\
TMP_DIR=$java_tmp \\
VALIDATION_STRINGENCY=LENIENT

# Clip overlapping reads to prevent inflated depth in variant calling.
java -Xmx${java_mem}g -XX:+AggressiveOpts -XX:+AggressiveHeap  -jar /scratch/DMP/DUDMP/TRANSGEN/lchen/ngs/9.Scripts/third_party_tools/UMIs/fgbio-0.6.1.jar ClipBam \\
    -i $analysis_dir/$samp_id.fq.aln.sorted.dedup.bam  \\
    -o $analysis_dir/$samp_id.sort.bam \\
    -m $analysis_dir/$samp_id.clipbam.metrics \\
    -r $ref_fasta \\
    --clip-overlapping-reads=true \\
    -c Hard 


END_OF_OUTPUT_SCRIPT_pt1

my $commandumi1 =<<END_OF_OUTPUT_SCRIPT_pt2;
#!/bin/bash

fq1_seq="`ls $fastq_dir/$samp_id\*_R1_\*.gz`";
fq2_seq="`ls $fastq_dir/$samp_id\*_R3_\*.gz`";
UMISEQ=$umi_seq

# Convert fastq to bam and annotate with both UMIs by supplying --read-structure: B=Barcode, T=Template, M=Molecular Identifier) 
# Validate read structures here: https://fulcrumgenomics.github.io/fgbio/validate-read-structure.html
# worked for IDT (UMI in index)

java -Xmx${java_mem}g -XX:+AggressiveOpts -XX:+AggressiveHeap  -jar /scratch/DMP/DUDMP/TRANSGEN/lchen/ngs/9.Scripts/third_party_tools/UMIs/fgbio-0.6.1.jar FastqToBam \\
--input \$fq1_seq \$fq2_seq \$UMISEQ \\
--output $analysis_dir/$samp_id.fq.unaln.umi.bam \\
--read-structures +T +T +M --sort true --umi-tag RX \\
--sample $samp_id --library $samp_id --read-group-id $samp_id


END_OF_OUTPUT_SCRIPT_pt2

my $commandumi1a =<<END_OF_OUTPUT_SCRIPT_pt2a;
#!/bin/bash

fq1_seq="`ls $fastq_dir/$samp_id\*_R1_\*.gz`";
fq2_seq="`ls $fastq_dir/$samp_id\*_R2_\*.gz`";

## It may work very well for Q33 but it won't handle IDT!
# Convert fastq to bam and annotate with both UMIs by supplying --read-structure: B=Barcode, T=Template, M=Molecular Identifier) 
# Validate read structures here: https://fulcrumgenomics.github.io/fgbio/validate-read-structure.html
java -Xmx${java_mem}g -XX:+AggressiveOpts -XX:+AggressiveHeap  -jar /scratch/DMP/DUDMP/TRANSGEN/lchen/ngs/9.Scripts/third_party_tools/UMIs/fgbio-0.6.1.jar FastqToBam \\
--input \$fq1_seq \$fq2_seq \\
--output $analysis_dir/$samp_id.fq.unaln.umi.precheck.bam \\
--read-structures $rstr --sort true --umi-tag RX \\
--sample $samp_id --library $samp_id --read-group-id $samp_id

samtools view -h -o $analysis_dir/$samp_id.fq.unaln.umi.precheck.sam $analysis_dir/$samp_id.fq.unaln.umi.precheck.bam
sed 's/RX:Z:\\([A-Z]*\\)-\\([A-Z]*\\)/RX:Z:\\1\\2/' $analysis_dir/$samp_id.fq.unaln.umi.precheck.sam > $analysis_dir/$samp_id.fq.unaln.umi.sam
samtools view -S -b $analysis_dir/$samp_id.fq.unaln.umi.sam > $analysis_dir/$samp_id.fq.unaln.umi.bam

END_OF_OUTPUT_SCRIPT_pt2a


my $commandumi2 =<<END_OF_OUTPUT_SCRIPT_pt3;

# Remove Illumina sequence adaptors (picard)
java -Xmx${java_mem}g -XX:ParallelGCThreads=$java_cpu -jar $picard MarkIlluminaAdapters MAX_RECORDS_IN_RAM=4000000 \\
INPUT=$analysis_dir/$samp_id.fq.unaln.umi.bam \\
OUTPUT=$analysis_dir/$samp_id.fq.unaln.marked.umi.bam \\
TMP_DIR=/scratch/DMP/DUDMP/TRANSGEN/transgen-mdx/ngs//temp_dir \\
M=$analysis_dir/$samp_id.adaptor.metrics

# Convert BAMs to Fastq format for alignment (picard)
java -Xmx${java_mem}g -XX:ParallelGCThreads=$java_cpu -jar $picard SamToFastq MAX_RECORDS_IN_RAM=4000000 \\
INPUT=$analysis_dir/$samp_id.fq.unaln.marked.umi.bam  \\
FASTQ=$analysis_dir/$samp_id.fq.unaln.marked.umi.fq \\
CLIPPING_ATTRIBUTE=XT CLIPPING_ACTION=2 INTERLEAVE=true NON_PF=true \\
TMP_DIR=/scratch/DMP/DUDMP/TRANSGEN/transgen-mdx/ngs//temp_dir

# Align (bwa)
bwa mem -M -c 1 -t $num_cpu -k $seed -p \\
$ref_fasta $analysis_dir/$samp_id.fq.unaln.marked.umi.fq \\
> $analysis_dir/$samp_id.fq.aln.umi.sam

# Merge unaligned BAM and aligned BAM to preserve Tags from unaligned file (picard)

java -Xmx${java_mem}g -XX:ParallelGCThreads=$java_cpu -jar $picard MergeBamAlignment MAX_RECORDS_IN_RAM=4000000 \\
R=$ref_fasta \\
UNMAPPED_BAM=$analysis_dir/$samp_id.fq.unaln.marked.umi.bam \\
ALIGNED_BAM=$analysis_dir/$samp_id.fq.aln.umi.sam \\
O=$analysis_dir/$samp_id.fq.aln.merged.umi.bam \\
CREATE_INDEX=true ADD_MATE_CIGAR=true CLIP_ADAPTERS=false CLIP_OVERLAPPING_READS=true INCLUDE_SECONDARY_ALIGNMENTS=true \\
MAX_INSERTIONS_OR_DELETIONS=-1 \\
PRIMARY_ALIGNMENT_STRATEGY=MostDistant \\
ATTRIBUTES_TO_RETAIN=XS \\
TMP_DIR=/scratch/DMP/DUDMP/TRANSGEN/transgen-mdx/ngs//temp_dir


# Group reads by UMI (fgbio)

java -Xmx${java_mem}g -XX:+AggressiveOpts -XX:+AggressiveHeap  -jar /scratch/DMP/DUDMP/TRANSGEN/lchen/ngs/9.Scripts/third_party_tools/UMIs/fgbio-0.6.1.jar GroupReadsByUmi \\
-i $analysis_dir/$samp_id.fq.aln.merged.umi.bam \\
-f $analysis_dir/$samp_id.umi.metrics -s adjacency -m 30 -t RX -T MI --min-umi-length $expected_umi_len \\
-o $analysis_dir/$samp_id.fq.aln.merged.umi.grp.bam

END_OF_OUTPUT_SCRIPT_pt3


my $commandumi3 =<<END_OF_OUTPUT_SCRIPT_pt4;

java -Xmx${java_mem}g -XX:ParallelGCThreads=$java_cpu -jar $picard SortSam MAX_RECORDS_IN_RAM=4000000 \\
SORT_ORDER=coordinate \\
INPUT=$analysis_dir/$samp_id.fq.aln.merged.umi.grp.bam \\
OUTPUT=$analysis_dir/$samp_id.fq.aln.merged.umi.grp.sort.bam \\
CREATE_INDEX=true \\
TMP_DIR=$java_tmp \\
VALIDATION_STRINGENCY=LENIENT

java -Xmx${java_mem}g -XX:ParallelGCThreads=$java_cpu -jar $picard MarkDuplicates MAX_RECORDS_IN_RAM=4000000 \\
REMOVE_DUPLICATES=$rm_dups BARCODE_TAG=MI \\
INPUT=$analysis_dir/$samp_id.fq.aln.merged.umi.grp.sort.bam \\
OUTPUT=$analysis_dir/$samp_id.sort.bam \\
METRICS_FILE=$analysis_dir/$samp_id.mark_dups.metrics \\
ASSUME_SORTED=true \\
CREATE_INDEX=true \\
TMP_DIR=$java_tmp \\
VALIDATION_STRINGENCY=LENIENT

    
## Clip overlapping reads to prevent inflated depth in variant calling.
#java -Xmx${java_mem}g -XX:+AggressiveOpts -XX:+AggressiveHeap  -jar /scratch/DMP/DUDMP/TRANSGEN/lchen/ngs/9.Scripts/third_party_tools/UMIs/fgbio-0.6.1.jar ClipBam \\
 #   -i $analysis_dir/$samp_id.fq.aln.merged.umi.grp.sort.dedup.bam \\
  #  -o $analysis_dir/$samp_id.sort.bam \\
   # -m $analysis_dir/$samp_id.clipbam.metrics \\
    #-r $ref_fasta \\
  #  --clip-overlapping-reads=true \\
  #  -c Hard 



END_OF_OUTPUT_SCRIPT_pt4

my $commandumicfdna =<<END_OF_OUTPUT_SCRIPT_pt5;

# The next step produces consensus reads, making it unnecessary to duplicate-mark the reads or perform further processing of the raw data.
#


# BAM needs to be coordinate sorted to call consensus reads
java -Xmx${java_mem}g  -XX:+AggressiveOpts -XX:+AggressiveHeap  -jar /scratch/DMP/DUDMP/TRANSGEN/lchen/ngs/9.Scripts/third_party_tools/UMIs/fgbio-0.6.1.jar SortBam \\
    -i $analysis_dir/$samp_id.fq.aln.merged.umi.grp.bam \\
    -o $analysis_dir/$samp_id.fq.aln.merged.umi.grp.sort.bam \\
    --sort-order TemplateCoordinate \\
    --max-records-in-ram 4000000


# Call consensus reads (fgbio)
java -Xmx${java_mem}g -XX:+AggressiveOpts -XX:+AggressiveHeap  -jar /scratch/DMP/DUDMP/TRANSGEN/lchen/ngs/9.Scripts/third_party_tools/UMIs/fgbio-0.6.1.jar CallMolecularConsensusReads \\
    -i $analysis_dir/$samp_id.fq.aln.merged.umi.grp.sort.bam \\
    -o $analysis_dir/$samp_id.fq.aln.merged.umi.grp.sort.consensus.bam \\
    --min-reads 3 \\
    --min-input-base-quality 30 \\
    --tag MI


# Filter consensus reads - Error correction (fgbio)
# 15/05/2020: min-base-qual 30 -> 40;--max-no-call-fraction 0.1 -> 0.05
java -Xmx${java_mem}g -XX:+AggressiveOpts -XX:+AggressiveHeap  -jar /scratch/DMP/DUDMP/TRANSGEN/lchen/ngs/9.Scripts/third_party_tools/UMIs/fgbio-0.6.1.jar FilterConsensusReads \\
    -i $analysis_dir/$samp_id.fq.aln.merged.umi.grp.sort.consensus.bam \\
    -o $analysis_dir/$samp_id.fq.aln.merged.umi.grp.sort.consensus.filtered.bam \\
    -r $ref_fasta \\
    --min-reads 3 --max-read-error-rate 0.05 --min-base-quality 40 --max-base-error-rate 0.1 --max-no-call-fraction 0.05 \\
    --reverse-per-base-tags true


# Convert BAMs back to FASTQ for alignment (picard)
java -Xmx${java_mem}g -XX:ParallelGCThreads=$java_cpu -jar $picard SamToFastq MAX_RECORDS_IN_RAM=4000000 \\
INPUT=$analysis_dir/$samp_id.fq.aln.merged.umi.grp.sort.consensus.filtered.bam \\
FASTQ=$analysis_dir/$samp_id.fq.aln.merged.umi.grp.sort.consensus.filtered.fq \\
CLIPPING_ATTRIBUTE=XT CLIPPING_ACTION=2 INTERLEAVE=true NON_PF=true \\
TMP_DIR=$java_tmp


# Align consensus sequences (bwa)
bwa mem -M -c 1 -t 4 -k $seed -p \\
$ref_fasta $analysis_dir/$samp_id.fq.aln.merged.umi.grp.sort.consensus.filtered.fq \\
> $analysis_dir/$samp_id.fq.aln.merged.umi.grp.sort.consensus.filtered.aln.sam


# Aligneed and unaligned bam files need to be sorted in the same order to be merged...
# Sort aligned bam
java -Xmx${java_mem}g -XX:ParallelGCThreads=$java_cpu -jar $picard SortSam MAX_RECORDS_IN_RAM=4000000 \\
    SORT_ORDER=queryname \\
    INPUT=$analysis_dir/$samp_id.fq.aln.merged.umi.grp.sort.consensus.filtered.aln.sam \\
    OUTPUT=$analysis_dir/$samp_id.fq.aln.merged.umi.grp.sort.consensus.filtered.aln.sort.bam \\
    CREATE_INDEX=true \\
    TMP_DIR=$java_tmp \\
    VALIDATION_STRINGENCY=LENIENT


# resort consensus.filtered.bam bam
java -Xmx${java_mem}g -XX:ParallelGCThreads=$java_cpu -jar $picard SortSam MAX_RECORDS_IN_RAM=4000000 \\
    SORT_ORDER=queryname \\
    INPUT= $analysis_dir/$samp_id.fq.aln.merged.umi.grp.sort.consensus.filtered.bam \\
    OUTPUT=$analysis_dir/$samp_id.fq.aln.merged.umi.grp.sort.consensus.filtered.resort.bam \\
    CREATE_INDEX=true \\
    TMP_DIR=$java_tmp \\
    VALIDATION_STRINGENCY=LENIENT


## Merge aligned and unaligned bam files to retain RX tags
java -Xmx${java_mem}g -XX:ParallelGCThreads=$java_cpu -jar $picard MergeBamAlignment MAX_RECORDS_IN_RAM=4000000 \\
    R=$ref_fasta \\
    UNMAPPED_BAM=$analysis_dir/$samp_id.fq.aln.merged.umi.grp.sort.consensus.filtered.resort.bam \\
    ALIGNED_BAM=$analysis_dir/$samp_id.fq.aln.merged.umi.grp.sort.consensus.filtered.aln.sort.bam \\
    O=$analysis_dir/$samp_id.fq.aln.merged.umi.grp.sort.consensus.filtered.aln.sort.merged.bam \\
    CREATE_INDEX=true ADD_MATE_CIGAR=true CLIP_ADAPTERS=true CLIP_OVERLAPPING_READS=true INCLUDE_SECONDARY_ALIGNMENTS=true SO=coordinate \\
    MAX_INSERTIONS_OR_DELETIONS=-1 \\
    PRIMARY_ALIGNMENT_STRATEGY=MostDistant \\
    ATTRIBUTES_TO_RETAIN=XS \\
    TMP_DIR=$java_tmp
    
    
# Clip overlapping reads to prevent inflated depth in variant calling.
java -Xmx${java_mem}g -XX:+AggressiveOpts -XX:+AggressiveHeap  -jar /scratch/DMP/DUDMP/TRANSGEN/lchen/ngs/9.Scripts/third_party_tools/UMIs/fgbio-0.6.1.jar ClipBam \\
    -i $analysis_dir/$samp_id.fq.aln.merged.umi.grp.sort.consensus.filtered.aln.sort.merged.bam \\
    -o $analysis_dir/$samp_id.sort.bam \\
    -m $analysis_dir/$samp_id.clipbam.metrics \\
    -r $ref_fasta \\
    --clip-overlapping-reads=true \\
    -c Hard 
    

END_OF_OUTPUT_SCRIPT_pt5


open OUT, "> $output_script" or die "can't write script to run bwa in $output_script: $!\n\n";

	## different scenario
	## 1. IDT-UMI only -- solved
	## 2. Q33-UMI only -- solved
	## 3. IDT-UMI and IDT-UDI mix pool -- solved
	## 4. Q33-UMI and IDT-UDI mix pool -- solved
	## 5. IDT-UMI and Q33-UMI mix pool -- needs UMI infor perl sample base! add in new $umi ($seq_type replace old $umi)
	## 6. IDT-UDI (or any none UMI runs) only -- solved
	## 7. Avenio and IDT-UDI -- solved
	## 8. Avenio with IDT-UMI -- needs UMI info per sample base! add in new $umi ($seq_type replace old $umi)
	## To Note, Q33 and Avenio have their own UDIs

if($seq_type ne 'umi' && $seq_type ne 'cfdna'){ ## non-UMI samples

  print OUT $command1;
  
} else{
	
	## add in UMI info into bam file -- IDT UMI is in the index

	
	if($umi eq 'IDT'){
	
		print OUT $commandumi1;
	
	} else {
	
		print OUT $commandumi1a;
	
	}
	
	## generate UMI group
	print OUT $commandumi2;
	
	## divert pipeline for cfDNA or Q33
	
	if($seq_type eq 'cfdna'){
		print OUT $commandumicfdna;
	
	} else {
	
		print OUT $commandumi3;
	}
	

}

close OUT;

sub usage() {
        my $usage =<<END;

bwa.pl

Usage: perl bwa.pl [options]

Options

 --num_cpu			Number of cpus per BWA job, default: 12
 --java_tmp						Path to Java tmp folder
 --ref_fasta			Path to reference file
 --analysis_dir			Path to analysis directory
 --fastq_dir			Path to fastq directory
 --samp_id			Sample ID
 --seq_plat			Sequencing platform, default: Illumina
 --barcode_seq			Sequencing index
 --rm_dups				remove PCR duplicates: true or false
 --help				Display this message and quit
 --output_script		Path to write script to
 --seed					BWA seed length (default is 50bp )

Updated: 05th June 2016 by Wei.Yuan\@icr.ac.uk
Updated: 18th Dec 2018 by Lina.Yuan\@icr.ac.uk

END
        print $usage;
}
