#!/usr/bin/perl -w ######################################
#                                                       #
# gatk.pl                                               #
# Perl script to call the gatk                          #
#                                                       #
#                                                       #
# 13-11-13 Lina.Chen@icr.ac.uk - v0.1                   #
# 28-11-13 Lina.Chen@icr.ac.uk - v0.2                   #
# 13-05-14 Lina.Chen@icr.ac.uk - v0.3                   #
# 04-09-18 Lina.Chen@icr.ac.uk - updated to GATK4                   #
#########################################################
use strict;
use Getopt::Long;

# variables needed:
my $output_script = '';									# path to write gatk script
my $gatk='';										 	#path gatk;
my $java_tmp = '';                                      # temp folder for java
my $picard = '';                             			# path picard
my $num_cpu = 1;                                        # number of cpus
my $mem = '14g';                                        # default memory (GB)
my $ref_fasta = '';                                     # path to reference fasta
my $dbsnp_vcf = '';										# path to dbSNP VCF file
my $indel_vcf = '';										# path to 1k indel file
my $analysis_dir = '';                                  # path to analysis folder
my $proj_id = '';                                       # project ID
my $gatk_grp = '';                                      # gatk group id
my $tumour = '';                                       # tumour ID
my $control = '';                                       # control ID
my $chr_input='';												# get chr intervals for variant calling
my $bamtools = '';                                      # path bamtools
my $seqtype='';
my $control_setting='';
my $PON='';
my $id=1;
my $vardict_filter='/scratch/DMP/DUDMP/TRANSGEN/transgen-mdx/ngs/9.scripts/1.pipeline_modules/3.variantcalling/filter_vardict.py'; ## temp fix, will move to resource file after testing

my $args = scalar(@ARGV);
my $help;
my $result = GetOptions (
  "output_script=s" => \$output_script,
  "gatk=s" => \$gatk,
  "java_tmp=s" => \$java_tmp,
  "picard=s" => \$picard,
  "num_cpus=i" => \$num_cpu,
  "mem=s" => \$mem,
  "ref_fasta=s" => \$ref_fasta,
  "dbsnp_vcf=s" => \$dbsnp_vcf,
  "indel_vcf=s" => \$indel_vcf,
  "PON=s" => \$PON,
  "analysis_dir=s" => \$analysis_dir,
  "tumour=s" => \$tumour,
  "control=s" => \$control,
  "chr=s" => \$chr_input,
  "seqtype=s" => \$seqtype,
  "id=s" => \$id,
  "help" => \$help,
  );

# print usage message if requested or no args supplied
if(defined($help) || !$args) {
  &usage;
  exit(0);
}

my $gatk_input = $tumour.'.sort.bam';
my $vardict_cmd='';
my $gatk_output=$analysis_dir;
if($control ne 'None'){

	$control_setting = "--normal-sample $control";
	$gatk_input="$control.sort.bam -I ".$gatk_input;
	$vardict_cmd=<<END_OF_OUTPUT_SCRIPT;

mv $gatk_output/$id.${tumour}.recalibrated.bam $gatk_output/$id.${tumour}.merged.recalibrated.bam
java -Xmx$mem -Djava.io.tmpdir=$java_tmp -jar $gatk SplitReads \\
  -I $gatk_output/$id.${tumour}.merged.recalibrated.bam \\
  -R $ref_fasta  --split-sample \\
  -O $gatk_output
mv $gatk_output/$id.${tumour}.merged.recalibrated.${tumour}.bam  $gatk_output/$id.${tumour}.recalibrated.bam
mv $gatk_output/$id.${tumour}.merged.recalibrated.${tumour}.bai  $gatk_output/$id.${tumour}.recalibrated.bai


END_OF_OUTPUT_SCRIPT


}

#my @input_bams = map {"$analysis_dir/$proj_id/Alignments/".$_.'.sort.bam'} @samples;

#my $gatk_input = join(" -I ", @input_bams);
my @chrs=split(",", $chr_input);
my $chr=join(" -L ", @chrs);
my $bed = $chr;
$bed =~s/\.intervals//;

my $dupfilter='--read-filter NotDuplicateReadFilter';
$dupfilter='' if($seqtype eq 'amplicon' || $seqtype eq 'umi'  || $seqtype eq 'cfdna');

my $vardict_output="$id.SNV.${tumour}.vardict.vcf";
$vardict_output="$id.SNV.${tumour}.vcf" if($seqtype eq 'cfdna');

my $mutect_output="$id.SNV.${tumour}.mutect.vcf";
$mutect_output="$id.SNV.${tumour}.vcf" if($seqtype ne 'cfdna');


my $commands =<<END_OF_OUTPUT_SCRIPT;
#!/bin/bash


java -Xmx$mem -Djava.io.tmpdir=$java_tmp -jar $gatk  BaseRecalibrator \\
  -L $chr  \\
  --interval-set-rule UNION \\
  --interval-merging-rule OVERLAPPING_ONLY \\
  --interval-padding 200 \\
  -I $gatk_input \\
  -R $ref_fasta  \\
  -O $gatk_output/$id.${tumour}.recalibrated.table  \\
  --known-sites $dbsnp_vcf \\
  --known-sites $indel_vcf \\
  --read-filter NotSecondaryAlignmentReadFilter \\
  --read-filter NotSupplementaryAlignmentReadFilter $dupfilter \\
  --read-filter OverclippedReadFilter --filter-too-short 60 --dont-require-soft-clips-both-ends true 

java -Xmx$mem -Djava.io.tmpdir=$java_tmp -jar $gatk  ApplyBQSR \\
  -L $chr   \\
  --bqsr-recal-file $gatk_output/$id.${tumour}.recalibrated.table  \\
  --interval-set-rule UNION \\
  --interval-merging-rule OVERLAPPING_ONLY \\
  --interval-padding 200 \\
  -I $gatk_input \\
  -R $ref_fasta  \\
  -O $gatk_output/$id.${tumour}.recalibrated.bam  \\
  --read-filter NotSecondaryAlignmentReadFilter \\
  --read-filter NotSupplementaryAlignmentReadFilter $dupfilter \\
  --read-filter OverclippedReadFilter --filter-too-short 60 --dont-require-soft-clips-both-ends true
  
java -Xmx$mem -Djava.io.tmpdir=$java_tmp -jar $gatk Mutect2 -L $chr  \\
  --tumor-sample $tumour $control_setting \\
  -R $ref_fasta  \\
  -I $gatk_output/$id.${tumour}.recalibrated.bam \\
  -O $gatk_output/$mutect_output \\
  --dont-trim-active-regions true \\
  --dont-use-soft-clipped-bases true \\
  --interval-set-rule UNION \\
  --interval-merging-rule OVERLAPPING_ONLY \\
  --interval-padding 200 \\
  --max-reads-per-alignment-start 0 \\
  --annotation-group StandardMutectAnnotation \\
  --panel-of-normals $PON \\
  --germline-resource $indel_vcf \\
  -bamout $gatk_output/$id.SNV.${tumour}.bamout.bam \\
  --read-filter NotSecondaryAlignmentReadFilter \\
  --read-filter NotSupplementaryAlignmentReadFilter $dupfilter \\
  --read-filter OverclippedReadFilter --filter-too-short 60 --dont-require-soft-clips-both-ends false
  
  
## modify VarDict output to decrease The minimum mean position of variants in the read from 8 to 1
$vardict_cmd

VarDict -3 -G $ref_fasta -N $tumour -f 0.0001 -K -U -I 100 \\
    -k 0 -b  $gatk_output/$id.${tumour}.recalibrated.bam \\
    -z -c 1 -S 2 -E 3 \\
    $bed | teststrandbias.R | var2vcf_valid.pl -N $tumour -E -f 0.0001 -p 1 > $gatk_output/$id.${tumour}.vardict.unfiltered.vcf

## VarDict filter

python $vardict_filter -i $gatk_output/$id.${tumour}.vardict.unfiltered.vcf -o $gatk_output/$id.${tumour}.vardict.filtered.vcf

bgzip -c $gatk_output/$id.${tumour}.vardict.filtered.vcf > $gatk_output/$id.${tumour}.vardict.filtered.vcf.gz

tabix -p vcf  $gatk_output/$id.${tumour}.vardict.filtered.vcf.gz

bcftools view -f 'PASS'  $gatk_output/$id.${tumour}.vardict.filtered.vcf.gz  > $gatk_output/$vardict_output



END_OF_OUTPUT_SCRIPT
#  --max-assembly-region-size 500 \\
#   --af-of-alleles-not-in-resource -1 \\
#--read-filter SampleReadFilter \\
# --sample $tumour \\


open OUT, "> $output_script" or die "can't write script to $output_script: $!\n\n";
print OUT $commands;
close OUT;


sub usage() {
        my $usage =<<END;

snv.pl

Usage: perl snv.pl [options]

Options
 --gatk							Path to GATK
 --java							Path to Java
 --java_tmp						Path to Java tmp folder
 --picard						Path to picard
 --dbsnp_vcf						Path to dbSNP file
 --indel_vcf						Path to indel file from 1k genome project
 --num_cpu						Number of cpus per gatk job, default: 24
 --mem							Default memory per gatk job: 20g
 --ref_fasta						Path to reference file
 --analysis_dir						Path to analysis directory
 --proj_id						Project ID
 --samp_id						Sample ID
 --gatk_grp						GATK group
 --help							Display this message and quit
 --output_script					Path to write script to

created: 03-05-2019 Lina.Chen\@icr.ac.uk

END
        print $usage;
}
