import subprocess
import argparse
import os
import sys


#def get_args():
#	
#	parser = argparse.ArgumentParser(description='FASTQ to BAM alignment script')
#
#	parser.add_argument('--sample_id', type=str, required=True, help='Sample ID')
#	parser.add_argument('--fq1', type=os.path.abspath, required=True, help='Path to FASTQ read 1')
#	parser.add_argument('--fq2', type=os.path.abspath, required=True, help='Path to FASTQ read 2')
#	parser.add_argument('--ref', type=os.path.abspath, required=True, help='Path to reference fasta')
#	parser.add_argument('--seq_type', type=str, required=True, help='Sequencing type taken from samplesheet "seq_type" column')
#	parser.add_argument('--umi', type=str, required=True, help='UMI setting taken from samplesheet "umi" columns')
#
#	args = parser.parse_args()
#
#	return args
#


def kill(args):
	
	for f in [sample_id, fq1, fq2, ref]:
		if os.path.exists(f):
			if os.path.isfile(f):
				pass
		else:
			print(file, 'not found')
			sys.exit()


def align(args):


	seq_plat = 'ILLUMINA'
	barcode_seq = 'NNNNNNN'
	seed = 50
	num_cpu = 4
	java_cpu = 2
	java_mem = 32
	picard_jar = '/opt/gridware/apps/picard-tools/2.8.1/picard.jar'
	picard = f'java -Xmx{java_mem}g -XX:ParallelGCThreads={java_cpu} -jar {picard_jar}'	
	RG=f'@RG\tID:{sample_id}\tPL:{seq_plat}\tPU:{barcode_seq}\tLB:{sample_id}\tSM:{samp_id}'

	capture_command = f'''
	
	#!/bin/bash
d
	bwa mem -M -c 1 -t {num_cpu} -k {seed} -R {RG} \
	{ref} {fq1} {fq2} > {samp_id}.aln.sam
	
	
	{picard} SortSam MAX_RECORDS_IN_RAM=4000000 \
	SORT_ORDER=coordinate \
	INPUT={samp_id}.aln.sam \
	OUTPUT=$samp_id.fq.aln.sorted.bam \\
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

	'''


print(f'''

Hello {name}

Test

''')	


	
