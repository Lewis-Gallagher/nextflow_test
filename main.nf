#!/usr/bin/env nextflow

params.threads = 4
params.read_pairs = "$baseDir/data/fastq/*_{R1,R2}*.fastq.gz"
params.reference = "$baseDir/data/genome.fa"
params.saveReference = "true"
params.outdir = "$baseDir/out"

//markdups metrics
params.rm_dups=true
params.java_tmp=""
params.java_cpu=2


//--------------------------------------------------


log.info """\
	------------
	RMH PIPELINE
	------------
	Read pairs:	${params.read_pairs}
	Reference: 	${params.reference}
	Output dir: 	${params.outdir}
	
 	
"""
.stripIndent()

logParams(params, "nextflow_parameters.txt")


def logParams(p, n) {
   File file = new File(n)
   file.write "Parameter:\tValue\n"

   for(s in p) {
     file << "${s.key}:\t${s.value}\n"
 }
}


//-------------------------------------------------
// 	Send fastq pairs to to reads channel
//-------------------------------------------------

Channel
    .fromFilePairs(params.read_pairs, flat: false, size:2)
    .ifEmpty { exit 1, "Read pairs could not be found: ${params.read_pairs}" }
    .into { reads_ch }


// Index channels

Channel
    .fromPath(params.reference)
    .ifEmpty { exit 1, "Reference fasta could not be found: ${params.reference}" }
    .into { fastaForBWAIndex_ch; fastaForSamtoolsIndex_ch }

//fasta_for_bwa_index = Channel
 //   .fromPath("${params.reference}")
//fasta_for_samtools_index = Channel
//    .fromPath("${params.reference}")

reference = file(params.reference)


    
  //--------------------------------------------\\
 //    		     WORKFLOW                    \\
//------------------------------------------------\\


process makeBWAindex {
    publishDir path: { params.saveReference ? "${params.outdir}/reference_genome" : params.outdir },
       	 saveAs: { params.saveReference ? it : null }, mode: 'copy'

    input:
    file fasta from fastaForBWAIndex_ch

    output:
    file "*.{amb,ann,bwt,pac,sa}" into bwaIndex_ch

    script:

    """
    bwa index $fasta

    """
}

process makeFASTAindex {
    publishDir path: { params.saveReference ? "${params.outdir}/reference_genome" : params.outdir },
       	saveAs: { params.saveReference ? it : null }, mode: 'copy'
    
    input:
    file fasta from fastaForSamtoolsIndex_ch

    output:
    file "*.fai" into faiIndex_ch

    script:

    """
    samtools faidx $fasta

    """
}


process bwa {
    tag "$sample_id"
    container 'lewis_bwa:v0.1'

    publishDir "${params.outdir}/align", mode: 'copy', pattern: "*.bam", overwrite: true

    input:
    set sample_id, file(reads) from reads_ch

    output:
    set sample_id, file "${sample_id}.aln.bam" into bwaOutput_ch

    script:
    
    """
    bwa mem -M -c 1 -k 50 -t ${params.threads} ${params.reference} ${reads} \
    | samtools sort -o ${sample_id}.aln.bam -
    
    """
}

// process markduplicates {
//   container 'broadinstitute/gatk:latest'

//   output:
//   stdout result

//   """
//   gatk MarkDuplicates -h
//   """
// }

// result.view { it.trim() }


process markDuplicates {
    tag "$sample_id"
    publishDir "${params.outdir}/markdups", mode: 'copy', overwrite: true

    container 'broadinstitute/gatk:latest'

    input:
    set sample_id, file infile_bam from bwaOutput_ch

    output:
    set sample_id, file outfile_bam into markDuplicatesOutput_BAM_ch
    set sample_id, file outfile_metrics into markDuplicatesOutput_QC_ch

    script:
    outfile_bam = sample_id + ".dedup.bam"
    outfile_metrics = sample_id + ".mark_dups.metrics"
  
    """
    gatk MarkDuplicates \
        --REMOVE_DUPLICATES=${params.rm_dups} \
        --INPUT=${infile_bam} \
        --OUTPUT=${outfile_bam} \
        --METRICS_FILE=${outfile_metrics} \
        --ASSUME_SORTED=true \
        --CREATE_INDEX=false \
        --TMP_DIR=${params.java_tmp} \
        --VALIDATION_STRINGENCY=LENIENT
 
    """
}


process clipBam {
    tag "$sample_id"
    publishDir "${params.outdir}/clipBam", mode: 'copy', overwrite: true

    container 'lewis_bwa:v0.1'

    input:
    set sample_id, file infile_bam from markDuplicatesOutput_BAM_ch
    file reference 

    output:
    set sample_id, file outfile_bam into clipBamOutput_BAM_ch
    set sample_id, file outfile_metrics into clipBamOutput_QC_ch

    script:
    outfile_bam = sample_id + ".sort.bam"
    outfile_metrics = sample_id + ".clipbam.metrics"

    """
    java -jar /usr/bin/fgbio.jar ClipBam \
        --input ${infile_bam}  \
        --output ${outfile_bam} \
        --metrics ${outfile_metrics} \
        --ref ${reference} \
        --clip-overlapping-reads=true \
        --clipping-mode Hard 
    
    """
}




Channel.fromPath('test_*.txt').set{ test_ch }

process foo {
	input:
	file x from test_ch
	
	script:
	extension = x.getExtension()
	"""
	 echo ${extension} > extension.txt
	"""
}


targets = [ file(params.intervals), file(params.bed_interval), file(params.hotspots) ]
targets_extensions = ['', '.targetbed', '.hotspots' ]

process calculateHsMetrics {

	publishDir "${params.outdir}/qc", mode: 'copy', overwrite: true

	input:
	file bwaBam from clipBamOutput_BAM_ch
	each target from targets
	each extension from targets_extensions
	

	output:
	set sample_id, file "*.{qc.metrics*, *qc.coverage*, *base.coverage*}" into hsMetricsOutput_ch

	script:
	
	"""
    gatk CalculateHsMetrics \
    	REFERENCE_SEQUENCE=${reference} \
    	INPUT=${bwaBam} \
    	OUTPUT=${sample_id}.qc.metrics${ext} \
    	BAIT_INTERVALS=${targets} \
    	TARGET_INTERVALS=${targets} \
    	BAIT_SET_NAME=${sample_id} \
    	METRIC_ACCUMULATION_LEVEL=ALL_READS \
    	TMP_DIR=${params.java_tmp} \
    	PER_TARGET_COVERAGE=${sample_id}.qc.coverage${ext} \
    	PER_BASE_COVERAGE=${sample_id}.base.coverage${ext} \
    	VALIDATION_STRINGENCY=LENIENT    

	"""

process collectTargedPcrMetrics {

    input:
    file bwaBam from clipBamOutput_BAM_ch
    file targets from target_ch
    file reference

    output:
    set sample_id, file outfile_hsmetrics into hsMetricsOutput_ch

    script:

    """
    CollectTargetedPcrMetrics \
	REFERENCE_SEQUENCE=${reference} \
	INPUT=${bwaBam} \
	OUTPUT=${sample_id}.qc.metrics \
	AMPLICON_INTERVALS=${targets} \
	TARGET_INTERVALS=${targets} \
	CUSTOM_AMPLICON_SET_NAME=${sample_id} \
	METRIC_ACCUMULATION_LEVEL=ALL_READS \
	TMP_DIR=${params.java_tmp{ \
	PER_TARGET_COVERAGE=${sample_id}.qc.coverage \
	PER_BASE_COVERAGE=${sample_id}.base.coverage \
	VALIDATION_STRINGENCY=LENIENT

    perl $getexoncoverage --dir $output_dir --samp $samp_id --tcontent $tcontent
    cp $output_dir/${samp_id}* $transfer_dir/.


    """
}



process baseRecalibartor {
    publishDir "${params.outdir}/snv", mode: 'copy', overwrite: true

    input:
    file(bam) from clipBamOutput_ch_BAM
    file(reference)
    
    output:
    set sample_id, file(recal_table) into baseRecalibratorOutput_ch

    script:
    
    """
    gatk BaseRecalibrator \
    	-L ??? \
  		--interval-set-rule UNION \
		--interval-merging-rule OVERLAPPING_ONLY \
 		--interval-padding 200 \
  		-I $gatk_input \
  		-R $ref_fasta  \
  		-O ${sample_id}.recalibrated.table \
  		--known-sites $dbsnp_vcf \
  		--known-sites $indel_vcf \
  		--read-filter NotSecondaryAlignmentReadFilter \
  		--read-filter NotSupplementaryAlignmentReadFilter $dupfilter \
  		--read-filter OverclippedReadFilter \
		--filter-too-short 60 \
		--dont-require-soft-clips-both-ends true

    """
    }


process applyBQSR {
    publishDir "${params.outdir}/snv", mode: 'copy', overwrite: true

	input:
	file(recal_table) from baseRecalibratorOutput_ch
	file(bam) from clipBamOutput_ch_BAM
	file(reference)

	output:
	set sample_id, file(recalibratedBam) into applyBQSROutput_ch

	script:

	"""
	gatk ApplyBQSR \
		-L ??? \
		--bqsr-recal-file ${recal_table} \
		--interval-set-rule UNION \
		--interval-merging-rule OVERLAPPING_ONLY \
		--interval-padding 200 \
		-I ${bam} \
		-R ${reference} \
		-O ${recalibratedBam} \
		--read-filter NotSecondaryAlignmentReadFilter \
 		--read-filter NotSupplementaryAlignmentReadFilter $dupfilter \
  		--read-filter OverclippedReadFilter \
		--filter-too-short 60 \
		--dont-require-soft-clips-both-ends true
	
	"""


    }





