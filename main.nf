#!/usr/bin/env nextflow

params.threads = 4
//params.fastqs = "$baseDir/data/fastq/*_{R1,R2,R3}_*.fastq.gz"
params.reference = "$baseDir/data/genome.fa"
params.saveReference = "true"
params.outdir = "$baseDir/out"
params.read_structures = "+T +T +M"
params.fgbio = "/apps/fgbio/0.6.1/fgbio-0.6.1.jar"
params.seq_type = "cfdna"
params.umi="IDT"
params.umi_family_size = 3
params.seed = 50
reference = file(params.reference)

//markdups metrics
params.rm_dups=true
params.java_tmp=""
params.java_cpu=2


//--------------------------------------------------


log.info """\
	------------
	RMH PIPELINE
	------------
	Read pairs:	${params.fastqs}
	Reference: 	${params.reference}
	Output dir: ${params.outdir}
	
 	
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
// 	Send fastq pairs/trios to to reads channel
//-------------------------------------------------

// Giving size: -1 should allow for any number of files in the group i.e. R1,R2 or R1,R2,R3 
Channel
    .fromFilePairs("$baseDir/data/fastq/*_{R1,R2,R3}_*.fastq.gz", size: -1)
    .ifEmpty { exit 1, "Input FASTQs could not be found." }
    .into { fastq_ch }


//// Index channels
//
//Channel
//    .fromPath(params.reference)
//    .ifEmpty { exit 1, "Reference fasta could not be found: ${params.reference}" }
//    .into { fastaForBWAIndex_ch; fastaForSamtoolsIndex_ch }
//
////fasta_for_bwa_index = Channel
// //   .fromPath("${params.reference}")
////fasta_for_samtools_index = Channel
////    .fromPath("${params.reference}")
//
//
//
//    
//  //--------------------------------------------\\
// //    		        WORKFLOW                   \\
////------------------------------------------------\\
//
//
//process makeBWAindex {
//    publishDir path: { params.saveReference ? "${params.outdir}/reference_genome" : params.outdir },
//       	 saveAs: { params.saveReference ? it : null }, mode: 'copy'
//
//    input:
//    file fasta from fastaForBWAIndex_ch
//
//    output:
//    file "*.{amb,ann,bwt,pac,sa}" into bwaIndex_ch
//
//    when:
//    params.seqtype != 'umi' || params.seqtype != 'cfdna'
//
//    script:
//
//    """
//    bwa index $fasta
//
//    """
//}
//
//process makeFASTAindex {
//    publishDir path: { params.saveReference ? "${params.outdir}/reference_genome" : params.outdir },
//       	saveAs: { params.saveReference ? it : null }, mode: 'copy'
//    
//    input:
//    file fasta from fastaForSamtoolsIndex_ch
//
//    output:
//    file "*.fai" into faiIndex_ch
//
//    when:
//    params.seq_type != 'umi' || params.seq_type != 'cfdna'
//
//    script:
//
//    """
//    samtools faidx $fasta
//
//    """
//}


process bwa {
    tag "$sample_id"
    container 'lewis_bwa:v0.1'

    publishDir "${params.outdir}/Alignments", mode: 'copy', pattern: "*.bam", overwrite: true

    input:
    set sample_id, file(reads) from bwaInput_ch

    output:
    set sample_id, file(outfileBam) into bwaOutput_ch

    when:
    params.seq_type != 'umi' || params.seq_type != 'cfdna'

    script:
    outfileBam = sample_id + '.aln.bam'

    """
    bwa mem -M -c 1 -k 50 -t ${params.threads} ${params.reference} ${reads} \
    | samtools sort -o $alnBam -
    
    """
}


//process sortSam {
//    tag "$sample_id"
//    publishDir "${params.outdir}/Alignments", mode: 'copy', overwrite: true
//    
//    input:
//    set sample_id, file(infileBam) bwaOutput_ch
//
//    output:
//    set sample_id, file(outfileBam) into sortSamOutput_ch
//
//    when:
//    params.seq_type != 'umi' || params.seq_type != 'cfdna'
//
//    script:
//
//    outfileBam = sample_id + '.fq.aln.sorted.bam'
//    
//    """
//
//    java -Xmx${params.java_mem}g -XX:ParallelGCThreads=$params.java_cpu -jar $params.picard SortSam MAX_RECORDS_IN_RAM=4000000 \
//        SORT_ORDER=coordinate \
//        INPUT=${infileBam} \
//        OUTPUT=${outfileBam} \
//        CREATE_INDEX=true \
//        TMP_DIR=${params.java_tmp} \
//        VALIDATION_STRINGENCY=LENIENT
//
//    """
//}


process markDuplicates {
    tag "$sample_id"
    publishDir "${params.outdir}/markdups", mode: 'copy', overwrite: true

    container 'broadinstitute/gatk:latest'

    input:
    set sample_id, file infile_bam from bwaOutput_ch

    output:
    set sample_id, file outfile_bam into markDuplicatesOutput_BAM_ch
    set sample_id, file outfile_metrics into markDuplicatesOutput_QC_ch

    when:
    params.seq_type != 'umi' || params.seq_type != 'cfdna'

    script:
    outfile_bam = sample_id + ".dedup.bam"
    outfile_metrics = sample_id + ".mark_dups.metrics"
  
    """
    java -Xmx${params.java_mem}g -XX:ParallelGCThreads=$params.java_cpu -jar picard.jar
        REMOVE_DUPLICATES=${params.rm_dups} \
        INPUT=${infile_bam} \
        OUTPUT=${outfile_bam} \
        METRICS_FILE=${outfile_metrics} \
        ASSUME_SORTED=true \
        CREATE_INDEX=false \
        TMP_DIR=${params.java_tmp} \
        VALIDATION_STRINGENCY=LENIENT
 
    """
}


process clipBam {
    tag "$sample_id"
    publishDir "${params.outdir}/Alignments", mode: 'copy', overwrite: true

    container 'lewis_bwa:v0.1'

    input:
    set sample_id, file(infile_bam) from markDuplicatesOutput_BAM_ch
    file(reference) 

    output:
    set sample_id, file outfile_bam into clipBamOutput_BAM_ch
    set sample_id, file outfile_metrics into clipBamOutput_QC_ch

    when:
    params.seq_type != 'umi' || params.seq_type != 'cfdna'

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


//  //----------------------------------------\\
//  \\              cfDNA analysis            //
//  //----------------------------------------\\



process FastqToBam {
    publishDir path: '${params.outdir}/Alignments', mode: 'copy', overwrite: true

    input:
    set sample_id, file(reads) from fastqs_ch

    output:
    file(outfileBam) into FastqToBam_ch

    when:
    params.seq_type == 'cfdna'


    script:
    // `--input $reads` will pass all reads matched in fastq glob i.e. R1,R2,R3 in IDT or just R1,R2 in Avenio.
    // read structres are defined with the global variable params.read_structures

    
    outfileBam = sample_id + 'unmapped.umi.bam'
    
    if ( params.umi == 'IDT')
        """
        java -Xmx${params.java_mem}g -XX:+AggressiveOpts -XX:+AggressiveHeap -jar $params.fgbio FastqToBam \
            --input $reads \
            --output ${sample_id}.fq.unaln.umi.bam \
            --read-structures $params.read_structures \
            --sort true \
            --umi-tag RX \
            --sample ${sample_id} \
            --library ${sample_id} \
            --read-group-id ${sample_id}

        """

    else if ( params.umi == 'AVENIO' )
        """
        java -jar params.fgbio FastqToBam \
            --input $reads \
            --output ${sample_id}.fq.unaln.umi.precheck.bam \
            --read-structures $params.read_structures \
            --sort true \
            --umi-tag RX \
            --sample ${sample_id} \
            --library ${sample_id} \
            --read-group-id ${sample_id}

        samtools view \
            -h -o ${sample_id}.fq.unaln.umi.precheck.sam \
            ${sample_id}.fq.unaln.umi.precheck.bam

        sed 's/RX:Z:\\([A-Z]*\\)-\\([A-Z]*\\)/RX:Z:\\1\\2/' \
        ${sample_id}.fq.unaln.umi.precheck.sam > \
        ${sample_id}.fq.unaln.umi.sam

        samtools view \
        -S -b ${sample_id}.fq.unaln.umi.sam > \
        ${outfileBam}

        """

}


process mergeBamAlignment {
    publishDir path: '${params.outdir}/Alignments', mode: 'copy', overwrite: true

    input:
    set sample_id, file(unalignedBam) from fastqToBam_ch
    file(alignedBam) from bwaOutput_ch
    file(reference)

    output:
    file(outfileBam) into mergedBam_ch

    when:
    params.seq_type == 'cfdna' || params.seq_type == 'AVENIO'

    script:
    outfileBam = sample_id + '.merged.bam'

    """
    java -Xmx${params.java_mem}g -XX:ParallelGCThreads=$java_cpu -jar $params.picard MergeBamAlignment MAX_RECORDS_IN_RAM=4000000 \
        R=${reference} \
        UNMAPPED_BAM=${unalignedBam} \
        ALIGNED_BAM=${alignedBam} \
        OUTPUT=${outfileBam} \
        CREATE_INDEX=true \
        ADD_MATE_CIGAR=true \
        CLIP_ADAPTERS=false \
        CLIP_OVERLAPPING_READS=true \
        INCLUDE_SECONDARY_ALIGNMENTS=true \
        MAX_INSERTIONS_OR_DELETIONS=-1 \
        PRIMARY_ALIGNMENT_STRATEGY=MostDistant \
        ATTRIBUTES_TO_RETAIN=XS \
        TMP_DIR=$params.java_tmp 

    """
}

process GroupReadsByUmi {

    publishDir path: '${params.outdir}/Alignments', mode: 'copy', overwrite: true

    input:
    set sample_id, file(mergedBam) from mergedBam_ch

    output:
    set sample_id, file(groupedBam) into groupedBam_ch
    set sample_id, file(groupedBam_metrics) into groupedBamMetrics_ch

    when:
    params.seq_type == 'cfdna' || params.seq_type == 'AVENIO'

    script:
    outfileBam = sample_id + '.grouped.bam'

    """
    java -Xmx${params.java_mem}g -XX:+AggressiveOpts -XX:+AggressiveHeap -jar $params.fgbio GroupReadsByUmi \
    --input ${sample_id}.fq.aln.merged.umi.bam \
    --family-size-histogram $groupedBam_metrics \
    --strategy adjacency \
    --min-map-q 30 \
    --raw-tag RX \
    --assign-tag MI \
    --min-umi-length $params.expected_umi_len \
    --output ${outfileBam}

    """
}


process callMolecularConsensusReads {
    publishDir path: '${params.outdir}/Alignments', mode: 'copy', overwrite: true

    input:
    set sample_id, file(infileBam) from groupedBam_ch

    output:
    file(outfileBam) into consensusBam_ch

    when:
    params.seq_type == 'cfdna' || params.seq_type == 'AVENIO'

    script:
    outfileBam = sample_id + '.consensus.bam'

    """
    java -Xmx${params.java_mem}g -XX:+AggressiveOpts -XX:+AggressiveHeap -jar $params.fgbio CallMolecularConsensusReads \
        --input ${infileBam} \
        --output ${outfileBam} \
        --min-reads ${params.umi_family_size} \
        --min-input-base-quality 30 \
        --tag MI
    """
}


process FilterConsensusReads {
    publishDir path: '${params.outdir}/Alignments', mode: 'copy', overwrite: true

    input:
    set sample_id, file(infileBam) from consensusBam_ch
    file(reference)

    output:
    file(outfileBam) into filterConsensusReads_ch

    when:
    params.seq_type == 'cfdna' || params.seq_type == 'AVENIO'

    script:
    outfileBam = sample_id + '.consensus.filtered.bam'
    
    """
    java -Xmx${params.java_mem}g -XX:+AggressiveOpts -XX:+AggressiveHeap -jar $params.fgbio FilterConsensusReads \
        --input ${infileBam} \
        --output ${outfileBam} \
        --ref ${reference} \
        --min-reads 3 \
        --max-read-error-rate 0.05 \
        --min-base-quality 40 \
        --max-base-error-rate 0.1 \
        --max-no-call-fraction 0.05 \
        --reverse-per-base-tags true
    """
}

process alignConsensusBam {
    publishDir path: '${params.outdir}/Alignments', mode: 'copy', overwrite: true

    input:
    set sample_id, file(infileBam) from filterConsensus_ch

    output:
    file(outfileBam) into alnConsensusBam_ch

    when:
    params.seq_type == 'cfdna' || params.seq_type == 'AVENIO'

    script:
    outfileBam = sample_id + '.consensus.aln.bam'

    """
    java -Xmx${params.java_mem}g -XX:ParallelGCThreads=$java_cpu -jar $params.picard FastqToSam MAX_RECORDS_IN_RAM=4000000 \
        INPUT=${infileBam} \
        FASTQ=/dev/stdout \
        CLIPPING_ATTRIBUTE=XT \
        CLIPPING_ACTION=2 \
        INTERLEAVE=true \
        NON_PF=true \
        TMP_DIR=${params.java_tmp} |
    bwa mem -M -c 1 -k $params.seed -t ${params.threads} -p /dev/stdin |
    java -Xmx${params.java_mem}g -XX:ParallelGCThreads=$java_cpu -jar $params.picard MergeBamAlignment MAX_RECORDS_IN_RAM=4000000 \
        SORT_ORDER=coordinate \
        INPUT=/dev/stdin \
        OUTPUT=${outfileBam} \
        CREATE_INDEX=true \
        TMP_DIR=${params.java_tmp} \
        ADD_MATE_CIGAR=true \
        CLIP_ADAPTERS=true \
        CLIP_OVERLAPPING_READS=true \
        INCLUDE_SECONDARY_ALIGNMENTS=true \
        VALIDATION_STRINGENCY=LENIENT

    """

}

process clipConsensusBam {
    publishDir path: '${params.outdir}/Alignments', mode: 'copy', overwrite: true

    input:
    set sample_id, file(infileBam) from alnConsensusBam_ch
    file(reference)

    output:
    file(outfileBam) into sortBam_ch
    file(outfileMetrics) into clipbamMetrics_ch

    when:
    params.seq_type == 'cfdna' || params.seq_type == 'AVENIO'

    script:
    outfileBam = sample_id + '.sort.bam'
    outfileMetrics = sample_id + '.clipbam.metrics'
    
    """
    java -Xmx${params.java_mem}g -XX:+AggressiveOpts -XX:+AggressiveHeap -jar $params.fgbio ClipBam \
        --input ${infileBam} \
        --ouitput ${outfileBam} \
        --metrics ${outfileMetrics} \
        --ref $reference \
        --clip-overlapping-reads=true \
        --clipping-mode Hard

    """

}



//process sortSam_cfDNA {
//    publishDir path: '${params.outdir}/Alignments', mode: 'copy', overwrite: true
//
//    input:
//    set sample_id, file(infileBam) from groupedBam_ch
//
//    output:
//    file(outfileBam) into sortSam_ch2
//
//    when:
//    params.seq_type = 'cfdna'
//
//    script:
//    outfileBam = sample_id + 'fq.aln.merged.umi.grp.sort.bam'
//
//    """
//    java -Xmx${params.java_mem}g -XX:ParallelGCThreads=$java_cpu -jar $params.picard SortSam MAX_RECORDS_IN_RAM=4000000 \
//        SORT_ORDER=coordinate \
//        INPUT=${infileBam} \
//        OUTPUT=${outfileBam} \
//        CREATE_INDEX=true \
//        TMP_DIR=$java_tmp \
//        VALIDATION_STRINGENCY=LENIENT
//
//    """
//}
//
//
//process MarkDuplicates {
//    publishDir path: '${params.outdir}/Alignments', mode: 'copy', overwrite: true
//
//    input:
//    set sample_id, file(infileBam) from sortSam_ch2
//
//    output:
//    file(outfileBam) into markDups_ch
//
//    when:
//    params.seq_type = 'cfdna'
//
//    script:
//    outfileBam = sample_id + '.sort.bam'
//
//}




//Channel.fromPath('test_*.txt').set{ test_ch }
//
//process foo {
//	input:
//	file x from test_ch
//	
//	script:
//	extension = x.getExtension()
//	"""
//	 echo ${extension} > extension.txt
//	"""
//}
//
//
//targets = [ file(params.intervals), file(params.bed_interval), file(params.hotspots) ]
//targets_extensions = ['', '.targetbed', '.hotspots' ]
//
//process calculateHsMetrics {
//
//	publishDir "${params.outdir}/qc", mode: 'copy', overwrite: true
//
//	input:
//	file bwaBam from clipBamOutput_BAM_ch
//	each target from targets
//	each extension from targets_extensions
//	
//
//	output:
//	set sample_id, file "*.{qc.metrics*, *qc.coverage*, *base.coverage*}" into hsMetricsOutput_ch
//
//	script:
//	
//	"""
//    gatk CalculateHsMetrics \
//    	REFERENCE_SEQUENCE=${reference} \
//    	INPUT=${bwaBam} \
//    	OUTPUT=${sample_id}.qc.metrics${ext} \
//    	BAIT_INTERVALS=${targets} \
//    	TARGET_INTERVALS=${targets} \
//    	BAIT_SET_NAME=${sample_id} \
//    	METRIC_ACCUMULATION_LEVEL=ALL_READS \
//    	TMP_DIR=${params.java_tmp} \
//    	PER_TARGET_COVERAGE=${sample_id}.qc.coverage${ext} \
//    	PER_BASE_COVERAGE=${sample_id}.base.coverage${ext} \
//    	VALIDATION_STRINGENCY=LENIENT    
//
//	"""
//
//process collectTargedPcrMetrics {
//
//    input:
//    file bwaBam from clipBamOutput_BAM_ch
//    file targets from target_ch
//    file reference
//
//    output:
//    set sample_id, file outfile_hsmetrics into hsMetricsOutput_ch
//
//    script:
//
//    """
//    CollectTargetedPcrMetrics \
//	REFERENCE_SEQUENCE=${reference} \
//	INPUT=${bwaBam} \
//	OUTPUT=${sample_id}.qc.metrics \
//	AMPLICON_INTERVALS=${targets} \
//	TARGET_INTERVALS=${targets} \
//	CUSTOM_AMPLICON_SET_NAME=${sample_id} \
//	METRIC_ACCUMULATION_LEVEL=ALL_READS \
//	TMP_DIR=${params.java_tmp{ \
//	PER_TARGET_COVERAGE=${sample_id}.qc.coverage \
//	PER_BASE_COVERAGE=${sample_id}.base.coverage \
//	VALIDATION_STRINGENCY=LENIENT
//
//    perl $getexoncoverage --dir $output_dir --samp $samp_id --tcontent $tcontent
//    cp $output_dir/${samp_id}* $transfer_dir/.
//
//
//    """
//}
//
//
//
//process baseRecalibartor {
//    publishDir "${params.outdir}/snv", mode: 'copy', overwrite: true
//
//    input:
//    file(bam) from clipBamOutput_ch_BAM
//    file(reference)
//    
//    output:
//    set sample_id, file(recal_table) into baseRecalibratorOutput_ch
//
//    script:
//    
//    """
//    gatk BaseRecalibrator \
//    	-L ??? \
//  		--interval-set-rule UNION \
//		--interval-merging-rule OVERLAPPING_ONLY \
// 		--interval-padding 200 \
//  		-I $gatk_input \
//  		-R $ref_fasta  \
//  		-O ${sample_id}.recalibrated.table \
//  		--known-sites $dbsnp_vcf \
//  		--known-sites $indel_vcf \
//  		--read-filter NotSecondaryAlignmentReadFilter \
//  		--read-filter NotSupplementaryAlignmentReadFilter $dupfilter \
//  		--read-filter OverclippedReadFilter \
//		--filter-too-short 60 \
//		--dont-require-soft-clips-both-ends true
//
//    """
//    }
//
//
//process applyBQSR {
//    publishDir "${params.outdir}/snv", mode: 'copy', overwrite: true
//
//	input:
//	file(recal_table) from baseRecalibratorOutput_ch
//	file(bam) from clipBamOutput_ch_BAM
//	file(reference)
//
//	output:
//	set sample_id, file(recalibratedBam) into applyBQSROutput_ch
//
//	script:
//
//	"""
//	gatk ApplyBQSR \
//		-L ??? \
//		--bqsr-recal-file ${recal_table} \
//		--interval-set-rule UNION \
//		--interval-merging-rule OVERLAPPING_ONLY \
//		--interval-padding 200 \
//		-I ${bam} \
//		-R ${reference} \
//		-O ${recalibratedBam} \
//		--read-filter NotSecondaryAlignmentReadFilter \
// 		--read-filter NotSupplementaryAlignmentReadFilter $dupfilter \
//  		--read-filter OverclippedReadFilter \
//		--filter-too-short 60 \
//		--dont-require-soft-clips-both-ends true
//	
//	"""
//
//
//    }





