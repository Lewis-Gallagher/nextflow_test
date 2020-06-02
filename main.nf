#!/usr/bin/env nextflow
params.pool_id = "Pool_TEST"
params.samplesheet="$baseDir/Analysis/${params.pool_id}.hpc.csv"
params.threads = 4
params.fastqs = "$baseDir/fastq/${params.pool_id}/*_{R1,R2,R3}_*.fastq.gz"
params.reference = "/scratch/DMP/DUDMP/TRANSGEN/transgen-mdx/ngs/7.resources/hg19/hg19.fa"
params.outdir = "$baseDir/Analysis/"
params.read_structures = "+T +M +T"
params.seq_type = "cfdna"
params.umi="IDT"
params.umi_family_size = 1
params.expected_umi_len = 9
params.seed = 50
params.rm_dups=true

params.java_tmp="tmp/"
params.java_cpu=9
params.java_mem=128
params.fgbio_jar = "/apps/fgbio/0.6.1/fgbio-0.6.1.jar"
params.picard_jar="/apps/picard-tools/2.8.2/picard-2.8.2.jar"
params.picard="java -Xmx${params.java_mem}g -XX:ParallelGCThreads=${params.java_cpu} -jar ${params.picard_jar}"
params.fgbio="java -Xmx${params.java_mem}g -XX:+AggressiveOpts -XX:+AggressiveHeap -jar ${params.fgbio_jar}"
reference = file(params.reference)

//--------------------------------------------------


log.info """

	------------
	RMH PIPELINE
	------------
	FASTQ:	        ${params.fastqs}
	Reference: 	${params.reference}
	Output dir:     ${params.outdir}
	
 	
"""
.stripIndent()


logParams(params, params.outdir+ "/nextflow_parameters.txt")

// Writes all run parameters to text file
def logParams(p, n) {
	
	File file = new File(n)
	file.write "Parameter:\tValue\n"

	for(s in p) {
	file << "${s.key}:\t${s.value}\n"
	}
}



//^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
//
//                     CHANNELS
//
//vvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvv


//Channel.fromPath( file(params.samplesheet) )
//        .splitCsv(header: true, sep: ',')
//        .map{row ->
//	    def pool_id = row['pool_id']
//          def sample_ID = row['samp_id']
//          def run_id = row['run_id']
//	    def index_id = row['index_id']
//	    def seq_plat = row['seq_plat']
//	    def seq_type = row['seq_type']
//	    def target_bed = row['target_bed']
//	    def gatk_grp = row['gatk_grp']
//	    def vpanel = row['vpanel']
//	    def umi = row['umi']
//	    def tag = row['tag']
//	    def tumour_type = row['tumour_type']
//	    def trial_id = row['tumour_type']
//	    def tumour_content = row['tumour_content']
//
//	    return [ 'Pool ID': pool_id, 
//		     'Sample ID': sample_ID, 
//		     'Run ID': run_id, 
//		     'Index ID': index_id, 
//		     'Sequencing Platform': seq_plat, 
//		     'Analysis type': seq_type, 
//		     'Target BED': target_bed, 
//		     'GATK group': gatk_grp, 
//		     'Virtual Panel': vpanel, 
//		     'UMI setting': umi, 
//		     'Tumour/Normal': tag, 
//		     'Tumour type': tumour_type, 
//		     'Trial ID': trial_id, 
//		     'Tumour content': tumour_content  ]
//        }
//	.flatMap()
//        .into { samplesheet_info }
//
//samplesheet_info.subscribe { println "samplesheet_info: ${it}" }




// Send fastq pairs/trios to to reads channel
// Giving size: -1 should allow for any number of files in the group i.e. R1,R2 or R1,R2,R3 
Channel
	.fromFilePairs( params.fastqs, size: -1)
	.ifEmpty { exit 1, "Input FASTQs could not be found." }
	.set { fastq_channel }



//^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
//					  	 
//			MAIN			 
//						 
//vvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvv




//  /~/----------------------------------------\~\
//  \~\             Capture analysis           /~/
//  /~/----------------------------------------\~\


process bwa {

	when:
	params.seq_type == 'capture'

	input:
	set sample_id, file(infileFastqs) from fastq_channel

	output:
	set sample_id, file(outfileBam) into bwaOutput
	
	script:
	outfileBam = sample_id + 'aln.bam'
	"""
	bwa mem -M -c 1 -t ${params.threads} -k ${params.seed} -p \
        ${params.reference} ${infileFastqs} > ${outfileBam}

	"""
}


process sortSam {

	when:
	params.seq_type == 'capture'

	input:
	set sample_id, file(infileBam) from bwaOutput

	output:
	set sample_id, file(outfileBam) into sortSamOutput

	script:
	outfileBam = sample_id + '.sorted.bam'
	"""
	${params.picard} SortSam \
                SORT_ORDER=queryname \
                INPUT=${infileBam} \
                OUTPUT=${outfileBam} \
                CREATE_INDEX=true \
                VALIDATION_STRINGENCY=LENIENT \
                TMP_DIR=${params.java_tmp}

	"""
}


process markDuplicates {

	when:
        params.seq_type == 'capture'

	input:
        set sample_id, file(infileBam) from sortSamOutput

        output:
        set sample_id, file(outfileBam) into markDuplicatesOutput
	set sample_id, file(outfileMetrics) into markDuplicatesOutput_metrics

        script:
        outfileBam = sample_id + '.marked.bam'
	outfileMetrics = sample_id + '.mark_dups.metrics'
        """
	${picard} MarkDuplicates \
		REMOVE_DUPLICATES=true \
		INPUT=${infileBam} \
		OUTPUT=${outfileBam} \
		METRICS_FILE=${sample_id}.mark_dups.metrics \
		ASSUME_SORTED=true \
		CREATE_INDEX=true \
		TMP_DIR=${params.java_tmp} \
		VALIDATION_STRINGENCY=LENIENT
	"""
}


process clipBam {
        publishDir path: "${params.outdir}/Alignments", mode: 'copy', overwrite: true

        when:
        params.seq_type == 'capture'

        input:
        set sample_id, file(infileBam) from markDuplicatesOutput

        output:
        set sample_id, file(outfileBam) into clipBamOutput
        set sample_id, file('*.clipbam.metrics') into clipBamOutput_metrics
 
        script:
        outfileBam = sample_id + '.sort.bam'
	outfileMetrics = sample_id + '.clipbam.metrics'
        """
	${params.fgbio} ClipBam \
                --input ${infileBam} \
                --output ${outfileBam} \
                --metrics ${outfileMetrics} \
                --ref ${params.reference} \
                --clip-overlapping-reads=true \
                --clipping-mode Hard

	"""
}


//  /~/----------------------------------------\~\
//  \~\               Q33 analysis             /~/
//  /~/----------------------------------------\~\


process fastqToBam {

	when:
	params.seq_type == 'umi'	

	input:
	set sample_id, file(infileFastqs) from fastqs_channel

	output:
	set sample_id, file(outfileBam) into fastqToBamOutput

	script:
	outfileBam = sample_id + '.bam'
	"""
	${params.fgbio} FastqToBam \
                --input ${infileFastqs} \
                --read-structures ${params.read_structures} \
                --output ${sample_id}.precheck.sam \
                --sort true \
                --umi-tag RX \
                --sample ${sample_id} \
                --library ${sample_id} \
                --read-group-id ${sample_id}
	
	samtools view -h -o ${sample_id}.precheck.sam ${sample_id}.precheck.bam
	sed 's/RX:Z:\([A-Z]*\)-\([A-Z]*\)/RX:Z:\1\2/' ${sample_id}.precheck.sam > ${sample_id}.sam
	samtools view -S -b ${sample_id}.sam > ${outfileBam}
	"""
}


process markIlluminaAdapters {

        when:
        params.seq_type == 'umi'

        input:
        set sample_id, file(infileBam) from fastqToBamOutput

        output:
        set sample_id, file(outfileBam) into markIlluminaAdaptersOutput
	set sample_id, file(outfileMetrics) into markIlluminaAdaptersOutput_metrics

        script:
        outfileBam = sample_id + '.marked.bam'

	"""
	 ${params.picard} MarkIlluminaAdapters \
                INPUT=${infileBam} \
                OUTPUT=${outfileBam} \
                METRICS=${outfileMetrics} \
		TMP_DIR=${params.java_tmp}
	"""
}


process samToFastq {

        when:
        params.seq_type == 'umi'

        input:
        set sample_id, file(infileBam) from markIlluminaAdaptersOutput

        output:
        set sample_id, file(outfileFastq) into Output samToFastqOutput

        script:
        outfileBam = sample_id + '.fastq'
	"""
	${params.picard} SamToFastq \
                INPUT=${infileBam} \
                FASTQ=${outfileFastq} \
                CLIPPING_ATTRIBUTE=XT \
                CLIPPING_ACTION=2 \
                INTERLEAVE=true \
                NON_PF=true \
                TMP_DIR=${params.java_tmp}
	"""
}


process bwa {

        when:
        params.seq_type == 'umi'

        input:
        set sample_id, file(infileBam) from markIlluminaAdaptersOutput

        output:
        set sample_id, file(outfileBam) into bwaOutput

        script:
        outfileBam = sample_id + '.aln.bam'
	"""
        bwa mem -M -c 1 -t ${params.threads} -k ${params.seed} -p \
        ${params.reference} ${infileFastq} > ${outfileBam}
        """
}


process mergeBamAlignment

        when:
        params.seq_type == 'umi'

        input:
        set sample_id, file(infileBam) from bwaOutput

        output:
        set sample_id, file(outfileBam) into mergeBamAlignmentOutput

        script:
        outfileBam = sample_id + '.merged.bam'
	"""
        ${params.picard} MergeBamAlignment \
                R=${params.reference} \
                UNMAPPED_BAM=${infileUnmappedBam} \
                ALIGNED_BAM=${infileMappedBam} \
                OUTPUT=${outfileBam} \
                CREATE_INDEX=true \
                ADD_MATE_CIGAR=true \
                CLIP_ADAPTERS=false \
                CLIP_OVERLAPPING_READS=true \
                INCLUDE_SECONDARY_ALIGNMENTS=true \
                MAX_INSERTIONS_OR_DELETIONS=-1 \
                PRIMARY_ALIGNMENT_STRATEGY=MostDistant \
                ATTRIBUTES_TO_RETAIN=XS \
                TMP_DIR=${params.java_tmp}
        """
}


process groupReadsByUmi

        when:
        params.seq_type == 'umi'

        input:
        set sample_id, file(infileBam) from mergeBamAlignmentOutput

        output:
        set sample_id, outfileBam into groupReadsByUmiOutput
	set sample_id, outfileMetrics into groupReadsByUmiOutput_metrics

        script:
        outfileBam = sample_id + '.bam'
	outfileMetrics = sample_id + '.umi.metrics'
	"""
	${params.fgbio} GroupReadsByUmi \
                --input ${infileBam} \
                --family-size-histogram ${outfileMetrics} \
                --strategy adjacency \
                --min-map-q 30 \
                --raw-tag RX \
                --assign-tag MI \
                --min-umi-length ${params.expected_umi_len} \
                --output ${outfileBam}
	"""	
}


process sortSam {

        when:
        params.seq_type == 'umi'

        input:
        set sample_id, file(infileBam) from groupReadsByUmiOutput

        output:
        set sample_id, file(outfileBam) into sorSamOutput

        script:
        outfileBam = sample_id + '.bam'
	"""
	${params.picard} SortSam \
                SORT_ORDER=queryname \
                INPUT=${infileBam} \
                OUTPUT=${outfileBam} \
                CREATE_INDEX=true \
                VALIDATION_STRINGENCY=LENIENT \
                TMP_DIR=${params.java_tmp}
	"""
}


process markDuplicates {
	publishDir path: "${params.outdir}/Alignments", mode: 'copy', overwrite: true

        when:
        params.seq_type == 'umi'

        input:
        set sample_id, file(infileBam) from sorSamOutput

        output:
        set sample_id, file(outfileBam) into markDuplicatesOutput

        script:
        outfileBam = sample_id + '.bam'
	"""
	${picard} MarkDuplicates \
                REMOVE_DUPLICATES=true \
                INPUT=${infileBam} \
                OUTPUT=${outfileBam} \
                METRICS_FILE=${sample_id}.mark_dups.metrics \
                ASSUME_SORTED=true \
                CREATE_INDEX=true \
                TMP_DIR=${params.java_tmp} \
                VALIDATION_STRINGENCY=LENIENT	
	"""
}



//  /~/----------------------------------------\~\
//  \~\              cfDNA analysis            /~/
//  /~/----------------------------------------\~\



process fastqToBam {
	
	when:
	params.seq_type == 'cfdna'

	input:
	set sample_id, file(fastqs) from fastq_channel

	output:
	set sample_id, file(outfileBam) into fastqToBamOutput

	script:
	outfileBam = sample_id + '.unaln.bam'
	"""
	${params.fgbio}	FastqToBam \
		--input ${fastqs} \
		--read-structures ${params.read_structures} \
		--output ${outfileBam} \
		--sort true \
		--umi-tag RX \
		--sample ${sample_id} \
		--library ${sample_id} \
		--read-group-id ${sample_id} 
	"""
}


process  markIlluminaAdapters {

	publishDir path: "${params.outdir}/Alignments", pattern: "*.metrics",  mode: 'copy', overwrite: true
	
	when:
        params.seq_type == 'cfdna'

	input:
	set sample_id, file(infileBam) from fastqToBamOutput

	output:
	set sample_id, file(outfileBam) into markIlluminaAdaptersOutput1, markIlluminaAdaptersOutput2
	set sample_id, file('*.metrics') into markIlluminaAdaptersOutput_metrics

	script:
	outfileBam = sample_id + '.marked.bam'
	"""
	${params.picard} MarkIlluminaAdapters \
		INPUT=${infileBam} \
		OUTPUT=${outfileBam} \
		METRICS=${sample_id}.adapter.metrics
	"""
}


process samToFastq {
	
	when:
        params.seq_type == 'cfdna'	
	
	input:
	set sample_id, file(infileBam) from markIlluminaAdaptersOutput1

	output:
	set sample_id, file(outfileFastq) into samToFastqOutput

	script:
	outfileFastq = sample_id + '.fastq'
	"""
	${params.picard} SamToFastq \
		INPUT=${infileBam} \
		FASTQ=${outfileFastq} \
		CLIPPING_ATTRIBUTE=XT \
		CLIPPING_ACTION=2 \
		INTERLEAVE=true \
		NON_PF=true \
		TMP_DIR=${params.java_tmp}
	"""
}


process bwaRawReads {

	when:
        params.seq_type == 'cfdna'

	input:
        set sample_id, file(infileFastq) from samToFastqOutput

	output:
	set sample_id, file(outfileBam) into bwaRawReadsOuput

	script:
	outfileBam = sample_id + '.aln.bam' 
	"""
	bwa mem -M -c 1 -t ${params.threads} -k ${params.seed} -p \
	${params.reference} ${infileFastq} > ${outfileBam}
	"""
}


process mergeBamAlignment {

	when:
        params.seq_type == 'cfdna'

        input:
        set sample_id, file(infileMappedBam) from bwaRawReadsOuput
	set sample_id, file(infileUnmappedBam) from markIlluminaAdaptersOutput2

        output:
        set sample_id, file(outfileBam) into mergeBamAlignmentOutput

        script:
        outfileBam = sample_id + '.merged.bam'
	"""
	${params.picard} MergeBamAlignment \
		R=${params.reference} \
		UNMAPPED_BAM=${infileUnmappedBam} \
		ALIGNED_BAM=${infileMappedBam} \
		OUTPUT=${outfileBam} \
		CREATE_INDEX=true \
		ADD_MATE_CIGAR=true \
		CLIP_ADAPTERS=false \
		CLIP_OVERLAPPING_READS=true \
		INCLUDE_SECONDARY_ALIGNMENTS=true \
		MAX_INSERTIONS_OR_DELETIONS=-1 \
		PRIMARY_ALIGNMENT_STRATEGY=MostDistant \
		ATTRIBUTES_TO_RETAIN=XS \
		TMP_DIR=${params.java_tmp}
	"""

}


process groupReadsByUmi {

	when:
        params.seq_type == 'cfdna' || params.seq_type == 'AVENIO'
        
        input:
        set sample_id, file(infileBam) from mergeBamAlignmentOutput

        output:
        set sample_id, file(outfileBam) into groupReadsByUmiOutput
        set sample_id, file('*.metrics') into groupReadsByUmiOutput_metrics

        when:

        script:
        outfileBam = sample_id + '.grouped.bam'
        outfileMetrics = sample_id + '.umi.metrics'
        """
        ${params.fgbio} GroupReadsByUmi \
                --input ${infileBam} \
                --family-size-histogram ${outfileMetrics} \
                --strategy adjacency \
                --min-map-q 30 \
                --raw-tag RX \
                --assign-tag MI \
                --min-umi-length ${params.expected_umi_len} \
                --output ${outfileBam}
        """
}


process callMolecularConsensusReads {

	when:
	params.seq_type == 'cfdna' || params.seq_type == 'AVENIO'
     
	input:
        set sample_id, file(infileBam) from groupReadsByUmiOutput

        output:
        set sample_id, file(outfileBam) into callMolecularConsensusReadsOutput

        script:
        outfileBam = sample_id + '.consensus.bam'
        
	"""
        ${params.fgbio} CallMolecularConsensusReads \
            --input ${infileBam} \
            --output ${outfileBam} \
            --min-reads ${params.umi_family_size} \
            --min-input-base-quality 30 \
            --tag MI
        """
}


process filterConsensusReads {
	
	when:
	params.seq_type == 'cfdna' || params.seq_type == 'AVENIO'

        input:
        set sample_id, file(infileBam) from callMolecularConsensusReadsOutput

        output:
        set sample_id, file('*.filtered.fastq') into filterConsensusReadsFastqOutput
	set sample_id, file('*.filtered.bam') into filterConsensusReadsBamOutput

        script:
        outfileFastq = sample_id + '.filtered.fastq'

        """
        ${params.fgbio} FilterConsensusReads \
        	--input ${infileBam} \
        	--output ${sample_id}.filtered.bam \
        	--ref ${params.reference} \
        	--min-reads 3 \
        	--max-read-error-rate 0.05 \
       		--min-base-quality 40 \
        	--max-base-error-rate 0.1 \
        	--max-no-call-fraction 0.05 \
        	--reverse-per-base-tags true
	
	 ${params.picard} SamToFastq \
               	INPUT=${sample_id}.filtered.bam \
               	FASTQ=${outfileFastq} \
               	CLIPPING_ATTRIBUTE=XT \
               	CLIPPING_ACTION=2 \
               	INTERLEAVE=true \
               	NON_PF=true \
               	TMP_DIR=${params.java_tmp}
        """
}
 

process consensusBwa {

	when:
        params.seq_type == 'cfdna' || params.seq_type == 'AVENIO'

        input:
        set sample_id, file(infileFastq) from filterConsensusReadsFastqOutput

        output:
        set sample_id, file(outfileBam) into consensusBwaOutput

        script:
        outfileBam = sample_id + 'consensus.aln.bam'

        """
	bwa mem -M -c 1 -t ${params.threads} -k ${params.seed} -p \
        ${params.reference} ${infileFastq} \
		| ${params.picard} SortSam \
			SORT_ORDER=queryname \
			INPUT=/dev/stdin \
			OUTPUT=${outfileBam} \
			CREATE_INDEX=true \
			VALIDATION_STRINGENCY=LENIENT \
			TMP_DIR=${params.java_tmp}	
	"""
}


process consensusSortSam {

	when:
        params.seq_type == 'cfdna' || params.seq_type == 'AVENIO'

        input:
        set sample_id, file(infileBam) from filterConsensusReadsBamOutput

        output:
        set sample_id, file(outfileBam) into consensusSortSamOutput

        script:
        outfileBam = sample_id + 'consensus.unaln.bam'

        """
        ${params.picard} SortSam \
        	SORT_ORDER=queryname \
                INPUT=${infileBam} \
                OUTPUT=${outfileBam} \
                CREATE_INDEX=true \
                VALIDATION_STRINGENCY=LENIENT \
                TMP_DIR=${params.java_tmp}
        """
}


process consensusMergeBamAlignment {

        when:
        params.seq_type == 'cfdna'

        input:
	set sample_id, file(infileMappedBam) from consensusBwaOutput
	set sample_id, file(infileUnmappedBam) from consensusSortSamOutput
        
	output:
	set sample_id, file(outfileBam) into consensusMergeBamAlignmentOutput

        script:
        outfileBam = sample_id + '.merged.bam'
        """
        ${params.picard} MergeBamAlignment \
                R=${params.reference} \
                UNMAPPED_BAM=${infileUnmappedBam} \
                ALIGNED_BAM=${infileMappedBam} \
                OUTPUT=${outfileBam} \
                CREATE_INDEX=true \
                ADD_MATE_CIGAR=true \
                CLIP_ADAPTERS=false \
                CLIP_OVERLAPPING_READS=true \
                INCLUDE_SECONDARY_ALIGNMENTS=true \
                MAX_INSERTIONS_OR_DELETIONS=-1 \
                PRIMARY_ALIGNMENT_STRATEGY=MostDistant \
                ATTRIBUTES_TO_RETAIN=XS \
                TMP_DIR=${params.java_tmp}
        """
}


process consensusClipBam {

	publishDir path: "${params.outdir}/Alignments", mode: 'copy', overwrite: true

        when:
        params.seq_type == 'cfdna'

        input:
        set sample_id, file(infileBam) from consensusMergeBamAlignmentOutput

        output:
        set sample_id, file(outfileBam) into consensusClipBamOutput
	set sample_id, file(outfileMetrics) into consensusClipBamOutput_metrics	

        script:
        outfileBam = sample_id + '.sort.bam'
	outfileMetrics = sample_id + '.clipbam.metrics'
	"""
	${params.fgbio} ClipBam \
		--input ${infileBam} \
		--output ${outfileBam} \
		--metrics ${outfileMetrics} \
		--ref ${params.reference} \
		--clip-overlapping-reads=true \
		--clipping-mode Hard	
	"""
}
