#!/usr/bin/env nextflow

params.threads = 4

params.read_pairs = "$baseDir/data/1901548-BIOMEDE_B169-T_S32_{R1,R2}*.fastq.gz"
params.reference = "$baseDir/data/genome.fa"
params.saveReference = "true"
params.outdir = "$baseDir/out"

//markdups metrics
params.rm_dups=true
params.java_tmp=""
params.java_cpu=2

log.info """\

	TEST PIPELINE
	-------------
	Reads: ${params.read_pairs}
	Reference: ${params.reference}
  """
  .stripIndent()


Channel
    .fromFilePairs(params.read_pairs, flat: true)
    .ifEmpty { exit 1, "Read pairs could not be found: ${params.read_pairs}" }
    .set { reads }
		.println()

fasta_for_bwa_index = Channel
    .fromPath("${params.reference}")
fasta_for_samtools_index = Channel
    .fromPath("${params.reference}")

reference = file(params.reference)


// ---------------------------------------------



process makeBWAindex {
		publishDir path: { params.saveReference ? "${params.outdir}/reference_genome" : params.outdir },
						 saveAs: { params.saveReference ? it : null }, mode: 'copy'
    input:
    file fasta from fasta_for_bwa_index

    output:
    file "*.{amb,ann,bwt,pac,sa}" into bwa_index

    script:
    """
    bwa index $fasta
    """
}

process makeFASTAindex {
		publishDir path: { params.saveReference ? "${params.outdir}/reference_genome" : params.outdir },
						 saveAs: { params.saveReference ? it : null }, mode: 'copy'
    input:
		file fasta from fasta_for_samtools_index

    output:
    file "*.fai" into fai_index

    script:
    """
    samtools faidx $fasta
    """
}


process bwa {

    container 'lewis_bwa:v0.1'

    publishDir "${params.outdir}/align", mode: 'copy', pattern: "*.bam", overwrite: true

		input:
		set sample_id, file(r1), file(r2) from reads
		reference

		output:
		set sample_id, file("${sample_id}.aln.bam") into bwaOutput

    script:
    """
    bwa mem -M -c 1 -k 50 -t ${params.threads} ${reference} ${r1} ${r2} \
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
  publishDir "${params.outdir}/markdups", mode: 'copy', overwrite: true

  container 'broadinstitute/gatk:latest'

  input:
  set sample_id, file(infile_bam) from bwaOutput 

  output:
  set sample_id, file(outfile_bam) into markDuplicatesOutput_BAM
  set sample_id, file(outfile_metrics) into markDuplicatesOutput_QC

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
  publishDir "${params.outdir}/clipBam", mode: 'copy', overwrite: true

  container 'lewis_bwa:v0.1'

  input:
  set sample_id, file(infile_bam) from markDuplicatesOutput_BAM
  file(reference)

  output:
  set sample_id, file(outfile_bam) into clipBamOutput_BAM
  set sample_id, file(outfile_metrics) into clipBamOutput_QC

  script:
  outfile_bam = sample_id + ".sort.bam"
  outfile_metrics = sample_id + ".clipbam.metrics"

  """
  java -jar /usr/bin/fgbio.jar ClipBam \
    -i ${infile_bam}  \
    -o ${outfile_bam} \
    -m ${outfile_metrics} \
    -r ${reference} \
    --clip-overlapping-reads=true \
    -c Hard 
  """

}