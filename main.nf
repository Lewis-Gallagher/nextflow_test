#!/usr/bin/env nextflow

params.threads = 4

params.read_pairs = "$baseDir/data/*_R{1,2}_*.fastq.gz"
params.reference = "$baseDir/data/chr15.fa"
params.saveReference = true
params.outdir = "$baseDir/out"

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
    .into { reads }
		.println()

fasta_for_bwa_index = Channel
    .fromPath("${params.reference}")
fasta_for_samtools_index = Channel
    .fromPath("${params.reference}")

reference = file(params.reference)

process makeBWAindex {
		publishDir path: { params.saveReference ? "${params.outdir}/reference_genome" : params.outdir },
						 saveAs: { params.saveReference ? it : null }, mode: 'copy'
    input:
    file fasta from fasta_for_bwa_index

    output:
    file "*.{amb,ann,bwt,pac,sa}" into bwa_index

    script:
    """
    ~/Software/bwa/bwa index $fasta
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
    publishDir "${params.outdir}/align", mode: 'copy', pattern: "*.sam", overwrite: true

		input:
		set id, file(r1), file(r2) from reads
		file(index) from bwa_index
		file reference

		output:
		set id, file("${id}.aln.sam")

    script:
    """
        ~/Software/bwa/bwa mem -M -c 1 -k 50 -t ${params.threads} $params.reference ${r1} ${r2} -o ${id}.aln.sam

    """
}
