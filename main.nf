#!/usr/bin/env nextflow
params.scriptDir = "${baseDir}/bin"
params.pool_id = "Pool_1261"
params.samplesheet="/$baseDir/data/Analysis/${params.pool_id}.hpc.csv"
params.bwa_cpus = 4
params.cfdna_cpus = 10
params.fastqs = "/$baseDir/data/fastq/${params.pool_id}/*_{R1,R2,R3}_*.fastq.gz"
params.reference = "/$baseDir/data/hg19/hg19.fa"
params.outdir = "/$baseDir/data/Analysis/${pool_id}"
params.resources = "/scratch/DMP/DUDMP/TRANSGEN/transgen-mdx/ngs/7.resources"
params.seed = 50
params.min_reads=3
params.java_tmp="${baseDir}/tmp/"
params.java_cpu=9
params.java_mem=128
params.fgbio_jar = "/opt/gridware/apps/fgbio/0.6.1/fgbio-0.6.1.jar"
params.picard_jar = "/opt/gridware/apps/picard-tools/2.8.2/picard-2.8.2.jar"
//params.fgbio_jar = "/usr/local/bin/fgbio.jar"
//params.picard_jar = "/usr/local/bin/picard.jar"
params.picard="java -Xmx${params.java_mem}g -XX:ParallelGCThreads=${params.java_cpu} -jar ${params.picard_jar}"
params.fgbio="java -Xmx${params.java_mem}g -XX:ParallelGCThreads=${params.java_cpu} -XX:+AggressiveOpts -XX:+AggressiveHeap -jar ${params.fgbio_jar}"
reference = file(params.reference)

params.mutect_cpus = 4
params.gatk_mem = "30g"
params.gatk_engine = "/apps/gatk/4.0.5.1/gatk-package-4.0.5.1-local.jar"
params.dbsnp = "${params.resources}/hg19/bedfiles/dbsnp_142.b37.20150102.snp.withchr.vcf"
params.indel = "${params.resources}/hg19/bedfiles/dbsnp_142.b37.20150102.indel.withchr.vcf"
params.PON_WGS = "${params.resources}/hg19/PanelOfNormal/panel_of_normal_WGS.vcf.gz"

//--------------------------------------------------


log.info """

        ------------
        RMH PIPELINE
        ------------
        FASTQ:          ${params.fastqs}
        Reference:      ${params.reference}
        Output dir:     ${params.outdir}


"""
.stripIndent()



// \=\--------------------------------------\=\
// /=/		     CHANNELS	  	    /=/	
// \=\--------------------------------------\=\



Channel.fromPath( file(params.samplesheet) )
        .splitCsv(header: true, sep: ',')
        .map{row ->[
           pool_id: row.pool_id,
           sample_id: row.samp_id,
           run_id: row.run_id,
           index_id: row.index_id,
           seq_plat: row.seq_plat,
           seq_type: row.seq_type,
           target_bed: row.target_bed,
           gatk_grp: row.gatk_grp,
           vpanel: row.vpanel,
           umi: row.umi,
           tag: row.tag,
           tumour_type: row.tumour_type,
           trial_id: row.tumour_type,
           tumour_content: row.tumour_content,
	   fastq_path:file("/${baseDir}/data/fastq/${row.pool_id}/${row.samp_id}*_R{1,2,3}_*.fastq.gz")]
	}
	.set { samplesheet }


process get_bed {

	input:
	val(info) from samplesheet

	output:
	val(info) into samplesheet_align
//	tuple val(sample_id), val(bedfile) into bedfile_path_qc, bedfile_path_snv 
//	tuple val(sample_id), val(info) into bedfile_path_qc, bedfile_path_snv

	script:
	sample_id = info['sample_id']
	panel = info['target_bed']

	if ( panel == "ABCBIO" )
		info["bedfile"] = params.ABCBIO
	if ( panel == "ABCBIOv1_1" )
		info["bedfile"] = params.ABCBIOv1_1
	if ( panel == "AFP5_amplicon" )
		info["bedfile"] = params.AFP5_amplicon
	if ( panel == "AgilentctDNA" )
		info["bedfile"] = params.AgilentctDNA
	if ( panel == "AVENIOE" )
		info["bedfile"] = params.AVENIOE
	if ( panel == "AVENIOT" )
		info["bedfile"] = params.AVENIOT
	if ( panel == "BRCA" )
		info["bedfile"] = params.BRCA
	if ( panel == "CRUK" )
		info["bedfile"] = params.CRUK
	if ( panel == "ct_BREAST" )
		info["bedfile"] = params.ct_BREAST
	if ( panel == "ct_GI" )
		info["bedfile"] = params.ct_GI
	if ( panel == "ct_PAED_Diagnostic" )
		info["bedfile"] = params.ct_PAED_Diagnostic
	if ( panel == "ct_PAED_Diagnostic2" )
		info["bedfile"] = params.ct_PAED_Diagnostic2
	if ( panel == "DALEK_TFord" )
		info["bedfile"] = params.DALEK_TFord
	if ( panel == "DDR" )
		info["bedfile"] = params.DDR
	if ( panel == "DDRv1" )
		info["bedfile"] = params.DDRv1
	if ( panel == "exome" )
		info["bedfile"] = params.exome
	if ( panel == "GI" )
		info["bedfile"] = params.GI
	if ( panel == "GI2" )
		info["bedfile"] = params.GI2
	if ( panel == "GSIC" )
		info["bedfile"] = params.GSIC
	if ( panel == "hg19_cds_bed" )
		info["bedfile"] = params.hg19_cds_bed
	if ( panel == "IDTExome" )
		info["bedfile"] = params.IDTExome
	if ( panel == "MEDEXOME" )
		info["bedfile"] = params.MEDEXOME
	if ( panel == "MolecularID" )
		info["bedfile"] = params.MolecularID
	if ( panel == "OGTHaem" )
		info["bedfile"] = params.OGTHaem
	if ( panel == "PAED" )
		info["bedfile"] = params.PAED
	if ( panel == "PAED2" )
		info["bedfile"] = params.PAED2
	if ( panel == "PAEDFUSION" ) 
		info["bedfile"] = params.PAEDFUSION
	if ( panel == "PrimeExome" )
		info["bedfile"] = params.PrimeExome
	if ( panel == "refflatgene" )
		info["bedfile"] = params.refflatgene
	if ( panel == "RMHQ33" )
		info["bedfile"] = params.RMHQ33
	if ( panel == "RMH200" )
		info["bedfile"] = params.RMH200
	if ( panel == "RMSfusion" ) 
		info["bedfile"] = params.RMSfusion
	if ( panel == "TST170" )
		info["bedfile"] = params.TST170
	if ( panel == "TwistExome" )
		info["bedfile"] = params.TwistExome
	if ( panel == "wgs" )
		info["bedfile"] = params.wgs
 


	"""
	"""

}


params.ABCBIO = "${params.resources}/hg19/bedfiles/ABCBIO_v1.0"
params.ABCBIOv1_1 = "${params.resources}/hg19/bedfiles/ABCBIO_v1.3"
params.AFP5_amplicon = "${params.resources}/hg19/bedfiles/AFP5_amplicon_track"
params.AgilentctDNA = "${params.resources}/hg19/bedfiles/panel_AgilentGI.v2"
params.AVENIOE = "${params.resources}/hg19/bedfiles/Avenio_Expanded_Panel"
params.AVENIOT = "${params.resources}/hg19/bedfiles/Avenio_targeted_panel"
params.BRCA = "${params.resources}/hg19/bedfiles/BRCA_amplicon"
params.CRUK = "${params.resources}/hg19/bedfiles/CRUK_SMP2_SNV_5bp_buffer_SMP2-02"
params.ct_BREAST = "${params.resources}/hg19/bedfiles/ct_BREAST"
params.ct_GI = "${params.resources}/hg19/bedfiles/ct_GI"
params.ct_PAED_Diagnostic = "${params.resources}/hg19/bedfiles/ct_PAED_Diagnostic"
params.ct_PAED_Diagnostic2 = "${params.resources}/hg19/bedfiles/ct_PAED_Diagnostic_v2"
params.DALEK_TFord = "${params.resources}/hg19/bedfiles/DALEK_TFord"
params.DDR = "${params.resources}/hg19/bedfiles/DDR"
params.DDRv1 = "${params.resources}/hg19/bedfiles/panel_DDR.v1.0"
params.exome = "/scratch/DCS/CANBIO/wyuan/ngs/resources/GRCh37/rod_files/SureSelect_V4/SureSelect_V4_All_Exon.anno.bed"
params.GI = "${params.resources}/hg19/bedfiles/format_w5bp_V1.0"
params.GI2 = "${params.resources}/hg19/bedfiles/format_v2.2"
params.GSIC = "${params.resources}/hg19/bedfiles/GSIC"
params.hg19_cds_bed = "${params.resources}/hg19/bedfiles/hg19_CDS"
params.IDTExome = "${params.resources}/hg19/bedfiles/IDTExome"
params.MEDEXOME = "${params.resources}/hg19/bedfiles/MEDEXOME"
params.MolecularID = "${params.resources}/hg19/bedfiles/MolecularID"
params.OGTHaem = "${params.resources}/hg19/bedfiles/panel_OGTHaem.v1.0"
params.PAED = "${params.resources}/hg19/bedfiles/panel_PAED_V5_w5.1"
params.PAED2 = "${params.resources}/hg19/bedfiles/panel_PAED_V2"
params.PAEDFUSION = "${params.resources}/hg19/bedfiles/161108_HG19_FusionV1_EZ_HX3_capture_targets"
params.PrimeExome = "${params.resources}/hg19/bedfiles/PrimeExome"
params.refflatgene = "/scratch/DCS/CANBIO/wyuan/ngs/resources/GRCh37/genomes/refFlat.gene.nochr.bed"
params.RMHQ33 = "${params.resources}/hg19/bedfiles/RMHQ33"
params.RMH200 = "${params.resources}/hg19/bedfiles/RMH200"
params.RMSfusion = "${params.resources}/hg19/bedfiles/RMSfusion"
params.TST170 = "${params.resources}/hg19/bedfiles/TST170.v1.0"
params.TwistExome = "${params.resources}/hg19/bedfiles/Twist_Exome"
params.wgs = "${params.resources}/hg19/bedfiles/ens61_human_chromosomes_and_MT.bedToolsGenomeFile"
params.WGS_gnomad = "${params.resources}/hg19/bedfiles/ens61_human_chromosomes_and_MT.gnomad.withaf"


//test.subscribe { println "value: $it" }

//^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
//
//                     PROCESSES
//
//vvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvv



//process get_sample_info_align {
//
//	input:
//	val info from samplesheet_align
//
//	output:
//	set sample_id, fastqs, seq_type, umi, read_structures, expected_umi_length into sample_info_align
//	
//	script:
//
//	sample_id = info['sample_id']
//        fastqs = info['fastq_path'].sort().join(' ')
//        seq_type = info['seq_type']
//        umi = info['umi']
//
//        if (umi == 'IDT') {
//                read_structures = '+T +M +T'
//                expected_umi_length = 9 }
//        if (umi == 'AVENIO') {
//                read_structures = '6M+T 6M+T'
//                expected_umi_length = 12 }
//        if (umi == 'Q33') {
//                read_structures = '+T 12M12S+T'
//                expected_umi_length = 12 }
//	
//	"""
//	"""
//}



process align {
	publishDir path: "${params.outdir}/Alignments", pattern: "${sample_id}.sort.ba?", mode: 'copy', overwrite: true
        publishDir path: "${params.outdir}/Stats", pattern: '*.metrics', mode: 'copy', overwrite: true
	publishDir path: "${params.outdir}/scripts", pattern: "map.${sample_id}.sh", mode: 'copy', overwrite: true
	
	tag "MAP_${sample_id}"

	input:
	val(info) from samplesheet_align

	output:
	file("bwa.${sample_id}.sh")
        set info, file("${sample_id}.sort.ba?") into alignOutput, alignOutput1
        set info, file('*.metrics') into alignOutput_metrics

	script:
	sample_id = info['sample_id']
        fastqs = info['fastq_path'].sort().join(' ')
        pool_id = info['pool_id']
        seq_type = info['seq_type']
	umi = info['umi']
	seq_plat = info['seq_plat']
	
	barcode_seq = 'NNNNNNN'
	fastq_dir = "${baseDir}/data/fastq/${pool_id}"
	rm_dups = seq_type == "amplicon" ? 'false' : 'true'
	num_cpus = seq_type == "cfdna" ? params.cfdna_cpus : params.bwa_cpus

	"""
	perl ${params.scriptDir}/bwa.pl \
		--output_script bwa.${sample_id}.sh \
		--java_tmp ${params.java_tmp} \
		--picard ${params.picard_jar} \
		--num_cpu ${num_cpus} \
		--ref_fasta ${params.reference} \
		--analysis_dir . \
		--fastq_dir ${fastq_dir} \
		--samp_id ${sample_id} \
		--seq_plat ${seq_plat} \
		--barcode_seq ${barcode_seq} \
		--umi ${umi} \
		--seq_type ${seq_type} \
		--run_type none \
		--rm_dups ${rm_dups}

	. bwa.${sample_id}.sh
	"""
}


process qc {
	publishDir path: "${params.outdir}/Stats", pattern: "${sample_id}.*", mode: 'copy', overwrite: true
	publishDir path: "${params.outdir}/scripts", pattern: "qc.${sample_id}.sh",  mode: 'copy', overwrite: true

	tag "QC_${sample_id}"

	input:
	set val(info), file(infileBam) from alignOutput

	output:	
	file("qc.${sample_id}.sh")
	tuple val(info), file("${sample_id}.*") into qcOutput
	
        script:
        sample_id = info['sample_id']
        seq_type = info['seq_type']
        tcontent = info['tumour_content']
	bedfile = info["bedfile"]	

	if ( tcontent == null ) {
		tcontent = 0.5 
	}
	
	if ( bedfile != null ) {
		interval = "${bedfile}.CoverageCalculator.bed.intervals"
	} else {
		interval = 'wgs'
	}

	"""
	perl ${params.scriptDir}/qc.pl \
		--output_script qc.${sample_id}.sh \
		--ref_fasta ${params.reference} \
		--java_tmp ${params.java_tmp} \
		--picard ${params.picard_jar} \
		--getexoncoverage ${params.scriptDir}/get_exon_coverage.pl \
		--analysis_dir . \
		--samp_id ${sample_id} \
		--interval_file ${interval} \
		--seq_type ${seq_type} \
		--tcontent ${tcontent} \
		--MID_file ${params.MolecularID}
	
	. qc.${sample_id}.sh
	"""
}


// SNV module currently cannot handle wgs sequencing in such a way that it will split jobs by chromosome
// SNV module currenlty cannot handle multiple tumour samples per control. Either tumour-only or 1 tumour per normal sample is supported

process snv {
	publishDir path: "${params.outdir}/Variants/${sample_id}", pattern: "{*${sample_id}.vcf*, 1.SNV.${sample_id}.bamout.ba?, 1.${sample_id}.vardict*}", mode: 'copy', overwrite: true
        publishDir path: "${params.outdir}/scripts", pattern: "snv.${sample_id}.sh",  mode: 'copy', overwrite: true
        tag "SNV_${sample_id}"

	when:
	seq_type != "lcwgs"

	input:
	tuple val(info), file(infileBam) from alignOutput1

	output:
	file("snv.${sample_id}.sh")
	file("{*.${sample_id}.vcf*, 1.SNV.${sample_id}.bamout.ba?, 1.${sample_id}.vardict*}") into snvOutput

	script:
	sample_id = info["sample_id"]
	gatk_grp = info["gatk_grp"]
        seq_type = info["seq_type"]
	tag = info["tag"]
	bedfile = info["bedfile"]
	
	count = 1
	control = "None"
	intervals = "${bedfile}.bed.intervals"

	if ( tag == "normal" || tag == "control" ) {
		control = gatk_grp
	} else {
		tumour = gatk_grp
	}
	if ( seq_type == 'cfdna' ) {
		mutect_cpus = 10
		gatk_mem = 55g
	} else {
		mutect_cpus = params.mutect_cpus
		gatk_mem = params.gatk_mem
	}

	"""
	perl ${params.scriptDir}/snv.pl \
		--output_script snv.${gatk_grp}.sh \
		--java_tmp ${params.java_tmp} \
		--picard ${params.picard_jar} \
		--gatk ${params.gatk_engine} \
		--ref_fasta ${params.reference} \
		--analysis_dir . \
		--control ${control} \
		--tumour ${tumour} \
		--dbsnp_vcf ${params.dbsnp} \
		--indel_vcf ${params.indel} \
		--PON ${params.PON_WGS} \
		--mem ${params.gatk_mem} \
		--num_cpus ${mutect_cpus} \
		--chr ${intervals} \
		--seqtype ${seq_type} \
		--id ${count}

	. snv.${sample_id}.sh
	"""

}



//process align {
//
//	tag {sample_id}
//        publishDir path: "${params.outdir}/Alignments", pattern: '*.sort.ba?', mode: 'copy', overwrite: true
//        publishDir path: "${params.outdir}/Stats", pattern: '*.metrics', mode: 'copy', overwrite: true
//
//	when:
//	params.align
//
//        input:
////        val info from sample_info_align
//	set sample_id, fastqs, seq_type, umi, read_structures, expected_umi_length from sample_info_align	
//
//        output:
//        set info, file('*.sort.ba?') into alignOutput
//        set info, file('*.metrics') into alignOutput_metrics
//
//        script:
////        sample_id = info['sample_id']
////        fastqs = info['fastq_path'].sort().join(' ')
////        pool_id = info['pool_id']
////        seq_type = info['seq_type']
////	umi = info['umi']
////
////	if (umi == 'IDT') {
////                read_structures = '+T +M +T'
////                expected_umi_length = 9 }
////        if (umi == 'AVENIO') {
////                read_structures = '6M+T 6M+T'
////                expected_umi_length = 12 }
////        if (umi == 'Q33') {
////                read_structures = '+T 12M12S+T'
////                expected_umi_length = 12 }
//
//
//	if ( seq_type == "capture" | seq_type == "exome" )
//		"""
//		bwa mem -M -c 1 -t ${params.threads} -k ${params.seed} -p \
//        	${params.reference} ${fastqs} > ${sample_id}.aln.bam
//
//        	${params.picard} SortSam \
//        	        SORT_ORDER=queryname \
//        	        INPUT=${sample_id}.aln.bam \
//        	        OUTPUT=${sample_id}.aln.sorted.bam \
//        	        CREATE_INDEX=true \
//        	        VALIDATION_STRINGENCY=LENIENT \
//        	        TMP_DIR=${params.java_tmp}
//
//        	${params.picard} MarkDuplicates \
//        	        REMOVE_DUPLICATES=true \
//        	        INPUT=${sample_id}.aln.sorted.bam \
//        	        OUTPUT=${sample_id}.aln.sorted.marked.bam \
//        	        METRICS_FILE=${sample_id}.mark_dups.metrics \
//        	        ASSUME_SORTED=true \
//        	        CREATE_INDEX=true \
//        	        TMP_DIR=${params.java_tmp} \
//        	        VALIDATION_STRINGENCY=LENIENT
//
//        	${params.fgbio} ClipBam \
//        	        --input ${sample_id}.aln.sorted.marked.bam \
//        	        --output ${sample_id}.sort.bam \
//        	        --metrics ${sample_id}.clipbam.metrics \
//        	        --ref ${params.reference} \
//        	        --clip-overlapping-reads=true \
//        	        --clipping-mode Hard	
//		"""
//
//	else if ( seq_type == "umi" )
//		"""
//		${params.fgbio} FastqToBam \
//               		--input ${fastqs} \
//               		--read-structures ${read_structures} \
//               		--output ${sample_id}.umi.unaln.bam \
//               		--sort true \
//               		--umi-tag RX \
//               		--sample ${sample_id} \
//               		--library ${sample_id} \
//               		--read-group-id ${sample_id}
//
//	        samtools index ${sample_id}.umi.unaln.bam
//	
//	        samtools view -h ${sample_id}.umi.unaln.bam > ${sample_id}.umi.unaln.inter.sam
//	        sed 's/RX:Z:\\([A-Z]*\\)-\\([A-Z]*\\)/RX:Z:\\1\\2/' ${sample_id}.umi.unaln.inter.sam > ${sample_id}.umi.unaln.sam
//	        samtools view -S -b ${sample_id}.umi.unaln.sam > ${sample_id}.umi.unaln.edit.bam
//	
//	        ${params.picard} MarkIlluminaAdapters \
//	                INPUT=${sample_id}.umi.unaln.edit.bam \
//	                OUTPUT=${sample_id}.umi.unaln.marked.bam \
//	                METRICS=${sample_id}.adapter.metrics
//	
//	        ${params.picard} SamToFastq \
//	                INPUT=${sample_id}.umi.unaln.marked.bam \
//	                FASTQ=${sample_id}.umi.unaln.marked.bam.fastq \
//	                CLIPPING_ATTRIBUTE=XT \
//	                CLIPPING_ACTION=2 \
//	                INTERLEAVE=true \
//	                NON_PF=true \
//	                TMP_DIR=${params.java_tmp}
//	
//	        bwa mem -M -c 1 -t ${params.threads} -k ${params.seed} -p \
//	                ${params.reference} ${sample_id}.umi.unaln.marked.bam.fastq > ${sample_id}.umi.marked.aln.bam
//	
//	        ${params.picard} MergeBamAlignment \
//	                R=${params.reference} \
//	                UNMAPPED_BAM=${sample_id}.umi.unaln.marked.bam \
//	                ALIGNED_BAM=${sample_id}.umi.marked.aln.bam \
//	                OUTPUT=${sample_id}.umi.marked.aln.merged.bam \
//	                CREATE_INDEX=true \
//	                ADD_MATE_CIGAR=true \
//	                CLIP_ADAPTERS=false \
//	                CLIP_OVERLAPPING_READS=true \
//	                INCLUDE_SECONDARY_ALIGNMENTS=true \
//	                MAX_INSERTIONS_OR_DELETIONS=-1 \
//	                PRIMARY_ALIGNMENT_STRATEGY=MostDistant \
//	                ATTRIBUTES_TO_RETAIN=XS \
//	                TMP_DIR=${params.java_tmp}
//	
//	        ${params.fgbio} GroupReadsByUmi \
//	                --input ${sample_id}.umi.marked.aln.merged.bam \
//	                --family-size-histogram ${sample_id}.umi.metrics \
//	                --strategy adjacency \
//	                --min-map-q 30 \
//	                --raw-tag RX \
//	                --assign-tag MI \
//	                --min-umi-length ${expected_umi_length} \
//	                --output ${sample_id}.grouped.bam
//
//	 	${params.picard} SortSam \
//         	       SORT_ORDER=coordinate \
//         	       INPUT=${sample_id}.grouped.bam \
//         	       OUTPUT=${sample_id}.umi.marked.aln.merged.sorted.bam \
//         	       CREATE_INDEX=true \
//         	       VALIDATION_STRINGENCY=LENIENT \
//         	       TMP_DIR=${params.java_tmp}
//
//        	${params.picard} MarkDuplicates \
//         	       REMOVE_DUPLICATES=true \
//         	       INPUT=${sample_id}.umi.marked.aln.merged.sorted.bam \
//         	       OUTPUT=${sample_id}.sort.bam \
//         	       METRICS_FILE=${sample_id}.mark_dups.metrics \
//         	       ASSUME_SORTED=true \
//         	       CREATE_INDEX=true \
//         	       TMP_DIR=${params.java_tmp} \
//         	       VALIDATION_STRINGENCY=LENIENT
//		"""	
//	
//	else if ( seq_type == "cfdna" )
//		"""
//		${params.fgbio} FastqToBam \
//                        --input ${fastqs} \
//                        --read-structures ${read_structures} \
//                        --output ${sample_id}.umi.unaln.bam \
//                        --sort true \
//                        --umi-tag RX \
//                        --sample ${sample_id} \
//                        --library ${sample_id} \
//                        --read-group-id ${sample_id}
//
//                samtools index ${sample_id}.umi.unaln.bam
//
//                samtools view -h ${sample_id}.umi.unaln.bam > ${sample_id}.umi.unaln.inter.sam
//                sed 's/RX:Z:\\([A-Z]*\\)-\\([A-Z]*\\)/RX:Z:\\1\\2/' ${sample_id}.umi.unaln.inter.sam > ${sample_id}.umi.unaln.sam
//                samtools view -S -b ${sample_id}.umi.unaln.sam > ${sample_id}.umi.unaln.edit.bam
//
//                ${params.picard} MarkIlluminaAdapters \
//                        INPUT=${sample_id}.umi.unaln.edit.bam \
//                        OUTPUT=${sample_id}.umi.unaln.marked.bam \
//                        METRICS=${sample_id}.adapter.metrics
//
//                ${params.picard} SamToFastq \
//                        INPUT=${sample_id}.umi.unaln.marked.bam \
//                        FASTQ=${sample_id}.umi.unaln.marked.bam.fastq \
//                        CLIPPING_ATTRIBUTE=XT \
//                        CLIPPING_ACTION=2 \
//                        INTERLEAVE=true \
//                        NON_PF=true \
//                        TMP_DIR=${params.java_tmp}
//
//                bwa mem -M -c 1 -t ${params.threads} -k ${params.seed} -p \
//                        ${params.reference} ${sample_id}.umi.unaln.marked.bam.fastq > ${sample_id}.umi.marked.aln.bam
//
//                ${params.picard} MergeBamAlignment \
//                        R=${params.reference} \
//                        UNMAPPED_BAM=${sample_id}.umi.unaln.marked.bam \
//                        ALIGNED_BAM=${sample_id}.umi.marked.aln.bam \
//                        OUTPUT=${sample_id}.umi.marked.aln.merged.bam \
//                        CREATE_INDEX=true \
//                        ADD_MATE_CIGAR=true \
//                        CLIP_ADAPTERS=false \
//                        CLIP_OVERLAPPING_READS=true \
//                        INCLUDE_SECONDARY_ALIGNMENTS=true \
//                        MAX_INSERTIONS_OR_DELETIONS=-1 \
//                        PRIMARY_ALIGNMENT_STRATEGY=MostDistant \
//                        ATTRIBUTES_TO_RETAIN=XS \
//                        TMP_DIR=${params.java_tmp}
//
//                ${params.fgbio} GroupReadsByUmi \
//                        --input ${sample_id}.umi.marked.aln.merged.bam \
//                        --family-size-histogram ${sample_id}.umi.metrics \
//                        --strategy adjacency \
//                        --min-map-q 30 \
//                        --raw-tag RX \
//                        --assign-tag MI \
//                        --min-umi-length ${expected_umi_length} \
//                        --output ${sample_id}.grouped.bam
//
//		${params.fgbio} CallMolecularConsensusReads \
//	                --input ${sample_id}.grouped.bam \
//	                --output ${sample_id}.umi.marked.aln.merged.sort.consensus.bam \
//	                --min-reads ${params.min_reads} \
//	                --min-input-base-quality 30 \
//	                --tag MI
//	
//	        ${params.fgbio} FilterConsensusReads \
//	                --input ${sample_id}.umi.marked.aln.merged.sort.consensus.bam \
//	                --output ${sample_id}.umi.marked.aln.merged.sort.consensus.filtered.bam \
//	                --ref ${params.reference} \
//	                --min-reads 3 \
//	                --max-read-error-rate 0.05 \
//	                --min-base-quality 40 \
//	                --max-base-error-rate 0.1 \
//	                --max-no-call-fraction 0.05 \
//	                --reverse-per-base-tags true
//	
//	        ${params.picard} SortSam \
//	                SORT_ORDER=queryname \
//	                INPUT=${sample_id}.umi.marked.aln.merged.sort.consensus.filtered.bam \
//	                OUTPUT=${sample_id}.umi.marked.aln.merged.sort.consensus.filtered.sorted.bam \
//	                CREATE_INDEX=true \
//	                VALIDATION_STRINGENCY=LENIENT \
//	                TMP_DIR=${params.java_tmp}
//	
//	        ${params.picard} SamToFastq \
//	                INPUT=${sample_id}.umi.marked.aln.merged.sort.consensus.filtered.sorted.bam \
//	                FASTQ=${sample_id}.umi.marked.aln.merged.sort.consensus.filtered.sorted.fastq \
//	                CLIPPING_ATTRIBUTE=XT \
//	                CLIPPING_ACTION=2 \
//	                INTERLEAVE=true \
//	                NON_PF=true \
//	                TMP_DIR=${params.java_tmp}
//	
//	        bwa mem -M -c 1 -t ${params.threads} -k ${params.seed} -p \
//	               ${params.reference} ${sample_id}.umi.marked.aln.merged.sort.consensus.filtered.sorted.fastq > ${sample_id}.consensus.aln.bam \
//	
//	        ${params.picard} SortSam \
//	               SORT_ORDER=queryname \
//	               INPUT=${sample_id}.consensus.aln.bam \
//	               OUTPUT=${sample_id}.consensus.aln.sorted.bam \
//	               CREATE_INDEX=true \
//	               VALIDATION_STRINGENCY=LENIENT \
//	               TMP_DIR=${params.java_tmp}
//	
//	        ${params.picard} MergeBamAlignment \
//	               R=${params.reference} \
//	               UNMAPPED_BAM=${sample_id}.umi.marked.aln.merged.sort.consensus.filtered.sort.bam \
//	               ALIGNED_BAM=${sample_id}.consensus.aln.sorted.bam \
//	               OUTPUT=${sample_id}.consensus.aln.sorted.merged.bam \
//	               CREATE_INDEX=true \
//	               ADD_MATE_CIGAR=true \
//	               CLIP_ADAPTERS=false \
//	               CLIP_OVERLAPPING_READS=true \
//	               INCLUDE_SECONDARY_ALIGNMENTS=true \
//	               MAX_INSERTIONS_OR_DELETIONS=-1 \
//	               PRIMARY_ALIGNMENT_STRATEGY=MostDistant \
//	               ATTRIBUTES_TO_RETAIN=XS \
//	               TMP_DIR=${params.java_tmp}
//	
//	        ${params.fgbio} ClipBam \
//	               --input ${sample_id}.consensus.aln.sorted.merged.bam \
//	               --output ${sample_id}.sort.bam \
//	               --metrics ${sample_id}.clipbam.metrics \
//	               --ref ${params.reference} \
//	               --clip-overlapping-reads=true \
//	               --clipping-mode Hard
//		"""
//
//
//}
//
//
//process align_umi1 {
//        publishDir path: "${params.outdir}/Alignments", pattern: '*.metrics', mode: 'copy', overwrite: true
//	publishDir path: "${params.outdir}/scripts/", pattern: ".command.sh", mode: 'copy', overwrite: true, saveAs: "align.${sample_id}.sh"
//	tag {sample_id}
//
//        when:
//        info['seq_type'] == 'cfdna' | info['seq_type'] == 'umi'
//
//        input:
//        val info from sample_info_channel_umi
//
//        output:
//        set sample_id, seq_type, umi, file('*.grouped.ba?') into alignUmiGroupOutput_amplicon, alignUmiGroupOutput_cfdna
//
//        script:
//        sample_id = info['sample_id']
//        fastqs = info['fastq_path'].sort().join(' ')
//        seq_type = info['seq_type']
//        umi = info['umi']
//
//	if (umi == 'IDT') {
//                read_structures = '+T +M +T' 
//		expected_umi_length = 9 }
//        if (umi == 'AVENIO') {
//                read_structures = '6M+T 6M+T' 
//		expected_umi_length = 12 }
//        if (umi == 'Q33') {
//                read_structures = '+T 12M12S+T' 
//		expected_umi_length = 12 }
//
//	
//	"""
//	${params.fgbio} FastqToBam \
//               --input ${fastqs} \
//               --read-structures ${read_structures} \
//               --output ${sample_id}.umi.unaln.bam \
//               --sort true \
//               --umi-tag RX \
//               --sample ${sample_id} \
//               --library ${sample_id} \
//               --read-group-id ${sample_id}
//
//	samtools index ${sample_id}.umi.unaln.bam
//	
//	samtools view -h ${sample_id}.umi.unaln.bam > ${sample_id}.umi.unaln.inter.sam
//        sed 's/RX:Z:\\([A-Z]*\\)-\\([A-Z]*\\)/RX:Z:\\1\\2/' ${sample_id}.umi.unaln.inter.sam > ${sample_id}.umi.unaln.sam
//	samtools view -S -b ${sample_id}.umi.unaln.sam > ${sample_id}.umi.unaln.edit.bam
//
//	${params.picard} MarkIlluminaAdapters \
//                INPUT=${sample_id}.umi.unaln.edit.bam \
//                OUTPUT=${sample_id}.umi.unaln.marked.bam \
//                METRICS=${sample_id}.adapter.metrics
//
//	${params.picard} SamToFastq \
//		INPUT=${sample_id}.umi.unaln.marked.bam \
//		FASTQ=${sample_id}.umi.unaln.marked.bam.fastq \
//		CLIPPING_ATTRIBUTE=XT \
//		CLIPPING_ACTION=2 \
//		INTERLEAVE=true \
//		NON_PF=true \
//		TMP_DIR=${params.java_tmp}
//
//	bwa mem -M -c 1 -t ${params.threads} -k ${params.seed} -p \
//		${params.reference} ${sample_id}.umi.unaln.marked.bam.fastq > ${sample_id}.umi.marked.aln.bam
//
//	${params.picard} MergeBamAlignment \
//	        R=${params.reference} \
//	        UNMAPPED_BAM=${sample_id}.umi.unaln.marked.bam \
//	        ALIGNED_BAM=${sample_id}.umi.marked.aln.bam \
//	        OUTPUT=${sample_id}.umi.marked.aln.merged.bam \
//	        CREATE_INDEX=true \
//	        ADD_MATE_CIGAR=true \
//	        CLIP_ADAPTERS=false \
//	        CLIP_OVERLAPPING_READS=true \
//	        INCLUDE_SECONDARY_ALIGNMENTS=true \
//	        MAX_INSERTIONS_OR_DELETIONS=-1 \
//	        PRIMARY_ALIGNMENT_STRATEGY=MostDistant \
//	        ATTRIBUTES_TO_RETAIN=XS \
//	        TMP_DIR=${params.java_tmp}
//
//	${params.fgbio} GroupReadsByUmi \
//		--input ${sample_id}.umi.marked.aln.merged.bam \
//		--family-size-histogram ${sample_id}.umi.metrics \
//		--strategy adjacency \
//		--min-map-q 30 \
//		--raw-tag RX \
//		--assign-tag MI \
//		--min-umi-length ${expected_umi_length} \
//		--output ${sample_id}.grouped.bam
//	"""
//}
//
//
//process align_umi2 {
//
//	tag {sample_id}
//	publishDir path: "${params.outdir}/scripts/", pattern: ".command.sh", mode: 'copy', overwrite: true, saveAs: "align.${sample_id}.sh"
//        publishDir path: "${params.outdir}/Alignments", pattern: '*.{sort.ba?,metricsi,.command.sh}', mode: 'copy', overwrite: true
//
//        when:
//        seq_type == 'umi'
//
//        input:
//        set sample_id, seq_type, umi, file(infileBam) from alignUmiGroupOutput_amplicon
//
//        output:
//        set sample_id, file('*.sort.ba?') into alignUmiOutput
//
//        script:
//	
//	if (umi == 'IDT') {
//                expected_umi_length = 9 }
//
//        if (umi == 'AVENIO') {
//                expected_umi_length = 12 }
//
//        if (umi == 'Q33') {
//                expected_umi_length = 12 }
//	
//	"""
//	${params.picard} SortSam \
//                SORT_ORDER=coordinate \
//                INPUT=${infileBam} \
//                OUTPUT=${sample_id}.umi.marked.aln.merged.sorted.bam \
//                CREATE_INDEX=true \
//                VALIDATION_STRINGENCY=LENIENT \
//                TMP_DIR=${params.java_tmp}
//
//        ${params.picard} MarkDuplicates \
//                REMOVE_DUPLICATES=true \
//                INPUT=${sample_id}.umi.marked.aln.merged.sorted.bam \
//                OUTPUT=${sample_id}.sort.bam \
//                METRICS_FILE=${sample_id}.mark_dups.metrics \
//                ASSUME_SORTED=true \
//                CREATE_INDEX=true \
//                TMP_DIR=${params.java_tmp} \
//                VALIDATION_STRINGENCY=LENIENT
//
//	"""
//}
//
//
//process align_consensus {
//	tag {sample_id}
//	publishDir path: "${params.outdir}/scripts/", pattern: ".command.sh", mode: 'copy', overwrite: true, saveAs: "align.${sample_id}.sh"
//        publishDir path: "${params.outdir}/Alignments", pattern: '*.{sort.ba?,metrics}', mode: 'copy', overwrite: true
//
//        when:
//        seq_type == 'cfdna'
//
//        input:
//        set sample_id, seq_type, umi, file(infileBam) from alignUmiGroupOutput_cfdna
//
//        output:
//        set info, file('*.sort.ba?') into alignConsensusOutput
//
//	script:
//	"""
//	${params.fgbio} CallMolecularConsensusReads \
//                --input ${infileBam} \
//                --output ${sample_id}.umi.marked.aln.merged.sort.consensus.bam \
//                --min-reads ${params.min_reads} \
//                --min-input-base-quality 30 \
//                --tag MI
//
//        ${params.fgbio} FilterConsensusReads \
//                --input ${sample_id}.umi.marked.aln.merged.sort.consensus.bam \
//                --output ${sample_id}.umi.marked.aln.merged.sort.consensus.filtered.bam \
//                --ref ${params.reference} \
//                --min-reads 3 \
//                --max-read-error-rate 0.05 \
//                --min-base-quality 40 \
//                --max-base-error-rate 0.1 \
//                --max-no-call-fraction 0.05 \
//                --reverse-per-base-tags true
//
//        ${params.picard} SortSam \
//                SORT_ORDER=queryname \
//                INPUT=${sample_id}.umi.marked.aln.merged.sort.consensus.filtered.bam \
//                OUTPUT=${sample_id}.umi.marked.aln.merged.sort.consensus.filtered.sorted.bam \
//                CREATE_INDEX=true \
//                VALIDATION_STRINGENCY=LENIENT \
//                TMP_DIR=${params.java_tmp}
//
//        ${params.picard} SamToFastq \
//                INPUT=${sample_id}.umi.marked.aln.merged.sort.consensus.filtered.sorted.bam \
//                FASTQ=${sample_id}.umi.marked.aln.merged.sort.consensus.filtered.sorted.fastq \
//                CLIPPING_ATTRIBUTE=XT \
//                CLIPPING_ACTION=2 \
//                INTERLEAVE=true \
//                NON_PF=true \
//                TMP_DIR=${params.java_tmp}
//
//        bwa mem -M -c 1 -t ${params.threads} -k ${params.seed} -p \
//               ${params.reference} ${sample_id}.umi.marked.aln.merged.sort.consensus.filtered.sorted.fastq > ${sample_id}.consensus.aln.bam \
//
//        ${params.picard} SortSam \
//               SORT_ORDER=queryname \
//               INPUT=${sample_id}.consensus.aln.bam \
//               OUTPUT=${sample_id}.consensus.aln.sorted.bam \
//               CREATE_INDEX=true \
//               VALIDATION_STRINGENCY=LENIENT \
//               TMP_DIR=${params.java_tmp}
//
//        ${params.picard} MergeBamAlignment \
//               R=${params.reference} \
//               UNMAPPED_BAM=${sample_id}.umi.marked.aln.merged.sort.consensus.filtered.sort.bam \
//               ALIGNED_BAM=${sample_id}.consensus.aln.sorted.bam \
//               OUTPUT=${sample_id}.consensus.aln.sorted.merged.bam \
//               CREATE_INDEX=true \
//               ADD_MATE_CIGAR=true \
//               CLIP_ADAPTERS=false \
//               CLIP_OVERLAPPING_READS=true \
//               INCLUDE_SECONDARY_ALIGNMENTS=true \
//               MAX_INSERTIONS_OR_DELETIONS=-1 \
//               PRIMARY_ALIGNMENT_STRATEGY=MostDistant \
//               ATTRIBUTES_TO_RETAIN=XS \
//               TMP_DIR=${params.java_tmp}
//
//        ${params.fgbio} ClipBam \
//               --input ${sample_id}.consensus.aln.sorted.merged.bam \
//               --output ${sample_id}.sort.bam \
//               --metrics ${sample_id}.clipbam.metrics \
//               --ref ${params.reference} \
//               --clip-overlapping-reads=true \
//               --clipping-mode Hard
//
//	"""
//}
//
//
//
//
//
//
//
//
//
//process qc {
//
//	publishDir path: "${params.outdir}/Stats", mode: 'move', overwrite: true
//	
//	input:
//	val info from sample_info_qc
//	val(bedfile) from bedfile_path
//	
//	when:
//	params.qc
//
//	output:
//	
//
//	script:
//        sample_id = info['sample_id']
//        seq_type = info['seq_type']
//        umi = info['umi']
//	target_bed = info['target_bed']
//	tcontent = info['tumour_content']
//	
//	infileBam = '${params.outdir}/Alignments/${sample_id}.sort.bam'
//	interval = "${bedfile}.CoverageCalculator.bed.intervals"
//	bed_intervls = "${bedfile}.bed.intervals"
//
//
//	if ( seq_type == "capture" | seq_type == "exome" ) {
//		"""
//	        $params.picard CalculateHsMetrics MAX_RECORDS_IN_RAM=500000 \
//	                REFERENCE_SEQUENCE=${params.reference} \
//	                INPUT=${infileBam} \
//	                OUTPUT=${sample_id}.qc.metrics \
//	                BAIT_INTERVALS=${interval} \
//	                TARGET_INTERVALS=${interval} \
//	                BAIT_SET_NAME=${sample_id} \
//	                METRIC_ACCUMULATION_LEVEL=ALL_READS \
//	                TMP_DIR=${params.java_tmp} \
//	                PER_TARGET_COVERAGE=${sample_id}.qc.coverage \
//	                PER_BASE_COVERAGE=${sample_id}.base.coverage \
//	                VALIDATION_STRINGENCY=LENIENT
//	
//	        $params.picard CalculateHsMetrics MAX_RECORDS_IN_RAM=500000 \
//	                REFERENCE_SEQUENCE=${params.reference} \
//	                INPUT=${infileBam} \
//	                OUTPUT=${sample_id}.qc.metrics.targetbed \
//	                BAIT_INTERVALS=${bed_interval} \
//	                TARGET_INTERVALS=${bed_interval} \
//	                BAIT_SET_NAME=${sample_id} \
//	                METRIC_ACCUMULATION_LEVEL=ALL_READS \
//	                TMP_DIR=${params.java_tmp} \
//	                PER_TARGET_COVERAGE=${sample_id}.qc.coverage.targetbed \
//	                PER_BASE_COVERAGE=${sample_id}.base.coverage.targetbed \
//	                VALIDATION_STRINGENCY=LENIENT
//
//		perl get_exon_coverage.pl --dir ${params.outdir}/Stats --samp ${sample_id} --tcontent ${tcontent}
//		"""	
//
//	} else if ( seq_type == "umi" | seq_type == "cfdna" | seq_type == "amplicon" ) {
//		"""
//		$params.picard CollectTargetedPcrMetrics MAX_RECORDS_IN_RAM=500000 \
//			REFERENCE_SEQUENCE=$params.reference \
//			INPUT=${sample_id}.sort.bam \
//			OUTPUT=${sample_id}.qc.metrics \
//			AMPLICON_INTERVALS=${interval} \
//			TARGET_INTERVALS=${interval} \
//			CUSTOM_AMPLICON_SET_NAME=${sample_id} \
//			METRIC_ACCUMULATION_LEVEL=ALL_READS \
//			TMP_DIR=${params.java_tmp} \
//			PER_TARGET_COVERAGE=${sample_id}.qc.coverage \
//			PER_BASE_COVERAGE=${sample_id}.base.coverage \
//			VALIDATION_STRINGENCY=LENIENT
//
//		perl get_exon_coverage.pl --dir ${params.outdir}/Stats --samp ${sample_id} --tcontent ${tcontent}
//		"""
//
//	else
//		"""
//		${params.picard} CollectWgsMetrics MAX_RECORDS_IN_RAM=500000 \
//			REFERENCE_SEQUENCE=${params.reference} \
//			INPUT=${sample_id}.sort.bam \
//			OUTPUT=${sample_id}.qc.metrics \
//			TMP_DIR=${params.java_tmp} \
//			VALIDATION_STRINGENCY=LENIENT
//
//		perl get_exon_coverage.pl --dir ${params.outdir}/Stats --samp ${sample_id} --tcontent ${tcontent}
//		"""
//}
