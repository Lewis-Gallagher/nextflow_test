#!/bin/bash

seed=50
java_tmp="tmp/"
java_cpu=9
java_mem=128
fgbio_jar="/apps/fgbio/0.6.1/fgbio-0.6.1.jar"
picard_jar="/apps/picard-tools/2.8.2/picard-2.8.2.jar"
picard="java -Xmx${java_mem}g -XX:ParallelGCThreads=${java_cpu} -jar ${picard_jar}"
fgbio="java -Xmx${java_mem}g -XX:+AggressiveOpts -XX:+AggressiveHeap -jar ${fgbio_jar}"

#####################################################################################

sample_id=$1
fq1=$2
fq2=$3
reference=$4

bwa mem -M -c 1 -t ${bwa_threads} -k ${seed} -p \
        ${reference} ${fq1} ${fq2} > ${sample_id}.aln.bam

${picard} SortSam \
	SORT_ORDER=queryname \
        INPUT=${sample_id}.aln.bam \
        OUTPUT=${sample_id}.aln.sort.bam \
        CREATE_INDEX=true \
        VALIDATION_STRINGENCY=LENIENT \
        TMP_DIR=${java_tmp}

${picard} MarkDuplicates \
        REMOVE_DUPLICATES=true \
        INPUT=${sample_id}.aln.sort.bam  \
        OUTPUT=${sample_id}.aln.sort.marked.bam  \
        METRICS_FILE=${sample_id}.mark_dups.metrics \
        ASSUME_SORTED=true \
        CREATE_INDEX=true \
        TMP_DIR=${java_tmp} \
        VALIDATION_STRINGENCY=LENIENT

${fgbio} ClipBam \
        --input ${sample_id}.aln.sort.marked.bam \
        --output ${sample_id}.aln.sort.marked.bam \
        --metrics ${sample_id}.clipbam.metrics  \
        --ref ${reference} \
        --clip-overlapping-reads=true \
        --clipping-mode Hard


