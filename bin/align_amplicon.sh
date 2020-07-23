#!/bin/bash

java_tmp="tmp/"
java_cpu=9
java_mem=128
fgbio_jar="/apps/fgbio/0.6.1/fgbio-0.6.1.jar"
picard_jar="/apps/picard-tools/2.8.2/picard-2.8.2.jar"
picard="java -Xmx${java_mem}g -XX:ParallelGCThreads=${java_cpu} -jar ${picard_jar}"
fgbio="java -Xmx${java_mem}g -XX:+AggressiveOpts -XX:+AggressiveHeap -jar ${fgbio_jar}"

read_structures="+T 12M12S+T"


#####################################################################################

sample_id=$1
fq1=$2
fq2=$3
reference=$4


${fgbio} FastqToBam \
	--input ${fq1} ${fq2} \
        --read-structures ${read_structures} \
        --output ${sample_id}.bam \
        --sort true \
        --umi-tag RX \
        --sample ${sample_id} \
        --library ${sample_id} \
        --read-group-id ${sample_id}


${picard} MarkIlluminaAdapters \
        INPUT=${sample_id}.bam \
        OUTPUT=${sample_id}.marked.bam \
        METRICS=${sample_id}.mark_dups.metrics \
	TMP_DIR=${java_tmp}


${picard} SamToFastq \
        INPUT=${sample_id}.marked.bam \
        FASTQ=${sample_id}.marked.bam.fastq \
        CLIPPING_ATTRIBUTE=XT \
        CLIPPING_ACTION=2 \
        INTERLEAVE=true \
        NON_PF=true \
        TMP_DIR=${java_tmp}


bwa mem -M -c 1 -t ${threads} -k ${seed} -p \
        ${reference} ${sample_id}.marked.bam.fastq  > ${sample_id}.aln.bam 


${picard} MergeBamAlignment \
        R=${reference} \
        UNMAPPED_BAM=${sample_id}.marked.bam \
        ALIGNED_BAM=${sample_id}.aln.bam \
        OUTPUT=${sample_id}.aln.merged.bam \
        CREATE_INDEX=true \
        ADD_MATE_CIGAR=true \
        CLIP_ADAPTERS=false \
        CLIP_OVERLAPPING_READS=true \
        INCLUDE_SECONDARY_ALIGNMENTS=true \
        MAX_INSERTIONS_OR_DELETIONS=-1 \
        PRIMARY_ALIGNMENT_STRATEGY=MostDistant \
        ATTRIBUTES_TO_RETAIN=XS \
        TMP_DIR=${java_tmp}


${fgbio} GroupReadsByUmi \
        --input ${sample_id}.aln.merged.bam \
        --family-size-histogram ${sample_id}.umi.metrics \
        --strategy adjacency \
        --min-map-q 30 \
        --raw-tag RX \
        --assign-tag MI \
        --min-umi-length 9 \
        --output ${sample_id}.aln.merged.group.bam


${picard} SortSam \
        SORT_ORDER=queryname \
        INPUT=${sample_id}.aln.merged.group.bam \
        OUTPUT=${sample_id}.aln.merged.group.sorted.bam \
        CREATE_INDEX=true \
        VALIDATION_STRINGENCY=LENIENT \
        TMP_DIR=${java_tmp}


${picard} MarkDuplicates \
        REMOVE_DUPLICATES=true \
        INPUT=${sample_id}.aln.merged.group.sorted.bam \
        OUTPUT=${sample_id}.aln.merged.group.sorted.marked.bam \
        METRICS_FILE=${sample_id}.mark_dups.metrics \
        ASSUME_SORTED=true \
        CREATE_INDEX=true \
        TMP_DIR=${java_tmp} \
        VALIDATION_STRINGENCY=LENIENT
