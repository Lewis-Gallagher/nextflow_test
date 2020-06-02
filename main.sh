#!/bin/bash

#BSUB -J "mdx"
#BSUB -P DMPMXHAAC
#BSUB -n 4
#BSUB -R "span[hosts=1]"
#BSUB -o output.%J
#BSUB -e errors.%J

module add java/sun8
module add nextflow
module add samtools
module add bwa

# Initiate NextFlow Job
nextflow ./main.nf

