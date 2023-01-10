#!/bin/bash

if [ $# -ne 5 ]
then
    echo ""
    echo "Usage: $0 <hg38.wgs|hg38.wes|hg38.haloplex> <genome.fa> <output prefix> <sample1.read1.fq.gz> <sample1.read2.fq.gz>"
    echo ""
    exit -1
fi

SCRIPT=$(readlink -f "$0")
BASEDIR=$(dirname "$SCRIPT")
export PATH=${BASEDIR}/../conda/bin:${PATH}
source activate align

# CMD parameters
ATYPE=${1}
GENOME=${2}
OUTP=${3}
FQ1=${4}
FQ2=${5}
THREADS=8

# BWA alignment
bwa mem -R "@RG\tID:${OUTP}\tSM:${OUTP}" -t ${THREADS} ${GENOME} ${FQ1} ${FQ2} | samtools fixmate -m - - | samtools sort -o ${OUTP}.srt.bam -
samtools index ${OUTP}.srt.bam

# Mark duplicates
if [[ ${ATYPE} = *"haloplex"* ]]; then
    mv ${OUTP}.srt.bam ${OUTP}.bam
    mv ${OUTP}.srt.bam.bai ${OUTP}.bam.bai
else
    samtools markdup ${OUTP}.srt.bam ${OUTP}.bam
    samtools index ${OUTP}.bam
    rm ${OUTP}.srt.bam ${OUTP}.srt.bam.bai
fi

# Run stats using unfiltered BAM
samtools idxstats -@ ${THREADS} ${OUTP}.bam > ${OUTP}.idxstats
samtools flagstat -@ ${THREADS} ${OUTP}.bam > ${OUTP}.flagstat

# Run alfred for BAM statistics
alfred qc -b ${BASEDIR}/../genome/${ATYPE}.bed.gz -r ${GENOME} -j ${OUTP}.alfred.json.gz -o ${OUTP}.alfred.tsv.gz ${OUTP}.bam
