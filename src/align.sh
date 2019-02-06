#!/bin/bash

if [ $# -ne 5 ]
then
    echo ""
    echo "Usage: $0 <hg19.wgs|hg19.wes|hg19.haloplex> <genome.fa> <output prefix> <sample1.read1.fq.gz> <sample1.read2.fq.gz>"
    echo ""
    exit -1
fi

SCRIPT=$(readlink -f "$0")
BASEDIR=$(dirname "$SCRIPT")

# Add all required binaries
export PATH=/g/funcgen/bin:${PATH}

# CMD parameters
ATYPE=${1}
GENOME=${2}
OUTP=${3}
FQ1=${4}
FQ2=${5}
THREADS=4

# Redirect tmp to scratch directory if available
if [ -n "${SCRATCHDIR}" ]
then
    export TMP=${SCRATCHDIR}
fi

# Generate IDs
FQ1ID=`echo ${OUTP} | sed 's/$/.1.fq/'`
FQ2ID=`echo ${OUTP} | sed 's/$/.2.fq/'`

# Fastqc
mkdir -p ${OUTP}_prefastqc/ && fastqc -t ${THREADS} -o ${OUTP}_prefastqc/ ${FQ1} && fastqc -t ${THREADS} -o ${OUTP}_prefastqc/ ${FQ2}
    
# BWA
bwa mem -t ${THREADS} -R "@RG\tID:${OUTP}\tSM:${OUTP}" ${GENOME} ${FQ1} ${FQ2} | samtools view -bT ${GENOME} - > ${OUTP}.raw.bam

# Sort & Index
samtools sort -@ ${THREADS} -m 1536M -o ${OUTP}.srt.bam ${OUTP}.raw.bam && rm ${OUTP}.raw.bam && samtools index ${OUTP}.srt.bam

# Mark duplicates
if [[ ${ATYPE} = *"haloplex"* ]]; then
    mv ${OUTP}.srt.bam ${OUTP}.rmdup.bam
    mv ${OUTP}.srt.bam.bai ${OUTP}.rmdup.bam.bai
else
    bammarkduplicates markthreads=${THREADS} tmpfile=${OUTP}_`date +'%H%M%S'` I=${OUTP}.srt.bam O=${OUTP}.rmdup.bam M=${OUTP}.metrics.tsv index=1 rmdup=0
    rm ${OUTP}.srt.bam ${OUTP}.srt.bam.bai
fi
    
# Run stats using unfiltered BAM
samtools idxstats ${OUTP}.rmdup.bam > ${OUTP}.idxstats
samtools flagstat ${OUTP}.rmdup.bam > ${OUTP}.flagstat

# Run alfred for BAM statistics
alfred qc -b ${BASEDIR}/../R/${ATYPE}.bed.gz -r ${GENOME} -j ${OUTP}.alfred.json.gz -o ${OUTP}.alfred.tsv.gz ${OUTP}.rmdup.bam

# Filter duplicates, unmapped reads, chrM and unplaced contigs
CHRS=`zcat ${BASEDIR}/../R/${ATYPE}.bed.gz | cut -f 1 | sort -k1,1V -k2,2n | uniq | tr '\n' ' '`
samtools view -F 1024 -b ${OUTP}.rmdup.bam ${CHRS} > ${OUTP}.bam
samtools index ${OUTP}.bam
rm ${OUTP}.rmdup.bam ${OUTP}.rmdup.bam.bai
