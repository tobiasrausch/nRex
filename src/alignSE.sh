#!/bin/bash

if [ $# -ne 4 ]
then
    echo ""
    echo "Usage: $0 <hg19.wgs|hg19.wes|hg19.haloplex> <genome.fa> <output prefix> <sample.read.fq.gz>"
    echo ""
    exit -1
fi

SCRIPT=$(readlink -f "$0")
BASEDIR=$(dirname "$SCRIPT")
BWAALIGN=1  # BWA

# Add all required binaries
export PATH=/g/funcgen/bin:${PATH}

# CMD parameters
ATYPE=${1}
GENOME=${2}
OUTP=${3}
FQ1=${4}
THREADS=4

# Redirect tmp to scratch directory if available
if [ -n "${SCRATCHDIR}" ]
then
    export TMP=${SCRATCHDIR}
fi

# Generate IDs
FQ1ID=`echo ${OUTP} | sed 's/$/.1.fq/'`

# Fastqc
mkdir -p ${OUTP}_prefastqc/ && fastqc -t ${THREADS} -o ${OUTP}_prefastqc/ ${FQ1}
    
# BWA or Bowtie
if [ ${BWAALIGN} -eq 1 ]
then
    bwa mem -t ${THREADS} -R "@RG\tID:${OUTP}\tSM:${OUTP}" ${GENOME} ${FQ1} | samtools view -bT ${GENOME} - > ${OUTP}.raw.bam
else
    bowtie2 --threads ${THREADS} --rg-id "${OUTP}" --rg "SM:${OUTP}" --very-sensitive -x ${GENOME} -U ${FQ1} | samtools view -bT ${GENOME} - > ${OUTP}.raw.bam
fi

# Sort & Index
samtools sort -@ ${THREADS} -m 1536M -o ${OUTP}.srt.bam ${OUTP}.raw.bam && rm ${OUTP}.raw.bam && samtools index ${OUTP}.srt.bam

# Mark duplicates or deactivate for transposase assays
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
