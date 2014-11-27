#!/bin/bash

if [ $# -lt 3 ]
then
    echo "**********************************************************************"
    echo "nRex: Somatic single-nucleotide variant calling."
    echo "This program comes with ABSOLUTELY NO WARRANTY."
    echo ""
    echo "nRex (Version: 0.0.1)"
    echo "Contact: Tobias Rausch (rausch@embl.de)"
    echo "**********************************************************************"
    echo ""
    echo "Usage: $0 <genome.fa> <tumor.bam> <normal1.bam> ..."
    echo ""
    exit -1
fi

SCRIPT=$(readlink -f "$0")
BASEDIR=$(dirname "$SCRIPT")

# Wrapper script for mpileup, freebayes and GATK (optional)
${BASEDIR}/src/paired.sh $@

# Filter somatic SNVs
SAMPLEID=`echo ${2} | sed 's/^.*\///' | sed 's/\..*$//'`
if [ -f ${SAMPLEID}/${SAMPLEID}.mpileup.vcf ]
then
    echo "Mpileup somatic filtering"
    python ${BASEDIR}/python/somaticSNV.py -v ${SAMPLEID}/${SAMPLEID}.mpileup.vcf -o ${SAMPLEID}/${SAMPLEID}.somatic.mpileup.vcf
    cat ${SAMPLEID}/${SAMPLEID}.somatic.mpileup.vcf | grep -v "^#" | cut -f 1-5 | sort | uniq > ${SAMPLEID}/${SAMPLEID}.somatic.mpileup.bed
fi
if [ -f ${SAMPLEID}/${SAMPLEID}.freebayes.vcf ]
then
    echo "Freebayes somatic filtering"
    python ${BASEDIR}/python/somaticSNV.py -v ${SAMPLEID}/${SAMPLEID}.freebayes.vcf -g 0 -o ${SAMPLEID}/${SAMPLEID}.somatic.freebayes.vcf
    cat ${SAMPLEID}/${SAMPLEID}.somatic.freebayes.vcf | grep -v "^#" | cut -f 1-5 | sort | uniq > ${SAMPLEID}/${SAMPLEID}.somatic.freebayes.bed
fi
if [ -f ${SAMPLEID}/${SAMPLEID}.gatk.vcf ]
then
    echo "GATK somatic filtering"
    python ${BASEDIR}/python/somaticSNV.py -v ${SAMPLEID}/${SAMPLEID}.gatk.vcf -o ${SAMPLEID}/${SAMPLEID}.somatic.gatk.vcf
    cat ${SAMPLEID}/${SAMPLEID}.somatic.gatk.vcf | grep -v "^#" | cut -f 1-5 | sort | uniq > ${SAMPLEID}/${SAMPLEID}.somatic.gatk.bed
else
    cp ${SAMPLEID}/${SAMPLEID}.somatic.mpileup.bed ${SAMPLEID}/${SAMPLEID}.somatic.gatk.bed
fi

# Take the consensus
cat ${SAMPLEID}/${SAMPLEID}.somatic.mpileup.vcf | head -n 1000 | grep "^#" | ${BASEDIR}/src/htslib/bgzip > ${SAMPLEID}/${SAMPLEID}.nRex.vcf.gz
sort ${SAMPLEID}/${SAMPLEID}.somatic.mpileup.bed ${SAMPLEID}/${SAMPLEID}.somatic.freebayes.bed | uniq -d | sort - ${SAMPLEID}/${SAMPLEID}.somatic.gatk.bed | uniq -d > ${SAMPLEID}/${SAMPLEID}.grepConsensus
cat ${SAMPLEID}/${SAMPLEID}.somatic.mpileup.vcf | grep -v "^#" | grep -w -Ff ${SAMPLEID}/${SAMPLEID}.grepConsensus | ${BASEDIR}/src/htslib/bgzip >> ${SAMPLEID}/${SAMPLEID}.nRex.vcf.gz

# Clean-up
rm ${SAMPLEID}/${SAMPLEID}.grepConsensus ${SAMPLEID}/${SAMPLEID}.somatic.gatk.* ${SAMPLEID}/${SAMPLEID}.somatic.freebayes.* ${SAMPLEID}/${SAMPLEID}.somatic.mpileup.*

