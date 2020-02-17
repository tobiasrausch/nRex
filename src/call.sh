#!/bin/bash

if [ $# -lt 4 ]
then
    echo ""
    echo "Usage: $0 <hg19.wgs|hg19.wes|hg19.haloplex> <genome.fa> <output prefix> <sample1.bam> <sample2.bam> ... <sampleN.bam>"
    echo ""
    exit -1
fi

SCRIPT=$(readlink -f "$0")
BASEDIR=$(dirname "$SCRIPT")
FBY=1 # Freebayes or Strelka

# Add all required binaries
export PATH=/g/funcgen/bin:${PATH}
STRELKAROOT=/g/solexa/bin/software/strelka-2.9.2.centos6_x86_64/

# CMD parameters
ATYPE=${1}
GENOME=${2}
OUTP=${3}
shift 3
THREADS=4

# Redirect tmp to scratch directory if available
if [ -n "${SCRATCHDIR}" ]
then
    export TMP=${SCRATCHDIR}
fi

# Freebayes or Strelka
if [ ${FBY} -eq 1 ]
then
    freebayes --no-partial-observations --min-repeat-entropy 1 --report-genotype-likelihood-max --min-alternate-fraction 0.15 --fasta-reference ${GENOME} --genotype-qualities $@ -v ${OUTP}.vcf
    #freebayes --no-partial-observations --min-repeat-entropy 1 --report-genotype-likelihood-max --min-alternate-fraction 0.01 --fasta-reference ${GENOME} --genotype-qualities $@ -v ${OUTP}.vcf
else
    if [[ ${ATYPE} = *"wes"* ]]; then
	${STRELKAROOT}/bin/configureStrelkaGermlineWorkflow.py --exome --bam `echo ${@} | sed 's/[ \t][ \t]*/ --bam /g'` --referenceFasta ${GENOME} --runDir ${OUTP}_strelkaWES
    else
	${STRELKAROOT}/bin/configureStrelkaGermlineWorkflow.py --bam `echo ${@} | sed 's/[ \t][ \t]*/ --bam /g'` --referenceFasta ${GENOME} --runDir ${OUTP}_strelkaWES
    fi
    ${OUTP}_strelkaWES/runWorkflow.py -m local -j ${THREADS}
    bcftools view -f PASS ${OUTP}_strelkaWES/results/variants/variants.vcf.gz > ${OUTP}.vcf
    rm -rf ${OUTP}_strelkaWES
fi
cat ${OUTP}.vcf | grep "^#" > ${OUTP}.vcf.tmp
cat ${OUTP}.vcf| grep -v "^#" | sort -k1,1V -k2,2n >> ${OUTP}.vcf.tmp
mv ${OUTP}.vcf.tmp ${OUTP}.vcf
bgzip ${OUTP}.vcf
tabix ${OUTP}.vcf.gz

# Normalize VCF
bcftools norm -O z -o ${OUTP}.norm.vcf.gz -f ${GENOME} -m -both ${OUTP}.vcf.gz
tabix ${OUTP}.norm.vcf.gz
rm ${OUTP}.vcf.gz ${OUTP}.vcf.gz.tbi

# Filtering
if [ ${FBY} -eq 1 ]
then
    # Fixed threshold filtering
    if [[ ${ATYPE} = *"haloplex"* ]]; then
	# Stringent parameters
	bcftools filter -O z -o ${OUTP}.norm.filtered.vcf.gz -e '%QUAL<=20 || %QUAL/INFO/AO<=2 || SAF<=2 || SAR<=2' ${OUTP}.norm.vcf.gz
	# Relaxed parameters
	#bcftools filter -O z -o ${OUTP}.norm.filtered.vcf.gz -e '%QUAL<=20 || %QUAL/INFO/AO<=2 || SAF<=1 || SAR<=1' ${OUTP}.norm.vcf.gz
    else
	# Stringent parameters
	bcftools filter -O z -o ${OUTP}.norm.filtered.vcf.gz -e '%QUAL<=20 || %QUAL/INFO/AO<=2 || SAF<=2 || SAR<=2 || RPR<=2 || RPL<=2' ${OUTP}.norm.vcf.gz
	# Relaxed parameters
	#bcftools filter -O z -o ${OUTP}.norm.filtered.vcf.gz -e '%QUAL<=20 || %QUAL/INFO/AO<=2 || SAF<=1 || SAR<=1 || RPR<=0 || RPL<=0' ${OUTP}.norm.vcf.gz
    fi
else 
    bcftools filter -O z -o ${OUTP}.norm.filtered.vcf.gz -e '%QUAL<=20' ${OUTP}.norm.vcf.gz
fi
tabix ${OUTP}.norm.filtered.vcf.gz
rm ${OUTP}.norm.vcf.gz ${OUTP}.norm.vcf.gz.tbi

# Subset to target regions
bcftools view -T <(zcat ${BASEDIR}/../R/${ATYPE}.bed.gz) ${OUTP}.norm.filtered.vcf.gz | bgzip > ${OUTP}.${ATYPE}.vcf.gz
tabix ${OUTP}.${ATYPE}.vcf.gz
