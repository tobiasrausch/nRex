#!/bin/bash

if [ $# -lt 3 ]
then
    echo "Usage: $0 <genome.fa> <tumor.bam> <normal1.bam> ..."
    exit -1
fi

SCRIPT=$(readlink -f "$0")
BASEDIR=$(dirname "$SCRIPT")
GENOME=${1}
shift

# Programs
PICARD=${BASEDIR}/picard/dist/picard.jar
SAM=${BASEDIR}/samtools/samtools
BCF=${BASEDIR}/bcftools/bcftools
BAMADDRG=${BASEDIR}/bamaddrg/bamaddrg
FREEBAYES=${BASEDIR}/freebayes/bin/freebayes
GATK=/g/solexa/bin/software/GenomeAnalysisTK-2.5-2-gf57256b/GenomeAnalysisTK.jar

# Tmp directory
DSTR=$(date +'%a_%y%m%d_%H%M')
export TMP=/tmp/tmp_nRex_${DSTR}
mkdir -p ${TMP}
JAVAOPT="-Xmx32g -Djava.io.tmpdir=${TMP}"
PICARDOPT="MAX_RECORDS_IN_RAM=5000000 TMP_DIR=${TMP} VALIDATION_STRINGENCY=SILENT"

# Merge and remove duplicates
SAMPLEARR=( $@ )
BAMLIST=""
for ((  i = 0 ;  i < ${#SAMPLEARR[@]};  i++  ))
do
    SAMPLEID=`echo ${SAMPLEARR[$i]} | sed 's/^.*\///' | sed 's/\..*$//'`
    mkdir -p ${SAMPLEID}

    # Have we processed this sample already?
    if [ ! -f ${SAMPLEID}/${SAMPLEID}.rmdup.bam ]
    then
        # Clean .bam file
	java ${JAVAOPT} -jar ${PICARD} CleanSam I=${SAMPLEARR[$i]} O=${SAMPLEID}/${SAMPLEID}.clean.bam ${PICARDOPT}

        # Adjusting .bam file read groups to match GATK requirements
	java ${JAVAOPT} -jar ${PICARD} AddOrReplaceReadGroups I=${SAMPLEID}/${SAMPLEID}.clean.bam O=${SAMPLEID}/${SAMPLEID}.rg.bam SO=coordinate ID=${SAMPLEID} LB=${SAMPLEID}LB RGPL=illumina PU=laneX SM=${SAMPLEID} CN=nRex ${PICARDOPT}
	rm ${SAMPLEID}/${SAMPLEID}.clean.bam

        # Reorder resulting .bam file according to the reference
	java ${JAVAOPT} -jar ${PICARD} ReorderSam I=${SAMPLEID}/${SAMPLEID}.rg.bam O=${SAMPLEID}/${SAMPLEID}.sort.bam REFERENCE=${GENOME} ${PICARDOPT}
	rm ${SAMPLEID}/${SAMPLEID}.rg.bam

        # Remove duplicates
	java ${JAVAOPT} -jar ${PICARD} MarkDuplicates I=${SAMPLEID}/${SAMPLEID}.sort.bam O=${SAMPLEID}/${SAMPLEID}.rmdup.bam  M=${SAMPLEID}/${SAMPLEID}.metric.log PG=null MAX_FILE_HANDLES_FOR_READ_ENDS_MAP=1000 ${PICARDOPT}
	rm ${SAMPLEID}/${SAMPLEID}.sort.bam

        # Index resulting, sorted .bam file
	${SAM} index ${SAMPLEID}/${SAMPLEID}.rmdup.bam
    fi

    # Append BAM file
    BAMLIST=${BAMLIST}","${SAMPLEID}/${SAMPLEID}.rmdup.bam
done

# Mpileup
SAMPLEID=`echo ${SAMPLEARR[0]} | sed 's/^.*\///' | sed 's/\..*$//'`
if [ ! -f ${SAMPLEID}/${SAMPLEID}.mpileup.vcf ]
then
    ${SAM} mpileup -B -t DP,DV,SP -ugf ${GENOME} `echo ${BAMLIST} | sed 's/,/ /g'` | ${BCF} call -cv - > ${SAMPLEID}/${SAMPLEID}.mpileup.vcf
fi

# Freebayes
if [ ! -f ${SAMPLEID}/${SAMPLEID}.freebayes.vcf ]
then
    ${FREEBAYES} --min-alternate-fraction 0.15 --fasta-reference ${GENOME} --genotype-qualities `echo ${BAMLIST} | sed 's/,/ -b /g'` -v ${SAMPLEID}/${SAMPLEID}.freebayes.vcf
fi

# Optional: GATK's HaplotypeCaller
if [ -f ${GATK} ]
then
    if [ ! -f ${SAMPLEID}/${SAMPLEID}.gatk.vcf ]
    then
	java ${JAVAOPT} -jar ${GATK} -T HaplotypeCaller -R ${GENOME} `echo ${BAMLIST} | sed 's/,/ -I /g'` -o ${SAMPLEID}/${SAMPLEID}.gatk.vcf 
    fi
fi
