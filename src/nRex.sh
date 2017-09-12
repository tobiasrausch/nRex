#!/bin/bash

if [ $# -ne 5 ]
then
    echo "**********************************************************************"
    echo "nRex: Single-nucleotide variant calling."
    echo "This program comes with ABSOLUTELY NO WARRANTY."
    echo ""
    echo "nRex (Version: 0.0.5)"
    echo "Contact: Tobias Rausch (rausch@embl.de)"
    echo "**********************************************************************"
    echo ""
    echo "Usage: $0 <wgs|wex|haloplex> <genome.fa> <output prefix> <sample1.read1.fq.gz> <sample1.read2.fq.gz>"
    echo ""
    exit -1
fi

SCRIPT=$(readlink -f "$0")
BASEDIR=$(dirname "$SCRIPT")

export PERL5LIB=${BASEDIR}/perl/lib/perl5/:${BASEDIR}/perl/lib/5.24.0/
export PATH=/g/funcgen/bin:${BASEDIR}/perl/bin:${PATH}

# CMD parameters
ATYPE=${1}
GENOME=${2}
OUTP=${3}
FQ1=${4}
FQ2=${5}
THREADS=4

# Programs
VEP=${BASEDIR}/vep/vep.pl
VEP_DATA=${BASEDIR}/vepcache

# Redirect tmp to scratch directory if available
if [ -n "${SCRATCHDIR}" ]
then
    export TMP=${SCRATCHDIR}
fi

# Generate IDs
FQ1ID=`echo ${OUTP} | sed 's/$/.1.fq/'`
FQ2ID=`echo ${OUTP} | sed 's/$/.2.fq/'`
BAMID=`echo ${OUTP} | sed "s/$/.align/"`

# Fastqc
mkdir -p ${OUTP}_prefastqc/ && fastqc -t ${THREADS} -o ${OUTP}_prefastqc/ ${FQ1} && fastqc -t ${THREADS} -o ${OUTP}_prefastqc/ ${FQ2}
    
# BWA
bwa mem -t ${THREADS} -R "@RG\tID:${OUTP}\tSM:${OUTP}" ${GENOME} ${FQ1} ${FQ2} | samtools view -bT ${GENOME} - > ${BAMID}.bam

# Sort & Index
samtools sort -@ ${THREADS} -o ${BAMID}.srt.bam ${BAMID}.bam && rm ${BAMID}.bam && samtools index ${BAMID}.srt.bam

# Mark duplicates
if [ ${ATYPE} == "haloplex" ]
then
    mv ${BAMID}.srt.bam ${BAMID}.rmdup.bam
    mv ${BAMID}.srt.bam.bai ${BAMID}.rmdup.bam.bai
else
    bammarkduplicates markthreads=${THREADS} I=${BAMID}.srt.bam O=${BAMID}.rmdup.bam M=${BAMID}.metrics.tsv index=1 rmdup=0
    rm ${BAMID}.srt.bam ${BAMID}.srt.bam.bai
fi
    
# Run stats using unfiltered BAM
samtools idxstats ${BAMID}.rmdup.bam > ${OUTP}.idxstats
samtools flagstat ${BAMID}.rmdup.bam > ${OUTP}.flagstat

# Filter duplicates, unmapped reads, chrM and unplaced contigs
if [ ${ATYPE} == "haloplex" ]
then
    CHRS=`cat ${BASEDIR}/../bed/haloplex.bed | cut -f 1 | sort -k1,1V -k2,2n | uniq | tr '\n' ' '`
else
    CHRS=`cat ${BASEDIR}/../bed/exome.bed | cut -f 1 | sort -k1,1V -k2,2n | uniq | tr '\n' ' '`
fi
samtools view -F 1024 -b ${BAMID}.rmdup.bam ${CHRS} > ${BAMID}.final.bam
samtools index ${BAMID}.final.bam
rm ${BAMID}.rmdup.bam ${BAMID}.rmdup.bam.bai

# Run stats using filtered BAM, TSS enrichment, error rates, etc.
if [ ${ATYPE} == "haloplex" ]
then
    alfred qc -b ${BASEDIR}/../bed/haloplex.bed -r ${GENOME} -o ${OUTP}.bamStats ${BAMID}.final.bam
else
    alfred qc -b ${BASEDIR}/../bed/exome.bed -r ${GENOME} -o ${OUTP}.bamStats ${BAMID}.final.bam
fi

# Freebayes
freebayes --no-partial-observations --min-repeat-entropy 1 --report-genotype-likelihood-max --min-alternate-fraction 0.15 --fasta-reference ${GENOME} --genotype-qualities ${BAMID}.final.bam -v ${BAMID}.vcf
bgzip ${BAMID}.vcf
tabix ${BAMID}.vcf.gz

# Normalize VCF
vt normalize ${BAMID}.vcf.gz -r ${GENOME} | vt decompose_blocksub - | vt decompose - | vt uniq - | bgzip > ${BAMID}.norm.vcf.gz
tabix ${BAMID}.norm.vcf.gz
rm ${BAMID}.vcf.gz ${BAMID}.vcf.gz.tbi

# Fixed threshold filtering
if [ ${ATYPE} == "haloplex" ]
then
    bcftools filter -O z -o ${BAMID}.norm.filtered.vcf.gz -e '%QUAL<=20 || %QUAL/AO<=2 || SAF<=2 || SAR<=2' ${BAMID}.norm.vcf.gz
else
    bcftools filter -O z -o ${BAMID}.norm.filtered.vcf.gz -e '%QUAL<=20 || %QUAL/AO<=2 || SAF<=2 || SAR<=2 || RPR<=2 || RPL<=2' ${BAMID}.norm.vcf.gz
fi
tabix ${BAMID}.norm.filtered.vcf.gz
rm ${BAMID}.norm.vcf.gz ${BAMID}.norm.vcf.gz.tbi
exit;

# VEP
perl ${VEP} --species homo_sapiens --assembly GRCh37 --offline --no_progress --no_stats --sift b --polyphen b --ccds --hgvs --symbol --numbers --regulatory --canonical --protein --biotype --tsl --appris --gene_phenotype --af --af_1kg --af_esp --af_exac --pubmed --variant_class --no_escape --vcf --minimal --dir ${VEP_DATA} --fasta ${VEP_DATA}/Homo_sapiens.GRCh37.75.dna.primary_assembly.fa --input_file ${OUTP}/${BAMID}.norm.filtered.vcf.gz --output_file ${OUTP}/${BAMID}.vep.vcf --plugin ExAC,${VEP_DATA}/ExAC.r0.3.1.sites.vep.vcf.gz --plugin Blosum62 --plugin CSN --plugin Downstream --plugin GO --plugin LoFtool,${BASEDIR}/../bed/LoFtool_scores.txt --plugin TSSDistance --plugin MaxEntScan,${VEP_DATA}/Plugins/maxentscan/
bgzip ${OUTP}/${BAMID}.vep.vcf
tabix ${OUTP}/${BAMID}.vep.vcf.gz

# Convert to BCF
bcftools view -O b -o ${OUTP}/${BAMID}.vep.bcf ${OUTP}/${BAMID}.vep.vcf.gz
bcftools index ${OUTP}/${BAMID}.vep.bcf
rm ${OUTP}/${BAMID}.vep.vcf.gz ${OUTP}/${BAMID}.vep.vcf.gz.tbi ${OUTP}/${BAMID}.norm.filtered.vcf.gz ${OUTP}/${BAMID}.norm.filtered.vcf.gz.tbi

