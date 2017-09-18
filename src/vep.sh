#!/bin/bash

if [ $# -ne 1 ]
then
    echo "**********************************************************************"
    echo "VEP annotation"
    echo "This program comes with ABSOLUTELY NO WARRANTY."
    echo ""
    echo "nRex (Version: 0.0.5)"
    echo "Contact: Tobias Rausch (rausch@embl.de)"
    echo "**********************************************************************"
    echo ""
    echo "Usage: $0 <input.vcf.gz>"
    echo ""
    exit -1
fi

SCRIPT=$(readlink -f "$0")
BASEDIR=$(dirname "$SCRIPT")

# Load VEP dependencies
module load Perl BioPerl DBD-mysql HTSlib
export PERL5LIB=${BASEDIR}/BioDbBuild/lib/perl5/x86_64-linux-thread-multi:${PERL5LIB}

# CMD parameters
VCF=${1}
ID=`echo ${VCF} | sed 's/.vcf.gz$//'`
FORK=8

# Programs
VEP=${BASEDIR}/vep/vep
VEP_DATA=${BASEDIR}/vepcache

# VEP
${VEP} --fork ${FORK} --species homo_sapiens --assembly GRCh37 --offline --everything --stats_file ${ID}.vep.html --format vcf --cache --dir ${VEP_DATA} --fasta ${VEP_DATA}/Homo_sapiens.GRCh37.75.dna.primary_assembly.fa --input_file ${VCF} --output_file ${ID}.vep.vcf --no_escape --custom ${VEP_DATA}/gnomad.genomes.r2.0.1.sites.noVEP.vcf.gz,gnomADg,vcf,exact,0,AF_AFR,AF_AMR,AF_ASJ,AF_EAS,AF_FIN,AF_NFE,AF_OTH --plugin Blosum62 --plugin CSN --plugin Downstream --plugin GO --plugin LoFtool,${VEP_DATA}/Plugins/LoFtool_scores.txt --plugin TSSDistance --plugin MaxEntScan,${VEP_DATA}/Plugins/maxentscan/
bgzip ${ID}.vep.vcf
tabix ${ID}.vep.vcf.gz

# Convert to BCF
bcftools view -O b -o ${ID}.vep.bcf ${ID}.vep.vcf.gz
bcftools index ${ID}.vep.bcf
rm ${ID}.vep.vcf.gz ${ID}.vep.vcf.gz.tbi

