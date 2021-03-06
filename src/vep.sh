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
module load HTSlib/1.5-foss-2016b BCFtools/1.5-foss-2016b Perl/5.24.1-foss-2016b BioPerl/1.7.1-foss-2016b-Perl-5.24.1 DBD-mysql/4.042-foss-2016b-Perl-5.24.1
export PERL5LIB=${BASEDIR}/BioDbBuild/lib/perl5/x86_64-linux-thread-multi:${PERL5LIB}

# CMD parameters
VCF=${1}.input.vcf
ID=`echo ${VCF} | sed 's/.vcf.gz$//'`
THREADS=4

# Programs
VEP=${BASEDIR}/vep/vep
VEP_DATA=${BASEDIR}/vepcache

# Rename chrs
rm -f rn.chrs
for CHR in 1 2 3 4 5 6 7 8 9 10 11 12 13 14 15 16 17 18 19 20 21 22 X Y
do
    echo chr${CHR} ${CHR} >> rn.chrs
done

# Decompose multi-allelics and normalize InDels
bcftools annotate --rename-chrs rn.chrs ${1} | bcftools norm --threads ${THREADS} -f ${VEP_DATA}/Homo_sapiens.GRCh37.75.dna.primary_assembly.fa -m -both - > ${VCF}
# Without chr renaming
# bcftools norm --threads ${THREADS} -f ${VEP_DATA}/Homo_sapiens.GRCh37.75.dna.primary_assembly.fa -m -both ${1} > ${VCF}
rm rn.chrs

# VEP
${VEP} --fork ${THREADS} --species homo_sapiens --assembly GRCh37 --offline --everything --no_stats --format vcf --cache --dir ${VEP_DATA} --fasta ${VEP_DATA}/Homo_sapiens.GRCh37.75.dna.primary_assembly.fa --input_file ${VCF} --vcf --output_file ${ID}.vep.vcf --no_escape --custom ${VEP_DATA}/gnomad.genomes.r2.0.1.sites.noVEP.vcf.gz,gnomADg,vcf,exact,0,AF_AFR,AF_AMR,AF_ASJ,AF_EAS,AF_FIN,AF_NFE,AF_OTH --plugin Blosum62 --plugin CSN --plugin Downstream --plugin GO --plugin LoFtool,${VEP_DATA}/Plugins/LoFtool_scores.txt --plugin TSSDistance --plugin MaxEntScan,${VEP_DATA}/Plugins/maxentscan/ --plugin MPC,${VEP_DATA}/Plugins/fordist_constraint_official_mpc_values_v2.txt.gz --plugin ExACpLI,${VEP_DATA}/Plugins/ExACpLI_values.txt
bgzip ${ID}.vep.vcf
tabix ${ID}.vep.vcf.gz
rm ${VCF}

# Convert to BCF
bcftools view -O b -o ${ID}.vep.bcf ${ID}.vep.vcf.gz
bcftools index ${ID}.vep.bcf
rm ${ID}.vep.vcf.gz ${ID}.vep.vcf.gz.tbi

