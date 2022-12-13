#!/bin/bash

for CHR in `seq 1 22` X
do
    wget ftp://ftp.1000genomes.ebi.ac.uk/vol1/ftp/data_collections/1000G_2504_high_coverage/working/20220422_3202_phased_SNV_INDEL_SV/1kGP_high_coverage_Illumina.chr${CHR}.filtered.SNV_INDEL_SV_phased_panel.vcf.gz
    wget ftp://ftp.1000genomes.ebi.ac.uk/vol1/ftp/data_collections/1000G_2504_high_coverage/working/20220422_3202_phased_SNV_INDEL_SV/1kGP_high_coverage_Illumina.chr${CHR}.filtered.SNV_INDEL_SV_phased_panel.vcf.gz.tbi
    FILE=`ls *chr${CHR}.*.vcf.gz`
    if [ -f ${FILE} ]
    then
	bcftools view -m2 -M2 -v snps,indels --no-version -O b -o chr${CHR}.bcf ${FILE}
	bcftools index chr${CHR}.bcf
    fi
done
