#!/bin/bash

#wget -O- ftp://ftp.ncbi.nlm.nih.gov/genomes/all/GCA/000/001/405/GCA_000001405.15_GRCh38/seqs_for_alignment_pipelines.ucsc_ids/GCA_000001405.15_GRCh38_no_alt_analysis_set.fna.gz | gzip -d > GCA_000001405.15_GRCh38_no_alt_analysis_set.fna  
#samtools faidx GCA_000001405.15_GRCh38_no_alt_analysis_set.fna  
#wget http://ftp.1000genomes.ebi.ac.uk/vol1/ftp/release/20130502/supporting/GRCh38_positions/ALL.chr{{1..22},X,Y}_GRCh38.genotypes.20170504.vcf.gz{,.tbi}
 
#for chr in {1..22} X Y; do  
#    (bcftools view --no-version -h ALL.chr${chr}_GRCh38.genotypes.20170504.vcf.gz | grep -v "^##contig=<ID=[GNh]" | sed 's/^##contig=<ID=MT/##contig=<ID=chrM/;s/^##contig=<ID=\([0-9XY]\)/##contig=<ID=chr\1/'; bcftools view --no-version -H -c 2 ALL.chr${chr}_GRCh38.genotypes.20170504.vcf.gz | sed 's/^/chr/') | bcftools norm --no-version -Ou -m -any | bcftools norm --no-version -Ob -o ALL.chr${chr}_GRCh38.genotypes.20170504.bcf -d none -f GCA_000001405.15_GRCh38_no_alt_analysis_set.fna && bcftools index -f ALL.chr${chr}_GRCh38.genotypes.20170504.bcf  
#done

# Rename
for F in *.bcf
do
    ID=`echo ${F} | sed 's/ALL.//' | sed 's/_.*$/.bcf/'`
    echo ${ID}
    #mv ${F} ${ID}
    #mv ${F}.csi ${ID}.csi
done

# Concatenate without chrY
bcftools concat chr1.bcf chr2.bcf chr3.bcf chr4.bcf chr5.bcf chr6.bcf chr7.bcf chr8.bcf chr9.bcf chr10.bcf chr11.bcf chr12.bcf chr13.bcf chr14.bcf chr15.bcf chr16.bcf chr17.bcf chr18.bcf chr19.bcf chr20.bcf chr21.bcf chr22.bcf chrX.bcf | cut -f 1-8 | bcftools view -O b -o sites.bcf -
bcftools index sites.bcf
