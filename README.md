nRex
====

nRex is a germline & somatic single-nucleotide variant calling method.

Installing nRex
---------------

`git clone --recursive https://github.com/tobiasrausch/nRex.git`

`cd nRex/`

`make all`


Build BED file with chromosomes and exonic regions for QC
---------------------------------------------------------

`cd R/ && Rscript exon.R`


Running nRex
------------

`./src/nRex.sh <hg19.wgs|hg19.wes|hg19.haloplex> <genome.fa> <output prefix> <sample1.read1.fq.gz> <sample1.read2.fq.gz>`

