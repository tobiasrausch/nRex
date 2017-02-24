nRex
====

nRex is a germline & somatic single-nucleotide variant calling method. It uses FreeBayes and vt. If you can improve the current workflow please let me know.

Installing nRex
---------------

To build nRex you need to do the following:

`git clone https://github.com/tobiasrausch/nRex.git`

`cd nRex/`

`make all`


Running nRex
------------

nRex needs the reference genome and one bam file each for the tumor and the normal sample.

`./src/nRex.sh <wgs|wex|haloplex> <genome.fa> <output prefix> <sample1.read1.fq.gz> <sample1.read2.fq.gz>`

