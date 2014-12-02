nRex
====

nRex is a somatic single-nucleotide variant calling meta-method. It uses mpileup, freebayes and optionally GATK. If you can improve the current workflow please let me know.

Installing nRex
---------------

To build nRex you need to do the following:

`git clone --recursive https://github.com/tobiasrausch/nRex.git`

`cd nRex/`

`make all`


Running nRex
------------

nRex needs the reference genome and one bam file each for the tumor and the normal sample.

`./nRex.sh <ref.fa> <tumor.bam> <normal.bam> ...`

There is also a dockerized nRex version available [here](https://registry.hub.docker.com/u/trausch/nRex/).
