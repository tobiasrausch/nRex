# nRex

nRex is a germline & somatic single-nucleotide variant calling pipeline

## Installing nRex

`git clone --recursive https://github.com/tobiasrausch/nRex.git`

`cd nRex/`

`make all`

## Download reference genome

`cd genome/ && ./prepare_genome.sh`

## Download haplotype reference panel

`cd refpanel/ && ./download.sh`

## Running nRex for GRCh38

`./src/nRex.sh <output prefix> <sample1.read1.fq.gz> <sample1.read2.fq.gz>`
