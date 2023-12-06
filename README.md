# nRex

nRex is a germline & somatic single-nucleotide variant calling pipeline for whole-genome human genomics sequencing data (>=30x).

## Installing nRex

`git clone --recursive https://github.com/tobiasrausch/nRex.git`

`cd nRex/`

`make all`

## Download reference genome

`cd genome/ && ./prepare_genome.sh`

## Download a haplotype reference panel

`cd refpanel/ && ./download.sh`

`cd maps/ && ./download.sh`

## Optional: VEP annotation of variants

To activate the annotation of variants using VEP, you need to download an offline annotation cache for GRCh38

`make .vep`

## Running nRex for GRCh38

There is a pipeline for short- and long-reads. For short-reads:

`./src/nRex.sh <output prefix> <sample1.read1.fq.gz> <sample1.read2.fq.gz>`

and for long-reads:

`./src/nRexLR.sh <output prefix> <sample1.read1.fq.gz> <sample1.read2.fq.gz>`

## Postprocessing the output of the pipeline

A few helper scripts to summarize the output of the various tools.

### Aggregating QC statistics across multiple samples

`./src/aggregate.sh table *.qc.summary`

### Summarizing VEP output

To generate a table of annotated variants, you can use

`python3 scripts/vep.py -v sample.vep.bcf`
