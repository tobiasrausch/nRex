# nRex

nRex is a germline & somatic single-nucleotide variant calling pipeline for whole-genome human genomics sequencing data (>=30x).

## Installing nRex

`git clone --recursive https://github.com/tobiasrausch/nRex.git`

`cd nRex/`

`make all`

## Download reference genome

`cd genome/ && ./prepare_genome.sh`

## Optional: VEP annotation of variants

To activate the annotation of variants using VEP, you need to download an offline annotation cache for GRCh38

`make .vep`

## Optional: Phasing of variants using a haplotype reference panel

To activate the phasing of variants, you need to download a haplotype reference panel

`cd refpanel/ && ./download.sh`

and install shapeit (see shapeit documentation)

`make .shapeit`

## Running nRex for GRCh38

`./src/nRex.sh <output prefix> <sample1.read1.fq.gz> <sample1.read2.fq.gz>`

## Postprocessing the output of the pipeline

### Aggregating QC statistics across multiple samples

`./src/aggregate.sh table *.qc.summary`

### Summarizing VEP output

To generate a table of annotated variants, you can use

`python3 scripts/vep.py -v sample.vep.bcf`
