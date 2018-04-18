#!/bin/bash

# ExAC mask
wget ftp://ftp.broadinstitute.org/pub/ExAC_release/release0.3.1/resources/exome_calling_regions.v1.interval_list
cat exome_calling_regions.v1.interval_list | grep -v "^@" | grep -v "^MT" | sed 's/^/chr/' | gzip -c > exac.hg19.bed.gz
rm exome_calling_regions.v1.interval_list

# SGDP regions to be filtered
wget https://github.com/lh3/sgdp-fermi/releases/download/v1/sgdp-263-hs37d5.tgz
tar -xzf sgdp-263-hs37d5.tgz
zcat sgdp-263-hs37d5/um75-hs37d5.bed.gz | grep -v "_" | grep -v "^GL" | grep -v "^MT" | grep -v "hs37d5" | sed 's/^/chr/' | gzip -c > sgdp.hg19.bed.gz
rm sgdp-263-hs37d5.tgz
rm -rf sgdp-263-hs37d5/
