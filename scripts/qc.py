#! /usr/bin/env python

import csv
import argparse
import collections
import os
import json
import gzip

# Parse command line
parser = argparse.ArgumentParser(description='Aggregate QC statistics')
parser.add_argument('-p', '--prefix', required=True, dest='prefix', help='file prefix (required)')
args = parser.parse_args()

qc = dict()
qc['Sample'] = args.prefix

# Trim galore
for rn in ["1", "2"]:
    filep = args.prefix + "." + rn + ".fq.gz_trimming_report.txt"
    if (os.path.exists(filep)) and (os.path.isfile(filep)):
        with open(filep) as f:
            for line in f:
                if line.startswith("Reads with adapters"):
                    qc['AdaptersRead' + rn] = line[(line.find('(')+1):line.find(')')]

# Parse alignment statistics
filep = args.prefix + ".alfred.tsv.gz"
if (os.path.exists(filep)) and (os.path.isfile(filep)):
    with gzip.open(filep, "rb") as f:
        columns = None
        for line in f.read().decode("ascii").splitlines():
            if line.startswith("ME"):
                if columns is None:
                    columns = line.strip().split('\t')
                else:
                    records = line.strip().split('\t')
                    qc['PercUnmapped'] = str(round(float(records[columns.index('UnmappedFraction')]) * 100, 2)) + "%"
                    qc['PercSeqError'] = str(round(float(records[columns.index('ErrorRate')]) * 100, 2)) + "%"
                    qc['MedianCoverage'] = records[columns.index('MedianCoverage')]
                    qc['SDCoverage'] = records[columns.index('SDCoverage')]
                    qc['MedianInsertSize'] = records[columns.index('MedianInsertSize')]
                    qc['PercDuplicate'] = str(round(float(records[columns.index('DuplicateFraction')]) * 100, 2)) + "%"
                    qc['PercMappedSameChr'] = str(round(float(records[columns.index('MappedSameChrFraction')]) * 100, 2)) + "%"

# Variants
filep = args.prefix + ".vcf.gz"
if (os.path.exists(filep)) and (os.path.isfile(filep)):
    lnum = 0
    with gzip.open(filep, 'rt') as f:
        for line in f:
            if not line.startswith("#"):
                lnum += 1
    qc['UnfilteredVariants'] = lnum
filep = args.prefix + ".hg38.wgs.vcf.gz"
if (os.path.exists(filep)) and (os.path.isfile(filep)):
    lnum = 0
    with gzip.open(filep, 'rt') as f:
        for line in f:
            if not line.startswith("#"):
                lnum += 1
    qc['FilteredVariants'] = lnum
    
                    
# Output QC dictionary
#print(qc.keys())
for key in ['Sample', 'AdaptersRead1', 'AdaptersRead2', 'PercUnmapped', 'PercSeqError', 'MedianCoverage', 'SDCoverage', 'MedianInsertSize', 'PercDuplicate', 'PercMappedSameChr', 'UnfilteredVariants', 'FilteredVariants']:
    if key in qc.keys():
        print(key, qc[key])
