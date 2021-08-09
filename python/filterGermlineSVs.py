#! /usr/bin/env python

from __future__ import print_function
import argparse
import sys
import cyvcf2


# Parse command line
parser = argparse.ArgumentParser(description='Filter germline SNPs using matched tumor(s)')
parser.add_argument('-v', '--vcf', metavar='input.vcf.gz', required=True, dest='vcf', help='input VCF file (required)')
parser.add_argument('-g', '--germ', metavar='blood', required=True, dest='germ', help='germline sample (required)')
parser.add_argument('-t', '--tumor', metavar='tumor', required=False, dest='tumor', help='tumor sample(s), separated by comma (required)')
parser.add_argument('-q', '--quality', metavar='20', required=False, dest='quality', help='min. quality (optional)')
parser.add_argument('-o', '--out', metavar='out.vcf', required=True, dest='out', help='output VCF file (required)')
args = parser.parse_args()

# Min. quality
minqual = 20
if args.quality:
    minqual = int(args.quality)

# Parse VCF
vcf = cyvcf2.VCF(args.vcf)
w = cyvcf2.Writer(args.out, vcf)

samples = list(vcf.samples)
tumidx = set()
if args.tumor:
    for t in args.tumor.split(','):
        if t not in samples:
            print(t, "does not exist in VCF!", sep=" ", file=sys.stderr)
            quit()
        tumidx.add(samples.index(t))
germidx = set()
for g in args.germ.split(','):
    if g not in samples:
        print(g, "does not exist in VCF!", sep=" ", file=sys.stderr)
        quit()
    germidx.add(samples.index(g))
if len(germidx) and len(tumidx):
    for record in vcf:
        # Ignore multi-allelics
        if len(record.ALT) > 1:
            continue

        # Any filter set?
        if record.FILTER is not None:
            continue
        
        # Ignore imprecise
        if record.INFO.get('PRECISE') is None:
            continue

        # Simple quality filter
        if record.QUAL <= minqual:
            continue

        # SR mapping quality filter
        if record.INFO.get("SRMAPQ") <= minqual:
            continue

        # Good consensus alignment
        if record.INFO.get("SRQ") <= 0.8:
            continue

        # At least 3 split reads
        if record.INFO.get("SR") < 3:
            continue

        # Germline carrier?
        keepSite = False
        for gIdx in germidx:
            if (record.gt_types[gIdx] != vcf.HOM_REF) and (record.gt_types[gIdx] != vcf.UNKNOWN):
                keepSite = True
        if not keepSite:
            continue

        # Tumor confirmed?
        tumorConfirmed = False
        rv = record.format('RV')
        for idx in tumidx:
            if rv[idx][0] >= 1:
                tumorConfirmed = True
        if not tumorConfirmed:
            continue

        # Keep site
        w.write_record(record)
# Close file
w.close()
