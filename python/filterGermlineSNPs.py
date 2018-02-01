#! /usr/bin/env python

from __future__ import print_function
import argparse
import sys
import cyvcf2


# Parse command line
parser = argparse.ArgumentParser(description='Filter germline SNPs using matched tumor(s)')
parser.add_argument('-v', '--vcf', metavar='input.vcf.gz', required=True, dest='vcf', help='input VCF file (required)')
parser.add_argument('-g', '--germ', metavar='blood', required=True, dest='germ', help='germline sample (required)')
parser.add_argument('-t', '--tumor', metavar='tumor', required=True, dest='tumor', help='tumor sample(s), separated by comma (required)')
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

        # Simple quality filter
        if record.QUAL <= minqual:
            continue

        # Germline carrier?
        keepSite = False
        for gIdx in germidx:
            if record.gt_types[gIdx] != vcf.HOM_REF:
                keepSite = True
        if not keepSite:
            continue

        # Tumor confirmed?
        tumorConfirmed = False
        ao = record.format('AO')
        for idx in tumidx:
            # Multi-allelics?
            for i in range(len(ao[idx])):
                if ao[idx][i] >= 2:
                    tumorConfirmed = True
        if not tumorConfirmed:
            continue

        # Keep site
        w.write_record(record)

    # Close file
    w.close()
