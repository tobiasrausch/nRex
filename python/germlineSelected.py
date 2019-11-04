#! /usr/bin/env python

from __future__ import print_function
import argparse
import sys
import cyvcf2


# Parse command line
parser = argparse.ArgumentParser(description='Filter germline SNPs that are selected for in the matched tumor(s)')
parser.add_argument('-v', '--vcf', metavar='input.vcf.gz', required=True, dest='vcf', help='input VCF file (required)')
parser.add_argument('-g', '--germ', metavar='blood', required=True, dest='germ', help='germline sample (required)')
parser.add_argument('-t', '--tumor', metavar='tumor', required=True, dest='tumor', help='tumor sample(s), separated by comma (required)')
parser.add_argument('-q', '--quality', metavar='20', required=False, dest='quality', help='min. quality (optional)')
args = parser.parse_args()

# Min. quality
minqual = 20
if args.quality:
    minqual = int(args.quality)

# Parse VCF
vcf = cyvcf2.VCF(args.vcf)
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

        # Compute AF
        ao = record.format('AO')
        ro = record.format('RO')
        af = dict()
        maxGerm = 0
        maxTum = 0
        for idx, s in enumerate(samples):
            if (record.gt_types[idx] != vcf.UNKNOWN) and (record.gt_types[idx] != vcf.HOM_REF):
                if ro[idx][0] + ao[idx][0] > 0:
                    afs = float(ao[idx][0]) / float(ro[idx][0] + ao[idx][0])
                    af[s] = afs
                    if (idx in germidx) and (afs > maxGerm):
                        maxGerm = afs
                    elif (idx in tumidx) and (afs > maxTum):
                        maxTum = afs
        # Germline variant?
        if maxGerm > 0:
            if (maxGerm * 1.2 < maxTum):
                idstr = "."
                if record.ID is not None:
                    idstr = record.ID
                print(record.CHROM, record.POS, idstr, record.REF, ','.join(record.ALT), end="", sep="\t")
                for k in af.keys():
                    print("\t" + k + "(" + str(af[k]) + ")", end="")
                print()
