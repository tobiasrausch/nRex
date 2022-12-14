#! /usr/bin/env python

from __future__ import print_function
import argparse
import sys
import cyvcf2


# Parse command line
parser = argparse.ArgumentParser(description='Filter somatic variants')
parser.add_argument('-v', '--vcf', metavar='input.vcf.gz', required=True, dest='vcf', help='input VCF file (required)')
parser.add_argument('-g', '--germ', metavar='blood', required=True, dest='germ', help='germline sample(s) (required)')
parser.add_argument('-t', '--tumor', metavar='tumor', required=True, dest='tumor', help='tumor sample (required)')
parser.add_argument('-q', '--quality', metavar='20', required=False, dest='quality', help='min. quality (optional)')
parser.add_argument('-c', '--coverage', metavar='10', required=False, dest='coverage', help='min. coverage (optional)')
parser.add_argument('-n', '--ncutoff', metavar='1', required=False, dest='ncutoff', help='remove variant if control has this ALT read support (optional)')
parser.add_argument('-o', '--out', metavar='out.vcf', required=True, dest='out', help='output VCF file (required)')
args = parser.parse_args()


# Min. quality
minqual = 20
if args.quality:
    minqual = int(args.quality)

# Min. coverage
mincov = 10
if args.coverage:
    mincov = int(args.coverage)

# Normal ALT support cutoff
ncutoff = 1
if args.ncutoff:
    ncutoff = int(args.ncutoff)
    
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

        # No missing genotypes
        unknown = False
        for gIdx in germidx:
            if record.gt_types[gIdx] == vcf.UNKNOWN:
                unknown = True
        for tIdx in tumidx:
            if record.gt_types[tIdx] == vcf.UNKNOWN:
                unknown = True
        if unknown:
            continue

        # Tumor carrier?
        ao = record.format('AO')
        ro = record.format('RO')
        keepSite = False
        for tIdx in tumidx:
            if record.gt_types[tIdx] != vcf.HOM_REF:
                for i in range(len(ao[tIdx])):
                    dp = ro[tIdx][0] + ao[tIdx][i]
                    if (dp >= mincov) and (ao[tIdx][i] >= 2):
                        keepSite = True
        if not keepSite:
            continue

        # Absent in germline
        germline = False
        validCov = False
        for gIdx in germidx:
            for i in range(len(ao[gIdx])):
                dp = ro[gIdx][0] + ao[gIdx][i]
                if dp >= mincov:
                    validCov = True
                if ao[gIdx][i] >= ncutoff:
                    germline = True
                    break
        if not validCov:
            germline = True
        if germline:
            continue

        # Keep site
        w.write_record(record)

    # Close file
    w.close()
