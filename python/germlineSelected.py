#! /usr/bin/env python

from __future__ import print_function
import argparse
import sys
import cyvcf2


# Parse command line
parser = argparse.ArgumentParser(description='Filter germline SNPs using matched tumor(s)')
parser.add_argument('-v', '--vcf', metavar='input.vcf.gz', required=True, dest='vcf', help='input VCF file (required)')
args = parser.parse_args()

# Parse VCF
vcf = cyvcf2.VCF(args.vcf)
samples = list(vcf.samples)
for record in vcf:
    # Ignore multi-allelics
    if len(record.ALT) > 1:
        continue

    # Compute AF
    ao = record.format('AO')
    ro = record.format('RO')
    af = dict()
    for idx, s in enumerate(samples):
        if (record.gt_types[idx] != vcf.UNKNOWN) and (record.gt_types[idx] != vcf.HOM_REF):
            if ro[idx][0] + ao[idx][0] > 0:
                afs = float(ao[idx][0]) / float(ro[idx][0] + ao[idx][0])
                af[s] = afs

    # Germline carrier?
    maxGerm = 0
    for k in af.keys():
        if k.find("REM/") != -1:
            if af[k] > maxGerm:
                maxGerm = af[k]
    maxTum = 0
    for k in af.keys():
        if k.find("REM/") == -1:
            if af[k] > maxTum:
                maxTum = af[k]
    if maxGerm > 0:
        if (maxGerm * 1.2 < maxTum):
            idstr = "."
            if record.ID is not None:
                idstr = record.ID
            print(record.CHROM, record.POS, idstr, record.REF, ','.join(record.ALT), end="", sep="\t")
            for k in af.keys():
                print("\t" + k + "(" + str(af[k]) + ")", end="")
            print()




