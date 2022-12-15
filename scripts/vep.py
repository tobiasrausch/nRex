#! /usr/bin/env python

from __future__ import print_function
import csv
import argparse
import sys
import collections
import json
import cyvcf2

# Parse command line
parser = argparse.ArgumentParser(description='Create variant table')
parser.add_argument('-v', '--variants', metavar='sample.vcf.gz', required=True, dest='variants', help='VCF file (required)')
args = parser.parse_args()

# Desired VEP columns
descols = ["Consequence", "IMPACT", "SYMBOL", "Gene", "STRAND", "Feature", "EXON", "HGVSc", "HGVSp", "cDNA_position", "CDS_position", "Protein_position", "Amino_acids", "Codons", "Existing_variation", "SIFT", "PolyPhen", "gnomADe_AF", "gnomADg_AF", "MAX_AF", "CLIN_SIG"]
selected = ['missense_variant', 'splice_acceptor_variant', 'splice_donor_variant', 'stop_gained', 'stop_lost', 'start_lost', 'frameshift_variant', 'inframe_insertion', 'inframe_deletion']

# VCF/BCF parsing
varstore = collections.defaultdict(collections.defaultdict)
vcf = cyvcf2.VCF(args.variants)
for header in vcf.header_iter():
    if header['HeaderType'] == "INFO":
        if header['ID'] == "CSQ":
            vepcols = header['Description'].replace('Consequence annotations from Ensembl VEP. Format: "', '')[:-1].split('|')
            
addr = dict()
for idx, val in enumerate(vepcols):
    addr[val] = idx

print("chromosome", "position", "REF", "ALT", "carrier", "AF", '\t'.join(descols), sep='\t')
for record in vcf:
    csq = record.INFO.get('CSQ')
    if csq is None:
        print(record.CHROM, record.POS, record.REF, ','.join(record.ALT), "unknown_variant", sep='\t')
        continue
    transcripts = csq.split(',')
    for tr in transcripts:
        fields = tr.split('|')
        consout = False
        for constype in fields[addr['Consequence']].split('&'):
            if constype in selected:
                consout = True
                break
        af_vcf = record.INFO.get('AF')
        if (consout) and (fields[addr['CANONICAL']] == "YES") and (fields[addr['BIOTYPE']] == "protein_coding"):
            if len(fields) != len(vepcols):
                print("Error: VEP annotation is corrupted!", file=sys.stderr)
                sys.exit(1)
            # Debug fields
            #for (a, b) in zip(vepcols, fields):
            #    print(a, b)
            carrierStr=""
            for spl, gt in zip(vcf.samples, record.gt_types):
                # gt_types is array of 0,1,2,3==HOM_REF, HET, UNKNOWN, HOM_ALT
                if (gt != 0) and (gt != 2):
                    carrierStr += spl + "(" + str(gt) + "),"
            fields[addr['EXON']] = fields[addr['EXON']].replace('/', ';')
            print(record.CHROM, record.POS, record.REF, ','.join(record.ALT), carrierStr, af_vcf, '\t'.join([fields[addr[cname]] for cname in descols]), sep='\t')
