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
parser.add_argument('-f', '--vaf', metavar='0.15', required=False, dest='vaf', help='min. VAF')
parser.add_argument('-q', '--qual', metavar='20', required=False, dest='qual', help='min. quality')
parser.add_argument('-r', '--report', metavar='0.05', required=False, dest='report', help='min. reporting VAF')
parser.add_argument('-a', '--ao', metavar='2', required=False, dest='minao', help='min. alternative observation')
parser.add_argument('-n', '--no-consequence-selection', dest='nocons', action='store_false')
parser.set_defaults(nocons=True)
args = parser.parse_args()

minao = 2
if args.minao:
    minao = (int) (args.minao)
vafth = 0.15
if args.vaf:
    vafth = (float) (args.vaf)
qual = 20
if args.qual:
    qual = (int) (args.qual)
reportth = 0.05
if args.report:
    reportth = (float) (args.report)
if vafth < reportth:
    reportth = vafth
    
# Desired VEP columns
descols = ["Consequence", "IMPACT", "SYMBOL", "Gene", "STRAND", "Feature", "EXON", "HGVSc", "HGVSp", "cDNA_position", "CDS_position", "Protein_position", "Amino_acids", "Codons", "Existing_variation", "SIFT", "PolyPhen", "gnomADe_AF", "gnomADg_AF", "MAX_AF", "CLIN_SIG"]
selected = ['missense_variant', 'splice_acceptor_variant', 'splice_donor_variant', 'stop_gained', 'stop_lost', 'start_lost', 'frameshift_variant', 'inframe_insertion', 'inframe_deletion']

# VCF/BCF parsing
varstore = collections.defaultdict(collections.defaultdict)
vcf = cyvcf2.VCF(args.variants)
amPresent = False
for header in vcf.header_iter():
    if header['HeaderType'] == "INFO":
        if header['ID'] == "CSQ":
            vepcols = header['Description'].replace('Consequence annotations from Ensembl VEP. Format: "', '')[:-1].split('|')
        if header['ID'] == "AMCLASS":
            amPresent = True
            
addr = dict()
for idx, val in enumerate(vepcols):
    addr[val] = idx

# Header
if amPresent:
    print("chromosome", "position", "REF", "ALT", "qual", "carrier", "ampath", "amclass", '\t'.join(descols), sep='\t')
else:
    print("chromosome", "position", "REF", "ALT", "qual", "carrier", '\t'.join(descols), sep='\t')

# Parse VCF
for record in vcf:
    if (record.FILTER is None) or (record.FILTER == ".") or (record.FILTER == "PASS"):
        pass
    else:
        continue
    if record.QUAL < qual:
        continue
    csq = record.INFO.get('CSQ')
    if amPresent:
        amclass = record.INFO.get('AMCLASS')
        ampath = record.INFO.get('AMPATH')
    if csq is None:
        print(record.CHROM, record.POS, record.REF, ','.join(record.ALT), "unknown_variant", sep='\t')
        continue
    transcripts = csq.split(',')
    for tr in transcripts:
        fields = tr.split('|')
        consout = False
        pathogenic = False
        if fields[addr['CLIN_SIG']] == "pathogenic":
            pathogenic = True
        for constype in fields[addr['Consequence']].split('&'):
            if (not args.nocons) or (constype in selected) or (pathogenic):
                consout = True
                break
        if (consout) and (fields[addr['CANONICAL']] == "YES") and (fields[addr['BIOTYPE']] == "protein_coding"):
            if len(fields) != len(vepcols):
                print("Error: VEP annotation is corrupted!", file=sys.stderr)
                sys.exit(1)
            # Debug fields
            #for (a, b) in zip(vepcols, fields):
            #    print(a, b)
            carrierStr=""
            validSite=False
            aoPresent=True
            try:
                ao = record.format('AO')
            except KeyError:
                aoPresent=False
            if aoPresent:
                for spl, gt, ao, ro in zip(vcf.samples, record.gt_types, record.format('AO'), record.format('RO')):
                    # gt_types is array of 0,1,2,3==HOM_REF, HET, UNKNOWN, HOM_ALT
                    if (gt != 2):
                        if (ao[0] + ro[0]) > 0:
                            af = (float) (ao[0]) / (float) (ao[0] + ro[0])
                            if (ao[0] >= minao) and (af > vafth):
                                validSite = True
                            if (ao[0] >= minao) and (af > reportth):
                                carrierStr += spl + "(VAF=" + str(round(af,2)) + "),"
            else:
                for spl, gt, dp, ad in zip(vcf.samples, record.gt_types, record.format('DP'), record.format('AD')):
                    # gt_types is array of 0,1,2,3==HOM_REF, HET, UNKNOWN, HOM_ALT
                    if (gt != 2):
                        if dp > 0:
                            af = (float) (ad) / (float) (dp)
                            if (ad >= minao) and (af > vafth):
                                validSite = True
                            if (ad >= minao) and (af > reportth):
                                carrierStr += spl + "(VAF=" + str(round(af,2)) + "),"
            fields[addr['EXON']] = fields[addr['EXON']].replace('/', ';')
            if validSite:
                if amPresent:
                    print(record.CHROM, record.POS, record.REF, ','.join(record.ALT), record.QUAL, carrierStr, ampath, amclass, '\t'.join([fields[addr[cname]] for cname in descols]), sep='\t')
                else:
                    print(record.CHROM, record.POS, record.REF, ','.join(record.ALT), record.QUAL, carrierStr, '\t'.join([fields[addr[cname]] for cname in descols]), sep='\t')
