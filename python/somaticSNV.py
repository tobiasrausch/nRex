#! /usr/bin/env python

from __future__ import print_function
import scipy.stats
import vcf
import argparse
import csv
import sys
import collections
import numpy
import re

# Parse command line
parser=argparse.ArgumentParser(description='Filter for somatic SVs.')
parser.add_argument('-v', '--vcf', metavar='variants.vcf', required=True, dest='vcfFile', help='input vcf file (required)')
parser.add_argument('-g', '--gq', metavar='25', required=False, dest='gqCut', help='GQ threshold (optional)')
parser.add_argument('-o', '--out', metavar='out.vcf', required=True, dest='outFile', help='output vcf file (required)')
parser.add_argument('-t', '--type', metavar='SNV', required=False, dest='variantType', help='variant type [SNV, INDEL] (optional)')
args = parser.parse_args()


# Command-line args
gqCut=25
qualCut=20
if args.gqCut:
    gqCut=int(args.gqCut)
variantType=0
if (args.variantType) and (args.variantType=="INDEL"):
    variantType=1

# Parameters
pvCut=0.01
bafCut=0.04

# Output vcf records
if args.vcfFile:
    vcf_reader=vcf.Reader(open(args.vcfFile), 'r')
    vcf_reader.infos['SOMATIC'] = vcf.parser._Info('SOMATIC', 0, 'Flag', 'Somatic SNV.')
    vcf_writer = vcf.Writer(open(args.outFile, 'w'), vcf_reader, lineterminator='\n')
    for record in vcf_reader:
        if ((variantType==0) and (record.is_snp)) or ((variantType==1) and (record.is_indel)):
            gqRef=[]
            gqAlt=[]
            for call in record.samples:
                if (((re.search(r"[Nn]ormal", call.sample)!=None) or (re.search(r"[Cc]ontrol", call.sample)!=None)) and (call.called) and (call.gt_type==0)):
                #if ((re.search(r"NDNA", call.sample)!=None) and (call.called) and (call.gt_type==0)):
                    dvVal=getattr(call.data, 'DV', None)
                    vafPass=False
                    dpVal=0
                    if (dvVal==None):
                        adVal=getattr(call.data, 'AD', None)
                        if (adVal==None):
                            # Freebayes
                            roVal=getattr(call.data, 'RO', None)
                            aoVal=getattr(call.data, 'AO', None)
                            if (roVal is not None) and (aoVal is not None) and (len(record.INFO['TYPE'])==1):
                                if ((variantType==1)  or ((variantType==0) and (record.INFO['TYPE'][0]=="snp"))):
                                    dpVal=float(roVal+aoVal)
                                    if (dpVal>0):
                                        baf=float(aoVal)/dpVal
                                        if (baf<=bafCut):
                                            vafPass=True
                        else:
                            #GATK
                            dpVal=float(adVal[0]+adVal[1])
                            dvVal=float(adVal[1])
                            if (dpVal>0):
                                baf=dvVal/dpVal
                                if (baf<=2*bafCut):  # GATK seems to be counting differently that's why 2*
                                    vafPass=True
                    else:
                        #Mpileup
                        dpVal=float(call['DP'])
                        dvVal=float(call['DV'])
                        if (dpVal>0):
                            baf=dvVal/dpVal
                            if (baf<=bafCut):
                                vafPass=True
                    if (vafPass):
                        gqRef.append(call['GQ'])
                if ((re.search(r"[Tt]umo", call.sample)!=None) and (call.called) and (call.gt_type!=0)):
                #if ((re.search(r"CLL", call.sample)!=None) and (call.called) and (call.gt_type!=0)):
                    gqAlt.append(call['GQ'])
            # Check p-values if present
            pv4Pass=True
            if ('PV4' in record.INFO):
                lowP=0
                for pval in record.INFO['PV4']:
                    if (float(pval)<pvCut):
                        lowP+=1
                if (lowP>1) or (float(record.INFO['PV4'][0])<pvCut/2) or (float(record.INFO['PV4'][2])<pvCut/2):
                    pv4Pass=False
            if ((len(gqRef)>0) and (len(gqAlt)>0) and (numpy.median(gqRef)>=gqCut) and (numpy.median(gqAlt)>=gqCut) and (record.QUAL>=qualCut) and (pv4Pass)):
                record.INFO['SOMATIC']=True;
                vcf_writer.write_record(record)
