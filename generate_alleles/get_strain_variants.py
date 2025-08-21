#!/usr/bin/env python3.6
"""
Avery Davis Bell updates to Corinne Simonti original script
This script generates strain_variants.txt (and others?)
"""

# Initially written by Corinne Simonti
## imports
import argparse
import os
import sys
import matplotlib as mpl
import matplotlib.pyplot as plt
import matplotlib.ticker as ticker
from itertools import cycle
import pandas as pd
import numpy as np
import os
import math
from matplotlib.backends.backend_pdf import PdfPages
import scipy
from scipy import stats
import subprocess
import gzip
import random
import datetime
import time
import getopt
import re
import vcf # update - no 'vcf' package but syntax looks like 'pyvcf' as far as I can tell


def main():
    #### Parse arguments
    argp = argparse.ArgumentParser(prog = "get_strain_variants.py",
                                   description = "Avery Davis Bell updates to Corinne Simonti original script. \
                                       This script generates strain_variants.txt (and others?)")
    argp.add_argument("-trnas", dest = "strConfOut", required = True, type = str,
                      metavar = "tRNAs-confidence-set.out",
                      help = "Path to tRNAscan-SE .out file for high confidence set, e.g. ce11-tRNAs-confidence-set.out")
    argp.add_argument("-vcf", dest = "strVCF", required = True, type = str,
                      metavar = "variants.vcf",
                      help = "Path to VCF file for species whos tRNAs are provided")
    argp.add_argument("-out", dest = "strOut", required = True, type = str,
                      metavar = "outfilestem",
                      help = "Prefix for output files")
    args = argp.parse_args()

    # Assign args to names Corinne uses (instead of updating whole script)
    CHAN_OUT_FILE = args.strConfOut
    VCF_FILE = args.strVCF

    ##### Corinne's whole script (wasn't in functions)
    STATEMENT = "Starting at %s..." % time.ctime()
    print(STATEMENT)

    ## Dictionaries and lists.
    LOG_LIST = []
    write_list = []

    CEHREP_LIST = []
    CEHSKIP_LIST = []

    codon2aa = {}

    chr2start = {}
    chr2end = {}

    ce2seq = {}
    ce2loc = {}    # chr, start, end
    loc2ce = {}
    ce2info = {}    # AA, codon
    ce2strand = {}

    loc2allele = {}
    loc2hr = {}
    loc2ha = {}
    loc2het = {}

    nuc2match = {}
    nuc2match["A"] = "T"
    nuc2match["C"] = "G"
    nuc2match["G"] = "C"
    nuc2match["T"] = "A"

    LOG_LIST.append(STATEMENT + '\n')

    STATEMENT = "Empty dictionaries and lists generated."
    print(STATEMENT)
    LOG_LIST.append(STATEMENT + '\n')
    ## Get tRNA info.
    for line in open(CHAN_OUT_FILE):
        if line[0] == "S" or line[0] == "N" or line[0] == "-":
            continue
        seg = line.split("\t")
        CHR = str(seg[0]).strip()
        NUM = str(seg[1])
        TRNA = CHR + ".trna" + NUM
        START = int(seg[2])
        END = int(seg[3])
        if START > END:
            START = int(seg[3])
            END = int(seg[2])
        AA = str(seg[4])
        # if AA == "Undet":
        #    continue
        CODON = str(seg[5])
        # if CODON == "NNN":
        #    continue
        # NOTE = str(seg[14]) # assigned but never used & was causing issues; note field moves around in tRNAscan-SE output
        if CODON not in codon2aa:
            codon2aa[CODON] = AA
        ce2loc[TRNA] = [CHR,START,END]
        loc2ce[(CHR,START,END)] = TRNA
        ce2info[TRNA] = [AA,CODON]
        if "chr" in CHR or "Chr" in CHR: # get rid of unnecessary/unmatched 'Chrs'
            CHROM = CHR.split('r')[1]
        else:
            CHROM = CHR
        if CHROM not in chr2start:
            chr2start[CHROM] = []
            chr2end[CHROM] = []
        chr2start[CHROM].append(START)
        chr2end[CHROM].append(END)

    STATEMENT = "Reference tRNA information gathered."
    print(STATEMENT)
    LOG_LIST.append(STATEMENT + '\n')

    ## Assemble variant info.
    vcf_reader = vcf.Reader(open(VCF_FILE, 'r'))

    SAMPLES = vcf_reader.samples

    i = 0

    CHR = ''

    for record in vcf_reader:
        CHROM = str(record.CHROM)
        POS = int(record.POS)
        REF = str(record.REF)
        ALT = str(record.ALT[0])
        if CHROM == "MtDNA":
            break
        if CHROM != CHR:
            CHR = CHROM
            i = 0
            chr2start[CHR].sort()
            chr2end[CHR].sort()
            STATEMENT = "Chromosome " + CHR + " beginning..."
            print(STATEMENT)
            LOG_LIST.append(STATEMENT + '\n')
        if i >= len(chr2start[CHR]):
            continue
        if POS < chr2start[CHR][i]:
            continue
        if POS <= chr2end[CHR][i] and POS >= chr2start[CHR][i]:
            loc2allele[(CHROM,POS)] = [REF,ALT]
            loc2hr[(CHROM,POS)] = []
            loc2ha[(CHROM,POS)] = []
            loc2het[(CHROM,POS)] = []
            for sample in SAMPLES:
                if record.genotype(sample)['GT'] == "0/0":
                    loc2hr[(CHROM,POS)].append(sample)
                elif record.genotype(sample)['GT'] == "1/1":
                    loc2ha[(CHROM,POS)].append(sample)
                else:
                    loc2het[(CHROM,POS)].append(sample)
                    STATEMENT = "%s\t%d\t%s\t%s" % (CHROM,POS, sample, record.genotype(sample)['GT'])
                    LOG_LIST.append(STATEMENT + '\n')
        if POS > chr2end[CHR][i]:
            i += 1
        if i >= len(chr2start[CHR]):
            continue
        elif POS > chr2end[CHR][i]:
            i += 1
        if i >= len(chr2start[CHR]):
            continue
        elif POS > chr2end[CHR][i]:
            i += 1
        if i >= len(chr2start[CHR]):
            continue
        elif POS > chr2end[CHR][i]:
            i += 1
        if i >= len(chr2start[CHR]):
            continue
        elif POS > chr2end[CHR][i]:
            i += 1
        if i >= len(chr2start[CHR]):
            continue
        elif POS > chr2end[CHR][i]:
            i += 1

    STATEMENT = "Variant information assembled for strains."
    print(STATEMENT)
    LOG_LIST.append(STATEMENT + '\n')
    ## Print outfiles.
    # Variant file.
    outfn = args.strOut + "_strain_variants.txt"
    outfile = open(outfn, 'w')

    PREAMBLE = "##File generated by [gitrepo]wormtrna/updatedinitial/get_strain_variants.py.\n"
    outfile.write(PREAMBLE)

    HEADER = "##Chr\tPos\tRef\tAlt\thomRef\thomAlt\thet\n"
    outfile.write(HEADER)

    for (CHR,POS) in loc2allele:
        line = CHR + "\t" + str(POS) + "\t"
        line += ("\t").join(loc2allele[(CHR,POS)]) + "\t"
        if len(loc2hr[(CHR,POS)]) > 0:
            line += (",").join(loc2hr[(CHR,POS)]) + "\t"
        else:
            line += "None" + "\t"
        if len(loc2ha[(CHR,POS)]) > 0:
            line += (",").join(loc2ha[(CHR,POS)]) + "\t"
        else:
            line += "None" + "\t"
        if len(loc2het[(CHR,POS)]) > 0:
            line += (",").join(loc2het[(CHR,POS)])
        else:
            line += "None"
        outfile.write(line + "\n")

    STATEMENT = "%s written." % (outfn)
    print(STATEMENT)
    LOG_LIST.append(STATEMENT + "\n")
    # Log file.
    outfn = args.strOut + "_get_strain_variants.log"
    outfile = open(outfn, 'w')

    PREAMBLE = "##File generated by [gitrepo]wormtrna/updatedinitial/get_strain_variants.py.\n"
    outfile.write(PREAMBLE)

    PREAMBLE = "##Written on %s.\n" % time.ctime()
    outfile.write(PREAMBLE)

    for line in LOG_LIST:
        outfile.write(line)

    outfile.close()

    STATEMENT = "%s written." % (outfn)
    print(STATEMENT)

if __name__=="__main__":
    main()
