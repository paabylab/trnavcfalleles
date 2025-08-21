#!/usr/bin/env python3.6
"""
Avery Davis Bell updates to Corinne Simonti original script
Gets per-allele-information - number strains with each allele, etc
Files generated: allele_info.txt

"""
# Initially written by Corinne Simonti
## imports
import argparse
import os
import sys
import matplotlib as mpl
#import matplotlib.pyplot as plt
#import matplotlib.ticker as ticker
#from itertools import cycle
#import seaborn as sns
import pandas as pd
import numpy as np
import math
#from matplotlib.backends.backend_pdf import PdfPages
import scipy
from scipy import stats
import subprocess
import gzip
import random
import datetime
import time
import getopt
import re
import vcf # called vcf, installed as pyvcf

def main():
    #### Parse arguments
    argp = argparse.ArgumentParser(prog = "allele_info.py",
                                   description = "Avery Davis Bell updates to Corinne Simonti original script \
                                   Gets per-allele-information - number strains with each allele, etc \
                                   Files generated: allele_info.txt")
    argp.add_argument("-outstraintrnas", dest = "strStrainOut", required = True, type = str,
                      metavar = "strain-filtered.out",
                      help = "Path to tRNAscan-SE .out file for strain-specific tRNAs")
    argp.add_argument("-straintrninfo", dest = "strStrainInfo", required = True, type = str,
                      metavar = "strain_trnas_info.txt",
                      help = "Path to *strain_trnas_info.txt output file of build_alt_sequences.py")
    argp.add_argument("-vcf", dest = "strVCF", required = True, type = str,
                      metavar = "file.vcf",
                      help = "Path to VCF. Only used here to get the initial/full sample list.")
    argp.add_argument("-out", dest = "strOut", required = True, type = str,
                      metavar = "outfilestem",
                      help = "Prefix for output files")
    args = argp.parse_args()

    # Assign args to names Corinne uses (instead of updating whole script)
    OUT_FILE = args.strStrainOut
    INFO_FILE = args.strStrainInfo

    ##### Corinne's whole script (wasn't in functions)
    STATEMENT = "Starting at %s..." % time.ctime()
    print(STATEMENT)

    ## Initialize Dictionaries and lists.
    LOG_LIST = []
    write_list = []

    STRAIN_LIST = []

    strain2allele = {}
    trna2aa = {}    # aa, anticodon
    trna2best = {}
    allele2count = {}
    allele2score = {}    # Inf, HMM, 2'Str, Iso....ADB ony has Inf and Iso scores
    allele2aa = {}    # aa, anticodon
    allele2iso = {}    # Iso CM, Iso score
    allele2note = {}
    codon2aa = {}

    strain2invar = {}
    strain2missing = {}
    strain2lost = {}
    strain2alter = {}
    strain2inf = {}
    strain2hmm = {}
    strain2struct = {}

    LOG_LIST.append(STATEMENT + '\n')

    STATEMENT = "initial Dictionaries and lists generated."
    print(STATEMENT)
    LOG_LIST.append(STATEMENT + '\n')
    
    ## get strain info from vcf (earlier # strains was hardcoded)
    vcf_reader = vcf.Reader(open(args.strVCF, 'r'))
    vcf_sample_list = vcf_reader.samples
    
    ## Get tRNA info.
    for line in open(OUT_FILE):
        if line[0] == "S" or line[0] == "N" or line[0] == "-":
            continue
        seg = line.split("\t")
        NAME = str(seg[0])
        TRNA = ("_").join(NAME.split("_")[:-2]) # TRNA = NAME.split("_")[0]
        STRAIN = NAME.split("_")[-2]
        NUM = int(seg[1])
        START = int(seg[2])
        END = int(seg[3])
        AA = str(seg[4])
        CODON = str(seg[5])
        INF = float(seg[8]) # correct
        # HMM = float(seg[9]) # I do not have HMM!!
        # STR = float(seg[10]) # not sure what this is ....
        CM = str(seg[9]) # CM = str(seg[11]) # Iso CM is field 9 now
        ISO = float(seg[10]) # ISO = float(seg[12]) # ISO score is field 10 now 
        NOTE = str(seg[-1])
        if CODON not in codon2aa:
            codon2aa[CODON] = AA
        allele2count[(TRNA,STRAIN)] = NAME.split('_')[2]
        allele2score[(TRNA,STRAIN)] = [INF] # [INF,HMM,STR]
        allele2aa[(TRNA,STRAIN)] = [AA,CODON]
        allele2iso[(TRNA,STRAIN)] = [CM,ISO]
        allele2note[(TRNA,STRAIN)] = NOTE.strip('\n')
        if TRNA not in trna2best:
            trna2best[TRNA] = [STRAIN,INF]
        else:
            if INF > trna2best[TRNA][1]:
                trna2best[TRNA] = [STRAIN,INF]
            elif INF == trna2best[TRNA][1]:
                if int(allele2count[(TRNA,STRAIN)]) > int(allele2count[(TRNA,trna2best[TRNA][0])]):
                    trna2best[TRNA] = [STRAIN,INF]
    
           
    ## Add tRNA info for seqs in the strain info file that didn't make full tRNAs found by tRNAscan-SE
    istm = open(INFO_FILE)
    for line in istm:
        if line[0] == "#":
            continue
        seg = line.split('\t')
        newTRNA = seg[3]
        STRAINS = seg[-1].rstrip('\n').split(',')
        present = 0
        if newTRNA not in trna2best: # means it can't be/have best - wasn't found by trnascan-SE
            trna2best[newTRNA] = ["NA",0] # new - check here if something chokes
        for STRAIN in STRAINS:
            for (TRNA, hSTRAIN) in allele2score:
                if newTRNA == TRNA and hSTRAIN==STRAIN:
                    present += 1
        # Add this strain-tRNA pair if they aren't included
        if present==0:
            allele2score[(newTRNA,STRAIN[0])] = [0]
            allele2aa[(newTRNA,STRAIN[0])] = ["Undet","NNN"]
            allele2iso[(newTRNA,STRAIN[0])] = ["Undet", 0]
            allele2note[(newTRNA,STRAIN[0])] = "undetermined isotype / not found by tRNAscan-SE"
            allele2count[(newTRNA,STRAIN[0])] = str(len(STRAINS))
    istm.close()

    STATEMENT = "tRNA information gathered."
    print(STATEMENT)
    LOG_LIST.append(STATEMENT + '\n')
    ## Print outfiles.
    # Allele file.
    outfn = args.strOut + "_allele_info.txt"
    outfile = open(outfn, 'w')

    PREAMBLE = "##File generated by [gitrepo]wormtrna/updatedinitial/allele_info.py.\n"
    outfile.write(PREAMBLE)

    HEADER = "##tRNA\t#strains\tAA\tCodon\tInfernal\tNote\tLost\talleleCM\n" # "##tRNA\t#strains\tAA\tCodon\tInfernal\tHMM\t2'Str\tNote\tLost\talleleCM\n"
    outfile.write(HEADER)

    for (TRNA,STRAIN) in allele2count:
        LOST = "Functional"
        if trna2best[TRNA][0] == STRAIN:
            LOST = "Best"
        if allele2count[(TRNA,STRAIN)] == str(len(vcf_sample_list) + 1): # need to figure out why this +1 is needed
            LOST = "Invariant"
        AA = allele2aa[(TRNA,STRAIN)][0]
        INF = allele2score[(TRNA,STRAIN)][0]
        # HMM = allele2score[(TRNA,STRAIN)][1]
        # STR = allele2score[(TRNA,STRAIN)][2]
        if "IPD" in allele2note[(TRNA,STRAIN)]:
            LOST = "Altered"
        if "secondary" in allele2note[(TRNA,STRAIN)] or "undetermined isotype" in allele2note[(TRNA,STRAIN)] or "pseudo" in allele2note[(TRNA,STRAIN)]:
            LOST = "Lost"
        write_line = TRNA + "\t" + allele2count[(TRNA,STRAIN)] + "\t" + ("\t").join(allele2aa[(TRNA,STRAIN)])
        if INF != "NA":
            write_line += "\t%.1f" % (INF) # "\t%.1f\t%.1f\t%.1f" % (INF,HMM,STR)
        else:
            write_line += "\t%s" % (INF) # "\t%s\t%s\t%s" % (INF,HMM,STR)
        if (TRNA,STRAIN) in allele2note:
            write_line += "\t" + allele2note[(TRNA,STRAIN)]
        else:
            write_line += "\tNA"
        write_line += "\t" + LOST
        write_line += "\t" + allele2iso[(TRNA,STRAIN)][0] + "\n"
        outfile.write(write_line)

    STATEMENT = "%s written." % (outfn)
    print(STATEMENT)
    LOG_LIST.append(STATEMENT + "\n")
    # Log file.
    outfn = args.strOut + "_allele_info.log"
    outfile = open(outfn, 'w')

    PREAMBLE = "##File generated by [gitrepo]wormtrna/updatedinitial/allele_info.py.\n"
    outfile.write(PREAMBLE)

    PREAMBLE = "##Written on %s.\n" % datetime.date.today().isoformat()
    outfile.write(PREAMBLE)

    for line in LOG_LIST:
        outfile.write(line)

    outfile.close()

    STATEMENT = "%s written." % (outfn)
    print(STATEMENT)



if __name__=="__main__":
    main()
