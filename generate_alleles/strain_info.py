#!/usr/bin/env python3.6
"""
Avery Davis Bell updates to Corinne Simonti original script
Gets per-strain information - number tRNAs with mutations, that sort of thing
Files generated: strain_summary.txt ; strain_summary_wrtRef.txt ; strain_variable_allele_info.txt; variant_uniqueallele.txt

"""
# Initially written by Corinne Simonti
## imports
import argparse
import sys
import matplotlib as mpl
#import matplotlib.pyplot as plt
#import matplotlib.ticker as ticker
#from itertools import cycle
#import seaborn as sns
import pandas as pd
import numpy as np
import os
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
    argp = argparse.ArgumentParser(prog = "strain_info.py",
                                   description = "Avery Davis Bell updates to Corinne Simonti original script \
                                   Gets per-strain information - number tRNAs with mutations, that sort of thing \
                                   Files generated: strain_summary.txt ; strain_summary_wrtRef.txt ; strain_variable_allele_info.txt ; variant_uniqueallele.txt")
    argp.add_argument("-outstraintrnas", dest = "strStrainOut", required = True, type = str,
                      metavar = "strain-filtered.out",
                      help = "Path to tRNAscan-SE .out file for strain-specific tRNAs")
    argp.add_argument("-straintrnass", dest = "strStrainSS", required = True, type = str,
                      metavar = "strain-filtered.ss",
                      help = "Path to tRNAscan-SE .ss file for strain-specific tRNAs")
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
    SS_FILE = args.strStrainSS
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
    allele2score = {}    # Inf, HMM, 2'Str, Iso.....# ADB only has Inf and Iso scroes!!
    allele2aa = {}    # aa, anticodon
    allele2iso = {}    # Iso CM, Iso score
    allele2note = {}
    codon2aa = {}

    strain2invar = {}
    strain2missing = {}
    strain2lost = {}
    strain2alter = {}
    strain2inf = {}
    # strain2hmm = {} # no HMM, 2' structure in my runs...
    # strain2struct = {}

    strain2best = {}
    strain2func = {}
    strain2alter2 = {}

    st_aa2count = {}
    st_codon2count = {}
    st_aa2mm = {}
    st_codon2mm = {}

    codon2seq = {}
    ua2info = {} # Name, AA, Codon, Iso

    LOG_LIST.append(STATEMENT + '\n')

    STATEMENT = "Dictionaries and lists generated."
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
        INF = float(seg[8])
        # HMM = float(seg[9])
        # STR = float(seg[10])
        INF = float(seg[8]) # correct
        # HMM = float(seg[9]) # I do not have HMM!!
        # STR = float(seg[10]) # not sure what this is ....
        CM = str(seg[9]) # CM = str(seg[11]) # Iso CM is field 9 now
        ISO = float(seg[10]) # ISO = float(seg[12]) # ISO score is field 10 now
        NOTE = str(seg[-1])
        if CODON not in codon2aa:
            codon2aa[CODON] = AA
        allele2count[(TRNA,STRAIN)] = NAME.split('_')[-1].strip()
        allele2score[(TRNA,STRAIN)] = [INF] # [INF,HMM,STR]
        allele2aa[(TRNA,STRAIN)] = [AA,CODON]
        allele2iso[(TRNA,STRAIN)] = [CM,ISO]
        allele2note[(TRNA,STRAIN)] = NOTE.strip('\n')
        if TRNA not in trna2best:
            trna2best[TRNA] = INF
        elif INF > trna2best[TRNA]:
            trna2best[TRNA] = INF
        if "secondary" in NOTE or "undetermined isotype" in NOTE or "pseudo" in NOTE or "unexpected" in NOTE:
            continue
        if AA == "iMet":
            CODON = "iCAT"
        st_aa2count.setdefault((STRAIN,AA),0)
        st_codon2count.setdefault((STRAIN,CODON),0)
        st_aa2mm.setdefault((STRAIN,AA),0)
        st_codon2mm.setdefault((STRAIN,CODON),0)
        st_aa2count[(STRAIN,AA)] += 1
        st_codon2count[(STRAIN,CODON)] += 1
        if "IPD" in NOTE:
            st_aa2mm[(STRAIN,AA)] += 1
            st_codon2mm[(STRAIN,CODON)] += 1


    STATEMENT = "tRNA information gathered."
    print(STATEMENT)
    LOG_LIST.append(STATEMENT + '\n')
    ## Get strain info.
    for line in open(INFO_FILE):
        if line[0] == "#":
            continue
        seg = line.split('\t')
        TRNA = seg[3]
        AA = str(seg[5])
        CODON = str(seg[6])
        COUNT = str(seg[7]) # newly added!
        trna2aa[TRNA] = [AA,CODON]
        STRAINS = seg[-1].rstrip('\n').split(',')
        if TRNA not in trna2best: # means it can't be/have best - wasn't found by trnascan-SE
            trna2best[TRNA] = [0] # new - check here if something chokes
        for STRAIN in STRAINS:
            strain2allele[(STRAIN,TRNA)] = STRAINS[0]
            if STRAIN not in STRAIN_LIST:
                STRAIN_LIST.append(STRAIN)
            if (TRNA,STRAIN) not in allele2score:
                """
                if (TRNA,STRAIN) == ('chrX.trna232', 'ECA744'):
                    allele2score[(TRNA,STRAIN)] = [0,0,0]
                    allele2aa[(TRNA,STRAIN)] = ["Undet","NNN"]
                    allele2iso[(TRNA,STRAIN)] = ["Undet",0]
                    allele2note[(TRNA,STRAIN)] = "undetermined isotype"
                else:
                """
                # Set up if in info file but not tRNAscan-SE output (e.g., really short seqs) - this is what the specific edge case above was dealing with
                if (TRNA,STRAINS[0]) not in allele2score:
                    allele2score[(TRNA,STRAIN)] = [0]
                    allele2aa[(TRNA,STRAIN)] = ["Undet","NNN"]
                    allele2iso[(TRNA,STRAIN)] = ["Undet",0]
                    allele2note[(TRNA,STRAIN)] = "undetermined isotype / not found by tRNAscan-SE" 
                    allele2count[(TRNA, STRAIN)] = COUNT # formerly count wasn't happening! was only happening for ones in trnascan-SE out!
                else: # Others (normal)
                    allele2score[(TRNA,STRAIN)] = allele2score[(TRNA,STRAINS[0])]
                    allele2aa[(TRNA,STRAIN)] = allele2aa[(TRNA,STRAINS[0])]
                    allele2iso[(TRNA,STRAIN)] = allele2iso[(TRNA,STRAINS[0])]
                    allele2note[(TRNA,STRAIN)] = allele2note[(TRNA,STRAINS[0])]
                    AA = allele2aa[(TRNA,STRAINS[0])][0]
                    CODON = allele2aa[(TRNA,STRAINS[0])][1]
                    NOTE = allele2note[(TRNA,STRAINS[0])]
                    if "secondary" in NOTE or "undetermined isotype" in NOTE or "pseudo" in NOTE or "unexpected" in NOTE:
                        continue
                    if AA == "iMet":
                        CODON = "iCAT"
                    st_aa2count.setdefault((STRAIN,AA),0)
                    st_codon2count.setdefault((STRAIN,CODON),0)
                    st_aa2mm.setdefault((STRAIN,AA),0)
                    st_codon2mm.setdefault((STRAIN,CODON),0)
                    if "IPD" not in NOTE:
                        st_aa2count[(STRAIN,AA)] += 1
                        st_codon2count[(STRAIN,CODON)] += 1
                    else:
                        st_aa2mm[(STRAIN,AA)] += 1
                        st_codon2mm[(STRAIN,CODON)] += 1

    STATEMENT = "Strain allele information gathered."
    print(STATEMENT)
    LOG_LIST.append(STATEMENT + '\n')
    ## Get sequence info.
    for line in open(SS_FILE):
        seg = line.split()
        if line == "\n":
            continue
        if "trna" in seg[0]: # updated to deal with cases where _ is in contig names. still won't work if _ is in strain name!
            NAME = (".").join(seg[0].split('.')[:-1])
            TRNA = ("_").join(NAME.split("_")[:-2]) # TRNA = NAME.split("_")[0]
            STRAIN = NAME.split("_")[-2]
        NOTE = allele2note[(TRNA,STRAIN)]
        # if "secondary" in NOTE or "undetermined isotype" in NOTE or "pseudo" in NOTE or "unexpected" in NOTE:
        #    continue
        AA = allele2aa[(TRNA,STRAIN)][0]
        CODON = allele2aa[(TRNA,STRAIN)][1]
        if AA == "iMet":
            CODON = "iCAT"
        if seg[0] == "Seq:":
            SEQ = seg[-1].strip('\n')
            codon2seq.setdefault(CODON,[]).append(SEQ)


    ## Print outfiles.
    # Allele file (by strain).
    outfn = args.strOut + "_strain_variable_allele_info.txt"
    outfile = open(outfn, 'w')

    PREAMBLE = "##File generated by [gitrepo]wormtrna/updatedinitial/strain_info.py.\n"
    outfile.write(PREAMBLE)

    HEADER = "##Strain\ttRNA\tAA\tCodon\tInfernal\tNote\tLost\talleleCM\n"
    outfile.write(HEADER)

    for STRAIN in STRAIN_LIST:
        strain2invar[STRAIN] = 0
        strain2missing[STRAIN] = 0
        strain2lost[STRAIN] = 0
        strain2alter[STRAIN] = 0
        strain2best[STRAIN] = 0
        strain2alter2[STRAIN] = 0
        strain2inf[STRAIN] = []
        # strain2hmm[STRAIN] = []
        # strain2struct[STRAIN] = []
        
        # The counting in here is not perfect! Missing some somewhere! May want to do my own processing...
        for TRNA in trna2aa:
            if allele2count[(TRNA,"Reference")] == str(len(vcf_sample_list) + 1): # need to figure out why this + 1 is needed
                continue
            if allele2count[(TRNA,"Reference")] == "1": # should this be 2, for 'Reference' and N2?
                continue
            LOST = "Not Lost"
            AA = "NA"
            INF = 'NA'
            # HMM = 'NA'
            # STR = 'NA'
            if (TRNA,STRAIN) in allele2score:
                if allele2score[(TRNA,"Reference")][0]==0: # happens occasionally when pseudos etc are allowed in
                    INF = "NA"
                else:
                    INF = (allele2score[(TRNA,STRAIN)][0] - allele2score[(TRNA,"Reference")][0]) /  allele2score[(TRNA,"Reference")][0]
                # HMM = (allele2score[(TRNA,STRAIN)][1] - allele2score[(TRNA,"Reference")][1]) /  allele2score[(TRNA,"Reference")][1]
                # STR = (allele2score[(TRNA,STRAIN)][2] - allele2score[(TRNA,"Reference")][2]) /  allele2score[(TRNA,"Reference")][2]
                if allele2score[(TRNA,STRAIN)][0] == trna2best[TRNA]:
                    strain2best[STRAIN] += 1
            if INF == 0: # and HMM == 0 and STR == 0:
                strain2invar[STRAIN] += 1
            elif INF == "NA":
                strain2missing[STRAIN] += 1
            else:
                strain2inf[STRAIN].append(INF)
                # strain2hmm[STRAIN].append(HMM)
                # strain2struct[STRAIN].append(STR)
            if (TRNA,STRAIN) in allele2note:
                NOTE = allele2note[(TRNA,STRAIN)]
                if "IPD" in NOTE and "pseudo" not in NOTE and "undetermined isotype" not in NOTE and "secondary" not in NOTE and "unexpected" not in NOTE:
                    strain2alter2[STRAIN] += 1
            if (TRNA,STRAIN) in allele2aa:
                AA = allele2aa[(TRNA,STRAIN)][0]
                if trna2aa[TRNA][0] != AA:
                    LOST = "Altered"
                if "secondary" in allele2note[(TRNA,STRAIN)] or "undetermined isotype" in allele2note[(TRNA,STRAIN)] or "pseudo" in allele2note[(TRNA,STRAIN)] or "unexpected" in allele2note[(TRNA,STRAIN)]:
                    LOST = "Lost"
                    strain2lost[STRAIN] += 1
            else:
                LOST = "Missing"
            if LOST == "Altered":
                strain2alter[STRAIN] += 1
            write_line = STRAIN + "\t" + TRNA + "\t" + ("\t").join(trna2aa[TRNA])
            if INF != "NA":
                write_line += "\t%.5f" % (INF) # "\t%.5f\t%.5f\t%.5f" % (INF,HMM,STR)
            else:
                write_line += "\t%s" % (INF) # "\t%s\t%s\t%s" % (INF,HMM,STR)
            if (TRNA,STRAIN) in allele2note:
                write_line += "\t" + allele2note[(TRNA,STRAIN)]
            else:
                write_line += "\tNA"
            write_line += "\t" + LOST
            write_line += "\t" + AA + "\n"
            outfile.write(write_line)

    STATEMENT = "%s written." % (outfn)
    print(STATEMENT)
    LOG_LIST.append(STATEMENT + "\n")
    # Summary file (w/rt Ref).
    outfn = args.strOut + "_strain_summary_wrtRef.txt"
    outfile = open(outfn, 'w')

    PREAMBLE = "##File generated by [gitrepo]wormtrna/updatedinitial/strain_info.py.\n"
    outfile.write(PREAMBLE)

    HEADER = "##Strain\tInvariant\tMissing\tLost\tAltered\tInf\n" # "##Strain\tInvariant\tMissing\tLost\tAltered\tInf\tHMM\t2'Str\n"
    outfile.write(HEADER)

    for STRAIN in STRAIN_LIST:
        if len(strain2inf[STRAIN]) > 0:
            write_line = "%s\t%d\t%d\t%d\t%d\t%.5f\n" % (STRAIN,strain2invar[STRAIN],strain2missing[STRAIN],strain2lost[STRAIN],
                strain2alter[STRAIN],np.average(strain2inf[STRAIN])) #, np.average(strain2hmm[STRAIN]),np.average(strain2struct[STRAIN]))
        else:
            write_line = "%s\t%d\t%d\t%d\t%d\t0\n" % (STRAIN,strain2invar[STRAIN],strain2missing[STRAIN],strain2lost[STRAIN],strain2alter[STRAIN])
        outfile.write(write_line)

    STATEMENT = "%s written." % (outfn)
    print(STATEMENT)
    LOG_LIST.append(STATEMENT + "\n")
    # Summary file (w/rt tRNA).
    outfn = args.strOut + "_strain_summary.txt"
    outfile = open(outfn, 'w')

    PREAMBLE = "##File generated by [gitrepo]wormtrna/updatedinitial/strain_info.py.\n"
    outfile.write(PREAMBLE)

    HEADER = "##Strain\tBest\tFunctional\tMissing\tLost\tAltered\n"
    outfile.write(HEADER)

    for STRAIN in STRAIN_LIST: ## ** confirm this is right # here; corinne had 371 where I've put len(vcf_sample_list)
        FUNC = (len(vcf_sample_list) + 1) - strain2best[STRAIN] - strain2missing[STRAIN] - strain2lost[STRAIN] - strain2alter2[STRAIN]
        write_line = "%s\t%d\t%d\t%d\t%d\t%d\n" % (STRAIN,strain2best[STRAIN],FUNC,strain2missing[STRAIN],strain2lost[STRAIN],
            strain2alter2[STRAIN])
        outfile.write(write_line)

    STATEMENT = "%s written." % (outfn)
    print(STATEMENT)
    LOG_LIST.append(STATEMENT + "\n")
    # Summary file (w/rt tRNA).
    codon2aa["iCAT"] = "iMet"
    codon2aa["CAT"] = "Met"

    outfn = args.strOut + "_variant_uniqueallele.txt"
    outfile = open(outfn, 'w')

    PREAMBLE = "##File generated by [gitrepo]wormtrna/updatedinitial/strain_info.py.\n"
    outfile.write(PREAMBLE)

    HEADER = "##AA\tAnticodon\tStrain\tAA_count\tAA_mm\tAC_count\tAC_mm\tUnique\n"
    outfile.write(HEADER)

    for STRAIN in STRAIN_LIST:
        for CODON in codon2aa:
            # if CODON == "NNN":
            #    continue
            AA = codon2aa[CODON]
            write_line = AA + "\t" + CODON + "\t" + STRAIN + "\t"
            write_line += str(st_aa2count.setdefault((STRAIN,AA),0)) + "\t" + str(st_aa2mm.setdefault((STRAIN,AA),0)) + "\t"
            write_line += str(st_codon2count.setdefault((STRAIN,CODON),0)) + "\t" + str(st_codon2mm.setdefault((STRAIN,CODON),0)) + "\t"
            write_line += str(len(set(codon2seq[CODON])))
            outfile.write(write_line + "\n")

    STATEMENT = "%s written." % (outfn)
    print(STATEMENT)
    LOG_LIST.append(STATEMENT + "\n")
    # Log file.
    outfn = args.strOut + "_strain_info.log"
    outfile = open(outfn, 'w')

    PREAMBLE = "##File generated by [gitrepo]wormtrna/updatedinitial/strain_info.py.\n"
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
