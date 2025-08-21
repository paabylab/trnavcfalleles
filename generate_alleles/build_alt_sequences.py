#!/usr/bin/env python3.6
"""
Avery Davis Bell updates to Corinne Simonti original script
This script generates file strain_trnas_gen.txt, strain_trnas.fa, strain_trnas_info.txt
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
import os
import math
#from matplotlib.backends.backend_pdf import PdfPages
import scipy
from scipy import stats
import subprocess
import gzip
import random
import subprocess
import datetime
import time
import getopt
import re
import vcf # called vcf, installed as pyvcf


def main():
    #### Parse arguments
    argp = argparse.ArgumentParser(prog = "build_alt_sequences.py",
                                   description = "Avery Davis Bell updates to Corinne Simonti original script. \
                                       This script generates strain_trnas_gen.txt (the strain-specific tRNA sequences, I think),\
                                           strain_trnas.fa, strain_trnas_info.txt")
    argp.add_argument("-trnasout", dest = "strConfOut", required = True, type = str,
                      metavar = "tRNAs-confidence-set.out",
                      help = "Path to tRNAscan-SE .out file for high confidence set, e.g. ce11-tRNAs-confidence-set.out")
    argp.add_argument("-trnass", dest = "strConfSS", required = True, type = str,
                      metavar = "tRNAs-confidence-set.ss",
                      help = "Path to tRNAscan-SE .ss file for high confidence set, e.g. ce11-tRNAs-confidence-set.ss")
    argp.add_argument("-trnabed", dest = "strBed", required = True, type = str,
                      metavar = "tRNAs.bed",
                      help = "Path to tRNAscan-SE .bed file for tRNAs, e.g. ce11-tRNAs.bed")
    argp.add_argument("-vcf", dest = "strVCF", required = True, type = str,
                      metavar = "variants.vcf",
                      help = "Path to VCF file for species whos tRNAs are provided (just used for sample info here)")
    argp.add_argument("-vars", dest = "strVars", required = True, type = str,
                      metavar = "strain_variants.txt",
                      help = "Path to *_strain_variants.txt output of get_strain_variants.py")
    argp.add_argument("-out", dest = "strOut", required = True, type = str,
                      metavar = "outfilestem",
                      help = "Prefix for output files")
    args = argp.parse_args()

    # Assign args to names Corinne uses (instead of updating whole script)
    CHAN_OUT_FILE = args.strConfOut
    CHAN_SS_FILE = args.strConfSS
    CHAN_BED_FILE = args.strBed
    VCF_FILE = args.strVCF
    VAR_FILE = args.strVars

    ##### Corinne's whole script (wasn't in functions)
    STATEMENT = "Starting at %s..." % time.ctime()
    print(STATEMENT)

    ## Set up Dictionaries and lists.
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

    strain2seq = {}

    nuc2match = {}
    nuc2match["A"] = "T"
    nuc2match["C"] = "G"
    nuc2match["G"] = "C"
    nuc2match["T"] = "A"
    nuc2match["."] = "."

    pos2seq = {}
    var2pos = {}
    strain_pos2seq = {}
    seq2strain = {}

    trna2seq = {}

    LOG_LIST.append(STATEMENT + '\n')

    STATEMENT = "Empty/initial dictionaries and lists generated."
    print(STATEMENT)
    LOG_LIST.append(STATEMENT + '\n')
    ## Get tRNA info.
    for line in open(CHAN_OUT_FILE):
        if line[0] == "S" or line[0] == "N" or line[0] == "-":
            continue
        seg = line.split("\t")
        CHR = str(seg[0]).strip()
        if "chr" in CHR or "Chr" in CHR: # get rid of unnecessary/unmatched 'Chrs'
            CHROM = CHR.split('r')[1]
        else:
            CHROM = CHR
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
        # NOTE = str(seg[14]) # assigned but not used and was causing issues
        if CODON not in codon2aa:
            codon2aa[CODON] = AA
        ce2loc[TRNA] = [CHROM,START,END]
        loc2ce[(CHROM,START,END)] = TRNA
        ce2info[TRNA] = [AA,CODON]


    for line in open(CHAN_SS_FILE):
        seg = line.split()
        if line == "\n":
            continue
        if "trna" in line: # if line[:3] == "chr": # fix this so it works if starts with chr OR NOT
            TRNA = seg[0]
        if seg[0] == "Seq:":
            SEQ = seg[-1].strip('\n')
            ce2seq[TRNA] = SEQ.upper()


    STATEMENT = "Reference tRNA information gathered."
    print(STATEMENT)
    LOG_LIST.append(STATEMENT + '\n')
    ## Get strand info.
    for line in open(CHAN_BED_FILE):
        seg = line.split('\t')
        CHR = str(seg[0])
        CHROM = CHR # CHR.split('r')[1] # this split was for when bed file had chr and vcf didn't; in workflow I've editted bed file!!
        START = int(seg[1]) + 1
        END = int(seg[2])
        STRAND = str(seg[5])
        if (CHROM, START, END) in loc2ce:
            ce2strand[loc2ce[(CHROM,START,END)]] = STRAND


    STATEMENT = "Strand information gathered."
    print(STATEMENT)
    LOG_LIST.append(STATEMENT + '\n')
    ## Assemble variant info.
    vcf_reader = vcf.Reader(open(VCF_FILE, 'r'))
    STRAIN_LIST = vcf_reader.samples

    for line in open(VAR_FILE):
        if line[0] == "#":
            continue
        seg = line.split('\t')
        CHR = str(seg[0])
        POS = int(seg[1])
        REF = str(seg[2])
        ALT = str(seg[3])
        var2pos[(CHR,POS)] = [REF,ALT]
        HR = seg[4].split(',')
        if len(HR) != 0:
            for strain in HR:
                if len(REF) == 1:
                    strain_pos2seq[(strain,CHR,POS)] = REF
                else:
                    for i in range(len(REF)):
                        strain_pos2seq[(strain,CHR,POS+i)] = REF[i]
                        var2pos[(CHR,POS+i)] = [REF[i],"-"]
        HA = seg[5].split(',')
        if len(HA) != 0:
            for strain in HA:
                if len(REF) == 1:
                    strain_pos2seq[(strain,CHR,POS)] = ALT
                else:
                    for i in range(len(REF)):
                        if i == 0:
                            strain_pos2seq[(strain,CHR,POS)] = ALT
                        else:
                            strain_pos2seq[(strain,CHR,POS+i)] = "-"
        HET = seg[6].strip('\n').split(',')
        if len(HET) != 0:
            for strain in HET:
                if len(REF) == 1:
                    strain_pos2seq[(strain,CHR,POS)] = "."
                else:
                    for i in range(len(REF)):
                        if i == 0:
                            strain_pos2seq[(strain,CHR,POS)] = "."
                        else:
                            strain_pos2seq[(strain,CHR,POS+i)] = "."

    STATEMENT = "Variant information assembled for strains."
    print(STATEMENT)
    LOG_LIST.append(STATEMENT + '\n')
    ## Map sequence to position.
    for trna in ce2loc:
        CHR = ce2loc[trna][0]
        if ce2strand[trna] == "+":
            i = 0
            POS = ce2loc[trna][1]
            END = ce2loc[trna][2]
            while POS <= END:
                pos2seq[(CHR,POS)] = ce2seq[trna][i]
                i += 1
                POS += 1
        else:
            i = -1
            POS = ce2loc[trna][1]
            END = ce2loc[trna][2]
            while POS <= END:
                pos2seq[(CHR,POS)] = ce2seq[trna][i]
                i -= 1
                POS += 1

    STATEMENT = "Sequences mapped to genomic positions."
    print(STATEMENT)
    LOG_LIST.append(STATEMENT + '\n')
    ## Build sequences for other strains.
    for trna in ce2loc:
        for strain in STRAIN_LIST:
            strain2seq[(strain,trna)] = ""
        CHR = ce2loc[trna][0]
        START = ce2loc[trna][1]
        END = ce2loc[trna][2]
        if ce2strand[trna] == "+":
            POS = START
            while POS <= END:
                if (CHR,POS) in var2pos:
                    REF = var2pos[(CHR,POS)][0]
                    ALT = var2pos[(CHR,POS)][1]
                    if REF != pos2seq[(CHR,POS)]:
                        STATEMENT = "ERROR: Reference (" + pos2seq[(CHR,POS)] + ") does not match strain ref: " + REF + " at " + CHR + ":" + str(POS)
                        print(STATEMENT)
                        LOG_LIST.append(STATEMENT + '\n')
                    for strain in STRAIN_LIST:
                        strain2seq[(strain,trna)] += strain_pos2seq[(strain,CHR,POS)]
                else:
                    for strain in STRAIN_LIST:
                        strain2seq[(strain,trna)] += pos2seq[(CHR,POS)]
                POS += 1
        else:
            POS = END
            while POS >= START:
                if (CHR,POS) in var2pos:
                    REF = var2pos[(CHR,POS)][0]
                    ALT = var2pos[(CHR,POS)][1]
                    if ALT != "-" and nuc2match[REF] != pos2seq[(CHR,POS)]:
                        STATEMENT = "ERROR: Reference (" + pos2seq[(CHR,POS)] + ") does not match strain ref: " + nuc2match[REF] + " at " + CHR + ":" + str(POS)
                        print(STATEMENT)
                        LOG_LIST.append(STATEMENT + '\n')
                    for strain in STRAIN_LIST:
                        if strain_pos2seq[(strain,CHR,POS)] == ALT:
                            if ALT != "-":
                                if len(ALT) == 1:
                                    strain2seq[(strain,trna)] += nuc2match[strain_pos2seq[(strain,CHR,POS)]]
                                elif len(ALT) > 1:
                                    i = -1
                                    POS1 = 1
                                    LEN = len(ALT)
                                    while POS1 <= LEN:
                                        strain2seq[(strain,trna)] += nuc2match[ALT[i]]
                                        i -= 1
                                        POS1 += 1
                                else:
                                    print("ERROR: alternate allele is empty string.")
                            else:
                                strain2seq[(strain,trna)] += strain_pos2seq[(strain,CHR,POS)]
                        else:
                            strain2seq[(strain,trna)] += nuc2match[strain_pos2seq[(strain,CHR,POS)]]
                else:
                    for strain in STRAIN_LIST:
                        strain2seq[(strain,trna)] += pos2seq[(CHR,POS)]
                POS -= 1


    STATEMENT = "Sequences built."
    print(STATEMENT)
    LOG_LIST.append(STATEMENT + '\n')
    ## Collate identical sequences.
    for trna in ce2info:
        N2 = ce2seq[trna] # I think OK to use N2 here - it's a variable not actually the strain (strain only works for elegans)
        seq2strain[(trna,N2)] = ["Reference"]
        trna2seq[trna] = [N2]
        for strain in STRAIN_LIST:
            SEQ = strain2seq[(strain,trna)]
            if "." in SEQ:
                continue
            if (trna,SEQ) in seq2strain:
                seq2strain[(trna,SEQ)].append(strain)
            else:
                trna2seq[trna].append(SEQ)
                seq2strain[(trna,SEQ)] = [strain]

    STATEMENT = "Unique sequences assembled."
    print(STATEMENT)
    LOG_LIST.append(STATEMENT + '\n')
    ## Print outfiles.
    # Sequence file.
    outfn = args.strOut + "_strain_trnas_gen.txt"
    outfile = open(outfn, 'w')

    PREAMBLE = "##File generated by [gitrepo]wormtrna/updatedinitial/build_alt_sequences.py.\n"
    outfile.write(PREAMBLE)

    HEADER = "##Chr\tStart\tEnd\ttRNA\tStrand\tAA\tCodon\tnumVersions\tnumStrains\tMissing\n"
    outfile.write(HEADER)

    for trna in ce2info:
        write_line = ce2loc[trna][0] + "\t" + str(ce2loc[trna][1]) + "\t" + str(ce2loc[trna][2]) + "\t"
        write_line += trna + "\t" + ce2strand[trna] + "\t"
        write_line += ("\t").join(ce2info[trna]) + "\t"
        write_line += str(len(trna2seq[trna])) + "\t"
        NUM_STRAINS = []
        STRNUM_STRAINS = []
        for seq in trna2seq[trna]:
            NUM_STRAINS.append(len(seq2strain[(trna,seq)]))
            STRNUM_STRAINS.append(str(len(seq2strain[(trna,seq)])))
        write_line += (',').join(STRNUM_STRAINS) + "\t"
        write_line += str(len(STRAIN_LIST) - sum(NUM_STRAINS))
        outfile.write(write_line + '\n')

    STATEMENT = "%s written." % (outfn)
    print(STATEMENT)
    LOG_LIST.append(STATEMENT + "\n")
    # Fasta file.
    outfn = args.strOut + "_strain_trnas.fa"
    outfile = open(outfn, 'w')

    for trna in ce2info:
        for seq in trna2seq[trna]:
            info_line = ">" + trna + "_" + seq2strain[(trna,seq)][0] + "_" + str(len(seq2strain[(trna,seq)]))
            outfile.write(info_line + "\n")
            fa_line = seq
        # NEW: deal with cases where seq has --- in it
            if "-" in fa_line:
                STATEMENT = "Sequence updated to exclude -s originally in it (presumably from indels): " + info_line + " (original sequence: " + fa_line + ")"
                print(STATEMENT)
                LOG_LIST.append(STATEMENT + "\n")
                new_line = ""
                for char in fa_line:
                    if char!="-":
                        new_line = new_line + char
                        fa_line = new_line
            outfile.write(fa_line + "\n")


    STATEMENT = "%s written." % (outfn)
    print(STATEMENT)
    LOG_LIST.append(STATEMENT + "\n")
    # Fasta file.
    outfn = args.strOut + "_strain_trnas_info.txt"
    outfile = open(outfn, 'w')

    PREAMBLE = "##File generated by [gitrepo]wormtrna/updatedinitial/build_alt_sequences.py.\n"
    outfile.write(PREAMBLE)

    HEADER = "##Chr\tStart\tEnd\ttRNA\tStrand\tAA\tCodon\tnumStrains\tStrains\n"
    outfile.write(HEADER)

    for trna in ce2info:
        for seq in trna2seq[trna]:
            write_line = ce2loc[trna][0] + "\t" + str(ce2loc[trna][1]) + "\t" + str(ce2loc[trna][2]) + "\t"
            write_line += trna + "\t" + ce2strand[trna] + "\t"
            write_line += ("\t").join(ce2info[trna]) + "\t"
            write_line += str(len(seq2strain[(trna,seq)])) + "\t"
            write_line += (",").join(seq2strain[(trna,seq)])
            outfile.write(write_line + '\n')


    STATEMENT = "%s written." % (outfn)
    print(STATEMENT)
    LOG_LIST.append(STATEMENT + "\n")
    # Log file.
    outfn = args.strOut + "_build_alt_seqs.log"
    outfile = open(outfn, 'w')

    PREAMBLE = "##File generated by [gitrepo]wormtrna/updatedinitial/build_alt_sequences.py.\n"
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
