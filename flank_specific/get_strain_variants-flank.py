#!/usr/bin/env python3.6
"""
Avery Davis Bell updates to Corinne Simonti original script
Gets variants in tRNA flanking regions from a VCF
Files generated: strain_variants-flank.txt, get_strain_variants-flank.log

"""

# Initially written by Corinne Simonti
## imports
import argparse
import os
import sys
import gzip
import random
import subprocess
import datetime
import time
import getopt
import re
import vcf

def main():
    STATEMENT = "Starting at %s..." % time.ctime()
    print(STATEMENT)
    
    #### Parse arguments
    argp = argparse.ArgumentParser(prog = "get_strain_variants-flank.py",
                                   description = "Avery Davis Bell updates to Corinne Simonti original script. \
                                       This script generates strain_variants-flank.txt, get_strain_variants-flank.log")
    argp.add_argument("-trnasout", dest = "strDetailedOut", required = True, type = str,
                      metavar = "tRNAs-detailed.out",
                      help = "Path to tRNAscan-SE .out file for full set, e.g. ce11-tRNAs-detailed.out")
    argp.add_argument("-vcf", dest = "strVCF", required = True, type = str,
                      metavar = "variants.vcf",
                      help = "Path to VCF file for species whos tRNAs are provided")
    argp.add_argument("-out", dest = "strOut", required = True, type = str,
                      metavar = "outfilestem",
                      help = "Prefix for output files")
    args = argp.parse_args()

    # Assign args to names Corinne uses (instead of updating whole script)
    CHAN_OUT_FILE = args.strDetailedOut
    VCF_FILE = args.strVCF

    ##### Corinne's whole script (wasn't in functions)
    ## Dictionaries and lists.
    LOG_LIST = []
    write_list = []

    CEHREP_LIST = []
    CEHSKIP_LIST = []

    codon2aa = {}

    chr2start = {}
    chr2end = {}
    chr2pos = {}
    chr_pos2trna = {}
    chr_pos2region = {}
    chr_pos2rel = {}
    trna2note = {}
    trna2strand = {}

    ce2seq = {}
    ce2loc = {}    # chr, start, end
    loc2ce = {}
    ce2info = {}    # AA, codon
    ce2strand = {}

    loc2info = {}
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

    STATEMENT = "Dictionaries and lists generated."
    print(STATEMENT)
    LOG_LIST.append(STATEMENT + '\n')
    ## Get tRNA info.
    for line in open(CHAN_OUT_FILE):
        if line[0] == "S" or line[0] == "N" or line[0] == "-":
            continue
        seg = line.split("\t")
        STRAND = "+"
        CHR = str(seg[0]).strip()
        if "chr" in CHR or "Chr" in CHR: # get rid of unnecessary/unmatched 'Chrs'
            CHR = CHR.split('r')[1]
        NUM = str(seg[1])
        TRNA = CHR + ".trna" + NUM
        START = int(seg[2])
        END = int(seg[3])
        if START > END:
            START = int(seg[3])
            END = int(seg[2])
            STRAND = "-"
        AA = str(seg[4])
        # if AA == "Undet":
        #     continue
        CODON = str(seg[5])
        # if CODON == "NNN":
        #     continue
        NOTE = str(seg[-1]).strip('\n') # this should work with whatever version
        trna2note[TRNA] = NOTE
        trna2strand[TRNA] = STRAND # new, save this for output (instead of just use here)
        if CODON not in codon2aa:
            codon2aa[CODON] = AA
        ce2loc[TRNA] = [CHR,START,END]
        loc2ce[(CHR,START,END)] = TRNA
        ce2info[TRNA] = [AA,CODON]
        CHROM = CHR
        if CHROM not in chr2start:
            chr2start[CHROM] = []
            chr2end[CHROM] = []
        chr2start[CHROM].append(START)
        chr2end[CHROM].append(END)
        chr2pos.setdefault(CHROM,[])
        POS = (START - 40) # start with 40 bp before start. **Start is on chr, so smaller number REGARDLESS of strand
        while POS <= START: # look at each flank bp
            chr2pos[CHROM].append(POS)
            chr_pos2trna.setdefault((CHROM,POS),TRNA)
            if (START - POS) > 20:
                chr_pos2region.setdefault((CHROM,POS),"outer5") # 5' is to the L of gene if strand is +
                if STRAND == "-":
                    chr_pos2region[(CHROM,POS)] = "outer3" # 5' is to the R of gene if strand is -, so this to the L is 3'
            else:
                chr_pos2region.setdefault((CHROM,POS),"inner5")
                if STRAND == "-":
                    chr_pos2region[(CHROM,POS)] = "inner3"
            chr_pos2rel[(CHROM,POS)] = (POS - START)
            POS += 1
        POS = END
        while POS <= (END + 40): # start with END and go to 40 bp after **end is on chr, the bigger position REGARDLESS of strand
            chr2pos[CHROM].append(POS)
            chr_pos2trna.setdefault((CHROM,POS),TRNA)
            if (POS - END) > 20:
                chr_pos2region.setdefault((CHROM,POS),"outer3")
                if STRAND == "-":
                    chr_pos2region[(CHROM,POS)] = "outer5"
            else:
                chr_pos2region.setdefault((CHROM,POS),"inner3")
                if STRAND == "-":
                    chr_pos2region[(CHROM,POS)] = "inner5"
            chr_pos2rel[(CHROM,POS)] = (POS - END)
            POS += 1


    STATEMENT = "Reference tRNA information gathered."
    print(STATEMENT)
    LOG_LIST.append(STATEMENT + '\n')
    ## Assemble variant info.
    vcf_reader = vcf.Reader(open(VCF_FILE, 'r'))

    SAMPLES = vcf_reader.samples

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
            chr2start[CHR].sort()
            chr2end[CHR].sort()
            STATEMENT = "Chromosome " + CHR + " beginning..."
            print(STATEMENT)
            LOG_LIST.append(STATEMENT + '\n')
        if POS in chr2pos[CHROM]:
            loc2allele[(CHROM,POS)] = [REF,ALT] # based on VCF
            loc2hr[(CHROM,POS)] = []
            loc2ha[(CHROM,POS)] = []
            loc2het[(CHROM,POS)] = []
            loc2info[(CHROM,POS)] = []
            for sample in SAMPLES:
                if record.genotype(sample)['GT'] == "0/0":
                    loc2hr[(CHROM,POS)].append(sample)
                elif record.genotype(sample)['GT'] == "1/1":
                    loc2ha[(CHROM,POS)].append(sample)
                else:
                    loc2het[(CHROM,POS)].append(sample)
                    STATEMENT = "%s\t%d\t%s\t%s" % (CHROM,POS, sample, record.genotype(sample)['GT'])
                    LOG_LIST.append(STATEMENT + '\n')
            if "ANN" in record.INFO:
                for ann in record.INFO["ANN"]:
                    loc2info[(CHROM,POS)].append(ann)
        else:
            continue
        line = CHROM + "\t" + str(POS) + "\t"
        line += ("\t").join(loc2allele[(CHR,POS)]) + "\t" # alleles
        line += chr_pos2trna[(CHR,POS)] + "\t"  # tRNA
        line += trna2strand[chr_pos2trna[(CHR,POS)]] + "\t" # strand - NEW
        line += chr_pos2region[(CHR,POS)] + "\t" # region [flank]
        line += str(chr_pos2rel[(CHR,POS)]) + "\t" # relative position (flank)
        line += trna2note[chr_pos2trna[(CHR,POS)]] + "\t"
        # NEW: adding COUNTS of each genotype
        line += str(len(loc2hr[(CHR,POS)])) + "\t" # nHomRef
        line += str(len(loc2ha[(CHR,POS)])) + "\t" # nHomAlt
        line += str(len(loc2hr[(CHR,POS)]) + len(loc2ha[(CHR,POS)])) + "\t" # nNotMissingHet
        line += str(len(loc2het[(CHR,POS)])) + "\t" # nmissingOrHet
        # Back to old - add who has which genotype
        if len(loc2hr[(CHR,POS)]) > 0:
            line += (",").join(loc2hr[(CHR,POS)]) + "\t"
        else:
            line += "None" + "\t"
        if len(loc2ha[(CHR,POS)]) > 0:
            line += (",").join(loc2ha[(CHR,POS)]) + "\t"
        else:
            line += "None" + "\t"
        if len(loc2het[(CHR,POS)]) > 0:
            line += (",").join(loc2het[(CHR,POS)]) + "\t"
        else:
            line += "None" + "\t"
        loc2info.setdefault((CHR,POS),[])
        if len(loc2info[(CHR,POS)]) > 0:
            line += (",").join(loc2info[(CHR,POS)])
        else:
            line += "NA"
        write_list.append(line)



    STATEMENT = "Variant information assembled for strains."
    print(STATEMENT)
    LOG_LIST.append(STATEMENT + '\n')
    ## Print outfiles.
    # Variant file.
    outfn = args.strOut + "_strain_variants-flank.txt"
    outfile = open(outfn, 'w')

    PREAMBLE = "##File generated by [gitrepo]wormtrna/updatedinitial/get_strain_variants-flank.py.\n"
    outfile.write(PREAMBLE)

    HEADER = "##Chr\tPos\tRef\tAlt\ttRNA\tStrand\tRegion\trelPos\tNote\tnHomRef\tnHomAlt\tnNotMissingHet\tnmissingOrHet\thomRef\thomAlt\tmissingOrHet\tinfo\n"
    outfile.write(HEADER)

    for line in write_list:
        outfile.write(line + "\n")

    outfile.close()

    STATEMENT = "%s written." % (outfn)
    print(STATEMENT)
    LOG_LIST.append(STATEMENT + "\n")
    # Log file.
    outfn = args.strOut + "_get_strain_variants-flank.log"
    outfile = open(outfn, 'w')

    PREAMBLE = "##File generated by [gitrepo]wormtrna/updatedinitial/get_strain_variants-flank.py.\n"
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