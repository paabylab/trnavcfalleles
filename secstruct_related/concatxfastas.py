#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Worker script to concatenate fastas based on seq ID
Like seqkit concat BUT can do xfastas/ concatenates each line other than ID separately
Plus some tRNA specific stuff

Created on Wed Mar 19 13:42:22 2025

@author: Avery Davis Bell
"""
import argparse
import gzip
import time
import sys
import glob
#%%
def readseqs(strPath, bUReplace):
    """
    Read in one fastx [or fasta] file as dict of seq names

    Parameters
    ----------
    strPath : str
        path to file to read
    bUReplace : bool
        True or False, replace any Us in seq with Ts

    Returns
    -------
    dictionary, keys are seq IDs, values are list: each line of sequence or similar after >seq_ID line

    """
    hSeqs = {}
    iSeqs = 0
    iLine = 0
    
    istm = open(strPath, "r")
    for strLine in istm:
        iLine += 1
        strLine = strLine.strip()
        if strLine[0]==">":
            strID = strLine[1:]
            iSeqs+=1
            iID = iLine
            hSeqs[strID] = []
        if iLine==(iID + 1):
            if bUReplace:
                newSeq = ""
                for let in strLine:
                    if let=="U":
                        newSeq = newSeq + "T"
                    else:
                        newSeq = newSeq + let
            else:
                newSeq = strLine
            hSeqs[strID].append(newSeq)
        elif iLine > iID + 1:
            hSeqs[strID].append(strLine)
            
    istm.close()
    
    print("...." + str(iSeqs) + " seqs read from " + strPath + "\n")
    return(hSeqs)
#%%
def replaceAnticodon(hToReplace, strCodFile, bUReplace):
    """
    Replaces Ns in first element of each aaToReplace with bases from strCodFile

    Parameters
    ----------
    hToReplace : dict
        dictionary with sequence IDs as keys, list of [sequence, any other lines from fasta] as value
    strCodFile : str
        path to fasta file that has the same sequence IDs as in hToReplace
    bUReplace: bool
        replace Us in anticodon seqs. Passed to readseqs

    Returns
    -------
    hToReplace but with values modified to have the sequences from strCodFile

    """
    hCods = readseqs(strCodFile, bUReplace)
    
    for id, aSeqs in hToReplace.items():
        if id in hCods.keys():
            iN = 0
            strOut = ""
            for let in aSeqs[0]:
                if let=="N":
                    newlet = hCods[id][0][iN]
                    iN+=1
                else:
                    newlet = let
                strOut = strOut + newlet
            hToReplace[id][0] = strOut
            
    return(hToReplace)

#%%
def rep(x, n):
    out = ""
    for i in range(n):
        out = out + x
    return(out)

#%%
def concatwrite(astrOrder, hhSeqs, bFull, strFill, strOut):
    """
    Combines seqs (and other lines) based on seq ID from multiple fastas read in already as dict
    Writes these out

    Parameters
    ----------
    astrOrder : list of strs
        ORDER sequences should be combined in - keys of hhSeqs!
    hhSeqs : dict of dicts
        seqs to combine. outer keys are sequence order as in astrOrder
        Inner keys are sequence names, values are lists of sequences to combine/write
    bFull : bool
        True or False: keep all sequences, like full/outer join [even if they don't show up in a file or two']
    strFill : str
        ill with N bases/residues for IDs missing in some files when using -full
    strOut : str
        Filepath for output. If ends with .gz, output will be gzipped

    Returns
    -------
    None.

    """
    # Meta info about length of input seqs - used if bFull and for below
    hLenFill = {}
    for fname, hSeq in hhSeqs.items():
      hLenFill[fname] = [len(hSeq[list(hSeq.keys())[0]][0]), len(hSeq[list(hSeq.keys())[0]])] # length of seq, number of seqs to write

    # save Meta info about length of input seqs - just in case
    ostminf = open(strOut + "_seq_length_info.txt", "w")
    HEADER = ["input_file_name", "sequence_length", "number_seq_lines"]
    ostminf.write("\t".join(HEADER) + "\n")
    for fname, aiLen in hLenFill.items():
        wline = [fname, str(aiLen[0]), str(aiLen[1])]
        ostminf.write("\t".join(wline) + "\n")
    ostminf.close()
    
    # Get sequences to combine based on bFull
    hIDs = {}
    for fname, hSeq in hhSeqs.items():
        for id, seqs in hSeq.items():
            if id in hIDs.keys():
                hIDs[id] += 1
            else:
                hIDs[id] = 1
    if bFull: # use all IDs regardless; figure out length fill in for each
        useIDs = list(hIDs.keys())
    else: # only use those in all hhSeqs
        useIDs = []
        for id, iN in hIDs:
            if iN==len(hhSeqs.keys()):
                useIDs.append(id)
    
    # Combine & write out
    if strOut.endswith("gz"):
        ostm = gzip.open(strOut, mode = 'wb')
    else:
        ostm = open(strOut, mode = "w")
    for id in useIDs:
        wrID = ">" + id
        astrWSeqs = [] 
        for i in range(hLenFill[list(hLenFill.keys())[0]][1]):
            astrWSeqs.append("")
        for strOrd in astrOrder:
            hSeq = hhSeqs[strOrd]
            if id in hSeq.keys(): # always true if not bFull; normal case
                for i in range(len(hSeq[id])):
                    astrWSeqs[i] = astrWSeqs[i] + hSeq[id][i]
                # NEed it to be length of as many fasta lines, but combine within those....
            else: # make fake out seqs for this one [fill in]
                for i in range(hLenFill[strOrd][1]):
                    astrWSeqs[i] = astrWSeqs[i] + rep(strFill, hLenFill[strOrd][0])
        wLines = [wrID + "\n"]
        for strW in astrWSeqs:
            wLines.append(strW + "\n")
        if strOut.endswith("gz"):
            for line in wLines:
                ostm.write(line.encode("utf-8"))
        else:
            for line in wLines:
                ostm.write(line)
                
    # close out
    ostm.close()

#%%
def main():
    #### Parse arguments
    argp = argparse.ArgumentParser(prog = "concatxfastas.py",
                                   description = "Worker script to concatenate fastas based on seq ID \
                                   Like seqkit concat BUT can do xfastas/ concatenates each line other than ID separately \
                                   Plus some tRNA specific stuff")
    argp.add_argument("-files", dest = "strFiles", required = True, type = str,
                      metavar = "file1,file2,file3",
                      help = "Comma-separated list of files [plain text] to concatenate by seq ID. \
                          They will be combined IN THIS ORDER. Globs fine if they're not in sym links!")
    argp.add_argument("-out", dest = "strOut", required = True, type = str,
                      metavar = "out.xfasta.gz",
                      help = "Output filepath. If ends in .gz, will be gzipped")
    argp.add_argument("-full", dest = "bFull", required = False, type = bool, default = True,
                      metavar = "[True, False]",
                      help = "True or False: keep all sequences, like full/outer join")
    argp.add_argument("-fill", dest = "strFill", required = False, type = str, default = "N",
                      metavar = "N",
                      help = "fill with N bases/residues for IDs missing in some files when using -full")
    argp.add_argument("-replaceU", dest = "bU", required = False, type = bool, default = True,
                      metavar = "[True, False]",
                      help = "True or False: replace any Us in seq with Ts")
    argp.add_argument("-replaceCodon", dest = "bCod", required = False, type = bool, default = True,
                      metavar = "[True, False]",
                      help = "Replace anticodons [in 'antarm' file!!] that are N'ed out with their sequence. If so, must provide codon fasta too")
    argp.add_argument("-codonfasta", dest = "strCodFile", required = False, type = str,
                      metavar = "anticodon.xfasta",
                      help = "If -replaceCodon is True, path to fasta file with codon seq and only codon seq. Will replace Ns in 'antarm' file")
    args = argp.parse_args()
    
    
    #### Set up logging
    statement = "Starting at %s..." % time.ctime()
    print(statement)
    log_list = []
    log_list.append(statement + '\n')
    
    #### Read in all seqs
    # Parse files
    astrFilesStart = args.strFiles.strip().split(",")
    astrFiles = []
    for strFile in astrFilesStart:
        cleanfile = glob.glob(strFile)[0] # DANGER WILL ROBINSON
        astrFiles.append(cleanfile)
    # Log
    statement = "Reading in sequences (at %s)..." % time.ctime()
    print(statement)
    log_list.append(statement)
    # Do the work
    hhSeqs = {}
    for strFile in astrFiles:
        hhSeqs[strFile] = readseqs(strFile, args.bU)

    #### Replace codon if desired
    if args.bCod:
        # log
        statement = "Replacing N codons in antarm fasta (at %s)..." % time.ctime()
        print(statement)
        log_list.append(statement)
        # Get full file name as needed
        strAnticF = glob.glob(args.strCodFile)[0]
        # Do replacement
        antarm = ""
        for strFile in astrFiles: # find which one
            if "antarm" in strFile:
                antarm = strFile
        hhSeqs[antarm] = replaceAnticodon(hhSeqs[antarm], strAnticF, args.bU)
    
    #### Combine and save out
    # Log
    statement = "Combining sequences (at %s)..." % time.ctime()
    print(statement)
    log_list.append(statement)
    # Do the work
    concatwrite(astrOrder = astrFiles, hhSeqs = hhSeqs, bFull = args.bFull, strFill = args.strFill, strOut = args.strOut)
    
    #### Done! Log
    statement = "Processing done! Writing log file and wrapping up (at %s)..." % time.ctime()
    print(statement)
    log_list.append(statement + '\n')
    logfilepath = "concatxfastas-" + time.strftime("%Y%m%d-%H%M%S") + ".log"
    logfile = open(logfilepath, "w")
    PREAMBLE = "##File generated by [gitrepo]wormtrna concatxfastas.py at %s.\n" % time.ctime()
    logfile.write(PREAMBLE)
    for line in log_list:
        logfile.write(line)
    logfile.close()

if __name__=="__main__":
    main()
