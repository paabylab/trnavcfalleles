#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Get information on variants in tRNA genes from VCF - similar to get_strain_variants.py, 
BUT gets relative position of each variant within its tRNA (for 'easy' mapping to structure from output of secstruct2pieces.py)
Created on Thu Mar 20 09:53:44 2025

@author: Avery Davis Bell
"""

import argparse
import os
import sys
import time
import vcf # update - no 'vcf'-titled package but syntax looks like 'pyvcf' as far as I can tell
import numpy as np
import gzip
#%%
# FUNCTIONS FROM SECSTRUCT2SEQPIECES.PY - SLIGHTLY UPDATED IF NEEDED
def readss(ssfile):
    """
    Reads in each tRNA from ssfile in format to use for downstream analysis

    Parameters
    ----------
    ssfile : character
        Path to tRNAscan-SE .SS file for tRNAs to process here

    Returns
    -------
    list of the following - 
    t2loc : DICT 
        Key: tRNA name. 
        Value: to its location - stranded (CHR, start of gene, end of gene - start > end for - strand)
    t2locabs : DICT
        Key: tRNA name. 
        Value: its location where start always comes first on chromosome [CHR, start, end]
    loc2t : DICT
        Key: tRNA location (CHR, abs start [L on chr], abs end [R on chr])
        value: tRNA name
    t2strand : DICT
        Key: tRNA name
        Value: strand
    t2len : DICT
        Key: tRNA name 
        Value: length (from SS file - not calculated or updated here)
    t2codon : DICT
        Key: tRNA name 
        Value: codon info: [type/aa, anticodon, [rel start, rel stop], [genome start of gene, genome end of gene], [lowest genome pos of gene - L bound regardless of strand, R bound regardless of strand]]
    t2seq : DICT
        Key: tRNA name 
        Value: list of sequence, 2ndary structure **introns removed if there was an intron**
    t2intron : DICT
        Key:  tRNA name for any with intron to intron info (NOT for all)
        Value: [[start, end], [start on chr, end on chr]] 1) list intron relative position before it was removed from seq/struct 2) list ABS GENOMIC POSITION in smaller number, bigger number format!!
    tpseud : LIST
        list of tRNA names if flagged as possible pseudogenes, to save which tRNAs are flagged as possible pseudogenes just in case

    """

    # Initialize dictionaries, lists
    t2loc = {} # tRNA name to its location - stranded
    t2locabs = {} # tRNA name to its location where start always comes first on chromosome
    loc2t = {} # tRNA location to name
    t2strand = {} # tRNA name to strand
    t2len = {} # tRNA name to length
    t2codon = {} # tRNA name to codon info [relative position AND info about codon]
    t2seq = {} # tRNA name to sequence, 2ndary structure
    t2intron = {} # tRNA name for any with intron to intron info
    tpseud = []  # save which tRNAs are flagged as possible pseudogenes just in case
    
    # Process SS file
    iLine = 0
    strName = "" # tRNA name
    istm = open(ssfile)
    for strLine in istm:
        iLine += 1
        if "trna" in strLine: # this line starts this tRNA, save about tRNA info
            astrHead = strLine.strip().split()
            strName = astrHead[0]
            CHR = strName.split(".")[0]
            if CHR.startswith("chr"): # sometimes chr sneaks through in some files but not others
                CHR = CHR[3:]
            astrPos = astrHead[1].strip("(|)").split("-")
            aiPos = []
            for pos in astrPos:
                aiPos.append(int(pos))
            t2loc[strName] = [CHR, aiPos[0], aiPos[1]] # ** tRNA-wise not L-R [diff based on strand]
            if(aiPos[0] > aiPos[1]): # - stranded
                t2strand[strName] = "-"
                t2locabs[strName] = [CHR, aiPos[1], aiPos[0]]
            else:
                t2strand[strName] = "+"
                t2locabs[strName] = t2loc[strName]
            t2len[strName] = int(astrHead[3])
            loc2t[(t2locabs[strName][0], t2locabs[strName][1], t2locabs[strName][2])] = strName # location in L-R/absolute looks up tRNA
        elif strLine.startswith("Type"): # info including anticodon; save that
            astrLine = strLine.strip().split()
            strType = astrLine[1]
            strAnti = astrLine[3]
            astrRelPos = astrLine[5].split("-")
            astrAbsPos = astrLine[6].strip("(|)").split("-")
            aiAbsPos = []
            for pos in astrAbsPos:
                aiAbsPos.append(int(pos))
            t2codon[strName] = [strType, strAnti, [int(astrRelPos[0]), int(astrRelPos[1])], aiAbsPos] # type, codon, relative position, absolute position diff based on strand
            if t2strand[strName] == "-": # add absolute position L then R, regardless of strand
                t2codon[strName].append([aiAbsPos[1], aiAbsPos[0]])
            elif t2strand[strName] == "+":
                t2codon[strName].append(aiAbsPos)
        elif strLine.startswith("Possible intron"): # get intron info
            astrInInfo = strLine.strip().split(":")[1].split()
            astrRelPos = astrInInfo[0].split("-") # **1 indexed **
            t2intron[strName] = [[int(astrRelPos[0]), int(astrRelPos[1])]] # RELATIVE position
            # ABSOLUTE position: start + this if +; start [rel] - this if -. -1 for INCLUSIVE of starts/ends
            if t2strand[strName]=="+":
                t2intron[strName].append([t2loc[strName][1] + (int(astrRelPos[0])-1), t2loc[strName][1] + (int(astrRelPos[1])-1)])
            elif t2strand[strName]=="-":
                t2intron[strName].append([t2loc[strName][1] - (int(astrRelPos[1])-1), t2loc[strName][1] - (int(astrRelPos[0])-1)])
        elif strLine.startswith("Possible pseudogene"): # flag pseudogene
            tpseud.append(strName)
        elif strLine.startswith("Seq:"): # save seq
            oneSeq = strLine.strip().split()[1]
            t2seq[strName] = [oneSeq]
        elif strLine.startswith("Str:"): # save structure
            oneStruct = strLine.strip().split()[1]
            t2seq[strName].append(oneStruct)
        elif strLine.startswith("\n"): # this breaks between tRNAs
            continue
        
    istm.close()
    
    # Process & remove introns as needed [index: start-1:end]
    if len(t2intron) > 0:
        for strT, aaiPos in t2intron.items():
            aiPos = aaiPos[0]
            # Olds
            oldseq = t2seq[strT][0]
            oldstruct = t2seq[strT][1]
            # News
            newseq = oldseq[0:(aiPos[0]-1)] + oldseq[aiPos[1]:]
            newstruct = oldstruct[0:(aiPos[0]-1)] + oldstruct[aiPos[1]:]
            # Replace
            t2seq[strT][0] = newseq
            t2seq[strT][1] = newstruct
  
    
    # Return
    return([t2loc, t2locabs, loc2t, t2strand, t2len, t2codon, t2seq, t2intron, tpseud])

#%%
#%%    
def onepos2s(loc, cod, seq):
    """
    Maps each position in seq to its predicted tRNA segment based on secondary structure, codon info

    Parameters
    ----------
    loc : list
        [chromosome, int start pos [can be > end], int end pos]
    cod : list
        codon info: [type/aa, anticodon, [rel start, rel stop], [genome start of gene, genome end of gene], [lowest genome pos of gene - L bound regardless of strand, R bound regardless of strand]]
    seq : list
        list of sequence, 2ndary structure **introns removed if there was an intron**.
        Each value in the strings matches exactly

    Returns
    -------
    If tRNA is well structured: list of :
        hstr2pos : DICT
            Keys: all possible parts of tRNA structures [short format]. Values: Keep BOUNDS of structure, in pythonic indexing. Empty if that structure not observed
        hFlags : DICT
            empty if no flags (if tRNA was particularly well behaved/easy to characterize)
            keys: flag, values: list of description, number if relevant
        hstr2pos : DICT
            Most important, potentially! Keys are 1-indexed positions in tRNA (one for each)
                Values are list of [overall structure name, sub structure part, index this structure is from 1 in tRNA, number from one this base is within this structure part]
    If not enough structures possible after finding codon, just one dict:
         {"CodonIssue": ["Codon position too far over to leave space for other structures", cod[2]]}
    """
    # --- Break sec struct into chunks
    chunks = []
    chpos = []
    prev = ""
    i = 0
    for char in seq[1]:
        i+=1
        if char==prev:
            chunks[len(chunks)-1].append(char)
            chpos[len(chpos)-1].append(i)
        else:
            chunks.append([char])
            chpos.append([i])
        prev = char
    
    # Numbers to np.array
    for i in range(len(chpos)):
        chpos[i] = np.array(chpos[i])
    
    # --- Set up output, trackers
    hPosType = {} # keys will be int positions in tRNA; 
    #   values will be [STRUCTURE, SUBstructure, index of this substruct - as in hPossibleStructs,
    #   relative position in substructure
    
    hPossibleStructs = {"acc_stem_l": [1, "Req"], 
                        "acc_d_link": [2, "Opt"],
                        "d_arm_l": [3, "Req"],
                        "d_loop": [4, "Req"],
                        "d_arm_r": [5, "Req"],
                        "d_ant_link": [6, "Opt"],
                        "ant_arm_l": [7, "Req"],
                        "ant_loop_l": [8, "Req"],
                        "anticodon": [9, "Req"],
                        "ant_loop_r": [10, "Req"],
                        "ant_arm_r": [11, "Req"],
                        "variable_loop": [12, "Req"],
                        "t_arm_l": [13, "Req"],
                        "t_loop": [14, "Req"],
                        "t_arm_r": [15, "Req"],
                        "t_acc_link": [16, "Opt"],
                        "acc_stem_r": [17, "Req"],
                        "acc_overhang": [18, "Req"]
                        }
    aUsedStructs = [] # track which ones have been assigned and therefore shouldn't be used!
    
    # --- Anticodon [arm, actual antic, loop]
    antloop = 0
    antpos = list(range(cod[2][0], cod[2][1] + 1))
    antposch = []
    
    for i in range(len(chunks)):
        inch = 0
        for j in antpos:
            for k in range(len(chpos[i])):
                if j==chpos[i][k]:
                    inch+=1
                    antposch.append(k)
        if inch == len(antpos): # This chunk has anticodon, process it
            antloop = i # Know which one anticodon loop is
            # ANTICODON ITSELF
            ## Save info about it
            for j in range(len(antpos)):
                hPosType[antpos[j]] = ["anticodon", "anticodon", hPossibleStructs["anticodon"][0], j + 1] # 1 indexed
            ## Delete it from original data
            chpos[i] = np.delete(chpos[i], antposch)
            chunks[i] = np.delete(chunks[i], antposch)
            aUsedStructs.append(hPossibleStructs["anticodon"][0]) # anticodon itself
            # ANTICODON LOOP BITS
            aAL = []
            aAR = []
            for j in range(len(chpos[i])):
                if chpos[i][j]<antpos[0]: ## Right anticodon loop
                    hPosType[chpos[i][j]] = ["anticodon", "loop_L", hPossibleStructs["ant_loop_l"][0], j + 1]
                    aAL.append(j)
                if chpos[i][j]>antpos[len(antpos) - 1]: ## Left anticodon loop
                    hPosType[chpos[i][j]] = ["anticodon", "loop_R", hPossibleStructs["ant_loop_r"][0], j + 1]
                    aAR.append(j)
            ## Delete from orig data
            chpos[i] = np.delete(chpos[i], aAL + aAR)
            chunks[i] = np.delete(chunks[i], aAL + aAR)
            ## Note
            aUsedStructs.append(hPossibleStructs["ant_loop_l"][0])
            aUsedStructs.append(hPossibleStructs["ant_loop_r"][0])
            
            break # no need to keep going - done with this chunk, will start again to do next one
    if antloop==0: # EDGE case esp when try to force pseudos
        hFlags = {"CodonIssue": "Anticodon wasn't in one structure, often because I was trying to force something"}
        return([hFlags])
    
    # Do anticodon arms
    iant1 = antloop - 1
    iant2 = antloop + 1
    iAntMatch = 0
    bAnt = False # Flag
    antLen = 5 # default 
    if len(chunks[iant1])==len(chunks[iant2]) & len(chunks[iant1])>=(antLen - 2): # also needs to be a reasonable size
        for i in range(len(chunks[iant1])):
            if((chunks[iant1][i]==">" and chunks[iant2][i]=="<") or (chunks[iant1][i]=="<" and chunks[iant2][i]==">")):
                iAntMatch+=1
        if iAntMatch==len(chunks[iant1]): # Confirmed they're same length and look sensible, do 'em!
            bAnt = True # Flag this looks nice
            antLen = iAntMatch
    
    ## Grab the ones for the pair, whether pretty or not 
    hiAntL = {} # key will be index of list, values the ones IN THERE to grab
    hiAntR = {}
    iCtAntL = 0 # when gets to antLen, have grabbed enough
    iCtAntR = 0 # when gets to 5, have grabbed enough
    # Get them into hPosType; pop them from inputs
    ## Left
    for i in range(iant1, 0, -1): # L counts down/ goes away from anticodon loop
        for j in reversed(range(len(chpos[i]))):
            if i in hiAntL.keys():
                hiAntL[i].append(j)
            else:
                hiAntL[i] = [j]
            iCtAntL += 1
            if iCtAntL == antLen:
                break
        if iCtAntL==antLen:
            break
    ## Right
    for i in range(iant2, len(chpos)): # R arm counts up from anticodon loop
        for j in range(len(chpos[i])):
            if i in hiAntR.keys():
                hiAntR[i].append(j)
            else:
                hiAntR[i] = [j]
            iCtAntR += 1
            if iCtAntR == antLen:
                break
        if iCtAntR==antLen:
            break
    # Get them into hPosType; pop them from inputs
    myI = antLen + 1
    for i, aiJ in hiAntL.items():
        for j in aiJ:
            myI-=1
            hPosType[chpos[i][j]] = ["anticodon", "arm_L", hPossibleStructs["ant_arm_l"][0], myI]
        chpos[i] = np.delete(chpos[i], aiJ)
        chunks[i] = np.delete(chunks[i], aiJ)
    myI = 0 # want to count UP here too
    for i, aiJ in hiAntR.items():
        for j in aiJ:
            myI += 1
            hPosType[chpos[i][j]] = ["anticodon", "arm_R", hPossibleStructs["ant_arm_r"][0], myI]
        chpos[i] = np.delete(chpos[i], aiJ)
        chunks[i] = np.delete(chunks[i], aiJ)
    aUsedStructs.append(hPossibleStructs["ant_arm_l"][0]) # document
    aUsedStructs.append(hPossibleStructs["ant_arm_r"][0]) # document    

    # --- Deal with SPECIAL CASE: anticodon is so far over there's not room for other structures
    if max(hiAntR.keys()) > (len(chunks) - 5):
        return([{"CodonIssue": ["Codon position too far over to leave space for other structures", cod[2]]}]) # list makes it work with other things possibly returned

    # --- Acceptor stem [this could go at any point - doesn't depend on assignment of other structures]
    iASov = len(chpos) -1
    for j in range(len(chpos[iASov])):
        hPosType[chpos[iASov][j]] = ["acceptor_stem", "overhang", hPossibleStructs["acc_overhang"][0], j + 1]
    chpos[iASov] = np.array([])
    chunks[iASov] = []
    aUsedStructs.append(hPossibleStructs["acc_overhang"][0]) # document
    ## paired stem IF PRETTY: just get length
    iaL = 0
    iaR = iASov - 1
    bA = False
    accLen = 7 # default
    if len(chpos[iaL])<=len(chpos[iaR]) and len(chpos[iaL]) >= accLen:
        if (chunks[iaL][0]==">" and chunks[iaR][(len(chunks[iaR]) - accLen):len(chunks[iaR])][0]=="<") or (chunks[iaL][0]=="<" and chunks[iaR][(len(chunks[iaR]) - accLen):len(chunks[iaR])][0]==">"):
            accLen = len(chpos[iaL]) # only let it be LONGER!!
            bA = True
        # ***** NOTE: fixed - ACC STEM OFTEN RIGHT NEXT TO T, SO THIS IS OVERPENALIZING
        # ***NOTE: not fixed - this doesn't work as I'd like if there's matching stem but one mismatch

    ## paired stem IF NOT PRETTY: grab the first 7, last 7 that aren't overhang
    #           simplify idea: this actually should work if pretty too - if so, want to set length differently. DONE!
    hiAcL = {} # key will be index of list, values the ones IN THERE to grab
    hiAcR = {}
    iCtAcL = 0 # when gets to accLen, have grabbed enough
    iCtAcR = 0 # when gets to 7, have grabbed enough
    # Find which ones to get
    ## Left
    for i in range(len(chpos)):
        for j in range(len(chpos[i])):
            if i in hiAcL.keys():
                hiAcL[i].append(j)
            else:
                hiAcL[i] = [j]
            iCtAcL += 1
            if iCtAcL == accLen:
                break
        if iCtAcL==accLen:
            break
    ## Right
    for i in reversed(range(len(chpos))):
        for j in reversed(range(len(chpos[i]))):
            if i in hiAcR.keys():
                hiAcR[i].append(j)
            else:
                hiAcR[i] = [j]
            iCtAcR += 1
            if iCtAcR == accLen:
                break
        if iCtAcR==accLen:
            break
    # Get them into hPosType; pop them from inputs
    myI = 0
    for i, aiJ in hiAcL.items():
        for j in aiJ:
            myI+=1
            hPosType[chpos[i][j]] = ["acceptor_stem", "arm_L", hPossibleStructs["acc_stem_l"][0], myI]
        chpos[i] = np.delete(chpos[i], aiJ)
        chunks[i] = np.delete(chunks[i], aiJ)
    myI = accLen + 1 # want to count down here
    for i, aiJ in hiAcR.items():
        for j in aiJ:
            myI -= 1
            hPosType[chpos[i][j]] = ["acceptor_stem", "arm_R", hPossibleStructs["acc_stem_r"][0], myI]
        chpos[i] = np.delete(chpos[i], aiJ)
        chunks[i] = np.delete(chunks[i], aiJ)
    aUsedStructs.append(hPossibleStructs["acc_stem_l"][0]) # document
    aUsedStructs.append(hPossibleStructs["acc_stem_r"][0]) # document
            
    # --- D arm
    ## Bounded by 2 that are equal in length with opposite sec structure [in theory], and must be before anticodon
    ## Updated to grab 4 bases if a sec struct issue
    idL = 0
    idR = 0
    idLoop = 0
    bD = False # flag
    dLen = 4 # default
    for i in range(min(hiAntL.keys())): # This loop figures out if nicely structured; goes up to first involved in anticodon arm
        if i==0 or len(chunks[i])==0:
            continue # first is always acceptor, has to start with at least 2nd; it didn't like trying to get into an empty one
        iNotDot = 0
        for s in chunks[i]:
            if s!=".":
                iNotDot +=1
        if iNotDot == len(chunks[i]) and iNotDot==len(chunks[i + 2]) and len(chunks[i])>=(dLen - 2): # make sure not just length of 1 or whatever
            if((chunks[i][0]==">" and chunks[i + 2][0]=="<") or (chunks[i][0]=="<" and chunks[i + 2][0]==">")):
                dLen = max([iNotDot, dLen]) # if there's 2 that match and 1 that doesn't, this WILL mess up. trying this way...
            # check structs are opposite, then assign idL as i, idR as i+2, idLoop as i+1, and break
                idL = i
                idR = i+2
                idLoop = i+1
                bD = True # flag we're good to assign this [do this elsewhere!!]
                break
    if not bD:
       # Figure out where loop is based on anticodon
       for i in range(min(hiAntL.keys()), 1, -1): # expect it to be close
           iNotDot = 0
           for s in chunks[i]:
               if s!=".":
                   iNotDot+=1
           if iNotDot >=1: # some pairing here, going to assume this is the ARM
               iDot = 0
               for j in range(i - 1 , 1, -1): # loop through to find all .s (enough to be a loop)
                   for s in chunks[j]:
                       if s==".":
                           iDot += 1
                   if iDot >= 5: # d loop needs to be big
                       idLoop = j
                       break
               if iDot>=5:
                   break # out of outer loop
    ## Grab loop - should know where it is regardless now
    if idLoop==0:
        # print("Big issue, don't have d loop, somethings gonna be weird, NOT FIXED " + loc[0] + ":" + str(loc[1]) + "-" + str(loc[2]))
        hFlags = {"D_loop" : ["Unable to identify D loop, VERY WEIRD & UNTESTED", ""]}
        if cod[0]=="Fake": # one way this can happen is when I try to force a pseudogene
            hFlags["CodonIssue"] = ["Codon position not defined - not processed despite trying & faking it", ""]
        return([hFlags])
    # ** DEAL WITH WEIRD CASE WHERE PAIRING ISN'T PREDICTED ***: CAN'T TAKE THE WHOLE LOOP FIRST IN AT LEAST ONE CASE
    #       need to find the arms before getting rid of the loop ##
    # doing that below
    else:
        for j in range(len(chpos[idLoop])):
            hPosType[chpos[idLoop][j]] = ["d", "loop", hPossibleStructs["d_loop"][0], j + 1]
        chpos[idLoop] = np.array([])
        chunks[idLoop] = []
        aUsedStructs.append(hPossibleStructs["d_loop"][0]) # document
    ## do paired stem
    # # Hash to grab which ones are which; update to deal with LOOP!!
    hiDL = {} # key will be index of list, values the ones IN THERE to grab
    hiDR = {}
    iCtDL = 0 # when gets to dLen, have grabbed enough
    iCtDR = 0 # when gets to dLen, have grabbed enough
    # Find which ones to get
    ## Left
    ## Deal with WEIRD CASE where there isn't an arm that can be found
    chEmpty = []
    bDArmIssue = False
    for i in range(len(chunks)):
        if len(chunks[i])==0:
            chEmpty.append(i)
    if (idLoop - 1) in chEmpty:
        bDArmIssue = True
        hiDL = {idLoop: []} # Make sure note d loop in hiDL instead
    ## keep going if not
    if not bDArmIssue:
        for i in range(idLoop -1, 1, -1):
            for j in reversed(range(len(chpos[i]))):
                if i in hiDL.keys():
                    hiDL[i].append(j)
                else:
                    hiDL[i] = [j]
                iCtDL += 1
                if iCtDL == dLen:
                    break
            if iCtDL==dLen:
                break
    ## Right: updated for if first is not paired ,need to go into loop
    for i in range(idLoop + 1, min(hiAntL.keys()) + 1): # better to be short than to go beyond anticodon!!
        for j in range(len(chpos[i])):
            if i in hiDR.keys():
                hiDR[i].append(j)
            else:
                hiDR[i] = [j]
            iCtDR += 1
            if iCtDR == dLen:
                break
        if iCtDR==dLen:
            break
    # Get them into hPosType; pop them from inputs
    myI = dLen + 1
    if not bDArmIssue:
        for i, aiJ in hiDL.items():
            for j in aiJ:
                myI-=1
                hPosType[chpos[i][j]] = ["d", "arm_L", hPossibleStructs["d_arm_l"][0], myI]
            chpos[i] = np.delete(chpos[i], aiJ)
            chunks[i] = np.delete(chunks[i], aiJ)
        aUsedStructs.append(hPossibleStructs["d_arm_l"][0]) # document
    myI = 0
    for i, aiJ in hiDR.items():
        for j in aiJ:
            myI += 1
            hPosType[chpos[i][j]] = ["d", "arm_R", hPossibleStructs["d_arm_r"][0], myI]
        chpos[i] = np.delete(chpos[i], aiJ)
        chunks[i] = np.delete(chunks[i], aiJ)
    aUsedStructs.append(hPossibleStructs["d_arm_r"][0]) # document

    # --- T arm
    ## do based on pos from end/acceptor arm
    itL = 0
    itR = 0
    itLoop = 0
    bT = False # flag
    tLen = 4 # default
    ## Figure out where they are if perfectly structured
    for i in range(min(hiAcR.keys()), min(hiAcR.keys()) - 3, -1): # It needs to either be RIGHT NEXT to or JUST ONE OVER From acc stem!
        if len(chunks[i])==0:
            continue
        iNotDot = 0
        for s in chunks[i]:
            if s!=".":
                iNotDot += 1
        if iNotDot == len(chunks[i]) and (iNotDot==len(chunks[i - 2])) and len(chunks[i]) >= (tLen - 1): # or (iNotDot + 1)==len(chunks[i-2])# Sometimes off by one here, need to allow for that - want to find FIRST pair
            if((chunks[i][0]==">" and chunks[i - 2][0]=="<") or (chunks[i][0]=="<" and chunks[i - 2][0]==">")):
                bT = True
                tLen = max(iNotDot, tLen) ## minimum in case of issues...fix?
                itL = i -2
                itR = i
                itLoop = i -1
                break
    ## Figure out where loop is if not nicely structured
    if not bT:
        for i in range(min(hiAcR.keys()), max(hiAntR.keys()) + 1, -1):
            if len(chunks[i])==0:
                continue
            iNotDot = 0
            for s in chunks[i]:
                if s!=".":
                    iNotDot += 1
            if iNotDot >=1: # some pairing here, first walking back from acc stem is t arm
                iDot = 0
                for j in range(i - 1, 1, -1 ): # loop through to find all .s (enough to be a loop)
                    for s in chunks[j]:
                        if s==".":
                            iDot +=1
                    if iDot>= 3: # t loop not sure size basis
                        itLoop = j
                        break
                if iDot>=3:
                    break # out of outer loop
    ### Grab loop - should know where it is regardless now
    if itLoop==0:
        # print("Big issue, don't have t loop, somethings gonna be weird, NOT FIXED " + loc[0] + ":" + str(loc[1]) + "-" + str(loc[2]))
        hFlags = {"T_loop" : ["Unable to identify T loop, VERY WEIRD & UNTESTED", ""]}
        if cod[0]=="Fake": # one way this can happen is when I try to force a pseudogene
            hFlags["CodonIssue"] = ["Codon position not defined - not processed despite trying & faking it", ""]
        return([hFlags])
    else:
        for j in range(len(chpos[itLoop])):
            hPosType[chpos[itLoop][j]] = ["t", "loop", hPossibleStructs["t_loop"][0], j + 1]
        chpos[itLoop] = np.array([])
        chunks[itLoop] = []
        aUsedStructs.append(hPossibleStructs["t_loop"][0]) # document
    ## do paired stem
    # # Hash to grab which ones are which; update to deal with LOOP!!
    hiTL = {} # key will be index of list, values the ones IN THERE to grab
    hiTR = {}
    iCtTL = 0 # when gets to tLen, have grabbed enough
    iCtTR = 0 # when gets to 7, have grabbed enough
    # Find which ones to get
    ## Left
    for i in range(itLoop -1, 1, -1):
        for j in reversed(range(len(chpos[i]))):
            if i in hiTL.keys():
                hiTL[i].append(j)
            else:
                hiTL[i] = [j]
            iCtTL += 1
            if iCtTL == tLen:
                break
        if iCtTL==tLen:
            break
    ## Right
    for i in range(itLoop + 1, len(chpos)):
        for j in range(len(chpos[i])):
            if i in hiTR.keys():
                hiTR[i].append(j)
            else:
                hiTR[i] = [j]
            iCtTR += 1
            if iCtTR == tLen:
                break
        if iCtTR==tLen:
            break
    # Get them into hPosType; pop them from inputs
    myI = tLen + 1
    for i, aiJ in hiTL.items():
        for j in aiJ:
            myI-=1
            hPosType[chpos[i][j]] = ["t", "arm_L", hPossibleStructs["t_arm_l"][0], myI]
        chpos[i] = np.delete(chpos[i], aiJ)
        chunks[i] = np.delete(chunks[i], aiJ)
    myI = 0
    for i, aiJ in hiTR.items():
        for j in aiJ:
            myI += 1
            hPosType[chpos[i][j]] = ["t", "arm_R", hPossibleStructs["t_arm_r"][0], myI]
        chpos[i] = np.delete(chpos[i], aiJ)
        chunks[i] = np.delete(chunks[i], aiJ)
    aUsedStructs.append(hPossibleStructs["t_arm_l"][0]) # document
    aUsedStructs.append(hPossibleStructs["t_arm_r"][0]) # document

    # ---- Variable loop: ANYTHING REMAINING between anticodon arm & t arm
    hVar = {} # record what to grab/pop
    for i in range(max(hiAntR.keys()), min(hiTL.keys()) + 1):
        if(len(chpos[i]) > 0):
            for j in range(len(chpos[i])):
                if i in hVar.keys():
                    hVar[i].append(j)
                else:
                    hVar[i] = [j]
    # Get into hPosType; pop from inputs
    myI = 0
    for i, aiJ in hVar.items():
        for j in aiJ:
            myI+=1
            hPosType[chpos[i][j]] = ["variable_loop", "variable_loop", hPossibleStructs["variable_loop"][0], myI]
        chpos[i] = np.delete(chpos[i], aiJ)
        chunks[i] = np.delete(chunks[i], aiJ)
    aUsedStructs.append(hPossibleStructs["variable_loop"][0]) # document
    
    # --- Clean up: any remaining little connector bits, name them somehow sensibly. [these are optional bits]
    # Collect
    hLeftover = {}
    for i in range(len(chpos)):
        if(len(chpos[i])) > 0:
            for j in range(len(chpos[i])):
                if i in hLeftover.keys():
                    hLeftover[i].append(j)
                else:
                    hLeftover[i] = [j]
    # Sort
    hAccDLink = {}
    hDAntLink = {}
    hTAccLink = {}
    for i, aiJ in hLeftover.items():
        if(i in range(max(hiAcL), min(hiDL) + 1)): # acc_d_link
            hAccDLink[i] = aiJ
            aUsedStructs.append(hPossibleStructs["acc_d_link"][0]) # go ahead and document
        elif(i in range(max(hiDR), min(hiAntL) + 1)): # d_ant_link
            hDAntLink[i] = aiJ
            aUsedStructs.append(hPossibleStructs["d_ant_link"][0]) # go ahead and document
        elif(i in range(max(hiTR), min(hiAcR) + 1)): # t_acc_link
            hTAccLink[i] = aiJ
            aUsedStructs.append(hPossibleStructs["t_acc_link"][0]) # go ahead and document
    
    # Save, pop any needed
    hList = [hAccDLink, hDAntLink, hTAccLink]
    hDescrip = ["acceptor_d_link", "d_anticodon_link", "t_acceptor_link"]
    hShort = ["acc_d_link", "d_ant_link", "t_acc_link"]
    for k in range(len(hList)): # easier/faster than an if/else
        myI = 0    
        myh = hList[k]
        for i, aiJ in myh.items():
            for j in aiJ:
                myI += 1
                hPosType[chpos[i][j]] = [hDescrip[k], hDescrip[k], hPossibleStructs[hShort[k]][0], myI]
            chpos[i] = np.delete(chpos[i], aiJ)
            chunks[i] = np.delete(chunks[i], aiJ)
        aUsedStructs.append(hPossibleStructs[hShort[k]][0]) # document
                
    
    # --- Confirm all required elements have been ID'ed; no un-IDed positions remain (flag or break if not?)
    iReqMissing = 0 
    for struct, info in hPossibleStructs.items():
        if info[1]=="Req" and info[0] not in aUsedStructs:
            iReqMissing+=1
    iRemain = 0
    for i in range(len(chpos)):
        iRemain += len(chpos[i])
    
    # --- Add absolute genomic position? IDEALLY YES [fn]
    #           [what do with INTRONS though]
    #       Also, do not HAVE that info for diff/strain-spec alleles. So that makes sense to RE-INTEGRATE later [or, better yet, just have known variants in terms of rel pos in tRNA as well as genomic pos]
    
    # --- Save/manage flags
    hFlags = {} # Description, then number if relevant.
    if iReqMissing > 0:
        hFlags["ReqMissing"] = ["Some required structures missing", iReqMissing]
    if iRemain > 0:
        hFlags["Remaining"] = ["Not all tRNA bases assigned to a structure", iRemain]
    if not bAnt:
        hFlags["Anticodon"] = ["Anticodon arm stem not as expected", ""]
    if not bA:
        hFlags["Acceptor"] = ["Acceptor stem not as expected", ""]
    if not bD:
        hFlags["D_arm"] = ["D arm not as expected", ""]
    if idLoop==0:
        hFlags["D_loop"] = ["Unable to identify D loop, VERY WEIRD & UNTESTED", ""]
    if not bT:
        hFlags["T_arm"] = ["T arm not as expected", ""]
    if itLoop==0:
        hFlags["T_loop"] = ["Unable to identify T loop, VERY WEIRD & UNTESTED", ""]
    
    # --- Get boundaries for EACH structure [can combine later for, e.g., FASTA output]
    hstr2pos = {} # Keep BOUNDS of structure, in pythonic indexing. doing at end for simplicity of thinking [would've been faster to do on own]
    for struct, metinfo in hPossibleStructs.items():
        if struct not in hstr2pos.keys(): # initialize for each one, even if not observed in this tRNA
            hstr2pos[struct] = []
        aiThis = []
        for i1, info in hPosType.items():
            if info[2]==metinfo[0]: # match on structure number, names are different
                aiThis.append(i1 - 1) # i1 are 1 indexed, want output to be 0 indexed
        if(len(aiThis)>0):
            hstr2pos[struct] = [min(aiThis), max(aiThis)]
                
    # --- Return
    # return([chunks, chpos, hPosType, aUsedStructs, hFlags]) # TESTING
    return([hPosType, hFlags, hstr2pos])
            ### fix return to be what I actually want
            # Possibly: start-stop for each position [reverse way]. Make sure have ALL though. (could put NAs here when don't have them)
            #       meh, that's SUPER easy downstream
#%%
def pos2struct(t2len, t2loc, t2codon, t2seq, iMinLn = 71, bFakeCods = True):
    """
    Runs onepos2s for all tRNas after filtering out those that are too short or don't have codon info specified

    Parameters [all are outputs of readss unless otherwise noted]
    ----------
    t2len : DICT
        Key: tRNA name. 
        Value: to its length
    t2loc : DICT 
        Key: tRNA name. 
        Value: to its location - stranded (CHR, start of gene, end of gene - start > end for - strand)
    t2codon : DICT
        Key: tRNA name 
        Value: codon info: [type/aa, anticodon, [rel start, rel stop], [genome start of gene, genome end of gene], [lowest genome pos of gene - L bound regardless of strand, R bound regardless of strand]]
    t2seq : DICT
        Key: tRNA name 
        Value: list of sequence, 2ndary structure **introns removed if there was an intron**
    iMinLn : INT
        value below which tRNA processing fn won't be called
    bFakeCods : bool
        True or False, if codon position is 0, 0 / NNN pseudo-gene, fake it out as being 34-36 and try anyway [if gene long enough]

    Returns
    -------
    hash with tRNA names as keys; values are returns from ht2PosStruct OR a hash flagging why that wasn't/coudn't be run

    """
    ht2PosStruct = {}
    for tname in t2loc.keys():
        if t2len[tname] < iMinLn:
            ht2PosStruct[tname] = [{"TooShort": ["TRNA sequence TOO SHORT - not processed", t2len[tname]]}] # just flag
        elif t2codon[tname][2][0]== 0: # don't process - it can't find codon
            if bFakeCods:
                fakecod = ['Fake', 'Fake', [34, 36], [34, 36], [34, 36]]
                ht2PosStruct[tname] = onepos2s(t2loc[tname], fakecod, t2seq[tname])
            else:
                ht2PosStruct[tname] = [{"CodonIssue": ["Codon position not defined - not processed", ""]}] # just flag
        else:
            ht2PosStruct[tname] = onepos2s(t2loc[tname], t2codon[tname], t2seq[tname])
    return(ht2PosStruct)
        
#%%
#### new functions this script ###
def pullvcfvars(loc2t, strVCF):
    """
    Gets info from VCF on variants that overlap with positions in loc2t

    Parameters
    ----------
    loc2t : dict
        keys: chr, start, end [in absolute positions - not stranded], value: tRNA this is
    strVCF : str
        Path to VCF to process

    Returns
    -------
    dictionary of variant info:
        for each variant internal to a tRNA, key is (CHROM, POS) and value is list of:
            [[tRNA name this overlaps with] , [variant's alleles], [samples with 0/0 gts], [samples with 1/1 gts], [samples with missing or het gts]]

    """
    # Get all tRNAs for easy sorting
    chr2pos = {}
    for loc, tname in loc2t.items():
        chrom, start, end = loc
        if chrom not in chr2pos.keys():
            chr2pos[chrom] = []
        chr2pos[chrom].append([start, end])
    
    # Begin VCF processing
    vcf_reader = vcf.Reader(open(strVCF, 'r'))
    SAMPLES = vcf_reader.samples
    
    var2info = {}
    i = 0 # counter
    CHR = '' # start nowhere
    
    # Process each record
    for record in vcf_reader:
        CHROM = str(record.CHROM)
        POS = int(record.POS)
        REF = str(record.REF)
        ALT = str(record.ALT[0])
        if CHROM == "MtDNA":
            continue
        if CHROM.startswith("chr"): # sometimes chr sneaks through in some files but not others
            CHROM = CHROM[3:]
        if CHROM!=CHR:
            if CHROM not in chr2pos.keys(): # just keep going if don't have tRNAs on this chromosome [unlikely - but happens for mtDNAs!]
                continue
            CHR = CHROM
            chr2pos[CHR].sort(key = lambda x: x[0])
            i = 0 # keep track of # tRNAs processed
            print("Chromosome " + CHR + " beginning....")
        if i >= len(chr2pos[CHR]):
            continue
        if POS >= chr2pos[CHR][i][0] and POS<=chr2pos[CHR][i][1]:
            refs = []
            alts = []
            hetmiss = []
            for sample in SAMPLES:
                if record.genotype(sample)['GT'] =="0/0":
                    refs.append(sample)
                elif record.genotype(sample)['GT'] == "1/1":
                    alts.append(sample)
                else:
                    hetmiss.append(sample)
            var2info[(CHROM, POS)] = [loc2t[(CHR, chr2pos[CHR][i][0], chr2pos[CHR][i][1])], [REF,ALT], refs, alts, hetmiss]
        if POS > chr2pos[CHR][i][1]:
            i += 1
        if i >= len(chr2pos[CHR]):
            continue
        elif POS > chr2pos[CHR][i][1]:
            i += 1
        if i >= len(chr2pos[CHR]):
            continue
        elif POS > chr2pos[CHR][i][1]:
            i += 1
        if i >= len(chr2pos[CHR]):
            continue
        elif POS > chr2pos[CHR][i][1]:
            i += 1
        if i >= len(chr2pos[CHR][i]):
            continue
        elif POS > chr2pos[CHR][i][1]:
            i += 1
        # .... I can't imagine all of these are necessary but I'm not fixing it now
            
    # Return
    return(var2info)
        
#%%
def rep(x, n, mytype = "list"):
    if mytype=="str":
        out = ""
        for i in range(n):
            out = out + x
    elif mytype=="list":
        out = []
        for i in range(n):
            out.append(x)
    return(out)
#%%
def vartrnainfo(var2info, ht2PosStruct, t2intron, t2strand, t2locabs):
    """
    Gets tRNA information - relative position, what that is in secondary structure, etc - for each variant

    Parameters
    ----------
    var2info : dict
        output of pullvcfvars: for each variant internal to a tRNA, key is (CHROM, POS) and value is list of:
            [[tRNA name this overlaps with] , [variant's alleles], [samples with 0/0 gts], [samples with 1/1 gts], [samples with missing or het gts]]
    ht2PosStruct : dict
        output of pos2struct - see that for details
    t2intron : dict
        mapping of tRNAs with introns: keys are tRNAs IF they have introns, values are [[rel position in tRNA start, end], [abs position in genome start, end]]
    t2strand : dict
        keys: tRNA name, value: strand ('+' or '-')
    t2locabs : dict
        keys: tRNA name, value: [chr, start, end on CHR - start less than end regardless of strand]

    Returns
    -------
    dict: Keys are (chr, pos) of variant
        var2trna, variant to trna info: list of [tRNA name, tRNA strand, var's pos in coding tRNA seq [relative],
                                                 structure var is in, sub structure var is in, int structure number is in tRNA,
                                                 relative position of nt/this var in its substructure [for comparing across tRNAs]]
        

    """
    var2trna = {}
    
    for vpos, info in var2info.items():
        CHR, POS = vpos
        tname, als, refs, alts, hetmiss = info

        # tRNA position info
        iRelPos, strStr, strSubStr, iStrN, iStrRelPos = rep("NA", 5)
        if len(ht2PosStruct[tname])==3: # not just flaggs
        # Determine relative position of variant in tRNA (where possible)
            if t2strand[tname]=="+":
                iRelPos = POS - t2locabs[tname][1] + 1
            elif t2strand[tname]=="-":
                iRelPos = t2locabs[tname][2] - POS + 1
            if tname in t2intron.keys(): # deal with intron - adjust relative position if in coding, NA out if not
                if iRelPos in range(t2intron[tname][0][0], t2intron[tname][0][1] + 1):
                    iRelPos = "NA"
                else:
                    if iRelPos > t2intron[tname][0][0]: # if it's below, can stay as it is; if above (here), need to correct by intron size
                        iRelPos = iRelPos - (len(range(t2intron[tname][0][0], t2intron[tname][0][1])) + 1)
        # Determine info about this position
            if iRelPos!="NA":
                hPosType, hFlags, hstr2pos = ht2PosStruct[tname]
                strStr, strSubStr, iStrN, iStrRelPos = hPosType[iRelPos]
        
        # Add everything to output dicts
        var2trna[vpos] = [tname, t2strand[tname], iRelPos, strStr, strSubStr, iStrN, iStrRelPos]
    
    return(var2trna)
# %%
def writevars(var2info, var2trna, strOut):
    """
    Writes out one line per variant overlapping tRNA - to strOut

    Parameters
    ----------
    var2info : dict
        (output by pullvcfvars)
    var2trna : dict
        (output by vartrnainfo)
    strOut : str
        Path to output file. If .gz, will be gzipped

    Returns
    -------
    None.
    
    Side effects
    ---------
    file written to strOut. Columns:
        chrom: chromosome 
        pos: position of **VCF variant**
        ref: ref allele of variant - as in VCF
        alt: alt allele of variant - as in VCF
        nHomRef: # samples in VCF with 0/0 gts
        nHomAlt: # samples in VCF with 1/1 gts
        nNotMissingHet: # samples in VCF with 0/0 or 1/1 gts (not 0/1, missing)
        tRNA: tRNA gene name
        tRNA_strand: strand of tRNA
        tRNA_pos: position of variant **relative to tRNA non intronic bases
        structure: structure of tRNA this base is in
        substructure: substructure of tRNA this base is in
        nSubStructure: count substructure is from L to R (of all substructures)
        substructure_pos: position of variant **within its substructure
        nMissingHet: # samples in VCF with 0/1 or missing gts
        homRef: comma-separated samples with 0/0 gts
        homAlt: comma-separated samples with 1/1 gts
        missingOrHettRNA: comma-separated samples with 1/1 or ./. gts

    """
    
    # Info to write
    PREAMBLE = "##File generated by [gitrepo]wormtrna get_strain_variants_relpos.py at %s \n" % time.ctime()
    HEADER = "\t".join(["chrom", "pos", "ref", "alt", "nHomRef", "nHomAlt", "nNotMissingHet", "nMissingHet",
                        "tRNA", "tRNA_strand", "tRNA_pos", "structure", "substructure", "nSubStructure", "substructure_pos",
                        "homRef", "homAlt", "missingOrHet"]) + "\n"
    wlines = [PREAMBLE, HEADER]
    for vpos, info in var2info.items():
        tname, als, refs, alts, hetmiss = info
        astrtinfo = []
        for el in var2trna[vpos]:
            astrtinfo.append(str(el))
        aline = [str(vpos[0]), str(vpos[1])] + als + [str(len(refs)), str(len(alts)), str(len(refs) + len(alts)), str(len(hetmiss))] + \
            astrtinfo + [",".join(refs), ",".join(alts), ",".join(hetmiss)] 
        wlines.append("\t".join(aline) + "\n")

    # Open output file
    if strOut.endswith("gz"):
        ostm = gzip.open(strOut, mode = 'wb')
        for line in wlines:
            ostm.write(line.encode("utf-8"))
    else:
        ostm = open(strOut, mode = "w")
        for line in wlines:
            ostm.write(line)
    
    # close
    ostm.close()

#%%
def main():
    #### Parse arguments
    argp = argparse.ArgumentParser(prog = "get_strain_variants_relpos.py",
                                   description = "Get information on variants in tRNA genes from VCF - similar to get_strain_variants.py, \
                                   BUT gets relative position of each variant within its tRNA (for 'easy' mapping to structure from output of secstruct2pieces.py)")
    # I/O
    argp.add_argument("-trnass", dest = "strSS", required = True, type = str,
                      metavar = "trnas.SS",
                      help = "Path to tRNAscan-SE .SS file for tRNAs to process here (reference tRNAs). Should have TRUE genomic coordinates to match those in VCF")
    argp.add_argument("-vcf", dest = "strVCF", required = True, type = str,
                      metavar = "variants.vcf",
                      help = "Path to VCF file for species whos tRNAs are provided")
    argp.add_argument("-out", dest = "strOut", required = True, type = str,
                      metavar = "outfile.txt.gz",
                      help = "Path to output file. If .gz, will be gzipped")
    # Related to 2ary structure processing
    argp.add_argument("-forcepseud", dest = "bFakeCods", required = False, type = bool,
                      default = True,
                      metavar = "[True, False]",
                      help = "During secondary structure processing: If gene is pseudo with codon listed as NNN as 0-0, try to force processing by trying codon at 34-36. (Default : True)")
    # Parse
    args = argp.parse_args()
    
    #### Set up logging
    statement = "....Starting at %s..." % time.ctime()
    print(statement)
    log_list = []
    log_list.append(statement + '\n')
    
    #### Read tRNA info
    statement = "....Reading in tRNA info (at %s)..." % time.ctime()
    print(statement)
    log_list.append(statement)
    t2loc, t2locabs, loc2t, t2strand, t2len, t2codon, t2seq, t2intron, tpseud = readss(args.strSS)
    
    #### Pull variants that are within tRNAs
    statement = "....Reading in variant info (at %s)..." % time.ctime()
    print(statement)
    log_list.append(statement + "\n")
    var2info = pullvcfvars(loc2t, args.strVCF)
    # log how many read in
    statement = "...Read in " + str(len(var2info)) + " variants that overlapped with provided tRNAs"
    print(statement)
    log_list.append(statement + "\n")
    
    #### Get tRNA secondary structure info
    statement = "...Mapping tRNA positions to tRNA structure for variant interpretation (at %s)..." % time.ctime()
    print(statement)
    log_list.append(statement + '\n')
    ht2PosStruct = pos2struct(t2len, t2loc, t2codon, t2seq, iMinLn = 71, bFakeCods = args.bFakeCods)
    
    
    #### Get relative position in tRNA of variants
    statement = "...Curating tRNA information such as relative position for each variant (at %s)..." % time.ctime()
    print(statement)
    log_list.append(statement + '\n')
    var2trna = vartrnainfo(var2info, ht2PosStruct, t2intron, t2strand, t2locabs)
    
    #### Save out
    statement = "...Writing tRNA-overlapping variant info (at %s)..." % time.ctime()
    print(statement)
    log_list.append(statement + '\n')
    writevars(var2info, var2trna, args.strOut)
    
    #### Done! Log
    statement = "Processing done! Writing log file and wrapping up (at %s)..." % time.ctime()
    print(statement)
    log_list.append(statement + '\n')
    logfilepath = "get_strain_variants_relpos-" + time.strftime("%Y%m%d-%H%M%S") + ".log"
    logfile = open(logfilepath, "w")
    PREAMBLE = "##File generated by [gitrepo]wormtrna sget_strain_variants_relpos.py at %s.\n" % time.ctime()
    logfile.write(PREAMBLE)
    for line in log_list:
        logfile.write(line)
    logfile.close()

if __name__=="__main__":
    main()
