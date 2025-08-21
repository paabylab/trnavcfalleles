#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Break tRNAs into their structures (t arm etc) \
    based on tRNAScan-SE secondary structure. For alignment.
    
Created on Wed Mar  5 09:53:41 2025

@author: abell65
"""
## imports
import argparse
import os
import subprocess
import time
# import vcf # update - no 'vcf' package but syntax looks like 'pyvcf' as far as I can tell
import re
import numpy as np
import gzip

#%%
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
        Value: [start, end] intron relative position before it was removed from seq/struct
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
            t2intron[strName] = [int(astrRelPos[0]), int(astrRelPos[1])]
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
        for strT, aiPos in t2intron.items():
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
    elif antloop<5: # EDGE case when try to force pseudos - it's really long and anticodon not found
        hFlags = {"CodonIssue": "Anticodon wasn't far enough into structure, often because I was trying to force something"}
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
    
def saveinfo(t2seq, t2intron, tpseud, ht2PosStruct, strOStem):
    """
    Writes 3 output files: 
        *_tRNA2struct_info.txt.gz: Information on tRNAs, one row per tRNA. NB doesn't include global stuff like position, strand as this doesn't make sense for strain-specific sequences that might be provided [has to come from elsewhere]
          Columns (DESCRIPTIONS of them can be found in summary file):
            tRNA: gene name
            PossiblePseudogene: T or F
            Intron: T or F
            TooShort: T or F
            CodonIssue: T or F
            ReqStructuresMissing: number of required tRNA cloverleaf structures missing in this one (should be 0) 
            RemainingBases: number of bases in tRNA not assigned to a structure (should be 0)
            AnticodonArmIssue: T or F. Anticodon arm stem not as expected (could be just one unpaired base, or worse)
            AccStemIssue: T or F. Acceptor stem not as expected (could be just one unpaired base, or worse)
            DArmIssue: T or F. D arm not as expected (could be just one unpaired base, or worse)
            DLoopIssue: T or F. Unable to identify D loop, VERY WEIRD & UNTESTED
            TArmIssue: T or F. T arm not as expected (could be just one unpaired base, or worse)
            TLoopIssue: T or F. Unable to identify T loop, VERY WEIRD & UNTESTED
        *_tRNA2struct_summary.txt: summary of number of tRNAs with various characteristics. Columns:
            Category: flag category (intron, codon issue, etc)
            Description: description of category (as defined in this script)
            N: number of tRNAs flagged for this category
        *_tRNA2struct.txt.gz: For each base in tRNA, what structure does it map to. 
                            For all tRNAs where this was doable (NOT necessarily all tRNAs), one line per base in non-intronic length of that tRNA [1 indexed].
          Columns:
            tRNA: gene name
            pos.tRNA: position in the tRNA (1-length of tRNA) this row describes
            Structure: overall structure this position is assigned to (e.g. 'D' for D arm)
            SubStructure: substructure this position is assigned to (e.g. 'arm_L' for L/first part of paired arm before D loop)
            num.structure: number this substructure is in the tRNA (with acceptor stem L starting at 1)
            pos.substructure: position this base is within its substructure [for comparing across tRNA]
            nucleotide: nt at this structure (from sequence in input SS file)
            secstruct: sec structure in . >< format at this structure (from sec structure in input SS file)

    Parameters
    ----------
    t2seq : DICT
        Key: tRNA name 
        Value: list of sequence, 2ndary structure **introns removed if there was an intron**
    t2intron : DICT
        Key:  tRNA name for any with intron to intron info (NOT for all)
        Value: [start, end] intron relative position before it was removed from seq/struct
    tpseud : LIST
        list of tRNA names if flagged as possible pseudogenes, to save which tRNAs are flagged as possible pseudogenes just in case
    ht2PosStruct : DICT
        output of pos2struct; see that for details. hash with tRNA names as keys; values are returns from ht2PosStruct OR a hash flagging why that wasn't/coudn't be run
    strOStem : str
        prefix for output filepath (including any path info).

    Returns
    -------
    None.

    """
    # --- Initialize summary trackers
    hInfo = {} # for tRNA info outfile. key is name
        # Structure: TooShort (T/F); CodonIssue (F or description); ReqMissing (# of req structures missing - usually 0);
        # RemainingBases (# bases unassigned to structure - usually 0), AnticodonArmIssue anticodon arm issue (T/F), acc stem issue (T/F),
        # DArmIssue d arm not as expected T/F, DLooIssuep couldn't ID D loop T/F, TArmIssue t arm not as expected T/F, TLoopIssue couldn't ID t Loop T/F        
    iNoFlags, iTooShort, iCodonIssue, iReqMissing, iRemaining, iAntArm, iAcc, iDArm, iDLoop, iTArm, iTLoop = [0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0]
    
    # --- Do major output along with tracking
    datoutf = strOStem + "_tRNA2struct.txt.gz"
    ostm = gzip.open(datoutf, mode = 'wb')
    
    ## Preamble
    PREAMBLE = "##File generated by [gitrepo]wormtrna secstruct2seqpieces.py on %s.\n" % time.ctime()
    ostm.write(PREAMBLE.encode("utf-8"))
    ## Header
    HEADER = "\t".join(["tRNA", "pos.tRNA", "Structure", "SubStructure", "num.structure", "pos.substructure", "nucleotide", "secstruct"]) + "\n"
    ostm.write(HEADER.encode("utf-8"))
    
    for tname, info in ht2PosStruct.items():
        # Set up
        if len(info)==3:
            hPosType, hFlags, hstr2pos = info
        elif len(info)==1:
            hFlags = info[0]
            hPosType = False # just for logic flow later
        # Identify any flags (& continue if needed)
        infNoFlags, infTooShort, infCodonIssue, infReq, infRemain, infAntArm, infAcc, infDArm, infDLoop, infTArm, infTLoop = ["F", "F", "F", "0", "0", "F", "F", "F", "F", "F", "F"]
        if hFlags=={}:
            iNoFlags += 1
            infNoFlags = "T"
        for flag, flinfo in hFlags.items():
            if flag=="TooShort":
                iTooShort += 1 # counter
                infTooShort = "T" # what output info will be
            elif flag=="CodonIssue":
                iCodonIssue += 1
                infCodonIssue = "T"
            elif flag=="ReqMissing":
                iReqMissing += 1
                infReq = str(flinfo[1])
            elif flag=="Remaining":
                iRemaining += 1
                infRemain = str(flinfo[1])
            elif flag=="Anticodon":
                iAntArm += 1
                infAntArm = "T"
            elif flag=="Acceptor":
                iAcc += 1
                infAcc = "T"
            elif flag=="D_arm":
                iDArm += 1
                infDArm = "T"
            elif flag=="D_loop":
                iDLoop += 1
                infDLoop = "T"
            elif flag=="T_arm":
                iTArm += 1
                infTArm = "T"
            elif flag=="T_loop":
                iTLoop += 1
                infTLoop = "T"
        hInfo[tname] = ("\t".join([infNoFlags, infTooShort, infCodonIssue, infReq, infRemain, infAntArm, infAcc, infDArm, infDLoop, infTArm, infTLoop]) + "\n")
        if hPosType==False:
            continue
        # Process, save actual tRNA info
        for i in sorted(list(hPosType.keys())): # do in order along tRNA; includes SEQUENCE
            wline = "\t".join([tname, str(i)] + hPosType[i][0:2] + [str(hPosType[i][2]), str(hPosType[i][3]),\
                                                                    t2seq[tname][0][i - 1], t2seq[tname][1][i - 1]]) + "\n"
            ostm.write(wline.encode("utf-8"))
    
    ostm.close()
    
    # --- Write tRNA info file
    infoutf = strOStem + "_tRNA2struct_info.txt.gz"
    ostm = gzip.open(infoutf, mode = 'wb')
    ## Preamble
    PREAMBLE = "##File generated by [gitrepo]wormtrna secstruct2seqpieces.py on %s.\n" % time.ctime()
    ostm.write(PREAMBLE.encode("utf-8"))
    ## Header
    HEADER = "\t".join(["tRNA", "PossiblePseudogene", "Intron", "NoFlags", "TooShort", "CodonIssue", "ReqStructuresMissing", "RemainingBases", "AnticodonArmIssue", "AccStemIssue", "DArmIssue", "DLoopIssue", "TArmIssue", "TLoopIssue"]) + "\n"
    ostm.write(HEADER.encode("utf-8"))
    # Data
    iPseud, iIntron = [0, 0]
    for tname, strInfo in hInfo.items():
        oPseud = "F"
        oIntron = "NA"
        if tname in tpseud:
            iPseud += 1
            oPseud = "T"
        if tname in t2intron.keys():
            iIntron += 1
            oIntron = str(t2intron[tname][0]) + "-" + str(t2intron[tname][1])
        wline = tname + "\t" + oPseud + "\t" + oIntron + "\t" + strInfo 
        ostm.write(wline.encode("utf-8"))
    
    ostm.close()
    
    # --- Write summary file
    sumoutf = strOStem + "_tRNA2struct_summary.txt"
    ostm = open(sumoutf, "w")
    
    ## Preamble
    PREAMBLE = "##File generated by [gitrepo]wormtrna secstruct2seqpieces.py on %s.\n" % time.ctime()
    ostm.write(PREAMBLE)
    ## Header
    HEADER = "\t".join(["Category", "Description", "N"]) + "\n"
    ostm.write(HEADER)
    
    # Data
    astrCats = ["PossiblePseudogene", "Intron", "NoFlags", "TooShort", "CodonIssue", "ReqStructuresMissing", "RemainingBases", "AnticodonArmIssue", "AccStemIssue", "DArmIssue", "DLoopIssue", "TArmIssue", "TLoopIssue"]
    astrDescrips = ["Flagged as possible pseudogene in SS input", "Flagged as having intron in SS input", "Not flagged at all - followed all rules perfectly as far as we can tell", "TRNA sequence TOO SHORT - not processed here",\
                    "Codon position not defined OR at a place that can't be made sense of with other structures - not processed", \
                        "At least 1 required sec structure missing", "Not all tRNA bases assigned to a structure", "Anticodon arm stem not as expected (could be just one unpaired base, or worse)", \
                            "Acceptor stem not as expected (could be just one unpaired base, or worse)", "D arm not as expected (could be just one unpaired base, or worse)", "Unable to identify D loop, VERY WEIRD & UNTESTED", \
                                "T arm not as expected (could be just one unpaired base, or worse)", "Unable to identify T loop, VERY WEIRD & UNTESTED"]
    aiNs = [iPseud, iIntron, iNoFlags, iTooShort, iCodonIssue, iReqMissing, iRemaining, iAntArm, iAcc, iDArm, iDLoop, iTArm, iTLoop]
    for i in range(len(astrCats)):
        wline = "\t".join([astrCats[i], astrDescrips[i], str(aiNs[i])]) + "\n"
        ostm.write(wline)
        
    ostm.close()

#%%
def writetfastas(fstem, ht2PosStruct, t2seq, regfa = True, x = False, nameprefix = "", maskantic = True):
    """
    Writes one fasta, x fasta, or both per tRNA secondary structure part. << recoded to (( etc

    Parameters
    ----------
    fstem : string
        Fully articulated path to file writing destination; what segment and file suffix (.fasta or .xfasta) will be added
    ht2PosStruct : DICT
        output of pos2struct; see that for details. hash with tRNA names as keys; values are returns from ht2PosStruct OR a hash flagging why that wasn't/coudn't be run
    t2seq : DICT
        Key: tRNA name 
        Value: list of sequence, 2ndary structure **introns removed if there was an intron**
    regfa : bool, optional
        write normal fasta?. The default is True (normal fasta written).
    x : bool, optional
        write x fasta?. The default is False (only normal fasta written).
    nameprefix : STR, optional
        Prefix to add before tRNA sequence name **including any delimiter**. The default is "".
    maskantic : bool, optional
        N-mask anticodon in antarm out file? Default is True

    Returns
    -------
    None.

    """
    # Set up sequence parts to combine for fastas
    #       keys: names for output fasta files; values: 1-based inds of this structure as in hPosType
    hFas = {"accstemL2d" : [1, 2], # L acceptor stem & any linker
            "darm2ant" : [3, 4, 5, 6], # D arm & any linker. WAS MISSING LINKER
            "antarm" : [7, 8, 9, 10, 11], # Anticodon arm. REMEMBER N OUT ACTUAL ANTICODON
            "anticodon" : [9], # Anticodon itself - NNN will be included in antarm out, actual seq here
            "varloop" : [12], # variable loop
            "tarm2acc" : [13, 14, 15, 16], # t arm and any linker with acceptor stem
            "accstemR2end" : [17, 18]} # acceptor stem final part & overhang
    
    # Iterate over fastas to write
    for struct, aiWhich in hFas.items():
        # Set up file(s)
        if regfa:
            wfastaf = fstem + "_" + struct + ".fasta"
            wfasta = open(wfastaf, "w")
        if x:
            xfastaf = fstem + "_" + struct +  ".xfasta"
            xfasta = open(xfastaf, "w")
        
        # Pull out ALL tRNAs with these structures & write them
        for tname, info in ht2PosStruct.items():
            if len(info)!=3:
                continue
            hPosType, hFlags, hstr2pos = info
            aiNts = []
            aiAnticodon = []
            for i1, posinfo in hPosType.items():
                if posinfo[2] in aiWhich:
                    aiNts.append(i1 - 1) # 0 index this for string
                if struct=="antarm": # ONLY if this, N out anticodon
                    if posinfo[2]==hFas["anticodon"][0] and maskantic: # only adds in here if one should mask it
                        aiAnticodon.append(i1 - 1)
            if len(aiNts) > 0: # If has this structure, write it!
                outName = nameprefix + tname # seq name
                # Seq info to write
                aiNts.sort()
                astrNts = []
                astrStructs = []
                for i in aiNts:
                    if i in aiAnticodon: # N out if appropriate
                        astrNts.append("N")
                    else:
                        astrNts.append(t2seq[tname][0][i]) # this looks like it should be fine but clearly...isn't
                    # Recode sec struct nicely
                    if t2seq[tname][1][i]==">":
                        astrStructs.append("(")
                    elif t2seq[tname][1][i]=="<":
                        astrStructs.append(")")
                    else:
                        astrStructs.append(t2seq[tname][1][i])
                        
                # Write fasta if desired
                if regfa:
                    wfasta.write(">" + outName + "\n")
                    wfasta.write("".join(astrNts) + "\n")
                # Write x fasta if desired
                if x:
                    xfasta.write(">" + outName + "\n")
                    xfasta.write("".join(astrNts) + "\n")
                    xfasta.write("".join(astrStructs) + "\n")
    
        # Close
        if regfa:
            wfasta.close()
        if x:
            xfasta.close()

#%%
def main():
    #### Parse arguments
    argp = argparse.ArgumentParser(prog = "secstruct2seqpieces.py",
                                   description = "Break tRNAs into their structures (t arm etc) \
                                       based on tRNAScan-SE secondary structure. For alignment.")
    argp.add_argument("-trnass", dest = "strSS", required = True, type = str,
                      metavar = "trnas.SS",
                      help = "Path to tRNAscan-SE .SS file for tRNAs to process here")
    argp.add_argument("-outdir", dest = "strOutDir", required = True, type = str,
                      metavar = "outdir",
                      help = "Output directory")
    argp.add_argument("-out", dest = "strOut", required = True, type = str,
                      metavar = "outfilestem",
                      help = "Prefix for output files")
    argp.add_argument("-forcepseud", dest = "bFakeCods", required = False, type = bool,
                      default = True,
                      metavar = "[True, False]",
                      help = "If gene is pseudo with codon listed as NNN as 0-0, try to force processing by trying codon at 34-36. (Default : True)")
    argp.add_argument("-fasta", dest = "bFasta", required = False, type = bool,
                      default = True,
                      metavar = "[True, False]",
                      help = "Write a fasta for each sensible secondary structure chunk of tRNAs.\
                          [default : True]")
    argp.add_argument("-xfasta", dest = "bXFasta", required = False, type = bool,
                      default = True,
                      metavar = "[True, False]",
                      help = "Write a X fasta (where line after sequence is secondary structure in (). format) for each sensible secondary structure chunk of tRNAs.\
                          [default : True]")
    argp.add_argument("-fanameprefix", dest = "strFastPrefix", required = False, type = str, 
                      default = "",
                      metavar = "species_prefix",
                      help = "Optional prefix for all sequences in output fasta/xfasta files - e.g., species identifier if files will be compined across multiple species.")
    argp.add_argument("-maskanticodon", dest = "bMask", required = False, type = bool,
                      default = True,
                      metavar = "[True, False]",
                      help = "N out anticodons in output fasta/X fasta? (they are saved in their own FASTA) [default : True]")
    args = argp.parse_args()
    
    #### Set up logging
    statement = "Starting at %s..." % time.ctime()
    print(statement)
    log_list = []
    log_list.append(statement + '\n')
    
    #### Read in tRNA info
    statement = "Reading in tRNA info (at %s)..." % time.ctime()
    print(statement)
    log_list.append(statement)
    t2loc, t2locabs, loc2t, t2strand, t2len, t2codon, t2seq, t2intron, tpseud = readss(args.strSS)
    
    #### Map position to structure
    statement = "Mapping tRNA positions to tRNA structure (at %s)..." % time.ctime()
    print(statement)
    log_list.append(statement + '\n')
    ht2PosStruct = pos2struct(t2len, t2loc, t2codon, t2seq, iMinLn = 71, bFakeCods = args.bFakeCods)
       
    #### Save & summarize
    statement = "Writing tRNA sec structure pieces info outputs (at %s)..." % time.ctime()
    print(statement)
    log_list.append(statement + '\n')
    strOutStem = args.strOutDir + "/" + args.strOut
    saveinfo(t2seq, t2intron, tpseud, ht2PosStruct, strOStem = strOutStem)
    
    #### Split seqs into FASTAs based on structure if desired [input option]
    strNamePref = ""
    if args.strFastPrefix!="":
        strNamePref = "".join(args.strFastPrefix.strip().split()) + "_"
        
    if args.bFasta or args.bXFasta:
        statement = "Writing tRNAs in sec structure chunks to fasta file(s) (at %s)..." % time.ctime()
        print(statement)
        log_list.append(statement + '\n')
        fastem = args.strOutDir + "/" + args.strOut
        writetfastas(fstem = fastem, ht2PosStruct = ht2PosStruct, t2seq = t2seq, regfa = args.bFasta, x = args.bXFasta, nameprefix = strNamePref, maskantic = args.bMask)
    
    #### Done! Log
    statement = "Processing done! Writing log file and wrapping up (at %s)..." % time.ctime()
    print(statement)
    log_list.append(statement + '\n')
    logfilepath = args.strOutDir + "/secstruct2seqpieces-" + time.strftime("%Y%m%d-%H%M%S") + ".log"
    logfile = open(logfilepath, "w")
    PREAMBLE = "##File generated by [gitrepo]wormtrna secstruct2seqpieces.py at %s.\n" % time.ctime()
    logfile.write(PREAMBLE)
    for line in log_list:
        logfile.write(line)
    logfile.close()

if __name__=="__main__":
    main()