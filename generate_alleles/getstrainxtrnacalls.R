#! /usr/bin/env/ Rscript
# Summarize tRNA alleles information from Python scripts, tRNAscan-SE, catalogmissingness_trnavcfscripts.R                                            
# by Avery Davis Bell, begun 2024.11.19

if(length(.libPaths())==1){
  mylibloc <- .libPaths()[1]
}else{ # Presumes on PACE. This is not the best...
  mylibloc <- .libPaths()[grep(R.Version()$platform, .libPaths())]
}
require(data.table, lib.loc = mylibloc) 
require(argparser)

#### Functions ####
getfanames<-function(fapath){
  # Gets character vector of seq names from a fasta 
  # In: fapath, path to fasta
  # Out: character vector of same length as seqs in fasta. Of seq names - doesn't keep anything else
  
  # Get just name lines
  ls<-readLines(fapath)
  ls<-ls[substr(ls, 1, 1)==">"]
  
  # Strip > and anything after any possible whitespace
  snames<-sapply(strsplit(ls, ">|\\s"), function(x) x[2])
  return(snames)
}

strainxtrna<-function(trnas, missinfo){
  # Gets each strain's allele for each tRNA! Then overlays missing info and NAs out where genotype call was missing
  # In: trnas, trna information data.table. Key columns are tRNA, allelename, Strains (comma-sep ID of strains called with this allele [or with missingness])
  #     missinfo, missing genotype call info data.table. Key columns are tRNA, missingOrHet (comma-sep ID of all strains called with this allele)
  # Out: data.table; rows are strains and columns (other than the first) are tRNAs.
  #         First column is strain, then all others contain tRNA allele ID called for that strain at that gene.
  #         If there was no data or known missingness at that gene in that strain, it's NA'ed out
  
  # Set up all strains, tRNAs to fill in info for
  strns<-unique(c(unlist(lapply(missinfo$missingOrHet, function(x) strsplit(x, ",")[[1]])),
                  unlist(lapply(trnas$Strains, function(x) strsplit(x, ",")[[1]]))))
  strns<-strns[strns!="None"] # just in case it slipped through
  
  gs<-sort(unique(trnas$tRNA))
  
  # Fill in info from trnas
  ## Data to extract from
  pert<-lapply(gs, function(g){
    tmp<-lapply(trnas[tRNA==g, Strains], function(x) strsplit(x, ",")[[1]])
    names(tmp)<-trnas[tRNA==g, allelename]
    return(tmp)
  })
  names(pert)<-gs
  ## Set up matrix
  sxt<-matrix(nrow = length(strns), ncol = length(gs))
  rownames(sxt)<-strns
  colnames(sxt)<-gs
  ## Fill in matrix
  for(g in gs){ # for loops fill matrices better, often
      oneg<-pert[[g]]
      for(myn in names(oneg)){
        mystrns<-oneg[[myn]]
        sxt[mystrns, g]<-myn
      }
  }
  
  # Overlay missingness information
  pertm<-lapply(missinfo$missingOrHet, function(x) strsplit(x, ",")[[1]])
  names(pertm)<-missinfo$tRNA
  lns<-sapply(pertm, length)
  pertm<-pertm[lns!=0]
  
  for(g in names(pertm)){
    sxt[pertm[[g]], g]<-NA # these are missing, so we're overlaying NAs
  }
  
  # Format (data.table?) & Return
  out<-data.table(sxt, keep.rownames = T)
  setnames(out, "rn", "strain")
  return(out)
}

getalleleinfo<-function(tsout, fanames){
  # Collects allele information for ALL observed alleles (including those not found by tRNAscan-SE)
  # In: tsout, tRNAscan-SE .out relevant output as a data.table. Columns: allelename, AA, Codon, Infernal, AlleleCM, IsotypeScore, Note
  #     fanames, vector of all tRNA allele names in the study (ie those provided as input to tRNAscan-SE)
  # Out: data.table with one row per allele in fanames. Columns:
  #     tRNA, tRNA gene ID
  # allelename, specific allele this row describes (as in input)
  # AA,  amino acid (as in input)
  # Codon, codon (as in input)
  # Infernal, infernal score (as in input)
  # AlleleCM, expected amino acid based on backbone [confirm] (as in input)
  # IsotypeScore, (as in input)
  # Note, (as in input)
  # VariableInPop - T or F: is more than one allele observed in population?
  # Classification - newly added key here. Is this allele Lost, Altered (isotype switch), Best (highest Infernal score), or Functional (lower but func infernal score)
  
  # Subfunctions
  classonetrna<-function(trnadt, functhresh = 20){
    # Classifies each allele in one tRNA as Altered/Lost/Functional/Best
    # In: trnadt, data.table for one tRNA. columns must include AA, AlleleCM, Infernal, Note
    #     functhresh, functional threshold for Infernal score [should be redundant!]
    # Out: data.table of classifications for each allele. Columns VariableInPop (T/F), Classification (Altered/Lost/Functional/Best)

    # Find any that are lost and altered; save their classifications and remove from processing
    out<-rep(NA, nrow(trnadt))
    if(length(out)==1){
      isvar<-F
    }else{
      isvar<-T
    }
    out[trnadt[, AA!=AlleleCM | grepl("IPD", Note)]] <- "Altered" # sometimes things are IPD and pseudo - want pseudo to override so it goes last
    out[trnadt[, grepl("pseudo", Note) | grepl("secondary", Note) | grepl("undetermined isotype", Note) | Codon=="NNN"]] <- "Lost" # Codon=="NNN" added later; big miss!
    
    
    # If any remain to classify, ID best vs functional
    if(sum(is.na(out))>0){
      toclass<-which(is.na(out))
      dat.toclass<-trnadt[toclass,]
      best<-which(dat.toclass$Infernal == max(dat.toclass$Infernal))
      func<-which(dat.toclass$Infernal < max(dat.toclass$Infernal) & dat.toclass$Infernal > functhresh)
      out[toclass][best]<-"Best"
      out[toclass][func]<-"Functional"
    }
    
    # Return
    return(data.table(VariableInPop = isvar, Classification = out))
  }
  
  # Data.table of any not included in tsout
  others<-data.table(tRNA = sapply(strsplit(fanames[!fanames%in%tsout$allelename], "_"), function(x) paste(x[1:(length(x) -2)], collapse = "_")), # this *should* dealwith cases where there's an underscore in chromosome ID
                     allelename = fanames[!fanames%in%tsout$allelename],
                     AA = "Undet", Codon = "NNN", Infernal = 0, AlleleCM = "Undet",
                     IsotypeScore = NaN, Note = "undetermined isotype / not found by tRNAscan-SE")
  
  # Add Altered/Invariant/Lost/Functional/Best annotations for alleles that are in tsout
  dat<-copy(tsout)
  # dat[, tRNA:=dat[, tstrsplit(allelename, "_")]$V1] # doesn't work when _ in chr ID
  dat[, tRNA:=sapply(strsplit(allelename, "_"), function(x) paste(x[1:(length(x) -2)], collapse = "_"))]
  setcolorder(dat, "tRNA")
  dat<-rbind(dat, others)
  setkey(dat, tRNA)
  dat<-data.table(dat, dat[, classonetrna(.SD), by = tRNA][,.(VariableInPop, Classification)])
  
  # Return
  return(dat)
}

allelesumm<-function(allinfo.tscan, sxtdat){
  # Gets counts of strains observing each allele for each tRNA
  # In: allinfo.tscan, allele info (output of getalleleinfo) - key columns are tRNA, allelename
  #     sxtdat, strains x tRNAs data.table as output from strainxtrna 
  # Out: allinfo.tscan with new columns added:
  #   n.allele, # strains with this allele
  #   n.called, # strains with allele called at this tRNA
  #   n.missing, # strains with no allele called at this tRNA
  
  # Subfunctions
  ctone<-function(trna.id, allinfo.tscan, sxtdat){
    # Gets 3 counts for one tRNA for all alleles
    # In: trna.id, tRNA to count
    #     allinfo.tscan, info on all alleles for this trna (possibly more - narrowing is done)
    #     sxtdat, data to count in
    # Out: data.table with columns allelename, n.allele, n.called, n.missing
    trnadat<-allinfo.tscan[tRNA==trna.id, ]
    thisdat<-sxtdat[, get(trna.id)]
    
    # count each allele
    nas<-sapply(trnadat$allelename, function(x) sum(thisdat==x, na.rm = T))
    
    # Count totals
    nc<-sum(!is.na(thisdat))
    nm<-sum(is.na(thisdat))
    
    # Format & return
    out<-data.table(allelename = trnadat$allelename, 
                    n.allele = nas,
                    n.called = nc,
                    n.missing = nm)
    return(out)
  }
  
  # Get all counts
  cts<-rbindlist(lapply(allinfo.tscan[, unique(tRNA)], ctone, allinfo.tscan, sxtdat)) 
  setkey(allinfo.tscan, allelename)
  setkey(cts, allelename)
  out<-allinfo.tscan[cts]
  setkey(out, tRNA)
  
  # Return
  return(out)
}

strainsumm<-function(allinfo.tscan, sxtdat, cs = c("Best", "Functional", "Altered", "Lost")){
  # Summarize number of observed alleles per strain per classification overall and per amino acid
  # In: allinfo.tscan, allele info (output of getalleleinfo) - key columns are allelename, AA, Classification
  #     sxtdat, strains x tRNAs data.table as output from strainxtrna 
  # Out: long format data with one row per strain, amino acid combination (and all together). Columns:
  #   strain,
  # AA, amino acid alleles/tRNAs are combined for - if 'all', that's top-level strain summary
  # <one for each classification - Altered, lost, best, functional> number of alleles called <this classification> for this amino acid in this strain
  # nNA, number of alleles not called due to genotype missingness
  #       [per AA is the ones actually coded for -- not the alleleCM. could be delivering wrong ones or...?]
  
  # Subfunctions
  onestrnsumm<-function(allinfo.tscan, strndat, cs){
    # Summarizes # alleles of which class per strain
    # In: allinfo.tscan, allele info (output of getalleleinfo) - key columns are allelename, AA, Classification. Keyed by allelename!!
    #     strndat, one-row data.table data for one strain; one column for each tRNA
    #     cs, all classifications in Classification column to quantify
    # Out: data.table with one row per AA (and all, combined). Columns AA, <all in cs>, nNA - all but first are counts for in this strain how many alleles were in each category
    
    # Ones with non-NA alleles
    theseals<-allinfo.tscan[unlist(strndat)] # almost works, but need the tRNA ID even if allele is NA....
    theseals<-theseals[!is.na(tRNA), .(tRNA, allelename, AA, AlleleCM, Classification)]
    
    # Fill in those with NA alleles
    thesena<-unique(allinfo.tscan[tRNA %in% names(strndat)[is.na(strndat)]][, .(tRNA, AA)]) ## but need to strip extra alleles - that's what 2nd part does
    thesena<-thesena[AA!="Undet", ] # sometimes alleles have multiple AAs - not useful here
    thesena[, `:=`(allelename = NA, AlleleCM = NA, Classification = NA)]
    setcolorder(thesena, names(theseals))
    theseals<-rbind(theseals,thesena) 
    
    # Count/summarize
    ## All
    allc<-theseals[ , sapply(cs, function(x) sum(Classification==x, na.rm = T))]
    allc.na<-theseals[, sum(is.na(Classification))]
    allc<-data.table(AA = "all", t(allc), allc.na)
    ## Per AA
    setkey(theseals, AA)
    peraa<-theseals[ , sapply(cs, function(x) data.table(t(sum(Classification==x, na.rm = T)))), by = AA]
    peraa.na<-theseals[, sum(is.na(Classification)), by = AA]
    peraa<-peraa[peraa.na]
    ## Combine
    outc<-rbindlist(list(allc, peraa), use.names = F)
    setnames(outc, c("AA", cs, "nNA"))
    
    # Return
    return(outc)
  }
  
  # Set up
  # cs<-allinfo.tscan[, unique(Classification)] # all possible classifications to count. DECIDED TO PROVIDE AS FN ARG SO SAME EACH TIME
  setkey(allinfo.tscan, allelename)
  setkey(sxtdat, strain)
  
  # Get summary
  perstrn<-sxtdat[, onestrnsumm(allinfo.tscan, .SD, cs), by = strain]
  return(perstrn)
}

geneupsumm<-function(al.cts, cs = c("Best", "Functional", "Altered", "Lost")){
  # Summarizes number total alleles, number with each classification, if any strains have missing genotype calls at gene at multiple levels:
  #     per tRNA gene
  #     per anticodon 
  #     per AA
  # In: al.cts, output of allelesumm: information for each allele. Columns as described in that fn
  #     cs, categories to count number of alleles of - in Classification column of input
  # Out: List of data.table summaries:
  #       $bygene, one row per tRNA gene. Columns:
  #           tRNA, gene ID
  #					AA, AA codon codes for
  #					Codon, actual codon in tRNA
  #					VariableInPop, is gene variable in population
  #					nAlleles, # alleles observed in population
  #					anyMissingCalls, T or F, any alleles have any missing genotype calls in sequence
  #					Best, number alleles classified as best
  #					Functional, number alleles classified as functional
  #					Altered, number alleles classified as altered
  #					Lost, number alleles classified as altered
  #       $bycodon, one row per triplet codon - the one that's actually CALLED, not the backbone!  Columns:
  #         AA, AA codon codes for
  #					Codon, actual codon in tRNA
  #					nGenes, # genes that have this codon
  #					VariableInPop, is ANY gene for this codon variable in population
  #					nAlleles, # alleles observed in population across any gene
  #					nGenesMissingCalls, how many genes have alleles with any missing genotype calls (./. in VCF)
  #					Best, number alleles classified as best across all genes
  #					Functional, number alleles classified as functional across all genes
  #					Altered, number alleles classified as altered across all genes
  #					Lost, number alleles classified as altered across all genes
  #       $byaa, one row per amino acid - the one that's actually CALLED, not the backbone! (from the codon)
  #         AA, AA summarized here
  #					nCodons, # codons for this AA in this genome
  #					nGenes, # genes that have this AA
  #					VariableInPop, is ANY gene for this AA variable in population
  #					nAlleles, # alleles observed in population across any gene for this AA
  #					nGenesMissingCalls, how many genes have alleles with any missing genotype calls (./. in VCF)
  #					Best, number alleles classified as best across all genes for this AA
  #					Functional, number alleles classified as functional across all genes for this AA
  #					Altered, number alleles classified as altered across all genes for this AA
  #					Lost, number alleles classified as altered across all genes for this AA
  
  # Subfunctions
  ctclass<-function(oneset, cs){
    # For input rows, counts how many alleles in each classification (# rows where Classification==this) and IF any strains have missingness here (n.missing > 0 any row)
    # In: oneset, al.cts subset to summarize
    #     cs, categories to count number of alleles of - in Classification column of input
    # Out: one-row data.table with columns;
    #       nAlleles, number of alleles observed
    #       anyMissingCalls, did any of these alleles have any strains with any missing genotypecalls
    #       < one for each in cs>
    
    myc<-oneset[ , sapply(cs, function(x) sum(Classification==x, na.rm = T))]
    myc.na<-oneset[, sum(n.missing>0) > 0]
    myc<-data.table(nAlleles = sum(myc), anyMissingCalls = myc.na, t(myc))
    return(myc)
  }
  
  # Per tRNA
  setkey(al.cts, tRNA)
  t.summ<-data.table(al.cts[, .(AA[1], Codon[1], VariableInPop[1]), by = tRNA],
                     al.cts[, ctclass(.SD, cs), by = tRNA][, .SD, .SDcols = -1])
  setnames(t.summ, c("V1", "V2", "V3"), c("AA", "Codon", "VariableInPop"))
  
  # Per codon
  setkey(t.summ, AA, Codon)
  c.summ<-data.table(t.summ[, .(.N, sum(VariableInPop, na.rm = T)>0, sum(nAlleles), sum(anyMissingCalls)), by = .(AA, Codon)], # adds n genes, totals
                     t.summ[, data.table(t(sapply(cs, function(x) sum(get(x))))), by = .(AA, Codon)][, .SD, .SDcols = -c(1, 2)]) # SUM of how many alleles in each class across ALL GENES here
  setnames(c.summ, c("N", "V2", "V3", "V4"), c("nGenes", "VariableInPop", "nAlleles", "nGenesMissingCalls")) # here it's are ANY of genes VariableInPop
  
  # Per AA
  setkey(c.summ, AA)
  a.summ<-data.table(c.summ[, .(.N, sum(nGenes), sum(VariableInPop, na.rm = T)>0, sum(nAlleles), sum(nGenesMissingCalls)), by = AA], # adding number codons, total n genes, etc
                     c.summ[,  data.table(t(sapply(cs, function(x) sum(get(x))))), by = AA][, .SD, .SDcols = -1])
  setnames(a.summ, c("N", "V2", "V3", "V4", "V5"), c("nCodons", "nGenes", "VariableInPop", "nAlleles", "nGenesMissingCalls"))
  
  # Return
  return(list(bygene = t.summ, bycodon = c.summ, byaa = a.summ))
}

#### Arguments & inputs ####
p<-arg_parser("Summarize tRNA alleles information from Python scripts, tRNAscan-SE, catalogmissingness_trnavcfscripts.R", 
              name = "getstrainxtrnacalls.", hide.opts = TRUE)

p<-add_argument(p, "--trnainfo",
                help = "Path to *_strain_trnas_info.txt output of build_alt_sequences.py",
                type = "character")
p<-add_argument(p, "--missinfo",
                help = "Path to *_trna_variant_info_wmissingness_pertrna.txt output of catalogmissingness_trnavcfscripts.R",
                type = "character")
p<-add_argument(p, "--trnascanout",
                help = "Path to strain-specific fasta tRNAscan-SE output (.out file)",
                type = "character")
p<-add_argument(p, "--trnafasta",
                help = "Path to strain_trnas.fa from build_alt_sequences.py (used to check allele names)",
                type = "character")
p<-add_argument(p, "--refstrain",
                help = "ID of reference strain(s) to strip out of final data; comma-separated if more than one",
                type = "character")
p<-add_argument(p, "--baseoutname",
                help = "Base name for all output files",
                type = "character",
                default = "out")
p<-add_argument(p, "--outdir",
                help = "Outer output directory. Sub-directories will be created internally. **NB: if you provide getwd() here (quote wrapped), current directory will be used",
                type = "character")

# Parse arguments
cat("....Parsing arguments....\n")
p<-parse_args(p)

# outdir
if(p$outdir=="getwd()"){
  p$outdir<-getwd()
}
if(!dir.exists(p$outdir)){dir.create(p$outdir, recursive = T)}
setwd(p$outdir)

# Ref strain
refstrains<-strsplit(p$refstrain, ",")[[1]]

# Read in data
cat("....Reading in data....\n")
trnas<-fread(p$trnainfo, header = T) # tRNA info: who has what, use to map onto allele names [add allele name here!!]
setnames(trnas, "##Chr", "Chr")

missinfo<-fread(p$missinfo, header = T) # Who's missing per tRNA

fanames<-getfanames(p$trnafasta)

tsout<-fread(p$trnascanout, skip = 3) # tRNAscan-SE out. **if format changes this might not work!!
setnames(tsout, c("allelename", "drop1", "drop2", "drop3", "AA", "Codon", "drop4", "drop5", "Infernal", "AlleleCM", "IsotypeScore", "Note"))
tsout<-tsout[, .(allelename, AA, Codon, Infernal, AlleleCM, IsotypeScore, Note)]

# alinfo<-fread(p$alleleinfo, header = T) # info for each allele - but without full 'allele name' [needed for summaries, not for generating data matrix]
# setnames(alinfo, c("##tRNA", "#strains"), c("tRNA", "nStrains"))

# *** ideally get the lost, functional, best chars into whatever is downstream of tsout

#### Generate data matrix ####
cat("....Generating data matrix....\n")
# Annotate trna info with allele names 
trnas[, allelename:=paste(tRNA, trnas[, tstrsplit(Strains, ",")]$V1, numStrains, sep = "_"),]

# Get matrix from trna info; add missingness on top
sxtdat<-strainxtrna(trnas, missinfo)
## Strip reference strain out if know who it is! (and save separately)
if(length(refstrains)>0){
  setkey(sxtdat, strain)
  refout<-as.data.table(t(sxtdat[refstrains, .SD, .SDcols = -1 ]), keep.rownames = T)
  setnames(refout, c("tRNA", refstrains))
  write.table(refout, file.path(p$outdir, paste0(p$baseoutname, "_referencestrain", p$refstrain, "_pertRNAalleles.txt")),
              sep = "\t", quote = F, row.names = F)
  
  sxtdat<-sxtdat[!strain%in%refstrains, ]
}

# Save
write.table(sxtdat, gzfile(file.path(p$outdir, paste0(p$baseoutname, "_strainsxtrnas_alleles_wmissing.txt.gz"))),
            sep = "\t", quote = F, row.names = F)

#### Generate summaries ####
cat("....Generating summaries....\n")
# Get full allele information (with allele names mapped on)
allinfo.tscan<-getalleleinfo(tsout, fanames)

# By allele (easiest: # strain with each allele)
al.cts<-allelesumm(allinfo.tscan, sxtdat)
write.table(al.cts, file.path(p$outdir, paste0(p$baseoutname, "_alleleinfo_counts_wmissing.txt")),
            sep = "\t", quote = F, row.names = F)

# By strain (per AA & all together)
bystrain<-strainsumm(allinfo.tscan, sxtdat)
write.table(bystrain, file.path(p$outdir, paste0(p$baseoutname, "_strain_alleletype_counts_wmissing.txt")),
            sep = "\t", quote = F, row.names = F)
bystrain.var<-strainsumm(allinfo.tscan[VariableInPop==T, ], sxtdat)
write.table(bystrain.var, file.path(p$outdir, paste0(p$baseoutname, "_strain_alleletype_counts_wmissing_variablegenesonly.txt")),
            sep = "\t", quote = F, row.names = F)

# By gene, codon, AA: this is categorized by the codon/AA actually in the sequence (can differ across alleles)
bytca<-geneupsumm(al.cts)
write.table(bytca$bygene, file.path(p$outdir, paste0(p$baseoutname, "_bygeneup_tRNA_counts_wmissing.txt")),
            sep = "\t", quote = F, row.names = F)
write.table(bytca$bycodon, file.path(p$outdir, paste0(p$baseoutname, "_bygeneup_codon_counts_wmissing.txt")),
            sep = "\t", quote = F, row.names = F)
write.table(bytca$byaa, file.path(p$outdir, paste0(p$baseoutname, "_bygeneup_AA_counts_wmissing.txt")),
            sep = "\t", quote = F, row.names = F)

# Do I have anywhere: total number genes that are variable [broken down by pseudo and not]? This may be next script territory...(data and plot)
#       Kinda goes with number alleles per, that sort of thing

#### Script completion message & session information ####
cat("....getstrainxtrnacalls.R processing complete! Session information:....\n")
sessionInfo()

