#! /usr/bin/env/ Rscript
# Investigate potential stem compensatory mutations                               
# by Avery Davis Bell, begun 2026.02.10

if(length(.libPaths())==1){
  mylibloc <- .libPaths()[1]
}else{ # Presumes on PACE. This is not the best...
  mylibloc <- .libPaths()[grep(R.Version()$platform, .libPaths())]
}
require(data.table, lib.loc = mylibloc) 
require(argparser)
require(ggplot2)
require(ggforce)

# --- plotting
myggtheme<-theme_bw() +
  theme(axis.title.x = element_text(size = 12), axis.text.x = element_text(size = 11), 
        axis.title.y = element_text(size = 12), axis.text.y = element_text(size = 11),
        legend.title = element_text(size = 12), legend.text = element_text(size=11), 
        strip.text.x = element_text(size = 11), strip.text.y = element_text(size = 11), title = element_text(size = 14), 
        strip.text.x.top = element_text(size = 11), strip.text.x.bottom = element_text(size = 11), strip.text.y.right = element_text(size = 11), strip.text.y.left = element_text(size = 11),
        plot.subtitle = element_text(size = 13))

#### Functions ####
# --- data wrangling
fread.mult<-function(sinfo, exampfile, myselect = NA, ...){
  # Reads in a file for each species in sinfo
  # In: sinfo, data.table with columns infilename (exactly how all files have this species in their name), displayname, shortname
  #     exampfile, path to file to process but wherever infilename is, should say SAMP instead. Should have header
  #     myselect, OPTIONAL list of columns to select when reading in (drops others)
  #     ... passed to fread
  # Out: data.table containing whatever's in exampfile for each species in sinfo. Columns added to the front are displayname and shortname
  
  dat<-rbindlist(
    lapply(1:nrow(sinfo), function(s){
      if(all(is.na(myselect))){
        dat1<-fread(file = gsub("SAMP", sinfo[s, infilename], exampfile), header = T, ...)
      }else{
        dat1<-fread(file = gsub("SAMP", sinfo[s, infilename], exampfile), header = T, select = myselect, ...)
      }
      dat1<-data.table(dat1, sinfo[s, .(displayname, shortname)])
      setcolorder(dat1, c("displayname", "shortname"))
    })
  )
  return(dat)
}

ckpair<-function(char1dt, char2dt, char1name = "homAlt", char2name = "homAlt",
                 keepseparate = c('pos', 'ref', 'alt', 'nHomRef', 'nHomAlt', 'nNotMissingHet', 'nMissingHet', 'tRNA_pos', 'substructure', 'nSubStructure', 'substructure_pos')){
  # Checks if same ID(s) occurs in char1 and char2, which are single character objects containing lists of IDs separated by commas
  # In: car1dt and char2dt, two ONE ROW data.tables, each of which has any columns to append to output 
  #             as well as the character objects of interest in the columns with char1name and car2name
  #         ***same structure otherwise!!
  #     char1name and char2name, column names of objs to check
  #     keepseparate, vector of column names to include in output for both _1 and _2 EVEN IF they have same value
  # Returns NULL if no IDs overlap between the two inputs of interest; ONE rowdata.table if they do with:
  #     All columns prior to charname columns from input data tables, suffixed with _1 and _2 where they had diff values or were in keepseparate
  #         and keeping the same in common where they had the same value
  #     nInCommon, number IDs shared between these two [same in each row]
  #     whichInCommon, comma-separated IDs shared between these two
  
  id1<-strsplit(char1dt[, get(char1name)], ",")[[1]]
  id2<-strsplit(char2dt[, get(char2name)], ",")[[1]]
  if(sum(id1%in%id2)>0){
    whshared<-id1[id1%in%id2]
    # Get columns to spit out
    c1out<-copy(char1dt[, .SD, .SDcols = -which(names(char1dt)==char1name)])
    c2out<-copy(char2dt[, .SD, .SDcols = -which(names(char2dt)==char2name)])
    ## Keep common ones once, remove from the others
    #commns<-names(c1out)[c1out==c2out]
    commns <- names(c1out)[vapply( # tricky when there's an NA or other mild weirdness, including NA = NA now
      names(c1out),
      function(nm) {
        a <- c1out[[nm]]; b <- c2out[[nm]]
        # treat factors as character to avoid spurious non-equality
        if (is.factor(a)) a <- as.character(a)
        if (is.factor(b)) b <- as.character(b)
        # TRUE when equal, including NA==NA
        isTRUE(all.equal(a, b, check.attributes = FALSE))
      },
      logical(1L)
    )]
    commns<-commns[!commns%in%keepseparate]
    
    commonout<-c1out[, commns, with = F]
    c1out<-c1out[, .SD, .SDcols = -commns]
    c2out<-c2out[, .SD, .SDcols = -commns]
    ## Name suffixes
    setnames(c1out, names(c1out), paste(names(c1out), "1", sep = "_"))
    setnames(c2out, names(c2out), paste(names(c2out), "2", sep = "_"))
    # Combine with sharing info
    out<-data.table(commonout,
                    c1out, 
                    c2out,
                    nInCommon = length(whshared),
                    whichInCommon = paste(whshared, collapse = ","))
  }else{
    out<-NULL
  }
  return(out)
}

readss<-function(ssfile, restructuress = F){
  # <copied from other scripts of mine>
  # Reads in tRNA secondary structures from tRNAscan-SE .ss output file; returns data.table with one row per tRNA
  # In: ssfile, path to tRNAscan-SE .ss output file
  #     restructuress, T or F: replace default << in structure seqs with (( and >> with ))
  # Out: data.table with one row per tRNA described in ssfile. Columns:
  #     #       name, name of tRNA from first row of this record  ssfile
  # length, length of tRNA from first row of this record ss file
  # type, protein coding type of amino acid (e.g. Leu)
  # anticodon, sequence of anticodon
  # anticodon.start, first position of anticodon in tRNA
  # anticodon.stop, last position of anticodon in tRNA (first + 2)
  # type.score, score from the type line of SS (isotype score?)
  # sequence, tRNA nucleotide sequence
  # secstruct, secondary structure - as in ssfile if restructuress=F; otherwise changed so (s used instead of <s
  
  ## Subfunctions
  ssvec2dt<-function(ssvec){
    # Processes the lines from a .ss file corresponding to one tRNA into one row of a data.table (keeping most but not all info)
    # In: ssvec, 6 element vector with each element corresponding to that row of a tRNA's info in .ss file
    # Out: one-row data.table, columns:
    #       name, name of tRNA from first row ssfile
    # length, length of tRNA from first row ss file
    # type, protein coding type of amino acid (e.g. Leu)
    # anticodon, sequence of anticodon
    # anticodon.start, first position of anticodon in tRNA
    # anticodon.stop, last position of anticodon in tRNA (first + 2)
    # type.score, score from the type line of SS (isotype score?)
    # sequence, tRNA nucleotide sequence
    # secstruct, secondary structure as in ssfile
    
    # Get info nicely formatted
    ssvec.split<-strsplit(ssvec, split = "\t")
    ## Row 1
    tname<-strsplit(ssvec.split[[1]][1], " ")[[1]][1]
    tlen<-strsplit(ssvec.split[[1]][2], " ")[[1]][2]
    ## 'Type' info
    tline<-which(sapply(ssvec.split, function(x) startsWith(x[1], "Type")))
    splittype<-strsplit(ssvec.split[[tline]], " ")
    type<-splittype[[1]][2]
    anticodon<-splittype[[2]][2]
    anticstart<-as.integer(strsplit(splittype[[2]][4], "-")[[1]][1])
    anticstop<-as.integer(strsplit(splittype[[2]][4], "-")[[1]][2])
    typescore<-as.numeric(splittype[[3]][2])
    ## Sequence
    seqline<-which(sapply(ssvec.split, function(x) startsWith(x[1], "Seq")))
    seq<-strsplit(ssvec.split[[seqline]], " ")[[1]][2]
    ## Structure
    strline<-which(sapply(ssvec.split, function(x) startsWith(x[1], "Str")))
    sstruct<-strsplit(ssvec.split[[strline]], " ")[[1]][2]
    
    # Combine into data.table; return
    out<-data.table(name = tname,
                    length = tlen,
                    type = type,
                    anticodon = anticodon,
                    anticodon.start = anticstart,
                    anticodon.stop = anticstop,
                    type.score = typescore,
                    sequence = seq,
                    secstruct = sstruct)
    return(out)
  }
  
  # Read in as giant string char
  allstr<-readLines(ssfile)
  if(sum(substr(allstr, 1, 2)=="#")> 0){ # remove any comment lines
    allstr<-allstr[-which(substr(allstr, 1, 2)=="#")]
  }
  
  # Get list, one vector per tRNA
  ss.list<-list()
  ss.one<-c()
  for(i in 1:length(allstr)){
    if(allstr[i]==""){
      ss.list[[length(ss.list) + 1]]<-ss.one
      ss.one<-c()
    }else{
      ss.one<-c(ss.one, allstr[i])
    }
  }
  
  # Format into datatable, one row per tRNA
  ssdt<-rbindlist(lapply(ss.list, ssvec2dt))
  
  # Optionally replace << with paren characters
  if(restructuress==T){
    ssdt[, secstruct:= ssdt[, gsub("<", ")", gsub(">", "(", secstruct))]] # this looks 'backwards' but needs to be this way for RNAplot
  }
  
  # Return
  return(ssdt)
}

# --- Analysis (moreso)
charpairals<-function(shpair, gvs, als){
  # Tags allele information with information about variant pairs of interest
  # In: shpair, variant pair info for one tRNA gene (might be multiple pairs)
  #     gvs, all variants in this tRNA gene as from input
  #     als, allele info for this gene as from input
  # Out: data.table with one row per allele. Each row contains critical input allele information and columns:
  # is.reference, is this allele the reference
  # is.pair, this allele is one with at least one putatively comp. pair
  # n.pairs.in, count of N potentially compensatory pairs this allele harbors
  # has.pair.mut, T/F potential comp pair mutation (any!) seen in this allele
  # has.nonpair.mut, T/F other mutations that aren't in comp pair
  # n.muts, total N non-ref muts in allele

  # ------ FN starts here
  # Subfunctions
  charoneal<-function(al, shpair.use, gvs){
    # Characterizes ONE allele
    # Input: one row of input als dt narrowed to columns displayname, tRNA, allelename, AA, Codon, AlleleCM, Infernal, IsotypeScore,
    #                  Classification, n.allele, n.called, n.missing, freq.ofallinclmissing, nAls
    #       , full other dts
    # Out: one-row data.table with all allele info from al input and columns:
    # is.reference, is this allele the reference
    # is.pair, this allele is one with at least one putatively comp. pair
    # n.pairs.in, count of N potentially compensatory pairs this allele harbors
    # has.pair.mut, T/F potential comp pair mutation (any!) seen in this allele
    # has.nonpair.mut, T/F other mutations that aren't in comp pair
    # n.muts, total N non-ref muts in allele
    
    alstrn.split<-al[, tstrsplit(allelename, "_")]
    alstrn<-unlist(alstrn.split[, ncol(alstrn.split) - 1, with = F])
    if(alstrn=="Reference"){ # don't need to process variants if it's ref strain
      is.ref<-T # is this allele ref
      is.pair<-F # this allele is one with comp pair
      n.pairs.in<-0 # count of N potentially compensatory pairs this allele harbors
      has.pair.mut<-F # T/F potential comp pair mutation (any!) seen in this allele
      has.nonpair.mut<-F # T/F other mutations that aren't in comp pair
      n.muts<-0 # total N non-ref muts in allele
    }else{ # process variant data
      is.ref<-F
      # Deal if it's in a pair
      strnin<-sapply(shpair.use$strains, function(x) alstrn%in%x)
      if(sum(strnin)>0){
        is.pair<-T
        n.pairs.in<-sum(strnin)
        has.pair.mut<-T # pair mut is in it or no?
      }else{
        # Deal if it's not in a pair
        is.pair<-F
        n.pairs.in<-0
      }
      # Find out total N vars in this allele ( whether they're in pair or not), record as n.muts
      almuts<-gvs[grepl(alstrn, homAlt)]
      n.muts<-nrow(almuts)
      ## set up for pair comparison
      almutinp<-unique(rbind(almuts[, .(pos,ref,alt)][shpair.use[, .(pos_1, ref_1, alt_1)], 
                                         on = .(pos = pos_1, ref = ref_1, alt = alt_1), nomatch = 0L],
                   almuts[, .(pos,ref,alt)][shpair.use[, .(pos_2, ref_2, alt_2)], 
                                            on = .(pos = pos_2, ref = ref_2, alt = alt_2), nomatch = 0L])) # merges on position & checks if there are rows in merge, checks both pair muts
      if(nrow(almutinp)>0){
        has.pair.mut<-T
        if(nrow(almuts)==nrow(almutinp)){
          has.nonpair.mut<-F
        }else if(nrow(almuts)>nrow(almutinp)){ # must have non-pair mutations
          has.nonpair.mut<-T
        }
      }else{
        has.pair.mut<-F
        has.nonpair.mut<-T # obligately has some
      }
    }
    
    # --- Format & return
    out<-data.table(al,
                    is.reference = as.logical(is.ref),
                    is.pair = as.logical(is.pair),
                    n.pairs.in = as.numeric(n.pairs.in),
                    has.pair.mut = as.logical(has.pair.mut),
                    has.nonpair.mut = as.logical(has.nonpair.mut),
                    n.muts = as.numeric(n.muts)
                    )
    return(out)
  }
  
  # --- Set up data
  # Narrow allele information to what I'll use here for ease TAKING AWAY WHAT IT IS KEYED ON THOUGH - FIX FOR LONGER TERM USE
  aluse<-als[, .( allelename, AA, Codon, AlleleCM, Infernal, IsotypeScore,
                 Classification, n.allele, n.called, n.missing, freq.ofallinclmissing, nAls)]
  # split strains as needed
  shpair.use<-data.table(shpair, strains = strsplit(shpair$whichInCommon, ","))
  
  # --- Process each allele
  altag<-aluse[, charoneal(.SD, shpair.use, gvs), by = .I]
  return(altag)
}

#### Arguments and inputs ####
p<-arg_parser("Plot/analyze tRNA gene body variation data (from get_strain_variants_relpos.py)", 
              name = "genebodyvariation.R", hide.opts = TRUE)
# Input files
p<-add_argument(p, "--speciesf",
                help = "File containing information on all species to process here. Columns infilename (exactly how all files have this species in their name),
                displayname (name that should be used for plot outputs etc), shortname (no-spaces name for ouptut files, sorting, etc - either shorter than or same as infilename, probably).
                In order you'd like plots to be in!",
                type = "character")
p<-add_argument(p, "--genevars",
                help = "EXAMPLE Path to file with tRNA gene body info output by get_strain_variants_relpos.py
                Where species ID/species specific info is, put SAMP instead",
                type = "character")
p<-add_argument(p, "--trnasecstruct",
                help = "Path to file with information on how tRNA sec structures in --genevars input are arranged, etc. Columns:
                structure (as in --genevars), substructure (as in --genevars), nSubStructure (as in --genevars), canonicalnbp - length of this structure in a 72/73 bp tRNA [for plotting],
                structure_plot and substructure_plot: categories & names prettified/simplified for how you'd like to plot them;
                structure_plot_level & substructure_plot_level: RANKINGS of unique structure_plot and substructure_plot for ordering in plot",
                type = "character")
p<-add_argument(p, "--trnainfo",
                help = "Path to *_genes_swfixedvarclassifs.txt tRNA gene-characterizing output of isotype switch variation script. For species combined",
                type = "character")
p<-add_argument(p, "--talleleinfo",
                help = "Path to tRNA allelic information (_alleles_swfixedvarclassifs.txt.gz)",
                type = "character")
p<-add_argument(p, "--ssfile",
                help = "Path to EXAMPLE file with tRNA secondary structures (as output by tRNAscan-SE). For all alleles INCLUDING ref.",
                type = "character")

# Output related
p<-add_argument(p, "--baseoutname",
                help = "Base name for all output files",
                default = "out")
p<-add_argument(p, "--outdir",
                help = "Output directory path. if getwd(), current will be used",
                default = "out")

# Parse arguments
cat("....Parsing arguments....\n")
p<-parse_args(p)

# outdir
if(p$outdir=="getwd()"){
  p$outdir<-getwd()
}
if(!dir.exists(p$outdir)){dir.create(p$outdir, recursive = T)}
setwd(p$outdir)

#### Read in and format data ####
cat("....Reading in data....\n")
sinfo<-fread(p$speciesf, header = T)
tstruct<-fread(p$trnasecstruct, header = T)

# var info
vdat.all<-fread.mult(sinfo = sinfo, exampfile = p$genevars, skip = 1)
## remove 'chr' if there
vdat.all[startsWith(tRNA, "chr"), tRNA:= sapply(tRNA, function(x) substr(x, 4, nchar(x)))]

# tRNA info
tinfo<-fread(p$trnainfo, header = T)
## remove 'chr' if there
tinfo[startsWith(tRNA, "chr"), tRNA:= sapply(tRNA, function(x) substr(x, 4, nchar(x)))]

# tRNA allele info
alinfo<-fread(p$talleleinfo, header = T)
## remove 'chr' if there
alinfo[startsWith(tRNA, "chr"), tRNA:= sapply(tRNA, function(x) substr(x, 4, nchar(x)))]

# tRNA secondary structures (not used til much later)
secs<-rbindlist(lapply(1:nrow(sinfo), function(i){
  data.table(displayname = sinfo[i, displayname],
    readss(ssfile = gsub("SAMP", sinfo[i, infilename], p$ssfile),
         restructuress = T))
}))
## remove 'chr' if there
secs[startsWith(name, "chr"), name:= sapply(name, function(x) substr(x, 4, nchar(x)))]


#### Look for paired mutations possibilities ####
cat("....Finding mutations in stem regions in the same strains....\n")
# --- Narrow to the stems of interest
mystems<-tstruct[structure%in%c("acceptor_stem", "d", "anticodon", "t") & substructure%in%c("arm_L", "arm_R")]
vdat<-vdat.all[structure%in%mystems$structure & substructure%in%mystems$substructure] 

# --- Keep cases where there is at least one variant in arm_L & arm_R (within same structure, species, and gene)
setkey(vdat, displayname, tRNA, structure)
# Count per group
arm_counts <- vdat[, .(
    n_L = sum(substructure == "arm_L", na.rm = TRUE),
    n_R = sum(substructure == "arm_R", na.rm = TRUE)
    ),
  by = .(displayname, tRNA, structure)
]

# Keep ones that have 1+ in each arm
ac_keep<-arm_counts[n_L>0 & n_R>0, .(displayname, tRNA, structure)]
setkey(ac_keep, displayname, tRNA, structure)
vposspair<-vdat[ac_keep][substructure%in%c("arm_L", "arm_R")]
# EXCLUDE INDELS
vposspair<-vposspair[nchar(ref)==1 & nchar(alt)==1]

# --- Narrow to keep cases where 2+ variants are in same STRAIN
# Get these
setkey(vposspair, displayname, tRNA, structure)
strnshare<-vposspair[, {
  L <- .SD[substructure=="arm_L"]
  R<- .SD[substructure=="arm_R"]
  
  if (nrow(L) == 0L || nrow(R) == 0L){
    NULL   # nothing to check in this group
  }else{
    
    # all L × R combinations (row indices)
    idx <- CJ(i = seq_len(nrow(L)), j = seq_len(nrow(R)), sorted = FALSE)
    
    # apply ckpair to each pair; this returns internally
    idx[ , ckpair(L[i], R[j]), ]
  }
},
  .SDcols = -which(names(vposspair)%in%c("homRef", "missingOrHet")),
  by = .(displayname, tRNA, structure)
]

# --- Integrate with tRNA info (so can have easy pseud filter, for example)
setkey(tinfo, displayname, tRNA)
setkey(strnshare, displayname, tRNA)
strnshare<-tinfo[strnshare]
# Save out
write.table(strnshare, file = file.path(p$outdir, paste0(p$baseoutname, "_strainSharingStemSNVPairs.txt")),
            sep = "\t", quote = F, row.names = F)

# Little summaries
## By tRNA classifications
setkey(strnshare, displayname, all.lost, gene.classif)
nbylost<-data.table(tested = "Pseud. vs. not (all lost alleles)", strnshare[, .N, by  = .(displayname, all.lost)])
nbyswcl<-data.table(tested = "Iso. sw. classification (of not all pseud)", strnshare[all.lost==F, .N, by = .(displayname,gene.classif)])

## By singleton SNP at one, both places?
strnshare[, nSingletonOfPair:=((nHomAlt_1 == 1) + (nHomAlt_2 == 1))]
setkey(strnshare, displayname, all.lost, nSingletonOfPair)
nbysing<-data.table(tested = "N SNVs in this pair that are singletons (split by pseud & not)",
                    strnshare[, .N, by  = .(displayname, all.lost, nSingletonOfPair)])

## By STRUCTURE
setkey(strnshare, displayname, all.lost, structure)
nbystruct<-data.table(tested = "tRNA structure (split by pseud & not)",
                      strnshare[, .N, by  = .(displayname, all.lost, structure)])

## by tRNA - not super useful, isn't counting tRNAs
setkey(strnshare, displayname, all.lost, tRNA)
nbygene<-data.table(tested = "tRNA gene (split by pseud & not)",
                      strnshare[, .N, by  = .(displayname, all.lost, tRNA)])

## N is tRNAs not pairs
ntrnas<-data.table(tested = "number of tRNA genes, split by pseud & not",
                   strnshare[,length(unique(tRNA)), by  = .(displayname, all.lost)])
setnames(ntrnas, "V1", "N_genes")

## Save out summaries
list.summ<-list(nbylost, nbyswcl, nbysing, nbystruct, nbygene, ntrnas)
names(list.summ)<-c("nbylost", "nbyswcl", "nbysing", "nbystruct", "nbygene", "ntrnas")
invisible(lapply(names(list.summ), function(x){
  write.table(list.summ[[x]], 
              file = file.path(p$outdir, paste0(p$baseoutname, "strainSharingStemSNVPossPairs_summary_", x, ".txt")),
              sep = "\t", quote = F, row.names = F)
}))

#### Get allele freq distribution of alleles with/without these pairs #### 
# for non pseud I think
# keep infernal info too
# classify allele by N mutations it has, too?

# thinking through...
# For each strain with the pair, what allele is it in
#       note that as a PAIR allele
#       note if it also has --other-- variants [go back to var dat]
# For each var in the pair, what OTHER alleles is it in (other strains not in the pair)
#     note that as a member 1 of pair allele, etc
#     note if it also has --other-- variants [go back to var dat]
# Note other alleles as other variants? [basically - these DO NOT have any pairs of interest. maybe count how many variants they do have.]

setkey(strnshare, displayname, tRNA)
setkey(vdat.all, displayname, tRNA)
setkey(alinfo, displayname, tRNA)
alpairtag<-strnshare[all.lost==F ,{
  shpair<- .SD
  gvs <- vdat.all[.BY, on = .(displayname, tRNA)]  
  als <- alinfo[.BY, on = .(displayname, tRNA)]
  # (optional) drop key columns from the other inputs to avoid duplicate-key columns downstream
  gvs <- gvs[, !c("displayname","tRNA"), with = FALSE]
  als<-als[ , !c("displayname","tRNA"), with = FALSE]
  
  charpairals(shpair, gvs, als)
} ,
by = .(displayname, tRNA)]

write.table(alpairtag, file.path(p$outdir, paste0(p$baseoutname, "_tRNAsNonPseudwstrainSharingStemSNVPairs_alleleinfo.txt")),
            sep = "\t", quote = F, row.names = F)

# After this, summarize/plot: highlight the ones that have the pair vs those that don't for example
## setting up

# how many category boundaries we could have at most (based on global tRNA levels)
n_lvls <- length(unique(alpairtag$tRNA))
boundary_df <- data.frame(
  xintercept = seq(1.5, n_lvls - 0.5, by = 1)
)

# actual plot
alplt<-ggplot(alpairtag) +
  ggforce::geom_sina(aes(tRNA, freq.ofallinclmissing, color = is.pair, shape = has.pair.mut), alpha = 0.4) +
  # add vertical boundary lines between categories
  geom_vline(
    data = boundary_df,
    aes(xintercept = xintercept),
    color = "grey80",
    linewidth = 0.4,
    inherit.aes = FALSE
  ) +
  scale_color_manual(values = c("FALSE" = "black", "TRUE" = "red")) + 
  xlab("") + ylab("Allele frequency") +
  labs(color = "Allele has potential\ncompensatory pair", shape = "Allele carries\na pair mutation") +
  facet_wrap(~displayname, scales = "free_x") +
  myggtheme + theme(legend.position = "bottom", axis.text.x = element_text(angle = 90), panel.grid = element_blank())

pdf(file.path(p$outdir, paste0(p$baseoutname, "_tRNAsNonPseudwstrainSharingStemSNVPairs_alleleinfo_dotplot.pdf")), 8, 5)
print(alplt)
dev.off()

#### Script completion message & session information ####
cat("....investigate_compmuts.R processing complete! Session information:....\n")
sessionInfo()