#! /usr/bin/env/ Rscript
# Begin tRNA gene/allele tree explorations
# By Avery Davis Bell, begun 2025.04.29
require(ggtree) # from Bioconductor
require(tidytree)
require(ape)
require(phangorn)
require(data.table)
require(argparser)
require(ggplot2)
require(ggnewscale) # lets do multiple scales
require(RColorBrewer)

#### Functions ####
xfasta2fasta<-function(xfastaf, ofilepath){
  # Writes info in xfasta out to normal fasta
  # In: xfastaf, path to x fasta
  #     ofilepath, file to write to

  # --- Subfunctions (borrowed from other scripts)
  readxfasta<-function(xfastaf){
    # Reads an X fasta into a data.table with columns name, sequence, secstruct
    
    mylines<-readLines(xfastaf)
    out<-data.table(name = sapply(strsplit(mylines[seq(1, length(mylines), 3)], ">"), function(x) x[2]),
                    sequence = mylines[seq(2, length(mylines), 3)],
                    secstruct = mylines[seq(3, length(mylines), 3)])
    return(out)
  }
  
  writefasta<-function(ssdt, ofilepath){
    # Writes seqs from readss out in fasta format
    # In: ssdt, data.table with columns name, sequence [at least]
    #     ofilepath, file to write to
    # Out: none
    
    nlines<-ssdt[, paste0(">", name)]
    seqlines<-ssdt[, sequence]
    wlines<-c()
    for(i in 1:length(nlines)){
      wlines<-c(wlines, nlines[i], seqlines[i])
    }
    writeLines(wlines, con = ofilepath, sep = "\n")
  }
  
  # --- Read & write
  dat<-readxfasta(xfastaf)
  writefasta(dat, ofilepath)
}

njtree<-function(dnaobj, mlmod="JC69"){
  # Builds NJ tree object, doesn't save intermediates
  # In: dnaobj to build tree from - from ape read.dna. All seqs need to be same length.
  #     mlmod, passed to model argument of dist.ml
  # Out: NJ tree object - NJ(dist.ml(phydat(init obj)
  
  phy<-phyDat(dnaobj, type = "DNA", levels = NULL)
  mldist<-dist.ml(phy, model = mlmod)
  njout<-NJ(mldist)
  
  return(njout)
  
  # # Build the tree object - Annalise's code\
  # trnas_phydat <- phyDat(trnas, type="DNA", levels=NULL) # from phangorn pkg
  # trnas_dist <- dist.ml(trnas_phydat, model="JC69") # in phangorn
  # trnas_NJ <- NJ(trnas_dist)
  # 
  # trnas_NJ
}

addtreeinfo<-function(infodt, tree, labinfocol = "allelename"){
  # Annotates tree with info from infodt **For the subset in tree**
  # In: infodt, data.table with all info to add
  #     tree, tree to annotate
  #     labinfocol, name of column in infodt that corresponds to tip labels in tree
  # Out: tree but with info annotated
  
  myinfo<-copy(infodt)
  setkeyv(myinfo, labinfocol)
  tmplabel<-data.table(label = tree$tip.label, 
                       myinfo[tree$tip.label, ])
  tree.out<-full_join(tree, tmplabel, by = "label")
  return(tree.out)
}

plotmytree<-function(tree.lab, aacols, multspecies = T, spcols, clado = F, aacolbackbone = T, alfreq = T, aalty = "dotted",
                     classifcols, offset.clad = 0.4, offset.br = 0.2, labtxtsize = 0.5,
                     treelinesz = 0.1, mylinesz = 0.1,
                     labelantic = F){
  # Plots tRNA tree with all the nice stuff. I'm not going to write all the req'd columns here
  # In: tree.lab, tree with all the columns/labels of interest
  #     aacols, colors for each amino acid 
  #     multspecies, if T, colors by species 
  #     spcols, named vector of colors to use for species
  #     clado, if T, does cladogram (no branch lengths)
  #     aacolbackbone, if T, colors by AlleleCM.orig - original what gene has; if F, colors by the amino acid 'grabbed'
  #     alfreq, if T, plots tip points scaled by allele freq
  #     aalty, line type for connecting tips to labels
  #     classifcols, color for values in labclassif column
  #     offset.clad, offset starter for labels if cladogram
  #     offset.br, offset starter for labels if branch length not none
  #     labtxtsize, size for text labels
  #     treelinesz, lwd for tree
  #     mylinesz, width of line connecting branch to outside
  #     labelantic, if TRUE, labels with AA | anticodon instead of just AA
  # Out: ggplot
  
  # Set up skeleton
  if(clado==T){
    plt<-ggtree(tree.lab, layout = "circular", lwd = treelinesz, branch.length = "none")
    myoff<-offset.clad
    myaaoff<-offset.clad
  }else{
    plt<-ggtree(tree.lab, layout = "circular", lwd = treelinesz) + xlim(-0.5, NA)
    myoff<-offset.br
    myaaoff<-0
  }
  
  # figure out aa label
  if(labelantic==F){
    aalab<-expr(AA)
  }else{
    aalab<-expr(paste(AA, Codon, sep = " | "))
  }
  
  # Add colors/labels for the amino acid (coded for or backbone)
  if(aacolbackbone==T){
    plt<-plt + geom_tiplab(aes(color = AlleleCM, label = eval(aalab)), # right now doing AlleleCM not AlleleCM.orig - the UNDET super variable stuff messes with that a lot
                           align = T, show.legend = F, size = labtxtsize, 
                           linetype = aalty,
                           linesize = mylinesz, 
                           offset = myaaoff,
                           geom = "text") +
      scale_color_manual(values = aacols) +
      new_scale_color()
  }else{
    plt<-plt + geom_tiplab(aes(color = AA, label = eval(aalab)),# FOR WITHIN ONE AlleleCM/amino acid, should color by WHAT IT CODES FOR instead: color & label are AA, because can say they're all Leus
                           align = T, show.legend = F, size = labtxtsize,
                           linetype = aalty,
                           linesize = mylinesz, 
                           offset = myaaoff,
                           geom = "text") +
      scale_color_manual(values = aacols) +
      new_scale_color()
  }
  # Add tip size for allele freq
  if(alfreq==T){
    plt<-plt + geom_tippoint(aes(size = freq.ofallinclmissing), alpha = 0.1) +
      labs(size = "Allele frequency\nin species")
  }
  
  # Add labels for classification
  plt<-plt + geom_tiplab(aes(color = labclassif, label = labclassif), # right now this changes the dotted colors...
                         align = T, show.legend = F, linetype=NULL, # linetype NULL keeps it from overwriting that
                         size = labtxtsize,
                         offset = myoff + myaaoff, # offset to get it to show up too; may depend on xlim stuff 
                         geom = "text") +
    scale_color_manual(values = classifcols)
  
  
  # Add species colors tip labs
  if(multspecies==T){
      # Other idea: give 3 circles, T/F if species is there
      plt<-plt + geom_tiplab(aes(label = ele), color = ele.col,
                             align = T, show.legend = F, linetype=NULL,
                             geom = "text",
                             offset = myaaoff + myoff*3,
                             size = labtxtsize*4) +
        geom_tiplab(aes(label = bri), color = bri.col,
                    align = T, show.legend = F, linetype=NULL,
                    geom = "text",
                    offset = myaaoff + myoff*3 + myoff/3,
                    size = labtxtsize*4) +
        geom_tiplab(aes(label = trop), color = trop.col,
                    align = T, show.legend = F, linetype=NULL,
                    geom = "text",
                    offset = myaaoff + myoff*3 + (myoff/3)*2,
                    size = labtxtsize*4) 
      # OK this works IF set label in data.table as you want and give it one color (instead of doing scale)
      # Internal to here - fix offset size so they're right next to each other, offset from other labels
      
      # UPDATE OFFSET HERE to pass to next label level when this works [want classification outside this if it exists but not spaced by it if not]
    }
  
 
  #...add tRNA gene name? other?
  
  # Species should probably be label out at end too...allele freq maybe more OK to leave as point lab....
  # ***OR HIGHLIGHT, TRY THAT **** [could be diff grayscale colors.....]....have to do based on NODES....
  # tree.tib<-as_tibble(tree.lab)
  # plt + geom_highlight(data = tree.tib, aes(node = node, fill = displayname)) + scale_fill_manual(values = spcols)
  #     [might want to rename the geom_hilight data]
  #     SUPER SLOW, let's see if it can work at all
  # Works BUT most nodes aren't assoc with a species because they're internal, so that doesn't work
  
  # Return
  return(plt)
}


#### Arguments & inputs ####
p<-arg_parser("Begin tRNA gene/allele tree explorations", 
              name = "trna_trees_init.R", hide.opts = TRUE)

# Inputs
p<-add_argument(p, "--speciesf",
                help = "File containing information on all species processed in preceeding script. Columns infilename (exactly how all files have this species in their name),
                displayname (name that should be used for plot outputs etc - in other input files), shortname (no-spaces name for ouptut files, sorting, etc - either shorter than or same as infilename, probably).
                In order you'd like plots to be in!",
                type = "character")
p<-add_argument(p, "--alleleplus",
                help = "*alleleinfo_counts_worigcodonetc.txt  output of show_mutational_variation.R, with all allele info including Codon.orig, AlleleCM.orig, VariableInPop, Classification, etc",
                default = "~/GaTech Dropbox/Avery Bell/PaabyLab/Projects/tRNA/tRNAVarsFromVCFs/InitAlleleandLocPlots/threespecies_early2024CaeNDR_alleleinfo_counts_worigcodonetc.txt")
p<-add_argument(p, "--pseudref",
                help = "*_pseudoinref_genes.txt output of show_mutational_variation.R, with info on if each gene is all pseud",
                default = "~/GaTech Dropbox/Avery Bell/PaabyLab/Projects/tRNA/tRNAVarsFromVCFs/InitAlleleandLocPlots/threespecies_early2024CaeNDR_pseudoinref_genes.txt")
p<-add_argument(p, "--xfasta",
                help = "Foursale alignment xfasta for ALL tRNA alleles in ALL species of interest",
                default = "~/GaTech Dropbox/Avery Bell/PaabyLab/Projects/tRNA/tRNAVarsFromVCFs/DataFromPACE/tRNAAlignmentsSecStruct/fullalignments/allalleles_3species_202312and202401caendr_foursale_aligned.xfasta.gz")

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

#### Read in data ####
cat("....Reading in data....\n")
sinfo<-fread(p$speciesf, header = T) # not sure I need this?
# --- Allele info
info<-fread(p$alleleplus)
# Add if they're always pseud to each allele
psinf<-fread(p$pseudref)
setkey(info, tRNA)
setkey(psinf, tRNA)
ps.annot<-psinf[n_strains_func==0, .(displayname, tRNA)]
info[paste(displayname, tRNA)%in%ps.annot[, paste(displayname, tRNA)], all.pseud:=TRUE]
info[is.na(all.pseud), all.pseud:=F]

# Add species to front of allelename to match with xfasta. ***NOT INPUT BASED; HARDCODED....should change ideally
info[displayname=="C. briggsae", allelename:=paste("cbrig", allelename, sep = "_")]
info[displayname=="C. elegans", allelename:=paste("cele", allelename, sep = "_")]
info[displayname=="C. tropicalis", allelename:=paste("ctrop", allelename, sep = "_")]
info[, allelename:=paste(allelename, "trna1", sep = ".")] # all in xfasta have .trna1 after

# --- Sequence/alignment info
# re-format as fasta - did not work as x fasta
tmpfa<-file.path(p$outdir, paste0("tmp", sample.int(1e04, 1), ".fasta"))
xfasta2fasta(p$xfasta, tmpfa)
# ape to read.dna
talls<-read.FASTA(tmpfa, type = "DNA") 
# Clean up
file.remove(tmpfa)

# order info in terms of talls
info<-info[match(names(talls), allelename)]

#### Generate tree objects [before plotting] ####
cat("....Generating tree objects....\n")
# --- All species, generally all genes
allnj<-njtree(talls) # Takes some time! Lots of seqs!
allnonpseud<-njtree(talls[info[all.pseud==F, allelename]]) # exclude any with only pseud alleles.
save.image(file.path(p$outdir, paste0(p$baseoutname, "_trees.RData"))) ## for ease of resuscitation

# Just reference alleles [quickest way to just include one copy per gene]
refnj<-njtree(talls[info[grepl("Reference", allelename), allelename]])
save.image(file.path(p$outdir, paste0(p$baseoutname, "_trees.RData"))) ## for ease of resuscitation
refnonpseud<-njtree(talls[info[grepl("Reference", allelename) & all.pseud==F & Classification!="Lost", allelename]])
save.image(file.path(p$outdir, paste0(p$baseoutname, "_trees.RData"))) ## for ease of resuscitation

# --- Individual species
specnjs<-lapply(1:nrow(sinfo), function(spind){
  # All alleles
  allals<-njtree(talls[info[displayname==sinfo[spind, displayname], allelename]])
  # All not all pseud
  nonpseudals<-njtree(talls[info[displayname==sinfo[spind, displayname] & all.pseud==F, allelename]])
  # Ref alleles
  refnj<-njtree(talls[info[displayname==sinfo[spind, displayname] & grepl("Reference", allelename), allelename]])
  # Ref alleles not pseud
  refnonpseud<-njtree(talls[info[displayname==sinfo[spind, displayname] & grepl("Reference", allelename) & all.pseud==F & Classification!="Lost", allelename]])
  
  return(list(allals = allals, nonpseudals = nonpseudals, refals = refnj, refnonpseud = refnonpseud))
})
names(specnjs)<-sinfo$displayname
save.image(file.path(p$outdir, paste0(p$baseoutname, "_trees.RData"))) ## for ease of resuscitation

# --- Individual aas - mostly focusing on ones where isotype switching looks interesting. These are QUICK!
# Leu - all with ORIGINAL coded as leu (and not all pseud)
leunj<-njtree(talls[info[AlleleCM.orig=="Leu" & all.pseud==F, allelename]])
save.image(file.path(p$outdir, paste0(p$baseoutname, "_trees.RData"))) ## for ease of resuscitation

# Ser - all with ORIGINAL coded as Ser (and not all pseud)
sernj<-njtree(talls[info[AlleleCM.orig=="Ser" & all.pseud==F, allelename]])
save.image(file.path(p$outdir, paste0(p$baseoutname, "_trees.RData"))) ## for ease of resuscitation

# Thr - all with ORIGINAL coded as Ser (and not all pseud)
thrnj<-njtree(talls[info[AlleleCM.orig=="Thr" & all.pseud==F, allelename]])
save.image(file.path(p$outdir, paste0(p$baseoutname, "_trees.RData"))) ## for ease of resuscitation
           
# Trp - all with ORIGINAL coded as Ser (and not all pseud)
trpnj<-njtree(talls[info[AlleleCM.orig=="Trp" & all.pseud==F, allelename]])
save.image(file.path(p$outdir, paste0(p$baseoutname, "_trees.RData"))) ## for ease of resuscitation

## NEW: all individual amino acids. Based on AlleleCM.orig. NOT for all pseuds
aas.do<-info[, sort(unique(AlleleCM.orig))]
aas.do<-as.character(aas.do[aas.do!="Undet"])
aanjs<-lapply(aas.do, function(aa){
  return(njtree(talls[info[AlleleCM.orig==aa & all.pseud==F, allelename]]))
})
names(aanjs)<-aas.do
save.image(file.path(p$outdir, paste0(p$baseoutname, "_trees.RData"))) ## for ease of resuscitation


#### Make some tree plots ####
cat("....Making tree plots....\n")
# --- Set up color scales
# Amino acids
aas<-c("Ala","Arg","Asn","Asp","Cys","Gln","Glu","Gly","His","Ile","iMet","Leu","Lys","Met","Phe","Pro","SeC","Ser","Thr","Trp","Tyr","Val","Undet")
aacols <- c("#8B0000","#FF0010","#FF5005","#FFA405","#FFE100","#005C31","#00998F","#2BCE48","#94FFB5","#00FFFF","#AFEEEE","#00BFFF","#87CEFA","#3283FE","#0000FF","#00008B","#AA0DFE","#DEA0FD","#FE00FA","#FFC0CB","#FF69B4","#FF1493","#D3D3D3")
names(aacols)<-aas
# Species
spcols<-c("black", "gray40", "gray80") # brewer.pal(nrow(sinfo), "Set2")
names(spcols)<-sinfo$displayname

# Mismatch/lost classifications?? should add to info as well
info[Classification=="Altered", labclassif:="Mismatch"]
info[Classification=="Lost", labclassif:="Lost"]
# *** UPDATING so that any with Undet AA is also called Lost - isn't necessarily at least in first iteration of input
info[AA=="Undet", labclassif:="Lost"]

classifcols<-c("red", "orange")
names(classifcols)<-c("Lost", "Mismatch")

# --- Merge in labels with tree info 
# FACTOR all things that should be factors
info[, AlleleCM:=factor(AlleleCM, levels = names(aacols))]
info[, AlleleCM.orig:=factor(AlleleCM.orig, levels = names(aacols))]
info[, AA:=factor(AA, levels = names(aacols))]
info[, labclassif:=factor(labclassif, levels = names(classifcols))]
info[, displayname:=factor(displayname, levels = names(spcols))]
# > info[displayname=="C. elegans", spshort:="|\n\n"]
# > info[displayname=="C. briggsae", spshort:="\n|\n"]
# > info[displayname=="C. tropicalis", spshort:="\n\n|"]
# > info[, spshort:=factor(spshort, levels = c("|\n\n", "\n|\n", "\n\n|"))]
info[, ele:=factor(ifelse(displayname=="C. elegans", "-", ""), levels = c("-", ""))] # factor part is critical
info[, bri:=factor(ifelse(displayname=="C. briggsae", "-", ""), levels = c("-", ""))]
info[, trop:=factor(ifelse(displayname=="C. tropicalis", "-", ""), levels = c("-", ""))]
ele.col<-spcols[1]
bri.col<-spcols[2]
trop.col<-spcols[3]
# ele.col<-c("#00000000", spcols[1])
# names(ele.col)[1]<-NA
# bri.col<-c("#00000000",spcols[2])
# names(bri.col)[1]<-NA
# trop.col<-c("#00000000",spcols[3])
# names(trop.col)[1]<-NA

# Merge in relevant subsets
## All together
allnj.lab<-addtreeinfo(info, allnj)
allnonpseud.lab<-addtreeinfo(info, allnonpseud)
refnj.lab<-addtreeinfo(info, refnj)
refnonpseud.lab<-addtreeinfo(info, refnonpseud)
## Species specific
specnjs.lab<-lapply(specnjs, function(onespec){
  allals<-addtreeinfo(info, onespec$allals)
  nonpseudals<-addtreeinfo(info, onespec$nonpseudals)
  refnj<-addtreeinfo(info, onespec$refals)
  refnonpseud<-addtreeinfo(info, onespec$refnonpseud)
  
  return(list(allals = allals, nonpseudals = nonpseudals, refals = refnj, refnonpseud = refnonpseud))
})
names(specnjs.lab)<-names(specnjs)

## Specific amino acids
leunj.lab<-addtreeinfo(info, leunj)
sernj.lab<-addtreeinfo(info, sernj)
thrnj.lab<-addtreeinfo(info, thrnj)
trpnj.lab<-addtreeinfo(info, trpnj)

## ALL
aanjs.lab<-lapply(aanjs, function(oneaa){
  return(addtreeinfo(info, oneaa))
})
names(aanjs.lab)<-names(aanjs)

# SAVE image
save.image(file.path(p$outdir, paste0(p$baseoutname, "_trees.RData"))) ## for ease of resuscitation

# --- Make plots that have all 3 species, all alleles (or without pseud at all alleles)
# Remember - AA: what it codes for. AlleleCM: what BACKBONE grabs. AlleleCM.orig: what backbone of MOST ALLELES grabs. ***

# --- Make plots with all species together
#    for all alleles, all without pseud at all alleles
#    with each allele labeled by its AlleleCM backbone and also version labeled by what its anticodon is
#    phylo, clado versions
aspinf<-data.table(dat = rep(c("allnj.lab", "allnonpseud.lab", "refnj.lab", "refnonpseud.lab"), each = 2),
                   datdescrip = rep(c("All alleles", "Alleles from all tRNA genes with at least one functional allele",
                                      "Reference strain alleles", "Reference strain alleles (genes)"), each = 2),
                   backbone = rep(c(T, F), 4),
                   backbonedescrip = rep(c("AA color: gene backbone", "AA color: anticodon"), 4))
# Phylos
pdf(file.path(p$outdir, paste0(p$baseoutname, "_phylotrees_allspeciesgenes.pdf")), 16, 16)
invisible(
lapply(1:nrow(aspinf), function(i){
  print(
    plotmytree(tree.lab = eval(as.name(aspinf[i, dat])), aacols, multspecies = T, spcols, clado = F,
               aacolbackbone = as.logical(aspinf[i, backbone]), alfreq = T, aalty = "solid",
               classifcols, offset.clad = NA, offset.br = 0.0125, labtxtsize = 0.4, 
               treelinesz = 0.1, mylinesz = 0.1) +
      ggtitle(aspinf[i, datdescrip], subtitle = aspinf[i, backbonedescrip])
  )
})
)
invisible(dev.off())

# Clados
pdf(file.path(p$outdir, paste0(p$baseoutname, "_cladotrees_allspeciesgenes.pdf")), 16, 16)
invisible(
  lapply(1:nrow(aspinf), function(i){
    print(
      plotmytree(tree.lab = eval(as.name(aspinf[i, dat])), aacols, multspecies = T, spcols, clado = T,
                 aacolbackbone = as.logical(aspinf[i, backbone]), alfreq = F, aalty = "solid",
                 classifcols, offset.clad = 2, labtxtsize = 0.4, 
                 treelinesz = 0.1, mylinesz = 0.1) +
        ggtitle(aspinf[i, datdescrip], subtitle = aspinf[i, backbonedescrip])
    )
  })
)
invisible(dev.off())

# --- Make individual species plots for all alleles, all without pseud at all alleles
#     with each allele labeled by its AlleleCM backbone and also version labeled by what its anticodon is
#     phylo, clado versions
spinf<-data.table(dat = rep(c("allals", "nonpseudals", "refals", "refnonpseud"), each = 2),
                  datdescrip = rep(c("All alleles", "Alleles from all tRNA genes with at least one functional allele",
                                     "Reference strain alleles", "Reference strain alleles (genes)"), each = 2),
                  backbone = rep(c(T, F), 4),
                  backbonedescrip = rep(c("AA color: gene backbone", "AA color: anticodon"), 4))

lapply(1:nrow(sinfo), function(spind){
  # Phylos
  pdf(file.path(p$outdir, paste0(p$baseoutname, "_phylotrees_", sinfo[spind, shortname],".pdf")), 16, 16)
  lapply(1:nrow(spinf), function(i){
    print(
      plotmytree(tree.lab = specnjs.lab[[sinfo[spind, displayname]]][[spinf[i, dat]]], aacols, multspecies = F, spcols, clado = F,
                 aacolbackbone = as.logical(spinf[i, backbone]), alfreq = T, aalty = "solid",
                 classifcols, offset.clad = NA, offset.br = 0.0125, labtxtsize = 0.4, 
                 treelinesz = 0.1, mylinesz = 0.1) +
        ggtitle(paste(sinfo[spind, displayname], "|", aspinf[i, datdescrip]), subtitle = aspinf[i, backbonedescrip])
    )
  })
  invisible(dev.off())
  
  # Clados
  pdf(file.path(p$outdir, paste0(p$baseoutname, "_cladotrees_", sinfo[spind, shortname],".pdf")), 16, 16)
  lapply(1:nrow(spinf), function(i){
    print(
      plotmytree(tree.lab = specnjs.lab[[sinfo[spind, displayname]]][[spinf[i, dat]]], aacols, multspecies = F, spcols, clado = T,
                 aacolbackbone = as.logical(spinf[i, backbone]), alfreq = F, aalty = "solid",
                 classifcols, offset.clad =2, labtxtsize = 0.4, 
                 treelinesz = 0.1, mylinesz = 0.1) +
        ggtitle(paste(sinfo[spind, displayname], "|", aspinf[i, datdescrip]), subtitle = aspinf[i, backbonedescrip])
    )
  })
  invisible(dev.off())
})

# --- Make individual AA plots - first few
aapinf<-data.table(dat = c("leunj.lab", "sernj.lab", "thrnj.lab", "trpnj.lab"),
                   datdescrip = c("Leucine backbone", "Serine backbone", "Threonine backbone", "Tryptophan backbone"),
                   datname = c("leu", "ser", "thr", "trp"))

## Phylo & clado in same PDF here
invisible(lapply(1:nrow(aapinf), function(i){
  pdf(file.path(p$outdir, paste0(p$baseoutname, "_trees_", aapinf[i, datname],".pdf")), 16, 16)
  # Phylo
  print(
    plotmytree(tree.lab = eval(as.name(aapinf[i, dat])), aacols, multspecies = T, spcols, clado = F,
               aacolbackbone = F, alfreq = T, aalty = "solid",
               classifcols, offset.clad = NA, offset.br = 0.05, labtxtsize = 2, 
               treelinesz = 0.5, mylinesz = 0.5) +
      ggtitle(aapinf[i, datdescrip], subtitle = "Alleles colored by their anticodon")
  )
  
  # Clado
  print(
    plotmytree(tree.lab = eval(as.name(aapinf[i, dat])), aacols, multspecies = T, spcols, clado = T,
               aacolbackbone = F, alfreq = T, aalty = "solid",
               classifcols, offset.clad = 2, labtxtsize = 2, 
               treelinesz = 0.5, mylinesz = 0.5) +
      ggtitle(aapinf[i, datdescrip], subtitle = "Alleles colored by their anticodon")
  )
  
  invisible(dev.off())
}))
# confirm want/don't want allele freq...
  
# --- Make individual AA plots: ALL
# phylos
pdf(file.path(p$outdir, paste0(p$baseoutname, "_phylotrees_eachAA.pdf")), 16, 16)
invisible(
lapply(names(aanjs.lab), function(oneaa){
  print(
    plotmytree(tree.lab = aanjs.lab[[oneaa]], aacols, multspecies = T, spcols, clado = F,
               aacolbackbone = F, alfreq = T, aalty = "solid",
               classifcols, offset.clad = NA, offset.br = 0.09, labtxtsize = 1.8, 
                 treelinesz = 0.3, mylinesz = 0.3, labelantic = T) +
      ggtitle(oneaa, subtitle = "Alleles colored by their anticodon-matching amino acid")
  )
})
)
invisible(dev.off())

# Clados
pdf(file.path(p$outdir, paste0(p$baseoutname, "_cladotrees_eachAA.pdf")), 16, 16)
invisible(
  lapply(names(aanjs.lab), function(oneaa){
    print(
      plotmytree(tree.lab = aanjs.lab[[oneaa]], aacols, multspecies = T, spcols, clado = T,
                 aacolbackbone = F, alfreq = T, aalty = "solid",
                 classifcols, offset.clad = 6, offset.br = NA, labtxtsize = 1.8, 
                 treelinesz = 0.3, mylinesz = 0.3, labelantic = T) +
        ggtitle(oneaa, subtitle = "Alleles colored by their anticodon-matching amino acid")
    )
  })
)
invisible(dev.off())

#### Save specific annotated trees for ease of use in other scripts...? ####
# DO this if want to

#### Script completion message & session information ####
cat("....trna_trees_init.R processing complete! Session information:....\n")
sessionInfo()