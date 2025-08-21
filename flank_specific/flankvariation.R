#! /usr/bin/env/ Rscript
# Plot/analyze tRNA flank regions variation data (from trnaflankvcfvars.nf)
# By Avery Davis Bell, begun 2025.02.25
require(argparser, quietly = T)
require(data.table, quietly = T)
require(ggplot2, quietly = T)
require(RColorBrewer, quietly = T)

#### Functions ####
# --- Data wrangling
dtpci<-function(x, n, suff = ""){
  # gets prop and 95% CIs as data.table
  if(n<=0){
    out<-data.table(p = NaN, low95ci = NaN, high95ci = NaN)
  }else{
    res<-binom.test(x, n)
    out<-data.table(p = res$estimate, low95ci = res$conf.int[1], high95ci = res$conf.int[2])
  }
  setnames(out, paste0(names(out), suff))
  return(out)
}

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
        dat1<-fread(gsub("SAMP", sinfo[s, infilename], exampfile), header = T, ...)
      }else{
        dat1<-fread(gsub("SAMP", sinfo[s, infilename], exampfile), header = T, select = myselect, ...)
      }
      dat1<-data.table(dat1, sinfo[s, .(displayname, shortname)])
      setcolorder(dat1, c("displayname", "shortname"))
    })
  )
  return(dat)
}

freadgetpertrna<-function(trnainfoalf, examptrnag, sinfo, stripchrpref = T){
  # Reads in file containing tRNA per allele information including original codon & allele as well as (potential) remolded one. - *alleleinfo_counts_worigcodonetc.txt output of show_mutational_variation.R
  # Gets relevant information (see output description) in per-tRNA [format for each species]
  # In: info, data.table with columns infilename (exactly how all files have this species in their name), displayname, shortname
  #     trnainfoalf, path to file containing tRNA per allele information including original codon & allele as well as (potential) remolded one. - *alleleinfo_counts_worigcodonetc.txt output of show_mutational_variation.R
  #     examptrnag, EXAMPLE Path to file containing tRNA per gene information including name & strand. - *strain_trnas_info.txt output of build_alt_sequences.py. Where species ID/species specific info is, put SAMP instead
  #     stripchrpref, T or F: strip "chr" if it occurs at beginning of tRNA name? [it is NOT in flank data]
  # Out: data.table with one row per tRNA per strain. Columns:
  ##  displayname, species display name
  ##  shortname, species short name
  ##  tRNA, gene name
  ##  Strand, gene strand
  ##  Codon.orig, 'original' codon (in ref seq, usually - corresponding with AlleleCM.orig) this gene encodes
  ##  AlleleCM.orig, 'original' AA (in ref seq, usually - corresponding with AlleleCM.orig) this gene encodes
  ##  all.altered, are all alleles at this gene altered/isotype switch?
  ##  any.lost, are any alleles at this gene lost/pseudogenized?
  ##  all.lost, are all alleles at this gene lost/pseudogenized?
  ##  Codon.coded, actual codon in this gene (possibly multiple if SOME but not ALL alleles are isotype switches)
  ##  AlCM.coded, actual AA coded for by this gene (possibly multiple if SOME but not ALL alleles are isotype switches)
  
  # Strand info from multi-species files
  stdat<-unique(fread.mult(sinfo, exampfile = examptrnag, myselect = c("tRNA", "Strand")))
  setkey(stdat, displayname, tRNA)
  
  # Info from alleles
  aldat<-fread(trnainfoalf)
  setkey(aldat, displayname, tRNA)
  
  tdat<-aldat[, .(Codon.orig = Codon.orig[1], AlleleCM.orig = AlleleCM.orig[1],
                  all.altered = unique(all.altered), any.lost = "Lost"%in%Classification, all.lost = all(Classification=="Lost"),
                  Codon.coded = paste(unique(Codon), collapse = ","), AA.coded = paste(unique(AA), collapse = ",")),
              by = .(displayname, tRNA)]
  
  # Return
  odat<-stdat[tdat]
  if(stripchrpref==T){
    odat[startsWith(tRNA, "chr"), tRNA:=substr(tRNA, 4, nchar(tRNA))]
  }
  setkey(odat, displayname, tRNA)
  return(odat)
}

calcflankvar<-function(onef, onet, annotdt = NULL, mutsplit,
                       lflank = c(-40, 0), rflank = c(0, 40)){
  # Calculates various nucleotide variation metrics for one set of tRNAs (i.e., all non-pseud tRNAs for a given species)
  # In: onef, tRNA flank variants for set of interest.
  #           Required columns: tRNA, relPos, Ref, Alt,nHomRef, nHomAlt, nNotMissingHet, nmissingOrHet, majorAl, minorAl, n.majorAl, n.minorAl,
  #           minorAlFreq
  #           <maybe the homRef etc columns with strain IDs in there>
  #     onet, information on ALL tRNAs (with or without variation) for this set of interest [freadgetpertrna output format]
  #     annotdt, optional data.table to prepend to each row of output - i.e., describing this set of tRNAs
  #     mutsplit, data.table with specific mutations (from nt, to nt) to break out calculations by. Function also does 'other' - any not captured here - and 'all.SNV' (all THESE combined)  and 'any' (all ANY combined)
  #           Required columns: from (nt), to (nt), mutname (from > to format!!)
  #     lflank, bounds of 5' flank to compute values within (what input spans) - c(min, max)
  #     rflank, bounds of 3' flank to compute values within (what input spans) - c(min, max)
  # Out: list of results:
  #     $pvartrna, proportion of tRNAs with variant at each location in flanks (no sliding windows; can do as average post-hoc for proportion).
  #           One row per mutation class *per position in tRNA flank*. Columns:
  #           <any in annotdt>
  #           ## varfreq, this is done for all variants (subset as other columns suggest) AND for only vars with allele freq < 0.05
  #           ## mutclass, all, any.SNV, other, or specific base change
  #           ##. alleles, - for multi-mutation classes, Ref > Alt or Major > Minor for specific mutations (summarize for both of these, sometimes will be same and sometimes different)
  #           ##.  relpos, relative position in flank
  #           ##   flank, 3' or 5'
  #           ##   n.tRNAs.var, # tRNAs with variation at this position
  #           ##   n.tRNAs.invar, # tRNAs withOUT variation at this position. NB missing not taken into account here - it's on a per tRNA not allele freq/per strain basis
  #           ##   p.tRNAs.var, proportion of tRNAs (missing excluded) that have variant at this position
  #           ##   low95ci.tRNAs.var, binomial 95% CI on proportion here - lower bound
  #           ##   high05ci.tRNAs.var, binomial 95% CI on proportion here - upper bound
  
  # SUBFUNCTIONS (operating on one mutation class/set)
  proptvar<-function(oneclass, onet, lflank, rflank){
    # Calculate proportion of tRNAs with variant at each location in flanks
    # In: oneclass, flank variants to process - narrowed to set of interest (all summarized together here). Must have column Strand
    #     onet, all tRNAs to compare this set against (denominator, basically)
    #     lflank, bounds of 5' flank to compute values within (what input spans) - c(min, max)
    #     rflank, bounds of 3' flank to compute values within (what input spans) - c(min, max)
    # Out: data.table with one row per position in lflank and rflank analyzed. Columns:
    ##.  relpos, relative position in flank
    ##   flank, 3' or 5'
    ##   n.tRNAs.var, # tRNAs with variation at this position
    ##   n.tRNAs.invar, # tRNAs withOUT variation at this position. NB missing not taken into account here - it's on a per tRNA not allele freq/per strain basis
    ##   p.tRNAs.var, proportion of tRNAs (missing excluded) that have variant at this position
    ##   low95ci.tRNAs.var, binomial 95% CI on proportion here - lower bound
    ##   high05ci.tRNAs.var, binomial 95% CI on proportion here - upper bound
    
    # Get each site
    out<-data.table(relpos = c(lflank[1]:lflank[2], rflank[1]:rflank[2]),
                    flank = c(rep("5'", length(lflank[1]:lflank[2])), # need flank info because 0 could otherwise be either
                              rep("3'", length(rflank[1]:rflank[2]))))
    # Get # tRNAs with mutations at EACH POSSIBLE position
    out<-data.table(out, rbindlist(lapply(1:nrow(out), function(i){
      vardat<-oneclass[relPos==out[i, relpos] & substr(Region, nchar(Region), nchar(Region)) == out[i, substr(flank, 1, 1)]]
      return(data.table(n.tRNAs.var = nrow(vardat)))
    })))
    out[, n.tRNAs.invar:= nrow(onet) - (n.tRNAs.var)]
    out<-data.table(out, rbindlist(lapply(1:nrow(out), 
                                          function(x) out[x, dtpci(x = n.tRNAs.var, n = n.tRNAs.var + n.tRNAs.invar, suff = ".tRNAs.var")])))
    
    # Return
    return(out)
  }
  
  # --- Break data into sub dts based on mutation class to operate over
  fdat.per<-append(list(all = onef,  ## ANY mutations
                 any.SNV = onef[paste(Ref, Alt, sep = " > ") %in% mutsplit[, mutname], ], ## ANY SNV mutation (in class in mutsplit)
                 other = onef[!paste(Ref, Alt, sep = " > ") %in% mutsplit[, mutname], ]), ## NOT any of the input SNV mutations
  append(
    lapply(1:nrow(mutsplit), function(x){ # Ref vs Alt
    onef[paste(Ref, Alt, sep = " > ")==mutsplit[x , mutname]]
  }),
    lapply(1:nrow(mutsplit), function(x){ # Major vs minor
    onef[paste(majorAl, minorAl, sep = " > ")==mutsplit[x , mutname]]
  })))
  # names(fdat.per)<-c("all", "any.SNV", "other", paste("Ref-Alt:", mutsplit$mutname), paste("Major-Minor:", mutsplit$mutname))
  metadat.per<-data.table(mutclass = c("all", "any.SNV", "other", mutsplit$mutname, mutsplit$mutname),
                          alleles = c("-", "-", "-", rep("Ref > Alt", nrow(mutsplit)), rep("Major > Minor", nrow(mutsplit)))) # names broken down as I want
  
  # --- Calculate proportion of tRNAs with variant at each location
  ## All variants
  pvartrna.all<-rbindlist(lapply(1:nrow(metadat.per), function(i){
    # **Flip relpos in input whenever strand is - strand! (Input has relative position based on start/end as L/R on chromosome, but has flank sensibly in terms of 5'/3' based on transcription direction)
    # just do here - don't want to re-do when do pvartrna.low below [which is what happened when I did this inside fn]
    # ** INPUT has relative position in terms of to the L or R of the gene - not corrected for strandedness - but flank call is about strandedness
    fdat.per[[i]][Strand=="-", relPos:=relPos*-1]
    
    # do calcs
    pt<-proptvar(oneclass = fdat.per[[i]], onet, lflank, rflank)
    pt<-data.table(metadat.per[i, ], pt)
    return(pt)
  }))
  ## Proportion with LOW FREQ SNP (as in Thornlow 2018) (MAF < 0.05)
  pvartrna.low<-rbindlist(lapply(1:nrow(metadat.per), function(i){
    pt<-proptvar(oneclass = fdat.per[[i]][minorAlFreq<0.05], onet, lflank, rflank)
    pt<-data.table(metadat.per[i, ], pt)
    return(pt)
  }))
  ## Combine
  pvartrna<-rbind(data.table(varfreq = "all", pvartrna.all),
                  data.table(varfreq = "MAF < 0.05", pvartrna.low))
  
  # ....above is PER SITE; add option to do SLIDING WINDOW??? ideally, yes....blarghhhhh. NO: for proportion, fine to do window as AVERAGE
  #         across the bases - possible I can do that in ggplot; can definitely do on top as needed [median in the window centered there...]

  
  # --- Calculate Pi...somehow....MAYBE, NOT DEFINITELY
  # Other nucleotide diversity calcs?
  # # lapply(1:nrow(metadat.per); fdat.per[[x]])

  # --- Return
  return(list(pvartrna = data.table(annotdt, pvartrna)))
}

# --- Plotting functions etc
myggtheme<-theme_bw() +
  theme(axis.title.x = element_text(size = 12), axis.text.x = element_text(size = 11), 
        axis.title.y = element_text(size = 12), axis.text.y = element_text(size = 11),
        legend.title = element_text(size = 12), legend.text = element_text(size=11), 
        strip.text.x = element_text(size = 11), strip.text.y = element_text(size = 11), title = element_text(size = 14), 
        strip.text.x.top = element_text(size = 11), strip.text.x.bottom = element_text(size = 11), strip.text.y.right = element_text(size = 11), strip.text.y.left = element_text(size = 11),
        plot.subtitle = element_text(size = 13))

flanklocplot<-function(dat, xcol = "relpos", ycol = "p.tRNAs.var", ymincol = "low95ci.tRNAs.var", ymaxcol = "high95ci.tRNAs.var",
                       colcol = "mutclass", colvec = NULL, addcis = F, xlab = "Position in tRNA gene region",
                       ylab = "Polymorphism (prop. tRNAs with variation)", gbodyspan = 20, grect = T,
                       labaxisby = 5, drawinnerflank = T, inner5 = -20, inner3 = 10, labcol = ""){
  # Makes location vs. variation-type plot for tRNA flank regions.
  # In: dat, data.table with all columns of interest **including flank col and all the other named ones
  #     xcol, name of column with relative position from tRNA
  #     ycol, name of column with numeric value to plot
  #     ymincol, name of column with lower 95% CI bounds (used for plotting if addcis = T)
  #     ymaxcol, name of column with upper 95% CI bounds (used for plotting if addcis = T)
  #     colcol, column to color lines by 
  #     colvec, optional - values for scale color manual (lines) - names are colcol values
  #     addcis, if T, 95% CI lines will be added
  #     gbodyspan, how much distance to put between L and R flanks
  #     grect, put rectangle where tRNA would be (between flanks)
  #     labaxisby, how often to put an x axis tick
  #     drawinnerflank, add lines where inner flank bounds are?
  #     inner5, inner3: where to add inner flank boundaries if drawinnerflank = T
  #     labcol, legend label for color
  
  # Set up data. KEEP ALL ORIGNAL INFO TOO FOR, EG, FACETING
  pdat<-data.table(x = dat[, get(xcol)], y = dat[,get(ycol)], ymin = dat[,get(ymincol)], ymax = dat[, get(ymaxcol)],
                   coldat = dat[,get(colcol)], flank = dat[, flank],
                   dat[, .SD, .SDcols = which(!names(dat)%in%c(xcol, ycol, ymincol, ymaxcol, colcol, "flank"))])
  if(!is.null(colvec)){
    pdat[, coldat:=factor(coldat, levels = names(colvec))]
  }
  ## Add space for gene body
  pdat[, xlab:=x]
  pdat[flank=="5'", x:=x-(gbodyspan/2)]
  pdat[flank=="3'", x:=x+(gbodyspan/2)]
  
  # Set up X axis
  mybreaks<-sort(unique(pdat[, x]))
  myxlabs<-sort(c(unique(pdat[,xlab]), 0))
  mysub<-which(myxlabs%%labaxisby==0)
  mybreaks<-mybreaks[mysub]
  myxlabs<-myxlabs[mysub]
  
  # Y axis
  myylim<-c(0, pdat[,max(y) + 0.05*max(y)])
  if(addcis==T){ # expand for ymax
    myylim<-c(0, pdat[,max(ymax) + 0.05*max(ymax)])
  }
  
  # Basic plot
  plt<-ggplot(mapping = aes(x, y)) + 
    geom_line(data = pdat[flank=="5'", ], aes(color = coldat)) +
    geom_line(data = pdat[flank=="3'", ], aes(color = coldat)) +
    scale_x_continuous(expand = c(0, 0), breaks = mybreaks, labels = myxlabs) + # need 0 for L and R flanks
    scale_y_continuous(expand = c(0,0), limits = myylim) +
    labs(color = labcol) +
    xlab(xlab) + ylab(ylab) + myggtheme 
  
  ## Add color scale if desired
  if(!is.null(colvec)){
    plt<- plt + scale_color_manual(values = colvec)
  }
  ## Add lines for 95%CI if desired
  if(addcis==T){
    plt<-plt + 
      geom_line(data = pdat[flank=="5'", ], aes(x = x, y = ymax, color = coldat), linetype = "dashed", linewidth = 0.2) +
      geom_line(data = pdat[flank=="3'", ], aes(x = x, y = ymax,color = coldat), linetype = "dashed", linewidth = 0.2) +
      geom_line(data = pdat[flank=="5'", ], aes(x = x, y = ymin, color = coldat), linetype = "dashed", linewidth = 0.2) +
      geom_line(data = pdat[flank=="3'", ], aes(x = x, y = ymin,color = coldat), linetype = "dashed", linewidth = 0.2)
  }
  ## Add rectangle where tRNA gene would be if desired
  if(grect==T){
    plt<-plt + geom_rect(data = pdat, xmin = pdat[flank=="5'", max(x)], xmax = pdat[flank=="3'", min(x)],
                         ymin = -Inf, ymax = Inf, color = "gray50")
  }
  ## Add inner flank boundaries
  if(drawinnerflank==T){
    plt<-plt + geom_vline(xintercept = c(inner5 - (gbodyspan/2), inner3 + gbodyspan/2)) # in theory should work....
      # *** WORKS IN ONE CASE, NEED TO CHECK IF CHANGE THINGS...couldn't get to work later...
  }
    
  # figure out
  #     colors
  #     inner, outer flank boundaries - lines or rects
  #     ....if do any smoothing across adjacent bases

  return(plt)
}

#### Arguments and inputs ####
p<-arg_parser("Plot/analyze tRNA flank regions variation data (from trnaflankvcfvars.nf)", 
              name = "flankvariation.R", hide.opts = TRUE)
# Input files
p<-add_argument(p, "--speciesf",
                help = "File containing information on all species to process here. Columns infilename (exactly how all files have this species in their name),
                displayname (name that should be used for plot outputs etc), shortname (no-spaces name for ouptut files, sorting, etc - either shorter than or same as infilename, probably).
                In order you'd like plots to be in!",
                type = "character")
p<-add_argument(p, "--flankvars",
                help = "EXAMPLE Path to file with tRNA flank variant info output by get_strain_variants-flank.py.
                Where species ID/species specific info is, put SAMP instead",
                type = "character")
p<-add_argument(p, "--trnainfoal",
                help = "Path to file containing tRNA per allele information including original codon & allele as well as (potential) remolded one. - *alleleinfo_counts_worigcodonetc.txt output of show_mutational_variation.R",
                type = "character")
p<-add_argument(p, "--trnainfog",
                help = "EXAMPLE Path to file containing tRNA per gene information including name & strand. - *strain_trnas_info.txt output of build_alt_sequences.py
                Where species ID/species specific info is, put SAMP instead",
                type = "character")
p<-add_argument(p, "--lflank",
                help = "Start and end of basepairs away from gene used in input for 5' flank, comma-separated. E.g -40,0",
                default = "-40,0")
p<-add_argument(p, "--rflank",
                help = "Start and end of basepairs away from gene used in input for 3' flank, comma-separated. E.g 0,40",
                default = "0,40")


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

# flank boundaries
lflank<-as.numeric(strsplit(p$lflank, ",")[[1]])
rflank<-as.numeric(strsplit(p$rflank, ",")[[1]])

#### Read in and format data ####
cat("....Reading in data....\n")
sinfo<-fread(p$speciesf, header = T)

# --- tRNA gene information
tdat<-freadgetpertrna(trnainfoalf = p$trnainfoal, examptrnag = p$trnainfog, sinfo = sinfo)

# --- flank information
# all/initial
fdat<-fread.mult(sinfo = sinfo, exampfile = p$flankvars, myselect = NA, skip = 1)
if("##Chr"%in%names(fdat)){
  setnames(fdat, "##Chr", "Chr")
}
# Annotate with major allele, minor allele info
fdat[, `:=`(majorAl = ifelse(nHomRef>=nHomAlt, Ref, Alt),
            minorAl = ifelse(nHomRef<nHomAlt, Ref, Alt),
            n.majorAl = ifelse(nHomRef>=nHomAlt, nHomRef, nHomAlt),
            n.minorAl = ifelse(nHomRef<nHomAlt, nHomRef, nHomAlt))]
fdat[, minorAlFreq:=n.minorAl/(n.minorAl + n.majorAl)] # missing excluded from num & denom
setcolorder(fdat, c("homRef", "homAlt", "missingOrHet"), after = ncol(fdat))

#### Calculate flank mutation metrics ####
cat("....Calculating flank mutation metrics....\n")

# --- Get other relevant info
# Mutation categories
mutsplit<-data.table(from = rep(c("A", "C", "G", "T"), each = 4),
                     to = rep(c("A", "C", "G", "T")))
mutsplit<-mutsplit[from!=to, ]
mutsplit[, mutname:=paste(from, to, sep = " > ")] # can further annotate here - TAM categories etc

# Splits/data to compute this in [test info dt]
## species: done in lapply loop below
## strand: 
strandinfo<-data.table(descrip = c("Either strand", "+ strand", "- strand"),
                       evaltext = c("Strand%in%c('+', '-')", "Strand=='+'", "Strand=='-'"))
## amino acid: MEH, not doing for now - too subdivided for now. Would be an option to come back to, could just add in loop or somewhere
gsplits<-data.table(descrip = c("All tRNA genes", "Active tRNA genes", "Inactive tRNA genes"),
                    evaltext = c("all.lost%in%c(T, F)", "all.lost==F", "all.lost==T"))

# SPLIT DATA by [get in super long format] -- notes leftover/ not yet incorp
#   Above including/excluding ones with Undet amino acid? (though that'll be pseud, probably)
#   Split BY AMINO ACID -- backbone/expected and possibly actual anticodon??

# --- Perform calculations
fvar<-lapply(sinfo$displayname, function(species){
  lapply(1:nrow(strandinfo), function(strandi){
    lapply(1:nrow(gsplits), function(gi){
      thisannot<-data.table(displayname = species, strand = strandinfo[strandi, descrip], genes = gsplits[gi, descrip])
      # narrow tdat AND fdat based on all of these things [tdat first, then fdat of ones still in tdat]
      onet<-tdat[displayname==species & eval(parse(text = strandinfo[strandi, evaltext])) & 
                   eval(parse(text = gsplits[gi, evaltext]))]
      onef<-fdat[displayname==species & tRNA%in%onet$tRNA, ]
      # Run calculations
      out<-calcflankvar(onef, onet, annotdt = thisannot, mutsplit = mutsplit, 
                        lflank = lflank, rflank = rflank)
      return(out)
    })
  })
})

# --- get in long format; SAVE data
pvartrna<-rbindlist(lapply(fvar, function(x) rbindlist(lapply(x, function(y) rbindlist(lapply(y, function(z) z$pvartrna))))))
write.table(pvartrna, gzfile(file.path(p$outdir, paste0(p$baseoutname, "_perpositiontRNAvarianceprop.txt.gz"))),
            sep = "\t", quote = F, row.names = F)

# For function that computes info per base (for any provided set of genes)
#   ** Deal with MISSINGNESS - denominator for proportion variation type calcs should be OF STRAINS WITH CALLS AT THAT SITE
#   Do per possible nt-to-nt. May need to do that in WINDOWS rather than per bp. [possible do everything in multiple sliding windows?]
#       ***NOT NECESSARILY REF-TO-ALT. MAYBE major allele-to-minor allele, presuming minors are newer?? ***
#   Couple different metrics. Pi, potentially; MAF [though...how to do over multiple genes/sites...]
#         just raw proportion of tRNAs that have variants at this position**
#   confidence intervals?? might be tough
#   ****include, probably the tRNAs that DO NOT have variants [as denominator/not altered]?!? will need their strand probably....

#### Make plots! ####
cat("....Making plots....\n")
# --- Set up plot info data
# Colors
onecol.all<-c("black")
names(onecol.all)<-"all"
onecol.snv<-c("black")
names(onecol.snv)<-"any.SNV"
onecol.oth<-c("black")
names(onecol.oth)<-"other"

tams<-c("C > T", "G > A") 
# orig had reverse comped but I think that was a strand vs flank assigning issue originally **I think this is right given everything was previously reverse comped I THINK (canonical top strand goes C>T, reverse G>A)

blpal<-colorRampPalette(c("lightblue", "darkblue"))
mutclasscol<-blpal(nrow(mutsplit)-2)
mutclasscol<-c("red", "darkorange", mutclasscol)
names(mutclasscol)<-c(tams, mutsplit[!mutname%in%tams,mutname])
mutclasscol<-rev(mutclasscol)

tamcol<-rep("gray30", nrow(mutsplit))
tamcol[which(mutsplit$mutname%in%tams)]<-c("red", "darkorange")
names(tamcol)<-mutsplit$mutname
tamcol<-tamcol[c(tams, mutsplit[!mutname%in%tams,mutname])]  
tamcol<-rev(tamcol)

# mutations
mutinfo<-data.table(descrip = c("All", "Any SNV", "Other (non-SNVs)", "Ref > Alt", "Major > Minor"),
                    seltext = c("mutclass=='all'", "mutclass=='any.SNV'", "mutclass=='other'", "alleles=='Ref > Alt'", "alleles=='Major > Minor'"),
                    colvec1 = c("onecol.all", "onecol.snv", "onecol.oth", rep("mutclasscol", 2)),
                    colvec2 = c(rep(NA, 3), rep("tamcol", 2)),
                    addcis = c(T, T, T, F, F))

# --- Location vs. proportion variant tRNAs
invisible(lapply(1:nrow(sinfo), function(speci){ # one PDF per species
  pdf(file.path(p$outdir, paste0(p$baseoutname, "_locVpropvartRNAs_", sinfo[speci, shortname], ".pdf")),
      13, 8)
  
  lapply(1:nrow(mutinfo), function(mutind){
    lapply(strandinfo$descrip, function(strandval){ # one page per strand
      # Facet by mutation MAF, gene set
      print(
      flanklocplot(dat = pvartrna[displayname==sinfo[speci, displayname] & strand==strandval &
                                    eval(parse(text = mutinfo[mutind, seltext])),],
                   colvec = eval(as.name(mutinfo[mutind, colvec1])), addcis = mutinfo[mutind, addcis],
                   labcol = "Specific mutation") +
        facet_grid(genes~varfreq) +
        ggtitle(sinfo[speci, displayname],
                subtitle = paste("Mutation class:", mutinfo[mutind, descrip], "|", "Strand:", strandval))
      )
      
      # IF there's another color vector (TAM), make another plot with that colorway
      if(!is.na(mutinfo[mutind, colvec2])){
        print(
        flanklocplot(dat = pvartrna[displayname==sinfo[speci, displayname] & strand==strandval &
                                      eval(parse(text = mutinfo[mutind, seltext])),],
                     colvec = eval(as.name(mutinfo[mutind, colvec2])), addcis = mutinfo[mutind, addcis],
                     labcol = "Specific mutation\n(second colorway)") +
          facet_grid(genes~varfreq) +
          ggtitle(sinfo[speci, displayname],
                  subtitle = paste("Mutation class:", mutinfo[mutind, descrip], "|", "Strand:", strandval))
        )
      }
    })
  })
  
  invisible(dev.off())
}))

# --- Distribution of each mutation across the various regions...

# --- proportions of ALL variants vs. low freq variants AGAINST each other (for each position) (could include CIs)

#### plots to make/ analyses to do list ####



#### Script completion message & session information ####
cat("....flankvariation.R processing complete! Session information:....\n")
sessionInfo()