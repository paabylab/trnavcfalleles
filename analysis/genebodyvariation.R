#! /usr/bin/env/ Rscript
# Plot/analyze tRNA gene body variation data (from get_strain_variants_relpos.py)
# By Avery Davis Bell, begun 2025.03.31
require(argparser, quietly = T)
require(data.table, quietly = T)
require(ggplot2, quietly = T)
require(RColorBrewer, quietly = T)
require(R.utils)
require(ggpattern)

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

freadgetpertrnaflags<-function(trnainfoalf, examptrnag, exampssflag, sinfo, stripchrpref = T){
  # Reads in file containing tRNA per allele information including original codon & allele as well as (potential) remolded one. - *alleleinfo_counts_worigcodonetc.txt output of show_mutational_variation.R
  # Gets relevant information (see output description) in per-tRNA [format for each species]; also includes info on whether sec struct decoding threw any flags
  # In: info, data.table with columns infilename (exactly how all files have this species in their name), displayname, shortname
  #     trnainfoalf, path to file containing tRNA per allele information including original codon & allele as well as (potential) remolded one. - *alleleinfo_counts_worigcodonetc.txt output of show_mutational_variation.R
  #     examptrnag, EXAMPLE Path to file containing tRNA per gene information including name & strand. - *strain_trnas_info.txt output of build_alt_sequences.py. Where species ID/species specific info is, put SAMP instead
  #     exampssflag, EXAMPLE Path to file containing tRNA per gene information including name & strand. - *strain_trnas_info.txt output of build_alt_sequences.pyWhere species ID/species specific info is, put SAMP instead
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
  ##  < All columns freom exampssflag - PossiblePseudogene, Intron, NoFlags, TooShort, etc etc - from sec struct analysis>
  
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
  
  # Info on sec structure stuff: Restricts to REFERENCE if alleles are included (vs. not having strain names)
  ssdat<-fread.mult(sinfo, exampfile = exampssflag)
  ssdat.new<-data.table(ssdat[, c(1:3, 5, 9 , 10), with = F], ssdat[, lapply(.SD, as.logical), .SDcols = -c(1:3, 5, 9 , 10)])
  setcolorder(ssdat.new, names(ssdat))
  ssdat<-ssdat.new
  if(ssdat[, sum(grepl("Reference", tRNA))] > 0){ # keep only reference and strip the allele info
    ssdat<-ssdat[grepl("Reference", tRNA)]
    ssdat[, tRNA:=ssdat[, tstrsplit(tRNA, "_Reference")]$V1]
  }
  setkey(ssdat, displayname, tRNA)
  
  # Format &  Return
  odat<-ssdat[stdat[tdat]]
  if(stripchrpref==T){
    odat[startsWith(tRNA, "chr"), tRNA:=substr(tRNA, 4, nchar(tRNA))]
  }
  setcolorder(odat, c(names(stdat), names(tdat)[!names(tdat)%in%c("displayname", "tRNA")]))
  odat[,i.shortname:=NULL]
  
  setkey(odat, displayname, tRNA)
  return(odat)
}

# --- computations
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

calctgvar<-function(onev, onet, tstruct, annotdt = NULL, mutsplit){
  # Calculates proportion of tRNAs in input set that have variation at each base/struct of 'canonical' tRNA
  #     based on relative position in onev input
  # In: onev, variant information. Key columns = tRNA, tRNA_pos, **structure, substructure, nSubStructure, substructure_pos!!
  #     onet, information on ALL tRNAs (with or without variation) for this set of interest [freadgetpertrna output format]
  #     tstruct, info on all tRNA sec struct pieces as in onev. Columns structure, substructure, nSubStructure, canonicalnbp required
  #     annotdt, optional data.table to prepend to each row of output - i.e., describing this set of tRNAs
  #     mutsplit, data.table with specific mutations (from nt, to nt) to break out calculations by. Function also does 'other' - any not captured here - and 'all.SNV' (all THESE combined)  and 'any' (all ANY combined)
  #           Required columns: from (nt), to (nt), mutname (from > to format!!)
  # Out: data.table with one row per mutation class x allele freq x tRNA sec. substructure combination. Columns:
  #   <any in annotdt>
  #   mutclass - name of mutclass here
  #   alleles - "-" if doesn't make a difference; Ref > Alt or Major > Minor for specific mutation classes
  #   structure, substructure, nSubStructure, canonicalnbp - from tstruct input
  #     n.tRNAs.var, # tRNAs with one or more variants in this substructure
  #     n.variants, absolute # variants observed in this substructure
  #     n.tRNAs.invar, # tRNAs with no variants recorded in this substructure
  #     <p, low95ci, high95ci>.tRNAs.var.any - proportion & 95% CI (binomial) for number variant tRNAs over all. Not corrected by structure length.
  #     <p, low95ci, high95ci>.varsPerBp - proportion & 95% CI (binomial) for number variants in region normalized to length of that region from input (may be imperfect, but closer to per-bp number). Multiple variants from same tRNA count where relevant.
  
  
  # --- Subfunctions
  propvareachstruct<-function(oneclass, onet,  tstruct){
    # For each substructure in tRNA, figure out how many tRNAs have mutations *within that structure* [doing per base gets really complicated with different lengths]
    #   get proportion of tRNAs of input set that have a mutation here
    # In: oneclass, tRNA gene body variants to process - narrowed to set of interest (all summarized together here).
    #             key columns: structure/substructure/nSubStructure
    #     onet, all tRNAs to compare this set against (denominator, basically)
    #     tstruct, info on all tRNA sec struct pieces as in onev. Columns structure, substructure, nSubStructure, canonicalnbp required
    # Out: data.table with one row per sub-secondary-structure in tRNA (from tstruct). Columns:
    #     <4 relevant ones from tstruct>
    #     n.tRNAs.var, # tRNAs with one or more variants in this substructure
    #     n.variants, absolute # variants observed in this substructure
    #     n.tRNAs.invar, # tRNAs with no variants recorded in this substructure
    #     <p, low95ci, high95ci>.tRNAs.var.any - proportion & 95% CI (binomial) for number variant tRNAs over all. Not corrected by structure length.
    #     <p, low95ci, high95ci>.varsPerBp - proportion & 95% CI (binomial) for number variants in region normalized to length of that region from input (may be imperfect, but closer to per-bp number). Multiple variants from same tRNA count where relevant.
    
    # Set up each substructure [...or could pass in...]
    out<-data.table(tstruct[, .(structure, substructure, nSubStructure, canonicalnbp)])
    
    # Get number tRNAs with & without mutations in each substructure
    out<-data.table(out, rbindlist(lapply(1:nrow(out), function(i){
      vodat<-oneclass[structure==out[i, structure] & substructure==out[i, substructure]]
      return(data.table(n.tRNAs.var = length(vodat[, unique(tRNA)]), n.variants = nrow(vodat)))
    })))
    out[, n.tRNAs.invar:=nrow(onet) - n.tRNAs.var]
    
    # Get proportions (using n variants?, n tRNAs with variants); raw AND divided by number of 'canonical' bp there
    out<-data.table(out, rbindlist(lapply(1:nrow(out), function(x){
      o1<-data.table(out[x, dtpci(x = n.tRNAs.var, n = n.tRNAs.var + n.tRNAs.invar, suff = ".tRNAs.var.any")], # proportion tRNAs with variant base at ANY of bases in substructure
                     out[x, dtpci(x = round(n.variants/canonicalnbp), n = n.tRNAs.var + n.tRNAs.invar, suff = ".varsPerBp")]) # proportion variants normalized to canonical length of structure; multiple variants can count multiply (would on a per-bp level)
      return(o1)
      })))
    
    # Return
    return(out)
  }
  
  # --- Break data into sub dts based on mutation class to operate over
  vdat.per<-append(list(all = onev,  ## ANY mutations
                        any.SNV = onev[paste(ref, alt, sep = " > ") %in% mutsplit[, mutname], ], ## ANY SNV mutation (in class in mutsplit)
                        other = onev[!paste(ref, alt, sep = " > ") %in% mutsplit[, mutname], ]), ## NOT any of the input SNV mutations
                   append(
                     lapply(1:nrow(mutsplit), function(x){ # Ref vs Alt
                       onev[paste(ref, alt, sep = " > ")==mutsplit[x , mutname]]
                     }),
                     lapply(1:nrow(mutsplit), function(x){ # Major vs minor
                       onev[paste(majorAl, minorAl, sep = " > ")==mutsplit[x , mutname]]
                     })))

  metadat.per<-data.table(mutclass = c("all", "any.SNV", "other", mutsplit$mutname, mutsplit$mutname),
                          alleles = c("-", "-", "-", rep("Ref > Alt", nrow(mutsplit)), rep("Major > Minor", nrow(mutsplit)))) # names broken down as I want

  # --- Calculate proportion of tRNAs with variant at each location
  ## All variants
  pvartrna.all<-rbindlist(lapply(1:nrow(metadat.per), function(i){
    pt<-propvareachstruct(oneclass = vdat.per[[i]], onet, tstruct)
    pt<-data.table(metadat.per[i, ], pt)
    return(pt)
  }))
  ## Proportion with LOW FREQ SNP (as in Thornlow 2018) (MAF < 0.05)
  pvartrna.low<-rbindlist(lapply(1:nrow(metadat.per), function(i){
    pt<-propvareachstruct(oneclass = vdat.per[[i]][minorAlFreq<0.05], onet, tstruct)
    pt<-data.table(metadat.per[i, ], pt)
    return(pt)
  }))
  ## Combine
  pvartrna<-rbind(data.table(varfreq = "all", pvartrna.all),
                  data.table(varfreq = "MAF < 0.05", pvartrna.low))

  # --- Return
  return(data.table(annotdt, pvartrna))
}

# --- plotting
myggtheme<-theme_bw() +
  theme(axis.title.x = element_text(size = 12), axis.text.x = element_text(size = 11), 
        axis.title.y = element_text(size = 12), axis.text.y = element_text(size = 11),
        legend.title = element_text(size = 12), legend.text = element_text(size=11), 
        strip.text.x = element_text(size = 11), strip.text.y = element_text(size = 11), title = element_text(size = 14), 
        strip.text.x.top = element_text(size = 11), strip.text.x.bottom = element_text(size = 11), strip.text.y.right = element_text(size = 11), strip.text.y.left = element_text(size = 11),
        plot.subtitle = element_text(size = 13))

trnalocplot<-function(dat, ycol = "p.tRNAs.var.any", ymincol = "low95ci.tRNAs.var.any", ymaxcol = "high95ci.tRNAs.var.any",
                      tstruct,
                      colcol = "mutclass", colvec = NULL, addcis = F, xlab = "Position in tRNA gene",
                      ylab = "Polymorphism (prop. tRNAs with variation)", drawsubstructs = T, labcol = "",
                      facs = NULL){
  # Plots location in gene body vs proportion variation - in substructure CHUNKS (lines) rather than per bp (given input data)
  # In: dat, data. Must have columns structure, substructure, nSubStructure, canonicalnbp - for figuring out x axis where to plot; ycol, ymincol, ymaxcol values
  #     ycol, name of column with numeric value to plot
  #     ymincol, name of column with lower 95% CI bounds (used for plotting if addcis = T)
  #     ymaxcol, name of column with upper 95% CI bounds (used for plotting if addcis = T)
  #     tstruct, structure info about plotting. Must have structure, substructure as in dat, also structure_plot and substructure_plot, structure_plot_level & substructure_plot_level (used here for naming and ORDERING for legend etc)
  #     colcol, column to color lines by 
  #     colvec, optional - values for scale color manual (lines) - names are colcol values
  #     addcis, if T, 95% CI lines will be added
  #     drawsubstructs, if T, substructures are highlighted/called out/something....
  #     labcol, legend label for color
  #     facs: if going to facet LATER (does NOT facet here), tell the name of those facets - x axis computation has to happen internally. NULL or char vec
  #           ALSO need to do this for columns that are split into diff colors or whatever
  
  # --- Get data for plotting
  # General
  pdat<-copy(dat)
  pdat[, y :=get(ycol)]
  pdat[, ymin :=get(ymincol)]
  pdat[, ymax:=get(ymaxcol)]
  pdat[, coldat := get(colcol)]
  if(!is.null(colvec)){
    pdat[, coldat:=factor(coldat, levels = names(colvec))]
  }
  # Structure stuff for plotting
  pdat<-merge(pdat, tstruct)
  # pdat[is.na(substructure_plot), substructure_plot:="other"] better to leave NA so it doesn't get assigned a pattern
  pdat[, structure_plot:=factor(structure_plot, levels = unique(tstruct[order(structure_plot_level), structure_plot]))]
  pdat[, substructure_plot:=factor(substructure_plot, levels = c(na.omit(unique(tstruct[order(substructure_plot_level), substructure_plot]))))]
  # X axis locations (start and end of line segment) ** WITHIN ONE FACET
  if(is.null(facs)){
    pdat<-pdat[order(nSubStructure)]
    pdat[,xend:=cumsum(canonicalnbp)]
    pdat[,xstart:= c(1, xend[1:(nrow(pdat) -1)] + 1)]
  }else{
    setkeyv(pdat, facs)
    pdat<-pdat[order(nSubStructure)]
    pdat[,xend:=cumsum(canonicalnbp), by = facs]
    pdat[,xstart:= c(1, xend[1:(.N -1)] + 1), by = facs]
  }
  pdat[, xseg:=xstart - 0.5]
  fakend<-pdat[nSubStructure==max(nSubStructure)]
  fakend[, xseg:=xend + 0.5]
  pdat<-rbind(pdat, fakend) # add an 'end point' for final [1bp] structure
  
  # --- Basic plot
  # Y axis
  myylim<-c(0, pdat[,max(y) + 0.05*max(y)])
  if(addcis==T){ # expand for ymax
    myylim<-c(0, pdat[,max(ymax) + 0.05*max(ymax)])
  }
  plt<-ggplot(pdat) + geom_step(aes(x = xseg, y = y, color = coldat)) +
    scale_x_continuous(expand = c(0, 0)) +
    scale_y_continuous(expand = c(0,0), limits = myylim) +
    labs(color = labcol) + 
    xlab(xlab) + ylab(ylab) + myggtheme 
  
  # --- Add ons
  ## Substruct marking if desired. it is SLOW esp w/ facets
  if(drawsubstructs==T){
    plt <- plt + geom_rect_pattern(data = pdat[coldat==coldat[1]], aes(xmin = xstart - 0.5, xmax = xend + 0.5, ymin = -Inf, ymax = Inf, 
                               fill = structure_plot, pattern_shape = substructure_plot, pattern_color = structure_plot), alpha = 0.4, color = "gray",
                               pattern = 'pch', pattern_density = 0.2, pattern_spacing = 0.03) + # pattern_frequency, pattern_scale doesn't seem to change anything
      geom_step(data = pdat, aes(x = xseg, y = y, color = coldat)) + 
      labs(color = labcol, fill = "tRNA structure component", pattern_shape = "Secondary structure sub-type") + # shape or pattern not changing that label
      theme(panel.grid = element_blank()) +
      guides(pattern_color = "none", fill = guide_legend(override.aes = list(pch = 0), order = 2), color = guide_legend(order = 1))
  }
  
  ## Add color scale if desired
  if(!is.null(colvec)){
    plt<- plt + scale_color_manual(values = colvec)
  }
  ## Add lines for 95%CI if desired [ FIX!!]
  if(addcis==T){
    plt<-plt + 
      geom_step(aes(x = xseg, y = ymin, color = coldat), linewidth = 0.2) +
      geom_step(aes(x = xseg, y = ymax, color = coldat), linewidth = 0.2) 
  }

  return(plt)
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
p<-add_argument(p, "--trnainfoal",
                help = "Path to file containing tRNA per allele information including original codon & allele as well as (potential) remolded one. - *alleleinfo_counts_worigcodonetc.txt output of show_mutational_variation.R",
                type = "character")
p<-add_argument(p, "--trnainfog",
                help = "EXAMPLE Path to file containing tRNA per gene information including name & strand. - *strain_trnas_info.txt output of build_alt_sequences.py
                Where species ID/species specific info is, put SAMP instead",
                type = "character")
p<-add_argument(p, "--ssflags",
                help = "Path to file containing per-allele information categorizing how well sec structure assignment went ('info' file output of secstruct2seqpieces.py)
                Where species ID/species specific info is, put SAMP instead",
                type = "character")
p<-add_argument(p, "--trnasecstruct",
                help = "Path to file with information on how tRNA sec structures in --genevars input are arranged, etc. Columns:
                structure (as in --genevars), substructure (as in --genevars), nSubStructure (as in --genevars), canonicalnbp - length of this structure in a 72/73 bp tRNA [for plotting],
                structure_plot and substructure_plot: categories & names prettified/simplified for how you'd like to plot them;
                structure_plot_level & substructure_plot_level: RANKINGS of unique structure_plot and substructure_plot for ordering in plot",
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

# --- tRNA gene information * including sec struct flags in ref alleles
tdat<-freadgetpertrnaflags(trnainfoalf = p$trnainfoal, examptrnag = p$trnainfog, exampssflag = p$ssflags, sinfo = sinfo,
                           stripchrpref = F) # here, 'chr' thing is kept in both elegans inputs...for some reason...

# --- variant information
# all/initial
vdat<-fread.mult(sinfo = sinfo, exampfile = p$genevars, myselect = c('chrom', 'pos', 'ref', 'alt', 'nHomRef', 'nHomAlt', 'nNotMissingHet', 'nMissingHet', 'tRNA', 'tRNA_strand', 'tRNA_pos', 'structure', 'substructure', 'nSubStructure', 'substructure_pos'), 
                 skip = 1) # drop the loooong lists of who is what

# Annotate with major allele, minor allele info
vdat[, `:=`(majorAl = ifelse(nHomRef>=nHomAlt, ref, alt),
            minorAl = ifelse(nHomRef<nHomAlt, ref, alt),
            n.majorAl = ifelse(nHomRef>=nHomAlt, nHomRef, nHomAlt),
            n.minorAl = ifelse(nHomRef<nHomAlt, nHomRef, nHomAlt))]
vdat[, minorAlFreq:=n.minorAl/(n.minorAl + n.majorAl)] # missing excluded from num & denom

#### Calculate mutation metrics ####
cat("....Calculating gene body mutation metrics....\n")

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
## active vs not:
gsplits<-data.table(descrip = c("All tRNA genes", "Active tRNA genes", "Inactive tRNA genes"),
                    evaltext = c("all.lost%in%c(T, F)", "all.lost==F", "all.lost==T"))
## Sec structure flag vs not:
ssinf<-data.table(descrip = c("Any sec. struct. calling", "No flags on sec. struct. calling"),
                  evaltext = c("NoFlags%in%c(T, F)", "NoFlags==T"))

# --- Perform calculations
varinfo<-lapply(sinfo$displayname, function(species){
  lapply(1:nrow(strandinfo), function(strandi){
    lapply(1:nrow(gsplits), function(gi){
      lapply(1:nrow(ssinf), function(ssi){
        thisannot<-data.table(displayname = species, strand = strandinfo[strandi, descrip], genes = gsplits[gi, descrip], secstructflags = ssinf[ssi, descrip])
        # narrow tdat AND vdat based on all of these things [tdat first, then vdat of ones still in tdat]
        onet<-tdat[displayname==species & eval(parse(text = strandinfo[strandi, evaltext])) & 
                     eval(parse(text = gsplits[gi, evaltext])) & eval(parse(text = ssinf[ssi, evaltext]))]
        onev<-vdat[displayname==species & tRNA%in%onet$tRNA, ]
        # Run calculations
        out<-calctgvar(onev, onet, tstruct = tstruct, annotdt = thisannot, mutsplit = mutsplit)
        return(out)
      })
    })
  })
})

# << before getting here, exclude tRNAs where I couldn't do the sec struct stuff [from denom]? or do in here? - does automatically based on ssinf>>>

# --- get in long format; SAVE data
pvartrna<-rbindlist(lapply(varinfo, function(x) rbindlist(lapply(x, function(y) rbindlist(lapply(y, function(z) rbindlist(lapply(z, function(a) return(a)))))))))
write.table(pvartrna, gzfile(file.path(p$outdir, paste0(p$baseoutname, "_persecstructtRNAvarianceprop.txt.gz"))),
            sep = "\t", quote = F, row.names = F)

#### Do some plots, baby! ####
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

# Norm/un norm
pinfo<-data.table(suffInDat = c("tRNAs.var.any", c("varsPerBp")),
                  ldescrip = c("Total proportion tRNAs with any variants in region,\nunnormalized for length",
                               "Proportion tRNAs with variants, normalized for length;\nsame tRNA can have multiple mutations count"), # for plot subtitles or whatever
                  sdescrip = c("anyVarsUncorr", "varsPerBP"), # for file name etc
                  yax = c("Polymorphism\n(prop. tRNAs with variation)", "Polymorphism\n(prop. tRNAs with variation per bp in structure)"))

# --- Actually do the plots
lapply(1:nrow(sinfo), function(spind){ # species
  lapply(1:nrow(pinfo), function(pind){ #   plot prop of tRNA norm & unnorm 
    pdf(file.path(p$outdir, paste0(p$baseoutname, "_locVpropvartRNAs_", sinfo[spind, shortname], "_", pinfo[pind, sdescrip] , ".pdf")),
        13, 8)
    lapply(1:nrow(mutinfo), function(mutind){ # mutations, one page per
      lapply(strandinfo$descrip, function(strandval){ # one page per strand
        lapply(c("all", "MAF < 0.05"), function(maf){ # one page per all/MAF < 0.05
          cat(paste(c("...Plotting: ", sinfo[spind, displayname], pinfo[pind, ldescrip], mutinfo[mutind, descrip], strandval, maf, "\n"), collapse = "; ")) ## debug logging, doesn't report mutation
          # Original plot for this mutation type, faceted by inactive/active genes + whether sec struct had any flags or not
          ## Sec structs not delineated - too slow to do all
          # print(
          #   trnalocplot(dat = pvartrna[displayname==sinfo[spind, displayname] & strand==strandval & varfreq==maf &
          #                                eval(parse(text = mutinfo[mutind, seltext])),],
          #               ycol = paste("p", pinfo[pind, suffInDat], sep = "."), 
          #               ymincol = paste("low95ci", pinfo[pind, suffInDat], sep = "."),
          #               ymaxcol = paste("high95ci", pinfo[pind, suffInDat], sep = "."),
          #               tstruct = tstruct, colcol = "mutclass", colvec = eval(as.name(mutinfo[mutind, colvec1])),
          #               addcis = mutinfo[mutind, addcis], xlab = "Position in tRNA gene (bp)",
          #               ylab = pinfo[pind, yax], drawsubstructs = F, labcol = "Specific mutation") +
          #     facet_grid(genes~secstructflags) +
          #     ggtitle(paste(sinfo[spind, displayname], " | ", pinfo[pind, ldescrip]),
          #             subtitle = paste("Mutation class:", mutinfo[mutind, descrip], " | ",
          #                              "Minor allele freq:", maf, " | ",
          #                              "Strand:", strandval))
          # )
          ## Sec structs delineated
          print(
            trnalocplot(dat = pvartrna[displayname==sinfo[spind, displayname] & strand==strandval & varfreq==maf &
                                         eval(parse(text = mutinfo[mutind, seltext])),],
                        ycol = paste("p", pinfo[pind, suffInDat], sep = "."), 
                        ymincol = paste("low95ci", pinfo[pind, suffInDat], sep = "."),
                        ymaxcol = paste("high95ci", pinfo[pind, suffInDat], sep = "."),
                        tstruct = tstruct, colcol = "mutclass", colvec = eval(as.name(mutinfo[mutind, colvec1])),
                        addcis = mutinfo[mutind, addcis], xlab = "Position in tRNA gene (bp)",
                        ylab = pinfo[pind, yax], drawsubstructs = T, labcol = "Specific mutation",
                        facs = c("genes", "secstructflags", "mutclass")) +
              facet_grid(genes~secstructflags, scales = "free_y") +
              ggtitle(paste(sinfo[spind, displayname], " | ", pinfo[pind, ldescrip]),
                      subtitle = paste("Mutation class:", mutinfo[mutind, descrip], " | ",
                                       "Minor allele freq:", maf, " | ",
                                       "Strand:", strandval))
          )
          
          # IF there's another color vector (TAM), make another plot with that colorway
          if(!is.na(mutinfo[mutind, colvec2])){
            ## Sec structs not delineated - too slow
            # print(
            #   trnalocplot(dat = pvartrna[displayname==sinfo[spind, displayname] & strand==strandval & varfreq==maf &
            #                                eval(parse(text = mutinfo[mutind, seltext])),],
            #               ycol = paste("p", pinfo[pind, suffInDat], sep = "."), 
            #               ymincol = paste("low95ci", pinfo[pind, suffInDat], sep = "."),
            #               ymaxcol = paste("high95ci", pinfo[pind, suffInDat], sep = "."),
            #               tstruct = tstruct, colcol = "mutclass", colvec = eval(as.name(mutinfo[mutind, colvec2])),
            #               addcis = mutinfo[mutind, addcis], xlab = "Position in tRNA gene (bp)",
            #               ylab = pinfo[pind, yax], drawsubstructs = F, labcol = "Specific mutation\n(second colorway)") +
            #     facet_grid(genes~secstructflags) +
            #     ggtitle(sinfo[spind, displayname],
            #             subtitle = paste("Mutation class:", mutinfo[mutind, descrip], " | ",
            #                              "Minor allele freq:", maf, " | ",
            #                              "Strand:", strandval))
            # )
            ## Sec structs delineated
            print(
              trnalocplot(dat = pvartrna[displayname==sinfo[spind, displayname] & strand==strandval & varfreq==maf &
                                           eval(parse(text = mutinfo[mutind, seltext])),],
                          ycol = paste("p", pinfo[pind, suffInDat], sep = "."), 
                          ymincol = paste("low95ci", pinfo[pind, suffInDat], sep = "."),
                          ymaxcol = paste("high95ci", pinfo[pind, suffInDat], sep = "."),
                          tstruct = tstruct, colcol = "mutclass", colvec = eval(as.name(mutinfo[mutind, colvec2])),
                          addcis = mutinfo[mutind, addcis], xlab = "Position in tRNA gene (bp)",
                          ylab = pinfo[pind, yax], drawsubstructs = T, labcol = "Specific mutation\n(second colorway)",
                          facs = c("genes", "secstructflags", "mutclass")) +
                facet_grid(genes~secstructflags, scales = "free_y") +
                ggtitle(paste(sinfo[spind, displayname], " | ", pinfo[pind, ldescrip]),
                        subtitle = paste("Mutation class:", mutinfo[mutind, descrip], " | ",
                                         "Minor allele freq:", maf, " | ",
                                         "Strand:", strandval))
            )
          } # end if another colorway
          
        }) # end lapply MAF
      }) # end lapply strand
    }) # end lappply mutations
    invisible(dev.off())
  }) # end plot prop of tRNA norm & unnorm 
}) # end lapply species

# plot w/ & w/o substructure fill/pattern labels
#   species: diff PDFs
#   plot prop of tRNA norm & unnorm [var.any, varsPerBp] - must be different plots: diff PDFs
#   MAFs - facets - ?
#   mut types - pages
# strands - strandinfo- dif pages
# active/inactive - gsplits - facets?
# with & without flags - ssinf - facets?

#### Script completion message & session information ####
cat("....genebodyvariation.R processing complete! Session information:....\n")
sessionInfo()

