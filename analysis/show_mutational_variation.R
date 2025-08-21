#! /usr/bin/env/ Rscript
# Visualize/analyze tRNA allele information (using outputs of getstrainxtrnacalls.R, possibly others)                                      
# by Avery Davis Bell, begun 2024.12.03

if(length(.libPaths())==1){
  mylibloc <- .libPaths()[1]
}else{ # Presumes on PACE. This is not the best...
  mylibloc <- .libPaths()[grep(R.Version()$platform, .libPaths())]
}
require(data.table, lib.loc = mylibloc) 
require(argparser)
require(ggplot2)
require(ggh4x)
require(scatterpie)
require(formattable)

#### plot theme etc ####
myggtheme<-theme_bw() +
  theme(axis.title.x = element_text(size = 12), axis.text.x = element_text(size = 11), 
        axis.title.y = element_text(size = 12), axis.text.y = element_text(size = 11),
        legend.title = element_text(size = 12), legend.text = element_text(size=11), 
        strip.text.x = element_text(size = 11), strip.text.y = element_text(size = 11), title = element_text(size = 14), 
        strip.text.x.top = element_text(size = 11), strip.text.x.bottom = element_text(size = 11), strip.text.y.right = element_text(size = 11), strip.text.y.left = element_text(size = 11),
        plot.subtitle = element_text(size = 13))

sci2carrot<-function(n){
  # converts n in scientific notation to a bquote-able version [text though]
  # e.g. 1.1e-14 goes to '1.1*10^-14'
  # vectorized - works on vector of ns
  
  # Get in e notation as character
  ne<-paste(formattable::scientific(n, digits = 1))
  
  # convert to funky carrot notation (character)
  ne.split<-strsplit(ne, "e")
  n.out<-sapply(ne.split, function(x){
    if(substr(x[2], 2, 2)=="0"){
      x[2]<-paste0("-", substr(x[2], 3, nchar(x[2])))
    }
    n.out<-paste0(x[1], "^", x[2])
    return(n.out)
  })
  
  # Return
  return(n.out)
}

#### Functions ####
fread.mult<-function(sinfo, exampfile){
  # Reads in a file for each species in sinfo
  # In: sinfo, data.table with columns infilename (exactly how all files have this species in their name), displayname, shortname
  #     exampfile, path to file to process but wherever infilename is, should say SAMP instead. Should have header
  # Out: data.table containing whatever's in exampfile for each species in sinfo. Columns added to the front are displayname and shortname
  
  dat<-rbindlist(
    lapply(1:nrow(sinfo), function(s){
      dat1<-fread(gsub("SAMP", sinfo[s, infilename], exampfile), header = T)
      dat1<-data.table(dat1, sinfo[s, .(displayname, shortname)])
      setcolorder(dat1, c("displayname", "shortname"))
    })
  )
  return(dat)
}

refpseud<-function(peral, sinfo){
  # Finds any genes with alleles that are pseudo in ref but has non-pseudo genes in other
  # In: peral, data.table with species and per-allele info (*alleleinfo_counts_wmissing.txt output of getstrainxtrnacalls.R)
  #     sinfo, species info for species in per allele. Columns infilename, displayname (in peral), shortname (in peral)
  # Out: List of data.tables:
  #       $pseudorefgs, one row per gene that is pseudo in ref, contains total number alleles and number that might be functional and strain numbers for different allele classifications
  #       $summ, per-species summary: n_pseudo_in_ref (# genes pseudo in ref), n_pseudo_in_ref_mult_alleles (# of these that are variable in population), n_pseudo_in_ref_has_fn_allele (# of these that have putatively functional allele)
  #       $pseudoandfunc, input (peral) info for each allele from genes that are pseudo in ref but also have putatively functional alleles
  
  # Narrow to alleles in genes that have pseudo ref allele
  psref<-peral[grepl("pseudo", Note) & grepl("Reference", allelename), ] # All pseud-in-ref alleles
  all.psref<-rbindlist(lapply(1:nrow(sinfo), function(s){ # get all the alleles for genes that are pseud-in-ref
    gs<-psref[shortname==sinfo[s, shortname], tRNA]
    out<-peral[shortname==sinfo[s, shortname] & tRNA%in%gs]
    return(out)
  }))
  
  # Get per gene?
  setkey(all.psref, displayname, tRNA)
  psrefg<-all.psref[, .(n_alleles = .N, n_pseudo_alleles = sum(grepl("pseudo", Note) | Classification=="Lost"), 
                        n_func_alleles = sum(Classification%in%c("Best", "Functional", "Altered")),
                        n_strains_pseudo = sum(n.allele[Classification=="Lost"]), 
                        n_strains_func = sum(n.allele[Classification%in%c("Best", "Functional", "Altered")])), 
                    by = .(displayname, tRNA)]
  
  # All alleles for genes that have putatively fn'al alleles despite pseudo in ref
  pseudandfunc<-rbindlist(lapply(1:nrow(sinfo), function(s){
    peral[shortname==sinfo[s, shortname] & tRNA%in%psrefg[displayname==sinfo[s, displayname] & n_func_alleles>0, tRNA], ]
  }))
  
  # High level counts
  setkey(psref, displayname)
  summ<- psref[, .N, by = displayname][psref[, sum(VariableInPop==T), by = displayname]][psrefg[, sum(n_func_alleles > 0), by = displayname]]
  setnames(summ, c("displayname", "n_pseudo_in_ref", "n_pseudo_in_ref_mult_alleles", "n_pseudo_in_ref_has_fn_allele"))
  
  # Return all levels
  return(list(pseudorefgs = psrefg, summ = summ, allpseudoandfunc = pseudandfunc))
}

codoncm<-function(peral){
  # Gets the codon each allele should have been - i.e., if allele's AA != AlleleCM, grab the codon from an allele of this gene where it is
  # Enables some plotting stuff
  # In: peral, data.table with per-allele info for one species (*alleleinfo_counts_wmissing.txt output of getstrainxtrnacalls.R) and species info.
  #         Will be keyed by displayname, tRNA; also needs columns AA, Codon, AlleleCM. All others kept
  # Out: input data.table (possibly reordered) with new columns added:
  #     Codon.orig, what the codon was at the original amino acid this tRNA codes for - usually same as Codon, but in case of isotype switches, this helps group it with AlleleCM
  #     all.altered, usually F; or T - all alleles were flagged as isotype switches here; useful to know about this. CAN be this AND have classification be 'Lost'
  
  onecodoncm<-function(onetrna){
    # Gets 'original' codon for all alleles for one tRNA. input needs allelename - gets everything back together afterwards
    # Flags cases where all alleles were 'Altered' (no match: 'original' codon is assigned to be the input/new codon - so this will be grouped with its others with same codon)
    outc<-rep(NA, nrow(onetrna))
    altflag<-F
    
    if(all(onetrna[, AA==AlleleCM])){ # no mismatch
      outc<-onetrna$Codon
    }else if(onetrna[, sum(AA==AlleleCM)] > 0){ # one but not all mismatch
      origc<-onetrna[AA==AlleleCM, Codon][1]
      outc[onetrna[, AA==AlleleCM]]<-onetrna[AA==AlleleCM, Codon]
      outc[onetrna[, AA!=AlleleCM]]<-origc
    }else if(all(onetrna[, AA!=AlleleCM])){ # all mismatch: keep 'old' mismatch codon
       altflag<-T
       outc<-onetrna$Codon
    }
   
    return(data.table(allelename = onetrna$allelename,
                      Codon.orig = outc, 
                      all.altered = altflag))
  }
  
  setkey(peral, displayname, tRNA, allelename)
  newdat<-peral[, onecodoncm(.SD), by = .(displayname, tRNA)]
  setkey(newdat, displayname, tRNA, allelename)
  out<-peral[newdat]
  return(out)
}

allelebar<-function(pdat, classifcols){
  # Plots one bar per gene showing allele frequency (bar sums to 100%); alleles are sections of bar and they're colored by allele classification
  # In: pdat, data for one plot. Columns include tRNA, Codon, freq.ofallinclmissing, Classification, ....
  #     classifcols, vector of length and order of possible classifications, with colors they should take
  # Out: plot
  
  # Data setup
  pdat<-copy(pdat)
  pdat[, Classification:=factor(Classification, levels = names(classifcols))]
  
  # Make plot
  bp<-ggplot(pdat, aes(tRNA, freq.ofallinclmissing)) + 
    geom_bar(aes(fill = Classification), stat = "identity", position = "stack", color = "black") +
    scale_fill_manual(values = classifcols) + 
    scale_y_continuous(limits = c(0,1), expand = c(0, 0)) + 
    xlab("tRNA genes") + ylab("Allele frequency in population") + 
    facet_nested(. ~ AA + Codon, scales = "free", space = "free") + # from ggh4x package! hurray!
    myggtheme + theme(panel.grid = element_blank(), legend.position = "bottom", axis.text.x = element_blank(), axis.ticks.x = element_blank())
  
  # need to organize by codon, AA, etc....and label that way
  # Honestly faceting by AA probably makes sense....but need widths to scale so bar widths stay the same:  + facet_grid(.~AA, scales = "free", space = "free")
  
  # Within bars, order by:
  # --Classification [done with factor]
  # --within classification, highest freq on bottom [...hmmm, need to do some other sort]
  return(bp)
}

allelescatterpie<-function(peral.one, classifcols, missingcol = "white", aa.order, lowpad = 0, percodonwd = 2,
                           myfacet = "both", labelnalleles = T){
  # Makes 'scatter pie' with each tRNA a pie chart all piled up and small based on amino acid *in AlleleCM* (should do AlleleCM but also doing codon messes with that...)
  #           data needs to have codon it 'should' be!! .... based on other alleles at gene .....
  # In: peral.one, data.table with per-allele info for one species (*alleleinfo_counts_wmissing.txt output of getstrainxtrnacalls.R)
  #               Columns: tRNA, Codon, Codon.orig ***one that matches rest of alleles at gene where possible - from codoncm function****
  #               , AlleleCM [or AA, still deciding], Classification, freq.ofallinclmissing, all.altered [from codoncm function, can flag where genes might should be with another aa]
  #     classifcols, vector of length and order of possible classifications, with colors they should take. 'Missing' is added here to end.
  #     missingcol, color of pie slice containing missing data [data is re-formatted to explicitly include that in this function]
  #     aa.order, order in which to plot amino acids
  #     lowpad, pad allele frequencies below this to equal this. **IMPLEMENT, HOW DO w/r/t total circle??**
  #     percodonwd, number of 'stacks' of pies per codon (or aa if myfacet == aaonly)
  #     myfacet, "both" or "aaonly": do AA and codon nested facets ("both") or AA only facets ("aaonly")
  #     labelnalleles, T or F: add label of numer of alleles for each gene with more than one allele?
  # Out: plot. If faceted only by AA, genes (pies) are sorted first by codon, then by n alleles. If faceted out by codon within AA, sorted by n alleles.
  
  # Subfunctions
  getbestaa<-function(oneg){
    # Gets the 'best' amino acid to classify each tRNA gene as: the AlleleCM of the majority of alleles/the ones that don't have isotype switch (if any)
    if(oneg[, all(all.altered)]){
      # Hard to say here, just pick one!
      myaa<-oneg[1, AlleleCM]
    }else{
      myaa<-oneg[AA==AlleleCM & freq.ofallinclmissing==max(freq.ofallinclmissing), AlleleCM][1]
    }
    return(myaa)
  }
  
  # ---- format data: keep tRNAs together first
  ## Assign each gene to its appropriate group(s)
  alldat<-copy(peral.one)
  setkey(alldat, tRNA)
  alldat<-alldat[, getbestaa(.SD), by = tRNA][alldat]
  setnames(alldat, "V1", "groupaa")
  
  ## Get data all together with x, y values etc -- based on the new group(s)
  pdat<-rbindlist(lapply(aa.order, function(aa){
    codons.one<-alldat[groupaa==aa, sort(unique(Codon.orig))]
    if("NNN"%in%codons.one){ # Happens if all alleles are mismatches and there's an NNN in here. Not clear if get-around-able...
      codons.one<-c(codons.one[codons.one!="NNN"], "NNN")
    }
    if(myfacet=="aaonly"){ # Do coordinates only within AA!
      dat.alls<-alldat[groupaa==aa, ]
      setkey(dat.alls, tRNA)
      dat<-dat.alls[, .N, by = tRNA][order(N, decreasing = T)][,.(tRNA)] ## Sort tRNA genes by N alleles (more -> first)
      
      # Add in x, y coords for plotting (based on order data's currently in)
      if(nrow(dat)>1){
        dat[, x:=rep(c(1:percodonwd), ceiling(nrow(dat)/percodonwd))[1:nrow(dat)]]
        maxy<-dat[, max(table(x))]
        dat[, y:=rep(1:maxy, each = percodonwd)[1:nrow(dat)]]
      }else if(nrow(dat)==1){ # only one gene: put in the middle at the bottom
        dat[, `:=`(x = median(1:percodonwd), y = 1)]
      }
      # Add alleles back in!!
      setkey(dat, tRNA)
      setkey(dat.alls, tRNA)
      odat<-dat.alls[dat]
      # Add number alleles count for easy access (on per-tRNA basis)
      odat[, NAlleles:=.N, by = tRNA]
      ## Sort alleles within genes so ones with same classification come next to each other (?)
      odat<-odat[order(Classification), .SD, by = tRNA] # make sure this works...
      # Return
      out<-odat
    }else if(myfacet=="both"){ # Do coordinates within AA and Codon.orig
      out<-rbindlist(lapply(codons.one, function(cod){
        dat.alls<-alldat[groupaa==aa & Codon.orig == cod]
        setkey(dat.alls, tRNA)
        dat<-dat.alls[, .N, by = tRNA][order(N, decreasing = T)][,.(tRNA)] ## Sort tRNA genes by N alleles (more -> first)
        
        # Add in x, y coords for plotting (based on order data's currently in)
        if(nrow(dat)>1){
          dat[, x:=rep(c(1:percodonwd), ceiling(nrow(dat)/percodonwd))[1:nrow(dat)]]
          maxy<-dat[, max(table(x))]
          dat[, y:=rep(1:maxy, each = percodonwd)[1:nrow(dat)]]
        }else if(nrow(dat)==1){ # only one gene: put in the middle at the bottom
          dat[, `:=`(x = median(1:percodonwd), y = 1)]
        }
        
        # Add alleles back in!!
        setkey(dat, tRNA)
        setkey(dat.alls, tRNA)
        odat<-dat.alls[dat]
        # Add number alleles count for easy access (on per-tRNA basis)
        odat[, NAlleles:=.N, by = tRNA]
        ## Sort alleles within genes so ones with same classification come next to each other (?)
        odat[order(Classification), .SD, by = tRNA] # make sure this works...
        # Return
        return(odat)
      }))
      return(out)
    }
    out[, Codon.orig:=factor(Codon.orig, levels = codons.one)] # order nicely
    return(out)
  }))
  
  ## For each tRNA with missingness, add an 'allele' with all tRNA info with Missing classification
  pdat[, value:=freq.ofallinclmissing] # **scatterpie needs value column; good to have this before missingness for org purposes
  formiss<-unique(pdat[, .(tRNA, groupaa, displayname, shortname, VariableInPop, n.called, n.missing, 
                           Codon.orig = Codon.orig[1], all.altered, x, y, NAlleles)])[n.missing>0] # these are the ones missing 'allele' needs to be created for
  formiss[, `:=`(allelename = paste(tRNA, "missing", sep = "_"), AA = NA, Codon = NA, Infernal = NA, AlleleCM = NA, IsotypeScore = NA, Note = NA,
                 Classification = "Missing",
                 n.allele = n.missing, freq.ofcalled = NA, freq.ofallinclmissing = NA, 
                 value = n.missing/(n.missing + n.called))]
  setcolorder(formiss, names(pdat)[names(pdat)%in%names(formiss)]) # sometimes I added more columns after writing this fn....
  pdat<-rbind(pdat[, which(names(pdat)%in%names(formiss)), with = F], formiss) # Now the 'missing' alleles are added in!
  setkey(pdat, tRNA)
  
  ## Re-level data as needed
  pdat[, groupaa:=factor(groupaa, levels = aa.order)]
  pdat[, AA:=factor(AA, levels = aa.order)]
  classifcols.wmiss<-c(classifcols, missingcol)
  names(classifcols.wmiss)<-c(names(classifcols), "Missing")
  pdat[, Classification:=factor(Classification, levels = names(classifcols.wmiss))] # added in missing here......
  
  
  ## Handle if lowpad given
  if(lowpad>0){
    pdat[value<lowpad, value:=lowpad]
  }
  ## Labels for N alleles: blank if one allele, number if more than one [for now]
  pdat[,lab.n:=ifelse(NAlleles>1, as.character(NAlleles), "") ]
  
  # ---- Plotting
  # Generate base plot: faceted by amino acid THEN codon
  plt<-ggplot() + geom_scatterpie(data = pdat, 
                                  aes(x = x, y = y, group = tRNA),
                                  pie_scale=11,
                                  size=0, 			# This changes thickness of outline
                                  cols="Classification",
                                  long_format = T) +
    scale_fill_manual(values=classifcols.wmiss) +
    coord_equal() +
    labs(fill = "Allele classification:") + 
    theme_void() +
    theme(	strip.text.x = element_text(size=16, margin = margin(b=4, t=4)),
           strip.background = element_rect(color="black", size=1, linetype="solid"),
           legend.position = "bottom", legend.text = element_text(size = 14), legend.title = element_text(size = 15),
           plot.title = element_text(size = 18), plot.subtitle = element_text(size = 16))
  ## Add facets based on function argument
  if(myfacet == "both"){
    ##  facet by AA then codon
    plt<-plt + facet_nested(. ~ groupaa + Codon.orig, space = "free", switch = "x")  # scales = "free", `coord_fixed()` doesn't support free scales. not needed I think?  
  }else if(myfacet == "aaonly"){
    ## facet by AA 
    plt<-plt + facet_grid(.~groupaa, switch="x")
  }
  ## Label N alleles if desired
  if(labelnalleles==T){
    # labdat based on how facet? [have more columns in it]
    if(myfacet=="aaonly"){
      labdat<-unique(pdat[, .(tRNA, groupaa, x, y, lab.n)])
    }else if(myfacet=="both"){
      labdat<-unique(pdat[, .(tRNA, groupaa, Codon.orig, x, y, lab.n)])
    }
    
    plt<-plt + geom_text(data = labdat, 
                          aes(x, y, label = lab.n),
                         nudge_x = 0.25)
    # basics is working! need to figure out space etc
  }
  
  return(plt)
  
  
  # -------- old notes
  # TO DO - ordering/data manipulation
  #       if AA only, order by dominant allele type rather than codon? tricky; currently know at least that codons are grouped together
  # ** pies need 'slice' for missing data. add one fake 'allele' per tRNA
    
  #         **make sure FLAG where all were altered, i.e. genes that maybe should 'live' with different amino acid *** (all.altered column)
  #         OR remove these....they show up weird places
  # ALSO:
  #     annotate with # alleles ** add label
  # **Add 'missing' as white space to LEGEND. check this is working first.  

  # Issues I'm seeing in plot:
  #     multiple 'best' alleles
  #     not seeing missingness showing up (white space in pie)
  
  #### DO ALL THE THINGS NEEDED ABOVE #####
  
  # Could modify so codons get stacked/labeled going up....probably? [rather than horizontally]
  
  # Return
      # optionally also return data?
  
  # ...call out codons that have isotype switches....? how will this work here.... OK based on tRNA or....???
  #         codon should be one in REF strain? (but sometimes that's a switch?)
}

getbestaa<-function(oneg){
  # Gets the 'best' amino acid to classify each tRNA gene as: the AlleleCM of the majority of alleles/the ones that don't have isotype switch (if any)
  if(oneg[, all(all.altered)]){
    # Hard to say here, just pick one!
    myaa<-oneg[1, AlleleCM]
  }else{
    myaa<-oneg[AA==AlleleCM & freq.ofallinclmissing==max(freq.ofallinclmissing), AlleleCM][1]
  }
  return(myaa)
}

posvsnalplt<-function(twloc, chrlens, sorder, aa.order, mytitle = "", myscales = "free_x",
                      addhist = T, highlightgs = NA, highlightdescrip = "Pseudo. in ref"){
  # Makes plot that has chromosomal location on X axis, number of alleles on the y axis. IGNORES chromosomes not in chr info
  # In: twloc, data to plot from. Columns displayname, Chr, pos, Codon, numVersions required
  #     chrlens, data.table with chromosome information for each species. Columns: displayname, Chr, Length - one row per chromosome per species
  #     sorder, order for species to be plotted/analyzed/etc in. Display names in char vector.
  #     aa.order, order to plot amino acids in
  #     mytitle, main title for plot
  #     myscales, passed to facet_grid
  #     addhist, on top of points, plot line showing number of genes?
  #     highlightgs, optional data.table of columns displayname, tRNA - to highlight on plot with an additional point overlay
  #     highlightdescrip, only used if highlightgs provided. What text should come next to the point in the legend showing what's highlighted?
  # Out: plot
  
  # --- Data set up
  ## Exclude chrs not in chr info
  setkey(twloc, displayname, Chr)
  setkey(chrlens, displayname, Chr)
  pdata<-twloc[chrlens]
  exinfo<-paste("Genes not on included chromosomes:",
                paste(sapply(sorder, function(s){
                  n<-twloc[displayname==s, sum(!Chr%in%chrlens[displayname==s, Chr])]
                  paste(n, "in", s)
                }), collapse = "; "))
  ## Add fake data to span length of chromosomes
  fakes<-rbind(data.table(displayname = chrlens$displayname,
                    shortname = NA,
                    Chr = chrlens$Chr,
                    Start = NA, End = NA, tRNA = NA, Strand = NA, AA = NA, Codon = NA, 
                    numVersions = NA,
                    numStrains = NA, Missing = NA, 
                    pos = 1,
                    Length = NA),
               data.table(displayname = chrlens$displayname,
                          shortname = NA,
                          Chr = chrlens$Chr,
                          Start = NA, End = NA, tRNA = NA, Strand = NA, AA = NA, Codon = NA, 
                          numVersions = NA,
                          numStrains = NA, Missing = NA, 
                          pos = chrlens$Length,
                          Length = NA)
  )
  
  pdata<-rbind(pdata, fakes)
  ## order things nicely
  pdata[, AA:=factor(AA, levels = aa.order)]
  pdata[, displayname:=factor(displayname, levels = sorder)]
  ## plotting pos
  pdata[, posmb:=pos/1e06]
  
  ## Bounds of possible plotting area beyond each sample's bounds
  maxln<-chrlens[, max(Length), by = Chr]
  setnames(maxln, "V1", "max.ln")
  setkey(maxln, Chr)
  setkey(chrlens, Chr)
  chrends<-maxln[chrlens]
  
  # --- Plot
  myylab<-"Number of alleles in population"
  plt<-ggplot(data = pdata) + 
    geom_rect(data = chrends, aes(xmin = Length/1e06, xmax = max.ln/1e06, ymin = -Inf, ymax = Inf),
              fill = "gray20") # Set up rectangles to 'block out' areas of plot that don't exist in given chromosomes
  if(addhist==T){ # Add histogram,
    plt<-plt + geom_histogram(aes(posmb), fill = "lightgray", 
                              )  # adds number genes in a given region
        # breaks = seq(0, chrends[, max.ln]) does NOT work, can't seem to easily do diff breaks in diff facets
    # Deal with adding an AXIS for this if desired!!
    myylab<-"Number of genes (gray)\nNumber of alleles (colors)"
  }
  plt<-plt + geom_point(aes(x = posmb, y = numVersions, color = AA), pch = 3, alpha = 0.4) + # adds number alleles for each gene
    xlab("Genomic position (Mb)") + ylab(myylab) +
    ggtitle(mytitle, subtitle = exinfo) + 
    facet_grid(displayname ~ Chr, scales = myscales, space = "free_x") +
    scale_x_continuous(expand = c(0, 0)) + expand_limits(y = 0) +
    labs(color = "Amino acid") +
    myggtheme + theme(legend.position = "bottom", panel.grid = element_blank(),
                      strip.text.y = element_text(face = "italic"))
  
  # Highlight points if desired
  if(is.data.table(highlightgs)){
    # Set up data
    setkey(highlightgs, tRNA, displayname)
    setkey(pdata, tRNA, displayname)
    hdata<-pdata[highlightgs][!is.na(Chr)]
    
    # Add to plot
    plt<-plt + geom_point(data = hdata, aes(posmb, numVersions, shape = highlightdescrip), size = 0.7) +
      labs(shape = "")
  }

  # Maybe pick better colors for AAs?
  
  # --- Return
  return(plt)
}

posbyaa<-function(twloc, chrlens, spec = "C. elegans", aa.order, mytitle = "", mysubt = "",
                  highlightgs = NA, highlightdescrip = "Pseudo. in ref"){
  # Plots where tRNA genes are in genome split by amino acid. **based on AA in the twloc input which comes from trnas.gen - NOT AlleleCM, just what was in ref etc
  # in 500kb bins currently!
  #  In: twloc, data to plot from. Columns displayname, Chr, pos, Codon, AA, numVersions required
  #     chrlens,  data.table with chromosome information for each species. Columns: displayname, Chr, Length - one row per chromosome per species
  #     spec, which species should this plot be for: data is narrowed to this
  #     aa.order, order amino acid facets should occur in
  #     mytitle, title for plot
  #     mysubt, subtitle for plot
  #     highlightgs, optional data.table of columns displayname, tRNA - to highlight on plot with an additional histogram overlay (semi-transparent)
  #     highlightdescrip, only used if highlightgs provided. What text should come next to the point in the legend showing what's highlighted?
  # Out: plot
  
  # ---- Format data; 
  # keep only genes actually on chromosomes included here
  pdat<-twloc[displayname==spec]
  chrdat<-chrlens[displayname==spec]
  setkey(pdat, Chr, displayname)
  setkey(chrdat, Chr, displayname)
  pdat<-pdat[chrdat]
  # Fake data to span length of chromosomes
  fakes<-rbind(data.table(displayname = spec,
                          shortname = NA,
                          Chr = chrdat$Chr, # rep(chrdat$Chr, each = length(aa.order)),
                          Start = NA, End = NA, tRNA = NA, Strand = NA, 
                          AA = rep(aa.order, nrow(chrdat)), 
                          Codon = NA, 
                          numVersions = NA,
                          numStrains = NA, Missing = NA, 
                          pos = 1,
                          Length = NA),
               data.table(displayname = spec,
                          shortname = NA,
                          Chr = chrdat$Chr, # rep(chrdat$Chr, each = length(aa.order)),
                          Start = NA, End = NA, tRNA = NA, Strand = NA, 
                          AA =  NA, #rep(aa.order, nrow(chrdat)), 
                          Codon = NA, 
                          numVersions = NA,
                          numStrains = NA, Missing = NA, 
                          pos = chrdat$Length, # rep(chrdat$Length, each = length(aa.order)),
                          Length = NA)
  )
  pdat<-rbind(pdat, fakes)
  ## order things nicely
  pdat[, AA:=factor(AA, levels = aa.order)]
  ## plotting pos
  pdat[, posmb:=pos/1e06]

  # Make plot(s)
  binwd<-500e3/1e06
  plt<-ggplot() + geom_histogram(data = pdat, aes(posmb, fill = "Genes"), binwidth = binwd) +
    xlab("Genomic position (Mb)") + ylab("Number of tRNA genes") +
    ggtitle(mytitle, subtitle = mysubt) +
    facet_grid(AA~Chr, scales = "free_x", space = "free_x") +
    scale_fill_manual(values = "darkblue") + # for legend to look nice
    labs(fill = "") +
    scale_x_continuous(expand = c(0, 0)) + expand_limits(y = 0) +
    myggtheme + theme(panel.grid = element_blank(), panel.spacing = unit(0.1, "cm", data = NULL),
                      legend.position = "bottom")
  
  
  # New hist with highlight gene set if desired
  if(is.data.table(highlightgs)){
    # Set up data
    setkey(highlightgs, tRNA, displayname)
    setkey(pdat, tRNA, displayname)
    hdata<-pdat[highlightgs][!is.na(Chr)][displayname==spec]
    
    # Add to plot
    plt<-plt + geom_histogram(data = hdata, aes(posmb, fill = highlightdescrip), color = "red", lwd = 0.1, binwidth = binwd, alpha = 0.3) +
      scale_fill_manual(values = c("blue", "red")) + # for legend to look nice
      labs(fill = "")
  }
  
  ## ***Ideally also make version that is prop OF TOTAL N OF THIS AA (within facet)....would be nice if that could be done automatically....
  ##            a bit tricky - can do within facet easily , but not for the whole ROW....unless I compute that myself....which I don't want to do right now
  
  # Return
  return(plt)
}

testlocdist<-function(twloc,  splitcol = "refpseud", splitdescrip = "Pseudo. in ref", chrlens,
                      sorder, mytitle = "", myscales = "free"){
  # Tests whether location distribution between two sets of tRNA genes is different with ks.test (per chromosome [and overall?]); also plots these
  #     FACETED by chr, species
  # In: twloc, data.table with tRNA location AND T/F column named splitcol. must have columns displayname, Chr, pos, <splitcol>
  #     splitcol, name of boolean (T/F) column to split data on and test position within
  #     splitdescrip, description of splitcol for outputs, plotting, etc
  #     chrlens, data.table with chromosome information for each species. Columns: displayname, Chr, Length - one row per chromosome per species
  #     sorder, order for species to be plotted/analyzed/etc in. Display names in char vector.
  #     mytitle, main title for plot
  #     myscales, scales for facets
  # Out: List of data and plots:
  #       $statres, KS test results for each chromosome and genomewide (points given genome-wide coordinates). Columns:
  #           displayname, species as in input
  #           Chr, all (for genome wide) or chromosome this result is for
  #           tested, description of what's tested - from splitdescrip input
  #           ks.D.stat, self-explanatory
  #           ks.p.value, self-explanatory
  #       $histplt, overlapping transparent histograms showing proportion of tRNA genes in each 500kb gene genome-wide, colored by splitcol. KS p-values annotated.
  #       $ptplt, plots where each gene is one 'tick' style point, colored and located by splitcol. KS p-values annotated.
  
  # --- Do KS tests
  ## Per chr
  ksres.chr<-rbindlist(lapply(sorder, function(sp){
    out<-rbindlist(lapply(chrlens[displayname==sp, Chr], function(mychr){
      res<-ks.test(twloc[displayname==sp & Chr==mychr & get(splitcol)==T, pos], 
                   twloc[displayname==sp & Chr==mychr & get(splitcol)==F, pos])
      out<-data.table(displayname = sp,
                      Chr = mychr,
                      tested = splitdescrip,
                      ks.D.stat = res$statistic,
                      ks.p.value = res$p.value)
      return(out)
    }))
    return(out)
  }))
  ## Genome-wide coords for power
  ksres.all<-rbindlist(lapply(sorder, function(sp){
    thisdat<-twloc[displayname==sp, ]
    chrdat<-chrlens[displayname==sp, ]
    
    # Get fake-out coords. ASSUMES CHRS ARE IN ORDER
    thisdat[Chr==chrdat[1, Chr], gpos:=pos]  
    for(i in 2:nrow(chrdat)){
      myadd<-chrdat[1:(i -1), sum(Length)]
      thisdat[Chr==chrdat[i, Chr], gpos:=pos + myadd]
    }
    
    # Do test
    res<-ks.test(thisdat[get(splitcol)==T, gpos],
                 thisdat[get(splitcol)==F, gpos])
    out<-data.table(displayname = sp,
                    Chr = "all_combined",
                    tested = splitdescrip,
                    ks.D.stat = res$statistic,
                    ks.p.value = res$p.value)
    return(out)
  }))
  
  # --- Format data for plotting
  pdat<-copy(twloc)
  pdat[, usecol:=get(splitcol)]
  pdat<-pdat[, .(displayname, Chr, pos, usecol)]
  ## Exclude chrs not in chr info
  setkey(pdat, displayname, Chr)
  setkey(chrlens, displayname, Chr)
  pdat<-pdat[chrlens]
  ## Add fake data to span length of chromosomes
  fakes<-rbind(data.table(displayname = chrlens$displayname,
                          Chr = chrlens$Chr, pos = 1,
                          usecol = NA, Length = NA),
               data.table(displayname = chrlens$displayname,
                          Chr = chrlens$Chr, pos = chrlens$Length,
                          usecol = NA, Length = NA)
  )
  
  pdat<-rbind(pdat, fakes)
  ## order things nicely
  pdat[, displayname:=factor(displayname, levels = sorder)]
  ## plotting pos
  pdat[, posmb:=pos/1e06]
  ## y position for tick plot for optimal changeability
  pdat[usecol==F, ticky:=0]
  pdat[usecol==T, ticky:= 0.1]
  
  ## Bounds of possible plotting area beyond each sample's bounds
  maxln<-chrlens[, max(Length), by = Chr]
  setnames(maxln, "V1", "max.ln")
  setkey(maxln, Chr)
  setkey(chrlens, Chr)
  chrends<-maxln[chrlens]
  
  # P-value labels
  ## per chr, for facets
  kslab<-copy(ksres.chr)
  kslab[, mylabel:=paste("KS", "~italic(p)==~", ifelse(ks.p.value<0.001, sci2carrot(ks.p.value), round(ks.p.value, digits = 3)))]
  ## overall, for subtitle
  subt<-paste0("Genome-wide KS p-values:\n", 
              ksres.all[, paste(displayname, ifelse(ks.p.value<0.001, paste(formattable::scientific(x = ks.p.value, digits = 1)), round(ks.p.value, 3)), 
                                sep = ": ", collapse = "; ")])
  #       scientific formatting NOT working in here; it's not doing the digits even though normally does and does when given one value...paste fixed it for mysterious reasons
  
  # --- Make plot(s)
  # Try 1: overlapping histogram
  bnsz<-500e03/1e06

  # hists need to be plotted separately for overlap to work nicely!! (fill = usecol and all data segments bars/colors segments of them, doesn't plot whole deal)
  ovhist<-ggplot() + 
    geom_rect(data = chrends, aes(xmin = Length/1e06, xmax = max.ln/1e06, ymin = -Inf, ymax = Inf),
              fill = "gray20") + # Set up rectangles to 'block out' areas of plot that don't exist in given chromosomes
    geom_histogram(data = pdat[usecol==F], aes(posmb, y = after_stat(density)*bnsz, fill = usecol), binwidth = bnsz) +
    geom_histogram(data = pdat[usecol==T], aes(posmb, y = after_stat(density)*bnsz, fill = usecol), binwidth = bnsz) +
    geom_histogram(data = pdat[is.na(usecol)], aes(posmb, y = after_stat(density)*bnsz, fill = usecol), binwidth = bnsz) +
    scale_fill_manual(values = c(rgb(0,0,1,0.5), rgb(1,0,0,0.5)), na.value = rgb(1,1,1,0), labels = c("No", "Yes", "")) +
    facet_grid(displayname ~ Chr, scales = myscales, space = "free_x") +
    scale_x_continuous(expand = c(0, 0)) + expand_limits(y = 0) +
    labs(fill = splitdescrip) +
    ylab("Proportion of tRNA genes") + xlab("Position in genome (Mb)") +
    ggtitle(mytitle, subtitle = subt) +
    myggtheme + theme(legend.position = "bottom", panel.grid = element_blank(),
                      strip.text.y = element_text(face = "italic"), panel.spacing = unit(0.1, "cm", data = NULL)) +
    geom_label(size = 3, data = kslab, mapping = aes(x = Inf, y = Inf, label = mylabel),
               hjust = 1.05, vjust = 1.5, parse = T)
  
  # Another version: points/vertical lines where actual data is
  tickplt<-ggplot() +
    geom_rect(data = chrends, aes(xmin = Length/1e06, xmax = max.ln/1e06, ymin = -Inf, ymax = Inf),
              fill = "gray20") + # Set up rectangles to 'block out' areas of plot that don't exist in given chromosomes
    geom_point(data = pdat, aes(x = posmb, y = ticky, color = usecol),  pch = "|", size = 4) +
    scale_color_manual(values = c(rgb(0, 0, 0, 0.5), rgb(1, 0, 0, 0.5)), na.value = rgb(1,1,1,0), labels = c("No", "Yes", "")) +
    facet_grid(displayname ~ Chr, scales = myscales, space = "free_x") +
    scale_x_continuous(expand = c(0, 0)) + 
    ylim(pdat[, min(ticky, na.rm = T)] - 3*pdat[, max(ticky, na.rm = T)], pdat[, max(ticky, na.rm = T)] + 3*pdat[, max(ticky, na.rm = T)]) + # if min stays 0 this in theory centers the thing
    labs(col = splitdescrip) +
    ylab("Proportion of tRNA genes") + xlab("Position in genome (Mb)") +
    ggtitle(mytitle, subtitle = subt) +
    myggtheme + theme(legend.position = "bottom", panel.grid = element_blank(),
                      strip.text.y = element_text(face = "italic"), panel.spacing = unit(0.1, "cm", data = NULL),
                      axis.ticks.y = element_blank(), axis.text.y = element_blank()) +
    geom_label(size = 3, data = kslab, mapping = aes(x = Inf, y = Inf, label = mylabel),
               hjust = 1.05, vjust = 1.5, parse = T)

  # --- Return (stats & plot)
  return(list(statres = rbind(ksres.all, ksres.chr),
              histplt = ovhist,
              ptplt = tickplt))
}

# alleletreemap<-function(pdat, classifcols){
#   # Trying treemap. IN DEVELOPMENT - DELETE/CHANGE IF NOT USING 
#   
#   # Data setup
#   pdat<-copy(pdat)
#   pdat[, Classification:=factor(Classification, levels = names(classifcols))]
#   
#   # Make plot
#   tm<-ggplot(pdat, aes(area = freq.ofallinclmissing, fill = Classification, subgroup = Classification)) + 
#     geom_treemap() +  geom_treemap_subgroup_border(color = NA) +
#     scale_fill_manual(values = classifcols) +
#     facet_nested_wrap(. ~ AA + Codon + tRNA, scales = "free") + # facet_nested(. ~ AA + Codon + tRNA, scales = "free", space = "free") +
#     myggtheme
#   
#   # Current issues: 
#   #       not split by gene. maybe facet can do that....It DOES! but then it's kinda just back to barchart issue
#   #             facet_nested_wrap makes this a bit better, but wraps all the facets (when I'd want it to keep everything within an amino acid together...)
#   #       **MISSING space not automatically added. Probably important. [could add it to data if needed, of course...]
#   #       classification that's most common goes first - lose factor ordering
# }

#### Arguments & inputs ####
p<-arg_parser("Initial Visualize/analyze tRNA allele information ", 
              name = "show_mutational_variation.R", hide.opts = TRUE)

# Input file related
p<-add_argument(p, "--speciesf",
                help = "File containing information on all species to process here. Columns infilename (exactly how all files have this species in their name),
                displayname (name that should be used for plot outputs etc), shortname (no-spaces name for ouptut files, sorting, etc - either shorter than or same as infilename, probably).
                In order you'd like plots to be in!",
                type = "character")
p<-add_argument(p, "--alleleinfo",
                help = "EXAMPLE path to *alleleinfo_counts_wmissing.txt output of getstrainxtrnacalls.R. Where species ID/species specific info is, put SAMP instead",
                type = "character")
p<-add_argument(p, "--trnasgen",
                help = "EXAMPLE path to *strain_trnas_gen.txt output of build_alt_sequences.py. Where species ID/species specific info is, put SAMP instead",
                type = "character")
p<-add_argument(p, "--trnainfo",
                help = "EXAMPLE path to *_bygeneup_tRNA_counts_wmissing.txt output of getstrainxtrnacalls.R. (Columns about tRNA; counts of alleles w/ various classifications.) 
                 Where species ID/species specific info is, put SAMP instead",
                type = "character")
p<-add_argument(p, "--chrlens",
                help = "File containing lengths of chromosomes to use for plotting; chromosomes/contigs not included from those analyses/plots. Columns displayname (matching species info in speciesf), Chr, Length",
                type = "character")

# Output file related
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

#### Read in data ####
cat("....Reading in data....\n")
sinfo<-fread(p$speciesf, header = T)

# Allele info
peral<-fread.mult(sinfo, exampfile = p$alleleinfo)
## Add proportions/allele frequencies
peral[, `:=`(freq.ofcalled = n.allele/n.called,
             freq.ofallinclmissing = n.allele/(n.called + n.missing))]

# tRNA summary info
pertrna<-fread.mult(sinfo, exampfile = p$trnainfo)

# tRNA location info
twloc<-fread.mult(sinfo, exampfile = p$trnasgen)
if("##Chr"%in%names(twloc)){
  setnames(twloc, "##Chr", "Chr")
}
twloc[, pos:=((Start + End)/2)] # add single midpoint loc for plotting

# chromosomes info
chrlens<-fread(p$chrlens, header = T)

#### Do some analyses ####
# Ref pseudogenes that are functional in other alleles?
refpseud.dat<-refpseud(peral, sinfo)
write.table(refpseud.dat$pseudorefgs, file.path(p$outdir, paste0(p$baseoutname, "_pseudoinref_genes.txt")),
            sep = "\t", quote = F, row.names = F)
write.table(refpseud.dat$summ, file.path(p$outdir, paste0(p$baseoutname, "_pseudoinref_summary.txt")),
            sep = "\t", quote = F, row.names = F)
write.table(refpseud.dat$allpseudoandfunc, file.path(p$outdir, paste0(p$baseoutname, "_pseudoinref_details_putativelyfunc.txt")),
            sep = "\t", quote = F, row.names = F)

#### Do some plotting: allele info only ####
cat("....Allele plots....\n")
# output dir
aplotdir<-file.path(p$outdir, "alleleinfoplots")
if(!dir.exists(aplotdir)){dir.create(aplotdir)}

# Plot allele info
classifcols<-c("skyblue", "skyblue4", "yellow3", "red4")
names(classifcols)<-c("Best", "Functional", "Altered", "Lost")

# Add the codon ones isotype switches *should've* been (to enable splitting plot by codon if want)       
#      AND adding orig AA here for the few where AlleleCM differs within gene (also done by allelescatterpie fn)
peral<-codoncm(peral)
setkey(peral, tRNA, shortname)
peral<-peral[, getbestaa(.SD), by = .(tRNA, shortname)][peral]
setnames(peral, "V1", "AlleleCM.orig")
setcolorder(peral, c("displayname", "shortname", "tRNA", "allelename", "AA", "Codon", "AlleleCM", "Codon.orig", "AlleleCM.orig"))

write.table(peral, file.path(p$outdir, paste0(p$baseoutname, "_alleleinfo_counts_worigcodonetc.txt")),
            sep = "\t", quote = F, row.names = F) # Save

## Doing species differently here because they have diff tRNA genes (did technically work with facet_nested(displayname ~ AA + Codon, scales = "free", space = "free"))
albars<-lapply(sinfo$displayname, function(x){
  allelebar(peral[displayname==x], classifcols) + ggtitle(x)
  })
pdf(file.path(aplotdir, paste0(p$baseoutname, "_tRNA_allele_barplots.pdf")), 64, 16)
invisible(lapply(albars, print))
invisible(dev.off())

## Stacked up pie charts
    # VERSIONS to make: [with titles]
    # species separately since will be weird with facets. X
    # pseudo in all strains included, excluded [refpseud.dat$pseudorefgs]
    # with alleles numbered and without X
    # split out by codon - width 2 (on top of just by AA - width 3) X
### Set up 
al.labs<-c(F, T)
myfacs<-data.table(myfac = c("aaonly", "both"),
                   mywd = c(3, 2))
aa.order<-peral[, sort(unique(AlleleCM))]
if("Undet"%in%aa.order){
  aa.order<-c(aa.order[aa.order!="Undet"], "Undet")
}

### Make plots
stpies<-lapply(1:nrow(sinfo), function(s){

  # All genes
  opdat<-peral[shortname==sinfo[s, shortname], ] # no gene narrowing here
  allps<-lapply(al.labs, function(al){
    labps<-lapply(1:nrow(myfacs), function(fi){
      plt<-allelescatterpie(peral.one = opdat, classifcols = classifcols, missingcol = "white", aa.order = aa.order,
                       lowpad = 0.03, percodonwd = myfacs[fi, mywd], myfacet = myfacs[fi, myfac],
                       labelnalleles = al) +
        ggtitle(sinfo[s, displayname], subtitle = "All tRNAs in reference - genes and pseudogenes")
      return(plt)
    })
    return(labps)
  })

  # Excluding genes pseud in ref
  opdat<-peral[shortname==sinfo[s, shortname] & !tRNA%in%refpseud.dat$pseudorefgs[displayname==sinfo[s, displayname], tRNA], ] # no gene narrowing here
  allps.nopseud<-lapply(al.labs, function(al){
    labps<-lapply(1:nrow(myfacs), function(fi){
      plt<-allelescatterpie(peral.one = opdat, classifcols = classifcols, missingcol = "white", aa.order = aa.order,
                            lowpad = 0.03, percodonwd = myfacs[fi, mywd], myfacet = myfacs[fi, myfac],
                            labelnalleles = al) +
        ggtitle(sinfo[s, displayname], subtitle = "tRNAs in reference, excluding pseudogenes")
      return(plt)
    })
    return(labps)
  })
  
  # Return plots?
  return(list(allps, allps.nopseud))
})

# Save in 2 PDFs per species [doing separately so can mess with PDF settings without re-making plots]
## Split by amino acid
myhts<-c(8, 12, 7) ## DANGER WILL ROBINSON: hard coded PDF heights
invisible(
lapply(1:nrow(sinfo), function(s){
  theseps<-stpies[[s]]
  # determine height of PDF based on max height of stacks...how know that?....N genes....NOT as simple right now, will consider later
  pdf(file.path(aplotdir, paste0(p$baseoutname, "_", sinfo[s, shortname], "_tRNA_allele_stackedpies_AASplit.pdf")),
      24, myhts[s])
  
  lapply(theseps, function(labps){
    lapply(labps, function(sps){
      print(sps[[1]])
      return(NULL)
    })
    return(NULL)
  })
  
  invisible(dev.off())
  return(NULL)
})
)

## Split by amino acid and codon [needs to be wider to see anything]
invisible(
  lapply(1:nrow(sinfo), function(s){
    theseps<-stpies[[s]]
    pdf(file.path(aplotdir, paste0(p$baseoutname, "_", sinfo[s, shortname], "_tRNA_allele_stackedpies_AAPlusCodonSplit.pdf")),
        48, 8)
    
    lapply(theseps, function(labps){
      lapply(labps, function(sps){
        print(sps[[2]])
        return(NULL)
      })
      return(NULL)
    })
    
    invisible(dev.off())
    return(NULL)
  })
)
# By codon being kinda weird with too many white circles?


#### Do some plotting: with location in genome information ####
cat("....Location plots....\n")
# output dir
lplotdir<-file.path(p$outdir, "genelocplusplots")
if(!dir.exists(lplotdir)){dir.create(lplotdir)}

# --- Genomic location vs. # tRNA alleles (anytime there's a dot or whatever that's where a tRNA gene is) ** could COLOR by amino acid! [get orig amino acid - from stacked pie fn potentially]
# Versions to make [in order]:
#     with and without gene count histogram
#     Normal and log scale y axis - add scale_y_log10() to existing output for ease
#     Pseudogenes highlighted and not
#     Free and fixed y axis scales - for seeing more of data and for comparing across species
## plot variable set up
pseudhighlight<-refpseud.dat$pseudorefgs[, .(displayname, tRNA)]
myscdt<-data.table(passval = c("free", "free_x"), # pass to plot
                   titleval = c("free", "fixed")) # describes y axis
maxln<-chrlens[, max(Length), by = Chr]
setnames(maxln, "V1", "max.ln")
setkey(maxln, Chr)
setkey(chrlens, Chr)
chrends<-maxln[chrlens]

## make plot
gtplts<-lapply(c(T, F), function(dohist){
  normax<-unlist(
    lapply(list(NA, pseudhighlight), function(myhighlights){
    plt.list<-
      lapply(1:nrow(myscdt), function(mysci){ # always want x to be free; y needs fixed
      usetitle<-paste0("tRNA gene location vs. number alleles",
                      ifelse(dohist, " and number of genes", ""),
                      "; ", myscdt[mysci, titleval], " scale y axes",
                       ifelse(is.data.table(myhighlights), ";\nreference pseudogenes highlighted", ";\n"))
      
      plt<-posvsnalplt(twloc, chrlens, sorder = sinfo$displayname, aa.order = aa.order, 
                  mytitle = usetitle, myscales = myscdt[mysci, passval], addhist = dohist, 
                  highlightgs = myhighlights, highlightdescrip = "Pseudo. in ref")
      return(plt)
    })
    return(plt.list) ## unlisting may not work here, need to check
  }), recursive = F)
  
  # ADD scale to each of these pre-generated plots....and add log y ax to title?
  #       And need to fix rectangle to work; just doing manually
  logyax<-lapply(normax, function(oneplt){
    mytitle<-paste0(oneplt$labels$title, "; log scale y axis")
    # y_limits_1<-ggplot_build(oneplt)$layout$panel_params[[1]]$y.range
    
    plt<-oneplt + scale_y_log10()
    # y_limits <- ggplot_build(plt)$layout$panel_params[[1]]$y.range
    plt<-plt +
      ggtitle(mytitle) + 
      geom_rect(data = chrends, aes(xmin = Length/1e06, xmax = max.ln/1e06, ymin = 0.9 , ymax = Inf), fill = "gray20") 
    return(plt)
  })
  return(list(normax = normax, logyax = logyax))
})

# Save [separate so can mess with formatting separate from plot generation...]
pdf(file.path(lplotdir, paste0(p$baseoutname, "_tRNA_locVsNAllelesEtc.pdf")), 10, 8)
invisible(lapply(gtplts, function(oneset){
  lapply(oneset, function(one4){
    lapply(one4, print)
  })
}))
invisible(dev.off())



# ---- tRNA location split by amino acid
# Versions:
#     with and without pseuodogenes seperately highlighted
#     with and without pseudogenes included at all
#     possibly...multiple definitions of pseudogenes? Like, any with any pseud allele (plus what already easily have - pseud in ref)

# New highlight data: genes with pseud alleles in any strain
anypseud<-pertrna[Lost>0, ]

# Generate plots
aalocplts<-lapply(sinfo$displayname, function(myspec){
  # With pseudogenes not highlighted
  allgs<-posbyaa(twloc, chrlens, spec = myspec, aa.order = aa.order, mytitle = myspec,
                 mysubt = "All genes, no highlights")
  
  nopdat<-twloc[displayname==myspec & !tRNA%in%pseudhighlight[displayname==myspec, tRNA]]
  nops<-posbyaa(nopdat, chrlens, spec = myspec, aa.order = aa.order, mytitle = myspec,
                mysubt = "Genes pseudo. in ref exlcuded")
  
  # Pseudogenes highlighted: 2 different ways
  ## pseudo in ref
  pref<-posbyaa(twloc, chrlens, spec = myspec, aa.order = aa.order, mytitle = myspec,
                mysubt = "All genes, pseudo in ref highlighted", highlightgs = pseudhighlight,
                highlightdescrip = "Pseudo. in ref")
  ## pseudo in ANY strain
  pany<-posbyaa(twloc, chrlens, spec = myspec, aa.order = aa.order, mytitle = myspec,
                 mysubt = "All genes, w/ pseudo alleles highlighted", highlightgs = anypseud,
                 highlightdescrip = "Pseudo. in any 1+ strains")
  
  # Return
  return(list(allgs, nops, pref, pany))
})
names(aalocplts)<-sinfo$displayname

# Save plots: all in one, with species titles on blank pages between species
pdf(file.path(lplotdir, paste0(p$baseoutname, "_tRNA_locByAA.pdf")), 8, 12)
invisible(lapply(sinfo$displayname, function(myspec){
  titleplt<-ggplot() + geom_text(aes(x = 1, y = 1, label = myspec), size = 16, fontface = "italic") + theme_void() # Title/blank page
  print(titleplt)
  
  lapply(aalocplts[[myspec]], print)
  
  return(NULL)
}))
invisible(dev.off())

#### Various splits of tRNAs and looking at whether location distributions differ based on these ####
# --- Add info to test about ---
for(sp in sinfo$displayname){
  # Pseudo in ref? 
  twloc[displayname==sp & tRNA%in%refpseud.dat$pseudorefgs[displayname==sp, tRNA], refpseud:=T]
  twloc[displayname==sp & !tRNA%in%refpseud.dat$pseudorefgs[displayname==sp, tRNA], refpseud:=F]
  
  # Any pseudo/lost alleles?
  twloc[displayname==sp & tRNA%in%anypseud[displayname==sp, tRNA], anylost:=T]
  twloc[displayname==sp & !tRNA%in%anypseud[displayname==sp, tRNA], anylost:=F]
  
  # Variable in population?
  twloc[displayname==sp & tRNA%in%pertrna[displayname==sp & VariableInPop==T, tRNA], variableInPop:=T]
  twloc[displayname==sp & !tRNA%in%pertrna[displayname==sp & VariableInPop==T, tRNA], variableInPop:=F]
  
  # Any strains with missing genotypes?
  twloc[displayname==sp & tRNA%in%pertrna[displayname==sp & anyMissingCalls==T, tRNA], anyMissingCalls:=T]
  twloc[displayname==sp & !tRNA%in%pertrna[displayname==sp & anyMissingCalls==T, tRNA], anyMissingCalls:=F]
  
  # More that 10 (?) strains with missing genotypes?
  #     possibly interesting but don't have that super easily right now
}

# --- Run tests, plots ---
# Set up test info
ksinfo<-data.table(colname = c("refpseud", "anylost", "variableInPop", "anyMissingCalls"),
                   descrip = c("Pseudo. in ref.?", "Any pseudo/lost alleles?",
                               "Variable in population?", "Any strains have any missing genotype calls?"))

# Run analyses
locsplitout<-lapply(1:nrow(ksinfo), function(i){
  testlocdist(twloc, splitcol = ksinfo[i, colname], splitdescrip = ksinfo[i, descrip],
              chrlens = chrlens, sorder = sinfo$displayname, mytitle = ksinfo[i, descrip],
              myscales = "free_x")
})

# Save data, plots
write.table(rbindlist(lapply(locsplitout, function(x) x$statres)),
            file.path(lplotdir, paste0(p$baseoutname, "_tRNA_groupings_location_KStests.txt")),
            sep ="\t", quote = F, row.names = F)

pdf(file.path(lplotdir, paste0(p$baseoutname, "_tRNA_groupings_location_KStests.pdf")),
    10, 7)
invisible(lapply(1:nrow(ksinfo), function(x){
  # Title page between tested categories
  titleplt<-ggplot() + geom_text(aes(x = 1, y = 1, label = ksinfo[x, descrip]), size = 16, fontface = "italic") + theme_void() # Title/blank page
  print(titleplt)
  
  # Actual plots
  print(locsplitout[[x]]$histplt)
  print(locsplitout[[x]]$ptplt)
}))
invisible(dev.off())

#### Misc notes ####
    #[make plots ## do things including & excluding tRNAs that are pseudogenized in all (only 'pseudo' alleles)]
# To facet by species, try facet_nested(displayname ~ AA + Codon, scales = "free", space = "free")

# Don't forget
## do things including & excluding tRNAs that are pseudogenized in all (only 'pseudo' alleles) AND 'all.altered' things - they also muddy
## do everything so can either facet of make multiple plots for multiple species
## show where tRNAs are in genome - will need other file(s)?

# Version where different alleles are just colors a la Annalise's plot - not saying anything about function?

# * plot ideas to make*
#     focus/call out isotype switching


#### Script completion message & session information ####
cat("....show_mutational_variation.R processing complete! Session information:....\n")
sessionInfo()
