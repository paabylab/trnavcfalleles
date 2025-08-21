#! /usr/bin/env/ Rscript
# Plot/analyze tRNA flank and gene body variation together (after flankvariation.R, genebodyvariation.R)
# By Avery Davis Bell, begun 2025.04.28
require(argparser, quietly = T)
require(data.table, quietly = T)
require(ggplot2, quietly = T)
require(RColorBrewer, quietly = T)
require(R.utils)
require(ggpattern)

#### Functions #### 
myggtheme<-theme_bw() +
  theme(axis.title.x = element_text(size = 12), axis.text.x = element_text(size = 11), 
        axis.title.y = element_text(size = 12), axis.text.y = element_text(size = 11),
        legend.title = element_text(size = 12), legend.text = element_text(size=11), 
        strip.text.x = element_text(size = 11), strip.text.y = element_text(size = 11), title = element_text(size = 14), 
        strip.text.x.top = element_text(size = 11), strip.text.x.bottom = element_text(size = 11), strip.text.y.right = element_text(size = 11), strip.text.y.left = element_text(size = 11),
        plot.subtitle = element_text(size = 13))

movingavg<-function(onef, n = 5, breakavg = "flank", colpos = "relpos", colsmooth = "p.tRNAs.var"){
  # Generates a moving average over n positions in colpos
  # CENTERS avg on each base pair - for first n/2 (and last), is off center - bp 1 is 1-5 average, as are 2 and 3; 4 is 2-6 avg
  # In: onef, data - all used to do moving average
  #     n, # positions to smooth over. Only tested if this is odd
  #     breakavg, name of column to only do avg within values - usually flank
  #     colpos, char name of column with positions to smooth over 
  #     colsmooth, char column name of column with values to do the average in
  # Out: data.table with one row per input. Columns:
  #     <breakavg>, as in input (keyed by this now)
  #     <colpos>, as in input (though may be re-sorted in some cases)
  #     avg, key output value - moving average of colpos (CENTERS avg on each base pair - for first n/2 (and last), is off center - bp 1 is 1-5 average, as are 2 and 3; 4 is 2-6 avg)
  
  # --- Subfunction
  oneavg<-function(tmpdat){
    tmpdat<-tmpdat[order(get(colpos)), ]
    tmpdat[, myn:=1:nrow(tmpdat)]
    # First and last n * ### doesn't need to be this many: just until can CENTER on it
    mysep<-floor(n/2)
    tmpdat[myn%in%c(1:mysep), avg:=tmpdat[1:n, mean(get(colsmooth))]] # mean is of first 5, index is just first 2
    tmpdat[myn%in%seq((max(myn)-mysep + 1), max(myn)), avg:=tmpdat[seq((max(myn)-n + 1), max(myn)), mean(get(colsmooth))]] # mean is of last 5, index is just last 2
    # All the others - avg centered on it
    avgs<-sapply((mysep+1):tmpdat[, max(myn)-mysep], function(i){
      return(tmpdat[(i - mysep):(i + mysep), mean(get(colsmooth))])
    })
    tmpdat[is.na(avg), avg:=avgs]
    
    out<-tmpdat[, .(get(colpos), avg)]
    setnames(out, c(colpos, "avg"))
    
    return(out)
  }
  
  # --- Data set up
  bdat<-copy(onef)
  setkeyv(bdat, breakavg)
  
  # --- Get averages
  out<-bdat[, oneavg(.SD), by = breakavg]
  
  # --- Return
  return(out)
}

flankgenelocplot<-function(onef, oneg, tstruct,
                           flankx = "relpos", flank = "flank", flanky = "p.tRNAs.var", flankymin = "low95ci.tRNAs.var", flankymax = "high95ci.tRNAs.var",
                           geney = "p.varsPerBp", geneymin = "low95ci.varsPerBp", geneymax = "high95ci.varsPerBp",
                           colcol = "mutclass", colvec = NULL, labcol = "", addcis = FALSE,
                           xlab = "Position in tRNA gene region",
                           ylab = "Polymorphism (prop. tRNAs with variation)",
                           drawinnerflank = T, inner5 = -20, inner3 = 10, 
                           drawsubstructs = T, 
                           facs = NULL,
                           flankavg = F, flankavgn = 5){
  # Makes plot that shows tRNA locus including flank on either side. Cobbles together plotting fns from flankvariation.R, genebodyvariation.R
  # In: onef, data on FLANKS to plot. For one gene set (or one + facets). Needs columns of values flankx, flank, flanky, flankymax, flankymin [generally defaults], colcol
  #     oneg, data on IN GENE BODY to plot. For one gene set (or one + facets). Needs columns of values geney, geneymin, geneymax, colcol
  #           plus structure, substructure, nSubStructure, canonicalnbp - for figuring out x axis where to plot
  #     tstruct, structure info about plotting. Must have structure, substructure as in dat, also structure_plot and substructure_plot, structure_plot_level & substructure_plot_level (used here for naming and ORDERING for legend etc)
  #     colcol, column to color lines by 
  #     colvec, optional - values for scale color manual (lines) - names are colcol values
  #     labcol, what label to give color in legend
  #     addcis, plot where 95% CI of each estimate is
  #     xlab, x axis label
  #     ylab, y axis label
  #     drawinnerflank, add lines where inner flank bounds are?
  #     inner5, inner3: where to add inner flank boundaries if drawinnerflank = T
  #     drawsubstructs, if T, substructures are highlighted/called out/something....
  #     facs: if going to facet LATER (does NOT facet here), tell the name of those facets - x axis computation has to happen internally. NULL or char vec
  #           ALSO need to do this for columns that are split into diff colors or whatever
  #     flankavg, if T, does flankavgn MOVING AVERAGE value for values in flanks instead of absolute values. **addcis must be F if this is used
  #     flankavgn, number bp to do flank moving average over
  
  # --- Format GENE BODY data for plotting
  # General
  gpdat<-copy(oneg)
  gpdat[, y :=get(geney)]
  gpdat[, ymin :=get(geneymin)]
  gpdat[, ymax:=get(geneymax)]
  gpdat[, coldat := get(colcol)]
  if(!is.null(colvec)){
    gpdat[, coldat:=factor(coldat, levels = names(colvec))]
  }
  # Structure stuff for plotting
  gpdat<-merge(gpdat, tstruct)
  # gpdat[is.na(substructure_plot), substructure_plot:="other"] better to leave NA so it doesn't get assigned a pattern
  gpdat[, structure_plot:=factor(structure_plot, levels = unique(tstruct[order(structure_plot_level), structure_plot]))]
  gpdat[, substructure_plot:=factor(substructure_plot, levels = c(na.omit(unique(tstruct[order(substructure_plot_level), substructure_plot]))))]
  # X axis locations (start and end of line segment) ** WITHIN ONE FACET/mutation type
  if(!is.null(facs)){
    usefacs<-c("coldat", facs)
  }else{
    usefacs<-"coldat"
  }
  setkeyv(gpdat, usefacs)
  gpdat<-gpdat[order(nSubStructure)]
  gpdat[,xend:=cumsum(canonicalnbp), by = usefacs]
  gpdat[,xstart:= c(1, xend[1:(.N -1)] + 1), by = usefacs]
  gpdat[, xseg:=xstart - 0.5]
  fakend<-gpdat[nSubStructure==max(nSubStructure)]
  fakend[, xseg:=xend + 0.5]
  gpdat<-rbind(gpdat, fakend) # add an 'end point' for final [1bp] structure
  
  # --- Format FLANK data for plotting
  fpdat<-data.table(x = onef[, get(flankx)], y = onef[,get(flanky)], ymin = onef[,get(flankymin)], ymax = onef[, get(flankymax)],
                   coldat = onef[,get(colcol)], flank = onef[, flank],
                   onef[, .SD, .SDcols = which(!names(onef)%in%c(flankx, flanky, flankymin, flankymax, colcol, "flank"))])
  if(!is.null(colvec)){
    fpdat[, coldat:=factor(coldat, levels = names(colvec))]
  }
  ## Add space for gene body
  fpdat[, xlab:=x]
  fpdat[flank=="3'", x:=x + gpdat[,max(xend)] + 1]
  ## Do moving avg if desired
  if(flankavg==T){
    if(addcis==T){
      stop("Can't do flank averaging when addcis==T, change flankavg or addcis to F")
    }
    setkeyv(fpdat, usefacs)
    # Get new (averaged) y values
    newy<-fpdat[,movingavg(onef = .SD, n = flankavgn, breakavg = "flank", colpos = "xlab", colsmooth = "y") , by = usefacs]
    # Replace old y values with these (merge in)
    setkeyv(newy, c(usefacs, "flank", "xlab"))
    setkeyv(fpdat, c(usefacs, "flank", "xlab"))
    fpdat<-merge(newy, fpdat)
    fpdat[, y:=avg]
  }
  
  # --- Connecting flank to gene body: need to do by facet, color split
  ## For each color split, facets if provided
  setkeyv(fpdat, usefacs)
  setkeyv(gpdat, usefacs)
  
  lconn<-rbind(fpdat[flank=="5'" & x == 0, x, by = usefacs][
    fpdat[flank=="5'" & x==0, y, by = usefacs]][
      fpdat[flank=="5'" & x==0, ymin, by = usefacs]][
        fpdat[flank=="5'" & x==0, ymax, by = usefacs]], # end xstart
    gpdat[xstart==1, .(x = xseg), by = usefacs][
      gpdat[xstart==1, y, by = usefacs]][
        gpdat[xstart==1, ymin, by = usefacs]][
          gpdat[xstart==1, ymax, by = usefacs]]
    )
  rconn<-rbind(gpdat[xstart==max(xstart), .(x = max(xseg)), by = usefacs][
    gpdat[xseg==max(xseg), y, by = usefacs]][
      gpdat[xseg==max(xseg), ymin, by = usefacs]][
        gpdat[xseg==max(xseg), ymax, by = usefacs]], # end x start
    fpdat[flank=="3'" & xlab==0, x, by = usefacs][
      fpdat[flank=="3'" & xlab==0, y, by = usefacs]][
        fpdat[flank=="3'" & xlab==0, ymin, by = usefacs]][
          fpdat[flank=="3'" & xlab==0, ymax, by = usefacs]]
    )
  
  # lconn<-data.table(x = c(fpdat[flank=="5'" & x == 0, x, by = usefacs], gpdat[xstart==1, xseg, by = usefacs]),
  #                   y = c(fpdat[flank=="5'" & x==0, y], gpdat[xstart==1, y]),
  #                   ymin = c(fpdat[flank=="5'" & x==0, ymin],  gpdat[xstart==1, ymin]),
  #                   ymax = c(fpdat[flank=="5'" & x==0, ymax],  gpdat[xstart==1, ymax]))
  # rconn<-data.table(x = c(gpdat[xstart==max(xstart), xseg], fpdat[flank=="3'" & xlab==0, x]),
  #                   y = c(gpdat[xstart==max(xstart), y], fpdat[flank=="3'" & xlab==0, y]),
  #                   ymin = c(gpdat[xstart==max(xstart), ymin], fpdat[flank=="3'" & xlab==0, ymin]),
  #                   ymax = c(gpdat[xstart==max(xstart), ymax], fpdat[flank=="3'" & xlab==0, ymax]))
  # ## If facets, need to do in there tooooooooooo ### can just expand by once I figure that out

  #--- Y axis
  myylim<-c(0, max(gpdat[,max(y) + 0.05*max(y)], fpdat[,max(y) + 0.05*max(y)]))
  if(addcis==T){ # expand for ymax
    myylim<-c(0, max(gpdat[,max(ymax) + 0.05*max(ymax)], fpdat[,max(ymax) + 0.05*max(ymax)]))
  }
  
  # --- X axis labels/breaks
  # Ticks every 10, labels...
  myxlabs<-c("-40", "-20", 1, 73, "+20", "+40")
  mybreaks<-c(-40, -20, 1, 73, fpdat[flank=="3'" & xlab==20, x][1], fpdat[flank=="3'" & xlab==40, x][1])
  # myminorbreaks<-seq(-40, 40, 10)[seq(-40, 40, 10)!=0]
  breaksec<-c(-30, -10, 36, fpdat[flank=="3'" & xlab==5, x][1], fpdat[flank=="3'" & xlab==25, x][1])
  xseclabs<-c("Outer 5'\nflank", "Inner 5'\nflank", "tRNA gene", "Inner 3'\nflank", "Outer 3'\nflank")
  
  # --- Begin plot: FLANK info
  plt<-ggplot()
  ## Add inner flank boundaries
  if(drawinnerflank==T){
    frectdat1<-fpdat[, .(xmin = inner5, xmax = 1, ymin = -Inf, ymax = Inf), by = usefacs]
    frectdat2<-fpdat[,.(xmin = x[flank=="3'" & xlab==0], xmax = x[flank=="3'" & xlab==inner3],
                        ymin = -Inf, ymax = Inf), by = usefacs ]
    
    plt<-plt + geom_rect(data = frectdat1, aes(xmin = xmin, xmax = xmax, ymin = ymin, ymax = ymax), 
                         fill = "lightgray", alpha = 0.5) +
      geom_rect(data = frectdat2, aes(xmin = xmin, xmax = xmax, ymin = ymin, ymax = ymax), 
                fill = "lightgray", alpha = 0.5)
  }
  ## draw flank mutation rate
  plt<-plt+ 
    geom_line(data = fpdat[flank=="5'", ], aes(x, y, color = coldat)) +
    geom_line(data = fpdat[flank=="3'", ], aes(x, y, color = coldat)) +
    scale_x_continuous(expand = c(0, 0), breaks = mybreaks, labels = myxlabs, 
                       sec.axis = sec_axis(~ ., breaks = breaksec, labels = xseclabs)) + 
    scale_y_continuous(expand = c(0,0), limits = myylim) +
    labs(color = labcol) +
    xlab(xlab) + ylab(ylab) + myggtheme + theme(panel.grid = element_blank())
  
  # --- Add gene body info
  plt<-plt + geom_step(data = gpdat, aes(x = xseg, y = y, color = coldat)) +
    geom_line(data = rconn, aes(x, y, color = coldat)) + geom_line(data = lconn, aes(x, y, color = coldat))
  ## Substruct marking if desired. it is SLOW esp w/ facets
  if(drawsubstructs==T){
    plt <- plt + geom_rect_pattern(data = gpdat[coldat==coldat[1]], aes(xmin = xstart - 0.5, xmax = xend + 0.5, ymin = -Inf, ymax = Inf, 
                                                                       fill = structure_plot, pattern_shape = substructure_plot, pattern_color = structure_plot), alpha = 0.4, color = "gray",
                                   pattern = 'pch', pattern_density = 0.2, pattern_spacing = 0.03) + # pattern_frequency, pattern_scale doesn't seem to change anything
      geom_step(data = gpdat, aes(x = xseg, y = y, color = coldat)) + 
      geom_line(data = rconn, aes(x, y, color = coldat)) + geom_line(data = lconn, aes(x, y, color = coldat)) + 
      labs(color = labcol, fill = "tRNA structure component", pattern_shape = "Secondary structure sub-type") + # shape or pattern not changing that label
      guides(pattern_color = "none", fill = guide_legend(override.aes = list(pch = 0), order = 2), color = guide_legend(order = 1))
  }
  
  # --- Color, 95% CIs, etc
  if(!is.null(colvec)){
    plt<- plt + scale_color_manual(values = colvec)
  }
  if(addcis==T){
    mylty <- "solid"
    plt<-plt + 
      geom_step(data= gpdat, aes(x = xseg, y = ymin, color = coldat), linetype = mylty, linewidth = 0.2) +
      geom_step(data = gpdat, aes(x = xseg, y = ymax, color = coldat), linetype = mylty, linewidth = 0.2) +
      geom_line(data = fpdat[flank=="5'", ], aes(x = x, y = ymax, color = coldat), linetype = mylty, linewidth = 0.2) +
      geom_line(data = fpdat[flank=="3'", ], aes(x = x, y = ymax,color = coldat), linetype = mylty, linewidth = 0.2) +
      geom_line(data = fpdat[flank=="5'", ], aes(x = x, y = ymin, color = coldat), linetype = mylty, linewidth = 0.2) +
      geom_line(data = fpdat[flank=="3'", ], aes(x = x, y = ymin,color = coldat), linetype = mylty, linewidth = 0.2) +
      geom_line(data = rconn, aes(x, y = ymin, color = coldat), linetype = mylty, linewidth = 0.2) + 
      geom_line(data = rconn, aes(x, y = ymax, color = coldat), linetype = mylty, linewidth = 0.2) +
      geom_line(data = lconn, aes(x, y = ymin, color = coldat), linetype = mylty, linewidth = 0.2) + 
      geom_line(data = lconn, aes(x, y = ymax, color = coldat), linetype = mylty, linewidth = 0.2)
  }
  
  # --- Return
  return(plt)
}

#### Arguments and inputs ####
p<-arg_parser("Plot/analyze tRNA flank and gene body variation together (after flankvariation.R, genebodyvariation.R)", 
              name = "genebodyvariation.R", hide.opts = TRUE)

# Inputs
p<-add_argument(p, "--flankpertrna",
                help = "Flank variant location, proportion of tRNAs info: *_perpositiontRNAvarianceprop.txt.gz output of flankvariation.R",
                type = "character")
p<-add_argument(p, "--genevarspertrna",
               help = "Within gene body variant info/proportion of tRNAs: *_persecstructtRNAvarianceprop.txt.gz output of genebodyvariation.R")
p<-add_argument(p, "--trnasecstruct",
                help = "Path to file with information on how tRNA sec structures are arranged as used in genebodyvariation.R call. Columns:
                structure (as in --genevars), substructure (as in --genevars), nSubStructure (as in --genevars), canonicalnbp - length of this structure in a 72/73 bp tRNA [for plotting],
                structure_plot and substructure_plot: categories & names prettified/simplified for how you'd like to plot them;
                structure_plot_level & substructure_plot_level: RANKINGS of unique structure_plot and substructure_plot for ordering in plot",
                type = "character")
p<-add_argument(p, "--speciesf",
                help = "File containing information on all species processed in preceeding script. Columns infilename (exactly how all files have this species in their name),
                displayname (name that should be used for plot outputs etc - in other input files), shortname (no-spaces name for ouptut files, sorting, etc - either shorter than or same as infilename, probably).
                In order you'd like plots to be in!",
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

#### Read in data ####
cat("....Reading in data....\n")
fdat<-fread(p$flankpertrna, header = T)
gdat<-fread(p$genevarspertrna, header = T)
tstruct<-fread(p$trnasecstruct, header = T)
sinfo<-fread(p$speciesf, header = T)
# consider doing averages across multiple bp for flanks, too....this is done internal to plot!

#### Make plots ####
cat("....Setting up plotting inf....\n")
# --- Set up plotting-relevant info
# Colors
onecol.all<-c("black")
names(onecol.all)<-"all"
onecol.snv<-c("black")
names(onecol.snv)<-"any.SNV"
onecol.oth<-c("black")
names(onecol.oth)<-"other"

tams<-c("C > T", "G > A") 
# orig had reverse comped but I think that was a strand vs flank assigning issue originally **I think this is right given everything was previously reverse comped I THINK (canonical top strand goes C>T, reverse G>A)

# Mutation categories
mutsplit<-data.table(from = rep(c("A", "C", "G", "T"), each = 4),
                     to = rep(c("A", "C", "G", "T")))
mutsplit<-mutsplit[from!=to, ]
mutsplit[, mutname:=paste(from, to, sep = " > ")] # can further annotate here - TAM categories etc

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
## strand: 
strandinfo<-data.table(descrip = c("Either strand", "+ strand", "- strand"),
                       fdescrip = c("eitherstrand", "plusstrand", "minusstrand"),
                       evaltext = c("Strand%in%c('+', '-')", "Strand=='+'", "Strand=='-'"))

# --- Make plots!
cat("....Making plots....\n")

# consider doing averages across multiple bp for flanks, too....Do FIRST, specifically for individual mutations X > Y

# Species [PDF]
# Strands [diff pages]
# Genes: **  + facet_wrap(~genes, nrow = 3, strip.position = "right") * ## FACTOR THIS FIRST
# Mutation types [diff pages...diff PDFs??]
#     Color ways
#     Only do averaged ones for specific mutation types, probably....but do both for specific mutation types?
#     'all' mutation with AND WITHOUT CIs....

# NOT: sec structure flags; not bp/length normalized in tRNA gene body

invisible(lapply(1:nrow(sinfo), function(spind){ # species
  invisible(lapply(1:nrow(strandinfo), function(stind){ # strand
    pdf(file.path(p$outdir, paste0(p$baseoutname, "_", sinfo[spind, shortname], "_", strandinfo[stind, fdescrip], ".pdf")), 10, 8)
    invisible(lapply(1:nrow(mutinfo), function(mind){ # mutation types
      invisible(lapply(c("all", "MAF < 0.05"), function(maf){ # MAF
        # 'Raw' data, with CIs where appropriate
        print(
        flankgenelocplot(onef = fdat[displayname==sinfo[spind, displayname] & strand==strandinfo[stind, descrip] & varfreq==maf & 
                                       eval(parse(text = mutinfo[mind, seltext]))],
                         oneg = gdat[ displayname==sinfo[spind, displayname] & strand==strandinfo[stind, descrip] & varfreq==maf & 
                                        eval(parse(text = mutinfo[mind, seltext])) & secstructflags=="Any sec. struct. calling"],
                         tstruct = tstruct,
                         colvec = eval(as.name(mutinfo[mind, colvec1])),
                         labcol = "Specific mutation",
                         addcis = mutinfo[mind, addcis],
                         facs = "genes",
                         flankavg = F) + 
          facet_wrap(~factor(genes, levels = c("All tRNA genes", "Active tRNA genes", "Inactive tRNA genes")), nrow = 3, strip.position = "right") +
          ggtitle(paste(c(sinfo[spind, displayname], strandinfo[stind, descrip], paste(maf, "frequencies")), collapse = " | "),
                  subtitle = paste(mutinfo[mind, descrip], "variants, flanks have per-bp rates (non averaged)"))
        )
        
        # 5 bp smoothed data - NEVER CIs
        print(
          flankgenelocplot(onef = fdat[displayname==sinfo[spind, displayname] & strand==strandinfo[stind, descrip] & varfreq==maf & 
                                         eval(parse(text = mutinfo[mind, seltext]))],
                           oneg = gdat[ displayname==sinfo[spind, displayname] & strand==strandinfo[stind, descrip] & varfreq==maf & 
                                          eval(parse(text = mutinfo[mind, seltext])) & secstructflags=="Any sec. struct. calling"],
                           tstruct = tstruct,
                           colvec = eval(as.name(mutinfo[mind, colvec1])),
                           labcol = "Specific mutation",
                           addcis = F,
                           facs = "genes",
                           flankavg = T) + 
            facet_wrap(~factor(genes, levels = c("All tRNA genes", "Active tRNA genes", "Inactive tRNA genes")), nrow = 3, strip.position = "right") +
            ggtitle(paste(c(sinfo[spind, displayname], strandinfo[stind, descrip], paste(maf, "frequencies")), collapse = " | "),
                    subtitle = paste(mutinfo[mind, descrip], "variants, flank rates averaged over 5bp"))
        )
        
        # Second colorway, if there
        if(!is.na(mutinfo[mind, colvec2])){
          # 'Raw' data, with CIs where appropriate
          print(
            flankgenelocplot(onef = fdat[displayname==sinfo[spind, displayname] & strand==strandinfo[stind, descrip] & varfreq==maf & 
                                           eval(parse(text = mutinfo[mind, seltext]))],
                             oneg = gdat[ displayname==sinfo[spind, displayname] & strand==strandinfo[stind, descrip] & varfreq==maf & 
                                            eval(parse(text = mutinfo[mind, seltext])) & secstructflags=="Any sec. struct. calling"],
                             tstruct = tstruct,
                             colvec = eval(as.name(mutinfo[mind, colvec2])),
                             labcol = "Specific mutation\n(second colorway)",
                             addcis = mutinfo[mind, addcis],
                             facs = "genes",
                             flankavg = F) + 
              facet_wrap(~factor(genes, levels = c("All tRNA genes", "Active tRNA genes", "Inactive tRNA genes")), nrow = 3, strip.position = "right") +
              ggtitle(paste(c(sinfo[spind, displayname], strandinfo[stind, descrip], paste(maf, "frequencies")), collapse = " | "),
                      subtitle = paste(mutinfo[mind, descrip], "variants, flanks have per-bp rates (non averaged); second colorway"))
          )
          
          # 5 bp smoothed data - NEVER CIs
          print(
            flankgenelocplot(onef = fdat[displayname==sinfo[spind, displayname] & strand==strandinfo[stind, descrip] & varfreq==maf & 
                                           eval(parse(text = mutinfo[mind, seltext]))],
                             oneg = gdat[ displayname==sinfo[spind, displayname] & strand==strandinfo[stind, descrip] & varfreq==maf & 
                                            eval(parse(text = mutinfo[mind, seltext])) & secstructflags=="Any sec. struct. calling"],
                             tstruct = tstruct,
                             colvec = eval(as.name(mutinfo[mind, colvec2])),
                             labcol = "Specific mutation\n(second colorway)",
                             addcis = F,
                             facs = "genes",
                             flankavg = T) + 
              facet_wrap(~factor(genes, levels = c("All tRNA genes", "Active tRNA genes", "Inactive tRNA genes")), nrow = 3, strip.position = "right") +
              ggtitle(paste(c(sinfo[spind, displayname], strandinfo[stind, descrip], paste(maf, "frequencies")), collapse = " | "),
                      subtitle = paste(mutinfo[mind, descrip], "variants, flank rates averaged over 5bp; second colorway"))
          )
          
        } # end if second colorway
      })) # end mutation types
      })) # end MAF
    invisible(dev.off())
  })) # end strand
})) # end lapply over species
# 
# # tplt[one of all parameters except genes] + facet_wrap(~genes, nrow = 3, strip.position = "right")
# pdf(file.path(p$outdir, "testsize.pdf"), 10, 8) # about right, could go slightly narrower if desired. Probably go a bit taller for title [or narrower]
# print(tplt + facet_wrap(~genes, nrow = 3, strip.position = "right") + ggtitle("I'm a main", sub = "I'm a sub"))
# dev.off()
# 
# # Test w/ moving average
# maplt<-flankgenelocplot(onef = fdat[displayname=="C. elegans" & strand=="Either strand" & varfreq=="all" & alleles=="Major > Minor"],
#                        oneg = gdat[ displayname=="C. elegans" & strand=="Either strand" & varfreq=="all" & alleles=="Major > Minor" & secstructflags=="Any sec. struct. calling"],
#                        tstruct = tstruct,
#                        colvec = mutclasscol,
#                        addcis = F,
#                        facs = "genes",
#                        flankavg = T) + facet_wrap(~factor(genes, levels = c("All tRNA genes", "Active tRNA genes", "Inactive tRNA genes")), nrow = 3, strip.position = "right")

#### Script completion message & session information ####
cat("....flankandgenebodyvar_visualize.R processing complete! Session information:....\n")
sessionInfo()