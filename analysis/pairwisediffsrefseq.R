#! /usr/bin/env/ Rscript
# Work on getting pairwise # differences between tRNA sequences (ref genome, presumably)
# by Avery Davis Bell, begun 2025.07.14

# Load packages
require(data.table)
require(argparser)
require(ggplot2)

#### Functions ####
myggtheme<-theme_bw() +
  theme(axis.title.x = element_text(size = 12), axis.text.x = element_text(size = 11), 
        axis.title.y = element_text(size = 12), axis.text.y = element_text(size = 11),
        legend.title = element_text(size = 12), legend.text = element_text(size=11), 
        strip.text.x = element_text(size = 11), strip.text.y = element_text(size = 11), title = element_text(size = 14), 
        strip.text.x.top = element_text(size = 11), strip.text.x.bottom = element_text(size = 11), strip.text.y.right = element_text(size = 11), strip.text.y.left = element_text(size = 11),
        plot.subtitle = element_text(size = 13))

aln2seq<-function(char){
  # little worker
  seqsp<-strsplit(char, "")[[1]]
  out<-paste(tolower(seqsp[seqsp!="-"]), collapse = "")
  return(out)
}

readxfasta<-function(xfastaf){
  # Reads an X fasta into a data.table with columns name, sequence, secstruct
  
  mylines<-readLines(xfastaf)
  out<-data.table(name = sapply(strsplit(mylines[seq(1, length(mylines), 3)], ">"), function(x) x[2]),
                  sequence = mylines[seq(2, length(mylines), 3)],
                  secstruct = mylines[seq(3, length(mylines), 3)])
  return(out)
}

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

alnseqdist<-function(forpairs){
  # Gets the number of sequence elements (bases, gaps) that do not match & that do match for all pairwise comparisons of input seqs
  # In: forpairs, data.table to get distance between 'sequence' column of. Names displayname, tRNA, seqlabel, Intron [NA if no intron, where it is if intron] are used
  # Out: data.table with one row per pairwise comparison done. Columns:
  #         genepair, tRNA genes used here, '-' separated 
  #         seqlabelpair, seq label of the genes compared here (so can match with other genes), '-' separated
  #         intron: count of genes in pair with predicted intron
  #         ntNoMatch, number of elements in alignment [nts] that were NOT identical
  #         ntMatch, number of elements in alignment [nts] that were identical
  #         stretchNoMatch, number of distinct STRETCHES of alignments that don't match - i.e., if elements are F F F T, this is 1 stretch non-match, 3 nts non match
  #         stretchMatch, number of distinct STRETCHES of alignments that match...not sure this will be used, but might as well have
  
  out<-rbindlist(lapply(1:(nrow(forpairs) - 1), function(i){
    rbindlist(lapply((i+1):nrow(forpairs), function(j){
      compvec<-(strsplit(forpairs[i, sequence], "")[[1]]!=strsplit(forpairs[j, sequence], "")[[1]])
      comprle<-rle(compvec)
      out<-data.table(genepair = paste(sort(unique(c(forpairs[i, tRNA], forpairs[j, tRNA]))), collapse = "-"),
                      seqlabelpair = paste(sort(unique(c(forpairs[i, seqlabel], forpairs[j, seqlabel]))), collapse = "-"),
                      intron = forpairs[i, sum(!is.na(Intron))] + forpairs[j, sum(!is.na(Intron))],
                      nNoMatch = sum(compvec),
                      nMatch = sum(!compvec),
                      stretchNoMatch = sum(comprle$values),
                      stretchMatch = sum(!comprle$values))
    }))
  }))
  
  return(out)
}

#### Arguments & inputs ####
p<-arg_parser("Work on getting pairwise # differences between tRNA sequences (ref genome, presumably)", 
              name = "pairwisediffsrefseq.R", hide.opts = TRUE)

# Inputs
p<-add_argument(p, "--speciesf",
                help = "File containing information on all species processed in preceeding script(s). Columns infilename (exactly how all files have this species in their name),
                displayname (name that should be used for plot outputs etc - in other input files), shortname (no-spaces name for ouptut files, sorting, etc - either shorter than or same as infilename, probably).
                In order you'd like plots to be in!",
                type = "character")
p<-add_argument(p, "--xfasta",
                help = "Foursale alignment xfasta for ALL tRNA alleles in ALL species of interest",
                default = "~/GaTech Dropbox/Avery Bell/PaabyLab/Projects/tRNA/tRNAVarsFromVCFs/DataFromPACE/tRNAAlignmentsSecStruct/fullalignments/allalleles_3species_202312and202401caendr_foursale_aligned.xfasta.gz")
p<-add_argument(p, "--introninfo",
                help = "EXAMPLE path to *_tRNA2struct_info.txt.gz output file of secstruct2seqpieces.py, used here to ID genes with suspected introns vs not.
                        where sample info goes, put SAMP instead (as in speciesf)",
                default = "~/Dropbox-GaTech//Avery Bell/PaabyLab/Projects/tRNA/tRNAVarsFromVCFs/202506CaeNDRAnalyses/DataFromPACE/tRNAAlignmentsSecStruct/allelesecstruct/SAMP/SAMP_tRNA2struct_info.txt.gz")
p<-add_argument(p, "--seqgroupinfo",
                help = "*_identicalseqgroupinfo.txt output of trna_location_analyses.R: information on all seqeunces from ref genome - columns:
                					displayname, species
					                tRNA, tRNAscan-SE tRNA sequence ID, with 'chr' stripped if it was present
					                seq, sequence of that tRNA
					                seqlabel, unique sequence number - same sequences have the same number
					                nwseqlabel, number tRNAs with this same sequence in ref genome",
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


#### Get in data: ref genome aligned seqs ####
cat("....Reading in data....\n")
sinfo<-fread(p$speciesf, header = T)
sgrps<-fread(p$seqgroupinfo, header = T)

# --- intron info
tinfo<-fread.mult(sinfo, p$introninfo)
# narrow to reference seq
tinfo<-tinfo[grepl("Reference", tRNA)]
# Fix name format
tinfo<-data.table(tinfo[, tstrsplit(tRNA, "_")], tinfo)
setnames(tinfo, c("V1", "tRNA"), c("tRNA", "longnameinfo"))
tinfo[substr(tRNA, 1, 3)=="chr", tRNA:=substr(tRNA, 4, nchar(tRNA))] 
tinfo.use<-tinfo[, .(displayname, tRNA, Intron)]
setkey(tinfo.use, displayname, tRNA)

# --- Aligned fasta with gaps
xf<-readxfasta(xfastaf =p$xfasta)[, .(name, sequence)]
# Narrow to just reference seqs
xf<-xf[grepl("Reference", name)]

# Fix name format & strip out species
sdat<-data.table(xf[, tstrsplit(name, "_")], xf)
setnames(sdat, c("shortname", "tRNA", "strain", "alextrainfo", "name", "sequence"))
sdat[substr(tRNA, 1, 3)=="chr", tRNA:=substr(tRNA, 4, nchar(tRNA))]
setkey(sinfo, shortname)
setkey(sdat, shortname)
sdat<-sinfo[sdat] # add full seq names
sdat[, `:=`(infilename=NULL, shortname=NULL, strain=NULL, alextrainfo=NULL)] # drop extraneous columns

# Add in intron info
setkey(sdat, displayname, tRNA)
sdat<-tinfo.use[sdat]

# Narrow to only UNIQUE sequences, annotating with seq group type info
## Add in seq grp info
setkey(sgrps, displayname, tRNA)
setkey(sdat, displayname, tRNA)
sdat<-sgrps[sdat]
## keep only first occurrence of each seqlabel [w/in species]
unqseqs<-sdat[0]
for(sp in sinfo$displayname){
  onedat<-sdat[displayname==sp, ]
  spout<-sdat[0]
  for(i in 1:nrow(onedat)){
    if(!onedat[i, seqlabel]%in%spout$seqlabel){
      spout<-rbind(spout, onedat[i, ])
    }
  }
  unqseqs<-rbind(unqseqs, spout)
}

#### Determine nt distance between seq ALIGNMENTS ####
cat("....Determining nt distance between seq alignments....\n")
sdists<-unqseqs[, alnseqdist(.SD), by = displayname]

write.table(sdists, gzfile(file.path(p$outdir, paste0(p$baseoutname, "_alignedseqdiffs_uniqueseqpairs.txt.gz"))),
            sep = "\t", quote = F, row.names = F)

#### Summary plots ####
# --- All genes
plt.nt<-ggplot(sdists) + geom_histogram(aes(nNoMatch)) +
  xlab("Number nucleotides/gaps in alignment that don't match") + ylab("Number of unique tRNA sequence pairs") +
  ggtitle("Number nts/gaps that don't match", subtitle = "All unique seq gene pairs") +
  myggtheme + facet_wrap(~displayname, nrow = 3)
plt.pnt<-ggplot(sdists) + geom_histogram(aes(nNoMatch/(nNoMatch + nMatch))) +
  xlab("Prop. nucleotides/gaps in alignment that don't match") + ylab("Prop. of unique tRNA sequence pairs") +
  ggtitle("Prop. nts/gaps that don't match", subtitle = "All unique seq gene pairs") +
  myggtheme + facet_wrap(~displayname, nrow = 3)
plt.ns<-ggplot(sdists) + geom_histogram(aes(stretchNoMatch)) +
  xlab("Number stretches of nucleotides/gaps in alignment that don't match") + ylab("Number of unique tRNA sequence pairs") +
  ggtitle("Number stretches of nts/gaps that don't match", subtitle = "All unique seq gene pairs") +
  myggtheme + facet_wrap(~displayname, nrow = 3)

# --- Pairs without any introns
plt.nointrons<-list(
  p1 = ggplot(sdists[intron==0]) + geom_histogram(aes(nNoMatch)) +
  xlab("Number nucleotides/gaps in alignment that don't match") + ylab("Number of unique tRNA sequence pairs") +
  ggtitle("Number nts/gaps that don't match", subtitle = "Unique seq gene pairs with no introns") +
  myggtheme + facet_wrap(~displayname, nrow = 3),
  
  p2 = ggplot(sdists[intron==0]) + geom_histogram(aes(nNoMatch/(nNoMatch + nMatch))) +
  xlab("Prop. nucleotides/gaps in alignment that don't match") + ylab("Prop. of unique tRNA sequence pairs") +
  ggtitle("Prop. nts/gaps that don't match", subtitle = "Unique seq gene pairs with no introns") +
  myggtheme + facet_wrap(~displayname, nrow = 3),
  
p3 =ggplot(sdists[intron==0]) + geom_histogram(aes(stretchNoMatch)) +
  xlab("Number stretches of nucleotides/gaps in alignment that don't match") + ylab("Number of unique tRNA sequence pairs") +
  ggtitle("Number stretches of nts/gaps that don't match", subtitle = "Unique seq gene pairs with no introns") +
  myggtheme + facet_wrap(~displayname, nrow = 3)
)

# --- Pairs with any introns
plt.wintrons<-list(
  p1 = ggplot(sdists[intron>0]) + geom_histogram(aes(nNoMatch)) +
    xlab("Number nucleotides/gaps in alignment that don't match") + ylab("Number of unique tRNA sequence pairs") +
    ggtitle("Number nts/gaps that don't match", subtitle = "Unique seq gene pairs WITH intron(s)") +
    myggtheme + facet_wrap(~displayname, nrow = 3),
  
  p2 = ggplot(sdists[intron>0]) + geom_histogram(aes(nNoMatch/(nNoMatch + nMatch))) +
    xlab("Prop. nucleotides/gaps in alignment that don't match") + ylab("Prop. of unique tRNA sequence pairs") +
    ggtitle("Prop. nts/gaps that don't match", subtitle = "Unique seq gene pairs WITH intron(s)") +
    myggtheme + facet_wrap(~displayname, nrow = 3),
  
  p3 =ggplot(sdists[intron>0]) + geom_histogram(aes(stretchNoMatch)) +
    xlab("Number stretches of nucleotides/gaps in alignment that don't match") + ylab("Number of unique tRNA sequence pairs") +
    ggtitle("Number stretches of nts/gaps that don't match", subtitle = "Unique seq gene pairs WITH intron(s)") +
    myggtheme + facet_wrap(~displayname, nrow = 3)
)

pdf(file.path(p$outdir, paste0(p$baseoutname, "_alignedseqdiffs_uniqueseqpairs_hists.pdf")), 8, 10)
# all genes
print(plt.nt)
print(plt.pnt)
print(plt.ns)
# no introns
print(plt.nointrons)
# any introns
print(plt.wintrons)
invisible(dev.off())

#### Script completion message & session information ####
cat("....pairwisediffsrefseq.R processing complete! Session information:....\n")
sessionInfo()
