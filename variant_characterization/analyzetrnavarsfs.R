#! /usr/bin/env/ Rscript
# Various analysis from genotype counts of variants in tRNAs, protein-coding exons          
# by Avery Davis Bell, begun 2025.08.26

if(length(.libPaths())==1){
  mylibloc <- .libPaths()[1]
}else{ # Presumes on PACE. This is not the best...
  mylibloc <- .libPaths()[grep(R.Version()$platform, .libPaths())]
}
require(data.table, lib.loc = mylibloc) 
require(argparser)
require(ggplot2)
require(ggforce) # sina ?

#### plot theme etc ####
myggtheme<-theme_bw() +
  theme(axis.title.x = element_text(size = 12), axis.text.x = element_text(size = 11), 
        axis.title.y = element_text(size = 12), axis.text.y = element_text(size = 11),
        legend.title = element_text(size = 12), legend.text = element_text(size=11), 
        strip.text.x = element_text(size = 11), strip.text.y = element_text(size = 11), title = element_text(size = 14), 
        strip.text.x.top = element_text(size = 11), strip.text.x.bottom = element_text(size = 11), strip.text.y.right = element_text(size = 11), strip.text.y.left = element_text(size = 11),
        plot.subtitle = element_text(size = 13))

#### Functions ####
fread.mult<-function(sinfo, exampfile, header = T){
  # Reads in a file for each species in sinfo
  # In: sinfo, data.table with columns infilename (exactly how all files have this species in their name), displayname, shortname
  #     exampfile, path to file to process but wherever infilename is, should say SAMP instead. Should have header
  # Out: data.table containing whatever's in exampfile for each species in sinfo. Columns added to the front are displayname and shortname
  
  dat<-rbindlist(
    lapply(1:nrow(sinfo), function(s){
      myf<-gsub("SAMP", sinfo[s, infilename], exampfile)
      if(file.exists(myf)){
        dat1<-fread(myf, header = header)
        dat1<-data.table(dat1, sinfo[s, .(displayname, shortname)])
        setcolorder(dat1, c("displayname", "shortname"))
      }else{
        cat(paste("\n**** ***** **** File", myf, "DOES NOT EXIST. Continuing without it. Make sure that's what you want! **** ***** ****\n")) 
      }
    })
  )
  return(dat)
}

findallbedovs<-function(allbed){
  # worker to ID overlapping exons, record & narrow down to non-overlapping ones....
  # In: allbed, data.table with columns "displayname" "shortname"   "chr"         "start"       "end"         "gene_id"     "geneclass"   "bp_length"
  # Out: list of data.tables:
  #     $allbedflagged, input bed but anotated with overlaps_any_row (T if this exon overlaps another) and overlaps_any_gene (T if any exon in this gene overlaps another) columns
  #     $mergedbed, input bed with overlapping exons merged within a gene, so start is now start of first one that overlapped and end is end of last one that overlapped and length reflects this. Has same annotations as above plus n_merged_row: NaN if wasn't merged, number exons merged together if was
  
  # 0) Safety: ensure start <= end in your data (I see rows like 14281439 > 14281436)
  mybed <- copy(allbed)[, id := .I]
  mybed[start > end, c("start","end") := .(end, start)]
  
  # work on flag for if any overlap, possibly DROP big overlaps & keep only one exon length!!
  inclusive<-F
  min_overlap<-1
  
  # 1) Self non-equi overlap join and explicitly carry both sides' start/end
  ov2 <- mybed[
    mybed,
    on = .(displayname, gene_id, chr, start <= end, end >= start),
    allow.cartesian = TRUE, nomatch = 0L,
    
    # j: explicitly retain both sides' ids and coordinates (plus any other fields you need)
    .(
      displayname, gene_id, chr,
      id_x = id, start_x = start, end_x = end,
      id_y = i.id, start_y = i.start, end_y = i.end,
      shortname_x  = shortname,
      shortname_y  = i.shortname,
      geneclass_x  = geneclass,
      geneclass_y  = i.geneclass,
      bp_length_x  = bp_length,
      bp_length_y  = i.bp_length
    )
  ][
    # 2) Exclude self-matches
    id_x != id_y
  ][
    # 3) Compute overlap length
    , overlap := {
      raw <- pmin(end_x, end_y) - pmax(start_x, start_y)
      if (inclusive) raw + 1L else raw
    }
  ][
    # 4) Keep only pairs with overlap > threshold
    overlap > min_overlap
  ]
  
  # 5) De-duplicate symmetric pairs (x,y) == (y,x)
  ov2[, pair := paste(pmin(id_x, id_y), pmax(id_x, id_y))]
  ov2_unique <- unique(ov2, by = "pair")
  
  #
  # Flag rows that participate in any overlap
  overlapping_ids <- unique(c(ov2_unique$id_x, ov2_unique$id_y))
  mybed[, overlaps_any_row := id %in% overlapping_ids]
  mybed[, overlaps_any_gene := (sum(overlaps_any_row) > 0), by = .(displayname, gene_id)]
  
  # --- MERGE overlapping exons for length-ing downstream
  ## find groups of exons to merge
  ovgroups<-lapply(mybed[, unique(displayname)], function(sp){
    ovgroups<-list()
    ov2_un_sp<-ov2_unique[displayname==sp, ]
    # SLOW, SHOULD'VE/COULD'VE DONE THIS WITHIN GENE ID TOO...ADD?
    for(i in 1:nrow(ov2_un_sp)){
      if(i==1){
        ovgroups[[1]]<-unlist(ov2_un_sp[ i, .(id_x, id_y)])
      }
      else{
        myids<-unlist(ov2_un_sp[ i, .(id_x, id_y)])
        flag<-F
        for(j in 1:length(ovgroups)){
          vec<-ovgroups[[j]]
          if(myids[1]%in%vec | myids[2]%in%vec){
            ovgroups[[j]]<-unique(c(ovgroups[[j]], myids))
            flag<-T
          }
        }
        if(flag==F){ # add as own group
          ovgroups[[length(ovgroups) + 1]]<-myids
        }
      }
    }
    return(ovgroups)
  })
  names(ovgroups)<-mybed[, unique(displayname)]
  
  ## Merge all the overlapping ones, save how many records were merged
  outbed.merge<-rbindlist(lapply(names(ovgroups), function(sp){
    onedat<-mybed[displayname==sp]
    out<-rbindlist(lapply(ovgroups[[sp]], function(grp){
      tomerge<-onedat[id%in%grp]
      out1<-tomerge[, .(displayname = unique(displayname), shortname = unique(shortname), chr = unique(chr), 
                        start = min(start), end = max(end), gene_id = unique(gene_id), geneclass = unique(geneclass),
                        overlaps_any_row = T, overlaps_any_gene = T, n_merged_row = .N,
                        id = paste(id, collapse = ","))]
      out1[, bp_length := end - start]
      setcolorder(out1, c(names(onedat), "n_merged_row"))
      return(out1)
    }))
    return(out)
  }))
  ## *** add in ones that DIDN'T get merged too!! n_merged_row = 0
  outbed.nomerge<-mybed[overlaps_any_row==F]
  outbed.nomerge[, n_merged_row := NaN]
  
  ## combine output
  outbed<-rbind(outbed.nomerge, outbed.merge)[, .SD, .SDcols = -c("id")]
  setkeyv(outbed, key(allbed))
  
  return(list(allbedflagged = mybed, 
              mergedbed = outbed))
}

mw.dtout<-function(x, y){
  # Runs wilcox.test, puts results in data.table
  res<-wilcox.test(x, y)
  out<-data.table(med1 = median(x),
                  med2 = median(y),
                  W = res$statistic, 
                  mw.p.value = res$p.value)
  return(out)
}

gtcounts2sfs<-function(gtc.one){
  # Takes a data.table with one row per variant of interest and returns site frequency spectrum (in terms of number singletons, number doubletons, etc)
  mincts<-gtc.one[, min(n_homref, n_homalt), by = .I]$V1
  osfs<-sapply(1:max(mincts), function(x){
    return(sum(mincts==x))
  })
  return(osfs)
}

gettajdsfs <- function(sfs){
  # Modified from https://jyanglab.com/slides/2022-agro932/w6lab.html#2 - Jinliang Yang
  #In:  sfs (site frequency spectrum): number of singletons, doubletons, ..., etc
  # Out: data.table with pi_uncorrected, w_theta_uncorrected, tajimas_d [NATURALLY length corrected since it's ratio]
  n <- length(sfs) + 1 # number of chromosomes      # ??? chromosomes?? mean samples? yeah probably
  ss <- sum(sfs) # number of segregating sites
  a1 <- sum(1 / seq_len(n-1)) 
  a2 <- sum(1 / seq_len(n-1)^2)
  b1 <- (n + 1) / (3 * (n - 1))
  b2 <- 2 * (n^2 + n + 3)/(9 * n * (n - 1))
  c1 <- b1 - 1/a1
  c2 <- b2 - (n + 2)/(a1 * n) + a2 / a1^2
  e1 <- c1 / a1
  e2 <- c2 / (a1^2 + a2)
  Vd <- e1 * ss + e2 * ss * (ss - 1) 
  theta_pi <- sum(2 * seq_len(n-1) * (n - seq_len(n-1)) * sfs)/(n*(n-1))
  theta_w <- ss / a1
  res <- (theta_pi - theta_w) / sqrt(Vd)
  
  return(data.table(n_seg_sites = ss, # number used here
                    pi_uncorrected = theta_pi,
                    w_theta_uncorrected = theta_w,
                    tajimas_d = res))
}

# SO FAR THIS IS NOT PER KB, ALSO NEED LENGTH INFORMATION!!
# may want to do the process per exon and per gene....exons might be more comparable to tRNAs....

#### Arguments & inputs ####
p<-arg_parser("tRNA gene location analyses", 
              name = "trna_location_analyses.R", hide.opts = TRUE)

# tRNA Input file related
p<-add_argument(p, "--speciesf",
                help = "File containing information on all species to process here. Columns infilename (exactly how all files have this species in their name),
                displayname (name that should be used for plot outputs etc), shortname (no-spaces name for ouptut files, sorting, etc - either shorter than or same as infilename, probably).
                In order you'd like plots to be in!. **If not all files exist for each species, only does analyses that it can for each species**",
                type = "character")
p<-add_argument(p, "--trnabed",
                help = "EXAMPLE path to bed files of tRNA gene locations (4th column is name), as used to generate allele counts input.
                Where species ID/species specific info is, put SAMP instead",
                type = "character")
p<-add_argument(p, "--protbed",
                help = "EXAMPLE path to bed files of protein-coding exon locations (4th column is name), as used to generate allele counts input.
                Where species ID/species specific info is, put SAMP instead",
                type = "character")
p<-add_argument(p, "--trnavarcts",
                help = "EXAMPLE path to file with one row per variant in any tRNA with genotype counts. Columns are CHROM, POS, ID, REF, ALT, gene_id, n_homref, n_homalt, n_het, n_miss
                Where species ID/species specific info is, put SAMP instead",
                type = "character")
p<-add_argument(p, "--protvarcts",
                help = "EXAMPLE path to file with one row per variant in any protein-coding exo with genotype counts. Columns are CHROM, POS, ID, REF, ALT, gene_id, n_homref, n_homalt, n_het, n_miss
                Where species ID/species specific info is, put SAMP instead",
                type = "character")


# Output file related
p<-add_argument(p, "--baseoutname",
                help = "Base name for all output files",
                type = "character",
                default = "out")
p<-add_argument(p, "--outdir",
                help = "Outer output directory. Sub-directories will be created internally as needed. **NB: if you provide getwd() here (quote wrapped), current directory will be used",
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
cat("....Reading in data and getting stats of interest....\n")
sinfo<-fread(p$speciesf, header = T)

# --- Read in bed files, get relevant info
tbed<-fread.mult(sinfo, p$trnabed, header = F)[, .(displayname, shortname, V1, V2, V3, V4)]
setnames(tbed, c("V1", "V2", "V3", "V4"), c("chr", "start", "end", "gene_id"))
tbed[, geneclass:="tRNA"]

pbed<-unique(fread.mult(sinfo, p$protbed, header = F))
setnames(pbed, c("V1", "V2", "V3", "V4"), c("chr", "start", "end", "gene_id"))
pbed[, geneclass:="protein"]

allbed<-rbind(tbed, pbed)
rm(tbed, pbed)

allbed[, bp_length:=abs(end - start)]
setkey(allbed, displayname, geneclass, gene_id)

# ** Get with merged exons
usebedinfo<-findallbedovs(allbed)

# get PER GENE length since only tagged vars with genes
gbed<-usebedinfo$mergedbed[, .(shortname = unique(shortname), geneclass = unique(geneclass), chr = unique(chr), start = min(start), end = max(end),
                               bp_length = sum(bp_length), n_exons = .N, overlaps_any_gene = (sum(overlaps_any_gene) > 0),
                               had_merged_exons = !all(is.na(n_merged_row))), by = .(displayname, gene_id)]

# save and document these
datdir<-file.path(p$outdir, "data")
if(!dir.exists(datdir)){dir.create(datdir, recursive = T)}

write.table(usebedinfo$allbedflagged, gzfile(file.path(datdir, paste0(p$baseoutname, "_bedrecords_lengthmergedinfo.txt.gz"))),
            quote = F, row.names = F, sep = "\t")
write.table(usebedinfo$mergedbed, gzfile(file.path(datdir, paste0(p$baseoutname, "_bedrecords_mergedoverlappingexons.txt.gz"))),
            quote = F, row.names = F, sep = "\t")
write.table(gbed, gzfile(file.path(datdir, paste0(p$baseoutname, "_bedrecords_combinedpergene.txt.gz"))),
            quote = F, row.names = F, sep = "\t")

# --- Read in allele count files, annotate with bed info of interest
tc<-data.table(fread.mult(sinfo, p$trnavarcts, header = T), geneclass="tRNA")
pc<-data.table(fread.mult(sinfo, p$protvarcts, header = T), geneclass = "protein")

gtcts<-rbind(tc, pc)
setkey(gtcts, displayname, shortname, geneclass, gene_id)
rm(tc, pc)

# get length of seq in here....blech depends on which exon....
setkey(gbed, displayname, shortname, geneclass, gene_id)
gtcts<-gbed[gtcts]

# --- Get pop gen metrics of interest: per VARIANT
gtcts[, nStrainsWithGT:=n_homref + n_homalt]
gtcts[, pStrainsWithGT:=nStrainsWithGT/(n_homref + n_homalt + n_het + n_miss)] # for filtering later
gtcts[, MAF_obsgts:=min(n_homref, n_homalt)/(n_homref + n_homalt), by = .I]
gtcts[, MAF_allstrains:=min(n_homref, n_homalt)/(n_homref + n_homalt + n_het + n_miss), by = .I]
# Save
write.table(gtcts, gzfile(file.path(datdir, paste0(p$baseoutname, "_pervariant_counts_maf.txt.gz"))),
            quote = F, row.names = F, sep = "\t")

# Summarize site coverage across gene sets, stat difference if any
scov<-gtcts[, .(nsites = .N,nsites95cov = sum(pStrainsWithGT>0.95), nsites90cov = sum(pStrainsWithGT>0.9), nsites80cov = sum(pStrainsWithGT>0.8)), by = .(displayname, geneclass)]
scov[, `:=`(psites95cov = nsites95cov/nsites, psites90cov = nsites90cov/nsites, psites80cov = nsites80cov/nsites)]
write.table(scov, file.path(datdir, paste0(p$baseoutname, "_genotypecoverage_nonmissing_perset.txt")),
            quote = F, row.names = F, sep = "\t")

# --- Get pop gen metrics of interest: per GENE
# .... actually better to POOL together!
# can do this if care about difference among genes, though....

# --- Get pop gen metrics over GENE SETS of interest
#       include or exclude indels
#       missingness is tricky, going to make SFS rare-weighted...Need to FILTER for this so it's not as big a deal
tajdsets<-rbindlist(list(
  data.table(gtypethreshold = "90%", sites = "INDELs included", gtcts[pStrainsWithGT>0.9, gettajdsfs(gtcounts2sfs(.SD)), by = .(displayname, geneclass)]),
  data.table(gtypethreshold = "90%", sites = "INDELs excluded", gtcts[pStrainsWithGT>0.9 & nchar(REF)==1 & nchar(ALT)==1 & nchar(ID) ==1, gettajdsfs(gtcounts2sfs(.SD)), by = .(displayname, geneclass)]), ## sometimes ID has the other allele for....reasons
  data.table(gtypethreshold = "80%", sites = "INDELs included", gtcts[pStrainsWithGT>0.8, gettajdsfs(gtcounts2sfs(.SD)), by = .(displayname, geneclass)]),
  data.table(gtypethreshold = "80%", sites = "INDELs excluded", gtcts[pStrainsWithGT>0.8 & nchar(REF)==1 & nchar(ALT)==1 & nchar(ID) ==1, gettajdsfs(gtcounts2sfs(.SD)), by = .(displayname, geneclass)]) ## sometimes ID has the other allele for....reasons
  
))
write.table(tajdsets, file.path(datdir, paste0(p$baseoutname, "_tajimasd_pergeneset_estimates.txt")),
            quote = F, row.names = F, sep = "\t")

# get numbers PER GENE, per exon
# INCLUDE AND EXCLUDE INDEL VARIANTS FOR TAJIMA'S D....

# get numbers concatenating genes?? [across analysis groups of interest - based on length; tRNA functional; etc - would need to get that set up sooner]
# exclude tRNAs with introns...
# exclude genes with overlapping exons [mult isoforms]....?


#### Analyze/compare across gene groups ####

# initial MAF stats
mafdir<-file.path(p$outdir, "mafrelated")
if(!dir.exists(mafdir)){dir.create(mafdir)}
## Mann-Whitney
mw.maf<-rbindlist(lapply(sinfo$displayname, function(sp){
  out<-rbind(
    data.table(displayname = sp, test = "MAF of all strains, tRNA vs protein",
                  mw.dtout(gtcts[displayname==sp & geneclass=="tRNA", MAF_allstrains], gtcts[displayname==sp & geneclass=="protein", MAF_allstrains])),
    data.table(displayname = sp, test = "MAF of gt observed, tRNA vs protein",
             mw.dtout(gtcts[displayname==sp & geneclass=="tRNA", MAF_obsgts], gtcts[displayname==sp & geneclass=="protein", MAF_obsgts]))
  )
  return(out)
}))
write.table(mw.maf, file.path(mafdir, paste0(p$baseoutname, "_mafMWres.txt")),
            quote = F, row.names = F, sep = "\t")
## Number, proportion that are singletons (observed one strain)
singprop<-gtcts[, .(n_variants =.N,n_singletons = sum(n_homref==1 | n_homalt==1)), by = .(displayname, geneclass)]
singprop[, prop_singletons:=n_singletons/n_variants]
chisq.sing<-singprop[, .(chisq = chisq.test(as.matrix(n_variants, n_singletons, nrow = 2))$statistic,
  df = chisq.test(as.matrix(n_variants, n_singletons, nrow = 2))$parameter,
  chisq.pval = chisq.test(as.matrix(n_variants, n_singletons, nrow = 2))$p.value), by = displayname]
setkey(singprop, displayname)
setkey(chisq.sing, displayname)
singprop.out<-singprop[chisq.sing]
write.table(singprop.out, file.path(mafdir, paste0(p$baseoutname, "_proportionsingletonvariants.txt")),
            quote = F, row.names = F, sep = "\t")

# Initial just plot all the MAFs NO FILTERS
mafp<-ggplot(gtcts) + geom_histogram(aes(MAF_allstrains, fill = geneclass), alpha = 0.4) +
  ylab("Variant count") + xlab("Minor allele frequency (of all strains)") +
  facet_wrap(~factor(displayname, levels = sinfo$displayname)) + myggtheme

bnsz<-0.01
mafp.prop<-ggplot(mapping = aes(MAF_allstrains, fill = geneclass)) +
  geom_histogram(data = gtcts[geneclass=="tRNA"], aes(y = after_stat(density)*bnsz), alpha = 0.6, binwidth = bnsz) +
  geom_histogram(data = gtcts[geneclass=="protein"], aes(y = after_stat(density)*bnsz), alpha = 0.6, binwidth = bnsz) +
  scale_fill_manual(values = c(rgb(0,0,1,0.5), rgb(1,0,0,0.5))) +
  ylab("Proportion of variants in each gene set") + xlab("Minor allele frequency (of all strains)") +
  facet_wrap(~factor(displayname, levels = sinfo$displayname)) + myggtheme

mafp.viol<-ggplot(gtcts) + geom_violin(aes(geneclass, MAF_allstrains, fill = geneclass)) + 
  geom_boxplot(aes(geneclass, MAF_allstrains), outliers = F, fill = NA, color = "darkgray") +
  ylab("Minor allele frequency (of all strains)") + xlab("Gene class") +
  facet_wrap(~factor(displayname, levels = sinfo$displayname)) + myggtheme

mafp.viol.logy<-ggplot(gtcts) + geom_violin(aes(geneclass, MAF_allstrains, fill = geneclass)) + 
  geom_boxplot(aes(geneclass, MAF_allstrains), outliers = F, fill = NA, color = "darkgray") +
  scale_y_log10() +
  ylab("Minor allele frequency (of all strains)") + xlab("Gene class") +
  facet_wrap(~factor(displayname, levels = sinfo$displayname)) + myggtheme

# Initial just plot all the MAFs NO FILTERS - of non-missing vars
mafp.obsgt<-ggplot(gtcts) + geom_histogram(aes(MAF_allstrains, fill = geneclass), alpha = 0.4) +
  ylab("Variant count") + xlab("Minor allele frequency (of strains with gts)") +
  facet_wrap(~factor(displayname, levels = sinfo$displayname)) + myggtheme

bnsz<-0.01
mafp.obsgt.prop<-ggplot(mapping = aes(MAF_obsgts, fill = geneclass)) +
  geom_histogram(data = gtcts[geneclass=="tRNA"], aes(y = after_stat(density)*bnsz), alpha = 0.6, binwidth = bnsz) +
  geom_histogram(data = gtcts[geneclass=="protein"], aes(y = after_stat(density)*bnsz), alpha = 0.6, binwidth = bnsz) +
  scale_fill_manual(values = c(rgb(0,0,1,0.5), rgb(1,0,0,0.5))) +
  ylab("Proportion of variants in each gene set") + xlab("Minor allele frequency (of strains with gts)") +
  facet_wrap(~factor(displayname, levels = sinfo$displayname)) + myggtheme

mafp.obsgt.viol<-ggplot(gtcts) + geom_violin(aes(geneclass, MAF_obsgts, fill = geneclass)) + 
  geom_boxplot(aes(geneclass, MAF_obsgts), outliers = F, fill = NA, color = "darkgray") +
  ylab("Minor allele frequency (of strains with gts)") + xlab("Gene class") +
  facet_wrap(~factor(displayname, levels = sinfo$displayname)) + myggtheme

mafp.obsgt.viol.logy<-ggplot(gtcts) + geom_violin(aes(geneclass, MAF_obsgts, fill = geneclass)) + 
  geom_boxplot(aes(geneclass, MAF_obsgts), outliers = F, fill = NA, color = "darkgray") +
  scale_y_log10() +
  ylab("Minor allele frequency (of strains with gts)") + xlab("Gene class") +
  facet_wrap(~factor(displayname, levels = sinfo$displayname)) + myggtheme


## Save plots
pdf(file.path(mafdir, paste0(p$baseoutname, "_mafdistplots.pdf")), 9, 4)
print(mafp + ggtitle("Raw counts - minor allele freq of all strains"))
print(mafp.prop + ggtitle("Prop. of variants; 0.01 bins"))
print(mafp.viol + ggtitle("MAF of all strains, absolute scale"))
print(mafp.viol.logy + ggtitle("MAF of all strains, log10 y"))

print(mafp.obsgt + ggtitle("Raw counts - minor allele freq of all strains"))
print(mafp.obsgt.prop + ggtitle("Prop. of variants; 0.01 bins, MAF obs gt strains"))
print(mafp.obsgt.viol + ggtitle("MAF of obs gt strains, absolute scale"))
print(mafp.obsgt.viol.logy + ggtitle("MAF of obs gt strains, log10 y"))
invisible(dev.off())


# Plots & stats
# FILTER to non-pseud tRNAs or something here? INTRONS EXCLUDE

#### Script completion message & session information ####
cat("....analyzetrnavarsfs.R processing complete! Session information:....\n")
sessionInfo()
