#! /usr/bin/env/ Rscript
# Process outputs of get_strain_variants.py and build_alt_sequences.py to better account for missingness in VCFs when analyzing tRNAs                                             
# by Avery Davis Bell, begun 2024.11.18

if(length(.libPaths())==1){
  mylibloc <- .libPaths()[1]
}else{ # Presumes on PACE. This is not the best...
  mylibloc <- .libPaths()[grep(R.Version()$platform, .libPaths())]
}
require(data.table, lib.loc = mylibloc) 
require(argparser)
require(ggplot2)
require(ggrepel)

#### Functions ####
ovlapvartrna<-function(vars, trnas){
  # Finds which tRNAs variants fall within (excludes them if they don't)
  # In: vars, data.table with columns Chr and Pos, Ref, Alt (any others will be carried through)
  #     trnas, data.table with columns tRNA, Chr, Start, End (any others will be carried through)
  # Out: data.table containing all entries that overlapped. Columns:
  #       Chr - chromosome of var and tRNA
  # Pos - position of var
  # Ref - reference allele of var
  # Alt - alt allele of var
  # tRNA - tRNA name/ID
  # Start - tRNA start
  # End - tRNA end
  # <any others from tRNA eg Strand, AA, then any others from vars eg homRef, homAlt>

  # Set up data
  ## Variants
  varscop<-copy(vars)
  varscop[, Start:=Pos]
  varscop[, End:=Pos]
  setkey(varscop, Chr, Start, End)
  ## tRNAS
  setkey(trnas, Chr, Start, End)
  
  # Do overlap
  out<-foverlaps(varscop, trnas, type = "any", nomatch = NULL) 

  # Format & return
  setcolorder(out, c("Chr", "Pos", "Ref", "Alt", "tRNA", "Start", "End"))
  out[, `:=`(i.Start = NULL, i.End = NULL)]
  return(out)
}

getbystrain<-function(varswgene, strains, totntrna = 721){
  # Summarizes number of tRNAs strain has missing-called variants in, number of distinct variants this is
  # In: varswgene, information on variants that overlap with tRNAs and which strains have what. Required columns:
  #                 Chr, Pos, Ref, Alt - about variant; tRNA - tRNA ID; missingOrHet - strains that are missing vs not
  #     strains, vector of all strains
  #     totntrna, total # tRNA (to get # without missingness, potentially)
  # Out: data.table with one row per strain. Columns:
  # strain, strain ID
  # nvars_missing, number of variants in tRNAs that had missing (or het) calls in this strain
  # ntrnas_missingvars, number of distinct tRNAs with one or more variants with missing (or het) calls in this strain
  # ntrnas_nomissingvars, number of tRNAs with no missing/het calls - just total number input minus previous number
  
  # Get row per strain obs
  usedat<-copy(varswgene[, .(Chr, Pos, Ref, Alt, tRNA, missingOrHet)])
  setkey(usedat, Chr, Pos, Ref, Alt)
  bystrain<-usedat[, .(tRNA, strain = strsplit(missingOrHet, ",")[[1]]), by = .(Chr, Pos, Ref, Alt)]
  bystrain<-bystrain[strain!="None"]
  
  # Count
  setkey(bystrain, strain)
  ctstrain<-bystrain[, .(nvars_missing = .N, ntrnas_missingvars = length(unique(tRNA))), by = strain]
  ctstrain<-rbind(ctstrain, data.table(strain = strains[!strains%in%ctstrain$strain],
                                       nvars_missing = 0, ntrnas_missingvars = 0)) # Add any strains not included
  setkey(ctstrain, strain)
  ctstrain[, ntrnas_nomissingvars:=totntrna - ntrnas_missingvars]
  
  # Return
  return(ctstrain)
}

#### Plotting related ####
myggtheme<-theme_bw() +
  theme(axis.title.x = element_text(size = 12), axis.text.x = element_text(size = 11), 
        axis.title.y = element_text(size = 12), axis.text.y = element_text(size = 11),
        legend.title = element_text(size = 12), legend.text = element_text(size=11), 
        strip.text.x = element_text(size = 11), strip.text.y = element_text(size = 11), title = element_text(size = 14), 
        strip.text.x.top = element_text(size = 11), strip.text.x.bottom = element_text(size = 11), strip.text.y.right = element_text(size = 11), strip.text.y.left = element_text(size = 11),
        plot.subtitle = element_text(size = 13))

#### Arguments & inputs ####
p<-arg_parser("Process outputs of get_strain_variants.py and build_alt_sequences.py to better account for missingness in VCFs when analyzing tRNAs", 
              name = "catalogmissingness_trnavcfscripts.R", hide.opts = TRUE)

p<-add_argument(p, "--varinfo",
                help = "Path to *_strain_variants.txt output of get_strain_variants.py. **strains are inferred from this",
                type = "character")
p<-add_argument(p, "--trnainfo",
                help = "Path to *_strain_trnas_info.txt output of build_alt_sequences.py",
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

#### Do the work ####
# --- Get in starting data
vars<-fread(p$varinfo, header = T)
setnames(vars, "##Chr", "Chr")

# Get total strain list from variants. Assuming all will be in first 2 rows
strains<-unique(c(unlist(vars[1:3, lapply(.(homRef, homAlt, het), function(x) strsplit(x, ",")[[1]])])))
if("None"%in%strains){ ## rmemove None if it shows up
  strains<-strains[strains!="None"]
}
    
# unique tRNA info
trnas<-unique(fread(p$trnainfo, header = T, select = c("##Chr", "Start", "End", "tRNA", "Strand", "AA", "Codon")))
setnames(trnas, "##Chr", "Chr")

# --- Annotate variants with which tRNAs they are in; nicer format
varswgene<-ovlapvartrna(vars, trnas)
setkey(varswgene, Chr, Pos)
varswgene[, nhomRef:=sum(strsplit(homRef, ",")[[1]]!="None"), by = .(Chr, Pos)]
varswgene[, nhomAlt:=sum(strsplit(homAlt, ",")[[1]]!="None"), by = .(Chr, Pos)]
varswgene[, nMissingOrHet:=sum(strsplit(het, ",")[[1]]!="None"), by = .(Chr, Pos)]
varswgene[, `:=`(homRef=NULL, homAlt=NULL)]
setnames(varswgene, "het", "missingOrHet")
setcolorder(varswgene, c('Chr', 'Pos', 'Ref', 'Alt', 'tRNA', 'Start', 'End', 'Strand', 'AA', 'Codon', 'nhomRef', 'nhomAlt', 'nMissingOrHet', 'missingOrHet'))
# Save!
write.table(varswgene, 
            file.path(p$outdir, paste0(p$baseoutname, "_trna_variant_info_wmissingness_pervariant.txt")),
            sep = "\t", quote = F, row.names = F)

# --- Summaries
# By tRNA
setkey(varswgene, tRNA)
bytrna<-varswgene[, .(Chr[1], Start[1], End[1], Strand[1], AA[1], Codon[1],
                      .N,
                      paste(unique(unlist(strsplit(missingOrHet, ",")))[unique(unlist(strsplit(missingOrHet, ",")))!="None"],
                            collapse = ",")), 
                  by = tRNA]
setnames(bytrna, c('tRNA', 'Chr', 'Start', 'End', 'Strand', 'AA', 'Codon', 'nVCFVars', 'missingOrHet'))
bytrna[, nMissingOrHet:=sum(strsplit(missingOrHet, ",")[[1]]!="None"), by = tRNA]
## Add back any that didn't have variants, too
bytrna<-rbind(bytrna,
              trnas[!tRNA%in%bytrna$tRNA, .(tRNA, Chr, Start, End, Strand, AA, Codon, nVCFVars = 0,
                                            missingOrHet = "", nMissingOrHet = 0)])
setkey(bytrna, tRNA)
setcolorder(bytrna, "missingOrHet", after = ncol(bytrna))
## Save
write.table(bytrna, 
            file.path(p$outdir, paste0(p$baseoutname, "_trna_variant_info_wmissingness_pertrna.txt")),
            sep = "\t", quote = F, row.names = F)

# By strain (# variants with missingness, # tRNAs with variants with missingness)
bystrain<-getbystrain(varswgene, strains, totntrna = nrow(trnas))
write.table(bystrain, 
            file.path(p$outdir, paste0(p$baseoutname, "_trna_variant_info_wmissingness_perstrain.txt")),
            sep = "\t", quote = F, row.names = F)

# --- Couple quick plots [I think I must...]
# By tRNA
bytplt<-ggplot(bytrna, aes(nMissingOrHet)) + geom_histogram(breaks = seq(-0.5, max(bytrna$nMissingOrHet)+0.5, 1)) +
  xlab("Number of strains with missing (or het) variant calls") + ylab("Number of tRNA genes") +
  myggtheme
bytplt.to10<-ggplot(bytrna, aes(nMissingOrHet)) + geom_histogram(breaks = seq(-0.5, max(bytrna$nMissingOrHet)+0.5, 1)) +
  xlab("Number of strains with missing (or het) variant calls") + ylab("Number of tRNA genes\ntruncated at 10") +
  ylim(c(0, 10)) + myggtheme

pdf(file.path(p$outdir, paste0(p$baseoutname, "_trna_variant_info_wmissingness_pertrna.pdf")), 7, 5.5)
print(bytplt)
print(bytplt.to10)
invisible(dev.off())

# by strain
bysplt<-ggplot(bystrain, aes(ntrnas_missingvars, nvars_missing)) + geom_point() +
  geom_label_repel(aes(ntrnas_missingvars, nvars_missing, label = strain)) + # aes(ntrnas_missingvars, nvars_missing, label = strain)
  xlab("Number of tRNAs with at least one variant called missing (KEY!)") + ylab("Number of variants in any tRNA called missing") +
  myggtheme
pdf(file.path(p$outdir, paste0(p$baseoutname, "_trna_variant_info_wmissingness_perstrain.pdf")), 7, 7)
print(bysplt)
invisible(dev.off())

#### Script completion message & session information ####
cat("....catalogmissingness_trnavcfscripts.R processing complete! Session information:....\n")
sessionInfo()