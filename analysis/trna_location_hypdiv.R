#! /usr/bin/env/ Rscript
# Examine tRNA-hyperdivergent region overlap                                    
# by Avery Davis Bell, begun 2026.02.10

if(length(.libPaths())==1){
  mylibloc <- .libPaths()[1]
}else{ # Presumes on PACE. This is not the best...
  mylibloc <- .libPaths()[grep(R.Version()$platform, .libPaths())]
}
require(data.table, lib.loc = mylibloc) 
require(argparser)
require(ggplot2)

#### Functions ####
combovbed<-function(bedf){
  # Combines any overlapping bed file regions
  # Out: bed file with non-overlapping regions as data.table. Columns chr, start, end
  bed<-fread(bedf)[ , .(V1, V2, V3)] # first columns have location, others don't
  setnames(bed, c("chr", "start", "end"))
  
  setkey(bed, chr, start, end)
  setorder(bed, chr, start ,end)
  
  # Do merge
  merged <- bed[, {
    out <- list()
    cs <- start[1]   # current_start
    ce <- end[1]     # current_end
    
    for (i in seq_len(.N)) {
      if (start[i] <= ce) {
        # still overlapping → extend interval
        ce <- max(ce, end[i])
      } else {
        # non-overlap → emit previous, start new
        out <- c(out, list(list(start = cs, end = ce)))
        cs <- start[i]
        ce <- end[i]
      }
    }
    # emit final interval
    out <- c(out, list(list(start = cs, end = ce)))
    rbindlist(out)
  }, by = chr]
  
  return(merged)
}

library(data.table)


tf_chisq2 <- function(dt1, col1, descrip1, dt2, col2, descrip2, na.rm = TRUE, correct = TRUE) {
  #' Chi-square comparison of TRUE/FALSE proportions across two data.tables
  #' originally by copilot, with adjustments
  #'
  #' @param dt1 data.table, first dataset
  #' @param col1 character(1), name of logical column in dt1
  #' @param descrip1 character(1), description of dataset 1 for output
  #' @param dt2 data.table, second dataset
  #' @param col2 character(1), name of logical column in dt2
  #' @param descrip2 character(1), description of dataset 2 for output
  #' @param na.rm logical; if TRUE (default), drop NAs before counting; 
  #'   if FALSE, NAs are treated as FALSE.
  #' @param correct logical; passed to stats::chisq.test (Yates correction for 2×2)
  #' @return one-row data.table with counts, proportions, and chi-square test results
  stopifnot(is.data.table(dt1), is.data.table(dt2))
  stopifnot(is.character(col1), length(col1) == 1L, col1 %in% names(dt1))
  stopifnot(is.character(col2), length(col2) == 1L, col2 %in% names(dt2))
  stopifnot(is.logical(dt1[[col1]]), is.logical(dt2[[col2]]))
  
  prep <- function(v) {
    if (na.rm) {
      v <- v[!is.na(v)]
    } else {
      v[is.na(v)] <- FALSE
    }
    v
  }
  
  v1 <- prep(dt1[[col1]])
  v2 <- prep(dt2[[col2]])
  
  n1_T <- sum(v1, na.rm = FALSE)
  n1   <- length(v1)
  n1_F <- n1 - n1_T
  
  n2_T <- sum(v2, na.rm = FALSE)
  n2   <- length(v2)
  n2_F <- n2 - n2_T
  
  # 2×2 table: rows = datasets (dt1, dt2); cols = (T, F)
  tab <- matrix(c(n1_T, n1_F, n2_T, n2_F),
                nrow = 2, byrow = TRUE,
                dimnames = list(c("dt1", "dt2"), c("T", "F")))
  
  # Chi-square test (Yates correction optional)
  ct <- suppressWarnings(stats::chisq.test(tab, correct = correct))
  
  data.table(
    descrip1 = descrip1,
    descrip2 = descrip2,
    col1 = col1,
    col2 = col2,
    n1 = n1,
    n1_T = n1_T,
    n1_F = n1_F,
    prop1_T = if (n1 > 0) n1_T / n1 else NA_real_,
    n2 = n2,
    n2_T = n2_T,
    n2_F = n2_F,
    prop2_T = if (n2 > 0) n2_T / n2 else NA_real_,
    chisq_stat = unname(ct$statistic),
    df = unname(ct$parameter),
    p_value = unname(ct$p.value),
    method = ct$method,
    correct = isTRUE(correct)
  )
}


#### Arguments & inputs ####
p<-arg_parser("tRNA gene location vs hypdiv regions analyses", 
              name = "trna_location_hypdiv.R", hide.opts = TRUE)

# tRNA Input file related
p<-add_argument(p, "--trnainfo",
                help = "Path to *_tRNAgeneinfo_combined.txt output of trna_location_analyses.R. For species combined",
                type = "character")
p<-add_argument(p, "--species",
                help = "Species that hypdiv info provided here is for (writing at least first version to do one species at a time). As specified in displayname column of tRNA info file.",
                default = "C. elegans")

# Hypdiv input file related
p<-add_argument(p, "--hypdivbed",
                help = "Path to hyperdivergent regions bed file from CaeNDR release matching the date, species used here",
                default = "20220216_c_elegans_divergent_regions_strain.bed")

# Prot coding input file related
p<-add_argument(p, "--pgenelocs",
                help = "Path to no-header, 3 column file specifying all PROTEIN CODING (non tRNA) genes' locations for this species.
                Columns: chr, start, end (NOT header'ed!)",
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

#### Read in & combine files of interest ####
cat("....Reading in and combining files....\n")

# --- Read in
tinfo<-fread(p$trnainfo)[displayname==p$species]
setnames(tinfo, c("Chr", "Start", "End"), c("chr", "start", "end"))

pinfo<-fread(p$pgenelocs)[!V1%in%c("MtDNA", "mtDNA")]
setnames(pinfo, c("chr", "start", "end"))

divs<-combovbed(p$hypdivbed) # combines any present in multiple strains/any overlaps

# --- Annotate genes with if they overlap hypdiv
setkey(divs, chr, start, end)

# tRNAs
setkey(tinfo, chr, start, end)
tinfo[, hypdiv.overlap:=!is.na(foverlaps(tinfo, divs, type = "any", mult = "first", which = T))]

# protein coding genes
setkey(pinfo, chr, start, end)
pinfo[, hypdiv.overlap:=!is.na(foverlaps(pinfo, divs, type = "any", mult = "first", which = T))]

# Save in case want long term
write.table(tinfo, file.path(p$outdir, paste0(p$baseoutname, "_trnainfo_whypdivoverlap.txt")),
            quote = F, row.names = F, sep = "\t")
write.table(tinfo, gzfile(file.path(p$outdir, paste0(p$baseoutname, "_protcgeneinfo_whypdivoverlap.txt.gz"))),
            quote = F, row.names = F, sep = "\t")

#### Analysis compare to distribution of protein coding genes ####
# --- Set up gene sets to do
# (all, nonpseud, pseud) vs protein coding
cgsets<-data.table(descrip = c("All tRNAs (pseud + not)", "Non-pseud tRNAs (not all alleles pseud)",
                               "Pseud tRNAs (all alleles pseud)"),
                   evaltext = c("nAlleles > 0", "nAlleles > Lost", "nAlleles==Lost"))

# --- do all the tests
phypdivchi<-rbindlist(lapply(1:nrow(cgsets), function(i){
  tf_chisq2(dt1 = tinfo[eval(parse(text = cgsets[i, evaltext])),],
            col1 = "hypdiv.overlap",
            descrip1 = cgsets[i, descrip],
            dt2 = pinfo,
            col2 = "hypdiv.overlap",
            descrip2 = "Protein-coding genes")
}))
phypdivchi[, p_bonf_corr:=p_value*nrow(phypdivchi)]
# Save
write.table(phypdivchi, file.path(p$outdir, paste0(p$baseoutname, "_tRNAsVsProtCod_hypdivregions_chisq.txt")),
            quote = F, row.names = F, sep = "\t")

#### Script completion message & session information ####
cat("....trna_location_hypdiv.R processing complete! Session information:....\n")
sessionInfo()

