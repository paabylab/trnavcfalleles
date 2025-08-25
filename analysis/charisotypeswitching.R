#! /usr/bin/env/ Rscript
# More focused investigation of isotype switching
# by Avery Davis Bell, begun 2025.05.13

# Load packages
require(data.table)
require(argparser)
require(ggplot2)
require(ggh4x)
require(formattable)
require(ggpattern)
require(cowplot)
require(circlize)

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

#### Functions: data ####
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

summbyaa<-function(alinfo, aas, trips){
  # Gets counts of isotype switches for each possible isotype switching - for comparing across amino acids, for example
  # In:
  # Out: Big list of data tables:
  #       $aafromto, Summary of n alleles etc for each amino acid from-to (allelecm-aa) combination. Columns (describing here, they're similar though not identical for all):
  #             #         displayname, species - from input
  #                 fromTo, from.val-to.val in case want already concatenated
  #                 from, from value (from input)
  #                 to, to value (from input)
  #                 nAlleles, # alleles that are from-to
  #                 nGenesThisSw, # genes with at least one from-to allele
  #                 nGenesThisSw.invariant, # genes with at least one from-to allele with only one allele (this one) in pop
  #                 nGenesThisSw.variant, # genes with at least one from-to allele with multiple alleles in pop
  #                 nGenesThisSw.allThisSw, # gnees with at least one from-to allele and all alleles in population are THIS isotype switch
  #                 nGenesThisSw.notallThisSw, # genes with at least one from-to allele and all alleles in pop are NOT THIS isotype sw
  #                 nAllelesThisFrom, # alleles in population that have this from value (to ignored) - like all the ones with leucine allelecm, for example
  #                 nAllelesThisFrom.notAllPseud, # alleles in population that have this from value (to ignored) - like all the ones with leucine allelecm, for example - excluding ones with all pseud alleles
  #                 nGenesThisFrom, # genes in population that have any alleles with this from value (to ignored) - like all the ones with leucine allelecm, for example 
  #                 nGenesThisFrom.notAllPseud, # genes in population that have any alleles with this from value (to ignored) - like all the ones with leucine allelecm, for example - excluding ones with all pseud alleles
  #                 nAllelesThisTo, # as with from, but for ones with the provided To value (e.g. leucine amino acid)
  #                 nAllelesThisTo.notAllPseud, # as with from, but for ones with the provided To value (e.g. leucine amino acid)
  #                 nGenesThisTo, # as with from, but for ones with the provided To value (e.g. leucine amino acid)
  #                 nGenesThisTo.notAllPseud # as with from, but for ones with the provided To value (e.g. leucine amino acid)
  #     $aafrom, summary of n alleles etc that mutate FROM the same backbone, per amino acid. First columns are counts for when they switch isotypes, last columns are overall counts of genes with this alleleCM (as in earlier output)
  #     $aato, as above, but for those that mutate TO the same called AA
  #     $tripfromto, as aafromto but for every pair of triplets/codons instead of amino acids. Amino acids also specified with columns. Includes codon pairs that wouldn't result in the amino acid changing.
  #     $trfrom.all, as in aafrom, but for each triplet/codon. Includes if switches to a codon that doesn't change amino acid.
  #     $trfrom.aasw, as trfrom.all but first columns ONLY count those that cause AA switch. (latter columns count all switched, invariant, not , etc)
  #     $trto.all, as in trpfrom.all but mutating TO specified triplet
  #     $trto.aasw, as in trpfrom.aasw but mutating TO specified triplet (for those that cause AA switch.)
  
  # --- Subfunction (do ONE)
  onesummby<-function(alinfo, from.col = "AlleleCM", to.col = "AA", from.val = "Glu", to.val = "Ala"){
    # Gets number of alleles, genes with specified switch-ness, as well as background/denominator numbers
    # In: alinfo, data.table of all alleles. Columns must include displayname, from.col and to.col values, switch, VariableInPop, all.altered, all.lost
    #     from.col, character name of column you want to count switching FROM
    #     to.col, character name of column you want to count switching TO
    #     from.val, value to look for in from.col
    #     to.val, value to look for in to.col
    # Out:  data.table with one row per displayname characterizing this switch. Columns: (doesn't have pseud stuff for w/ this sw bc can't be all pseud w/ sw)
    #         displayname, species - from input
    #         fromTo, from.val-to.val in case want already concatenated
    #         from, from value (from input)
    #         to, to value (from input)
    #         nAlleles, # alleles that are from-to
    #         nGenesThisSw, # genes with at least one from-to allele
    #         nGenesThisSw.invariant, # genes with at least one from-to allele with only one allele (this one) in pop
    #         nGenesThisSw.variant, # genes with at least one from-to allele with multiple alleles in pop
    #         nGenesThisSw.allThisSw, # gnees with at least one from-to allele and all alleles in population are THIS isotype switch
    #         nGenesThisSw.allThisSw.variant, # genes with ALL same from-to alleles AND multiple alleles in pop
    #         nGenesThisSw.notallThisSw, # genes with at least one from-to allele and all alleles in pop are NOT isotype switched
    #         nAllelesThisFrom, # alleles in population that have this from value (to ignored) - like all the ones with leucine allelecm, for example
    #         nAllelesThisFrom.notAllPseud, # alleles in population that have this from value (to ignored) - like all the ones with leucine allelecm, for example - excluding ones with all pseud alleles
    #         nGenesThisFrom, # genes in population that have any alleles with this from value (to ignored) - like all the ones with leucine allelecm, for example 
    #         nGenesThisFrom.notAllPseud, # genes in population that have any alleles with this from value (to ignored) - like all the ones with leucine allelecm, for example - excluding ones with all pseud alleles
    #         nGenesThisFrom.noSwNoPseud.invariant, # genes with this from value that have NO isotype switches and NO allelic variation (and all non-pseud alleles)
    #         nGenesThisFrom.noSwNoPseud.variant, # genes with this from value that have NO isotype switches and have allelicvariation (and all non-pseud alleles)
    #         nAllelesThisTo, # as with from, but for ones with the provided To value (e.g. leucine amino acid)
    #         nAllelesThisTo.notAllPseud, # as with from, but for ones with the provided To value (e.g. leucine amino acid)
    #         nGenesThisTo, # as with from, but for ones with the provided To value (e.g. leucine amino acid)
    #         nGenesThisTo.notAllPseud # as with from, but for ones with the provided To value (e.g. leucine amino acid)
    #         nGenesThisTo.noSwNoPseud.invariant, # genes with this from value that have NO isotype switches and NO allelic variation (and all non-pseud alleles)
    #         nGenesThisTo.noSwNoPseud.variant, # genes with this to value that have NO isotype switches and have allelicvariation (and all non-pseud alleles)
    
    # Data set up
    udat<-copy(alinfo) # just for cleanliness
    udat[, `:=`(from.col=get(from.col), to.col = get(to.col))]
    setkey(udat, displayname)
    
    # Get numbers
    out<-udat[, .(from = from.val, to = to.val,
                  nAlleles = sum(switch==T & from.col==from.val & to.col==to.val, na.rm = T),
                  nGenesThisSw = length(unique(tRNA[switch==T & from.col==from.val & to.col==to.val])),
                  nGenesThisSw.invariant = length(unique(tRNA[switch==T & from.col==from.val & to.col==to.val & VariableInPop==F])), 
                  nGenesThisSw.variant = length(unique(tRNA[switch==T & from.col==from.val & to.col==to.val & VariableInPop==T])), # could be all altered or not all altered
                  nGenesThisSw.allThisSw = length(unique(tRNA[switch==T & from.col==from.val & to.col==to.val & all.altered==T])),
                  nGenesThisSw.allThisSw.variant = length(unique(tRNA[switch==T & from.col==from.val & to.col==to.val & all.altered==T & VariableInPop==T])),
                  nGenesThisSw.notallThisSw = length(unique(tRNA[switch==T & from.col==from.val & to.col==to.val & all.altered==F & VariableInPop==T])),
                  nAllelesThisFrom = sum(from.col==from.val, na.rm = T),
                  nAllelesThisFrom.notAllPseud = sum(from.col==from.val & all.lost==F, na.rm = T), # alleles from genes that are not all pseud
                  nGenesThisFrom = length(unique(tRNA[from.col==from.val])), # n genes with this amino acid backbone, for example (plus subsets thereof)
                  nGenesThisFrom.notAllPseud = length(unique(tRNA[from.col==from.val & all.lost==F])),
                  nGenesThisFrom.noSwNoPseud.invariant = length(unique(tRNA[from.col==from.val & all.lost==F & any.switch==F & VariableInPop==F])), ####
                  nGenesThisFrom.noSwNoPseud.variant = length(unique(tRNA[from.col==from.val & all.lost==F & any.switch==F & VariableInPop==T])), ####
                  nAllelesThisTo = sum(to.col==to.val, na.rm = T),
                  nAllelesThisTo.notAllPseud = sum(to.col==to.val & all.lost==F, na.rm = T), # alleles from genes that are not all pseud
                  nGenesThisTo = length(unique(tRNA[to.col==to.val])), # n genes with this amino acid backbone, for example (plus subsets thereof)
                  nGenesThisTo.notAllPseud = length(unique(tRNA[to.col==to.val & all.lost==F])),
                  nGenesThisTo.noSwNoPseud.invariant = length(unique(tRNA[to.col==to.val & all.lost==F & any.switch==F & VariableInPop==F])), ####
                  nGenesThisTo.noSwNoPseud.variant = length(unique(tRNA[to.col==to.val & all.lost==F & any.switch==F & VariableInPop==T]))),
              by = displayname]
    
    out[, fromTo:=paste(from.val, to.val, sep = "-")]
    setcolorder(out, c("displayname", "fromTo"))
    
    # Return
    return(out)
  }
  
  # --- By amino acid
  # Set up all possible from-to combinations (for actual aas observed in these populations)
  aado<-data.table(aa.from = rep(aas, each = length(aas)),
                   aa.to = rep(aas, length(aas)))
  aado<-aado[aa.from!=aa.to] # not same AA 
  # Do counts for each of these that are NOT same AA 
  nsumm.aafromto<-rbindlist(lapply(1:nrow(aado), function(x){
    onesummby(alinfo, from.col = "AlleleCM", to.col = "AA", 
              from.val = aado[x, aa.from], to.val = aado[x, aa.to])
  }))
  
  
  # --- By triplet
  # Set up all possible from-to combinations
  trdo<-data.table(trip.from = rep(trips$Codon, each = length(trips$Codon)),
                   trip.to = rep(trips$Codon, length(trips$Codon)),
                   aa.from = rep(trips$AA, each = length(trips$AA)),
                   aa.to = rep(trips$AA, length(trips$AA)))
  trdo<-unique(trdo[trip.from!=trip.to, .(trip.from, trip.to)]) # not same AA or codon...NO!! including those for future analysis...
  
  # Do counts of each of these. INCLUDES codon switches that don't change amino acids for later use....
  nsumm.trfromto<-rbindlist(lapply(1:nrow(trdo), function(x){
    out<-onesummby(alinfo, from.col = "Codon.orig", to.col = "Codon",
              from.val = trdo[x, trip.from], to.val = trdo[x, trip.to])
    out[, `:=`(AAFromCod = trips[Codon==trdo[x, trip.from], paste(AA, collapse = ",")], 
               AAToCod = trips[Codon==trdo[x, trip.to], paste(AA, collapse = ",")])] # includes multiple where appropriate (iMet etc)
    setcolorder(out, c("displayname", "fromTo", "from", "to", "AAFromCod", "AAToCod"))
    return(out)
  }))
 
  # --- Add up for all w switches *from* certain value (and *to* separately). the TOTALS (from-to ending columns) shouldn't be added but rather copied 
  ## AA level
  setkey(nsumm.aafromto, displayname, from)
  nsumm.aafrom<-nsumm.aafromto[,.(nAlleles = sum(nAlleles), 
                                  nGenesSwFromThis = sum(nGenesThisSw),
                                  nGenesSwFromThis.invariant = sum(nGenesThisSw.invariant),
                                  nGenesSwFromThis.variant = sum(nGenesThisSw.variant),
                                  nGenesFromThis.AllSameSw = sum(nGenesThisSw.allThisSw),
                                  nGenesFromThis.allThisSw.variant = sum(nGenesThisSw.allThisSw.variant),
                                  nGenesFromThis.notallThisSw = sum(nGenesThisSw.notallThisSw),
                                  nAllelesThisFrom = nAllelesThisFrom[1],  # this is ANY with AlleleCM this value, not just if they have switch
                                  nAllelesThisFrom.notAllPseud = nAllelesThisFrom.notAllPseud[1],
                                  nGenesThisFrom = nGenesThisFrom[1], 
                                  nGenesThisFrom.notAllPseud = nGenesThisFrom.notAllPseud[1],
                                  nGenesThisFrom.noSwNoPseud.invariant = nGenesThisFrom.noSwNoPseud.invariant[1],
                                  nGenesThisFrom.noSwNoPseud.variant = nGenesThisFrom.noSwNoPseud.variant[1]) ,
                               by = .(displayname, from)]
  setkey(nsumm.aafromto, displayname, to)
  nsumm.aato<-nsumm.aafromto[, .(nAlleles = sum(nAlleles), 
                                 nGenesSwToThis = sum(nGenesThisSw),
                                 nGenesSwToThis.invariant = sum(nGenesThisSw.invariant),
                                 nGenesSwToThis.variant = sum(nGenesThisSw.variant),
                                 nGenesToThis.AllSameSw = sum(nGenesThisSw.allThisSw),
                                 nGenesToThis.allThisSw.variant = sum(nGenesThisSw.allThisSw.variant),
                                 nGenesToThis.notallThisSw = sum(nGenesThisSw.notallThisSw),
                                 nAllelesThisTo = nAllelesThisTo[1],
                                 nAllelesThisTo.notAllPseud = nAllelesThisTo.notAllPseud[1],
                                 nGenesThisTo = nGenesThisTo[1],
                                 nGenesThisTo.notAllPseud = nGenesThisTo.notAllPseud[1],
                                 nGenesThisTo.noSwNoPseud.invariant = nGenesThisTo.noSwNoPseud.invariant[1],
                                 nGenesThisTo.noSwNoPseud.variant = nGenesThisTo.noSwNoPseud.variant[1]),
                             by = .(displayname, to)]
  
  # do at triplet/codon level, keep in mind have included codon switches that don't change amino acid here for later use...
  #     Show, for ex, if there's like hyper mutable codon
  setkey(nsumm.trfromto, displayname, from)
  ## ANY mutation
  trfrom.all<-nsumm.trfromto[, .(AAFromCod = AAFromCod[1],
                                   nAlleles = sum(nAlleles), 
                                   nGenesSwFromThis = sum(nGenesThisSw),
                                   nGenesSwFromThis.invariant = sum(nGenesThisSw.invariant),
                                   nGenesSwFromThis.variant = sum(nGenesThisSw.variant),
                                   nGenesFromThis.AllSameSw = sum(nGenesThisSw.allThisSw),
                                   nGenesFromThis.allThisSw.variant = sum(nGenesThisSw.allThisSw.variant),
                                   nGenesFromThis.notallThisSw = sum(nGenesThisSw.notallThisSw),
                                   nAllelesThisFrom = nAllelesThisFrom[1],  # this is ANY with AlleleCM this value, not just if they have switch
                                   nAllelesThisFrom.notAllPseud = nAllelesThisFrom.notAllPseud[1],
                                   nGenesThisFrom = nGenesThisFrom[1], 
                                   nGenesThisFrom.notAllPseud = nGenesThisFrom.notAllPseud[1],
                                 nGenesThisFrom.noSwNoPseud.invariant = nGenesThisFrom.noSwNoPseud.invariant[1],
                                 nGenesThisFrom.noSwNoPseud.variant = nGenesThisFrom.noSwNoPseud.variant[1]),
                               by = .(displayname, from)]
  ## Mutation that causes amino acid switch, still get totals
  trfrom.aasw<-nsumm.trfromto[, .(AAFromCod = AAFromCod[1],
                                  nAlleles = sum(nAlleles[AAFromCod!=AAToCod]), 
                                  nGenesSwFromThis = sum(nGenesThisSw[AAFromCod!=AAToCod]),
                                  nGenesSwFromThis.invariant = sum(nGenesThisSw.invariant[AAFromCod!=AAToCod]),
                                  nGenesSwFromThis.variant = sum(nGenesThisSw.variant[AAFromCod!=AAToCod]),
                                  nGenesFromThis.AllSameSw = sum(nGenesThisSw.allThisSw[AAFromCod!=AAToCod]),
                                  nGenesFromThis.allThisSw.variant = sum(nGenesThisSw.allThisSw.variant[AAFromCod!=AAToCod]),
                                  nGenesFromThis.notallThisSw = sum(nGenesThisSw.notallThisSw[AAFromCod!=AAToCod]),
                                  nAllelesThisFrom = nAllelesThisFrom[1],  # this is ANY with AlleleCM this value, not just if they have switch
                                  nAllelesThisFrom.notAllPseud = nAllelesThisFrom.notAllPseud[1],
                                  nGenesThisFrom = nGenesThisFrom[1], 
                                  nGenesThisFrom.notAllPseud = nGenesThisFrom.notAllPseud[1],
                                  nGenesThisFrom.noSwNoPseud.invariant = nGenesThisFrom.noSwNoPseud.invariant[1],
                                  nGenesThisFrom.noSwNoPseud.variant = nGenesThisFrom.noSwNoPseud.variant[1]),
                              by = .(displayname, from)]
  
  # TO specific codons (not sure I care but whatever)
  setkey(nsumm.trfromto, displayname, to)
  ## ANY mutation
  trto.all<-nsumm.trfromto[, .(AAToCod = AAToCod[1],
                               nAlleles = sum(nAlleles), 
                                 nGenesSwToThis = sum(nGenesThisSw),
                                 nGenesSwToThis.invariant = sum(nGenesThisSw.invariant),
                                 nGenesSwToThis.variant = sum(nGenesThisSw.variant),
                                 nGenesToThis.AllSameSw = sum(nGenesThisSw.allThisSw),
                               nGenesToThis.allThisSw.variant = sum(nGenesThisSw.allThisSw.variant),
                               nGenesToThis.notallThisSw = sum(nGenesThisSw.notallThisSw),
                                 nAllelesThisTo = nAllelesThisTo[1],
                                 nAllelesThisTo.notAllPseud = nAllelesThisTo.notAllPseud[1],
                                 nGenesThisTo = nGenesThisTo[1],
                                 nGenesThisTo.notAllPseud = nGenesThisTo.notAllPseud[1],
                               nGenesThisTo.noSwNoPseud.invariant = nGenesThisTo.noSwNoPseud.invariant[1],
                               nGenesThisTo.noSwNoPseud.variant = nGenesThisTo.noSwNoPseud.variant[1]),
                             by = .(displayname, to)]
  ##  Mutation that causes amino acid switch
  trto.aasw<-nsumm.trfromto[, .(AAToCod = AAToCod[1],
                                nAlleles = sum(nAlleles[AAFromCod!=AAToCod]), 
                                nGenesSwToThis = sum(nGenesThisSw[AAFromCod!=AAToCod]),
                                nGenesSwToThis.invariant = sum(nGenesThisSw.invariant[AAFromCod!=AAToCod]),
                                nGenesSwToThis.variant = sum(nGenesThisSw.variant[AAFromCod!=AAToCod]),
                                nGenesToThis.AllSameSw = sum(nGenesThisSw.allThisSw[AAFromCod!=AAToCod]),
                                nGenesToThis.allThisSw.variant = sum(nGenesThisSw.allThisSw.variant[AAFromCod!=AAToCod]),
                                nGenesToThis.notallThisSw = sum(nGenesThisSw.notallThisSw[AAFromCod!=AAToCod]),
                               nAllelesThisTo = nAllelesThisTo[1],
                               nAllelesThisTo.notAllPseud = nAllelesThisTo.notAllPseud[1],
                               nGenesThisTo = nGenesThisTo[1],
                               nGenesThisTo.notAllPseud = nGenesThisTo.notAllPseud[1],
                               nGenesThisTo.noSwNoPseud.invariant = nGenesThisTo.noSwNoPseud.invariant[1],
                               nGenesThisTo.noSwNoPseud.variant = nGenesThisTo.noSwNoPseud.variant[1]),
                           by = .(displayname, to)]
  
  # --- Return
  return(list(aafromto = nsumm.aafromto,
              aafrom = nsumm.aafrom,
              aato = nsumm.aato,
              tripfromto = nsumm.trfromto, # sw of any codon, even if same AA
              trfrom.all = trfrom.all,
              trfrom.aasw = trfrom.aasw,
              trto.all = trto.all,
              trto.aasw = trto.aasw
              ))
}

ksmwtests<-function(x, y){
  resmw<-wilcox.test(x, y)
  resks<-ks.test(x, y)
  
  out<-data.table(test = c("Mann-Whitney", "Kolmogorov-Smirnov"),
                  testshort = c("MW", "KS"),
                  median.1 = c(median(x, na.rm = T), median(x, na.rm = T)),
                  median.2 = c(median(y, na.rm = T), median(y, na.rm = T)),
                  statisticname = c("W", "D"),
                  statisticvalue = c(resmw$statistic, resks$statistic),
                  p.value = c(resmw$p.value, resks$p.value))
  return(out)
}

propcis<-function(x, n){
  if(x>=0 & n>=1){
    res<-binom.test(x, n)
  }else{
    res<-list(estimate = NaN, conf.int = c(NaN, NaN), p.value = NaN)
  }
 
  return(data.table(p = res$estimate,
                    p.low95ci = res$conf.int[1],
                    p.high95ci = res$conf.int[2],
                    binom.pval = res$p.value))
}

ntcomp<-function(nt){
  # Gives complement of ONE nucleotide
  ntdt<-data.table(orig = c("A", "C", "G", "T"),
                   comp = c("T", "G", "C", "A"))
  setkey(ntdt, orig)
  return(ntdt[nt, comp])
}

chartheoranticmuts<-function(swant, ac2aa){
  # Determines all the theoretical ways the alleleCM for each allele could get to the observed anticodon
  # In: swant, alleles that have anticodon switches to characterize. Must have columns displayname, tRNA, allelename, AA, Codon, AlleleCM
  #     ac2aa, canonical anticodon (as in tRNA DNA gene sequence) to amino acid mapping data.table, with all possibilities (Sup/SeC etc). Columns dnaAnticodon, AminoAcid
  # Out: list of data.tables:
  #     mutpaths: ALL possible mutation pathways for each input allele (codon -> AlleleCM). Multiple rows per input alleles. Possibly more than needed! Columns:
  #       observed.antic, input observed anticodon
  #       initialAA, input AA from
  #       proposed.init.antic, which anticodon for initialAA under consideration here
  #       nMutsFromInit, number mutations actual anticodon is from this theoretical ancestral one
  #       whichMutsFromInit, comma-collapsed positions in anticodon the mutations would be to get from proposed original to observed
  #       displayname - from input
  #       tRNA - gene name this is an allele of (from input)
  #       allelename - allele this is from (from input)
  #       AA - AA this observed codon codes for (from input)
  #     mpathsumm: SUMMARY of all possible mutation pathways. Columns:
  #       observed.antic, input observed anticodon
  #       n1MutPath, number of ways to get from observed to initial with 1 mutation
  #       which1MutPath, comma separated initial anticodons that have a 1 mutation path to observed mutation. Blank if none (as are all following 'which' columns)
  #       which1MutPathChange, comma separated places 1 mutation could be for the 1 mutation paths
  #       n2MutPath, number of ways to get from observed to initial with 2 mutations
  #       which2MutPath, comma separated initial anticodons that have a 2 mutation path to observed mutation
  #       which2MutPathChange,  comma separated places 2 mutations could be for the 2 mutation paths (separate paths ; separated)
  #       n3MutPath, number of ways to get from observed to initial with 3 mutations
  #       which3MutPath,  comma separated initial anticodons that have a 2 mutation path to observed mutation
  #       which3MutPathChange,  comma separated places 3 mutations could be for the 3 mutation paths (always 1,2,3)
  #       displayname - from input
  #       tRNA - gene name this is an allele of (from input)
  #       allelename - allele this is from (from input)
  #       AA - AA this observed codon codes for (from input)
  
  oneantimut<-function(obscod, fromAA, ac2aa){
    # Determines theoretical ways to get from all fromAA codons to obscod (observed codon)
    # In: obscod, anticodon sequence observed
    #     fromAA, amino acid this is hypothesized to have mutated from
    #     ac2aa, canonical anticodon (as in tRNA DNA gene sequence) to amino acid mapping data.table, with all possibilities (Sup/SeC etc). Columns dnaAnticodon, AminoAcid
    # Out: list of data.tables:
    #   mutpaths, one row for each comparison of actual anticodon with any that could code for fromAA. columns:
    #       observed.antic, input observed anticodon
    #       initialAA, input AA from
    #       proposed.init.antic, which anticodon for initialAA under consideration here
    #       nMutsFromInit, number mutations actual anticodon is from this theoretical ancestral one
    #       whichMutsFromInit, comma-collapsed positions in anticodon the mutations would be to get from proposed original to observed
    #   mutpathsumm,  SUMMARY of all possible mutation pathways. Columns:
    #       observed.antic, input observed anticodon
    #       n1MutPath, number of ways to get from observed to initial with 1 mutation
    #       which1MutPath, comma separated initial anticodons that have a 1 mutation path to observed mutation. Blank if none (as are all following 'which' columns)
    #       which1MutPathChange, comma separated places 1 mutation could be for the 1 mutation paths
    #       n2MutPath, number of ways to get from observed to initial with 2 mutations
    #       which2MutPath, comma separated initial anticodons that have a 2 mutation path to observed mutation
    #       which2MutPathChange,  comma separated places 2 mutations could be for the 2 mutation paths (separate paths ; separated)
    #       n3MutPath, number of ways to get from observed to initial with 3 mutations
    #       which3MutPath,  comma separated initial anticodons that have a 2 mutation path to observed mutation
    #       which3MutPathChange,  comma separated places 3 mutations could be for the 3 mutation paths (always 1,2,3)
    
    # Set up
    ocseq<-strsplit(obscod, "")[[1]]
    startcods<-ac2aa[AminoAcid==fromAA, dnaAnticodon] # all the possible codon starting points
    
    # Get possible mutation paths
    mutpath<-rbindlist(lapply(startcods, function(cod){
      cseq<-strsplit(cod, "")[[1]]
      compcod<-(ocseq==cseq)
      out<-data.table(observed.antic = obscod,
                      initialAA = fromAA, 
                      proposed.init.antic = cod,
                      nMutsFromInit = sum(!compcod),
                      whichMutsFromInit = paste(which(!compcod), collapse = ","))
    }))
    
    # Format/summarize/return
    msumm<-data.table(observed.antic = obscod, initialAA = fromAA,
                      n1MutPath = mutpath[, sum(nMutsFromInit==1)],
                      which1MutPath = mutpath[nMutsFromInit==1, paste(proposed.init.antic, collapse = ",")],
                      which1MutPathChange = mutpath[nMutsFromInit==1, paste(unique(whichMutsFromInit), collapse = ",")],
                      n2MutPath = mutpath[, sum(nMutsFromInit==2)],
                      which2MutPath = mutpath[nMutsFromInit==2, paste(proposed.init.antic, collapse = ",")],
                      which2MutPathChange = mutpath[nMutsFromInit==2, paste(unique(whichMutsFromInit), collapse = ";")], # differentiate the possibilities
                      n3MutPath = mutpath[, sum(nMutsFromInit==3)],
                      which3MutPath = mutpath[nMutsFromInit==3, paste(proposed.init.antic, collapse = ",")],
                      which3MutPathChange = mutpath[nMutsFromInit==3, paste(unique(whichMutsFromInit), collapse = ",")]) # always has to be 1,2,3
    
    return(list(mutpaths = mutpath,
                mutpathsumm = msumm))
  }
  
  # Get each existing anticodon - AlleleCM pair
  all.pairs<-unique(swant[, .(obscod = Codon, fromAA = AlleleCM)])
  
  # Get theoretical mutation paths for each
  mutall<-lapply(1:nrow(all.pairs), function(x){
    oneantimut(obscod = all.pairs[x, obscod], fromAA = all.pairs[x, fromAA], ac2aa)
  })
  mutpaths<-rbindlist(lapply(mutall, function(x) x$mutpaths))
  mpathsumm<-rbindlist(lapply(mutall, function(x) x$mutpathsumm))
  
  # Merge in with allele info
  setkey(mutpaths, observed.antic, initialAA)
  setkey(mpathsumm, observed.antic, initialAA)
  setkey(swant, Codon, AlleleCM)
  
  mutpaths<-mutpaths[swant[, .(displayname, tRNA, allelename, AA, Codon, AlleleCM)], allow.cartesian = T]
  mpathsumm<-mpathsumm[swant[, .(displayname, tRNA, allelename, AA, Codon, AlleleCM)]] # now Codon called observed.antic; AlleleCM is now initialAA
  
  # Return
  return(list(mutpaths = mutpaths,
              mpathsumm = mpathsumm))
}

realanticmuts<-function(mutpaths, varinfo){
  # For alleles with anticodon mutation switches, takes the possible mutations and combines with observed in pop variants (where they exist)
  #     to predict mutation paths. If NO variants (e.g., fixed switch), picks most likely path(s) and notes this
  # In: mutpaths, $mutpaths output of chartheoranticmuts
  #     varinfo, info for variants in tRNAs with relative position/structure info. From Path to file with tRNA gene body info output by get_strain_variants_relpos.py
  # Out: data.table with one row per input ALLELE (collapsed across mutpaths). Columns:
  #         #     displayname, from input
  #     tRNA, from input
  #     allelename, from input
  #     AA, from input (amino acid this allele's anticodon matches)
  #     initialAA, from input (amino acid backbone)
  #     observed.antic, actual observed anticodon in allele
  #     proposed.init.antic, anticodon this one mutated FROM (multiple if can't tell apart) - predicted or observed when possible
  #     nMutsFromInit, number of mutations in proposed initial anticodon to get to eventual observed anticodon
  #     whichMutsFromInit, which spot(s) in anticodon had to mutate to get to observed (; separated if multiple possibilities)
  #     varObsInData, was this variant observed in population variation data, T or F?
  #     FLAG.anticodVarButNoMatch : if T, there WAS an anticodon mutation flagged but it didn't match data - have seen this when, e.g, indels in anticodon arms
  #     nPossibleMutations, count of # mutations that could lead to what we see (can't narrow further)
  #     varObsPos, if anticodon variant observed - its genomic position
  #     varObs.tRNApos, if anticodon variant observed - its position in tRNA
  #     varObs.ref, if anticodon variant observed  - its ref allele
  #     varObs.alt, if anticodon variant observed  - its alt allele
  #     varObs.anticPos, if anticodon variant observed - its position in anticodon
  #     varObs.nHomAlt, if anticodon variant observed - n strains homozygous alt allele here
  #     nNonAnticVariants, if this allele has variants NOT in anticodon, how many
  #     structureNonAnticVariants, if this allele has variants NOT in anticodon, where are they (,-separated)
  #     posNonAnticVariants, if this allele has variants NOT in anticodon, genomic position (,-separated)
  #     consistentAnticMutAcrossAlleles, T if all alleles of this gene in input have same predicted [or obs] anticodon mutations; F otherwise
  
  # Subfunctions
  realoneg<-function(mutpathsone, onevarinfo){
    # Gets best mutation paths for ONE tRNA's alleles
    # In: mutpathsone, whole rows for one observed tRNA. One row per switch allele for that tRNA.
    #     onevarinfo, variants in this tRNA. 0 or more rows
    # Out: data.table with one row per allele (I think) 
    #     displayname, from input
    #     tRNA, from input
    #     allelename, from input
    #     AA, from input (amino acid this allele's anticodon matches)
    #     initialAA, from input (amino acid backbone)
    #     observed.antic, actual observed anticodon in allele
    #     proposed.init.antic, anticodon this one mutated FROM (multiple if can't tell apart) - predicted or observed when possible
    #     nMutsFromInit, number of mutations in proposed initial anticodon to get to eventual observed anticodon
    #     whichMutsFromInit, which spot(s) in anticodon had to mutate to get to observed (; separated if multiple possibilities)
    #     varObsInData, was this variant observed in population variation data, T or F?
    #     FLAG.anticodVarButNoMatch : if T, there WAS an anticodon mutation flagged but it didn't match data - have seen this when, e.g, indels in anticodon arms
    #     nPossibleMutations, count of # mutations that could lead to what we see (can't narrow further)
    #     varObsPos, if anticodon variant observed - its genomic position
    #     varObs.tRNApos, if anticodon variant observed - its position in tRNA
    #     varObs.ref, if anticodon variant observed  - its ref allele
    #     varObs.alt, if anticodon variant observed  - its alt allele
    #     varObs.anticPos, if anticodon variant observed - its position in anticodon
    #     varObs.nHomAlt, if anticodon variant observed - n strains homozygous alt allele here
    #     nNonAnticVariants, if this allele has variants NOT in anticodon, how many
    #     structureNonAnticVariants, if this allele has variants NOT in anticodon, where are they (,-separated)
    #     posNonAnticVariants, if this allele has variants NOT in anticodon, genomic position (,-separated)
    
    setkey(mutpathsone, allelename)
    # --- If no observed anticodon variants, switch is presumably fixed: do prediction of mutation path
    if(!"anticodon"%in%onevarinfo$substructure){
      # Predicted anticodon switch
      mutpathsone[, isMin:=nMutsFromInit==min(nMutsFromInit), by = allelename]
      predsw<-mutpathsone[isMin==T, .SD, .SDcols = -c(ncol(mutpathsone))]
      ## If more than one predicted one, collapse to one row per allele
      predsw.o<-predsw[, .(nPossibleMutations = .N, observed.antic = observed.antic[1], proposed.init.antic = paste(proposed.init.antic, collapse = ","), 
                           nMutsFromInit = paste(nMutsFromInit, collapse = ","), whichMutsFromInit = paste(whichMutsFromInit, collapse = ";"), 
                           displayname = displayname[1], tRNA = tRNA[1], initialAA = initialAA[1], AA = AA[1]), by = allelename]
       # add more columns to match where there IS an observation here
      predsw.o[, `:=`(varObsInData=F, FLAG.anticodVarButNoMatch = F, varObsPos = "", varObs.tRNApos = "", varObs.ref = "", varObs.alt = "", varObs.anticPos = "", varObs.nHomAlt = "")]
      
      # Are there observed mutations NOT in anticodon? If yes, note them per allele
      if(nrow(onevarinfo)>0){
        # Get into per allele (instead of per variant) format. Here, only care about the ones where someone is homAlt. If strain is Reference, by default has no variants vs. ref
        predsw.o[, sname:=predsw.o[, tstrsplit(allelename, "_")]$V2]
        varsch<-rbindlist(lapply(predsw.o$sname, function(strn){
          if(strn=="Reference"){
            myvars<-onevarinfo[0]
          }else{
            myvars<-onevarinfo[grepl(strn, homAlt), .SD, .SDcols = -c((ncol(onevarinfo)-2): ncol(onevarinfo))]
          }
          out<-myvars[,.(sname = strn, nNonAnticVariants = .N, structureNonAnticVariants = paste(structure, substructure, sep = ":", collapse = ","), posNonAnticVariants = paste(pos, collapse = ","))]
          return(out)
        }))
        setkey(varsch, sname)
        setkey(predsw.o, sname)
        
        # combine, format for return
        predsw.o<-predsw.o[varsch]
        predsw.o[, sname:=NULL]
      }else{ # # If not, say not
        predsw.o[, `:=`(nNonAnticVariants = 0, structureNonAnticVariants = "", posNonAnticVariants = "")]
      }
    }else{ # --- At least one Variant in anticodon
      # Deal with anticodon variant(s)
      antvar<-onevarinfo[substructure=="anticodon", ]
      setkey(mutpathsone, allelename) # do for each allele in case of multiple
      predsw.o<-rbindlist(lapply(mutpathsone[, unique(allelename)], function(alname){
        onemut<-mutpathsone[alname]
        onemut[, sname:=onemut[,tstrsplit(allelename, "_")]$V2]
        # Anticodon variants for this allele
        oneant<-antvar[grepl(onemut[1, sname], homAlt) & nchar(ref)==1 & nchar(alt)==1, # nchars are for screening out indel
                         .SD, .SDcols = -c((ncol(onevarinfo)-2): ncol(onevarinfo))]
        posspred<-onemut[nMutsFromInit==nrow(oneant)]
        if(nrow(posspred)==0){ # flag if they don't match OR if indel
          # this CAN happen - have an anticodon mutation but it doesn't match. going to FLAG this and then go on. (e.g., anticodon arm indels change where anticodon is)
          # alternatively can't deal with anticodon mutation that is indel
          onemut[, isMin:=nMutsFromInit==min(nMutsFromInit), by = allelename]
          predsw<-onemut[isMin==T, .SD, .SDcols = -c(ncol(onemut))]
          predsw.o<-predsw[, .(nPossibleMutations = .N, observed.antic = observed.antic[1], proposed.init.antic = paste(proposed.init.antic, collapse = ","), 
                               nMutsFromInit = paste(nMutsFromInit, collapse = ","), whichMutsFromInit = paste(whichMutsFromInit, collapse = ";"), 
                               displayname = displayname[1], tRNA = tRNA[1], initialAA = initialAA[1], AA = AA[1]), by = allelename]
          predsw.o[, `:=`(varObsInData=F, FLAG.anticodVarButNoMatch = T,
                         varObsPos = "", varObs.tRNApos = "", varObs.ref = "", varObs.alt = "", varObs.anticPos = "", varObs.nHomAlt = "")]
        }else{
          ## Annotate observed mutations 
          oneant[tRNA_strand=="+", `:=`(nt.from = ref, nt.to = alt)]
          oneant[tRNA_strand=="-", `:=`(nt.from = ntcomp(ref), nt.to = ntcomp(alt))]
          ## Annotate predicted mutations - expand to multiple rows if relevant
          posspredmuts<-posspred[, .(observed.antic = observed.antic, initialAA = initialAA, proposed.init.antic = proposed.init.antic,
                                     nMutsFromInit = nMutsFromInit, whichMutsFromInit = as.integer(strsplit(whichMutsFromInit, ",")[[1]]),
                                     displayname = displayname, tRNA = tRNA, allelename = allelename, AA = AA,
            nt.from = as.vector(sapply(as.integer(strsplit(whichMutsFromInit, ",")[[1]]), function(x) substr(proposed.init.antic, x, x))),
            nt.to = as.vector(sapply(as.integer(strsplit(whichMutsFromInit, ",")[[1]]), function(x) substr(observed.antic, x, x))))]
          
          ## Pick overlapping ones
          setkey(posspredmuts, whichMutsFromInit, nt.from, nt.to)
          setkey(oneant, substructure_pos, nt.from, nt.to)
          predsw<-merge(posspredmuts, oneant)
          if(nrow(predsw)==0){ # Sometimes there ARE mutations but they DON'T make sense with what is called
           #  cat(paste("For allele", alname, "there are observed anticodon mutation(s) that DO NOT match the predicted ones, FYI...\n"))
            predsw.o<-posspredmuts[, .(nPossibleMutations = .N, observed.antic = observed.antic[1], proposed.init.antic = paste(proposed.init.antic, collapse = ","), 
                                 nMutsFromInit = paste(nMutsFromInit, collapse = ","), whichMutsFromInit = paste(whichMutsFromInit, collapse = ";"), 
                                 displayname = displayname[1], tRNA = tRNA[1], initialAA = initialAA[1], AA = AA[1]), by = allelename]
            predsw.o[, `:=`(varObsInData=F, FLAG.anticodVarButNoMatch = T,
                            varObsPos = "", varObs.tRNApos = "", varObs.ref = "", varObs.alt = "", varObs.anticPos = "", varObs.nHomAlt = "")]
          }else{
            # Format for return [collapse-y stuff did above, note some stuff about actual observed mutation]
            predsw.o<-predsw[, .(# predicted stuff
              allelename = allelename, nPossibleMutations = .N, observed.antic = observed.antic[1], proposed.init.antic = paste(proposed.init.antic, collapse = ","), 
              nMutsFromInit = paste(nMutsFromInit, collapse = ","), whichMutsFromInit = paste(whichMutsFromInit, collapse = ";"), 
              displayname = displayname.x[1], tRNA = tRNA.x[1], initialAA = initialAA[1], AA = AA[1],
              # observed stuff
              varObsInData = T, FLAG.anticodVarButNoMatch = F,
              varObsPos = paste(pos, collapse = ","), varObs.tRNApos = paste(tRNA_pos, collapse = ","),
              varObs.ref = paste(ref, collapse = ","), varObs.alt = paste(alt, collapse = ","),
              varObs.anticPos = paste(substructure_pos, collapse = ","), varObs.nHomAlt = paste(nHomAlt, collapse = ","))]
            
            # (if they DO NOT match, have a way to NOTE THIS. prob return predicted + new columns??) ************** COME BACK & IMPLEMENT THIS
          }
        }
        return(predsw.o)
      }))
      
      # Deal with any other variant(s) [should mirror what's done above when no anticodon mutations]
      nonantvar<-onevarinfo[substructure!="anticodon", ]
      alsnonant<-unique(mutpathsone[, .(displayname, tRNA, allelename)])
      if(nrow(nonantvar)==0){ # no other variants, note this
        alsnonant[, `:=`(nNonAnticVariants = 0, structureNonAnticVariants = "", posNonAnticVariants = "")]
      }else{ # there are other variants, note/catalog this [like did before]
        # ... just need to do per allele .... then merge back in with other per allele
        alsnonant[, sname:=alsnonant[, tstrsplit(allelename, "_")]$V2]
        varsch<-rbindlist(lapply(alsnonant$sname, function(strn){
          if(strn=="Reference"){
            myvars<-nonantvar[0]
          }else{
            myvars<-nonantvar[grepl(strn, homAlt), .SD, .SDcols = -c((ncol(nonantvar)-2): ncol(nonantvar))]
          }
          out<-myvars[,.(sname = strn, nNonAnticVariants = .N, structureNonAnticVariants = paste(structure, substructure, sep = ":", collapse = ","), posNonAnticVariants = paste(pos, collapse = ","))]
          return(out)
        }))
        setkey(varsch, sname)
        setkey(alsnonant, sname)
        alsnonant<-alsnonant[varsch]
        alsnonant[, sname:=NULL]
      }
      setkey(predsw.o, displayname, tRNA, allelename)
      setkey(alsnonant, displayname, tRNA, allelename)
      predsw.o<-predsw.o[alsnonant]
    }
    
    # --- Return
    setcolorder(predsw.o, c("displayname", "tRNA", "allelename", "AA", "initialAA", "observed.antic", "proposed.init.antic",
                            "nMutsFromInit", "whichMutsFromInit", "varObsInData", "FLAG.anticodVarButNoMatch","nPossibleMutations", "varObsPos", "varObs.tRNApos", "varObs.ref", "varObs.alt",
                            "varObs.anticPos", "varObs.nHomAlt", "nNonAnticVariants", "structureNonAnticVariants", "posNonAnticVariants"))
    return(predsw.o)
  }
  
  # Run this per tRNA
  bestmuts<-rbindlist(lapply(unique(mutpaths$displayname), function(sp){
    rbindlist(lapply(unique(mutpaths[displayname==sp, tRNA]), function(t){
      # print(c(sp, t))
      out<-realoneg(mutpathsone = mutpaths[displayname==sp & tRNA==t],
                    onevarinfo = varinfo[displayname==sp & tRNA==t])
      return(out)
    }))
  }))
  
  # Add if mutation is consistent across alleles (to best of knowledge) **for ANTICODON mutations
  setkey(bestmuts, displayname, tRNA)
  toctmut<-unique(bestmuts[, .(proposed.init.antic, nMutsFromInit, whichMutsFromInit, varObsInData, FLAG.anticodVarButNoMatch, 
                        nPossibleMutations, varObsPos, varObs.tRNApos), by = .(displayname, tRNA)]) # more fields exist but this should be plenty
  setkey(toctmut, displayname, tRNA)
  ctmutposs<-toctmut[, .N, by = .(displayname, tRNA)] # more fields exist but this should be plenty
  
  bestmuts<-bestmuts[ctmutposs]
  bestmuts[, consistentAnticMutAcrossAlleles:=ifelse(N==1, T, F)]
  bestmuts[, N:=NULL]
  
  # Return!
  return(bestmuts)
}

backbonemuts<-function(alsbackbone, varinfo){
  # Characterizes the mutations that make up cases where tRNAscan-SE predicts backbone mutations cause isotype switching
  # In: alsbackbone, all tRNA alleles predicted to have a backbone mutation. Columns must include:
  #       displayname, tRNA, allelename, <AA descriptive stuff>, n.allele
  #     varinfo, info for variants in tRNAs with relative position/structure info. From Path to file with tRNA gene body info output by get_strain_variants_relpos.py
  # Out: List of data.tables:
  #     vardetail, details of all variants observed for each allele. One row per allele per variant. Columns:
  #         allele descriptors from input (same for each row of an allele): displayname, shortname, tRNA, allelename, AA, Codon, AlleleCM, Codon.orig, AlleleCM.orig, Infernal, IsotypeScore, Note, VariableInPop, Classification, n.allele, n.called, n.missing, freq.ofcalled, freq.ofallinclmissing, all.altered, all.lost, switch, any.switch, nAlCM, whatMut, 
  #         variant descriptors from input (different for each row of an allele): chrom, pos, ref, alt, nHomRef, nHomAlt, nNotMissingHet, nMissingHet, tRNA_strand, tRNA_pos, structure, substructure, nSubStructure, substructure_pos
  #     varsumm, summary of variants per allele, one row per allele. Columns:
  #       allele descriptors from input: displayname, tRNA, allelename, AA, Codon, AlleleCM, Infernal, IsotypeScore, Note, VariableInPop, Classification, n.allele, n.called, n.missing, freq.ofcalled, freq.ofallinclmissing, all.altered, all.lost, switch, any.switch, nAlCM, whatMut
  #       nVars, number variants (from VCF) observed in this allele
  #       nStructuresWithVars, number structures variants come from
  #       nSubstructsWithVars, number structure-substructure pairs variants come from
  #       nAcceptorVars, number variants from acceptor stem
  #       nDVars, number variants from d arm
  #       nTVars, number variants from t arm
  #       nVLVars, number variants from variable loop
  #       nAntArmVars, number variants from anticodon arm (arm, loop, actual anticodon)
  #       nOtherStructVars, number variants from any other called tRNA structure
  #       nNotCalledInStructVars, count where there are variants that I couldn't assign to a tRNA substruct
  #       varSubstructs, comma-separated list of structure-substructure pair variants in this allele are in
  #       flagNoObsVar, flag: T if all of these are 0s/ no observed variants
  
  # Subfunctions
  bboneoneg<-function(oneal, onevarinfo){
    # Summarizes backbone mutations for ONE allele
    # In: oneal, one row of alsbackbone, oneal 
    #     onevarinfo, variants falling in this gene
    # Out: List of data.tables:
    #     vardetail, details of all variants in this allele. One row per variant. Columns:
    #         allele descriptors from input (same for each row; may be more than this): displayname, shortname, tRNA, allelename, AA, Codon, AlleleCM, Codon.orig, AlleleCM.orig, Infernal, IsotypeScore, Note, VariableInPop, Classification, n.allele, n.called, n.missing, freq.ofcalled, freq.ofallinclmissing, all.altered, all.lost, switch, any.switch, nAlCM, whatMut, 
    #         variant descriptors from input (different for each row): chrom, pos, ref, alt, nHomRef, nHomAlt, nNotMissingHet, nMissingHet, tRNA_strand, tRNA_pos, structure, substructure, nSubStructure, substructure_pos
   #              NB even if no variants seen will have one row, with NAs
    #     varsumm, summary of variants in this allele. Columns:
    #       allele descriptors from input (may be more than this): displayname, tRNA, allelename, AA, Codon, AlleleCM, Infernal, IsotypeScore, Note, VariableInPop, Classification, n.allele, n.called, n.missing, freq.ofcalled, freq.ofallinclmissing, all.altered, all.lost, switch, any.switch, nAlCM, whatMut
    #       nVars, number variants (from VCF) observed in this allele
    #       nStructuresWithVars, number structures variants come from
    #       nSubstructsWithVars, number structure-substructure pairs variants come from
    #       nAcceptorVars, number variants from acceptor stem
    #       nDVars, number variants from d arm
    #       nTVars, number variants from t arm
    #       nVLVars, number variants from variable loop
    #       nAntArmVars, number variants from anticodon arm (arm, loop, actual anticodon)
    #       nOtherStructVars, number variants from any other called tRNA structure
    #       nNotCalledInStructVars, count where there are variants that I couldn't assign to a tRNA substruct
    #       varSubstructs, comma-separated list of structure-substructure pair variants in this allele are in
    #       flagNoObsVar, flag: T if all of these are 0s/ no observed variants
    
    # Data set up for easy merging
    setkeyv(oneal, names(oneal)[names(oneal)%in%names(onevarinfo)]) # for easy merging
    setkeyv(onevarinfo, names(oneal)[names(oneal)%in%names(onevarinfo)])
    
    # Get actual mutations for this allele, plan to return
    oneal[, strain:=oneal[, tstrsplit(allelename, "_")]$V2]
    varsal<-merge(
      oneal[, .SD, .SDcols = -ncol(oneal)], # allele info
      onevarinfo[grepl(oneal$strain, homAlt), (.SD), .SDcols = -c( (ncol(onevarinfo)-2):ncol(onevarinfo) )] # variant info
      )
    
    # Summarize mutations (number, what structure they're in, etc)
    varsumm<-varsal[, .( # Housekeeping/keep allele info
                        displayname, tRNA, allelename, AA, Codon, AlleleCM, Infernal, IsotypeScore, Note, VariableInPop, Classification, n.allele,
                        n.called, n.missing, freq.ofcalled, freq.ofallinclmissing, all.altered, all.lost, switch, any.switch, nAlCM, whatMut,
                        gene.classif,
                        # Variant SUMMARY
                        nVars = .N, nStructuresWithVars = length(unique(structure)), nSubstructsWithVars = length(unique(paste(structure, substructure, sep = "-"))),
                        nAcceptorVars = sum(structure=="acceptor_stem", na.rm = T), nDVars = sum(structure=="d", na.rm = T), 
                        nTVars = sum(structure=="t", na.rm = T), nVLVars = sum(structure=="variable_loop", na.rm = T),
                        nAntArmVars = sum(structure=="anticodon", na.rm = T),
                        nOtherStructVars = sum(!structure%in%c("acceptor_stem", "d", "t", "variable_loop", "anticodon", NA)), 
                        nNotCalledInStructVars = sum(is.na(structure)), # KEY - sometimes its vars that I couldn't classify
                        varSubstructs = paste(structure, substructure, sep = "-", collapse = ","),
                        # Flag if no vars observed
                        flagNoObsVar = is.na(pos[1]))] # no obs mutations, flag this [ don't think it'll happen but prepping for edge case ]
    
    # Return (both)
    return(list(vardetail = varsal,
                varsumm = varsumm))
  }
  
  # Get var info for all alleles
  l.peral<-lapply(unique(alsbackbone$displayname), function(sp){
    lapply(alsbackbone[displayname==sp, unique(allelename)], function(al){
      return(bboneoneg(oneal = alsbackbone[displayname==sp & allelename==al], 
                       onevarinfo = varinfo[tRNA==alsbackbone[displayname==sp & allelename==al, tRNA]]))
    })
  })
  
  # Format & return
  vardetail<-rbindlist(lapply(l.peral, function(x) rbindlist(lapply(x, function(y) y$vardetail))))
  varsumm<-rbindlist(lapply(l.peral, function(x) rbindlist(lapply(x, function(y) y$varsumm))))
  return(list(vardetail = vardetail, 
              varsumm = varsumm))
}

#### Functions: plotting ####
stackedbar<-function(pinhclass, mycolors, xcol = "strain", stackcol = "propName", stacknumcol = "prop", xorder = NA,
                     legendlabel = "", myxlab = "Strain", myylab = "Proportion genes", mytitle = "", mysubt = "",
                     myleglabels = c("", "")){
  # Modified from chrlocenrichment_asederpim.R + used many other places
  # Makes stacked bar chart, one stacked bar per category in x col. Copied from other scripts.
  # In: pinhclass, data.table with all data to plot. Must have columns named with values of xcol, stackcol, stacknumcol
  #     mycolors, colors for bars: vector of length of categories in data; also used as ORDER for bars!! Names are categories, values are colors
  #     xcol, column containing category for X axis of bar chart
  #     stackcol, column containing category to split bars into
  #     stacknumcol, column containing numbers to plot as stacked bar
  #     xorder, OPTIONAL order in which to arrange categories on x axis (in xcol)
  #     legendlabel
  #     myxlab, myylab, mytitle, mysubt: plot labels as intuitively named
  #     myleglabels, provide legend label names (**in order they're plotted **) - in exorder!
  # Out: ggplot2 barplot object. Can facet this externally!   
  
  # Format data
  pdata<-copy(pinhclass)
  pdata[,mystackcol:=factor(pdata[,get(stackcol)], levels = names(mycolors))] # relevel factor
  if(!is.na(xorder[1])){
    # Re-level factor 
    pdata[,xcol:=factor(pdata[,get(xcol)], levels = xorder)]
  }else{pdata[,xcol:=get(xcol)]}
  
  plt<-ggplot(pdata, aes(xcol, eval(as.name(stacknumcol)))) +
    geom_bar(aes(fill = mystackcol), stat = "identity", position = "stack") + labs(fill = legendlabel) +
    xlab(myxlab) + ylab(myylab) + 
    scale_fill_manual(values = mycolors, labels = myleglabels) + # provided colors
    scale_y_continuous(expand = c(0, 0)) +
    ggtitle(mytitle, subtitle = mysubt) +
    theme_bw() + myggtheme
  
  return(plt)
}

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
                help = "Path to *alleleinfo_counts_worigcodonetc.txt output of show_mutational_variation.R. For all species in species info file combined.",
                type = "character")
p<-add_argument(p, "--genevars",
                help = "EXAMPLE Path to file with tRNA gene body info output by get_strain_variants_relpos.py Where species ID/species specific info is, put SAMP instead",
                type = "character")
p<-add_argument(p, "--anticodontable",
                help = "Path to anticodon table with columns dnaAnticodon and AminoAcid. Should map directly to tRNA anticodons here, including optional readings of CAT as Met or iMet and of TCA, TTA, CTA as SeC",
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
alinfo<-fread(p$alleleinfo, header = T)

# Add for if all lost/pseud
setkey(alinfo, displayname, tRNA)
alinfo<-alinfo[alinfo[, .(all.lost = ifelse(sum(Classification=="Lost")==.N, T, F)), by = .(displayname, tRNA)]]

# Variants info (keep only ones in genes that have any isotype switches...that's added later?)
varinfo<-fread.mult(sinfo, exampfile = p$genevars)

#### Narrowing/analyzing: Get numbers, proportions ####
isws<-alinfo[AA!=AlleleCM & Codon!="NNN" & Classification!="Lost"] # isotype switches that are NOT pseudos
alinfo[, switch:=(AA!=AlleleCM & Codon!="NNN" & Classification!="Lost")] # so have in whole thing; added later
alinfo[, any.switch:=ifelse(sum(switch)>0, T, F), by = tRNA] # Flag if ANY allele at this gene has a switch
setkey(alinfo, displayname)

# --- Just a little count summary: overall
setkey(isws, displayname)
nsumm<-isws[, .(nSwitchAls = .N, 
                nGenesWithSwitches = length(unique(tRNA)),
                nGenesWSw.invariant = length(unique(tRNA[VariableInPop==F])),
                nGenesWSw.variant = length(unique(tRNA[VariableInPop==T])),
                nGenesWSw.allSw = length(unique(tRNA[all.altered==T])),
                nGenesWSw.notallSw = length(unique(tRNA[all.altered==F])),
                nGenesWSw.allSw.variant = length(unique(tRNA[all.altered==T & VariableInPop==T]))),
            by = displayname]
nsumm[, `:=`(nAlleles = alinfo[, .N, by = displayname]$N,
             nAllelesNotAllPseud = alinfo[all.lost==F, .N, by = displayname]$N,
             nAllelesNonSwitch = alinfo[, .N, by = displayname]$N - nSwitchAls,
             nAllelesNonSwitchNotAllPseud = alinfo[all.lost==F, .N, by = displayname]$N - nSwitchAls,
             nGenes = alinfo[, length(unique(tRNA)), by = displayname]$V1,
             nGenesNotAllPseud = alinfo[all.lost==F, length(unique(tRNA)), by = displayname]$V1,
             nGenesNoSwitches = alinfo[, length(unique(tRNA)), by = displayname]$V1 - nGenesWithSwitches,
             nGenesNoSwitchesNonPseud = alinfo[all.lost==F, length(unique(tRNA)), by = displayname]$V1 - nGenesWithSwitches,
             nGenesNoSwitchesNonPseud.variant = alinfo[all.lost==F & any.switch==F & VariableInPop==T, length(unique(tRNA)), by = displayname]$V1,
             nGenesNoSwitchesNonPseud.invariant = alinfo[all.lost==F & any.switch==F & VariableInPop==F, length(unique(tRNA)), by = displayname]$V1)]
# Save
write.table(nsumm, file.path(p$outdir, paste0(p$baseoutname, "_nisotypeswitchsummary_all.txt")),
            quote = F, row.names = F, sep = "\t")

# --- Same count summary: by AA and triplet
# Set up all amino acids, triplets that are observed in this dataset that we want to investigate
aas<-alinfo[, sort(unique(c(AlleleCM, AA)))]
aas<-aas[aas!="Undet"] #excluding pseuds from this particular analysis of switching
trips<-unique(alinfo[, .(AA, Codon)])[AA%in%aas & Codon!="NNN", ] # exclude pseud
      # NOTE do have cases where same codon gets two amino acid - SeC/Sup and Met/iMet...
trips<-trips[order(AA)]

# Get and save counts
nsumm.aacod<-summbyaa(alinfo, aas, trips)

aadir<-file.path(p$outdir, "byaacodon")
if(!dir.exists(aadir)){dir.create(aadir, recursive = T)}
invisible(lapply(names(nsumm.aacod), function(dn){
  write.table(nsumm.aacod[[dn]],
              gzfile(file.path(aadir, paste0(p$baseoutname, "_", dn, ".txt.gz"))),
              sep = "\t", quote = F, row.names = F)
}))

# Summarize: by various types, which switches observed in which number/sets of SPECIES (kinda age-of-switch related)

# Summarize: for triplet switches, switching within amino acid class vs. between them as best can detect....
#       just one off don't see any in the minority of cases where you'd stay coding for same amino acid

# --- proportions and CIs of all of this this, like this and in long format?: FROM
# Focusing on key: proportion of all (non-all-pseud) genes broken down into if they have multiple detected alleles; if these have sw & not
## Set up this information for posterity; probably just as easy to not set it up but w/e
alcatnames<-data.table(longname = c("all switch - fixed",
              "all switch - variable",
              "some switches, some not",
              "no switch - fixed",
              "no switch - variable"),
              shortname = c("alcat1", "alcat2", "alcat3", "alcat4", "alcat5")) # alcat1.. etc in below
propinfo<-data.table(dataset = c("all", "per AA", "per codon"),
                     datname = c(expression(nsumm), expression(nsumm.aacod$aafrom), expression(nsumm.aacod$trfrom.all)), 
                     datby = c(expression("displayname"), expression(c("displayname", "from")), expression(c("displayname", "from"))), # data should already be keyed by this
                     alcat1 = c("nGenesWSw.invariant", "nGenesSwFromThis.invariant", "nGenesSwFromThis.invariant"),
                     alcat2 = c("nGenesWSw.allSw.variant", "nGenesFromThis.allThisSw.variant", "nGenesFromThis.allThisSw.variant"),
                     alcat3 = c("nGenesWSw.notallSw", "nGenesFromThis.notallThisSw", "nGenesFromThis.notallThisSw"),
                     alcat4 = c("nGenesNoSwitchesNonPseud.invariant", "nGenesThisFrom.noSwNoPseud.invariant", "nGenesThisFrom.noSwNoPseud.invariant"),
                     alcat5 = c("nGenesNoSwitchesNonPseud.variant", "nGenesThisFrom.noSwNoPseud.variant", "nGenesThisFrom.noSwNoPseud.variant"))
## Get these nums and proportions
longnps<-rbindlist(lapply(1:nrow(propinfo), function(i){
  # Get ns
  dat<-eval(propinfo[i, datname])
  out<-rbindlist(lapply(1:nrow(alcatnames), function(x){
    data.table(whatn = alcatnames[x, longname],
               dat[, get(propinfo[i, get(alcatnames[x, shortname])]), by = eval(eval(propinfo[i, datby]))])
  }))
  setnames(out, c("whatn",eval(propinfo[i, datby]), "n"))
  ## Add SUM of all (denominator)
  setkeyv(out, eval(propinfo[i, datby]))
  denom<-out[, sum(n), by = eval(eval(propinfo[i, datby]))]
  setnames(denom, "V1", "totalN")
  setkeyv(denom, eval(propinfo[i, datby]))
  out<-out[denom]
  
  # Add in ps
  out<-data.table(out, 
                  rbindlist(lapply(1:nrow(out), function(x){
                    propcis(x = out[x, n], n = out[x, totalN])[, .(p, p.low95ci, p.high95ci)]
                  })))
  
  # Names - from should be 'all' if not included
  if(!"from"%in%names(out)){
    out[, from:="all"]
  }
  setcolorder(out, "from")
  
  # Return [make sure same annotation info in all]
  return(out)
}))
## Save
write.table(longnps, file.path(p$outdir, paste0(p$baseoutname, "_longformNsPs_from_genes.txt")),
            sep = "\t", quote = F)

# --- proportions and CIs of all 'TO' amino acids, in long format
# i.e., for all codons that end up saying X amino acid, how many are switches
propinfo.to<-data.table(dataset = c("all", "per AA", "per codon"),
                     datname = c(expression(nsumm), expression(nsumm.aacod$aato), expression(nsumm.aacod$trto.all)), 
                     datby = c(expression("displayname"), expression(c("displayname", "to")), expression(c("displayname", "to"))), # data should already be keyed by this
                     alcat1 = c("nGenesWSw.invariant", "nGenesSwToThis.invariant", "nGenesSwToThis.invariant"),
                     alcat2 = c("nGenesWSw.allSw.variant", "nGenesToThis.allThisSw.variant", "nGenesToThis.allThisSw.variant"),
                     alcat3 = c("nGenesWSw.notallSw", "nGenesToThis.notallThisSw", "nGenesToThis.notallThisSw"),
                     alcat4 = c("nGenesNoSwitchesNonPseud.invariant", "nGenesThisTo.noSwNoPseud.invariant", "nGenesThisTo.noSwNoPseud.invariant"),
                     alcat5 = c("nGenesNoSwitchesNonPseud.variant", "nGenesThisTo.noSwNoPseud.variant", "nGenesThisTo.noSwNoPseud.variant"))

longto<-rbindlist(lapply(1:nrow(propinfo.to), function(i){
  # Get ns
  dat<-eval(propinfo.to[i, datname])
  out<-rbindlist(lapply(1:nrow(alcatnames), function(x){
    data.table(whatn = alcatnames[x, longname],
               dat[, get(propinfo.to[i, get(alcatnames[x, shortname])]), by = eval(eval(propinfo.to[i, datby]))])
  }))
  setnames(out, c("whatn",eval(propinfo.to[i, datby]), "n"))
  ## Combine sups and SeCs here - Sups are isotype switching Secs [codons should work on their own]
  if(propinfo.to[i, dataset]=="per AA"){
    tocomb<-out[to%in%c("SeC", "Sup")]
    setkey(tocomb, whatn, displayname) # MESSES UP IF NOT KEYED BY DISPLAYNAME, FIX LATER
    scomb<-tocomb[, .(n = sum(n)), by = .(whatn, displayname)]
    scomb[, to:="SeC"]
    setcolorder(scomb, names(out))
    out<-rbind(out[!to%in%c("SeC", "Sup")], scomb)
  }
  ## Add SUM of all (denominator)
  setkeyv(out, eval(propinfo.to[i, datby]))
  denom<-out[, sum(n), by = eval(eval(propinfo.to[i, datby]))]
  setnames(denom, "V1", "totalN")
  setkeyv(denom, eval(propinfo.to[i, datby]))
  out<-out[denom]
  
  # Add in ps
  out<-data.table(out, 
                  rbindlist(lapply(1:nrow(out), function(x){
                    propcis(x = out[x, n], n = out[x, totalN])[, .(p, p.low95ci, p.high95ci)]
                  })))
  
  # Names - from should be 'all' if not included
  if(!"to"%in%names(out)){
    out[, to:="all"]
  }
  setcolorder(out, "to")
  
  # Return [make sure same annotation info in all]
  return(out)
}))
## Save
write.table(longto, file.path(p$outdir, paste0(p$baseoutname, "_longformNsPs_to_genes.txt")),
            sep = "\t", quote = F)

# Number where there are multiple alleles in a gene that are altered to DIFF 1) AAs, 2) Codons - doesn't look like see this here

#### Annotating alleles, genes with which category they are (like in longfromNsPs counts) ####
# --- CLASSIFICATIONS:
#     no_switch_fixed: no isotype switching, totally fixed in pop
#     no_switch_variable: no isotype switching, multiple alleles in pop
#     all_switch_fixed: all are SAME isotype switch, totally fixed allele in pop 
#     all_switch_variable: all are SAME isotype switch, allelic variation in pop
#     all_switch_diff_variable: all are isotype switches but there can be DIFF isotype switches (to 2 new AAs, for example)
#     variable_switch: some alleles are isotype switches, some are not (could be multiple iso sw types)
#     NA - other, like all pseud or whatever
# --- Add classifications
# Classify each tRNA
nforclass<-alinfo[, .(AA = paste(unique(AA), collapse = ","), Codon = paste(unique(AA), collapse = ","), AlleleCM = paste(unique(AlleleCM), collapse = ","), 
                      all.lost = all.lost[1], nAls = .N, 
                      nSwNotSwClearAlleles = sum((AA!=AlleleCM | AA==AlleleCM) &  AA!="Undet" & AlleleCM!="Undet"), 
                      nSw = sum(AA!=AlleleCM & AA!="Undet" & AlleleCM!="Undet"), nNotSw = sum(AA==AlleleCM & AA!="Undet" & AlleleCM!="Undet"), # done switch manually here - so cases where there's a switch but its LOST not classified as variable_switch....
                      nTypesSw = length(unique(paste(AlleleCM, AA)[Classification!="Lost"]))),
                  by = .(displayname, tRNA)]
nforclass[all.lost==F & nAls==1 & nAls==nNotSw, gene.classif:="no_switch_fixed"]
nforclass[all.lost==F & nAls==1 & nAls==nSw, gene.classif:="all_switch_fixed"]
nforclass[all.lost==F & nAls>1 & nSwNotSwClearAlleles==nNotSw, gene.classif:="no_switch_variable"]
nforclass[all.lost==F & nAls>1 & nTypesSw==1 & nSwNotSwClearAlleles==nSw, gene.classif:="all_switch_variable"]
nforclass[all.lost==F & nAls>1 & nTypesSw > 1 & nSwNotSwClearAlleles==nSw, gene.classif:="all_switch_diff_variable"]
nforclass[all.lost==F & nSw>0 & nNotSw>0, gene.classif:="variable_switch"]

# Merge this in with allele info
setkey(nforclass, displayname, tRNA)
setkey(alinfo, displayname, tRNA)
alinfo<-alinfo[nforclass[, .(displayname, tRNA, nAls, nSw, nNotSw, nTypesSw, gene.classif)]]

# Add SECONDARY classifier about frequency - sometimes the switch is the main one, so maybe should consider it with all_switch_variable instead
ts.varnotrare<-unique(alinfo[gene.classif=="variable_switch" & switch==T & freq.ofallinclmissing > 0.05, .(displayname, tRNA)])
nforclass[ts.varnotrare, flagSwNotRare :=T]

# --- save this full annotated data in case want for other uses
write.table(nforclass, file.path(p$outdir, paste0(p$baseoutname, "_genes_swfixedvarclassifs.txt")), 
            quote = F, row.names = F, sep = "\t")
write.table(alinfo, gzfile(file.path(p$outdir, paste0(p$baseoutname, "_alleles_swfixedvarclassifs.txt.gz"))),
            quote = F, row.names = F, sep = "\t")

#### How many mutations does it take to get to each allele, each switch...etc ####
mutdir<-file.path(p$outdir, "mutpaths")
if(!dir.exists(mutdir)){dir.create(mutdir, recursive = T)}
# --- Get genes to characterize
# Have at least one switch; EXCLUDING all-pseud alleles; KEEPING cases where there are multiple AlleleCMs at gene [even tho hard to deal with]; NNN original codons
setkey(alinfo, displayname, tRNA)
alinfo[, nAlCM:=length(unique(AlleleCM[Classification!="Lost"])), by = .(displayname, tRNA)] # do NOT count pseudo ones!!
alschar<-alinfo[all.lost==F & switch==T & Codon!="NNN" & Codon.orig!="NNN"]
gschar<-alschar[, unique(tRNA), by = displayname]; setnames(gschar, "V1", "tRNA")

# Annotate whether what seems likeliest is codon mutated (usually) or backbone mutated (based on other alleles having different AlleleCMs)
alschar[, whatMut:=ifelse(nAlCM>1, "backbone", "anticodon")]

# --- Get number genes anticodon vs. backbone per gene.classif
setkey(gschar, displayname, tRNA)
setkey(alschar, displayname, tRNA)
gschar<-unique(gschar[alschar[, .(displayname, tRNA, gene.classif, whatMut)]])

setkey(gschar, displayname, gene.classif)
nanticbb<-rbind(
  # All gene classifications combined
  gschar[, .(gene.classif = "all", 
             nSw = .N,
             nAnticodonMut = sum(whatMut=="anticodon"),
             nBackboneMut = sum(whatMut=="backbone")),
         by = displayname],
  
  # by gene classification (more useful)
  gschar[, .(nSw = .N,
             nAnticodonMut = sum(whatMut=="anticodon"),
             nBackboneMut = sum(whatMut=="backbone")),
          by = .(displayname, gene.classif)]
)

write.table(nanticbb, 
            file.path(mutdir, paste0(p$baseoutname, "_ns_anticodonbackbonemutgenes.txt")), 
            quote = F, row.names = F, sep = "\t")

# --- ANTICODON mutations
# THEORETICAL
ac2aa<-fread(p$anticodontable, header = T)
# For each existing anticodon - AlleleCM pair
#     where did it come from
#     what are the codon possibilities for this
#     what mutations would have to occur to get there
antimutpaths<-chartheoranticmuts(swant = alschar[whatMut=="anticodon"], ac2aa)
gpathstemp<- antimutpaths$mpathsumm[, sum(n1MutPath>0), by = .(displayname, tRNA)]

# Overlay with ACTUAL mutation path *that happened* if evidence  
predantmuts<-realanticmuts(antimutpaths$mutpaths, varinfo) # I actually think this is working!
## Merge in other allele info of interest - especially gene classif I think
setkey(predantmuts, displayname, allelename)
setkey(alinfo, displayname, allelename)
predantmuts<-data.table(predantmuts, alinfo[predantmuts[, .(displayname, allelename)], .(freq.ofallinclmissing, gene.classif)])
setcolorder(predantmuts, c("freq.ofallinclmissing", "gene.classif"), after = 3)

# Save data
write.table(predantmuts, file.path(mutdir, paste0(p$baseoutname, "_anticodonmutations_allelepredictedpaths.txt")), 
            quote = F, row.names = F, sep = "\t" )

# Summarize
setkey(predantmuts, displayname, tRNA)
# intermediate - per gene counts [save this?]. generally counts n alleles per gene with each characteristic
predantg<-predantmuts[, .(gene.classif = gene.classif[1],
                          nAlsAnticMut = .N, nAlsVarObsInData = sum(varObsInData),
                          consistentAnticMutAcrossAlleles = consistentAnticMutAcrossAlleles[1],
                          nAlsNonAnticVars = sum(nNonAnticVariants > 0),
                          nAls1MutPath = sum("1"%in%strsplit(nMutsFromInit, ",")[[1]]),
                          nAls2MutPath = sum("2"%in%strsplit(nMutsFromInit, ",")[[1]]),
                          nAls3MutPath = sum("3"%in%strsplit(nMutsFromInit, ",")[[1]]),
                          nAlsOnePossMut = sum(nPossibleMutations==1), # one possible OR observed
                          nAlsMoreOnePossMut = sum(nPossibleMutations>1),
                          nAlsFLAG.anticodVarButNoMatch = sum(FLAG.anticodVarButNoMatch),
                          nAlsNonAnticVariants = sum(nNonAnticVariants > 0)
                          ), 
                      by = .(displayname, tRNA)]
write.table(predantg, file.path(mutdir, paste0(p$baseoutname, "_anticodonmutations_pergenesummaries.txt")), 
            quote = F, row.names = F, sep = "\t" )

## Summary to use
setkey(predantmuts, displayname, gene.classif)
setkey(predantg, displayname, gene.classif)
pantmuts.summ<-rbind(
  ## overall
  merge(
    # N alleles
    predantmuts[, .(gene.classif = "all",
                    nAllelesWAntMut = .N,
                    nAlsVarObsInDat = sum(varObsInData),
                    nAlsFLAG.anticodVarButNoMatch = sum(FLAG.anticodVarButNoMatch),
                    nAlsConsistentAnticMutAcrossAlleles = sum(consistentAnticMutAcrossAlleles),
                    nAls1MutPath = sum("1"%in%strsplit(nMutsFromInit, ",")[[1]]),
                    nAls2MutPath = sum("2"%in%strsplit(nMutsFromInit, ",")[[1]]),
                    nAls3MutPath = sum("3"%in%strsplit(nMutsFromInit, ",")[[1]]),
                    nAlsOnePossMut = sum(nPossibleMutations==1), # one possible OR observed
                    nAlsMoreOnePossMut = sum(nPossibleMutations>1),
                    nAlsNonAnticVars = sum(nNonAnticVariants > 0)
                    # same stuff for n genes...
    ),
    by = displayname],
    
    # N genes
    # per GENE stuff
    predantg[, .(
              nGenesWAntMut = .N,
                 nGeneswAllelewVarObsInDat = sum(nAlsVarObsInData > 0), # might not be ALL alleles have var Obs in data
                 nGenesAllAlleleVarObsInDat = sum(nAlsAnticMut==nAlsVarObsInData), # ALL alleles have var Obs in data
                 nGenesAnyAlleleFLAG.anticodVarButNoMatch = sum(nAlsFLAG.anticodVarButNoMatch > 0),
                 nGenesConsistentAnticMutAcrossAlleles = sum(consistentAnticMutAcrossAlleles),
                 nGenesAllAl1MutPath = sum(nAlsAnticMut==nAls1MutPath),
                 nGenesAny1MutPath = sum(nAls1MutPath > 0),
                 nGenesAllAl2MutPath = sum(nAlsAnticMut==nAls2MutPath),
                 nGenesAnyAl2MutPath = sum(nAls2MutPath > 0),
                 nGenesAllAlOnePossMut = sum(nAlsAnticMut==nAlsOnePossMut),
                 nGenesAnyAlOnePossMut = sum(nAlsOnePossMut > 0),
                 nGenesAllAlNonAnticVars = sum(nAlsNonAnticVariants==nAlsAnticMut),
                 nGenesAnyAlNonAnticVars = sum(nAlsNonAnticVariants > 0)
    ),
    by = displayname]
  )
  
  ,
  ## within gene classifications!! 
  merge(
    # N alleles
    predantmuts[, .(nAllelesWAntMut = .N,
                    nAlsVarObsInDat = sum(varObsInData),
                    nAlsFLAG.anticodVarButNoMatch = sum(FLAG.anticodVarButNoMatch),
                    nAlsConsistentAnticMutAcrossAlleles = sum(consistentAnticMutAcrossAlleles),
                    nAls1MutPath = sum("1"%in%strsplit(nMutsFromInit, ",")[[1]]),
                    nAls2MutPath = sum("2"%in%strsplit(nMutsFromInit, ",")[[1]]),
                    nAls3MutPath = sum("3"%in%strsplit(nMutsFromInit, ",")[[1]]),
                    nAlsOnePossMut = sum(nPossibleMutations==1), # one possible OR observed
                    nAlsMoreOnePossMut = sum(nPossibleMutations>1),
                    nAlsNonAnticVars = sum(nNonAnticVariants > 0)
                    # same stuff for n genes...
    ),
    by = .(displayname, gene.classif)],
    
    # N genes
    # per GENE stuff
    predantg[, .(nGenesWAntMut = .N,
                 nGeneswAllelewVarObsInDat = sum(nAlsVarObsInData > 0), # might not be ALL alleles have var Obs in data
                 nGenesAllAlleleVarObsInDat = sum(nAlsAnticMut==nAlsVarObsInData), # ALL alleles have var Obs in data
                 nGenesAnyAlleleFLAG.anticodVarButNoMatch = sum(nAlsFLAG.anticodVarButNoMatch > 0),
                 nGenesConsistentAnticMutAcrossAlleles = sum(consistentAnticMutAcrossAlleles),
                 nGenesAllAl1MutPath = sum(nAlsAnticMut==nAls1MutPath),
                 nGenesAny1MutPath = sum(nAls1MutPath > 0),
                 nGenesAllAl2MutPath = sum(nAlsAnticMut==nAls2MutPath),
                 nGenesAnyAl2MutPath = sum(nAls2MutPath > 0),
                 nGenesAllAlOnePossMut = sum(nAlsAnticMut==nAlsOnePossMut),
                 nGenesAnyAlOnePossMut = sum(nAlsOnePossMut > 0),
                 nGenesAllAlNonAnticVars = sum(nAlsNonAnticVariants==nAlsAnticMut),
                 nGenesAnyAlNonAnticVars = sum(nAlsNonAnticVariants > 0)
    ),
    by = .(displayname, gene.classif)]
  )
)

write.table(pantmuts.summ, file.path(mutdir, paste0(p$baseoutname, "_anticodonmutations_perspeciessummaries.txt")), 
            quote = F, row.names = F, sep = "\t" )

# --- BACKBONE mutations
bbvars<-backbonemuts(alsbackbone = alschar[whatMut=="backbone"], varinfo) # list of data.tables

# Save data
write.table(bbvars$vardetail,
            file.path(mutdir, paste0(p$baseoutname, "_backbonemutations_allperallele.txt")), 
            quote = F, row.names = F, sep = "\t")
write.table(bbvars$varsumm,
            file.path(mutdir, paste0(p$baseoutname, "_backbonemutations_summaryperallele.txt")), 
            quote = F, row.names = F, sep = "\t")

# Summarize
## Intermediate - get per gene
setkey(bbvars$varsumm, displayname, tRNA)
bbvarg<-bbvars$varsumm[, .(gene.classif = gene.classif[1],
                           nAlsBbMut = .N, 
                           nAls1BbMut = sum(nVars==1),
                           nAlsMultBbMut = sum(nVars > 1),
                           nAlsAcceptorVars = sum(nAcceptorVars > 0),
                           nAlsDVars = sum(nDVars > 0),
                           nAlsTVars = sum(nTVars > 0),
                           nAlsVLVars = sum(nVLVars > 0),
                           nAlsAntArmVars = sum(nAntArmVars > 0),
                           nAlsOtherStructVars = sum(nOtherStructVars > 0), 
                           nAlsNotCalledInStructVars = sum(nNotCalledInStructVars),
                           nAlsFlagNoObsVar = sum(flagNoObsVar),
                           uniqueAlsVarSubstructs = paste(unique(varSubstructs), collapse = ";")),
                       by = .(displayname, tRNA)]
write.table(bbvarg, 
            file.path(mutdir, paste0(p$baseoutname, "_backbonemutations_pergenesummaries.txt")), 
            quote = F, row.names = F, sep = "\t")

## Summary to use
setkey(bbvars$varsumm, displayname, gene.classif)
setkey(bbvarg, displayname, gene.classif)

bbvar.summ<-rbind(
  ## overall
  bbvarg[, .(
    gene.classif = "all",
    nGenesBBMut = .N,
    nGenesAllAlFlagNoObsVar = sum(nAlsFlagNoObsVar==nAlsBbMut),
    nGenesAnyAlFlagNoObsVar = sum(nAlsFlagNoObsVar > 0),
    nGenesMultBBMutAlleles = sum(nAlsBbMut>1), 
    nGenesAllAl1BbMut = sum(nAlsBbMut==nAls1BbMut),
    nGenesAnyAl1BbMut = sum(nAls1BbMut > 0),
    nGenesAllAlMultBbMut = sum(nAlsMultBbMut==nAlsBbMut),
    nGenesAnyAlMultBbMut = sum(nAlsMultBbMut>0),
    nGenesAllAlAcceptorVars = sum(nAlsAcceptorVars==nAlsBbMut),
    nGenesAnyAlAcceptorVars = sum(nAlsAcceptorVars>0),
    nGenesAllAlDVars = sum(nAlsDVars==nAlsBbMut),
    nGenesAnyAlDVars = sum(nAlsDVars>0),
    nGenesAllAlTVars = sum(nAlsTVars==nAlsBbMut),
    nGenesAnyAlTVars = sum(nAlsTVars>0),
    nGenesAllAlVLVars = sum(nAlsVLVars==nAlsBbMut),
    nGenesAnyAlVLVars = sum(nAlsVLVars>0),
    nGenesAllAlAntArmVars = sum(nAlsAntArmVars==nAlsBbMut),
    nGenesAnyAlAntArmVars = sum(nAlsAntArmVars>0),
    nGenesAllAlOtherStructVars = sum(nAlsOtherStructVars==nAlsBbMut),
    nGenesAnyAlOtherStructVars = sum(nAlsOtherStructVars>0),
    nGenesAllAlNotCalledInStructVars = sum(nAlsNotCalledInStructVars==nAlsBbMut),
    nGenesAnyAlNotCalledInStructVars = sum(nAlsNotCalledInStructVars>0)
  ), by = displayname],
  
  ## within gene classifications (more key)
  bbvarg[, .(
    nGenesBBMut = .N,
    nGenesAllAlFlagNoObsVar = sum(nAlsFlagNoObsVar==nAlsBbMut),
    nGenesAnyAlFlagNoObsVar = sum(nAlsFlagNoObsVar > 0),
    nGenesMultBBMutAlleles = sum(nAlsBbMut>1), 
    nGenesAllAl1BbMut = sum(nAlsBbMut==nAls1BbMut),
    nGenesAnyAl1BbMut = sum(nAls1BbMut > 0),
    nGenesAllAlMultBbMut = sum(nAlsMultBbMut==nAlsBbMut),
    nGenesAnyAlMultBbMut = sum(nAlsMultBbMut>0),
    nGenesAllAlAcceptorVars = sum(nAlsAcceptorVars==nAlsBbMut),
    nGenesAnyAlAcceptorVars = sum(nAlsAcceptorVars>0),
    nGenesAllAlDVars = sum(nAlsDVars==nAlsBbMut),
    nGenesAnyAlDVars = sum(nAlsDVars>0),
    nGenesAllAlTVars = sum(nAlsTVars==nAlsBbMut),
    nGenesAnyAlTVars = sum(nAlsTVars>0),
    nGenesAllAlVLVars = sum(nAlsVLVars==nAlsBbMut),
    nGenesAnyAlVLVars = sum(nAlsVLVars>0),
    nGenesAllAlAntArmVars = sum(nAlsAntArmVars==nAlsBbMut),
    nGenesAnyAlAntArmVars = sum(nAlsAntArmVars>0),
    nGenesAllAlOtherStructVars = sum(nAlsOtherStructVars==nAlsBbMut),
    nGenesAnyAlOtherStructVars = sum(nAlsOtherStructVars>0),
    nGenesAllAlNotCalledInStructVars = sum(nAlsNotCalledInStructVars==nAlsBbMut),
    nGenesAnyAlNotCalledInStructVars = sum(nAlsNotCalledInStructVars>0)
  ), by = .(displayname, gene.classif)]
)
write.table(bbvar.summ, 
            file.path(mutdir, paste0(p$baseoutname, "_backbonemutations_perspeciessummaries.txt")), 
            quote = F, row.names = F, sep = "\t")

#### Plot summaries of any of this ####
cat("....Making plots....\n")

# --- Allele frequencies of the various switching alleles (histogram, then split by various things?)
#     compare to non-switch alleles?? yep [overlap or facet]

# Set up data
alinfo[, displayname:=factor(displayname, levels = sinfo$displayname)]
tfcol<-c("blue", "red")

# Stats results
aftestres<-alinfo[all.lost==F, ksmwtests(freq.ofallinclmissing[switch==T],
                                         freq.ofallinclmissing[switch==F]), by = displayname]
## Save
write.table(aftestres, file.path(p$outdir, paste0(p$baseoutname, "_allelefreq_switchVnot_statres.txt")),
            quote = F, row.names = F, sep = "\t")
## Format to label plots with?
aftestres[, mylabel:=paste0(testshort, "~italic(p)==~", ifelse(p.value<0.01, sci2carrot(p.value), round(p.value, digits = 2)))]
aftestlab<-data.table(displayname = aftestres[testshort=="MW", displayname],
                      mylabel = paste("atop(", aftestres[testshort=="MW", mylabel], 
                            ",", aftestres[testshort=="KS", mylabel], ")"))
aftestlab[, displayname:=factor(displayname, levels = sinfo$displayname)]

# Make plots
## Counts: excludes genes with all alleles pseud
paf1<-ggplot(alinfo[all.lost==F, ], aes(freq.ofallinclmissing)) + 
  geom_histogram(aes(fill = switch), alpha = 0.4, breaks = seq(-0.005, 1.005, 0.01)) + 
  scale_fill_manual(values = tfcol) +  
  xlab("Allele frequency (n strains this allele/n total strains,\nincluding those with missing genotypes") +
  ylab("Number of alleles") + labs(fill = "Isotype switch") +
  ggtitle("Allele frequency vs. count of alleles") + 
  facet_wrap(~displayname, nrow = 3, scales = "free") + myggtheme + theme(panel.grid = element_blank()) +
  geom_label(size = 3, data = aftestlab, mapping = aes(x = 0.5, y = Inf, label = mylabel),
             parse = T, vjust = 1.0)
## proportion (of alleles in that group ideally)
  ## Many thanks to: https://github.com/donaldtmcknight/Proportional-histograms-and-density-plots-in-ggplot/blob/main/Tutorial_git.md
paf2<-ggplot(alinfo[all.lost==F, ], aes(freq.ofallinclmissing)) + 
  geom_histogram(position = "identity", aes( fill = switch, y = stat(width*density)), alpha = 0.4, breaks = seq(-0.005, 1.005, 0.01)) + 
  scale_fill_manual(values = tfcol) +  
  xlab("Allele frequency (n strains this allele/n total strains,\nincluding those with missing genotypes") +
  ylab("Proportion") + labs(fill = "Isotype switch") +
  ggtitle("Allele frequency vs prop. alleles with that frequency") + 
  facet_wrap(~displayname, nrow = 3, scales = "free") + myggtheme + theme(panel.grid = element_blank()) +
  geom_label(size = 3, data = aftestlab, mapping = aes(x = 0.5, y = Inf, label = mylabel),
             parse = T, vjust = 1.0)

# Save plots (& stat res? or just put on there)
pdf(file.path(p$outdir, paste0(p$baseoutname, "_allelefreq_switchVnot_hists.pdf")),
              6, 6.5)
print(paf1)
print(paf2)
invisible(dev.off())

# --- codon/AA on X axis, prop of various things on y [prop of genes with any altered alleles; prop of STRAINS with any altered alleles; other...?]
# Get data set up
aa.ord.from<-c("all", aas[aas%in%alinfo[,unique(AlleleCM)]])
paadat<-copy(longnps[from%in%aa.ord.from])
catcols<-c("red", "darkred", "purple", "darkblue", "blue")
names(catcols)<-alcatnames[c(1:3, 5, 4), longname]
paadat[, whatn:=factor(whatn, levels = names(catcols))]
paadat[, displayname:=factor(displayname, levels = sinfo$displayname)]

paadat[,isall:=ifelse(from=="all", "All", "Individual amino acids")]

# Number plot: separate out all...
nplt<-stackedbar(paadat, catcols, xcol = "from", stackcol = "whatn", stacknumcol = "n",
                 xorder = aa.ord.from, myxlab = "Original amino acid", myylab = "N tRNA genes",
                 legendlabel = "Allelic makeup of gene",
                 mytitle = "Number of tRNA genes with each allelic makeup", myleglabels = names(catcols)) +
  ggh4x::facet_grid2(displayname~isall, scales = "free", space = "free_x", independent = "y") +
  theme(axis.text.x = element_text(angle=45, vjust=1, hjust=1), strip.text.y = element_text(face = "italic"))

# Proportion plot
pplt<-stackedbar(paadat, catcols, xcol = "from", stackcol = "whatn", stacknumcol = "p",
                 xorder = aa.ord.from, myxlab = "Original amino acid", myylab = "Proportion this (original) amino acid's\ntRNA genes",
                 legendlabel = "Allelic makeup of gene",
                 mytitle = "Proportion of tRNA genes within each (from) amino acid", mysubt = "with each allelic makeup", myleglabels = names(catcols)) +
  ggh4x::facet_grid2(displayname~isall, scales = "free", space = "free_x") +
  theme(axis.text.x = element_text(angle=45, vjust=1, hjust=1),  strip.text.y = element_text(face = "italic")) 

# Save in same PDF
pdf(file.path(aadir, paste0(p$baseoutname, "_SwNotVarNotbyAAFrom.pdf")), 9, 8)
print(nplt)
print(pplt)
print(pplt + geom_text(data = paadat[whatn=="all switch - fixed"], aes(x = from, y = 0.1, label = totalN), color = "darkgray") +
        ggtitle( "Proportion of tRNA genes within each (from) amino acid", sub = "with each allelic makeup; total N genes annotated")) ## with total N annotated
print(pplt + 
        geom_segment(data = paadat[whatn=="no switch - fixed"], aes(x = from, y = p.low95ci, yend = p.high95ci), color = "gray") +
        ggtitle("Proportion of tRNA genes within each (from) amino acid", sub = "with each allelic makeup; first prop. 95% CI annotated")) ## with one set of CIs annotated?
invisible(dev.off())

# --- Can I plot for mutating TO something? Proportion with a given AA call that started vs. switched to that AA call?
ptoaa<-copy(longto[to%in%aa.ord.from])
ptoaa[, whatn:=factor(whatn, levels = names(catcols))]
ptoaa[, displayname:=factor(displayname, levels = sinfo$displayname)]
ptoaa[,isall:=ifelse(to=="all", "All", "Individual amino acids")]

# Number plot: separate out all...
tonplt<-stackedbar(ptoaa, catcols, xcol = "to", stackcol = "whatn", stacknumcol = "n",
                 xorder = aa.ord.from, myxlab = "Codon-matching amino acid", myylab = "N tRNA genes",
                 legendlabel = "Allelic makeup of gene",
                 mytitle = "Number of tRNA genes with each allelic makeup", myleglabels = names(catcols)) +
  ggh4x::facet_grid2(displayname~isall, scales = "free", space = "free_x", independent = "y") +
  theme(axis.text.x = element_text(angle=45, vjust=1, hjust=1), strip.text.y = element_text(face = "italic"))

# Proportion plot
topplt<-stackedbar(ptoaa, catcols, xcol = "to", stackcol = "whatn", stacknumcol = "p",
                 xorder = aa.ord.from, myxlab = "Codon-matching amino acid", myylab = "Proportion this (codon-matcing) amino acid's\ntRNA genes",
                 legendlabel = "Allelic makeup of gene",
                 mytitle = "Proportion of tRNA genes within each (to) amino acid", mysubt = "with each allelic makeup", myleglabels = names(catcols)) +
  ggh4x::facet_grid2(displayname~isall, scales = "free", space = "free_x") +
  theme(axis.text.x = element_text(angle=45, vjust=1, hjust=1),  strip.text.y = element_text(face = "italic")) 

pdf(file.path(aadir, paste0(p$baseoutname, "_SwNotVarNotbyAATo.pdf")), 9, 8)
print(tonplt)
print(topplt)
print(topplt + geom_text(data = ptoaa[whatn=="all switch - fixed"], aes(x = to, y = 0.1, label = totalN), color = "darkgray") +
      ggtitle("Proportion of tRNA genes within each (to) amino acid", sub = "with each allelic makeup; total N genes annotated"))
print(topplt + geom_segment(data = ptoaa[whatn=="no switch - fixed"], aes(x = to, y = p.low95ci, yend = p.high95ci), color = "gray") +
        ggtitle("Proportion of tRNA genes within each (to) amino acid", sub = "with each allelic makeup; first prop. 95% CI annotated"))
invisible(dev.off())

# --- Heat map or similar for from->to pairs?
# Basic to try out idea: box shaded by NUMBER of genes where this is observed (at 1 or more strains/alleles)
phdat<-copy(nsumm.aacod$aafromto[, .(displayname, from, to, nAlleles, nGenesThisSw, nGenesThisFrom.notAllPseud, nGenesThisTo.notAllPseud)])
## combine Sup/Sec
supsec<-phdat[to%in%c("Sup", "SeC")& !from=="Sup"]
setkey(supsec, displayname, from)
comb<-supsec[, .(to = "SeC", sum(nAlleles), sum(nGenesThisSw), nGenesThisFrom.notAllPseud[1], sum(nGenesThisTo.notAllPseud)), by = .(displayname, from)]
setnames(comb, names(supsec))
phdat<-rbind(phdat[!to%in%c("Sup", "SeC") & !from=="Sup"], comb)

aa.ph<-aas[!aas=="Sup"]
phdat[, displayname:=factor(displayname, levels = sinfo$displayname)]
phdat[, from:=factor(from, levels = aa.ph)]
phdat[, to:=factor(to, levels = aa.ph)]

# NUMBERS
hnum<-ggplot(phdat, aes(from, to)) + geom_tile(aes(fill = nGenesThisSw)) +
  scale_fill_gradient(low = "gray90", high = "darkblue") +
  labs(fill = "Number") +
  ggtitle("Number of genes with 1+ allele with given switch") +
  myggtheme + facet_wrap(~displayname) +
  theme(axis.text.x = element_text(angle=45, vjust=1, hjust=1), strip.text.x = element_text(face = "italic"))

# RESCALE to total N genes within each from/to separately?? [do also otherwise have this...]
phdat[, `:=`(fromProp = nGenesThisSw/nGenesThisFrom.notAllPseud,
             toProp = nGenesThisSw/nGenesThisTo.notAllPseud)]

hnum.pfrom<-ggplot(phdat, aes(from, to)) + geom_tile(aes(fill = fromProp)) +
  scale_fill_gradient(low = "gray90", high = "darkblue") +
  labs(fill = "Proportion") +
  ggtitle("Proportion of genes with 'from' amino acid \nwith 1+ allele with given switch") +
  myggtheme + facet_wrap(~displayname) +
  theme(axis.text.x = element_text(angle=45, vjust=1, hjust=1), strip.text.x = element_text(face = "italic"))
hnum.pto<-ggplot(phdat, aes(from, to)) + geom_tile(aes(fill = toProp)) +
  scale_fill_gradient(low = "gray90", high = "darkblue") +
  labs(fill = "Proportion") +
  ggtitle("Proportion of genes with 'to' amino acid \nwith 1+ allele with given switch") +
  myggtheme + facet_wrap(~displayname) +
  theme(axis.text.x = element_text(angle=45, vjust=1, hjust=1), strip.text.x = element_text(face = "italic"))

# Save
pdf(file.path(aadir, paste0(p$baseoutname, "_AAFromToHeatMaps.pdf")),
    13, 5)
print(hnum)
print(hnum.pfrom)
print(hnum.pto)
invisible(dev.off())

# --- Try a couple ways of visualizing # mutations/mutation paths/etc
# directory: mut dir

# Basic: show number of genes with switches, whether they're anticodon or backbone [this actually required no mutation info]
## [split by gene.classif]
## Get data set up
mtype<-data.table(indat = c("nAnticodonMut", "nBackboneMut"),
                  outname = c("Anticodon", "Backbone"))
long.nanticbb<-rbindlist(lapply(unique(nanticbb$displayname), function(spec){
  rbindlist(lapply(unique(nanticbb$gene.classif), function(gc){
    rbindlist(lapply(1:nrow(mtype), function(i){
      data.table(displayname = spec,
                 gene.classif = gc,
                 whatn = mtype[i, outname],
                 n = nanticbb[displayname==spec & gene.classif==gc, get(mtype[i, indat])]
      )
    }))
  }))
}))
long.nanticbb[, displayname:=factor(displayname, levels = sinfo$displayname)]
long.nanticbb[, gene.classif:=factor(gene.classif, levels = c("all", "all_switch_fixed", "all_switch_variable", "all_switch_diff_variable", "variable_switch"))]
## Cross tabulate: same thing but with where I had observed mutations that explained vs where it remains predicted*** [intersect other numbers]
## will use geom_bar_battern
obsant<-pantmuts.summ[gene.classif=="variable_switch", .(displayname, gene.classif, nObs = nGeneswAllelewVarObsInDat)] # this is number observed for anticodon stuff [more stringent]
obsbb<-bbvar.summ[gene.classif=="variable_switch", .(displayname, gene.classif, nObs = nGenesBBMut - nGenesAllAlFlagNoObsVar)] # number observed for backbone mutation stuff [less stringent]
totant<-long.nanticbb[gene.classif=="variable_switch" & whatn=="Anticodon"]
totbb<-long.nanticbb[gene.classif=="variable_switch" & whatn=="Backbone"]
setkey(obsant, displayname)
setkey(totant, displayname)
setkey(obsbb, displayname)
setkey(totbb, displayname)
typeandobs<-rbind(totant[obsant], totbb[obsbb])
setnames(typeandobs, "n", "nTotal")
typexobs<-rbindlist(lapply(unique(typeandobs$displayname), function(spec){
  rbindlist(
    lapply(c("Anticodon", "Backbone"), function(what){
      return(
        rbind(
          # Observed
          data.table(displayname = spec,
                     gene.classif = typeandobs[displayname==spec & whatn==what, gene.classif],
                     muttype = what,
                     observed_var = TRUE,
                     group = paste(what, "observed_var", sep = "_"),
                     n = typeandobs[displayname==spec & whatn==what, nObs]),
          # NOT observed
          data.table(displayname = spec,
                     gene.classif = typeandobs[displayname==spec & whatn==what, gene.classif],
                     muttype = what,
                     observed_var = FALSE,
                     group = paste(what, "not_observed_var", sep = "_"),
                     n = typeandobs[displayname==spec & whatn==what, nTotal - nObs])
        )
      ) 
  }))
}))
## SAVE this in case want it again
write.table(typexobs, file.path(mutdir, paste0(p$baseoutname, "_ns_anticodonbackbonemutgenes_obsVnot_varswitch.txt")), 
            quote = F, row.names = F, sep = "\t")

## Set up colors
mtcols<-c("purple", "goldenrod")
names(mtcols)<-mtype$outname

## Plot: species x, faceted by gene classif
p.nanticbb<-stackedbar(pinhclass = long.nanticbb, mycolors = rev(mtcols),
                       xcol = "displayname", stackcol = "whatn", stacknumcol = "n",
                       xorder = sinfo$displayname,
                       legendlabel = "Switch-underlying mutation\n(predicted)", myxlab = "Species", myylab = "Number tRNA genes",
                       myleglabels = rev(names(mtcols))) +
  facet_wrap(~gene.classif, nrow = 5, scales = "free_y") + theme(axis.text.x = element_text(face = "italic"))
## Plot: bene classif x, species facet
p.nanticbb.specfac<-stackedbar(pinhclass = long.nanticbb, mycolors = rev(mtcols),
                               xcol = "gene.classif", stackcol = "whatn", stacknumcol = "n",
                               xorder = c("all", "all_switch_fixed", "all_switch_variable", "all_switch_diff_variable", "variable_switch"),
                               legendlabel = "Switch-underlying mutation\n(predicted)", myxlab = "Gene classification", myylab = "Number tRNA genes",
                               myleglabels = rev(names(mtcols))) +
  facet_wrap(~displayname, nrow = 3, scales = "free_y") + 
  theme(strip.text = element_text(face = "italic"), axis.text.x = element_text(angle = 45, hjust = 1.05, vjust = 1))

## Plot: which ones have obs vars, only for variable_switch
typexobsnames<-c("Anticodon_observed_var", "Anticodon_not_observed_var", "Backbone_observed_var", "Backbone_not_observed_var")
typexobscols<-c("purple", "purple", "goldenrod", "goldenrod")
names(typexobscols)<-typexobsnames
typexobspats<-c("none", "stripe", "none", "stripe")
names(typexobspats)<-typexobsnames
typexobs[, group:=factor(group, levels = typexobsnames)]
p.nanticbb.obs<-ggplot(typexobs, aes(displayname, n)) +
  geom_bar_pattern(aes(fill = group, pattern = group),  stat = "identity",
                   color = "black", 
                   pattern_fill = "black",
                   pattern_angle = 45,
                   pattern_density = 0.1,
                   pattern_spacing = 0.025,
                   pattern_key_scale_factor = 0.6) + 
  scale_fill_manual(values = rev(typexobscols), breaks = rev(typexobsnames[c(1,3)]), labels = rev(mtype$outname)) + 
  scale_pattern_manual(values = typexobspats, breaks = typexobsnames[c(1, 2)], labels = c("Observed", "Not observed")) +
  xlab("") + ylab("Number of tRNA genes\nwith anticodon switches") + 
  labs(fill = "Mutation type predicted to be:", pattern = "Explanatory variant was:") +
  ggtitle("tRNA genes with some isotype switching alleles") +
  myggtheme + theme(axis.text.x = element_text(face = "italic"))

## Save plots
pdf(file.path(mutdir, paste0(p$baseoutname, "_mutTypeGeneClassifBarPlots.pdf")), 6, 8)
print(p.nanticbb)
print(p.nanticbb.specfac)
plot_grid(NULL, p.nanticbb.obs, NULL, labels = NULL, nrow = 3, rel_heights = c(1, 2, 1))
invisible(dev.off())


# For backbone muts - show where they are? Maybe? (not so useful though, across diff allelecms etc)

# # --- consider try Circos plot with lines colored by 'type' of gene?? ### IN DEV, NOT SAVED OUT
# # Initial try: circlize package chord diagram from all allele info
# lapply(1:nrow(sinfo), function(i){
#   # ***should maybe make this a function instead of doing here **
#   
#   # Set up data
#   pdat<-copy(alinfo[displayname==sinfo[i, displayname]])
#   pdat[AA=="Sup", AA:="SeC"]
#   setnames(pdat, c("AlleleCM", "AA"), c("from", "to"))    
#   ## sort by from AA and then tRNA so things get drawn in correct order. **COULD custom sort tRNAs here if want altered ones first or whatever
#   setkey(pdat, from, tRNA)
#   ## gene 'blocks' showing where genes start and end within each
#   blocks.ln<-pdat[, .N, by = .(from, tRNA)]
#   setkey(blocks.ln, tRNA)
#   blocks<-blocks.ln[, .(tRNA, end = cumsum(N)), by = from]
#   blocks[, start:=blocks[, c(1,  end[1:(.N -1)] + 1), by = from]$V1[1:nrow(blocks)]]
#   setcolorder(blocks, c("from", "tRNA", "start", "end"))
#   ## Re-order these blocks by AA before go forward
#   pdat<-pdat[c(aa.ph, "Undet")]
#   setkey(blocks, from, tRNA)
#   blocks<-blocks[c(aa.ph, "Undet")]
#   
#   # Clear previous plots
#   circos.clear()
#   
#   # Make plot
#   chordDiagram(x = pdat[, .(from, to, freq.ofallinclmissing)],
#                transparency = 0.5)  # order = c(aa.ph, "Undet"), 
#   
#   circos.trackPlotRegion(
#     ylim = c(0, 1),
#     track.height = 0.05,
#     bg.border = NA,
#     panel.fun = function(x, y) {
#       sector = CELL_META$sector.index
#       sector_blocks <- subset(blocks, from == sector)
#       # xrange <- get.cell.meta.data("xrange")
#       xlim <- CELL_META$xlim
#       for (i in 1:nrow(sector_blocks)) {
#         circos.rect(
#           xleft = xlim[1] + (xlim[2] - xlim[1]) * sector_blocks$start[i],
#           xright = xlim[1] + (xlim[2] - xlim[1]) * sector_blocks$end[i],
#           ybottom = 0,
#           ytop = 1,
#           col = NA,
#           border = "black"
#         )
#       }
#     }
#   )
#       # This block thing isn't working/aligning well yet
#   # Ala goes to only 80 but we have 99 distinct allelic connections....but WIDTH not taken into account there. Try a new way to do it
#       
#   # Give color based on other AA colors I've used (for tree??)
#   #col = XXCOLVEC[pdat$[whatever]],
#   
#   # Add color for pseud, probably: track OR lty/line border/etc
#   # somehow group together within some gene??
#   
#   # can't tel to-from if change line color....
#   # can't tell what's same gene and what's different.... 
#   # CAN add things outside to color code for other stuff
# })



#### Script completion message & session information ####
cat("....charisotypeswitching.R processing complete! Session information:....\n")
sessionInfo()