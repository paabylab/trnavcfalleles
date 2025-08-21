#! /usr/bin/env/ Rscript
# Somewhat in depth analyses of tRNA gene locations (moreso than show_mutational_variation.R) (using outputs of getstrainxtrnacalls.R, possibly others)                                      
# by Avery Davis Bell, begun 2024.12.18, added to 2025.07.09

if(length(.libPaths())==1){
  mylibloc <- .libPaths()[1]
}else{ # Presumes on PACE. This is not the best...
  mylibloc <- .libPaths()[grep(R.Version()$platform, .libPaths())]
}
require(data.table, lib.loc = mylibloc) 
require(argparser)
require(ggplot2)
require(seqinr)
require(ggnewscale)
require(hexbin)

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

freadfastseqs<-function(reftrnafa, sinfo, stripchr = T){
  # Reads in fasta files, returns in one data.table
  # In: example path - Reference genome tRNA fastas (for determining which have identical sequences).  Where species ID/species specific info is, put SAMP instead
  #         (soft link them if you don't want to F up actual name)
  #     sinfo, data.table with columns infilename (exactly how all files have this species in their name), displayname, shortname
  #     stripchr, T or F: if "chr" or "Chr" begins tRNA name (tRNAscan-SE ID), strip it out
  # Out: data.table with columns:
  #   displayname, species
  #   tRNA, tRNA seq name
  #   seq, char dna seq
  
  f.list<-lapply(1:nrow(sinfo), function(i){
    read.fasta(file = gsub("SAMP", replacement = sinfo[i, infilename], x = reftrnafa),
               seqtype = "DNA",
               as.string = T)
  })
  
  fas.dt<-rbindlist(lapply(1:nrow(sinfo), function(i){
    # want tRNAscan-SE names, not any other ones
    if(grepl("tRNAscan-SE ID", getAnnot(f.list[[i]][[1]]), fixed = T)){
      ninfo<-getAnnot(f.list[[i]])
      tnames<-sapply(sapply(sapply(ninfo, strsplit, split = "tRNAscan-SE ID: "), function(x) strsplit(x[[2]], ")")), function(x) x[1])
    }else{
      tnames<-names(f.list[[i]])
    }
    
    if(stripchr==T){
      tnames["chr"%in%substr(tnames, 1, 3) | "Chr" %in%substr(tnames, 1, 3)]<-sapply(tnames["chr"%in%substr(tnames, 1, 3) | "Chr" %in%substr(tnames, 1, 3)],
                                                                                     function(x) substr(x, 4, nchar(x)))
    }
    
    data.table(displayname = sinfo[i, displayname],
               tRNA = unlist(tnames, recursive = T),
               seq = sapply(f.list[[i]], function(x) x[1]))
  }))
  
  return(fas.dt)
}

pandci<-function(x, n){
  # Output is one-row data.table with columns p, p_low95ci, p_high95ci
  # Deals with cases where n is 0
  myp<-x/n
  if(n>0){
    bres<-binom.test(x, n)
    out<-data.table(p = myp,  p_low95ci = bres$conf.int[1], p_high95ci = bres$conf.int[2])
  }else{
    out<-data.table(p = NaN,  p_low95ci = NaN, p_high95ci = NaN)
  }
  
  return(out)
}

chrdomnums<-function(ts, allgs, chrdoms, descrip){
  # Gets numbers, proportions of genes in each chromosomal domain (and summed across domains). Also does chisq test
  # SHODDY LATER UPDATE to only have all be non tRNA (consistent with input), not separately compute that
  # In: ts, tRNA genes with location information. Required columns: displayname (species), Chr, domain (arm, center, etc), subdomain (right/left)
  #     allgs, one row per EACH gene (background -tRNAs NOT included). Req columns: displayname (species), Chr, domain (arm, center, etc), subdomain (right/left)
  #     chrdoms, chromosome domain information - used for computing gene density/kb, for example. Columns:  displayname (species), chr, domain, subdomain, start, end
  #     descrip, character description to assign to all outputs as a column (e.g., description of gene set)
  # Out: List of data tables:
  #     $persubd, counts of genes in each subdomain as well as genes-per-kb for each subdomain and geneset (all [not tRNA], tRNA)
  #     $cts.dup, counts of genes in domains summed within chromosomes and overall all chromosomes
  #     $trnaprop.subd, proportion of genes in each subdomain that are tRNAs
  #     $trnaprop.dup, proportion of genes that are tRNAs - in domains summed within chromosomes and overall across all chromosomes
  #     #all.chisq # chi-sq results for testing within this class across domains (added across chrs)
  #     $perchr.chisq = # chi-sq results for testing within this class across domains (added across chrs)
  # NB: REMOVES genes not categorized into a known bin
  
  # First, remove ones without associated domain [removes SPECIES without data too]
  ts<-ts[!is.na(domain), ]
  allgs<-allgs[!is.na(domain), ]
  
  # --- Compute numbers
  setkey(ts, displayname, Chr, domain, subdomain)
  setkey(allgs, displayname, Chr, domain, subdomain)
  # Each subdomain
  subd<-rbind(ts[, .(level = "subdomain", genes = "tRNA", .N), by = .(displayname, Chr, domain, subdomain)],
              allgs[, .(level = "subdomain", genes = "all", .N), by = .(displayname, Chr, domain, subdomain)])
  ## Add any 0s. Couldn've done counting like this but its sooo ugly
  for(sp in chrdoms[, unique(displayname)]){
    for(mychr in chrdoms[displayname==sp, unique(chr)]){
      for(mydom in chrdoms[displayname==sp & chr==mychr, unique(domain)]){
        for(mysub in chrdoms[displayname==sp & chr==mychr & domain==mydom, unique(subdomain)]){
          for(myg in c("tRNA", "all")){
            if(nrow(subd[displayname==sp & Chr==mychr & domain==mydom & subdomain==mysub & genes==myg])==0){
              subd<-rbind(subd,
                          data.table(displayname = sp, Chr = mychr, domain = mydom, 
                                     subdomain = mysub, level = "subdomain", genes = myg, N = 0)) 
          }
          }
        }
      }
    }
  }
  
  # Each domain per chromosome
  dpc<-rbind(ts[, .(subdomain = "all", level = "domain", genes = "tRNA", .N), by = .(displayname, Chr, domain)],
             allgs[, .(subdomain = "all", level = "domain", genes = "all", .N), by = .(displayname, Chr, domain)])
  ## Add any 0s. Couldn've done counting like this but its sooo ugly
  for(sp in chrdoms[, unique(displayname)]){
    for(mychr in chrdoms[displayname==sp, unique(chr)]){
      for(mydom in chrdoms[displayname==sp & chr==mychr, unique(domain)]){
          for(myg in c("tRNA", "all")){
            if(nrow(dpc[displayname==sp & Chr==mychr & domain==mydom & genes==myg])==0){
              dpc<-rbind(dpc,
                          data.table(displayname = sp, Chr = mychr, domain = mydom, subdomain = "all",
                                     level = "domain", genes = myg, N = 0)) 
            }
          }
        }
      }
  }
  
  # Total for same domain across genome. Living on the edge and assuming no 0s here.
  dall<-rbind(ts[, .(Chr = "all", subdomain = "all", level = "domain", genes = "tRNA", .N), by = .(displayname, domain)],
              allgs[, .(Chr = "all",  subdomain = "all", level = "domain", genes = "all", .N), by = .(displayname, domain)])
  
  # Add in estimated NOT tRNA number for each of these...
  # setkey(subd, displayname, Chr, domain, subdomain)
  # subd<-rbind(subd, 
  #             subd[genes=="all"][subd[genes=="tRNA"]][, .(displayname, Chr, domain, subdomain, level, genes = "all_minus_tRNA", N = N - i.N)])
  # 
  # setkey(dpc, displayname, Chr, domain)
  # dpc<-rbind(dpc,
  #            dpc[genes=="all"][dpc[genes=="tRNA"]][, .(displayname, Chr, domain, subdomain = "all", level, genes = "all_minus_tRNA", N = N - i.N)])
  # 
  # setkey(dall, displayname, domain)
  # dall<-rbind(dall, 
  #             dall[genes=="all"][dall[genes=="tRNA"]][, .(displayname, Chr = "all", domain, subdomain = "all", level, genes = "all_minus_tRNA", N = N - i.N)])
  # 
  # --- Chi-sq testing
  # Overall
  cs.all<-rbindlist(lapply(
    dall[, unique(displayname)], function(sp){
    mat<-matrix(c(dall[displayname==sp & genes=="tRNA", N],
                  dall[displayname==sp & genes=="all", N]), nrow = 2, byrow = T)
    res<-chisq.test(mat)
    out<-data.table(displayname = sp,
                    chr = "all",
                    level = "domain",
                    chisq.stat = res$statistic,
                    chisq.df = res$parameter,
                    chisq.pvalue = res$p.value)
    return(out)
  }))

  
  # Per chromosome - two ways
  cs.perc<-rbindlist(lapply(dpc[, unique(displayname)], function(sp){
    out<-rbindlist(
      lapply(dpc[displayname==sp, unique(Chr)], function(mychr){
      # Domains grouped together
      mat<-matrix(c(dpc[displayname==sp & Chr==mychr & genes=="tRNA", N],
                    dpc[displayname==sp & Chr==mychr & genes=="all", N]), nrow = 2, byrow = T)
      res.d<-chisq.test(mat)
      
      # Split subdomains
      mat<-matrix(c(subd[displayname==sp & Chr==mychr & genes=="tRNA", N],
                    subd[displayname==sp & Chr==mychr & genes=="all", N]), nrow = 2, byrow = T)
      res.s<-chisq.test(mat)
      
      # Combine & return
      out<-data.table(displayname = sp,
                      chr = mychr,
                      level = c("domain", "subdomain"),
                      chisq.stat = c(res.d$statistic, res.s$statistic),
                      chisq.df = c(res.d$parameter, res.s$parameter),
                      chisq.pvalue = c(res.d$p.value, res.s$p.value))
      return(out)
    })
    )
    return(out)
  }))
  
  # --- compute Various proportions? with binomial 95% CIs where appropriate...
  # genes per KB each region
  chrdoms[, length.kb := (end - start)/1e03]
  setkey(chrdoms, displayname, chr, domain, subdomain)
  setkey(subd, displayname, Chr, domain, subdomain)
  pkb<-chrdoms[subd]
  pkb[, N_per_kb:=N/length.kb]
  out.subd<-pkb[, .(displayname, shortname, chr, domain, subdomain, level, genes, N, N_per_kb )]  # this is to RETURN. maybe DON'T combine with other Ns.
  
  # For each region, proportion of genes there that are tRNAs (this would let you see genome info AND easily compare)
  ## All
  setkey(dall, displayname, domain)
  out.ptofall.overall<-dall[genes=="tRNA"][dall[genes=="all"]][, pandci(N, i.N), by = .(displayname, domain)]
  out.ptofall.overall[, `:=`(Chr = "all", domain = "all")]
  ## Domains per chromosome....not that useful
  setkey(dpc, displayname, Chr, domain)
  out.ptofall.domain<-dpc[genes=="tRNA"][dpc[genes=="all"]][, pandci(N, i.N), by = .(displayname, Chr, domain)]
  setcolorder(out.ptofall.overall, names(out.ptofall.domain))
  ## Subdomains (most useful for PLOTTING)
  setkey(subd, displayname, Chr, domain, subdomain)
  out.ptofall.subd<-subd[genes=="tRNA"][subd[genes=="all"]][, pandci(N, i.N), by = .(displayname, Chr, domain, subdomain)]
  
  #     ## other? Not currently

  # --- Format for return : counts
  out.subd[, description:=descrip]
  setcolorder(out.subd, "description")
  
  cts.dup<-rbind(dall, dpc) # Counts summed to one level or all levels
  cts.dup[, description:=descrip]
  setcolorder(cts.dup, "description")
  
  out.ptofall.subd[, description:=descrip]
  setcolorder(out.ptofall.subd, "description")
  
  trnaprop.dup<-rbind(out.ptofall.overall, out.ptofall.domain) # props for domains summed within chromosome or across genome
  trnaprop.dup[, description:=descrip]
  setcolorder(trnaprop.dup, "description")
  
  # --- Return
  return(list(persubd = out.subd, # counts and perkb
              cts.dup = cts.dup,
              trnaprop.subd = out.ptofall.subd,
              trnaprop.dup = trnaprop.dup,
              all.chisq = data.table(description = descrip, cs.all), # chi-sq results for testing within this class across domains (added across chrs)
              perchr.chisq = data.table(description = descrip, cs.perc) # chi-sq results for testing within this class across domains (added across chrs)
              ))
}

chitsnot<-function(all.cts.one){
  # Does chi-sq comparing tRNAs and NOT tRNAs genome wide and per chr. UPDATED to use 'all' as non-tRNA # when input updated
  # N: COUNTS of tRNA and all_minus_tRNA genes in each domain (arm, center, tip) for all chromosomes combined and across chromosomes
  # In: all.cts.one, data.table with columns domain, Chr, subdomain, level, genes (tRNA, all_minus_tRNA used), N

  out<-rbindlist(lapply(all.cts.one[, unique(Chr)], function(mychr){
    dat<-all.cts.one[Chr==mychr & level=="domain"]
    testdat<-as.matrix(c(dat[genes=="tRNA", N]), dat[genes=="all", N], nrow = 2, byrow = 2)
    res<-chisq.test(testdat)
    return(data.table(chr = mychr,
                      level = "domain",
                      chisq.stat = res$statistic,
                      chisq.df = res$parameter,
                      chisq.pvalue = res$p.value))
  }))
  
  return(out)
}

gsperdist<-function(ts, allgs, binsz, chrlens.one, descrip){
  # Gets PROPORTION of tRNAs, all genes, non-tRNAs in each chromosomal bin of provided binsz
  # UPDATED so all is the non-tRNA set, just split into tRNAs and not
  #     if last bin is < 1/2 binsz, it's combined with previous one
  # In: ts, tRNA genes with location information (for one species/set/whatever). Required columns: Chr, pos** (midpoint gene to use)
  #     allgs, one row per EACH gene (background - not just tRNAs but tRNAs presumably included). Req columns: displayname (species), Chr, pos (midpoint gene to use as position)
  #     chrlens, Chr and Length of each chr for this species
  #     descrip, character description to assign to all outputs as a column (e.g., description of gene set)
  # Out: data.table with one row per bin per set of genes (all, tRNA, all minus tRNA). Columns:
  # description, from input
  # genes, tRNA, all,
  # chr, bin chr
  # start, bin start
  # end, bin end
  # bin.mid, bin midpoint
  # n, number genes of this category in this bin
  # p, proportion of genes in this category in this bin
  # p_low95ci, lower 95% binomial CI bound on prop
  # p_high95ci, upper 95% binomial CI bound on prop
  
  # Set up bins
  mybs<-rbindlist(lapply(1:nrow(chrlens.one), function(i){
    out<-data.table(chr = chrlens.one[i, Chr],
                    start = seq(1, chrlens.one[i, Length], binsz))
    out[, end:=start + binsz - 1]
    
    # Deal with last bin
    out[ nrow(out), end:=chrlens.one[i, Length]]
    if(out[nrow(out), (end - start) < binsz/2]){
      out<-out[1:(nrow(out) - 1)]
      out[ nrow(out), end:=chrlens.one[i, Length]]
    }
    
    # Add mid point of bin for possible midpoint (rather than span) plotting
    out[, bin.mid:=start + (end - start)/2]
  }))
  
  # Overlap bins & genes, raw ns
  ts.use<-ts[, .(displayname, Chr, start = pos, end = pos) ]
  setnames(ts.use, "Chr", "chr")
  allgs.use<-allgs[, .(displayname, Chr, start = pos, end = pos)]
  setnames(allgs.use, "Chr", "chr")
  setkey(mybs, chr, start, end)
  setkey(ts.use, chr, start, end)
  setkey(allgs.use, chr, start, end)
  
  n.ts.bins<-foverlaps(mybs, ts.use)[ , .(bin.mid = bin.mid[1], n = sum(!is.na(start))), by = .(chr, i.start, i.end)]
  n.all.bins<-foverlaps(mybs, allgs.use)[ , .(bin.mid = bin.mid[1], n = sum(!is.na(start))), by = .(chr, i.start, i.end)]
  setnames(n.ts.bins, c("i.start", "i.end"), c("start", "end"))
  setnames(n.all.bins, c("i.start", "i.end"), c("start", "end"))
  
  # Get counts, proportion of genes per bin [diff gene classes, long form]
  
  ## Get props
  n.ts.bins<-data.table(n.ts.bins, 
                        rbindlist(lapply(1:nrow(n.ts.bins), function(i){
                          pandci(x = n.ts.bins[i, n], n = n.ts.bins[, sum(n)])
                        })))
  n.all.bins<-data.table(n.all.bins, 
                        rbindlist(lapply(1:nrow(n.all.bins), function(i){
                          pandci(x = n.all.bins[i, n], n = n.all.bins[, sum(n)])
                        })))
 
  ## Add descriptions, gene info
  n.ts.bins[, `:=`(description = descrip, genes = "tRNA")]
  setcolorder(n.ts.bins, c("description", "genes"))
  n.all.bins[, `:=`(description = descrip, genes = "all")]
  setcolorder(n.all.bins, c("description", "genes"))
  
  ## Combine
  out<-rbindlist(list(n.ts.bins, n.all.bins))
  setkey(out, chr, start, end)
  
  # Return
  return(out)
}

findidentgroups<-function(charvec){
  # Gets number of unique elements in charvec, how many of each there are, which original elements are identical
  # In: charvec, character vector
  # Out: data.table with one row per observed seq in charvec. Columns:
  #         seqlabel, number ID for this sequence 
  #         seq, actual sequence
  #         nobs, how many times this sequence was observed in charvec
  #         wh.obs, list column - vector of which elements of charvec match this sequence
  
  charvec.orig<-charvec
  
  charvec.un<-unique(charvec)
  
  out<-data.table(seqlabel = 1:length(charvec.un),
                  seq = charvec.un,
                  nobs = sapply(charvec.un, function(x) sum(charvec.orig==x)),
                  wh.obs = sapply(charvec.un, function(x) which(charvec.orig==x)))
}

getpairdist<-function(forpairs, diffchr = 30e06){
  # Gets all pairwise distances between 'pos' column of the provided data.table, provides TYPE of pair this was
  # In: forpairs, data.table to get distance between. Must have chromosomes Chr and pos; tRNA; AA, Codon, seqlabel, nwseqlabel, nAlleles, Lost
  #     diffchr, arbitrary number to assign as distance between pairs on different chromosome.
  # Out: data.table with one row per pair. Columns:
  #   genepair, tRNA ID of first pair member -tRNAID of second pair member[for annotating with info I get down the line] (sorted)
  #   pseud, 0, 1, or 2 - how many of the genes here had nAlleles==Lost (for downstream filtering)
  #   AA, amino acid pair (from AA column of input) [sorted]
  #   Codon, codon pair (from AA column of input) [sorted]
  #   isoacceptor, T if pair is same AA, F if not
  #   isodecoder, T if pair is same Codon, F if not
  #   sameseq, T if pair has same sequence (seqlabel), F if not
  #   lonelyseq, 0, 1, or 2 - how many of the genes are the only one of their classification that has it's seq ID (nwseqlabel = 1)
  #   samechr, T or F: was pair on same chromosome
  #   distance, distance in bp between genes
  
  out<-rbindlist(lapply(1:(nrow(forpairs) - 1), function(i){
    rbindlist(lapply((i+1):nrow(forpairs), function(j){
      checkchr = (forpairs[i, Chr]==forpairs[j, Chr])
      out<-data.table(genepair = paste(sort(unique(c(forpairs[i, tRNA], forpairs[j, tRNA]))), collapse = "-"),
                      pseud = sum(forpairs[i, nAlleles==Lost], forpairs[j, nAlleles==Lost]),
                      AA = paste(sort(c(forpairs[i, AA], forpairs[j, AA])), collapse = ","),
                      Codon = paste(sort(c(forpairs[i, Codon], forpairs[j, Codon])), collapse = ","),
                      chr = paste(sort(unique(c(forpairs[i, Chr], forpairs[j, Chr]))), collapse = ","),
                      isoacceptor = (forpairs[i, AA]==forpairs[j, AA]),
                      isodecoder = (forpairs[i, Codon]==forpairs[j, Codon]),
                      sameseq = (forpairs[i, seqlabel]==forpairs[j, seqlabel]),
                      lonelyseq = sum(forpairs[i, nwseqlabel==1],forpairs[j, nwseqlabel==1]),
                      samechr = checkchr,
                      distance = ifelse(checkchr, abs(forpairs[i, pos] - forpairs[j, pos]), diffchr))
    }))
  }))
  
  return(out)
}

pair2waystats<-function(set1, set2){
  # Does multiple stat tests comparing the distances between set1 and set2 (outputs of getpairdist split in some way of interest)
  # In: set1 and set2, data.tables to compare. Columns used are samechr and distance
  # Out: list of data.tables with stat results:
  #         $propsamechr, prop.test results comparing proportion gene pairs on same chromosome
  #         $mwdist, several Mann-Whitney/Wilcox test results: all chromosomes together but excluding pairs spanning chromosomes;
  #               pairs on each chromosome tested separately; all chromosomes together including pairs spanning chromosomes
  
  # --- Compare proportions of pairs on same vs diff chromosomes
  ptcres<-prop.test(x = c(set1[, sum(samechr)], set2[, sum(samechr)]), 
                    n = c(nrow(set1), nrow(set2)), 
                    alternative = "two.sided")
  ptcres.dt<-data.table(prop.samechr.set1 = ptcres$estimate[1],
                        prop.samechr.set2 = ptcres$estimate[2],
                        diff.prop.set2MinusSet1 = ptcres$estimate[2] - ptcres$estimate[1],
                        n.set1 = nrow(set1),
                        n.set2 = nrow(set2),
                        chisq = ptcres$statistic,
                        chisq.df = ptcres$parameter, 
                        chisq.pvalue = ptcres$p.value)
  
  # --- Compare distance distributions only for those on the same chromosome (Mann-Whitney), combining across chromosomes AND on individual chromsomes
  # All together
  mw.samechr<-wilcox.test(set1[samechr==T, distance], set2[samechr==T, distance], alternative = "two.sided")
  mw.samechr.dt<-data.table(chr = "all together, cross-chromosome comparisons excluded",
                            median.dist.set1 = set1[samechr==T, median(distance)],
                            median.dist.set2 = set2[samechr==T, median(distance)],
                            diff.median.set2MinusSet1 = set2[samechr==T, median(distance)] - set1[samechr==T, median(distance)],
                            n.set1 = nrow(set1[samechr==T]),
                            n.set2 = nrow(set2[samechr==T]),
                            W = mw.samechr$statistic,
                            mw.pvalue = mw.samechr$p.value)
  # Individual chrs
  mychrs<-sort(unique(c(set1[samechr==T, chr], set2[samechr==T, chr])))
  mw.samechr.indiv<-rbindlist(lapply(mychrs, function(thischr){
    mw.samechr<-wilcox.test(set1[samechr==T & chr==thischr, distance], set2[samechr==T & chr==thischr, distance], alternative = "two.sided")
    mw.samechr.dt<-data.table(chr = thischr,
                              median.dist.set1 = set1[samechr==T & chr==thischr, median(distance)],
                              median.dist.set2 = set2[samechr==T & chr==thischr, median(distance)],
                              diff.median.set2MinusSet1 = set2[samechr==T & chr==thischr, median(distance)] - set1[samechr==T, median(distance)],
                              n.set1 = nrow(set1[samechr==T & chr==thischr]),
                              n.set2 = nrow(set2[samechr==T & chr==thischr]),
                              W = mw.samechr$statistic,
                              mw.pvalue = mw.samechr$p.value)
    return(mw.samechr.dt)
  }))
  
  # Combined
  mw.samechr.dt<-rbind(mw.samechr.dt, mw.samechr.indiv)
  
  # --- Compare distance distributions grouping all [uses 'dummy distance' getpairdist put in for diff chrs] [interestingly this just puts median up at diff chr often]
  mw.diffchr<-wilcox.test(set1[, distance], set2[, distance], alternative = "two.sided")
  mw.diffchr.dt<-data.table(chr = "all together, cross-chromosome comparisons included",
                            median.dist.set1 = set1[, median(distance)],
                            median.dist.set2 = set2[, median(distance)],
                            diff.median.set2MinusSet1 = set2[, median(distance)] - set1[, median(distance)],
                            n.set1 = nrow(set1),
                            n.set2 = nrow(set2),
                            W = mw.diffchr$statistic,
                            mw.pvalue = mw.diffchr$p.value)
  
  # --- Format & return
  return(list(propsamechr = ptcres.dt,
              mwdist = rbind(mw.samechr.dt, mw.diffchr.dt)))
}

pair2seqdiff<-function(sp, t1, t2, fas.seqgrps, seqdiffs){
  # for a set of tRNAs in sp, pulls their info from seqdiffs *even when those genes aren't specifically in there -uses seqlabel, seqlabelpair
  # In: sp, species/display name
  #     t1, first tRNA in pair
  #     t2, second tRNA in pair
  #     fas.seqgrps, seq grp info. Columns displayname, tRNA, seqlabel are required
  #     seqdiffs, seq diff info. Row is returned (minus displayname, genepair)
  # Out: row of seqdiffs for pair matching these seq groups, minus displayname & genepair info
  
  sgp<-paste(sort(fas.seqgrps[displayname==sp & tRNA%in%c(t1, t2), seqlabel]), collapse = "-")
  orow<-seqdiffs[displayname==sp & seqlabelpair==sgp]
  if(nrow(orow)==0){ # add NAs if wasn't in data
    cl.orow<-lapply(orow, class)
    out<-as.data.table(matrix(rep(NA, ncol(orow)), nrow = 1))
    setnames(out, names(orow))
    orow<-out
    lapply(names(cl.orow), function(x){
      myfn<-paste0("as.", cl.orow[x])
      orow[, eval(x):=eval(str2expression(paste0(myfn, "(", NA, ")")))]
    })
  }
  return(orow[, .SD, .SDcols = -c(1,2)])
}

#### Arguments & inputs ####
p<-arg_parser("tRNA gene location analyses", 
              name = "trna_location_analyses.R", hide.opts = TRUE)

# tRNA Input file related
p<-add_argument(p, "--speciesf",
                help = "File containing information on all species to process here. Columns infilename (exactly how all files have this species in their name),
                displayname (name that should be used for plot outputs etc), shortname (no-spaces name for ouptut files, sorting, etc - either shorter than or same as infilename, probably).
                In order you'd like plots to be in!. **If not all files exist for each species, only does analyses that it can for each species**",
                type = "character")
p<-add_argument(p, "--trnasgen",
                help = "EXAMPLE path to *strain_trnas_gen.txt output of build_alt_sequences.py. Where species ID/species specific info is, put SAMP instead",
                type = "character")
p<-add_argument(p, "--trnainfo",
                help = "EXAMPLE path to *_bygeneup_tRNA_counts_wmissing.txt output of getstrainxtrnacalls.R. (Columns about tRNA; counts of alleles w/ various classifications.) 
                 Where species ID/species specific info is, put SAMP instead",
                type = "character")

# Other input files
p<-add_argument(p, "--genelocs",
                help = "EXAMPLE path to no-header, 3 column file specifying all PROTEIN CODING (non tRNA) genes' locations for a species.
                Columns: chr, start, end (NOT header'ed!)
                Where species ID/species specific info is, put SAMP instead")
p<-add_argument(p, "--chrdomains",
                help = "EXAMPLE path to chromosome domain info files (tip, arm, center) - rename them if you need to for this convention to work.
                Should have columns chr, domain, subdomain, start, end.
                Where species ID/species specific info is, put SAMP instead")
p<-add_argument(p, "--chrlens",
                help = "File containing lengths of chromosomes to use for plotting; chromosomes/contigs not included from those analyses/plots. Columns displayname (matching species info in speciesf), Chr, Length",
                type = "character")
p<-add_argument(p, "--reftrnafa",
                help = "Reference genome tRNA fastas (for determining which have identical sequences).  Where species ID/species specific info is, put SAMP instead
                (soft link them if you don't want to F up actual name)",
                type = "character")
p<-add_argument(p, "--alnseqdiffs",
                help = "Path to *_alignedseqdiffs_uniqueseqpairs.txt.gz output of pairwisediffsrefseq.R generated with same data as here/same seq labels etc. Nucleotide-type differences between reference tRNAs.",
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

# tRNA summary info
pertrna<-fread.mult(sinfo, exampfile = p$trnainfo)

# tRNA location info
twloc<-fread.mult(sinfo, exampfile = p$trnasgen)
if("##Chr"%in%names(twloc)){
  setnames(twloc, "##Chr", "Chr")
}
twloc[, pos:=((Start + End)/2)] # add single midpoint loc for plotting

## Make sure tRNA names match: strip 'chr' if ever there
pertrna[substr(tRNA, 1, 3)=="chr", tRNA:=sapply(tRNA, function(x) substr(x, 4, nchar(x)))]
twloc[substr(tRNA, 1, 3)=="chr", tRNA:=sapply(tRNA, function(x) substr(x, 4, nchar(x)))]

# chromosomes info
chrlens<-fread(p$chrlens, header = T)
chrdoms<-fread.mult(sinfo, exampfile = p$chrdomains)

# All gene locations
glocs<-fread.mult(sinfo, exampfile = p$genelocs, header = F)
setnames(glocs, c("V1", "V2", "V3"), c("Chr", "Start", "End"))
glocs[, pos:=((Start + End)/2)] # plotting midpoint


#### Location distribution arms/centers vs. all genes ####
cat("....Running chromosome region distribution vs. all genes analyses....\n") 
# --- Outdir
ddir<-file.path(p$outdir, "domain_analyses")
if(!dir.exists(ddir)){
  dir.create(ddir, recursive = T)
}

# --- Annotate each gene with its chromosomal region
# tRNAs
setkey(twloc, displayname, Chr, Start, End)
setkey(chrdoms, displayname, chr, start, end)
twloc<-foverlaps(twloc, chrdoms, by.x = c("displayname", "Chr", "Start", "End"),
                     by.y = c("displayname", "chr", "start", "end"))[, .(displayname, shortname, Chr, Start, End, domain, subdomain, tRNA, 
                                                                         Strand, AA, Codon, numVersions, numStrains, pos)]
# All genes
setkey(glocs, displayname, Chr, Start, End)
glocs<-foverlaps(glocs, chrdoms, by.x = c("displayname", "Chr", "Start", "End"),
                 by.y = c("displayname", "chr", "start", "end"))[, .(displayname, shortname, Chr, Start, End, domain, subdomain, pos)]


# --- Get numbers & Run chi-sq based on numbers. Per chromosome AND overall [summed across chromosomes]
# Gene sets to use
gsets.domain<-data.table(descrip = c("All tRNAs", "Functional tRNAs (at least one non-pseud allele)",
                                     "Pseudo-tRNAs (no non-pseud alleles)"),
                         inpertrna = c("Lost>=0", "nAlleles > Lost", "nAlleles==Lost"))

# Compute Ns
setkey(twloc, displayname, tRNA)
ndomain.list<-lapply(1:nrow(gsets.domain), function(i){
  gs<-pertrna[eval(parse(text = gsets.domain[i, inpertrna])), ]
  setkey(gs, displayname, tRNA)
  thisdat<-twloc[gs]
  
  out<-chrdomnums(ts = thisdat, allgs = glocs, chrdoms = chrdoms, descrip = gsets.domain[i, descrip])
  return(out)
})

# Combine across gene sets for ease of saving, plotting
all.persubd<-rbindlist(lapply(ndomain.list, function(x) x$persubd))
all.cts.dup<-rbindlist(lapply(ndomain.list, function(x) x$cts.dup))
all.trnaprop.subd<-rbindlist(lapply(ndomain.list, function(x) x$trnaprop.subd))
all.trnaprop.dup<-rbindlist(lapply(ndomain.list, function(x) x$trnaprop.dup))

all.chisqwithin<-rbindlist(lapply(ndomain.list, function(x) x$all.chisq))
all.chisqwithin.perc<-rbindlist(lapply(ndomain.list, function(x) x$perchr.chisq))

# Save data
write.table(all.persubd, file.path(ddir, paste0(p$baseoutname, "_geneCounts_subdomains_withperkb.txt")),
            sep = "\t", quote = F, row.names = F)
write.table(all.cts.dup,  file.path(ddir, paste0(p$baseoutname, "_geneCounts_domainsup_withperkb.txt")),
            sep = "\t", quote = F, row.names = F)
write.table(all.trnaprop.subd, file.path(ddir, paste0(p$baseoutname, "_tRNAGeneProportion_subdomains.txt")),
            sep = "\t", quote = F, row.names = F)
write.table(all.trnaprop.dup, file.path(ddir, paste0(p$baseoutname, "_tRNAGeneProportion_domainsup.txt")),
            sep = "\t", quote = F, row.names = F)

write.table(all.chisqwithin, file.path(ddir, paste0(p$baseoutname, "_withingeneset_regionchisq_allchrs.txt")),
            sep = "\t", quote = F, row.names = F)
write.table(all.chisqwithin.perc, file.path(ddir, paste0(p$baseoutname, "_withingeneset_regionchisq_perchr.txt")),
            sep = "\t", quote = F, row.names = F)

# --- Test gene sets against each other w/r/t proportion of genes in domains genome wide
setkey(all.cts.dup, displayname, description)
all.chisqbtwn<-all.cts.dup[, chitsnot(.SD), by = .(displayname, description)]
write.table(all.chisqbtwn, file.path(ddir, paste0(p$baseoutname, "_betweengeneset_regionchisq.txt")),
            sep = "\t", quote = F, row.names = F)

# --- data to generate - PROPORTION of each gene set in each XXkb region of genome -
nper500kb<-rbindlist(lapply(sinfo$displayname, function(sp){
  rbindlist(lapply(1:nrow(gsets.domain), function(i){
    tsdo<-pertrna[displayname==sp & eval(parse(text = gsets.domain[i, inpertrna])), tRNA]
    
    out<-data.table(displayname = sp,
               gsperdist(ts = twloc[displayname==sp & tRNA%in%tsdo],
                         allgs = glocs[displayname==sp],
                         chrlens.one = chrlens[displayname==sp],
                         binsz = 500e03,
                         descrip = gsets.domain[i, descrip])
    )
    return(out)
  }))
}))

write.table(nper500kb, file.path(p$outdir, paste0(p$baseoutname, "_genecounts_per500kb_trnasandall.txt")),
            sep = "\t", quote = F, row.names = F)

#### Make general location/domain-related plots ####
cat("Making some general location/domain-related plots....\n")
# --- Plot that shows where genes are (500kb bins) + where chromosomal domains are
# First: one line per gene set of interest (all tRNA genes, non-pseud tRNA genes, pseud tRNA genes, all not tRNA genes)
nper500kb[description=="All tRNAs" & genes=="tRNA", plot.set:="tRNAs and pseudo-tRNAs"]
nper500kb[description=="All tRNAs" & genes=="all", plot.set:="Protein-coding (?) genes"]
nper500kb[description=="Functional tRNAs (at least one non-pseud allele)" & genes=="tRNA", plot.set:="Functional tRNAs (at least one non-pseud allele)"]
nper500kb[description=="Pseudo-tRNAs (no non-pseud alleles)" & genes=="tRNA", plot.set:="Pseudo-tRNAs (no non-pseud alleles)"]
p.nper500kb<-nper500kb[!is.na(plot.set)]
# set up chr area plot
## Bounds of possible plotting area beyond each sample's bounds
maxln<-chrlens[, max(Length), by = Chr]
setnames(maxln, "V1", "max.ln")
setkey(maxln, Chr)
setkey(chrlens, Chr)
chrends<-maxln[chrlens]
setnames(chrends, "Chr", "chr")

# make plot
p500kb<-ggplot() + 
  geom_rect(data = chrends, aes(xmin = Length/1e06, xmax = max.ln/1e06, ymin = -Inf, ymax = Inf),
            fill = "gray20") + # Set up rectangles to 'block out' areas of plot that don't exist in given chromosomes
  geom_rect(data = chrdoms, aes(xmin = start/1e06, xmax = end/1e06, 
                                ymin = -Inf, ymax = Inf, fill = domain), alpha = 0.3 ) + # put chromosomal domains on this
  scale_fill_manual(values = c("arm" = "gray", "center" = "white", "tip" = "gray30")) +
  geom_line(data = p.nper500kb, aes(bin.mid/1e06, p, color = plot.set)) +
  ylab("Proportion of genes in set") + xlab("Position on chromosome (Mb)") +
  scale_x_continuous(expand = c(0, 0)) + expand_limits(y = 0) +
  facet_grid(displayname~chr, scales = "free") + myggtheme

# version with 95% CIs added as thin lines....
p500kb.ci<-p500kb + 
  geom_line(data = p.nper500kb, aes(bin.mid/1e06, p_low95ci, color = plot.set), linewidth = 0.1 ) +
  geom_line(data = p.nper500kb, aes(bin.mid/1e06, p_high95ci, color = plot.set), linewidth = 0.1) +
  ggtitle("Multiple gene sets, 95% binomial CI shown in fine lines")

# Version just comparing tRNA and not tRNA genes (easier to see)
p500kb.2genesets<-ggplot() + 
  geom_rect(data = chrends, aes(xmin = Length/1e06, xmax = max.ln/1e06, ymin = -Inf, ymax = Inf),
                              fill = "gray20") + # Set up rectangles to 'block out' areas of plot that don't exist in given chromosomes
  geom_rect(data = chrdoms, aes(xmin = start/1e06, xmax = end/1e06, 
                                ymin = -Inf, ymax = Inf, fill = domain), alpha = 0.3 ) + # put chromosomal domains on this
  scale_fill_manual(values = c("arm" = "gray", "center" = "white", "tip" = "gray30")) +
  geom_line(data = p.nper500kb[plot.set%in%c("tRNAs and pseudo-tRNAs", "Protein-coding (?) genes" )], aes(bin.mid/1e06, p, color = plot.set)) +
  ylab("Proportion of genes in set") + xlab("Position on chromosome (Mb)") +
  scale_x_continuous(expand = c(0, 0)) + expand_limits(y = 0) +
  facet_grid(displayname~chr, scales = "free") + myggtheme

p500kb.2genesets.ci<-p500kb.2genesets + 
  geom_line(data = p.nper500kb[plot.set%in%c("tRNAs and pseudo-tRNAs", "Protein-coding (?) genes" )], aes(bin.mid/1e06, p_low95ci, color = plot.set), linewidth = 0.1 ) +
  geom_line(data = p.nper500kb[plot.set%in%c("tRNAs and pseudo-tRNAs", "Protein-coding (?) genes" )], aes(bin.mid/1e06, p_high95ci, color = plot.set), linewidth = 0.1) +
  ggtitle("Two gene sets, 95% binomial CI shown in fine lines")


# SAVE
pdf(file.path(p$outdir, paste0(p$baseoutname, "_allGsvstRNA_per500kb_plots.pdf")), 12, 8)
print(p500kb)
print(p500kb.ci)
print(p500kb.2genesets)
print(p500kb.2genesets.ci)
invisible(dev.off())

# --- show genes per kb for diff categories in domains
# DATA
all.persubd[description=="All tRNAs" & genes=="tRNA", plot.set:="tRNAs and pseudo-tRNAs"]
all.persubd[description=="All tRNAs" & genes=="all", plot.set:="Protein-coding (?) genes"]

# add domain info
setkey(all.persubd, displayname, chr, domain, subdomain)
setkey(chrdoms, displayname, chr, domain, subdomain)
p.all.persubd<-chrdoms[all.persubd][!is.na(plot.set)]

# first plot: raw numbers
pperdom<-ggplot() + geom_rect(data = chrends, aes(xmin = Length/1e06, xmax = max.ln/1e06, ymin = -Inf, ymax = Inf),
                              fill = "gray20") + # Set up rectangles to 'block out' areas of plot that don't exist in given chromosomes
  geom_rect(data = chrdoms, aes(xmin = start/1e06, xmax = end/1e06, 
                                ymin = -Inf, ymax = Inf, fill = domain), alpha = 0.3 ) + # put chromosomal domains on this
  scale_fill_manual(values = c("arm" = "gray", "center" = "white", "tip" = "gray30")) +
  geom_segment(data = p.all.persubd, aes(x = start/1e06, y = N_per_kb, xend = end/1e06, yend = N_per_kb, color = plot.set)) +
  ylab("Genes per kb") + xlab("Position on chromosome (Mb)") +
  scale_x_continuous(expand = c(0, 0)) + expand_limits(y = 0) +
  facet_grid(displayname~chr, scales = "free") + myggtheme
  
# Next: normalize per-kb rate to total n genes: divided by median
setkey(p.all.persubd, displayname, plot.set)
p.all.persubd[, perkb_norm:=N_per_kb/median(N) * 100, by = .(displayname, plot.set)]

pperdom.norm<-ggplot() + geom_rect(data = chrends, aes(xmin = Length/1e06, xmax = max.ln/1e06, ymin = -Inf, ymax = Inf),
                                   fill = "gray20") + # Set up rectangles to 'block out' areas of plot that don't exist in given chromosomes
  geom_rect(data = chrdoms, aes(xmin = start/1e06, xmax = end/1e06, 
                                ymin = -Inf, ymax = Inf, fill = domain), alpha = 0.3 ) + # put chromosomal domains on this
  scale_fill_manual(values = c("arm" = "gray", "center" = "white", "tip" = "gray30")) +
  geom_segment(data = p.all.persubd, aes(x = start/1e06, y = perkb_norm, xend = end/1e06, yend = perkb_norm, color = plot.set)) +
  ylab("Genes density\nnormalized to number of genes") + xlab("Position on chromosome (Mb)") +
  scale_x_continuous(expand = c(0, 0)) + expand_limits(y = 0) +
  facet_grid(displayname~chr, scales = "free") + myggtheme

# --- within each chromosome, what proportion of genes of that class on that chromosome are in each domain?
# data setup
setkey(p.all.persubd, displayname, chr, plot.set)
p.all.persubd[,pOnChr:=N/sum(N) , by = .(displayname, chr, plot.set)]

# Make plot
pperdom.winchr<-ggplot() + geom_rect(data = chrends, aes(xmin = Length/1e06, xmax = max.ln/1e06, ymin = -Inf, ymax = Inf),
                                     fill = "gray20") + # Set up rectangles to 'block out' areas of plot that don't exist in given chromosomes
  geom_rect(data = chrdoms, aes(xmin = start/1e06, xmax = end/1e06, 
                                ymin = -Inf, ymax = Inf, fill = domain), alpha = 0.3 ) +
  scale_fill_manual(values = c("arm" = "gray", "center" = "white", "tip" = "gray30")) +
  new_scale_fill() +
  geom_rect(data = p.all.persubd, aes(xmin = start/1e06, ymin = -Inf, xmax = end/1e06, ymax = pOnChr, fill = plot.set), alpha = 0.6) +
  scale_fill_manual(values = c("Protein-coding (?) genes" = "blue", "tRNAs and pseudo-tRNAs" = "red")) +
  ylab("Proportion of genes\non each chromosome") + xlab("Position on chromosome (Mb)") +
  scale_x_continuous(expand = c(0, 0)) + expand_limits(y = 0) +
  facet_grid(displayname~chr, scales = "free") + myggtheme

# Save
pdf(file.path(ddir, paste0(p$baseoutname, "_genestrnas_perkbprops_indomains.pdf")), 12, 8)
print(pperdom)
print(pperdom.norm)
print(pperdom.winchr)
invisible(dev.off())

# ....still should come up with better way to show this...
#         count number of times middle has more, fewer than arms (or one arm)...?

#### Investigate distance between genes ####
# --- Read in FASTA files, determine which sequences are identical
# Read
fas<-freadfastseqs(p$reftrnafa, sinfo, stripchr = T)

# Figure out identical groups
setkey(fas, displayname)
idgs<-fas[, findidentgroups(seq), by = displayname]
fas.seqgrps<-rbindlist(lapply(1:nrow(idgs), function(i){
  dat<-fas[displayname==idgs[i, displayname],]
  dat[idgs[i, wh.obs][[1]], seqlabel:=idgs[i, seqlabel]]
  out<-dat[!is.na(seqlabel)]
  out[, nwseqlabel:=nrow(out)]
  return(out)
}))

## Save this
write.table(fas.seqgrps, file.path(p$outdir, paste0(p$baseoutname, "_identicalseqgroupinfo.txt")),
            quote = F, row.names = F, sep = "\t")

# Annotate tRNA location info with identical info
setkey(twloc, displayname, tRNA)
setkey(pertrna, displayname, tRNA)
setkey(fas.seqgrps, displayname, tRNA)
tallinfo<-fas.seqgrps[pertrna][twloc]
tallinfo[, `:=`(seq=NULL, shortname = NULL, i.shortname = NULL, AA=NULL, Codon=NULL, numVersions = NULL)] # clean up - remove ugly & dup columns
setnames(tallinfo, c("i.AA", "i.Codon"), c("AA", "Codon") )
setcolorder(tallinfo, c("displayname", "tRNA", "seqlabel", "nwseqlabel", "AA", "Codon", "Chr", "Start", "End", "Strand", "pos", "domain", "subdomain", 
                         "VariableInPop", "nAlleles", "numStrains", "anyMissingCalls", "Best", "Functional", "Altered", "Lost"))
## Save this - good to have. REFERENCE info.
write.table(tallinfo, file.path(p$outdir, paste0(p$baseoutname, "_tRNAgeneinfo_combined.txt")),
            quote = F, row.names = F, sep = "\t")

## Save little SUMMARY of this. Useful!
setkey(tallinfo, displayname)
gsumm<-tallinfo[,.(nGenesInclPseud = .N,
                   nGenesWRefFastaSeqs = sum(!is.na(seqlabel)),
                   nUniqueGSeqs = length(unique(seqlabel, na.rm = T)),
                   nGenesNotAllPseud = sum(nAlleles > Lost),
                   nGeneswRefFastaSeqs.notAllPseud = sum(nAlleles>Lost & !is.na(seqlabel)),
                   nUniqueGSeqs.notAllPseud = length(unique(seqlabel[nAlleles>Lost])),
                   nGenesNoPseudAls = sum(Lost==0),
                   nGenesWRefFastaSeqs.noPseudAls = sum(!is.na(seqlabel) & Lost==0),
                   nUniqueGSeqs.noPseudAls = length(unique(seqlabel[Lost==0]))
                   ),
                by = displayname]
write.table(gsumm, file.path(p$outdir, paste0(p$baseoutname, "_tRNAgenecounts_refsequniqueness.txt")),
            quote = F, row.names = F, sep = "\t")

# --- Get among-gene distances
# Directory
pdir<-file.path(p$outdir, "pairwisedistance")
if(!dir.exists(pdir)){dir.create(pdir, recursive = T)}

# Get data: for ALL pairs, with classifications including if any in pair were pseud (all pseud alleles) included
cat("...Getting all pairwise position differences between tRNAs, this takes a bit of time.....\n") 
setkey(tallinfo, displayname)
pairdists<-tallinfo[, getpairdist(.SD, diffchr = 30e06), by = displayname]

# --- Intersect with differences between sequences [alignments]
seqdiffs<-fread(p$alnseqdiffs, header = T)
setkey(seqdiffs, displayname, genepair)
setkey(pairdists, displayname, genepair)
pairdists<-seqdiffs[pairdists] # KEEP ALL here even though will only have these new numbers for genes that I could do sec structure assignment on

#     can tell from info here if seqs are identical - that's already in the pair thing. If so, assign 0 to pairwise differences
pairdists[sameseq==T & is.na(seqlabelpair), `:=`(nNoMatch=0, nMatch=pairdists[1, nNoMatch + nMatch], stretchNoMatch = 0)] # add in 0s here....


# Add in 'proper' seq differences for the (many) pairs where they're not the 'unique' chosen one: do based on seq label group
toaddsdiff<-pairdists[is.na(seqlabelpair) & sameseq==F, pair2seqdiff(sp = displayname, t1 = strsplit(genepair, "-")[[1]][1],
                                                                     t2 = strsplit(genepair, "-")[[1]][2],
                                                                     fas.seqgrps, seqdiffs), 
                      by = .(displayname, genepair)] # this is SLOW; if works, need to save it
pairdists<-rbind(pairdists[!is.na(seqlabelpair) | sameseq==T],
                  toaddsdiff[pairdists[is.na(seqlabelpair) & sameseq==F, .(displayname, genepair, pseud, AA, Codon, chr, isoacceptor, isodecoder, sameseq, lonelyseq, samechr, distance)]])

## Save this
write.table(pairdists, gzfile(file.path(pdir, paste0(p$baseoutname, "_tRNApairdistances.txt.gz"))),
            quote = F, row.names = F, sep = "\t")


# ***For number summ below, add in with/without a 'seq difference' number [0 or otherwise]...total, with, without...
#       intersect with intron stuff too, probably....

# --- Get NUMBER of pairs each class type...
## Set up classes
pdgclasses<-data.table(descrip = c("Different isoacceptor", "Same isoacceptor", "Same isoacceptor, different isodecoder",
                                   "Same isodecoder", "Same isodecoder, different sequence", "Same sequence"),
                       toeval = c("isoacceptor==F", "isoacceptor==T", "isoacceptor==T & isodecoder==F", 
                                  "isodecoder==T", "isodecoder==T & sameseq==F", "sameseq==T"))
  # diff isoacceptor necessarily is also different isodecoder, sequence
sdclasses<-data.table(descrip = c("Any seq aln difference, including NA or none", "Measurable seq aln difference (not NA, > 0)",
                                  "Seq alignments are identical [except perhaps intron] (not NA, ==0)"),
                      toeval = c("is.na(nNoMatch) | !is.na(nNoMatch)", "!is.na(nNoMatch) & nNoMatch>0", "!is.na(nNoMatch) & nNoMatch==0"))

## N summary  
gsets.pair<-data.table(descrip = c("All tRNAs (pseud & not)", "Functional tRNAs (at least one non-pseud allele in both pair members)"),
                       toeval = c("pseud>=0", "pseud==0"))

npairsumm<-rbindlist(lapply(1:nrow(gsets.pair), function(i){
  rbindlist(lapply(1:nrow(sdclasses), function(z){
    rbindlist(lapply(1:nrow(pdgclasses), function(j){
      pairdists[, .(Genes = gsets.pair[i, descrip],
                    pairclass = pdgclasses[j, descrip],
                    seqdiffclass = sdclasses[z, descrip], 
                    N = sum(eval(parse(text = pdgclasses[j, toeval])) & eval(parse(text = gsets.pair[i, toeval])) &
                      eval(parse(text = sdclasses[z, toeval])))),
                by = displayname]
    }))
  }))
}))
### Save this
write.table(npairsumm, file.path(pdir, paste0(p$baseoutname, "_tRNApairtypes_nsumm.txt")),
            quote = F, row.names = F, sep = "\t")
    

# --- Stats for appropriate ones of these
# Set up all 2-way comparisons to do
testinfo<-data.table(set1name = c("Different isoacceptor", "Same isoacceptor, different isodecoder", "Same isodecoder, different sequence", "Different sequence"),
                     set2name = c("Same isoacceptor", "Same isodecoder", "Same sequence", "Same sequence"),
                     set1eval = c("isoacceptor==F", "isoacceptor==T & isodecoder==F", "isodecoder==T & sameseq==F", "sameseq==F"),
                     set2eval = c("isoacceptor==T", "isodecoder==T", "sameseq==T", "sameseq==T"))

# Run for this and excluding pseud vs. not
testpds<-lapply(sinfo$displayname, function(sp){
  lapply(1:nrow(gsets.pair), function(i){
    lapply(1:nrow(testinfo), function(j){
      testres<-pair2waystats(set1 = pairdists[displayname==sp & eval(parse(text = gsets.pair[i, toeval])) & eval(parse(text = testinfo[j, set1eval]))],
                             set2 =  pairdists[displayname==sp & eval(parse(text = gsets.pair[i, toeval])) & eval(parse(text = testinfo[j, set2eval]))])
      mwdist<-data.table(displayname = sp,
                Genes = gsets.pair[i, descrip],
                 Set1 = testinfo[j, set1name],
                 Set2 = testinfo[j, set2name],
                testres$mwdist)
      propsamechr<-data.table(displayname = sp,
                              Genes = gsets.pair[i, descrip],
                              Set1 = testinfo[j, set1name],
                              Set2 = testinfo[j, set2name],
                              testres$propsamechr)
      return(list(mwdist = mwdist,
                  propsamechr = propsamechr))
    })
  })
})

## Combine across gsets
mwdist<-rbindlist(lapply(testpds, function(x) rbindlist(lapply(x, function(y) rbindlist(lapply(y, function(z) z$mwdist))))))
ptestsamechr<-rbindlist(lapply(testpds, function(x) rbindlist(lapply(x, function(y) rbindlist(lapply(y, function(z) z$propsamechr))))))

## Save
write.table(ptestsamechr, file.path(pdir, paste0(p$baseoutname, "_tRNApairdistance_propsamechr_stats.txt")),
            quote = F, row.names = F, sep = "\t")
write.table(mwdist, file.path(pdir, paste0(p$baseoutname, "_tRNApairdistance_MWcompare_stats.txt")),
            quote = F, row.names = F, sep = "\t")

# **include & exclude pseud for all of these? OR just exlc pseud? Somewhere above, automate doing this +/- pseud-incl pairs...

# --- GLM for ones on same chr
# Set up models to do: NEW
mymodinfo<-data.table(seqdiffdescrip = rep(c("Number of alignment bases/gaps that differ", 
                                             "Number of distinct stretches of alignment bases/gaps that differ"), each = 6), 
                   seqdiffmeasure = rep(c("nNoMatch", "stretchNoMatch"), each = 6), # ...do for stretches AND counts?
                   description = rep(c("Null: distance ~ 1 [no all-pseud genes]", "Chr. only: distance ~ chromosome [no all-pseud genes]",
                                       "Seq diffs, chromosome: distance ~ chromosome + seq. diff. [no all-pseud genes]",
                                       "Null: distance ~ 1 [with all-pseud genes]", "Chr. and pseud only: distance ~ chromosome + pseud",
                                       "Seq diffs, chrom, pseud: distance ~ chromosome + pseud + seq. diff"), 2),
                   data.slice= rep(c("samechr==T & pseud==0", "samechr==T & pseud==0", "samechr==T & pseud==0",
                                     "samechr==T", "samechr==T", "samechr==T"), 2),
                   use.formula = c("distance ~ 1", "distance ~ chr", "distance ~ chr + DIFF", "distance ~ 1", "distance ~ chr + pseud",
                                   "distance ~ chr + pseud + DIFF"))
tofix<-mymodinfo[, grep("DIFF", use.formula)]
for(i in tofix){
  mymodinfo[i, use.formula:=gsub("DIFF", seqdiffmeasure, use.formula)]
}

# Do the bunch of models PER SPECIES
mods<-lapply(sinfo$displayname, function(sp){
  out<-lapply(1:nrow(mymodinfo), function(i){
    out<-glm(formula(mymodinfo[i, use.formula]),
               data = pairdists[displayname==sp & eval(parse(text = mymodinfo[i, data.slice])), ], 
               family = Gamma(link = "log"))
    return(out)
  })
  names(out)<-sapply(1:nrow(mymodinfo), function(i) paste(unlist(mymodinfo[i, ]), collapse = "; "))
  return(out)
})
names(mods)<-sinfo$displayname
## SAVE just in case need later
save(mods, file = file.path(pdir, paste0(p$baseoutname, "_tRNApairdistance_seqdiffs_glms_gammalog.RData")))

# Summarize things that are ONE PER model: AICs, dispersions, seq diff beta & how much change per seq diff [may be same, multiplicative? check if need un-log etc]
modsumm<-rbindlist(lapply(sinfo$displayname, function(sp){
  out<-data.table(displayname = sp, mymodinfo)
  out[, `:=`(model.AIC = sapply(mods[[sp]], function(x) x$aic),
              model.dispersion = sapply(mods[[sp]], function(x){
                sum(residuals(x, type = "deviance")^2) / x$df.residual #  value near 1 suggests good fit. # Much greater than 1 → overdispersion (consider alternative models or robust SEs).
              })
  )]
  ## Add stuff just for ones with seq diff measure included. ASSUMES IT IS LAST, IF NOT THIS WILL BREAK*
  out[tofix, `:=`(seqdiff.beta = sapply(mods[[sp]][tofix], function(x) summary(x)$coefficients[length(x$coefficients), "Estimate"]),
                        seqdiff.beta.se = sapply(mods[[sp]][tofix], function(x) summary(x)$coefficients[length(x$coefficients), "Std. Error"]),
                        seqdiff.tval = sapply(mods[[sp]][tofix], function(x) summary(x)$coefficients[length(x$coefficients), "t value"]),
                        seqdiff.pval =  sapply(mods[[sp]][tofix], function(x) summary(x)$coefficients[length(x$coefficients), "Pr(>|t|)"]),
                        seqdiff.foldchange = exp(sapply(mods[[sp]][tofix], function(x) summary(x)$coefficients[length(x$coefficients), "Estimate"])),
                        seqdiff.percentchange = paste0(round((exp(sapply(mods[[sp]][tofix], function(x) summary(x)$coefficients[length(x$coefficients), "Estimate"]))-1)*100, 2), "%"))]
}))

#     betas are in log. exp(coef(model))  # Multiplicative effect - e.g. 2 equals doubling

write.table(modsumm, file.path(pdir, paste0(p$baseoutname, "_tRNApairdistance_seqdiffs_glmsummary_gammalog.txt")),
            sep = "\t", quote = F, row.names = F)

# Plots
    # Residuals vs Fitted: Look for patterns (should be random scatter).
    # Normal Q-Q: Deviance residuals should roughly follow a straight line.
    # Scale-Location: Checks homoscedasticity.
    # Residuals vs Leverage: Identifies influential points.
lapply(1:nrow(sinfo), function(si){
  pdf(file.path(pdir, paste0(p$baseoutname, "_tRNApairdistance_seqdiffs_glm_gammalog_plots_", sinfo[si, shortname], ".pdf")),
      12, 12)
  par(mfrow = c(2, 2))
  lapply(1:nrow(mymodinfo), function(i){
    plot(mods[[si]][[i]], main = mymodinfo[i, description],
         sub = mymodinfo[i, seqdiffmeasure])
  })
  par(mfrow = c(1, 1))
  invisible(dev.off())
})


# Model comparisons
## Set up
mymodcomps<-data.table(description = c("With seq diff (nNoMatch) vs without (chrom only; no all-pseud genes)",
                                       "With seq diff (stretchNoMatch) vs without (chrom only; no all-pseud genes)",
                                       "With seq diff (nNoMatch) vs without (chrom + pseud info)",
                                       "With seq diff (stretchNoMatch) vs without (chrom + pseud info)",
                                       "Chrom. only vs null (no all pseud genes)",
                                       "Chrom + pseud only vs null"),
                       m1 = c(2, 8, 5, 11, 1, 4), # null-er
                       m2 =  c(3, 9, 6, 12, 2, 5)# test-er
                         )
## McFadden's R2s * for seq diff *
mcr2<-rbindlist(lapply(names(mods), function(sp){
  rbindlist(lapply(1:nrow(mymodcomps), function(i){
    data.table(displayname = sp,
               model.comparison = mymodcomps[i, description],
               loglik.nuller = logLik(mods[[sp]][[mymodcomps[i, m1]]]),
               loglik.tester = logLik(mods[[sp]][[mymodcomps[i, m2]]]),
               mcfaddenR2 = 1 - (logLik(mods[[sp]][[mymodcomps[i, m2]]])/logLik(mods[[sp]][[mymodcomps[i, m1]]])))
  }))
}))
write.table(mcr2, file.path(pdir, paste0(p$baseoutname, "_tRNApairdistance_seqdiffs_glmgammalog_mcfaddensr2.txt")),
            sep = "\t", quote = F, row.names = F)

## ANOVAs
aovs<-lapply(mods, function(x){
  out<-lapply(tofix, function(i){
    anova(x[[i]], test = "Chi")
  })
  names(out)<-sapply(tofix, function(i) paste(unlist(mymodinfo[i, ]), collapse = "; "))
  return(out)
})
names(aovs)<-names(mods)
  # save just in case
save(aovs, file = file.path(pdir, paste0(p$baseoutname, "_tRNApairdistance_seqdiffs_glms_gammalog_anovas.RData")))

## Format to save out relevant info 
aov.summ<-rbindlist(lapply(names(aovs), function(sp){
  rbindlist(lapply(1:length(tofix), function(i){
    oav<-aovs[[sp]][[i]]
    out<-data.table(displayname = sp,
                    mymodinfo[tofix[i],.(seqdiffdescrip, description)],
                    as.data.table(oav, keep.rownames = T))
    setnames(out, "rn", "model_term")
    pdev<-oav$Deviance/sum(oav$Deviance, na.rm = T) # compare just model terms
    ptotdev<-oav$Deviance/sum(oav$`Resid. Dev`[1])
    out[,`:=`(prop.modeltermdeviance = pdev,
              percent.modeltermdeviance = round(pdev*100, 2),
              prop.totaldeviance = ptotdev,
              percent.totaldeviance = round(ptotdev*100,2))]
    return(out)
  }))
}))

write.table(aov.summ, file.path(pdir, paste0(p$baseoutname, "_tRNApairdistance_seqdiffs_glmgamma_anovapdev.txt")),
            sep = "\t", quote = F, row.names = F)


# TO DO
## OOOH MIGHT BE NICE TO ADD IN IF KNOWING THEY'RE THE SAME ACCEPTOR OR DECODER OR WHATEVER ACTUALLY ADDS INFO ON TOP OF THIS?? ##
##    pseud can come in AFTER
##    think more about doing w/w/o chr...
  
# --- PLOTS, BABY
# Proportion on same chromosome - different categories
## Set up data so all props are together
pdat.samechr<-rbind(ptestsamechr[, .(displayname, Genes, GeneSubset = Set1, prop = prop.samechr.set1, n = n.set1)],
                    ptestsamechr[Set1!="Different sequence", .(displayname, Genes, GeneSubset = Set2, prop = prop.samechr.set2, n = n.set2)])
set.ord<-c("Different isoacceptor", "Different sequence", "Same isoacceptor", "Same isoacceptor, different isodecoder", 
           "Same isodecoder", "Same isodecoder, different sequence", "Same sequence")
pdat.samechr[, `:=`(displayname = factor(displayname, levels= sinfo$displayname), GeneSubset = factor(GeneSubset, levels = set.ord))]
## Make plot
psamechr<-ggplot(pdat.samechr, aes(GeneSubset, prop)) + geom_bar(stat = "identity") +
  xlab("Gene subset") + ylab("Proportion of genes of this type that are on the same chromosome") +
  myggtheme + theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust = 1)) + 
  facet_grid(Genes~displayname)

pdf(file.path(pdir, paste0(p$baseoutname, "_tRNApairdistance_propsamechr_plot.pdf")), 12, 8)
print(psamechr)
invisible(dev.off())

# All pairwise tests done - distance distributions. EXCL PSEUDS.
pdisthists<-lapply(1:nrow(testinfo), function(i){
  # Set up data
  oned<-pairdists[pseud==0 & eval(parse(text = testinfo[i, set1eval])) | eval(parse(text = testinfo[i, set2eval])), ]
  oned[eval(parse(text = testinfo[i, set1eval])), classif:=testinfo[i, paste(set1name)]]
  oned[eval(parse(text = testinfo[i, set2eval])), classif:=testinfo[i, paste(set2name)]]
  
  mytitle<-testinfo[i, paste(set1name, "vs.", set2name)]
  
  # Make plot including those on diff chrs
  bnsz <- 0.5
  plt<-ggplot(mapping = aes(distance/1e06)) + geom_histogram(data = oned[classif==testinfo[i, paste(set1name)]],
                                                   aes(y = after_stat(density)*bnsz, fill = classif), alpha = 0.6, binwidth = bnsz) +
    geom_histogram(data = oned[classif==testinfo[i, paste(set2name)]],
                   aes(y = after_stat(density)*bnsz, fill = classif), alpha = 0.6, binwidth = bnsz) +
    scale_fill_manual(values = c(rgb(0,0,1,0.5), rgb(1,0,0,0.5))) +
    ylab("Proportion of gene pairs") + xlab("Distance between genes in pair (Mb)") +
    labs(fill = "Gene set") + ggtitle(mytitle, subtitle = "Includes genes on diff. chrs (30Mb); \nexcludes genes with all pseud. alleles") +
    myggtheme +
    facet_wrap(~displayname, nrow = 3)
  
  # Make plot excluding those on diff chrs
  bnsz<-0.25
  plt.samechr<-ggplot(mapping = aes(distance/1e06)) + geom_histogram(data = oned[samechr==T & classif==testinfo[i, paste(set1name)]],
                                                                     aes(y = after_stat(density)*bnsz, fill = classif), alpha = 0.6, binwidth = bnsz) +
    geom_histogram(data = oned[samechr==T & classif==testinfo[i, paste(set2name)]],
                   aes(y = after_stat(density)*bnsz, fill = classif), alpha = 0.6, binwidth = bnsz) +
    scale_fill_manual(values = c(rgb(0,0,1,0.5), rgb(1,0,0,0.5))) +
    ylab("Proportion of gene pairs") + xlab("Distance between genes in pair (Mb)") +
    labs(fill = "Gene set") + ggtitle(mytitle, subtitle = "For all tRNAs on same chr;\nexcludes genes with all pseud. alleles") +
    myggtheme +
    facet_wrap(~displayname, nrow = 3)
  
  # Return plots
  return(list(allplt = plt, samechrplt = plt.samechr))
})

pdf(file.path(pdir, paste0(p$baseoutname, "_tRNApairdistance_2waydistcompare_plot.pdf")), 6, 8)
lapply(pdisthists, function(x){
  print(x$allplt)
  print(x$samechrplt)
})
invisible(dev.off())
# NEED TO DOCUMENT THIS, FIGURE OUT WHY DUPLICATING


# Plot vs SEQ divergence!! [before do model in case that informs]
## Number stretches
dvdivp<-ggplot(pairdists[pseud==0 & samechr==T], aes(stretchNoMatch, distance/1e06)) + geom_hex() + geom_smooth(color = "red") +
  myggtheme + facet_wrap(~displayname) + ggtitle("Pairs with pseuds excluded") + ylab("Distance (Mb)") +
  xlab("Seq divergence - n non-matching stretches")
dvdivp.log<-ggplot(pairdists[pseud==0 & samechr==T], aes(stretchNoMatch, distance)) + geom_hex() + geom_smooth(color = "red") +
  myggtheme + facet_wrap(~displayname) + ggtitle("Pairs with pseuds excluded, log scale y") + ylab("Distance (bp)") + 
  xlab("Seq divergence - n non-matching stretches") +
  scale_y_log10()  # x logging doesn't make sense: lots of zeros
### plus chr facets
dvdivp.chr<-ggplot(pairdists[pseud==0 & samechr==T], aes(stretchNoMatch, distance/1e06)) + geom_hex() + geom_smooth(color = "red") +
  myggtheme + facet_grid(displayname~chr) + ggtitle("Pairs with pseuds excluded") + ylab("Distance (Mb)") + 
  xlab("Seq divergence - n non-matching stretches")
dvdivp.log.chr<-ggplot(pairdists[pseud==0 & samechr==T], aes(stretchNoMatch, distance)) + geom_hex() + geom_smooth(color = "red") +
  myggtheme + facet_grid(displayname~chr) + ggtitle("Pairs with pseuds excluded, log scale y") + ylab("Distance (bp)") + 
  xlab("Seq divergence - n non-matching stretches") +
  scale_y_log10()  # x logging doesn't make sense: lots of zeros

## Number nt diffs
dvdivp_n<-ggplot(pairdists[pseud==0 & samechr==T], aes(nNoMatch/nMatch, distance/1e06)) + geom_hex() + geom_smooth(color = "red") +
  myggtheme + facet_wrap(~displayname) + ggtitle("Pairs with pseuds excluded") + ylab("Distance (Mb)") +
  xlab("Seq divergence - prop nt diffs")
dvdivp.log_n<-ggplot(pairdists[pseud==0 & samechr==T], aes(nNoMatch/nMatch, distance)) + geom_hex() + geom_smooth(color = "red") +
  myggtheme + facet_wrap(~displayname) + ggtitle("Pairs with pseuds excluded, log scale y") + ylab("Distance (bp)") + 
  xlab("Seq divergence - prop nt diffs") +
  scale_y_log10()  # x logging doesn't make sense: lots of zeros
### plus chr facets
dvdivp.chr_n<-ggplot(pairdists[pseud==0 & samechr==T], aes(nNoMatch/nMatch, distance/1e06)) + geom_hex() + geom_smooth(color = "red") +
  myggtheme + facet_grid(displayname~chr) + ggtitle("Pairs with pseuds excluded") + ylab("Distance (Mb)") + 
  xlab("Seq divergence - prop nt diffs")
dvdivp.log.chr_n<-ggplot(pairdists[pseud==0 & samechr==T], aes(nNoMatch/nMatch, distance)) + geom_hex() + geom_smooth(color = "red") +
  myggtheme + facet_grid(displayname~chr) + ggtitle("Pairs with pseuds excluded, log scale y") + ylab("Distance (bp)") + 
  xlab("Seq divergence - prop nt diffs") +
  scale_y_log10()  # x logging doesn't make sense: lots of zeros

## Save
pdf(file.path(pdir, paste0(p$baseoutname, "_tRNApairdistance_vsseqdiv.pdf")), 10, 6)
print(dvdivp)
print(dvdivp.chr)
print(dvdivp.log)
print(dvdivp.log.chr)
print(dvdivp_n)
print(dvdivp.chr_n)
print(dvdivp.log_n)
print(dvdivp.log.chr_n)
invisible(dev.off())


#### Script completion message & session information ####
cat("....trna_location_analyses.R processing complete! Session information:....\n")
sessionInfo()

