#!/usr/bin/env nextflow
nextflow.enable.dsl=2

/*
By Avery Davis Bell
Gets VCF records overlapping tRNA and protein-coding-exon regions, counts genotypes for these
*/

/*
#### Set up input parameters & defaults ####
*/
// Inputs: data
params.runinfo = "" // tab-delimited file, one row per sample/species to run. Columns:
// id, sample/species ID (used for output file naming)
// vcf, path to whole genome vcf to subset
// vcftbi, path to .tbi index of above
// trnabed, path to bed file for all tRNA genes of interest, name in 4th column
// gtf, path to GTF file containing (at least) all protein-coding genes' exons
params.out = "" // output directory (multiple files created internally)

// Housekeeping:  create output directories
outdir = file(params.out)
outdir.mkdirs()

/*
#### Processes ####
*/

include { bedsubvcf as bedsubvcf_trna } from './modules/bedsubvcf.nf'
include { bedsubvcf as bedsubvcf_exon } from './modules/bedsubvcf.nf'


process gtf2exonbed{
    // Gets an exon bed with gene names in 4th column from GTF file
    publishDir outdir, mode: 'copy', overwrite: true, pattern: "*.bed"

    input:
    tuple val(samp), path(vcf), path(vcftbi), path(gtf)

    output:
    tuple val(samp), path(vcf), path(vcftbi), path("*.bed")

    """
    # Run as sh script to get around awk-inside-nextflow issues
    sh gtf2protexonsbed.sh ${gtf} ${samp}

    """
}

/*
#### Workflow ####
*/

workflow{
    // Set up channels
    formatChannel = Channel.fromPath(params.runinfo)
    .splitCsv(header: true, sep: '\t', strip: true)
    .map{row->
      return [row.id, row.vcf, row.vcftbi, row.trnabed]
    }
    .set{runInfotRNA}

    formatChannel = Channel.fromPath(params.runinfo)
    .splitCsv(header: true, sep: '\t', strip: true)
    .map{row->
      return [row.id, row.vcf, row.vcftbi, row.gtf]
    }
    .set{runInfoExon}

    // Do the work (fixed because can't call same process 2x separately)
    bedsubvcf_trna(runInfotRNA, "tRNA") // tRNA counting
    gtf2exonbed(runInfoExon) // exon formatting
    bedsubvcf_exon(gtf2exonbed.out, "protcodexonic") // exon counting

}
