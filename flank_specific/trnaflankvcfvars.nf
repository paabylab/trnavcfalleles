#!/usr/bin/env nextflow
nextflow.enable.dsl=2

/*
By Avery Davis Bell
Runs characterizations of mutations in VCF that fall within (reference) tRNA sequence FLANKS
*/

/*
#### Set up input parameters & defaults ####
*/

// Inputs: data
params.runinfo = "" // tab-delimited file, one row per sample/species to run. Columns:
  // id, sample/species ID (used for output file naming)
  // trnasout, tRNAscan-SE .out filepath
  // trnabed, tRNA bed file path
  // vcf, path to VCF for this species - used to find all the tRNA variants
  // vcftbi, path to VCF .tbi index
params.out = "" // Outer/parent output directory (one per sample/species ID created internally)
params.totalflank = 40 // flank length to get for tRNA flanking bed. **Just for output bed file, this does NOT change how Corinne's script works
params.innerflank = 20 // flank length for inner tRNA flank for output bed file (outer is this to total).  **Just for output bed file, this does NOT change how Corinne's script works

// Inputs: Software related
params.vcfconda = '/storage/home/hcoda1/2/abell65/.conda/envs/trnavcf' // path to conda environment set up to run scripts that import pyvcf
params.wormtrnarepo = "/storage/coda1/p-apaaby3/0/shared/labmembers/abell65/scripts/wormtrna" // Directory containing gitrepo python scripts run here -
  // updatedinitial/get_strain_variants-flank.py and utilityscripts/bed2flanks.py

// Housekeeping:  create output directories
outdir = file(params.out)
outdir.mkdirs()

/*
#### Processes ####
*/

process bed4vcf{
  // Generate bed file that includes tRNAs and generous flanking regions

  input:
  tuple val(myid), path(trnaout), path(trnabed), path(vcf), path(vcftbi)

  output:
  tuple val(myid), path(trnaout), path(vcf), path(vcftbi), path("*_flank_generous.bed")

  """
  # Strip any errant "chrs" from bed files [careful, this might not work for all species/VCFs]
  sed 's/chr//g' ${trnabed} > reformatted.bed

  # Get +/- 1kb (well more than enough)
  awk '{print \$1"\t"\$2-1000"\t"\$3+1000}' reformatted.bed > newpos.bed
  cut -f4- reformatted.bed | paste newpos.bed - > ${myid}_trnas_flank_generous.bed
  """
}

process subsetvcf{
  // slice VCF to be only tRNA + flanking regions;
  // Edits VCF so python vcf reader can work with it (unzips; removes "malformed filter lines", etc)

  input:
  tuple val(myid), path(trnaout), path(vcf), path(vcftbi), path(trnaflankbed)

  output:
  tuple val(myid), path(trnaout), path("*edited.vcf")

  """
  # Subset vcf, save un-compressed
  bcftools view -H \
  --regions-file ${trnaflankbed} \
  -O v \
  -o ${myid}_bedregions.vcf \
  ${vcf}

  # Get unzipped header & clean it
  bcftools view -h ${vcf} > initheader.txt
  grep -v "##FILTER" initheader.txt > cleanheader.txt

  # Combine into 'final file'
  cat cleanheader.txt ${myid}_bedregions.vcf > ${myid}_subset_edited.vcf
  """
}

process strainflankvars{
  // Run Corinne's get_strain_variants-flank.py

  errorStrategy 'finish'
  conda params.vcfconda
  publishDir outdir, mode: 'copy', overwrite: true, saveAs: { filename -> "$myid/$filename" }, pattern: "*variants-flank.txt"
  publishDir outdir, mode: 'copy', overwrite: true, saveAs: { filename -> "$myid/$filename" }, pattern: "*variants-flank.log"

  input:
  tuple val(myid), path(trnaout), path(smallvcf)

  output:
  tuple val(myid), path("*variants-flank.txt"), emit: data
  path("*.log"), emit: log

  """
  python ${params.wormtrnarepo}/updatedinitial/get_strain_variants-flank.py \
  -trnasout ${trnaout} \
  -vcf ${smallvcf} \
  -out ${myid}
  """
}

process bedflanks{
  // Makes bed file with exact flank regions of interest (split into inner, outer, 5', 3' and named as such)
  // Generated from tRNA bed file

  conda params.vcfconda
  publishDir outdir, mode: 'copy', overwrite: true, saveAs: { filename -> "$myid/$filename" }, pattern: "*trna_flanks.bed"

  input:
  tuple val(myid), path(trnaout), path(trnabed), path(vcf), path(vcftbi)

  output:
  tuple val(myid), path("*trna_flanks.bed")

  """
  python ${params.wormtrnarepo}/utilityscripts/bed2flanks.py \
  -bed ${trnabed} \
  -nflank ${params.totalflank} \
  -ninner ${params.innerflank} \
  -out ${myid}_trna_flanks.bed
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
      return [row.id, row.trnasout, row.trnabed, row.vcf, row.vcftbi]
    }
    .set{runInfo}

  // Process through get_strain_variants-flank.py
  bed4vcf(runInfo) | subsetvcf | strainflankvars

  // (Separately) get bed file of nicely-delineated flank regions
  bedflanks(runInfo)
}
