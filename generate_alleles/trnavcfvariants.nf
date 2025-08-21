#!/usr/bin/env nextflow
nextflow.enable.dsl=2

/*
By Avery Davis Bell
Runs characterizations of mutations in VCF that fall within (reference) tRNA sequence(s)
*/

/*
#### Set up input parameters & defaults ####
*/

// Inputs: data
params.runinfo = "" // tab-delimited file, one row per sample/species to run. Columns:
  // id, sample/species ID (used for output file naming)
  // trnasout, highest-confidence tRNAscan-SE .out filepath OR full - just keep track of whether you're including pseudogenes or not!
  // trnass, RNAscan-SE secondary structure (usually .SS) out filepath
  // trnabed, tRNA bed file path
  // vcf, path to VCF for this species - used to find all the tRNA variants
  // vcftbi, path to VCF .tbi index
  // refstrain, ref strain to strip out before data summarization
params.out = "" // Outer/parent output directory (one per sample/species ID created internally)

// Inputs: Software related
params.vcfconda = '/storage/home/hcoda1/2/abell65/.conda/envs/trnavcf' // path to conda environment set up to run scripts that import pyvcf
params.pyscriptdir = "/storage/coda1/p-apaaby3/0/shared/labmembers/abell65/scripts/wormtrna/updatedinitial" // Directory containing python scripts run here - get_strain_variants.py etc
params.rscriptdir = "/storage/coda1/p-apaaby3/0/shared/labmembers/abell65/scripts/wormtrna/trnavcfvars" // Directory containing R scripts run here (e.g. getstrainxtrnacalls.R)

// Housekeeping:  create output directories
outdir = file(params.out)
outdir.mkdirs()

/*
#### Processes ####
*/

process subsetvcf{
  // Slices vcf to only be tRNA regions (otherwise collectvarstrnas is quite slow)
  // AND Edits VCF so python vcf reader can work with it (unzips; removes "malformed filter lines", etc)

  input:
    tuple val(myid), path(inittrnasout), path(trnass), path(trnabed), path(vcf), path(vcftbi), val(refstrain)

  output:
    tuple val(myid), path(inittrnasout), path(trnass), val(refstrain), path("*reformatted.bed"), path("*edited.vcf")

  """
  # Strip any errant "chrs" from bed files [careful, this might not work for all species/VCFs]
  sed 's/chr//g' ${trnabed} > ${myid}_reformatted.bed

  # Subset vcf, save un-compressed
  bcftools view -H \
  --regions-file ${myid}_reformatted.bed \
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

process collectvarstrnas{
  // Runs get_strain_variants.py (Collects variants in tRNAs (important preprocessing step for next scripts))

  errorStrategy 'finish'
  conda params.vcfconda
  publishDir outdir, mode: 'copy', saveAs: { filename -> "$myid/data/$filename" }, pattern: "*variants.txt"
  publishDir outdir, mode: 'copy', saveAs: { filename -> "$myid/data/$filename" }, pattern: "*variants.log"

  input:
    tuple val(myid), path(inittrnasout), path(trnass), val(refstrain), path(trnabed), path(vcf)
    // GET ALL ONES NEEDED FOR INPUTS, won't use all of them but easier to keep all together

  output:
    tuple val(myid), path(inittrnasout), path(trnass), val(refstrain), path(trnabed), path(vcf), path("*_strain_variants.txt"), emit: data
    path("*.log"), emit: log

  """
  python ${params.pyscriptdir}/get_strain_variants.py \
  -trnas ${inittrnasout} \
  -vcf ${vcf} \
  -out ${myid}
  """
}

process buildaltseqs{
 // Runs build_alt_sequences.py (Generates strain-specific tRNA sequences)

 errorStrategy 'finish'
 conda params.vcfconda
 publishDir outdir, mode: 'copy', saveAs: { filename -> "$myid/data/$filename" }, pattern: "*gen.txt"
 publishDir outdir, mode: 'copy', saveAs: { filename -> "$myid/data/$filename" }, pattern: "*strain_trnas.fa"
 publishDir outdir, mode: 'copy', saveAs: { filename -> "$myid/data/$filename" }, pattern: "*strain_trnas_info.txt"
 publishDir outdir, mode: 'copy', saveAs: { filename -> "$myid/data/$filename" }, pattern: "*seqs.log"

 input:
  tuple val(myid), path(inittrnasout), path(trnass), val(refstrain), path(trnabed), path(vcf), path(strainvars)

 output:
  tuple val(myid), val(refstrain), path(vcf), path(strainvars), path("*strain_trnas_gen.txt"), path("*strain_trnas.fa"), path("*strain_trnas_info.txt"), emit: data
  path("*.log"), emit: log
  // _strain_trnas_gen.txt, _strain_trnas.fa, _strain_trnas_info.txt
  // with all - think it's too confusing but can reinstate if needed: tuple val(myid), path(inittrnasout), path(trnass), path(trnabed), path(vcf), path(strainvars), path("*strain_trnas_gen.txt"), path("*strain_trnas.fa"), path("*strain_trnas_info.txt"), emit: data

 """
 python ${params.pyscriptdir}/build_alt_sequences.py \
 -trnasout ${inittrnasout} \
 -trnass ${trnass} \
 -trnabed ${trnabed}\
 -vcf ${vcf} \
 -vars ${strainvars} \
 -out ${myid}
 """
}

process straintrnascanse{
  // Run tRNAscan-SE on strain-specific FASTA

  errorStrategy 'finish'
  publishDir outdir, mode: 'copy', saveAs: { filename -> "$myid/data/trnascan/$filename" }, pattern: "*strainspec_trnas*"

  input:
    tuple val(myid), val(refstrain), path(vcf), path(strainvars), path(strain_trnas_gen), path(strain_trnas_fa), path(strain_trnas_info)

  output:
    tuple val(myid), path("*strainspec_trnas.out"), path("*strainspec_trnas.SS"), path("*strainspec_trnas.stats"), path("*strainspec_trnas_isospecific.out") // paths of what **it** makes

  """
  # **** CONSIDER SETTING TMPDIR HERE

  tRNAscan-SE \
  -o ${myid}_strainspec_trnas.out \
  -f ${myid}_strainspec_trnas.SS \
  -m ${myid}_strainspec_trnas.stats \
  --log ${myid}_strainspec_trnas.log \
  --progress \
  --detail \
  --isospecific ${myid}_strainspec_trnas_isospecific.out \
  ${strain_trnas_fa}
  """
}

process missingness{
  // Run catalogmissingness_trnavcfscripts.R to get record of where missing VCF genotypes are

  publishDir outdir, mode: 'copy', saveAs: { filename -> "$myid/data/$filename" }, pattern: "*per*.txt"
  publishDir outdir, mode: 'copy', saveAs: { filename -> "$myid/summaries/$filename" }, pattern: "*.pdf"

  input:
    tuple val(myid), val(refstrain), path(vcf), path(strainvars), path(strain_trnas_gen), path(strain_trnas_fa), path(strain_trnas_info)

  output:
    tuple val(myid), val(refstrain), path(strain_trnas_fa), path(strain_trnas_info), path("*pervariant.txt"), path("*pertrna.txt"), path("*perstrain.txt"), emit: data
    tuple path("*_pertrna.pdf"), path("*_perstrain.pdf"), emit: plots

  """
  module load r/4.2.1-tidy # obviously not ideal, but for some reason R doesn't load well through the config file

  Rscript ${params.rscriptdir}/catalogmissingness_trnavcfscripts.R \
  --varinfo ${strainvars}\
  --trnainfo ${strain_trnas_info} \
  --baseoutname ${myid} \
  --outdir "getwd()"
  """
}

process strainxtrna{
  // Run getstrainxtrnacalls.R to get final data and summaries

  publishDir outdir, mode: 'copy', saveAs: { filename -> "$myid/summaries/$filename" } // make sure this grabs just newly generated files

  input:
    tuple val(myid), val(refstrain), path(strain_trnas_fa), path(strain_trnas_info), path(mpervar), path(mpert), path(mpers), path(strainspec_out), path(strainspec_SS), path(strainspec_stats), path(strainspec_isospecific)

  output:
    path("*")

  """
  module load r/4.2.1-tidy # obviously not ideal, but for some reason R doesn't load well through the config file

  Rscript ${params.rscriptdir}/getstrainxtrnacalls.R \
  --trnainfo ${strain_trnas_info} \
  --missinfo ${mpert} \
  --trnafasta ${strain_trnas_fa} \
  --trnascanout ${strainspec_out} \
  --refstrain ${refstrain} \
  --baseoutname ${myid} \
  --outdir "getwd()"
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
      return [row.id, row.trnasout, row.trnass, row.trnabed, row.vcf, row.vcftbi, row.refstrain]
    }
    .set{runInfo}

  // Generate data
  subsetvcf(runInfo) | collectvarstrnas
  buildaltseqs(collectvarstrnas.out.data)
  missingness(buildaltseqs.out.data)
  straintrnascanse(buildaltseqs.out.data)

  // summarize & collect data
  missingness.out.data | combine(straintrnascanse.out, by: 0) | set{ combodat }
  strainxtrna(combodat)
}
