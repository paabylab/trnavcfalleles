#!/usr/bin/env nextflow
nextflow.enable.dsl=2

/*
By Avery Davis Bell
Runs secstruct2seqpieces.py to break tRNAs into pieces,
then does alignments of pieces and recombines.
*/

/*
#### Set up input parameters & defaults ####
*/

// Inputs: data
params.runinfo = "" // tab-delimited file, one row per sample/species to run. Columns:
// id, sample/species ID (used for output file naming)
// trnass, tRNAscan-SE .ss filepath for ALL alleles
// specprefix, species prefix (no spaces) for FASTA seq naming for combining across species
params.out = "" // Outer/parent output directory (multiple directories created internally)
params.alignmentout = "all_species_alleles_foursale_aligned.xfasta.gz" // File name for final alignment file

// Inputs: software parameter Related
params.gapopen = 0 // -gapopen passed to foursale/clustal22
params.gapextend = 0 // -gapextend passed to foursale/clustalw2
params.combfull = "True" // -full for fastax concat'ing: True or False: keep all sequences, like full/outer join
params.combfullfill = "N" // -fill for fastax concat'ing: "fill with N bases/residues for IDs missing in some files when using -full"
params.replaceU = "True" // -replaceU for fastax concat'ing: True or False: replace any Us in seq with Ts

// Inputs: Software related
params.pyscriptdir = "/storage/coda1/p-apaaby3/0/shared/labmembers/abell65/scripts/wormtrna" // directory containing script secstructrelated/secstruct2seqpieces.py; utilityscripts/concatxfastas.py
params.foursalerunjar = "-cp .:/storage/coda1/p-apaaby3/0/shared/software/foursale-1.3/foursale.jar:/storage/coda1/p-apaaby3/0/shared/software/java-classes/commons-cli-1.9.0/commons-cli-1.9.0.jar:/storage/coda1/p-apaaby3/0/shared/software/java-classes/commons-lang3-3.17.0/commons-lang3-3.17.0.jar:/storage/coda1/p-apaaby3/0/shared/software/java-classes/commons-io-2.18.0/commons-io-2.18.0.jar de.biozentrum.bioinformatik.foursale_cmd.Main" // String to put after java to get foursale to run. Should include all cp dependencies.
params.clustalw2 = "/storage/coda1/p-apaaby3/0/shared/software/clustalw-2.1-linux-x86_64-libcppstatic/clustalw2" // fully articulated path to clustalw directory that runs with foursale (passed to foursale)
params.vcfconda = '/storage/home/hcoda1/2/abell65/.conda/envs/trnavcf' // path to conda environment set up to run scripts that import pyvcf (don't actually need pyvcf, but this one has other dependencies I use)

// Housekeeping:  create output directories
outdir = file(params.out)
outdir.mkdirs()

outsspiecedir = file(params.out + "/allelesecstruct")
outsspiecedir.mkdirs()

outalndir = file(params.out + "/fullalignments")
outalndir.mkdirs()

/*
#### Processes ####
*/

process secstruct2seqpieces{
  // runs secstruct2seqpieces.py

  conda params.vcfconda
  publishDir outsspiecedir, mode: 'copy', overwrite: true, pattern: "*"

  input:
  tuple val(myid), path(ss), val(sname)

  output:
  tuple val(myid), val(sname), path("*"), emit: all
  path("*"), emit: dirpaths

  """
  mkdir -p ${myid}
  python ${params.pyscriptdir}/secstructrelated/secstruct2seqpieces.py \
  -trnass ${ss} \
  -outdir ${myid} \
  -out ${myid} \
  -fasta True \
  -xfasta True \
  -fanameprefix ${sname} \
  -maskanticodon True
  """
}

process combinexfastas{
  // Combines X fastas across species (within secondary structure pieces)
  input:
  path spdir
  val ssname

  output:
  tuple val(ssname), path("*combined.xfasta"), emit: all
  path("*combined.xfasta"), emit: fpaths
  """
  cat */*${ssname}.xfasta > ${ssname}_combined.xfasta
  """
}


process foursale{
  // Runs foursale for each secondary structure alignment piece

  errorStrategy 'finish'

  input:
  tuple val(ssname), path(xfasta)

  output:
  tuple val(ssname), path("*foursale_aligned.xfasta"), emit: all
  path("*foursale_aligned.xfasta"), emit: xfasta

  """
  java ${params.foursalerunjar} \
  -clustalw2 ${params.clustalw2} \
  -in ${xfasta} \
  -out ${ssname}_foursale_aligned.xfasta \
  -gapextend ${params.gapextend} \
  -gapopen ${params.gapopen}
  """
}

process combinealignments{
  // runs concatxfastas.py
  // Combines aligned outputs across secondary structures; [ideally subs in codon seq, too, here??]; and maybe keep length of each input alignment just in case want it

  conda params.vcfconda
  publishDir outalndir, mode: 'copy', overwrite: true, pattern: "*"

  input:
  path alnstructs // all the ones to combine in here...
  path(combfast) // anticodon fasta in here...

  output:
  path("*")

  """
  python ${params.pyscriptdir}/utilityscripts/concatxfastas.py \
  -files accstemL2d_foursale_aligned.xfasta,darm2ant_foursale_aligned.xfasta,antarm_foursale_aligned.xfasta,varloop_foursale_aligned.xfasta,tarm2acc_foursale_aligned.xfasta,accstemR2end_foursale_aligned.xfasta \
  -full ${params.combfull} \
  -fill ${params.combfullfill} \
  -replaceU ${params.replaceU} \
  -replaceCodon True \
  -codonfasta anticodon_combined.xfasta \
  -out ${params.alignmentout}
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
      return [row.id, row.trnass, row.specprefix]
    }
    .set{runInfo}
  sschan = channel.fromList( ["accstemL2d", "darm2ant","antarm","anticodon","varloop", "tarm2acc", "accstemR2end"] ) // NB would need to update if changed secstruct2seqpieces.py!
    // need anticodon one to have it for replacement - DON'T need to align that one but probably easiest to just do...

  // By-species preprocessing
  secstruct2seqpieces(runInfo)
  allspec = secstruct2seqpieces.out.dirpaths.collect()

  // By-structure processing
  combinexfastas(allspec, sschan)
  foursale(combinexfastas.out.all)

  // Combine across structures
  combinealignments(foursale.out.xfasta.collect(), combinexfastas.out.fpaths.collect())
}
