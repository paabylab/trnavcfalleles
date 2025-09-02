#!/usr/bin/env nextflow
nextflow.enable.dsl=2

process bedsubvcf{
    // Subset VCF based on bed file, keeping region name, then get per-genotype counts
    label 'counts'

    input:
    tuple val(samp), path(vcf), path(vcftbi), path(bed)
    val(descrip)

    output:
    tuple val(samp), path("*.txt.gz")

    """
    # Make header info
    printf '##INFO=<ID=REGION,Number=.,Type=String,Description="Overlapping BED name(s)">\n' > header_region.info

    # bcftools subset VCF
    bcftools view -R ${bed} -Ou  ${vcf} \
    | bcftools annotate -a ${bed} -c CHROM,FROM,TO,INFO/REGION -h header_region.info -Ou \
    | bcftools view -m2 -M2 -Oz -o subset.annot.vcf.gz
    bcftools index -f -t subset.annot.vcf.gz # this one may not be needed here

    # query
    bcftools query -f '%CHROM\t%POS\t%ID\t%REF\t%ALT\t%INFO/REGION\t%N_PASS(FMT/GT="0/0" | FMT/GT="0|0")\t%N_PASS(FMT/GT="1/1" | FMT/GT="1|1")\t%N_PASS(FMT/GT="0/1" | FMT/GT="1/0" | FMT/GT="0|1" | FMT/GT="1|0")\t%N_PASS(FMT/GT="." | FMT/GT="./.")\n'   subset.annot.vcf.gz > ${samp}_${descrip}_gt_counts.txt

    # Format & zip
    echo -e "CHROM\tPOS\tREF\tALT\tID\tgene_id\tn_homref\tn_homalt\tn_het\tn_miss" > header.tmp
    cat header.tmp ${samp}_${descrip}_gt_counts.txt | gzip > ${samp}_${descrip}_gt_counts.txt.gz
    """
}
