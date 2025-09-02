#!/usr/bin/env bas
set -eu

# Inputs
gtf=$1
samp=$2

#### Script ####

zcat -f ${gtf} | awk -F'\t' '
BEGIN { OFS = "\t" }
# Skip comments
/^#/ { next }
# Keep only exon features and exclude MtDNA/mtDNA chromosomes
$3 != "exon" || $1 ~ /^(MtDNA|mtDNA)$/ { next }
{
chrom = $1; start = $4; end = $5; attrs = $9

# Initialize attributes
gene_id = tx_biotype = gene_biotype = product = ""

# Split attribute column on semicolons
n = split(attrs, parts, /;[ \t]*/)
for (i = 1; i <= n; i++) {
    s = parts[i]
    # trim
    gsub(/^[ \t]+|[ \t]+$/, "", s)
    if (s == "") continue

    # Handle key=value and key "value"
    key = s; val = ""
    if (index(s, "=")) {
        key = s; sub(/=.*/, "", key)
        val = s; sub(/^[^=]*=/, "", val)
    } else if (index(s, " ")) {
        key = s; sub(/ .*/, "", key)
        val = s; sub(/^[^ ]* /, "", val)
    }

    # Remove surrounding quotes if present
    gsub(/^"/, "", val); gsub(/"$/, "", val)

    if (key == "gene_id")           gene_id = val
    else if (key == "transcript_biotype") tx_biotype = val
    else if (key == "gene_biotype") gene_biotype = val
    else if (key == "product")      product = val
}

# Keep if protein-coding or hypothetical protein product
accept = 0
low_tx   = tolower(tx_biotype)
low_gene = tolower(gene_biotype)
low_prod = tolower(product)

if (low_tx == "protein_coding" || low_gene == "protein_coding") accept = 1
if (low_prod ~ /hypothetical protein/) accept = 1
if (!accept) next

# BED is 0-based: start-1
print chrom, start-1, end, gene_id
}
' > ${samp}_protcodexons.bed
