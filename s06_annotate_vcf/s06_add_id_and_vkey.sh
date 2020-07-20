#!/bin/bash

# s06_add_id_and_vkey.sh

# Fill ID column with Chrom-Pos-Ref-Alt data and add VariantKey, which represents the same data represented in some hexadecimal coding. 
# Assuming that no lines with identical chr, pos, ref and alt are present in VCF,
# these ID-s may be helpful for matching variants between R and VCF, later if needed.

# References about Variant Key:
# https://github.com/samtools/bcftools/pull/831
# https://github.com/Genomicsplc/variantkey
# https://doi.org/10.1101/473744

# Started: Alexey  Larionov, 23Apr2020
# Last updated: Rofaida Desoki, 05Jun2020

#BSUB -q inter
#BSUB -P re_gecip_inherited_cancer_predisposition
#BSUB -o s06_add_variant_key_%J.out
#BSUB -e s06_add_variant_key_%J.err
#BSUB -n 4
#BSUB -R "rusage[mem=8000]"
#BSUB -R "span[hosts=1]"
##BSUB -W 01:00

# Start message
echo "Add ID and Variant Key"
date
echo ""

# Stop at runtime errors
set -e

# Load bcftools module (note version of bcftools)
module load bcftools/1.10.1

# Files and folders
base_folder="/re_gecip/inherited_cancer_predisposition/rofaida/Rproject"
data_folder="${base_folder}/data/s06_annotate_vcf"
scripts_folder="${base_folder}/scripts/s06_annotate_vcf"
cd "${scripts_folder}" 

source_vcf="${data_folder}/s05_sort_tag_svep_clinvar.vcf.gz"

id_vcf="${data_folder}/s06_sort_tag_svep_clinvar_id.vcf.gz"
vkey_vcf="${data_folder}/s06_sort_tag_svep_clinvar_id_vkey.vcf.gz"

id_log="${data_folder}/s06_add_id.log"
vkey_log="${data_folder}/s06_add_vkey.log"

# Progress report
echo "source_vcf: ${source_vcf}"
echo "id_vcf: ${id_vcf}"
echo "vkey_vcf: ${vkey_vcf}"
echo ""

# Write ID
echo "Writing to ID column ..."
bcftools annotate \
  --set-id '%CHROM\_%POS\_%REF\_%ALT' \
  --output "${id_vcf}" \
  --output-type z \
  --threads 4 \
  "${source_vcf}" \
  &> "${id_log}"

# Index interim vcf
bcftools index "${id_vcf}"

# Annotate with bcftools add-variantkey plugin
echo "Adding VariantKey ..."
bcftools  +add-variantkey \
  --output "${vkey_vcf}" \
  --output-type z \
  --threads 4 \
  "${id_vcf}" \
  &> "${vkey_log}"

# Index output vcf
bcftools index "${vkey_vcf}"

# Remove interim VCF and index
rm "${id_vcf}" "${id_vcf}.csi" 

# Progress report
echo "Done"
echo ""
echo "source_vcf:"
echo ""
bcftools view "${source_vcf}" | tail | awk '{print $1,$2,$3,$4,$5,$6,$7,$8,$9,$10}'
echo ""
echo "output_vcf:"
echo ""
bcftools view "${vkey_vcf}" | tail | awk '{print $1,$2,$3,$4,$5,$6,$7,$8,$9,$10}'
echo ""

# Explore IDs and VariantKeys
id_vkey_tmp="${data_folder}/s06_id_vkey_tmp.txt"
bcftools query -f '%ID\t%VKX\n' "${vkey_vcf}" > "${id_vkey_tmp}"

echo "Esamples of generated ID and VariantKeys:"
head "${id_vkey_tmp}"
echo ""

echo "Number of generated variants:"
cat "${id_vkey_tmp}" | wc -l 
echo ""

echo "Number of unique IDs:"
cut -f1 "${id_vkey_tmp}" | sort | uniq | wc -l
echo ""

echo "Number of unique VariantKeys:"
cut -f2 "${id_vkey_tmp}" | sort | uniq | wc -l
echo ""

rm "${id_vkey_tmp}"

# Completion message
echo "Completed script"
date
