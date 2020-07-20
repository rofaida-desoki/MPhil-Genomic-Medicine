#!/bin/bash

# s01_fill_tags.sh

# Summary:
# Set INFO tags: AF, AC, AC_Hemi, AC_Hom, AC_Het, AN, ExcHet, HWE, MAF, NS

# Started: Alexey Larionov 22Apr2020
# Last updated: Rofaida Desoki 23Apr2020

#BSUB -q inter
#BSUB -P re_gecip_inherited_cancer_predisposition
#BSUB -o s01_fill_tags_%J.out
#BSUB -e s01_fill_tags_%J.err
#BSUB -n 4
#BSUB -R "rusage[mem=8000]"
#BSUB -R "span[hosts=1]"
##BSUB -W 01:00

# Start message
echo "Fill tags in VCF file"
date
echo ""

# Stop at runtime errors
set -e

# Load bcftools module
module load bcftools/1.9

# Files and folders
base_folder="/re_gecip/inherited_cancer_predisposition/rofaida/Rproject"
data_folder="${base_folder}/data/s06_annotate_vcf"
scripts_folder="${base_folder}/scripts/s06_annotate_vcf"

rm -fr "${data_folder}" # remove results folder, if existed
mkdir -p "${data_folder}"
cd "${scripts_folder}" 

source_bcf="${base_folder}/data/s05_make_vcf/s04_sorted.bcf"
output_vcf="${data_folder}/s01_sort_tag.vcf.gz"
log="${data_folder}/s01_sort_tag.log"

# Progress report
echo "source_bcf: ${source_bcf}"
echo "output_vcf: ${output_vcf}"
echo ""
echo "Started ..."

# Annotate with bcftools fill-tags plugin (note the output file format)
bcftools plugin \
  fill-tags \
  --output "${output_vcf}" \
  --output-type z \
  --threads 4 \
  "${source_bcf}" \
  &> "${log}"

# Index output vcf
bcftools index "${output_vcf}"

# Progress report
echo "Done"
echo ""
echo "source_bcf:"
echo ""
bcftools view "${source_bcf}" | tail  -n 1 | awk '{print $1,$2,$3,$4,$5,$6,$7,$8,$9,$10}'
echo ""
echo "output_vcf:"
echo ""
bcftools view "${output_vcf}" | tail  -n 1 | awk '{print $1,$2,$3,$4,$5,$6,$7,$8,$9,$10}'
echo ""

# Completion message
echo "Completed script"
date
