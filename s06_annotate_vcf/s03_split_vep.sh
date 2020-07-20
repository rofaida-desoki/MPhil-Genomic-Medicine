#!/bin/bash

# Split VEP
# s03_split_vep.sh

# Started: Alexey Larionov, 03May2020
# Last updated: Rofaida Desoki, 05Jun2020

# Intended use:
# bsub < s03_split_vep.sh

# Note:
# VEP writes annotations in a single CSQ filed within INFO column
# split-vep creates separate fields in the INFO column for each/selected VEP annotations

#BSUB -q inter
#BSUB -P re_gecip_inherited_cancer_predisposition
#BSUB -o s03_split_vep_%J.stdout
#BSUB -e s03_split_vep_%J.stderr
#BSUB -R "rusage[mem=6000]"
#BSUB -n 1
##BSUB -W 01:00

# Stop at runtime errors
set -e

# Start message
echo "Split VEP"
date
echo ""

# Load modules
module load bcftools/1.10.1

# Folders
base_folder="/re_gecip/inherited_cancer_predisposition/rofaida/Rproject"
data_folder="${base_folder}/data/s06_annotate_vcf"
scripts_folder="${base_folder}/scripts/s06_annotate_vcf"
cd "${scripts_folder}"

source_vcf="${data_folder}/s02_sort_tag_vep.vcf.gz"
tmp_vcf="${data_folder}/vep_split_tmp.vcf.gz"
output_vcf="${data_folder}/s03_sort_tag_svep.vcf.gz"

log="${data_folder}/s03_split_vep.log"

# Progress report
echo "--- Input and output files ---"
echo ""
echo "source_vcf: ${source_vcf}"
echo "tmp_vcf: ${tmp_vcf}"
echo "output_vcf: ${output_vcf}"
echo "log: ${log}"
echo ""
echo "--- List of detected VEP files ---"
echo ""
bcftools +split-vep -l "${source_vcf}"
echo ""
echo "Splitting VEP fields ..."

# Get number of VEP annotations
n=$(bcftools +split-vep -l "${source_vcf}" | wc -l)
n=$(( $n - 1 ))

# Split VEP field
bcftools +split-vep \
  -c 0-$n \
  --annot-prefix vep_ \
  --output "${tmp_vcf}" \
  --output-type z \
  "${source_vcf}" \
  &> "${log}"

# Index interim VCF
bcftools index "${tmp_vcf}"

# Remove unsplit CSQ column
bcftools annotate \
  --remove INFO/CSQ \
  --output "${output_vcf}" \
  --output-type z \
  "${tmp_vcf}" \
  &> "${log}"

# Index output VCF
bcftools index "${output_vcf}"

# Remove interim VCF and index
rm "${tmp_vcf}" "${tmp_vcf}.csi"

# Summary of INFO fields
echo "Number of INFO fields in split vep vcf:"
bcftools view -h "${output_vcf}" | grep ^##INFO | wc -l
echo ""

echo "List of INFO fields in split vep vcf:"
bcftools view -h "${output_vcf}" | grep ^##INFO
echo ""

# Completion message
echo "Done"
date
