#!/bin/bash

# s04_concatenate_gene_VCFs.sh

# Summary
# Concatenate VCF files for each gene, sort and index the concatenated VCF

# Started: Alexey Larionov 27Mar2020
#Last updated: Rofaida Desoki 20Apr2020

# Use on cluster:
# bsub < s04_concatenate_gene_VCFs.sh

# BSUB instructions:

#BSUB -q short
#BSUB -P re_gecip_inherited_cancer_predisposition
#BSUB -o s04_concatenate_gene_VCFs_%J.out
#BSUB -e s04_concatenate_gene_VCFs_%J.err
#BSUB -n 2
#BSUB -R "rusage[mem=8000]"
#BSUB -R "span[hosts=1]"
#BSUB -W 04:00

# Start message
echo "Concatenate VCFs, sort and index the concatenated file"
date
echo ""

# Stop at runtime errors
set -e

# Load bcftools module
module load bcftools/1.9

# Files and folders
base_folder="/re_gecip/inherited_cancer_predisposition/rofaida/Rproject"

data_folder="${base_folder}/data/s05_make_vcf"
list_of_files="${data_folder}/s03_list_of_files.txt"
concatenated_bcf="${data_folder}/s04_concatenated.bcf"
sorted_bcf="${data_folder}/s04_sorted.bcf"

scripts_folder="${base_folder}/scripts/s05_make_vcf"
cd "${scripts_folder}" 

# Progress report
echo "list_of_files: ${list_of_files}"
echo "sorted_bcf: ${sorted_bcf}"
echo ""
echo "Started concatenating, sorting and indexing ..."

# Concatenate
bcftools concat \
    --file-list "${list_of_files}" \
    --allow-overlaps \
    --remove-duplicates \
    --output-type b \
     > "${concatenated_bcf}"

# Sort
bcftools sort  \
  "${concatenated_bcf}" \
  --output-type b \
  --output-file "${sorted_bcf}"

# Iindex
bcftools index "${sorted_bcf}" 

# Remove intermediate file
rm "${concatenated_bcf}"

# Progress report
echo "Done"
echo ""

# Get basic info about content of the BCF file
echo "Counts for the output file"
echo ""
bcftools plugin counts "${sorted_bcf}" 
echo ""

# Completion message
echo "Done all tasks"
date
