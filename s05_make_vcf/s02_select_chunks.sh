#!/bin/bash

# s02_select_chunks.sh

# Summary
# GEL split the multi-sample VCF into ~1000 chunks to keep reasonable file sizes.  
# Use bedtools to find chunks that overlap with selected gene coordinates.  
# Write the chunks suffixes to a text file for downstream processing.  
# Add 1000 bp padding to the gene coordinates on both sides  

# Reference
# https://bedtools.readthedocs.io/en/latest/content/bedtools-suite.html 

# Started: Alexey Larionov 26Mar2020
# Last updated: Rofaida Desoki 07Apr2020

# Use on cluster:
# bsub < s02_select_chunks.sh

# BSUB instructions:

#BSUB -q inter
#BSUB -P re_gecip_inherited_cancer_predisposition
#BSUB -o s02_select_chunks_%J.out
#BSUB -e s02_select_chunks_%J.err
#BSUB -n 1
#BSUB -R "rusage[mem=4000]"
#BSUB -R "span[hosts=1]"
##BSUB -W 01:00

# Start message
echo "Select VCF chunks containing data for required genes"
date
echo ""

# Stop at runtime errors
set -e

# Load bcftools module
module load bedtools/2.26.0

# Files and folders
chunks_folder="/gel_data_resources/main_programme/aggregated_illumina_gvcf/GRCH38/20190228/docs"
chunks_file="${chunks_folder}/chunks_2019_02_28_sorted.bed"

base_folder="/re_gecip/inherited_cancer_predisposition/rofaida/Rproject"
genes_chr_file="${base_folder}/data/s04_genes_coordinates/genes_chr.bed"

output_folder="${base_folder}/data/s05_make_vcf"
tmp_file_1="${output_folder}/s02_tmp_1.txt" # a good practice of making temporary files uses mktemp command
tmp_file_2="${output_folder}/s02_tmp_2.txt" 
output_file="${output_folder}/s02_selected_chunks.txt"

# Progress report
echo "chunks_file: ${chunks_file}"
echo "genes_chr_file: ${genes_chr_file}"
echo "output_file: ${output_file}"
echo ""
echo "Started ..."

# Make the output folder
rm -fr "${output_folder}" # remove if existed
mkdir -p "${output_folder}"

# Intesect as left outer join
bedtools intersect -loj -a "${chunks_file}" -b "${genes_chr_file}" > "${tmp_file_1}"

# Keep only lines where the intersect was found (no need for BEGIN option since we aren't modifying anything
awk '$6 != "." {print}' "${tmp_file_1}" > "${tmp_file_2}"

# Keep only required columns: chr, start, end, gene and file suffix.
# Add 1000 bp padding to start and end
awk 'BEGIN{OFS="\t"} {print $1, $7-1000, $8+1000, $9, $5}' "${tmp_file_2}" > "${output_file}"
# This awk could be combined with the previous awk.
# It is kept separate just for clarity and for the opportunity to explore intermediate files.  

# Completion message
echo "Done"
date
