#!/bin/bash

# s01_check_file_formats.sh

# Summary 
# - the genes bed file generated from EnsDb encodes chromosomes in the notation of numbers: 1,2,3 ... etc  
# - the chunks bed file and bcf files use chr1,chr2,chr3 ... notation  
# Thus the genes bed file is updated to the "chr" notation, to keep the file formats consistent.  

# Started by: Alexey Larionov 26Mar2020
#Last updated: Rofaida Desoki 06Apr2020

# Use on cluster:
# bsub < s01_check_file_formats.sh

# BSUB instructions:

#BSUB -q inter
#BSUB -P re_gecip_inherited_cancer_predisposition
#BSUB -o s01_check_file_formats_%J.out
#BSUB -e s01_check_file_formats_%J.err
#BSUB -n 1
#BSUB -R "rusage[mem=4000]"
#BSUB -R "span[hosts=1]"
##BSUB -W 01:00

# Start message
echo "Explore formats of VCF and BED files"
date
echo ""

# Stop at runtime errors
set -e

# Load bcftools module
module load bcftools/1.9

# Files and folders
base_folder="/re_gecip/inherited_cancer_predisposition/rofaida/Rproject"
genes_file="${base_folder}/data/s04_genes_coordinates/genes.bed"
genes_chr_file="${base_folder}/data/s04_genes_coordinates/genes_chr.bed"

chunks_bed_folder="/gel_data_resources/main_programme/aggregated_illumina_gvcf/GRCH38/20190228/docs"
chunks_file="${chunks_bed_folder}/chunks_2019_02_28_sorted.bed"

vcf_folder="/gel_data_resources/main_programme/aggregated_illumina_gvcf/GRCH38/20190228/data"
bcf_file="${vcf_folder}/60k_GRCH38_germline_mergedgVCF_chr10_10001_3111060.bcf"

# Genes file
echo "genes_file: ${genes_file}"
echo ""
cat "${genes_file}"
echo ""

# Chunks bed file
echo "chunks_file: ${chunks_file}"
echo ""
head "${chunks_file}"
echo ""

# BCF file
echo "bcf_file: ${bcf_file}"
echo ""
bcftools view -h "${bcf_file}" | grep contig | head 
echo ""

# Update genes bed file
awk 'BEGIN{OFS="\t"} $1="chr"$1 {print}' "${genes_file}" > "${genes_chr_file}"

# Genes file
echo "genes_chr_file: ${genes_chr_file}"
echo ""
cat "${genes_chr_file}"
echo ""

# Completion message
echo "Done"
date
