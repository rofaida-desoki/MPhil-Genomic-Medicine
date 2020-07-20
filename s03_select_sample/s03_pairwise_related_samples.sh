#!/bin/bash

# s03_select_related_pairs.sh

# Last updated: Rofaida Desoki 30Mar2020

# Use on cluster:
# bsub < s03_pairwise_related_samples.sh

# BSUB instructions:

#BSUB -q inter
#BSUB -P re_gecip_inherited_cancer_predisposition
#BSUB -o s03_pairwise_related_samples_%J.out
#BSUB -e s03_pairwise_related_samples_%J.err
#BSUB -n 3
#BSUB -R "rusage[mem=4000]"
#BSUB -R "span[hosts=1]"
##BSUB -W 04:00

# Start message
echo "Select pairs of related participants"
date
echo ""

# Stop at runtime errors
set -e

# Folders (/home/rdesoki)
source_folder="/gel_data_resources/main_programme/aggregated_illumina_gvcf/GRCH38/20190228/principal_components_and_relatedness/relatedness"
working_folder="/re_gecip/inherited_cancer_predisposition/rofaida/Rproject/scripts/s03_select_sample"
cd "${working_folder}"

# Setting source and output files
source_file="${source_folder}/60k_HWE_30k_random_pairwise_kinship_full.txt"
output_file="${working_folder}/s03_pairwise_related_samples.txt"

# Progress report
echo "working_folder: ${working_folder}"
echo "source_file: ${source_file}"
echo "output_file: ${output_file}"
echo ""
echo "Top lines the source file:"
head "${source_file}"
echo ""
echo "Number of lines in the source file:"
cat "${source_file}" | wc -l  
echo ""

# Get the line with samples
echo "Selecting related pairs ..."
echo ""
awk  '$3 > 0.04419417 {print}' "${source_file}" > "${output_file}"

# Results report
echo "Top lines the output file:"
head "${output_file}"
echo ""
echo "Number of lines in the output file:"
cat "${output_file}" | wc -l  
echo ""

# Completion message
echo "Done"
date
