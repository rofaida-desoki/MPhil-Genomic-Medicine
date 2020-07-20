#!/bin/bash

# s07_explore_annotations.sh

# Started: Alexey Larionov, 11May2020
# Updated: Alexey 18Jun2020
# Last updated: Rofaida Desoki, 25Jun2020

# Running many bcftools queries about separate fields

#BSUB -q inter
#BSUB -P re_gecip_inherited_cancer_predisposition
#BSUB -o s07_explore_annotations_%J.stdout
#BSUB -e s07_explore_annotations_%J.stderr
#BSUB -R "rusage[mem=6000]"
#BSUB -n 1
## BSUB -W 01:00

# Stop at runtime errors
set -e

# Start message
echo "Extact and explore selected annotations from VCF file"
date
echo ""

# Load modules
#module load bcftools/1.10.2
module load bio/BCFtools/1.10.2-foss-2018b

# Files and folders
base_folder="/re_gecip/inherited_cancer_predisposition/rofaida/Rproject"
data_folder="${base_folder}/data/s06_annotate_vcf"

scripts_folder="${base_folder}/scripts/s06_annotate_vcf"
cd "${scripts_folder}" 

vcf="${data_folder}/s06_sort_tag_svep_clinvar_id_vkey.vcf.gz"
summary="${vcf%.vcf.gz}_annotations_summary.txt"

# Progress report
echo "vcf: ${vcf}"
echo "summary: ${summary}"
echo ""

# Load function and explore annotations
echo "Extracting annotations and preparing summary ..."
echo ""

source "${scripts_folder}/f01_explore_annotations.sh"
explore_annotations "${vcf}" > "${summary}"

# Completion message
echo "Done" 
date
