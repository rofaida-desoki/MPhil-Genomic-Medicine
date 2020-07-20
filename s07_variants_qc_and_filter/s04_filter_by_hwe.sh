#!/bin/bash

# s04_filter_by_hwe.sh

# Started: Alexey  Larionov 18Jun2020
# Last updated: Rofaida Desoki 06Jul2020

# Filtering by HWE is included here just for illustration: all actually analysed variants will be rare, so they will have HWE ~1.  
# If we analysed common variants, then the typicaly used HWE thresholds are around 10^-6 

#BSUB -q inter
#BSUB -P re_gecip_inherited_cancer_predisposition
#BSUB -o s04_filter_by_hwe_%J.out
#BSUB -e s04_filter_by_hwe_%J.err
#BSUB -n 4
#BSUB -R "rusage[mem=8000]"
#BSUB -R "span[hosts=1]"
##BSUB -W 01:00

# Start message
echo "Filter vcf by HWE"
date
echo ""

# Stop at runtime errors
set -e

# Load bcftools module
#module load bcftools/1.10.2
module load bio/BCFtools/1.10.2-foss-2018b

# Files and folders
base_folder="/re_gecip/inherited_cancer_predisposition/rofaida/Rproject"
data_folder="${base_folder}/data/s07_variants_qc_and_filter"
scripts_folder="${base_folder}/scripts/s07_variants_qc_and_filter"
cd "${scripts_folder}" 

source_vcf="${data_folder}/s02_sort_tag_svep_clinvar_id_vkey_fix_qual.vcf.gz"
output_vcf="${data_folder}/s04_sort_tag_svep_clinvar_id_vkey_fix_qual_hwe.vcf.gz"
log="${data_folder}/s04_sort_tag_svep_clinvar_id_vkey_fix_qual_hwe.log"

# Progress report
echo "source_vcf: ${source_vcf}"
echo "output_vcf: ${output_vcf}"
echo ""

# Filter by HWE
bcftools view \
  "${source_vcf}" \
  --exclude "HWE <= 0.001" \
  --output-type z \
  --threads 4 \
  --output-file "${output_vcf}" \
  &> "${log}"

bcftools index "${output_vcf}"

# Check result
source "${scripts_folder}/f01_explore_annotations.sh"
summary="${output_vcf%.vcf.gz}_annotations_summary.txt" 
explore_annotations "${output_vcf}" > "${summary}"

# Completion message
echo "Done"
date
