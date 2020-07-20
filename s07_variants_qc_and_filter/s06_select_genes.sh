#!/bin/bash

# s06_select_genes.sh

# Started: Alexey  Larionov 18Jun2020
# Last updated: Rofaida Desoki 08Jul2020

#BSUB -q inter
#BSUB -P re_gecip_inherited_cancer_predisposition
#BSUB -o s06_select_genes_%J.out
#BSUB -e s06_select_genes_%J.err
#BSUB -n 4
#BSUB -R "rusage[mem=8000]"
#BSUB -R "span[hosts=1]"
##BSUB -W 01:00

# Start message
echo "Select genes"
date
echo ""

# Stop at runtime errors
set -e

# Load bcftools module
#module load bcftools/1.10.2
module load bio/BCFtools/1.10.2-foss-2018b

# Folders
base_folder="/re_gecip/inherited_cancer_predisposition/rofaida/Rproject"
data_folder="${base_folder}/data/s07_variants_qc_and_filter"
scripts_folder="${base_folder}/scripts/s07_variants_qc_and_filter"
cd "${scripts_folder}" 

# Files
source_vcf="${data_folder}/s05_sort_tag_svep_clinvar_id_vkey_fix_qual_hwe_gnomad.vcf.gz" 
output_vcf="${data_folder}/s06_genes.vcf.gz" 
log="${data_folder}/s06_genes.log"

# Progress report
echo "source_vcf: ${source_vcf}"
echo "output_vcf: ${output_vcf}"
echo ""

# Select genes
bcftools view \
  "${source_vcf}" \
  --include 'vep_SYMBOL="MLH1" | vep_SYMBOL="MSH2" | vep_SYMBOL="MSH6" | vep_SYMBOL="PMS2" | vep_SYMBOL="MSH3" | vep_SYMBOL="PMS1" | vep_SYMBOL="EPCAM" | vep_SYMBOL="CDH1" | vep_SYMBOL="MAP3K6" | vep_SYMBOL="MYD88" | vep_SYMBOL="CTNNA1"' \
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
