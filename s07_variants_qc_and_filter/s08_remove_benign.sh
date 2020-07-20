#!/bin/bash

# s08_remove_benign.sh

# Known benign variants could be between predicted LoF and FIMs

# Started: Alexey  Larionov 18Jun2020
# Last updated: Rofaida Desoki 09Jul2020

#BSUB -q inter
#BSUB -P re_gecip_inherited_cancer_predisposition
#BSUB -o s08_remove_benign_%J.out
#BSUB -e s08_remove_benign_%J.err
#BSUB -n 4
#BSUB -R "rusage[mem=8000]"
#BSUB -R "span[hosts=1]"
##BSUB -W 01:00

# Start message
echo "Remove known benign variants"
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
source_vcf="${data_folder}/s07_genes_ClnVar_LoF_FIM.vcf.gz" 
output_vcf="${data_folder}/s08_genes_ClnVar_LoF_FIM_nonBenign.vcf.gz" 
log="${data_folder}/s08_genes_ClnVar_LoF_FIM_nonBenign.log"

# Progress report
echo "source_vcf: ${source_vcf}"
echo "output_vcf: ${output_vcf}"
echo ""

# Filterig criteria: benign + good revision status in ClinSig
bn_CLNSIG="( CLNSIG~'Benign' | CLNSIG~'Likely_benign' )"
good_CLNREVSTAT="( CLNREVSTAT='reviewed_by_expert_panel' | CLNREVSTAT='criteria_provided\,_multiple_submitters\,_no_conflicts' | CLNREVSTAT='criteria_provided\,_single_submitter' )"
ClinSig_benign="$bn_CLNSIG & $good_CLNREVSTAT"

# Progress report
echo "ClinSig_benign: ${ClinSig_benign}"
echo ""

# Selecting variants

bcftools view "${source_vcf}" \
  --exclude "${ClinSig_benign}" \
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
