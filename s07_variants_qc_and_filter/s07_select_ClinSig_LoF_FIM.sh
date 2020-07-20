#!/bin/bash

# s07_select_ClinSig_LoF_FIM.sh

# Started: Alexey Larionov 21Jun2020
# Last updated: Rofaida Desoki 09Jul2020

#BSUB -q inter
#BSUB -P re_gecip_inherited_cancer_predisposition
#BSUB -o s07_select_ClinSig_LoF_FIM_%J.out
#BSUB -e s07_select_ClinSig_LoF_FIM_%J.err
#BSUB -n 4
#BSUB -R "rusage[mem=8000]"
#BSUB -R "span[hosts=1]"
##BSUB -W 01:00

# Start message
echo "Select pathogenic in ClinSig + predicted LoF and FIM"
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
source_vcf="${data_folder}/s06_genes.vcf.gz" 
output_vcf="${data_folder}/s07_genes_ClnVar_LoF_FIM.vcf.gz" 
log="${data_folder}/s07_genes_ClnVar_LoF_FIM.log"

echo "source_vcf: ${source_vcf}"
echo "output_vcf: ${output_vcf}"
echo ""

# --- ClinSig --- #

# Pathogenic in ClinSig
path_CLNSIG="( CLNSIG~'Pathogenic' | CLNSIG~'Likely_pathogenic' | CLNSIG~'risk_factor' )"

# Good revision status in ClinSig
good_CLNREVSTAT="( CLNREVSTAT='reviewed_by_expert_panel' | CLNREVSTAT='criteria_provided\,_multiple_submitters\,_no_conflicts' | CLNREVSTAT='criteria_provided\,_single_submitter' )"

# Pathogenic in ClinSig + Good revision status in ClinSig
ClinSig_pathogenic="( $path_CLNSIG & $good_CLNREVSTAT )"

# --- LoF --- #

# Predicted Loss of Fuaction
vep_LoF="vep_IMPACT='HIGH'" 

# --- FIM --- #

# Functionally important missense variants
# SIFT + PolyPen + CADD
SIFT_PolyPhen_CADD="( vep_SIFT~'deleterious(' & vep_PolyPhen~'probably_damaging(' & vep_CADD_PHRED>=20 )"

# --- Inframe indels --- #

# inframe_deletion (~10) inframe_insertion (~4)
inframe_CADD="( ( vep_Consequence~'inframe_deletion' | vep_Consequence~'inframe_insertion' ) & vep_CADD_PHRED>=20 )"

# --- Combine filters --- #
combined_filters="$ClinSig_pathogenic | $vep_LoF | $SIFT_PolyPhen_CADD | $inframe_CADD"

# Progress report
echo "combined_filters: ${combined_filters}"
echo ""

# --- Select variants --- #
bcftools view "${source_vcf}" \
  --include "${combined_filters}" \
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
