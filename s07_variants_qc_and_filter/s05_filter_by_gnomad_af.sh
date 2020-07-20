#!/bin/bash

# s05_filter_by_gnomad_af.sh

# Started: Alexey  Larionov 18Jun2020
# Last updated: Rofaida Desoki 08Jul2020

# Filtering by gnomad AF is included here mainly for illustration: the later filtering by known/predicted effects (ClinVar,Lof<FIM) 
# will actually select variants that will be much more rare than AF 0.05  
# If this filtering was meaningful, then the typicaly would have used AF thresholds around 0.01-0.001 

#BSUB -q inter
#BSUB -P re_gecip_inherited_cancer_predisposition
#BSUB -o s05_filter_by_gnomad_af_%J.out
#BSUB -e s05_filter_by_gnomad_af_%J.err
#BSUB -n 4
#BSUB -R "rusage[mem=8000]"
#BSUB -R "span[hosts=1]"
##BSUB -W 01:00

# Start message
echo "Filter vcf by vep_gnomAD_NFE_AF"
date
echo ""

# Stop at runtime errors
set -e

# Load bcftools module
#'module load bcftools/1.10.2
module load bio/BCFtools/1.10.2-foss-2018b

# Files and folders
base_folder="/re_gecip/inherited_cancer_predisposition/rofaida/Rproject"
data_folder="${base_folder}/data/s07_variants_qc_and_filter"
scripts_folder="${base_folder}/scripts/s07_variants_qc_and_filter"
cd "${scripts_folder}" 

source_vcf="${data_folder}/s04_sort_tag_svep_clinvar_id_vkey_fix_qual_hwe.vcf.gz"
output_vcf="${data_folder}/s05_sort_tag_svep_clinvar_id_vkey_fix_qual_hwe_gnomad.vcf.gz"
log="${data_folder}/s05_sort_tag_svep_clinvar_id_vkey_fix_qual_hwe_gnomad.log"

# Progress report
echo "source_vcf: ${source_vcf}"
echo "output_vcf: ${output_vcf}"
echo ""

# Filter by vep_gnomAD_NFE_AF 
# Note that MAF is not AF
# Beware of many missed gnomad AF-s
bcftools view \
  "${source_vcf}" \
  --exclude "vep_gnomAD_NFE_AF > 0.05 & vep_gnomAD_NFE_AF < 0.95 " \
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
