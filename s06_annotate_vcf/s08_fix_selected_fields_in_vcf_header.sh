#!/bin/bash

# s08_fix_selected_fields_in_vcf_header.sh
# Fix the types of data in the columns (vep_CADD_PHRED) and (vep_gnomAD_NFE_AF) to be float instead of string.

# Started: Alexey Larionov, 18Jun2020
# Last updated: Rofaida Desoki, 26Jun2020

#BSUB -q inter
#BSUB -P re_gecip_inherited_cancer_predisposition
#BSUB -o s08_fix_selected_fields_in_vcf_header_%J.stdout
#BSUB -e s08_fix_selected_fields_in_vcf_header_%J.stderr
#BSUB -R "rusage[mem=6000]"
#BSUB -n 1
##BSUB -W 01:00

# Stop at runtime errors
set -e

# Start message
echo "Fix types of selected annotations"
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

source_vcf="${data_folder}/s06_sort_tag_svep_clinvar_id_vkey.vcf.gz"
output_vcf="${data_folder}/s08_sort_tag_svep_clinvar_id_vkey_fix.vcf.gz"
header_tmp="${data_folder}/s08_header.tmp"

# Progress report
echo "source_vcf: ${source_vcf}"
echo "output_vcf: ${output_vcf}"
echo ""

# Extract header
bcftools view -h "${source_vcf}" > "${header_tmp}"

# Fix selected numeric fields recorded as "String"
##INFO=<ID=vep_CADD_PHRED,Number=.,Type=String,Description="The CADD_PHRED field from INFO/CSQ">
##INFO=<ID=vep_gnomAD_NFE_AF,Number=.,Type=String,Description="The gnomAD_NFE_AF field from INFO/CSQ">

sed -i 's/ID=vep_CADD_PHRED,Number=.,Type=String/ID=vep_CADD_PHRED,Number=.,Type=Float/' "${header_tmp}"
sed -i 's/ID=vep_gnomAD_NFE_AF,Number=.,Type=String/ID=vep_gnomAD_NFE_AF,Number=.,Type=Float/' "${header_tmp}"

# Update header in the file
bcftools reheader -h "${header_tmp}" "${source_vcf}" > "${output_vcf}"
bcftools index "${output_vcf}"

# Check result
summary="${output_vcf%.vcf.gz}_annnotations_summary.txt"
source "${scripts_folder}/f01_explore_annotations.sh"
explore_annotations "${output_vcf}" > "${summary}"

echo "PS: note Type in vep_CADD_PHRED and vep_gnomAD_NFE_AF" >>  "${summary}"

# Remove temporary file
rm "${header_tmp}"

# Completion message
echo "Done" 
date
