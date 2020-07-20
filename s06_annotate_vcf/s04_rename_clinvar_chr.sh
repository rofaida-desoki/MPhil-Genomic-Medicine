#!/bin/bash

# s04_rename_clinvar_chr.sh

# Summary:
# Rename chromosomes in ClinVar using bcftools --rename-chrs function. 

# The used bcftools --rename-chrs function is not well documented, so it is not immediately clear what happens, 
# for instance, with some of the contigs not included in the list.  The headers' check suggests that data for such contigs are removed 
# (at least from the header).  It should be OK in our case, because we don't have variants outside of the canonical chromosomes.  

# Started: Alexey Larionov 22Apr2020
# Last updated: Rofaida Desoki 06Jun2020

#BSUB -q inter
#BSUB -P re_gecip_inherited_cancer_predisposition
#BSUB -o s04_rename_clinvar_chr_%J.out
#BSUB -e s04_rename_clinvar_chr_%J.err
#BSUB -n 4
#BSUB -R "rusage[mem=8000]"
#BSUB -R "span[hosts=1]"
##BSUB -W 01:00

# Stop at runtime errors
set -e

# Start message
echo "Rename chromosomes in ClinVar"
date
echo ""
echo "The script changes 1,2,3... style to chr1,chr2,chr3... style"
echo "This is necessary for compatibility of clinvar vcf with aggregated vcf provided by GEL"
echo ""

# Load required module
module load bcftools/1.9

# Files and folders
base_folder="/re_gecip/inherited_cancer_predisposition/rofaida/Rproject"
data_folder="${base_folder}/data/s06_annotate_vcf"
scripts_folder="${base_folder}/scripts/s06_annotate_vcf"
cd "${scripts_folder}"

public_clinvar_folder="/public_data_resources/clinvar/20200604/clinvar/vcf_GRCh38"
public_clinvar_vcf="${public_clinvar_folder}/clinvar_20200602.vcf.gz" # Note inconsistency in dates

custom_clinvar_folder="${base_folder}/resources/clinvar/20200604/GRCh38"
custom_clinvar_vcf="${custom_clinvar_folder}/clinvar_20200602.vcf.gz" 

chromosomes_translation_file="${scripts_folder}/s04_clinvar_chr_translation.txt"
log="${data_folder}/s04_rename_clinvar_chr.log"

gel_vcf="${data_folder}/s03_sort_tag_svep.vcf.gz" #to view and compare

# Make updated clinvar folder
mkdir -p "${custom_clinvar_folder}"

# Progress report
echo "--- Input and output files ---"
echo ""
echo "public_clinvar_vcf: ${public_clinvar_vcf}"
echo "custom_clinvar_vcf: ${custom_clinvar_vcf}"
echo "chromosomes_translation_file: ${chromosomes_translation_file}"
echo "gel_vcf: ${gel_vcf}"
echo ""

# Check headers from clinvar and gel vcf-s

echo "--- Contigs in public ClinVar VCF header ---"
echo ""
bcftools view -h "${public_clinvar_vcf}" | grep "^##contig"
echo ""

echo "--- Number of records in public ClinVar VCF ---"
echo "" #H = no header
bcftools view -H "${public_clinvar_vcf}" | wc -l
echo ""

echo "--- Contigs in GEL VCF ---"
echo ""
bcftools view -h "${gel_vcf}" | grep "^##contig" | head -n 27
echo "..."
echo ""

echo "--- Translation file ---"
echo ""
cat "${chromosomes_translation_file}"
echo ""

# Update ClinVar VCF

echo "Updating ClinVar VCF ..."

bcftools annotate \
  --rename-chrs "${chromosomes_translation_file}" \
  --output "${custom_clinvar_vcf}" \
  --output-type z \
  --threads 4 \
  "${public_clinvar_vcf}" \
  &> "${log}"

bcftools index "${custom_clinvar_vcf}"

# Progress report
echo "Done"
echo ""

echo "--- Contigs in the updated ClinVar VCF header ---"
echo ""
bcftools view -h "${custom_clinvar_vcf}" | grep "^##contig"
echo ""
echo "--- Number of records in updated ClinVar VCF ---"
echo ""
bcftools view -H "${custom_clinvar_vcf}" | wc -l
echo ""

# Completion message
echo "Completed script"
date
echo ""
