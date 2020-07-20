#!/bin/bash

# s02_filter_by_qual.sh

# Started: Alexey  Larionov 24May2020
# Updated: Alexey Larionov 19Jun2020
# Last updated: Rofaida Desoki 03Jul2020

# The SNP's QUAL threshold of 641 approximately corresponds to Ts/Tv~2 
#(as shown in the bcftools-stat ts-tv by QUAL plot)
# Arbitrarily, the same threshold (641) is selected for INDELs

# Filtering by call rate is NOT applied at this stage because it may slightly change after selecting specific cases & controls.
# So, this filter (call rate > 0.85) will be applied immediately before SKAT analysis.  

#BSUB -q inter
#BSUB -P re_gecip_inherited_cancer_predisposition
#BSUB -o s02_filter_by_qual_%J.out
#BSUB -e s02_filter_by_qual_%J.err
#BSUB -n 4
#BSUB -R "rusage[mem=8000]"
#BSUB -R "span[hosts=1]"
##BSUB -W 01:00

# Start message
echo "Filter SNPs and INDELs in VCF file by QUAL"
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

source_vcf="${base_folder}/data/s06_annotate_vcf/s08_sort_tag_svep_clinvar_id_vkey_fix.vcf.gz"
output_vcf="${data_folder}/s02_sort_tag_svep_clinvar_id_vkey_fix_qual.vcf.gz"

# Folder for temporary files
tmp_folder="${data_folder}/tmp"
rm -fr "${tmp_folder}" # remove folder if existed
mkdir "${tmp_folder}"

# Temporary files and logs 
snp_vcf="${tmp_folder}/s02_snp.vcf.gz"
nonsnp_vcf="${tmp_folder}/s02_nonsnp.vcf.gz"

snp_qual_vcf="${tmp_folder}/s02_snp_qual.vcf.gz"
nonsnp_qual_vcf="${tmp_folder}/s02_nonsnp_qual.vcf.gz"

concat_vcf="${tmp_folder}/s02_concat.vcf.gz"

snp_log="${tmp_folder}/s02_snp.log"
nonsnp_log="${tmp_folder}/s02_nonsnp.log"

snp_qual_log="${tmp_folder}/s02_snp_qual.log"
nonsnp_qual_log="${tmp_folder}/s02_nonsnp_qual.log"

concat_log="${tmp_folder}/s02_concat.log"
sort_log="${tmp_folder}/s02_sort.log"

# Progress report
echo "source_vcf: ${source_vcf}"
echo "output_vcf: ${output_vcf}"
echo ""
echo "Counts in source vcf:"
echo ""
bcftools +counts "${source_vcf}"
echo ""

# Select SNPs

bcftools view \
  "${source_vcf}" \
  --include "TYPE='snp'" \
  --output-type z \
  --threads 4 \
  --output-file "${snp_vcf}" \
  &> "${snp_log}"

# Make index
bcftools index "${snp_vcf}"

echo "SNP vcf counts:"
echo ""
bcftools +counts "${snp_vcf}"
echo ""

# Select non-SNPs

bcftools view \
  "${source_vcf}" \
  --include "TYPE!='snp'" \
  --output-type z \
  --threads 4 \
  --output-file "${nonsnp_vcf}" \
  &> "${nonsnp_log}"

# Make index
bcftools index "${nonsnp_vcf}"

echo "Non-SNP vcf counts:"
echo ""
bcftools +counts "${nonsnp_vcf}"
echo ""

# Filter SNP-s

bcftools view \
  "${snp_vcf}" \
  --include "QUAL >= 641" \
  --output-type z \
  --threads 4 \
  --output-file "${snp_qual_vcf}" \
  &> "${snp_qual_log}"

bcftools index "${snp_qual_vcf}"

echo "Filtered SNP counts:"
echo ""
bcftools +counts "${snp_qual_vcf}"
echo ""

# Filter INDEL-s : arbitrary choosing the same QUAL threshold as for SNP-s !!!

bcftools view \
  "${nonsnp_vcf}" \
  --include "QUAL >= 641" \
  --output-type z \
  --threads 4 \
  --output-file "${nonsnp_qual_vcf}" \
  &> "${nonsnp_qual_log}"

bcftools index "${nonsnp_qual_vcf}"

echo "Filtered non-SNP counts:"
echo ""
bcftools +counts "${nonsnp_qual_vcf}"
echo ""

# Concatenate non-SNPs with filtered SNPs

bcftools concat \
    "${snp_qual_vcf}" "${nonsnp_qual_vcf}" \
    --allow-overlaps \
    --output-type z \
    --output "${concat_vcf}" \
    &> "${concat_log}"

bcftools sort  \
  "${concat_vcf}" \
  --output-type z \
  --output-file "${output_vcf}" \
  &> "${sort_log}"

bcftools index "${output_vcf}" 

echo "Output vcf counts:"
echo ""
bcftools +counts "${output_vcf}"
echo ""

# Get detailed annnotations summary
source "${scripts_folder}/f01_explore_annotations.sh"
summary="${output_vcf%.vcf.gz}_annotations_summary.txt" 
explore_annotations "${output_vcf}" > "${summary}"

# Remove temporary files
rm -fr "${tmp_folder}"

# Progress report
echo "Done"
date
