#!/bin/bash

# Annotate with VEP
# s02_annotate_with_vep.sh

# Started: Alexey Larionov, 22Apr2020
# Last updated: Rofaida Desoki, 04Jun2020

# Intended use:
# bsub < s02_annotate_with_vep.sh

# Note:
# A template provided by GEL here:
# /gel_data_resources/example_scripts/annotate_variants_with_vep

#BSUB -q inter
#BSUB -P re_gecip_inherited_cancer_predisposition
#BSUB -o s02_annotate_with_vep_%J.stdout
#BSUB -e s02_annotate_with_vep_%J.stderr
#BSUB -R "rusage[mem=16000]"
#BSUB -R "span[hosts=1]"
#BSUB -n 4
##BSUB -W 01:00

# Stop at runtime errors
set -e

# Start message
echo "Annotate with VEP"
date
echo ""

# Load required module (note bcftools version)
module load vep/98
module load bcftools/1.10.1
echo ""

# Folders
base_folder="/re_gecip/inherited_cancer_predisposition/rofaida/Rproject"
data_folder="${base_folder}/data/s06_annotate_vcf"
scripts_folder="${base_folder}/scripts/s06_annotate_vcf"
cd "${scripts_folder}"

# Reference genome
b38_fasta="/public_data_resources/reference/GRCh38/GCA_000001405.15_GRCh38_no_alt_analysis_set.fna"

# Source and ouptput (note the file formats)
source_vcf="${data_folder}/s01_sort_tag.vcf.gz"
vep_vcf="${data_folder}/s02_sort_tag_vep.vcf.gz"

# Log and report
vep_log="${data_folder}/s02_vep_annotation.log"
vep_report="${data_folder}/s02_vep_report.html"

# VEP cache version and folder
cache_version="98"
cache_folder="/tools/apps/vep/98/ensembl-vep/.vep"

# Plugins folder
plugins_folder="/tools/apps/vep/plugin/VEP_plugins"

# CADD plugin data #these are sources for cadd annotation plugin to use
cadd_data_folder="/public_data_resources/CADD/v1.5/GRCh38"
cadd_snv_data="${cadd_data_folder}/whole_genome_SNVs.tsv.gz"
cadd_indels_data="${cadd_data_folder}/InDels.tsv.gz"

# Num of threads
n_threads="4"

# File to copy information about VEP cache
vep_cache_info="${data_folder}/s02_vep_cache_info.txt"

# Progress report
echo "--- Settings ---"
echo ""
echo "source_vcf: ${source_vcf}"
echo "vep_vcf: ${vep_vcf}"
echo "vep_log: ${vep_log}"
echo "vep_report: ${vep_report}"
echo ""
echo "n_threads: ${n_threads}"
echo ""
echo "--- Data sources ---"
echo ""
echo "cache_folder: ${cache_folder}"
echo "cache_version: ${cache_version}"
echo ""
echo "b38_fasta: ${b38_fasta}"
echo ""
echo "See used VEP cache description in the following file:"
echo "${vep_cache_info}"
echo ""
echo "CADD annotation files:"
echo "${cadd_snv_data}"
echo "${cadd_indels_data}"
echo ""
echo "Working..."

# Keep information about used VEP cache #copy used versions of different program into the new cache_info file
cp "${cache_folder}/homo_sapiens/98_GRCh38/info.txt" "${vep_cache_info}"

# Annotate VCF
vep \
  --input_file "${source_vcf}" \
  --output_file "${vep_vcf}" \
  --vcf \
  --force_overwrite \
  --compress_output bgzip \
  --stats_file "${vep_report}" \
  --fork "${n_threads}" \
  --offline \
  --cache \
  --dir_cache "${cache_folder}" \
  --cache_version "${cache_version}" \
  --species homo_sapiens \
  --assembly GRCh38 \
  --fasta "${b38_fasta}" \
  --check_ref \
  --everything \
  --total_length \
  --check_existing \
  --exclude_null_alleles \
  --pick \
  --gencode_basic \
  --dir_plugins "${plugins_folder}" \
  --plugin CADD,"${cadd_snv_data}","${cadd_indels_data}" \
  &> "${vep_log}"

#  --nearest symbol \
# Can not be used yet: needs cache update at the first run

# Index annotated vcf
bcftools index "${vep_vcf}" 

# List added VEP annotations
echo ""
echo "Added VEP annotations:"
echo ""
bcftools +split-vep -l "${vep_vcf}"
echo ""

# Summary of INFO fields
echo "Number of INFO fields in vep_vcf:"
bcftools view -h "${vep_vcf}" | grep ^##INFO | wc -l
echo ""

echo "List of INFO fields in vep_vcf:"
bcftools view -h "${vep_vcf}" | grep ^##INFO
echo ""

# Progress report
echo "Example of record in source_vcf:"
echo ""
bcftools view "${source_vcf}" | tail  -n 1 | awk '{print $1,$2,$3,$4,$5,$6,$7,$8,$9,$10}'
echo ""
echo "Example of record in vep_vcf:"
echo ""
bcftools view "${vep_vcf}" | tail  -n 1 | awk '{print $1,$2,$3,$4,$5,$6,$7,$8,$9,$10}'
echo ""

# Completion message
echo "Done"
date
