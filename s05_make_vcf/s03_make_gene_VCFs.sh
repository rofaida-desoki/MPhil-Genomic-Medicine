#!/bin/bash

# s03_make_gene_VCFs.sh

# Summary
# Loop through the file with genes and chunks:
# - for each chunk select the gene region for all 60k samples, 
# - then keep only 13k selcted samples, 
# - then keep only variants PASS-ed filters and present in the selected samples (AC>0)
#
# It could be possible to combine several bcftools calls in one.  However, I split them in steps for clarity 
# (although, sometimes selecting smples does not work well when combined with filtering)

# Started: Alexey Larionov 29Mar2020
# Last updated: Rofaida Desoki 16Apr2020

# Use on cluster:
# bsub < s03_make_gene_VCFs.sh

# BSUB instructions:

#BSUB -q long
#BSUB -P re_gecip_inherited_cancer_predisposition
#BSUB -o s03_make_gene_VCFs_%J.out
#BSUB -e s03_make_gene_VCFs_%J.err
#BSUB -n 4
#BSUB -R "rusage[mem=8000]"
#BSUB -R "span[hosts=1]"
#BSUB -W 99:00

# Start message
echo "Extract data for required genes and samples from the source VCFs"
date
echo ""

# Stop at runtime errors
set -e

# Load bcftools module
module load bcftools/1.9

# Files and folders
base_folder="/re_gecip/inherited_cancer_predisposition/rofaida/Rproject"
samples_file="${base_folder}/data/s03_select_sample/samples.txt"

multisample_vcf_folder="/gel_data_resources/main_programme/aggregated_illumina_gvcf/GRCH38/20190228/data"

output_folder="${base_folder}/data/s05_make_vcf"
genes_file="${output_folder}/s02_selected_chunks.txt"

list_of_files="${output_folder}/s03_list_of_files.txt"
> "${list_of_files}" # delete file content as if any #one arrow ">" to overwrite

# Progress report
echo "genes_file: ${genes_file}"
echo "samples_file: ${samples_file}"
echo ""
echo "Started extracting data ..."
echo ""

# Read file with selected chunks line by line  
while  read chr start end gene file_suffix
do
  
  # Source chunk
  source_file="${multisample_vcf_folder}/60k_GRCH38_germline_mergedgVCF_${file_suffix}.bcf"
    
  # Select gene region for all 60k samples
  region="${chr}:${start}-${end}"
  gene_60k_raw_bcf="${output_folder}/s03_${gene}_60k_raw.bcf"
  gene_60k_raw_log="${output_folder}/s03_${gene}_60k_raw.log"
  
  bcftools view \
    "${source_file}" \
    --regions "${region}" \
    --output-file "${gene_60k_raw_bcf}" \
    --output-type b \
    --threads 4 \
    &> "${gene_60k_raw_log}"
  
  bcftools index "${gene_60k_raw_bcf}"
  
  # Filter samples: keep only selected 13k samples
  gene_13k_raw_bcf="${output_folder}/s03_${gene}_13k_raw.bcf"
  gene_13k_raw_log="${output_folder}/s03_${gene}_13k_raw.log"
  
  bcftools view \
    "${gene_60k_raw_bcf}" \
    --samples-file "${samples_file}" \
    --force-samples \
    --output-file "${gene_13k_raw_bcf}" \
    --output-type b \
    --threads 4 \
    &> "${gene_13k_raw_log}"
  
  bcftools index "${gene_13k_raw_bcf}"
  
  # Update AC and AN in INFO field (just in case ...)
  # Note that some bcftools functions use "--output" instead of "--output-file"
  gene_13k_acan_bcf="${output_folder}/s03_${gene}_13k_acan.bcf"
  gene_13k_acan_log="${output_folder}/s03_${gene}_13k_acan.log"
    
  bcftools plugin fill-AN-AC \
    "${gene_13k_raw_bcf}" \
    --output "${gene_13k_acan_bcf}" \
    --output-type b \
    --threads 4 \
    &> "${gene_13k_acan_log}"

  bcftools index "${gene_13k_acan_bcf}"

  # Filter variants: keep only PASS variants with AC>0 in selected samples
  gene_13k_filt_bcf="${output_folder}/s03_${gene}_13k_filt.bcf"
  gene_13k_filt_log="${output_folder}/s03_${gene}_13k_filt.log"
  
  bcftools view \
    "${gene_13k_acan_bcf}" \
    --exclude 'AC=0' \
    --apply-filters 'PASS' \
    --output-file "${gene_13k_filt_bcf}" \
    --output-type b \
    --threads 4 \
    &> "${gene_13k_filt_log}"
  
  bcftools index "${gene_13k_filt_bcf}"
  
  # Progress report
  echo "----- ${gene} -----"
  echo ""
  date
  echo "region: ${region}"
  echo "source_file: ${source_file}"
  echo "gene_13k_filt_bcf: ${gene_13k_filt_bcf}"
  echo ""
  bcftools view "${gene_13k_filt_bcf}" | tail | awk '{print $1,$2,$3,$4,$5,$6,$7,$8,$9,$10}'
  echo ""
  
  # Add file to the list of bcf files
  echo "${gene_13k_filt_bcf}" >> "${list_of_files}"
  
  # Remove intermediate files #"${gene_60k_raw_log}" & "${gene_13k_raw_log}" & "${gene_13k_filt_log}" kept to view
  rm "${gene_60k_raw_bcf}" "${gene_60k_raw_bcf}.csi" \
    "${gene_13k_raw_bcf}" "${gene_13k_raw_bcf}.csi" \
    "${gene_13k_acan_bcf}" "${gene_13k_acan_log}" "${gene_13k_acan_bcf}.csi" 
    
  
done < "${genes_file}" # Next gene

# Completion message
echo "Done all genes"
date
