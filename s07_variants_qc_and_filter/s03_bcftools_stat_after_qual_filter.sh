#!/bin/bash

# s03_bcftools_stat_after_filters.sh

# Evaluate VCF file with bcftools stat

# Started: Alexey  Larionov 24May2020
# Updated: Alexey 18Jun2020
# Last updated: Rofaida Desoki 06Jul2020

#BSUB -q inter
#BSUB -P re_gecip_inherited_cancer_predisposition
#BSUB -o s03_bcftools_stat_after_qual_filter_%J.out
#BSUB -e s03_bcftools_stat_after_qual_filter_%J.err
#BSUB -n 4
#BSUB -R "rusage[mem=8000]"
#BSUB -R "span[hosts=1]"
##BSUB -W 01:00

# Start message
echo "Evaluate VCF after filtering by QUAL"
date
echo ""

# Stop at runtime errors
set -e

# Load modules
#module load bcftools/1.10.1
#module load texlive/2018 # for pdf output

module load bio/BCFtools/1.10.2-foss-2018b
module load vis/matplotlib/3.0.2-foss-2018b-Python-3.6.6
# texlive ? ?

# Files and folders
base_folder="/re_gecip/inherited_cancer_predisposition/rofaida/Rproject"
data_folder="${base_folder}/data/s07_variants_qc_and_filter/"
stats_folder="${data_folder}/s03_bcfstats_after_qual_filter"
scripts_folder="${base_folder}/scripts/s07_variants_qc_and_filter"
cd "${scripts_folder}" 

rm -fr "${stats_folder}" # remove folder if existed
mkdir -p "${stats_folder}"

vcf_file="${data_folder}/s02_sort_tag_svep_clinvar_id_vkey_fix_qual.vcf.gz"
stats_file="${stats_folder}/bcfstats.vchk"

# Progress report
echo "vcf_file: ${vcf_file}"
echo "stats_folder: ${stats_folder}"
echo ""

# Basic vcf counts
echo "VCF counts:"
echo ""
bcftools plugin counts "${vcf_file}" 
echo ""

# More vcf stats
echo "Calculating bcfstats..."
bcftools stats -s- "${vcf_file}" > "${stats_file}" 

# Plot the stats - requires texlive
echo "Making plots..."
#plot-vcfstats -p "${stats_folder}" "${stats_file}"
plot-vcfstats --no-PDF -p "${stats_folder}" "${stats_file}" # for Helix until texlive module is available
echo ""
   
# Progress report
echo "Done"
date

