#!/bin/bash

# s01_bcftools_stat_before_filters.sh

# Evaluate annotated VCF file with bcftools stat before filtering

# Started: Alexey  Larionov 24May2020
# Updated: Alexey Larionov 18Jun2020
# Last updated: Rofaida Desoki 29Jun2020

#BSUB -q inter
#BSUB -P re_gecip_inherited_cancer_predisposition
#BSUB -o s01_bcftools_stat_before_filters_%J.out
#BSUB -e s01_bcftools_stat_before_filters_%J.err
#BSUB -n 4
#BSUB -R "rusage[mem=8000]"
#BSUB -R "span[hosts=1]"
##BSUB -W 01:00

# Start message
echo "Evaluate annotated VCF file before filtering"
date
echo ""

# Stop at runtime errors
set -e

# Load modules
#module load bcftools/1.10.2
#module load texlive/2018 # for pdf output

module load bio/BCFtools/1.10.2-foss-2018b
module load vis/matplotlib/3.0.2-foss-2018b-Python-3.6.6
# texlive ? ?

# Files and folders
base_folder="/re_gecip/inherited_cancer_predisposition/rofaida/Rproject"
stats_folder="${base_folder}/data/s07_variants_qc_and_filter/s01_stats_before_filters"
scripts_folder="${base_folder}/scripts/s07_variants_qc_and_filter"
cd "${scripts_folder}" 

rm -fr "${stats_folder}" # remove folder if existed
mkdir -p "${stats_folder}"

vcf="${base_folder}/data/s06_annotate_vcf/s06_sort_tag_svep_clinvar_id_vkey_fix.vcf.gz"
stats="${stats_folder}/bcfstats.vchk"

# Progress report
echo "vcf: ${vcf}"
echo "stats_folder: ${stats_folder}"
echo ""

# Calculate bcfstats
echo "Calculating bcfstats..."
echo ""
bcftools stats -s- "${vcf}" > "{stats}"

# More vcf stats
echo "Calculating stats..."
bcftools stats -s- "${vcf_file}" > "${stats_file}" 

# Plot the stats ( PDF requires texlive)
echo "Making plots..."
#plot-vcfstats -p "${stats_folder}" "${stats_file}"
plot-vcfstats --no-PDF -p "${stats_folder}" "${stats}" #for Helix (until texlive module is available)
echo ""
   
# Progress report
echo "Done"
date

