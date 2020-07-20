#!/bin/bash

# s03_make_plink_inuvika.sh

# Started: Alexey Larionov 14Jul2020
# Last updated: Rofaida Desoki 17Jul2020

# Note:
# On 14Jun2020 this script worked on Inuvika, but didnt work on Helix  

# Description:
# Make plink (bed-bim-fam) from vcf and fam

# Use:
# chmod +x s03_make_plink_inuvika.sh # Inuvika
# ./s03_make_plink_inuvika.sh &> s03_make_plink_inuvika.log # Inuvika
# bsub < s03_make_plink.sh # Helix

# BSUB instructions:

#BSUB -q inter
#BSUB -P re_gecip_inherited_cancer_predisposition
#BSUB -o s03_make_plink_%J.out
#BSUB -e s03_make_plink_%J.err
#BSUB -n 1
#BSUB -R "rusage[mem=1000]"
#BSUB -R "span[hosts=1]"
#BSUB -W 01:00

# Start message
echo "Make plink (bed-bim-fam) from vcf and fam"
date
echo ""

# Stop at runtime errors
set -e

# Load PLINK2 module
module load plink/2.0 # Inuvika
#module load bio/PLINK/2.00-devel-20200409-x86_64 # Helix
 
# Files and folders
base_folder="/home/rdesoki/re_gecip/inherited_cancer_predisposition/rofaida/Rproject/data/s09_sporadic_colon_cancer" # Inuvika
#base_folder="/re_gecip/inherited_cancer_predisposition/tischkowitz/users/alexey/ra_2020/data/s09_sporadic_colon_cancer" # Helix

vcf="${base_folder}/vcf_fam/selected.vcf.gz"
fam="${base_folder}/vcf_fam/selected.fam"

mkdir -p "${base_folder}/plink"
plink="${base_folder}/plink/dataset"

# Progress report
echo "vcf: ${vcf}"
echo "fam: ${fam}"
echo "plink: ${plink}"
echo ""

# Convert data
echo "Converting data ..."

plink2 \
  --vcf "${vcf}" \
  --fam "${fam}" \
  --vcf-half-call "m" \
  --double-id \
  --make-bed \
  --out "${plink}"

echo ""

# Progress report
echo "Done"
date
