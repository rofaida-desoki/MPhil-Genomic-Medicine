#!/bin/bash

# s05_annotate_with_clinvar.sh

# Summary:
# Annotate with ClinVar using bcftools --annotate function. 

# Started: Alexey Larionov 22Apr2020
# Last updated: Rofaida Desoki 06Jun2020

#BSUB -q inter
#BSUB -P re_gecip_inherited_cancer_predisposition
#BSUB -o s05_annotate_with_clinvar_%J.out
#BSUB -e s05_annotate_with_clinvar_%J.err
#BSUB -n 4
#BSUB -R "rusage[mem=8000]"
#BSUB -R "span[hosts=1]"
##BSUB -W 01:00

# Stop at runtime errors
set -e

# Start message
echo "Annotate with ClinVar"
date
echo ""

# Load required module
module load bcftools/1.9
echo ""

# Files and folders
base_folder="/re_gecip/inherited_cancer_predisposition/rofaida/Rproject"
data_folder="${base_folder}/data/s06_annotate_vcf"
scripts_folder="${base_folder}/scripts/s06_annotate_vcf"
cd "${scripts_folder}"

source_vcf="${data_folder}/s03_sort_tag_svep.vcf.gz"
output_vcf="${data_folder}/s05_sort_tag_svep_clinvar.vcf.gz"

log="${data_folder}/s05_annotate_with_clinvar.log"

# ClinVar data (customised / updated contigs format: chr1 instead of 1)
clinvar_folder="${base_folder}/resources/clinvar/20200604/GRCh38"
clinvar_vcf="${clinvar_folder}/clinvar_20200602.vcf.gz" 

# Progress report
echo "--- Input and output files ---"
echo ""
echo "source_vcf: ${source_vcf}"
echo "output_vcf: ${output_vcf}"
echo "clinvar_vcf: ${clinvar_vcf}"
echo ""
echo "Working..."

# Annotate with bcftools
# http://www.htslib.org/doc/bcftools.html#annotate
# Adds INFO fields from clinvar
# (there was no fields with such name in the initial file) 
bcftools annotate \
  --annotations "${clinvar_vcf}" \
  --columns INFO \
  --output "${output_vcf}" \
  --output-type z \
  --threads 4 \
  "${source_vcf}" \
  &> "${log}"

# Index output vcf
bcftools index "${output_vcf}"

# Summary of INFO fields
echo ""
echo "Number of INFO fields in ClinVar-annotated vcf:"
bcftools view -h "${output_vcf}" | grep ^##INFO | wc -l
echo ""

echo "List of INFO fields in ClinVar-annotated vcf:"
bcftools view -h "${output_vcf}" | grep ^##INFO
echo ""

# Progress report
echo ""
echo "selected lines from output_vcf:"
echo "" # print some lines by --query to extract fields from bcf/vcf into user-frienlty format #-i include #-f format
bcftools query \
  -i 'ALLELEID != "."' \
  -f '%CHROM\t%POS\t%REF\t%ALT\t%ALLELEID\t%CLNSIG\t%CLNDN\n' "${output_vcf}" | head
echo ""

# Completion message
echo "Done"
date
echo ""

# --- Description of the added ClinVar annotations --- #

##ID=<Description="ClinVar Variation ID">

##INFO=<ID=ALLELEID,Number=1,Type=Integer,Description="the ClinVar Allele ID">
##INFO=<ID=CLNREVSTAT,Number=.,Type=String,Description="ClinVar review status for the Variation ID">

##INFO=<ID=CLNSIG,Number=.,Type=String,Description="Clinical significance for this single variant">
##INFO=<ID=CLNDN,Number=.,Type=String,Description="ClinVar's preferred disease name for the concept specified by disease identifiers in CLNDISDB">
##INFO=<ID=CLNDISDB,Number=.,Type=String,Description="Tag-value pairs of disease database name and identifier, e.g. OMIM:NNNNNN">

##INFO=<ID=CLNSIGCONF,Number=.,Type=String,Description="Conflicting clinical significance for this single variant">
##INFO=<ID=SSR,Number=1,Type=Integer,Description="Variant Suspect Reason Codes. One or more of the following values may be added: 0 - unspecified, 1 - Paralog, 2 - byEST, 4 - oldAlign, 8 - Para_EST, 16 - 1kg_failed, 1024 - other">

##INFO=<ID=CLNVI,Number=.,Type=String,Description="the variant's clinical sources reported as tag-value pairs of database and variant identifier">
##INFO=<ID=ORIGIN,Number=.,Type=String,Description="Allele origin. One or more of the following values may be added: 0 - unknown; 1 - germline; 2 - somatic; 4 - inherited; 8 - paternal; 16 - maternal; 32 - de-novo; 64 - biparental; 128 - uniparental; 256 - not-tested; 512 - tested-inconclusive; 1073741824 - other">

##INFO=<ID=GENEINFO,Number=1,Type=String,Description="Gene(s) for the variant reported as gene symbol:gene id. The gene symbol and id are delimited by a colon (:) and each pair is delimited by a vertical bar (|)">
##INFO=<ID=CLNVCSO,Number=1,Type=String,Description="Sequence Ontology id for variant type">
##INFO=<ID=CLNVC,Number=1,Type=String,Description="Variant type">
##INFO=<ID=MC,Number=.,Type=String,Description="comma separated list of molecular consequence in the form of Sequence Ontology ID|molecular_consequence">

##INFO=<ID=RS,Number=.,Type=String,Description="dbSNP ID (i.e. rs number)">
##INFO=<ID=CLNHGVS,Number=.,Type=String,Description="Top-level (primary assembly, alt, or patch) HGVS expression.">
##INFO=<ID=DBVARID,Number=.,Type=String,Description="nsv accessions from dbVar for the variant">

##INFO=<ID=AF_ESP,Number=1,Type=Float,Description="allele frequencies from GO-ESP">
##INFO=<ID=AF_EXAC,Number=1,Type=Float,Description="allele frequencies from ExAC">
##INFO=<ID=AF_TGP,Number=1,Type=Float,Description="allele frequencies from TGP">

##INFO=<ID=CLNSIGINCL,Number=.,Type=String,Description="Clinical significance for a haplotype or genotype that includes this variant. Reported as pairs of VariationID:clinical significance.">
##INFO=<ID=CLNDNINCL,Number=.,Type=String,Description="For included Variant : ClinVar's preferred disease name for the concept specified by disease identifiers in CLNDISDB">
##INFO=<ID=CLNDISDBINCL,Number=.,Type=String,Description="For included Variant: Tag-value pairs of disease database name and identifier, e.g. OMIM:NNNNNN">
