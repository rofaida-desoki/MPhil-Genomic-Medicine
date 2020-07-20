# f01_explore_annotations.sh

# Function to print information about selected annotation fields in VCF file

# Started: Alexey Larionov 18Jun2020
# Last updated: Rofaida Desoki 25Jun2020

# Notes:
# 1) Requires bcftools and expects the name of VCF file as the parameter
# 2) The function does not check the existance of VCF file, availability of annotations etc

explore_annotations(){

	# Read parameter(s)
	local vcf="${1}"
	
	# Start message
	echo "VCF counts and annotations summary"
	date
	echo ""
	echo "${vcf}"
	echo ""

	# Print counts
	echo "--- Counts ---"
	echo ""
	bcftools +counts "${vcf}"
	echo ""

	# Extract selected annotations
	local annotations="${vcf%.vcf.gz}_annotations.txt"
	bcftools query -H -f '%ID\t%vep_SYMBOL\t%vep_Consequence\t%vep_IMPACT\t%CLNSIG\t%CLNREVSTAT\t%CLNDN\t%vep_SIFT\t%vep_PolyPhen\t%vep_CADD_PHRED\t%vep_gnomAD_NFE_AF\t%HWE\t%ALLELEID\t%vep_Existing_variation\t%QUAL\n' "${vcf}" > "${annotations}"

	echo "--- vep_SYMBOL ---"
	echo ""
	awk 'NR>1 {print $2}' "${annotations}" | sort |  uniq -c | sort -nr
	echo ""

	echo "--- CLNSIG CLNREVSTAT CLNDN ---"
	echo ""
	echo "Split CLNSIG"
	echo ""
	awk 'NR>1 {print $5}' "${annotations}" |  tr "|" "\n" | sort |  uniq -c | sort -nr
	echo ""
	echo "CLNREVSTAT"
	echo ""
	awk 'NR>1 {print $6}' "${annotations}" | sort |  uniq -c | sort -nr
	echo ""
	echo "CLNDN"
	echo ""
	echo "Count of distinct diagnoses:"
	awk 'NR>1 {print $7}' "${annotations}" | sort |  uniq -c | wc -l
	echo ""
	awk 'NR>1 {print $7}' "${annotations}" | sort |  uniq -c | sort -nr | head
	echo "..."
	echo ""
	echo "--- vep_IMPACT_Conesquence ---"
	echo ""
	echo "vep_IMPACT"
	echo ""
	awk 'NR>1 {print $4}' "${annotations}" | sort |  uniq -c | sort -nr
	echo ""
	echo "Split vep_Consequence"
	echo ""
	awk 'NR>1 {print $3}' "${annotations}" |  tr "@" "\n" | sort |  uniq -c | sort -nr
	echo ""

	echo "--- vep_SIFT_PplyPhen ---"
	echo ""
	echo "vep_SIFT"
	echo ""
	awk 'NR>1 {print $8}' "${annotations}" | sort |  uniq -c | sort -nr | head
	echo "..."
	echo ""
	echo "Split vep_SIFT" # after removal of the values in brackets
	echo ""
	bcftools query -f '%vep_SIFT\n' "${vcf}" | sed 's/(.*)//' | sort |  uniq -c | sort -nr
	echo ""
	echo "vep_PolyPhen"
	echo ""
	awk 'NR>1 {print $9}' "${annotations}" | sort |  uniq -c | sort -nr | head
	echo "..."
	echo ""
	echo "Split vep_PolyPhen" # after removal of the values in brackets
	echo ""
	bcftools query -f '%vep_PolyPhen\n' "${vcf}" | sed 's/(.*)//' | sort |  uniq -c | sort -nr
	echo ""

	echo "--- vep_CADD_PHRED ---"
	echo ""
	echo "Total count of variants:"
	awk 'NR>1 {print}' "${annotations}" | wc -l
	echo ""
	echo "Missed CADD_PHRED scores:"
	awk 'NR>1 {print $10}' "${annotations}" | grep -xc "\."
	echo ""
	echo "Non-Missed CADD_PHRED scores:"
	awk 'NR>1 {print $10}' "${annotations}" | grep -vxc "\."
	echo ""
	awk 'NR>1 {print $10}' "${annotations}" | grep -vx "\."  | sort -g | head
	echo "..."
	awk 'NR>1 {print $10}' "${annotations}" | grep -vx "\." | sort -g | tail
	echo ""

	echo "--- vep_gnomAD_NFE_AF ---"
	echo ""
	echo "Total count of variants:"
	awk 'NR>1 {print}' "${annotations}" | wc -l
	echo ""
	echo "Missed NFE AFs:"
	awk 'NR>1 {print $11}' "${annotations}" | grep -xc "\."
	echo ""
	echo "Zero NFE AFs:"
	awk 'NR>1 {print $11}' "${annotations}" | grep -xc 0
	echo ""
	echo "Non-missed non-zero NFE AFs:"
	awk 'NR>1 {print $11}' "${annotations}" | grep -vxc [0\.]
	echo ""
	awk 'NR>1 {print $11}' "${annotations}" | grep -vx [0\.]  | sort -g | head
	echo "..."
	awk 'NR>1 {print $11}' "${annotations}" | grep -vx [0\.] | sort -g | tail
	echo ""

	echo "--- INFO fields ---"
	echo ""
	echo "Count of INFO fields:"
	bcftools view -h "${vcf}" | grep ^##INFO | wc -l
	echo ""
	echo "List of INFO fields:"
	bcftools view -h "${vcf}" | grep ^##INFO
	echo ""

	# Completion message
	echo "Done"
	date
}
