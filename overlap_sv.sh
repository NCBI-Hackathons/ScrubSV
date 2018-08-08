#!/usr/bin/env bash

################ Requirements: Installation of vcftools, R, and GenomicRanges (R package)

###### On DNANexus worker
### vcftools:   aps install vcftools
### R is installed already, just installed "GenomicRanges" 
###             source("https://bioconductor.org/biocLite.R")
###             biocLite("GenomicRanges")


################ Description ################
# This script detects overlays in SV falls and flags them in the output tabular


################  The script ################

min_args=1

if [ $# -ne $min_args ]; then
	echo "Specify only the SV VCF file"
	exit 1
else
        echo "Intersects SVs with SNP for subprocessing"
        echo "File path of SVs VCF.GZ file"
	echo "File path of SNP VCF.GZ file"
	echo "Output: Reports a tab formated SNP file splitted into DUP, DEL, INV."
fi

sv_vcf="${1}"
sv_vcf_name= "${sv_vcf%%.vcf}"
sv_vcf_name="${sv_vcf%%.vcf.gz}"
sv_vcf_name="${sv_vcf##*/}"
tabular_filename="${sv_vcf_name}.gz"

vcf-query -f '%ID\t%CHROM\t%POS\t%INFO/END\t%INFO/SVTYPE\n' "${sv_vcf}" > "${tabular_filename}"
Rscript overlap_sv.R --arg1="${tabular_filename}" --arg2=DEL
Rscript overlap_sv.R --arg1="${tabular_filename}" --arg2=DUP
