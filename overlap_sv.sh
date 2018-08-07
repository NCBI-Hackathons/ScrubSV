#!/usr/bin/env bash

module load vcftools
module load R

#first file
vcf-query -f '%ID\t%CHROM\t%POS\t%INFO/END\t%INFO/SVTYPE\n' data/H002.bwamem.ill.mapped.filter.sort.bam.lumpy_noalts.vcf > HG002_delly.tab
Rscript overlap_sv.R --arg1=HG002_delly.tab --arg2=DEL
Rscript overlap_sv.R --arg1=HG002_delly.tab --arg2=DUP
