#!/bin/bash

##Requires bedtools + SURVIVOR 

min_args=2

if [ $# -eq $min_args ]; then        
        if [ ! -f $1 ]; then
                echo "file $1 not exists"
        elif [ ! -f $2 ]; then
                echo "file $2 not exists"
	else
		svs=$1
		snp=$2
		zcat $svs > $svs'.tmp.vcf'
		~/workspace/SURVIVOR/Debug/SURVIVOR vcftobed $svs'.tmp.vcf' -1 -1 $svs'.bedpe'
		cat $svs'.bedpe' | grep -v 'TRA' | cut -f 1,2,5,7,11 > $svs'.bed'
		rm  $svs'.bedpe'
		bedtools intersect -wo -a $snp -b $svs'.bed'  > $svs'_intersectSNP.tab'
       fi

else
        echo "Intersects SVs with SNP for subprocessing"
        echo "File path of SVs VCF.GZ file"
	echo "File path of SNP VCF.GZ file"
fi


