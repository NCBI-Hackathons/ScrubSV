#!/bin/bash

##Requires bedtools + SURVIVOR 

min_args=3

if [ $# -eq $min_args ]; then        
        if [ ! -f $1 ]; then
                echo "file $1 not exists"
        elif [ ! -f $2 ]; then
                echo "file $2 not exists"
	else
		svs=$1
		snp=$2
		reads=$3
		zcat $svs > $svs'.tmp.vcf'
		~/workspace/SURVIVOR/Debug/SURVIVOR vcftobed $svs'.tmp.vcf' -1 -1 $svs'.bedpe'
		cat $svs'.bedpe' | grep -v 'TRA' | cut -f 1,2,5,7,11 > $svs'.bed'
		rm  $svs'.bedpe'
		bedtools intersect -wo -a $snp -b $svs'.bed'  > $svs'_intersectSNP.tab'
		##awk '{print $1":"$2"-"$2+1}'  $svs'_intersectSNP.tab' | uniq  >  $svs'_intersectSNP.reg'
		samtools depth -q 20 -a -b $svs'.bed'  $reads | gzip - > $svs'_cov.txt.gz'
		rm  $svs'.bed'
		
		grep 'DUP' $svs'_intersectSNP.tab' > $svs'_intersectSNP.DUP.tab'
		grep 'DEL' $svs'_intersectSNP.tab' > $svs'_intersectSNP.DEL.tab'
		grep 'INV' $svs'_intersectSNP.tab' > $svs'_intersectSNP.INV.tab'
		rm  $svs'_intersectSNP.tab'
       fi

else
        echo "Intersects SVs with SNP for subprocessing"
        echo "File path of SVs VCF.GZ file"
	echo "File path of SNP VCF.GZ file"
	echo "Output: Reports a tab formated SNP file splitted into DUP, DEL, INV."
fi


