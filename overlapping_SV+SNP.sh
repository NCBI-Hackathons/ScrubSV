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
		##for each CHR defined in the header!
		
		if [ -f $svs'_intersectSNP.DUP+cov.tab' ]; then
			echo "Res DUP file exists and is going to be deleted"
			rm $svs'_intersectSNP.DUP+cov.tab'
		fi

		 if [ -f $svs'_intersectSNP.INV+cov.tab' ]; then
                        echo "Res INV file exists and is going to be deleted"
                        rm $svs'_intersectSNP.INV+cov.tab'
                fi

		 if [ -f $svs'_intersectSNP.DEL.tab' ]; then
                        echo "Res file exists and is going to be deleted"
                        rm $svs'_intersectSNP.DEL.tab'
                fi

		for i in $(zcat $1 | grep 'ID=' |  cut -f 3 -d '='| grep 'length' | cut -f 1 -d ',' | less) 
		do
			echo "Processing chr: $i"
			echo "	Converting SVs vcf to bed file:"
			zcat $svs | grep "^$i" > $svs'.tmp.vcf'
			~/workspace/SURVIVOR/Debug/SURVIVOR vcftobed $svs'.tmp.vcf' -1 -1 $svs'.bedpe'
			cat $svs'.bedpe' | grep -v 'TRA' | cut -f 1,2,5,7,11 > $svs'.bed'
			rm  $svs'.bedpe'
		
			echo "	Split up into types"
			grep 'DUP' $svs'.bed' > $svs'.DUP.bed'
			grep 'INV' $svs'.bed' > $svs'.INV.bed'
			grep 'DEL' $svs'.bed' > $svs'.DEL.bed'
		
			echo "	Extract SNPs overlapping SVs and reformat into bed file for DUP and INV"
			bedtools intersect -wo -a $snp -b $svs'.DUP.bed' | awk '{print $1,$2,$2+1,$10,$11,$12,$13,$14,$15}' | sed 's/ /	/g' >  $svs'_intersectSNP.DUP.bed' &
			P1=$!
			bedtools intersect -wo -a $snp -b $svs'.INV.bed' | awk '{print $1,$2,$2+1,$10,$11,$12,$13,$14,$15}' | sed 's/ /	/g' >  $svs'_intersectSNP.INV.bed' &
			P2=$!
			bedtools intersect -wo -a $snp -b $svs'.DEL.bed' >>  $svs'_intersectSNP.DEL.tab' &
			P3=$!
			wait $P1 $P2 $P3

			echo "	Compute coverage and parse it over into a bed file"
			sambamba depth base -L  $svs'.DUP.bed' $reads | grep -v 'REF' | awk '{print $1,$2,$2+1,$3}' | sed 's/ /	/g' |  gzip - > $svs'_DUP_cov.txt.gz' &
			P1=$!
			sambamba depth base -L $svs'.INV.bed'  $reads | grep -v 'REF' |awk '{print $1,$2,$2+1,$3}' | sed 's/ /	/g' | gzip - > $svs'_INV_cov.txt.gz' &
			P2=$!
			wait $P1 $P2
		
			echo "	Intersect SV+SNP with coverage"
			bedtools intersect -wo -a $svs'_intersectSNP.DUP.bed' -b $svs'_DUP_cov.txt.gz' |  cut -f 1,2,4,5,6,8,9,13  >> $svs'_intersectSNP.DUP+cov.tab'
			bedtools intersect -wo -a $svs'_intersectSNP.INV.bed' -b $svs'_INV_cov.txt.gz' |  cut -f 1,2,4,5,6,8,9,13  >> $svs'_intersectSNP.INV+cov.tab'
		
			echo "	Cleaning up tmp files:"
			rm $svs'.DUP.bed'
			rm $svs'.INV.bed'
			rm $svs'.DEL.bed'
			rm $svs'_intersectSNP.DUP.bed'
			rm $svs'_intersectSNP.INV.bed'
			rm $svs'_DUP_cov.txt.gz'
			rm $svs'_INV_cov.txt.gz'
		done;
       fi

else
        echo "Intersects SVs with SNP for subprocessing"
        echo "File path of SVs VCF.GZ file"
	echo "File path of SNP VCF.GZ file"
	echo "File path of BAM file"
	echo "Output: Reports a tab formated SNP file splitted into DUP, DEL, INV."
fi


