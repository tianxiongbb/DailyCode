#!/bin/bash

if [ $# -lt 3 ];then
	echo0 1 $0" in.bed in.chrom.size out.prefix seperate_strand(0 or 1; default 0)"
	exit 1
fi

STRAND=0
if [ $# -gt 3 ];then
	STRAND=$4
fi

if [ $STRAND -eq 0 ];then
	factor=(`wc -l $1`)
	echo0 2 "calculate depth......"
	sort -k1,1 -k2,2n $1 > $3.temp.bed
	bedtools genomecov -i $3.temp.bed -g $2 -bg | awk -v factor=$factor 'BEGIN{FS=OFS="\t"} {print $1,$2,$3,$4/factor*1000000}' | sort -k1,1 -k2,2n > $3.bdg
	echo0 2 "make bigWig......"
	bedGraphToBigWig $3.bdg $2 $3.bw
	rm $3.temp.bed
else
	factor=(`wc -l $1`)
	echo0 2 "divide strands......"
	awk '$6=="+"' $1 | sort -k1,1 -k2,2n > $3.temp.watson.bed
	awk '$6=="-"' $1 | sort -k1,1 -k2,2n > $3.temp.crick.bed
	echo0 2 "calculate depth......"
	bedtools genomecov -i $3.temp.watson.bed -g $2 -bg | awk -v factor=$factor 'BEGIN{FS=OFS="\t"} {print $1,$2,$3,$4/factor*1000000}' | sort -k1,1 -k2,2n > $3.watson.bdg
	bedtools genomecov -i $3.temp.crick.bed -g $2 -bg | awk -v factor=$factor 'BEGIN{FS=OFS="\t"} {print $1,$2,$3,$4/factor*1000000}' | sort -k1,1 -k2,2n > $3.crick.bdg
	echo0 2 "make bigWig......"
	bedGraphToBigWig $3.watson.bdg $2 $3.watson.bw
	bedGraphToBigWig $3.crick.bdg $2 $3.crick.bw
	rm $3.temp.watson.bed
	rm $3.temp.crick.bed
fi
	


