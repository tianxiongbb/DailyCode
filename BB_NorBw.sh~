#!/bin/bash

help_info(){
	echo "usage:"
	echo "sh BB_NorBw.sh in.bw out.bw nor_factor (eg: Unique Mapped Reads / 1000000) chrom.size plus/minus"
	echo ""
}

if [ $# -lt 1 ];then
	help_info
	exit 1
fi

SN=`basename ${1}`
bigWigToBedGraph ${1} temp_${SN}.bedGraph

sort -k1,1 -k2,2n temp_${SN}.bedGraph > temp_${SN}_sorted.bedGraph

FACTOR=${3}

if [ "$5" == "minus" ];then
	awk -v fac="$FACTOR" '{FS="\t";OFS="\t"} {print $1,$2,$3,-$4/fac}' temp_${SN}_sorted.bedGraph > temp_${SN}_sorted_nor.bedGraph
else
	awk -v fac="$FACTOR" '{FS="\t";OFS="\t"} {print $1,$2,$3,$4/fac}' temp_${SN}_sorted.bedGraph > temp_${SN}_sorted_nor.bedGraph
fi


bedGraphToBigWig temp_${SN}_sorted_nor.bedGraph ${4} ${2}

rm temp_${SN}*
