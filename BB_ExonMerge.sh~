#!/bin/bash

if [ $# -lt 1 ];then
	echo $0" in.bed out.bed"
	exit 1
fi

awk 'BEGIN{FS=OFS="\t"} {if(NR==1){chrom=$1;start=$2;ender=$3;name=$4;strand=$6}else if($4==name){ender=$3}else{print chrom,start,ender,name,0,strand;chrom=$1;start=$2;ender=$3;name=$4;strand=$6}} END{print chrom,start,ender,name,0,strand}' $1 > $2


