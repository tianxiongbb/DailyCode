#!/bin/bash
sort -k4,4 -k1,1n $1 | awk 'BEGIN{FS=OFS="\t"} {if(NR==1){chrom=$1;start=$2;ender=$3;name=$4;strand=$6}else if($4==name){ender=$3}else{print chrom,start,ender,name,255,strand;chrom=$1;start=$2;ender=$3;name=$4;strand=$6}} END{print chrom,start,ender,name,255,strand}' > $2


