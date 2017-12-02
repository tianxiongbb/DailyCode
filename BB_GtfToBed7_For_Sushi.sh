#!/bin/bash

if [ $# -lt 2 ];then
	echo0 1 $0" in.gtf out.bed7"
	exit 1
fi

awk 'BEGIN{FS=OFS="\t"} {split($9,a,"transcript_id \"");split(a[2],b,"\";");if($7=="+"){strand=1}else{strand=-1};if(b[1]){if($3=="exon"){print $1,$4,$5,b[1],255,strand,"exon"};if($3~/utr/){print $1,$4,$5,b[1],255,strand,"utr"}}}' $1 > $2
