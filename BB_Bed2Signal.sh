#!/bin/bash

if [ $# -lt 3 ];then
	echo0 1 $0" in.bed2 picluster.bed signal.prefix factor(default: unique reads / 1e6)"
	exit 1
fi

echo0 2 "intersectBed......"
intersectBed -wo -a $1 -b $2 > $3.temp
echo0 2 "calculate factor......"
if [ $# -lt 4 ];then
	FACTOR=`awk '{sum=sum+$4/$5} END{print sum}' $1`
else
	FACTOR=$4
fi
echo0 2 "calculate rpm......"
awk -v factor=$FACTOR 'BEGIN{FS=OFS="\t"} {if(NR==FNR){a[$11]=a[$11]+$4/$5}else{if(a[$4]){print $4,a[$4]*1000000/factor}else{print $4,0}}}' $3.temp $2 > $3.rpm
rm $3.temp
