#!/bin/bash

if [ $# -lt 2 ];then
	echo0 1 $0" in.bed Chain1(aToB) Chain2(bToA) out.bed minMatch"
	exit 1
fi

liftOver $1 $2 ${4}_atob /dev/null -minMatch=$5 -multiple
awk 'BEGIN{FS=OFS="\t"} {print $1,$2,$3,$4,$3-$2+1,$6}' ${4}_atob > t && sort -k4,4 -k5,5nr t > ${4}_atob
rm t
awk '!a[$4]++' ${4}_atob > ${4}_atob_only
liftOver ${4}_atob_only $3 ${4}_btoa /dev/null -minMatch=$5 -multiple
awk 'BEGIN{FS=OFS="\t"} {print $1,$2,$3,$4,$3-$2+1,$6}' ${4}_btoa > t && sort -k4,4 -k5,5nr t > ${4}_btoa
rm t
awk '!a[$4]++' ${4}_btoa > ${4}_btoa_only
intersectBed -wo -s -a $1 -b ${4}_btoa_only > ${4}_intersect
awk '$4==$10' ${4}_intersect | cut -f 4 > ${4}_pass_list
awk 'BEGIN{FS=OFS="\t"} {if(NR==FNR){a[$1]=1}else{if(a[$4]){print $0}}}' ${4}_pass_list ${4}_atob_only > $4.reciprocal.bed
mv ${4}_atob_only ${4}.bed
rm ${4}_*

