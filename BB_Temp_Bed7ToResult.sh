#!/bin/bash

if [ $# -lt 2 ];then
	echo0 1 $0" in.bed7 TE.size"
	exit 1
fi

awk 'BEGIN{FS=OFS="\t"} {if(NR==FNR){a[$4]["all"]+=$7;a[$4][$5]+=$7}else{print $1,a[$1]["all"]?a[$1]["all"]:0,a[$1]["1p1"]?a[$1]["1p1"]:0,a[$1]["2p"]?a[$1]["2p"]:0,a[$1]["singleton"]?a[$1]["singleton"]:0}}' $1 $2 > ${1%.bed7}.result 


