#1/bin/bash

if [ $# -lt 3 ];then
	echo0 1 "BB_trackUrl.sh in.bed out.txt hgsid"
	exit 0
fi

awk -v hgsid=$3 'BEGIN{FS=OFS="\t"} {leng=$3-$2;if(leng>10000){start=$2-leng/2;ender=$3+leng/2}else{start=$2-5000;ender=$3+5000};start=int(start);ender=int(ender);print "http://genome.cse.ucsc.edu/cgi-bin/hgTracks?db=mm10&lastVirtModeType=default&lastVirtModeExtraState=&virtModeType=default&virtMode=0&nonVirtPosition=&position="$1"%3A"start"-"ender"&hgsid="hgsid"&hgt.psOutput=on"}' $1 > $2


