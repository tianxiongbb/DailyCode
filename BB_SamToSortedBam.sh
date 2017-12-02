#!/bin/bash

if [ $# -lt 2 ];then
	echo0 1 BB_SamToSortedBam.sh" in.sam out.prefix(default:./temp)"
	exit 0
fi

samtools view -bhS $1 > $2.bam
samtools sort $2.bam $2.sort
rm $2.bam
filesize=`ls -l $2.sort | awk '{print $5}'`
if [ $filesize -gt 1024 ];then
	echo0 2 "Awesome!!!"
else
	echo0 0 "Godsh!!!"
fi


