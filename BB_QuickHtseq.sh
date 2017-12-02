#!/bin/bash

help_info(){
	echo "usage:"
	echo "BB_QuickHtseq bam gtf out.dir&prefix no/yes/reverse gene_id/gene_name"
}

if [ $# -lt 2 ];then
	help_info
	exit 1
fi

htseq-count -m intersection-strict -f bam -s ${4} -t exon -q \
-i ${5} ${1} ${2} > ${3}.sig
BB_NorHtseqRPM.py ${3}.sig ${3}.rpm 

