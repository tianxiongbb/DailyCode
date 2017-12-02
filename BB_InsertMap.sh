#!/bin/bash

#This pipeline is used for small RNA-seq, based on Pip_SRNA_Mouse results
#Update: Jul 1, 2016
#Author: BB~~~Tianxiong Yu

usage () {
cat << EOF

This pipeline is used for small RNA-seq, especially for analysis of piRNA
Author: BB~~~Tianxiong Yu

usage:
	BB_Pip_SRNA.sh in.insert Index out.bed2

EOF
}

if [ $# -lt 1 ]; then
	usage
	exit 1
fi

#################
### procedure ###
#################
###map to picluster
bowtie -r -v 1 -a --best --strata -S -p 8 ${2} \
${1} > ${1}.sam
samtools view -uS -F0x4 $1.sam 2>/dev/null | \
bedtools bamtobed -i - > $1.bed && \
insertBed_to_bed2 ${1} $1.bed > ${3} && \
rm -rf $1.sam $1.bed








