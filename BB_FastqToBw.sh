#!/bin/bash

#help
help_info(){
	echo "usage:"
	echo "sh BB_FastqToBam.sh in.fastq out.bam.dir out.prefix genome"
	echo ""
	echo "Caution:"
	echo "Only used for mouse ChIP data, if other organsim, please modify the program"
	echo "The program will use these tools:"
	echo "1. Fastqc"
	echo "2. Bowtie"
	echo "3. Samtools"
}

if [ $# -lt 3 ];then
	help_info
	exit 1
fi

#
#/home/tongji2/piRNA/Software/FastQC/fastqc -O ${2} ${1}
bowtie -S -p 2 -v 2 -k 1 -q --best /data/tusers/yutianx/tongji2/Annotation/Index/${4}_bowtie/genome ${1} \
> ${2}${3}.sam 2> /data/tusers/yutianx/tongji2/piRNA/Output/Bowtie_log/${3}.log
samtools view -bhS ${2}${3}.sam > ${2}${3}.bam
samtools sort ${2}${3}.bam ${2}${3}_sort
rm ${2}${3}.sam
rm ${2}${3}.bam


