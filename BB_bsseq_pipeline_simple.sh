#!/bin/bash
cd $4

if [ $# -lt 4 ];then
	echo -e "\033[1;33;40m"$0" SRR_ID/no(already downloaded) sample_name CPU outdir"
	echo -e "only use this for paired-end BSseq data\033[0m"
	exit 1
fi

if [ "$1" != "no" ];then
	BB_download_GEO.sh $1 paired
	mv ${1}_1.fastq ${2}_1.fastq
	mv ${1}_2.fastq ${2}_2.fastq
fi

trim_galore --paired --trim1 ${2}_1.fastq ${2}_2.fastq
bismark --bowtie2 -B ${2} --genome_folder /data/tusers/yutianx/tongji2/Annotation/Index/mm10_bismark/ -1 ${2}_1_val_1.fq -2 ${2}_2_val_2.fq -p $3
bismark_methylation_extractor --buffer_size 15G -p --bedGraph --no_overlap ${2}_pe.bam

