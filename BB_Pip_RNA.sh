#!/bin/bash

##----INTRO-----------##
# Name=BB_Pip_RNA
# Date=Apr20 ,2016
# Update=Apr 20, 2016
# Update information:

########################
# Purpose
#This file is for doing some pre-analysis for each RNA-seq samples
#1.Mapping: bowtie2 + STAR
#2.Count Feature Signal: htseq-count
#3.Plot Results: to do

#set python path

#######--Arguments--#######
help_info(){
	echo "usage:"
	echo "bash BB_Pip_RNA <option>* [-l input_left.fq] [-r input_right.fq] [-S star_index] [-B bowtie_index] [-c CPU] [-p out.prefix] [-o output_dir] [-L yes] [-M number] [-m number] [-R Y] [-H genome.gtf]"
	echo "Please use bash BB_Pip_RNA"
	echo ""
	echo "Arguments:"
	echo "-l left RNA-seq fastq file."
	echo "-r right RNA-seq fastq file. [ set to N or not set if single-end ]"
	echo "-S star index for mapping"
	echo "-B bowtie2 index for remove rRNA"
	echo "-c CPU used.---default: 1"
	echo "-p prefix output name.---default: Out"
	echo "-o output directory.---default: ./"
	echo "-L library type, default is dUTP, if yes, opossite"
	echo "-M max allowed multiple mapped reads for STAR map.---default: 100 (-1 mean not allowed)"
	echo "-m max mismatch allowed for STAR map"
	echo "-R set to N if don't remove duplicates---default: removed"
	echo "-H gtf file for htseq count, htseq will use id for feature name\n set to N if do not calculate feature signal; set to P if calculate gene and picluster signal"
	echo ""
	echo "Caution:"
	echo "The program only can deal with paired-end fr-firststrand library. If not, please modify the program"
}


if [ $# -lt 2 ];then
	help_info
	exit 1
fi

CPU=1
OUTPATH=./
GENOME=mm10
PREFIX=Out
LIBRARY=reverse
RIGHT=N
MN=100
MM=2
REMOVE=Y
HTSEQ=N

while getopts ":l:r:S:B:p:c:o:L:M:m:R:H:" Arg
do
	case $Arg in
		l)	LEFT=$OPTARG;;
		r)  RIGHT=$OPTARG;;
		S)	SI=$OPTARG;;
		B)  BI=$OPTARG;;
		p)	PREFIX=$OPTARG;;
		c)	CPU=$OPTARG;;
		o)	OUTPATH=$OPTARG;;
		L)	LIBRARY=$OPTARG;;
		M)  MN=$OPTARG;;
		m)  MM=$OPTARG;;
		R)  REMOVE=$OPTARG;;
		H)  HTSEQ=$OPTARG;;
		?)	echo "Wrong parameter!!!"
			exit 1;;
	esac
done

###########################

echo "Start!"
echo "Configure Parameters!"
echo "-l "${LEFT}
echo "-r "${RIGHT}
echo "-S "${SI}
echo "-B "${BI}
echo "-p "${PREFIX}
echo "-c "${CPU}
echo "-o "${OUTPATH}
echo "-L "${LIBRARY}
echo "-M "${MN}
echo "-m "${MM}
echo "-r "${REMOVE}
echo "-H "${HTSEQ}

###--mkdir--###
mkdir $OUTPATH
cd $OUTPATH

###--Process--###
#################
#2.Bowtie2+STAR+SAM/BAM
###RNAseq normal mapping and junction mapping, and assemble into transcript.
#-Bowtie2-#
mkdir RNAseq-rRNA
mkdir ./star
if [ ${RIGHT} == "N" ];then
#remove rRNA from whole genome RNA-seq data
# bowtie2 -rRNA
	echo "###############"
	echo "single-end mode"
	echo "###############"
	echo "###############"
	echo "bowtie2: remove rRNA"
	echo "###############"
	bowtie2 -N 1 \
		-p ${CPU} \
		-q \
		-x ${BI} \
		-U ${LEFT} \
		-S ./RNAseq-rRNA/${PREFIX}_rRNA.sam \
		--un ./RNAseq-rRNA/${PREFIX}_RNAseq-rRNA.fq
#-STAR-#
#run reference based STAR for all reads from sample.fq
#index file has been get in the path /data/tongji2/index/mm10_star_gtf 
#if use another organism, please run STAR --genomeGenerate again
#use 16 CPU, 2 mismatch, 100 max multi-map, , write all sam attributes, write unmapped reads
	echo "###############"
	echo "STAR: mapping to genome"
	echo "###############"
	STAR 	--genomeDir ${SI} \
		--runThreadN ${CPU} \
		--readFilesIn ./RNAseq-rRNA/${PREFIX}_RNAseq-rRNA.fq \
		--outFileNamePrefix ./star/${PREFIX} \
		--outFilterMismatchNmax ${MM} \
		--outFilterMultimapNmax ${MN} \
		--outSAMattributes All \
		--outFilterIntronMotifs RemoveNoncanonicalUnannotated
else
###pair-end mode
	echo "###############"
	echo "paired-end mode"
	echo "###############"
	echo "###############"
	echo "bowtie2: remove rRNA"
	echo "###############"
	bowtie2 -N 1 \
		-p ${CPU} \
		-q \
		-x ${BI} \
		-1 ${LEFT} \
		-2 ${RIGHT} \
		-S ./RNAseq-rRNA/${PREFIX}_rRNA.sam \
		--un-conc ./RNAseq-rRNA/${PREFIX}_RNAseq-rRNA.fq
	echo "###############"
	echo "STAR: mapping to genome"
	echo "###############"
	STAR 	--genomeDir ${SI} \
		--runThreadN ${CPU} \
		--readFilesIn ./RNAseq-rRNA/${PREFIX}_RNAseq-rRNA.1.fq ./RNAseq-rRNA/${PREFIX}_RNAseq-rRNA.2.fq \
		--outFileNamePrefix ./star/${PREFIX} \
		--outFilterMismatchNmax ${MM} \
		--outFilterMultimapNmax ${MN} \
		--outSAMattributes All \
		--outFilterIntronMotifs RemoveNoncanonicalUnannotated
fi

#sam to rmdup.bam
echo "###############"
echo "samtools: sam to sorted or rmdup bam"
echo "###############"
rm -rf RNAseq-rRNA
samtools view -bhS -o ./star/${PREFIX}.bam ./star/${PREFIX}Aligned.out.sam
rm -rf ./star/*.sam
samtools sort ./star/${PREFIX}.bam ./star/${PREFIX}.sort
#rm -rf ./star/${PREFIX}.bam
if [ ${REMOVE} == "Y" ];then
	samtools rmdup ./star/${PREFIX}.sort.bam ./star/${PREFIX}.sort.rmdup.bam
	rm -rf ./star/${PREFIX}.sort.bam
	BAM=./star/${PREFIX}.sort.rmdup.bam
else
	BAM=./star/${PREFIX}.sort.bam
fi

htseq-count
if [ ${HTSEQ} != "N" && ${HTSEQ} != "P" ];then
	echo "###############"
	echo "htseq: calculate signal of each feature"
	echo "###############"
	mkdir htseq
	htseq-count -m union -f bam -s ${LIBRARY} -t exon -i gene_id -q ${BAM} \
	${HTSEQ} > ./htseq/${PREFIX}.sig
	BB_NorHtseqRPM.py ./htseq/${PREFIX}.sig ./htseq/${PREFIX}.rpm
fi

if [ ${HTSEQ} == "P" ];then
	echo "###############"
	echo "htseq: calculate signal of each feature"
	echo "###############"
	mkdir htseq
	htseq-count -m union -f bam -s ${LIBRARY} -t exon -i gene_id -q ${BAM} \
	/data/tongji2/piRNA/Output/Other/picluster.gtf > ./htseq/${PREFIX}_piRNA.sig
	BB_NorHtseqRPM.py ./htseq/${PREFIX}_piRNA.sig ./htseq/${PREFIX}_piRNA.rpm
	htseq-count -m union -f bam -s ${LIBRARY} -t exon -i gene_id -q ${BAM} \
	/data/tongji2/InputForRunBT/Reference/Mus_musculus.GRCm38.82_seleno_filtered_final3.gtf > ./htseq/${PREFIX}_gene.sig
	BB_NorHtseqRPM.py ./htseq/${PREFIX}_gene.sig ./htseq/${PREFIX}_gene.rpm
fi

#mv ./star/${PREFIX}Log.final.out ./
#bam to bigwig

echo "Well done!"
echo "^_^"
echo "Finished!!!"

#########################
