#!/bin/bash

##----INTRO-----------##
# Name=piFinder
# Date=Nov15 ,2016
# Update=Jan09, 2017
# Update information:
# Remove small RNA reads mapped to protein-coding genes in antisense and in intron
# Trim 5' and 3' more stricted though it has significant boundary or it is a exon

########################
# Purpose
# This pipeline can define pachytene piRNA gene via small-RNA-seq and/or RNA-seq

#######--Arguments--#######
help_info(){
	echo ""
	echo -e "\033[32m  =========================================================================================================================\033[0m"
	echo "   Author: Tianxiong & Kaili"
	echo ""
	echo "   Usage:"
	echo "	bash BB_smallRNA_map.sh <option>* [-q srna.fq or srna.insert] [-g genome] [-G reference.gtf]"
	echo ""
	echo "   Optional arguments:"
	echo "	-p prefix name for output files. --default: ./result"
	echo "	-c CPU number used for bowtie mapping. --default:1"
	echo ""
	echo -e "\033[31m  !!! Caution: \033[0m"
	echo -e "	The program only can deal with single-end ff-firststrand library. If not, please modify the program."
	echo "	Program needed: bowtie fastq_to_insert insertBed_to_bed2 samtools bedtools"
	echo "	Make sure all the program is contained by environment viariable PATH."
	echo ""
	echo "    (￣(工)￣) Enjoy yourself~~~"
	echo -e "\033[32m  =========================================================================================================================\033[0m"
	echo ""
}
##TODO(BB): all index file needed, add to the help info. 

if [ $# -lt 3 ];then
	help_info
	exit 1
fi

#############################
# ARGS reading and checking #
#############################
OUTPUT=./result.bed2
CPU=1
INDEXPATH=/data/tusers/yutianx/tongji2/Annotation/Index
CHROMPATH=/data/tusers/yutianx/tongji2/Annotation/ChromSize
GENOMEPATH=/data/tusers/yutianx/tongji2/Annotation/Fasta

while getopts "hvq:c:p:g:G:" OPTION; do
	case $OPTION in
		h)	help_info && exit 0 ;;
		q)	INPUT_FASTQ=${OPTARG} ;;
		p)	PREFIX=${OPTARG} ;;
		c)	CPU=${OPTARG} ;;
		v)	echo "BB_smallRNA_map VERSION: Beta 1.0" && exit 0 ;;
		g)	GENOME=${OPTARG};;
		G)	GTF=${OPTARG};;
		*)	usage && exit 1 ;;
	esac
done

ANNO_PATH=/data/tusers/yutianx/tongji2/Software/piPipes/common/${GENOME}



INDEXPATH=`dirname $INDEXPATH`"/"`basename $INDEXPATH`
CHROMPATH=`dirname $CHROMPATH`"/"`basename $CHROMPATH`
GENOMEPATH=`dirname $GENOMEPATH`"/"`basename $GENOMEPATH`

###########
# process #
###########
#check small RNA-seq data for mapping
if [ ! -f ${INPUT_FASTQ} ];then
	echo -e "\033[40;31m\033[1mthere is no file in "${INPUT_FASTQ}". Exit\033[0m"
	exit 1
fi

#fastq to insert for space saving
DATE=`date --date="-24 hour"`
echo -e "\033[31mcreate insert format file for saving space\t"$DATE"\033[0m"
SUFFIX=${INPUT_FASTQ##*.}
if [ "$SUFFIX" = "fastq" -o "$SUFFIX" = "fq" ];then
	fastq_to_insert ${INPUT_FASTQ} ${PREFIX}.insert
else
	cp ${INPUT_FASTQ} ${PREFIX}.insert
fi

#map to rRNA
DATE=`date --date="-24 hour"`
echo -e "\033[31mmap to rRNA\t"$DATE"\033[0m"
bowtie -r -v 1 -a --best --strata -S -p ${CPU} ${ANNO_PATH}/BowtieIndex/rRNA ${PREFIX}.insert --un ${PREFIX}_rRNA.insert --al /dev/null > ${PREFIX}.rRNA.sam
samtools view -uS -F0x4 ${PREFIX}.rRNA.sam 2>/dev/null | bedtools bamtobed -i - > ${PREFIX}.rRNA.bed && insertBed_to_bed2 ${PREFIX}.insert ${PREFIX}.rRNA.bed > ${PREFIX}.rRNA.bed2
rm ${PREFIX}.rRNA.bed ${PREFIX}.rRNA.sam 
rm ${PREFIX}.insert
awk '{a[$7]=$4} END{for(i in a){l[length(i)]+=a[i]};for(i=1;i<=50;i++){print i"\t"(l[i]?l[i]:0)}}' ${PREFIX}.rRNA.bed2 > ${PREFIX}.rRNA.lendis

#map to miRNA hairpin
DATE=`date --date="-24 hour"`
echo -e "\033[33mmap to hairpin\t"$DATE"\033[0m"
bowtie -r -v 1 -m 1 -S --best --strata -p ${CPU} ${ANNO_PATH}/BowtieIndex/hairpin ${PREFIX}_rRNA.insert --un ${PREFIX}_rRNA_miRNA.insert --al /dev/null > ${PREFIX}.hairpin.sam
samtools view -uS -F0x4 ${PREFIX}.hairpin.sam 2>/dev/null | bedtools bamtobed -i - > ${PREFIX}.hairpin.bed && insertBed_to_bed2 ${PREFIX}_rRNA.insert ${PREFIX}.hairpin.bed > ${PREFIX}.hairpin.bed2
rm ${PREFIX}.hairpin.bed ${PREFIX}.hairpin.sam 
rm ${PREFIX}_rRNA.insert
HAIRPIN_READS=`awk '{sum+=$4/$5} END{print sum}' ${PREFIX}.hairpin.bed2`
awk '{a[$7]=$4} END{for(i in a){l[length(i)]+=a[i]};for(i=1;i<=50;i++){print i"\t"(l[i]?l[i]:0)}}' ${PREFIX}.hairpin.bed2 > ${PREFIX}.hairpin.lendis

#map to other ncRNA like snRNA,snoRNA,tRNA,processed_transcript......
#DATE=`date --date="-24 hour"`
#echo -e "\033[32mmap to ncRNA\t"$DATE"\033[0m"
#bowtie -r -v 1 -a -S --best --strata -p ${CPU} ${INDEXPATH}/${GENOME}_bowtie/ncRNA ${PREFIX}_rRNA_miRNA.insert --un ${PREFIX}_rRNA_miRNA_ncRNA.insert --al /dev/null > ${PREFIX}.ncRNA.sam
#samtools view -uS -F0x4 ${PREFIX}.ncRNA.sam 2>/dev/null | bedtools bamtobed -i - > ${PREFIX}.ncRNA.bed && insertBed_to_bed2 ${PREFIX}_rRNA_miRNA.insert ${PREFIX}.ncRNA.bed > ${PREFIX}.ncRNA.bed2
#rm ${PREFIX}.ncRNA.bed ${PREFIX}.ncRNA.sam 
#rm ${PREFIX}_rRNA_miRNA.insert
#awk '{a[$7]=$4} END{for(i in a){l[length(i)]+=a[i]};for(i=1;i<=50;i++){print i"\t"(l[i]?l[i]:0)}}' ${PREFIX}.ncRNA.bed2 > ${PREFIX}.ncRNA.lendis

#length filtering to 24-32bp
DATE=`date --date="-24 hour"`
echo -e "\033[31mfilter length\t"$DATE"\033[0m"
awk '{FS=OFS="\t"} {if(length($1)>23 && length($1)<33){print $0}}' ${PREFIX}_rRNA_miRNA.insert > ${PREFIX}.pilikeRNA.insert
rm ${PREFIX}_rRNA_miRNA.insert

#map to genome
DATE=`date --date="-24 hour"`
echo -e "\033[33mmap to genome\t"$DATE"\033[0m"
bowtie -v 1 -r -a -S --best --strata -p ${CPU} ${INDEXPATH}/${GENOME}_bowtie/genome ${PREFIX}.pilikeRNA.insert > ${PREFIX}.piRNA.sam

samtools view -uS -F0x4 ${PREFIX}.piRNA.sam 2>/dev/null | bedtools bamtobed -i - > ${PREFIX}.piRNA.insert.bed && insertBed_to_bed2 ${PREFIX}.pilikeRNA.insert ${PREFIX}.piRNA.insert.bed > ${PREFIX}.piRNA.bed2
rm -rf ${PREFIX}.piRNA.sam ${PREFIX}.piRNA.insert.bed ${PREFIX}.piRNA.insert
echo "finished!"

