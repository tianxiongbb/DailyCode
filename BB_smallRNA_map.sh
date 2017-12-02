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
	echo "	bash BB_smallRNA_map.sh <option>* [-q srna.fq] [-g genome] [-G reference.gtf]"
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
fastq_to_insert ${INPUT_FASTQ} ${PREFIX}.insert

#map to rRNA
DATE=`date --date="-24 hour"`
echo -e "\033[31mmap to rRNA\t"$DATE"\033[0m"
bowtie -r -v 1 -a --best --strata -S -p ${CPU} ${INDEXPATH}/${GENOME}_bowtie/rRNA ${PREFIX}.insert --un ${PREFIX}_rRNA.insert --al /dev/null > ${PREFIX}.rRNA.sam
samtools view -uS -F0x4 ${PREFIX}.rRNA.sam 2>/dev/null | bedtools bamtobed -i - > ${PREFIX}.rRNA.bed && insertBed_to_bed2 ${PREFIX}.insert ${PREFIX}.rRNA.bed > ${PREFIX}.rRNA.bed2
rm ${PREFIX}.rRNA.bed ${PREFIX}.rRNA.sam 
rm ${PREFIX}.insert

#map to miRNA hairpin
DATE=`date --date="-24 hour"`
echo -e "\033[33mmap to hairpin\t"$DATE"\033[0m"
bowtie -r -v 1 -m 1 -S --best --strata -p ${CPU} ${INDEXPATH}/${GENOME}_bowtie/hairpin ${PREFIX}_rRNA.insert --un ${PREFIX}_rRNA_miRNA.insert --al /dev/null > ${PREFIX}.hairpin.sam
samtools view -uS -F0x4 ${PREFIX}.hairpin.sam 2>/dev/null | bedtools bamtobed -i - > ${PREFIX}.hairpin.bed && insertBed_to_bed2 ${PREFIX}_rRNA.insert ${PREFIX}.hairpin.bed > ${PREFIX}.hairpin.bed2
rm ${PREFIX}.hairpin.bed ${PREFIX}.hairpin.sam 
rm ${PREFIX}_rRNA.insert

#map to other ncRNA like snRNA,snoRNA,tRNA,processed_transcript......
DATE=`date --date="-24 hour"`
echo -e "\033[32mmap to ncRNA\t"$DATE"\033[0m"
bowtie -r -v 1 -a -S --best --strata -p ${CPU} ${INDEXPATH}/${GENOME}_bowtie/ncRNA ${PREFIX}_rRNA_miRNA.insert --un ${PREFIX}_rRNA_miRNA_ncRNA.insert --al /dev/null > ${PREFIX}.ncRNA.sam
samtools view -uS -F0x4 ${PREFIX}.ncRNA.sam 2>/dev/null | bedtools bamtobed -i - > ${PREFIX}.ncRNA.bed && insertBed_to_bed2 ${PREFIX}_rRNA_miRNA.insert ${PREFIX}.ncRNA.bed > ${PREFIX}.ncRNA.bed2
rm ${PREFIX}.ncRNA.bed ${PREFIX}.ncRNA.sam 
rm ${PREFIX}_rRNA_miRNA.insert

#length filtering to 24-32bp
DATE=`date --date="-24 hour"`
echo -e "\033[31mfilter length\t"$DATE"\033[0m"
awk '{FS=OFS="\t"} {if(length($1)>23 && length($1)<33){print $0}}' ${PREFIX}_rRNA_miRNA_ncRNA.insert > ${PREFIX}.pilikeRNA.insert
rm ${PREFIX}_rRNA_miRNA_ncRNA.insert

#map to genome
DATE=`date --date="-24 hour"`
echo -e "\033[33mmap to genome\t"$DATE"\033[0m"
bowtie -v 1 -r -a -S --best --strata -p ${CPU} ${INDEXPATH}/${GENOME}_bowtie/genome ${PREFIX}.pilikeRNA.insert > ${PREFIX}.pilikeRNA.sam

samtools view -uS -F0x4 ${PREFIX}.pilikeRNA.sam 2>/dev/null | bedtools bamtobed -i - > ${PREFIX}.pilikeRNA.insert.bed && insertBed_to_bed2 ${PREFIX}.pilikeRNA.insert ${PREFIX}.pilikeRNA.insert.bed > ${PREFIX}.pilikeRNA.bed2
rm -rf ${PREFIX}.pilikeRNA.sam ${PREFIX}.pilikeRNA.insert.bed ${PREFIX}.pilikeRNA.insert

###remove ncRNA reads from ensembl gtf file
DATE=`date --date="-24 hour"`
echo -e "\033[31mremove ncRNA reads\t"$DATE"\033[0m"
#snoRNA scaRNA snRNA miRNA rRNA Mt_tRNA Mt_rRNA
sort -k1,1 -k2,2n ${PREFIX}.pilikeRNA.bed2 > t_${PREFIX} && mv t_${PREFIX} ${PREFIX}.pilikeRNA.bed2 
grep "gene_biotype \"snoRNA\"" ${GTF} > t${PREFIX}.gtf && NL=(`wc -l t${PREFIX}.gtf`)
if [ $NL -gt 0 ];then
	BB_GtfToExonForEachGene.py t${PREFIX}.gtf ${PREFIX}.snoRNA.bed name && cat ${PREFIX}.snoRNA.bed >> temp_${PREFIX}_ncRNA.bed && rm ${PREFIX}.snoRNA.bed 
fi
grep "gene_biotype \"scaRNA\"" ${GTF} > t${PREFIX}.gtf && NL=(`wc -l t${PREFIX}.gtf`)
if [ $NL -gt 0 ];then
	BB_GtfToExonForEachGene.py t${PREFIX}.gtf ${PREFIX}.scaRNA.bed name && cat ${PREFIX}.scaRNA.bed >> temp_${PREFIX}_ncRNA.bed && rm ${PREFIX}.scaRNA.bed 
fi
grep "gene_biotype \"snRNA\"" ${GTF} > t${PREFIX}.gtf && NL=(`wc -l t${PREFIX}.gtf`)
if [ $NL -gt 0 ];then
	BB_GtfToExonForEachGene.py t${PREFIX}.gtf ${PREFIX}.snRNA.bed name && cat ${PREFIX}.snRNA.bed >> temp_${PREFIX}_ncRNA.bed && rm ${PREFIX}.snRNA.bed 
fi
grep "gene_biotype \"miRNA\"" ${GTF} > t${PREFIX}.gtf && NL=(`wc -l t${PREFIX}.gtf`)
if [ $NL -gt 0 ];then
	BB_GtfToExonForEachGene.py t${PREFIX}.gtf ${PREFIX}.miRNA.bed name && cat ${PREFIX}.miRNA.bed >> temp_${PREFIX}_ncRNA.bed && rm ${PREFIX}.miRNA.bed 
fi
grep "gene_biotype \"rRNA\"" ${GTF} > t${PREFIX}.gtf && NL=(`wc -l t${PREFIX}.gtf`)
if [ $NL -gt 0 ];then
	BB_GtfToExonForEachGene.py t${PREFIX}.gtf ${PREFIX}.rRNA.bed name && cat ${PREFIX}.rRNA.bed >> temp_${PREFIX}_ncRNA.bed && rm ${PREFIX}.rRNA.bed 
fi
grep "gene_biotype \"Mt_tRNA\"" ${GTF} > t${PREFIX}.gtf && NL=(`wc -l t${PREFIX}.gtf`)
if [ $NL -gt 0 ];then
	BB_GtfToExonForEachGene.py t${PREFIX}.gtf ${PREFIX}.Mt_tRNA.bed name && cat ${PREFIX}.Mt_tRNA.bed >> temp_${PREFIX}_ncRNA.bed && rm ${PREFIX}.Mt_tRNA.bed 
fi
grep "gene_biotype \"Mt_rRNA\"" ${GTF} > t${PREFIX}.gtf && NL=(`wc -l t${PREFIX}.gtf`)
if [ $NL -gt 0 ];then
	BB_GtfToExonForEachGene.py t${PREFIX}.gtf ${PREFIX}.Mt_rRNA.bed name && cat ${PREFIX}.Mt_rRNA.bed >> temp_${PREFIX}_ncRNA.bed && rm ${PREFIX}.Mt_rRNA.bed 
fi
sort -k1,1 -k2,2n temp_ncRNA.bed > t${PREFIX} && mv t${PREFIX} temp_ncRNA.bed
bedtools intersect -v -sorted -s -f 0.5 -wa -a ${PREFIX}.pilikeRNA.bed2 -b temp_${PREFIX}_ncRNA.bed | sort -k1,1 -k2,2n > ${PREFIX}.piRNA.bed2
rm temp_${PREFIX}_ncRNA.bed t${PREFIX}.gtf

