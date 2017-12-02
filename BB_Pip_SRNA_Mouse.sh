#!/bin/bash

#This pipeline is used for small RNA-seq, especially for analysis of piRNA
#Update: Apr 16, 2016
#Author: BB~~~Tianxiong Yu

#########
# USAGE #
#########
usage () {
cat << EOF

This pipeline is used for small RNA-seq, especially for analysis of piRNA
Author: BB~~~Tianxiong Yu

usage:
	BB_Pip_SRNA.sh
		-i Input small RNA-seq file (.fq or .fq.gz)
		-g Genome assemble used (eg. mm10, hg38)
		-c CPU used for the pipeline (default: 1)
		-o Output directory (default: ./)
		-p Prefix name for output

options:
	-h Help Information
	-v Version Information
EOF
}

#############################
# ARGS reading and checking #
#############################
CPU=1

while getopts "hvi:c:o:g:p:" OPTION; do
	case $OPTION in
		h)	usage && exit 0 ;;
		i)	INPUT_FASTQ=${OPTARG} ;;
		o)	OUTDIR=${OPTARG} ;;
		c)	CPU=${OPTARG} ;;
		v)	echo "BB_Pip_SRNA VERSION: Beta 1.0" && exit 0 ;;
		g)	GENOME=${OPTARG};;
		p)	PREFIX=${OPTARG};;
		*)	usage && exit 1 ;;
	esac
done

if [ $# -lt 1 ]; then
	usage
	exit 1
fi

[[ -z "${INPUT_FASTQ}" ]] && usage && echo "Missing option -i for input fastq file or file does not exist" "error"
[[ -z ${GENOME} ]]  && usage && echo "Missing option -g for specifying which genome assembly to use" "error"
[ ! -f "${INPUT_FASTQ}" ] && echo "Cannot find input file $INPUT_FASTQ" "error"

##################
## Path setting ##
##################
piPipe_Path=/data/tongji2/piPipes/bin/

#################
### procedure ###
#################
mkdir ${OUTDIR}
cd ${OUTDIR}

###mapping
# ${piPipe_Path}piPipes_fastq_to_insert ${INPUT_FASTQ} ${PREFIX}.insert

# bowtie2 -r -N 0 -a -p ${CPU} -x /data/tongji2/InputForRunBT/Index/${GENOME}_bowtie2/ \
# -U ${PREFIX}.insert -S ${PREFIX}.sam 2>/dev/null

# samtools view -uS -F0x4 ${PREFIX}.sam 2>/dev/null | \
# bedtools bamtobed -i - > ${PREFIX}.insert.bed && \
# ${piPipe_Path}piPipes_insertBed_to_bed2 ${PREFIX}.insert ${PREFIX}.insert.bed > ${PREFIX}.bed2 && \
# rm -rf ${PREFIX}.sam ${PREFIX}.insert.bed

# ###seperate uniq map and multiple map
# awk '{FS="\t";OFS="\t"} {if($5>1){print $0}}' ${PREFIX}.bed2 > ${PREFIX}_multi.bed2
# awk '{FS="\t";OFS="\t"} {if($5<2){print $0}}' ${PREFIX}.bed2 > ${PREFIX}_uniq.bed2

# ###get length distribution and strand distribution
# python /data/tongji2/piRNA/Code/Tianxiong/Apr_2016/Apr12_GetMultipleHumanLength.py \
# ${PREFIX}_multi.bed2 ${PREFIX}_multi_length.txt ${PREFIX}_multi_strand.txt

# python /data/tongji2/piRNA/Code/Tianxiong/Apr_2016/Apr12_GetMultipleHumanLength.py \
# ${PREFIX}_uniq.bed2 ${PREFIX}_uniq_length.txt ${PREFIX}_uniq_strand.txt

# python /data/tongji2/piRNA/Code/Tianxiong/Apr_2016/Apr12_GetMultipleHumanLength.py \
# ${PREFIX}.bed2 ${PREFIX}_length.txt ${PREFIX}_strand.txt

###map summary
TOTAL_READS=`awk 'BEGIN{sum=0} {sum=sum+$2} END{print sum}' ${PREFIX}.insert`
TOTAL_MAP_READS=`awk 'BEGIN{sum=0} {sum=sum+$4/$5} END{print sum}' ${PREFIX}.bed2`
UNIQ_MAP_READS=`awk 'BEGIN{sum=0} {sum=sum+$4/$5} END{print sum}' ${PREFIX}_uniq.bed2`
MULTI_MAP_READS=`awk 'BEGIN{sum=0} {sum=sum+$4/$5} END{print sum}' ${PREFIX}_multi.bed2`

###intersect with each feature
mkdir Intersect_Result
for RM in LINE SINE LTR DNA Satellite Simple_repeat
do
	intersectBed -s -F 0.5 -wo -a /data/tongji2/InputForRunBT/RepeatMask/${GENOME}_${RM}.bed -b \
	${PREFIX}.bed2 > ./Intersect_Result/All_${RM}.txt
done

for RM in LINE SINE LTR DNA Satellite Simple_repeat
do
	intersectBed -s -F 0.5 -wo -a /data/tongji2/InputForRunBT/RepeatMask/${GENOME}_${RM}.bed -b \
	${PREFIX}_uniq.bed2 > ./Intersect_Result/Uniq_${RM}.txt
done

for RM in LINE SINE LTR DNA Satellite Simple_repeat
do
	intersectBed -s -F 0.5 -wo -a /data/tongji2/InputForRunBT/RepeatMask/${GENOME}_${RM}.bed -b \
	${PREFIX}_multi.bed2 > ./Intersect_Result/Multi_${RM}.txt
done

for RM in LINE SINE LTR DNA Satellite Simple_repeat
do
	intersectBed -S -F 0.5 -wo -a /data/tongji2/InputForRunBT/RepeatMask/${GENOME}_${RM}.bed -b \
	${PREFIX}.bed2 > ./Intersect_Result/All_${RM}_RVS.txt
done

for RM in LINE SINE LTR DNA Satellite Simple_repeat
do
	intersectBed -S -F 0.5 -wo -a /data/tongji2/InputForRunBT/RepeatMask/${GENOME}_${RM}.bed -b \
	${PREFIX}_uniq.bed2 > ./Intersect_Result/Uniq_${RM}_RVS.txt
done

for RM in LINE SINE LTR DNA Satellite Simple_repeat
do
	intersectBed -S -F 0.5 -wo -a /data/tongji2/InputForRunBT/RepeatMask/${GENOME}_${RM}.bed -b \
	${PREFIX}_multi.bed2 > ./Intersect_Result/Multi_${RM}_RVS.txt
done

ALL_LINE_PACHY_RVS_READS=`intersectBed -S -f 0.5 -wa -a ${PREFIX}.bed2 -b /data/tongji2/InputForRunBT/RepeatMask/${GENOME}_LINE.bed | intersectBed -s -f 0.5 -wa -a - -b /data/tongji2/piRNA/Output/PiRNA_Loci/LiXin_Pachytene.bed | awk 'BEGIN{sum=0} {sum=sum+($4/$5)} END{print sum}' -`
ALL_SINE_PACHY_RVS_READS=`intersectBed -S -f 0.5 -wa -a ${PREFIX}.bed2 -b /data/tongji2/InputForRunBT/RepeatMask/${GENOME}_SINE.bed | intersectBed -s -f 0.5 -wa -a - -b /data/tongji2/piRNA/Output/PiRNA_Loci/LiXin_Pachytene.bed | awk 'BEGIN{sum=0} {sum=sum+($4/$5)} END{print sum}' -`
ALL_LTR_PACHY_RVS_READS=`intersectBed -S -f 0.5 -wa -a ${PREFIX}.bed2 -b /data/tongji2/InputForRunBT/RepeatMask/${GENOME}_LTR.bed | intersectBed -s -f 0.5 -wa -a - -b /data/tongji2/piRNA/Output/PiRNA_Loci/LiXin_Pachytene.bed | awk 'BEGIN{sum=0} {sum=sum+($4/$5)} END{print sum}' -`
ALL_LINE_PACHY_READS=`intersectBed -s -f 0.5 -wa -a ${PREFIX}.bed2 -b /data/tongji2/InputForRunBT/RepeatMask/${GENOME}_LINE.bed | intersectBed -s -f 0.5 -wa -a - -b /data/tongji2/piRNA/Output/PiRNA_Loci/LiXin_Pachytene.bed | awk 'BEGIN{sum=0} {sum=sum+($4/$5)} END{print sum}' -`
ALL_SINE_PACHY_READS=`intersectBed -s -f 0.5 -wa -a ${PREFIX}.bed2 -b /data/tongji2/InputForRunBT/RepeatMask/${GENOME}_SINE.bed | intersectBed -s -f 0.5 -wa -a - -b /data/tongji2/piRNA/Output/PiRNA_Loci/LiXin_Pachytene.bed | awk 'BEGIN{sum=0} {sum=sum+($4/$5)} END{print sum}' -`
ALL_LTR_PACHY_READS=`intersectBed -s -f 0.5 -wa -a ${PREFIX}.bed2 -b /data/tongji2/InputForRunBT/RepeatMask/${GENOME}_LTR.bed | intersectBed -s -f 0.5 -wa -a - -b /data/tongji2/piRNA/Output/PiRNA_Loci/LiXin_Pachytene.bed | awk 'BEGIN{sum=0} {sum=sum+($4/$5)} END{print sum}' -`
ALL_LINE_PACHY_RVS_READS_RVS=`intersectBed -S -f 0.5 -wa -a ${PREFIX}.bed2 -b /data/tongji2/InputForRunBT/RepeatMask/${GENOME}_LINE.bed | intersectBed -S -f 0.5 -wa -a - -b /data/tongji2/piRNA/Output/PiRNA_Loci/LiXin_Pachytene.bed | awk 'BEGIN{sum=0} {sum=sum+($4/$5)} END{print sum}' -`
ALL_SINE_PACHY_RVS_READS_RVS=`intersectBed -S -f 0.5 -wa -a ${PREFIX}.bed2 -b /data/tongji2/InputForRunBT/RepeatMask/${GENOME}_SINE.bed | intersectBed -S -f 0.5 -wa -a - -b /data/tongji2/piRNA/Output/PiRNA_Loci/LiXin_Pachytene.bed | awk 'BEGIN{sum=0} {sum=sum+($4/$5)} END{print sum}' -`
ALL_LTR_PACHY_RVS_READS_RVS=`intersectBed -S -f 0.5 -wa -a ${PREFIX}.bed2 -b /data/tongji2/InputForRunBT/RepeatMask/${GENOME}_LTR.bed | intersectBed -S -f 0.5 -wa -a - -b /data/tongji2/piRNA/Output/PiRNA_Loci/LiXin_Pachytene.bed | awk 'BEGIN{sum=0} {sum=sum+($4/$5)} END{print sum}' -`
ALL_LINE_PACHY_READS_RVS=`intersectBed -s -f 0.5 -wa -a ${PREFIX}.bed2 -b /data/tongji2/InputForRunBT/RepeatMask/${GENOME}_LINE.bed | intersectBed -S -f 0.5 -wa -a - -b /data/tongji2/piRNA/Output/PiRNA_Loci/LiXin_Pachytene.bed | awk 'BEGIN{sum=0} {sum=sum+($4/$5)} END{print sum}' -`
ALL_SINE_PACHY_READS_RVS=`intersectBed -s -f 0.5 -wa -a ${PREFIX}.bed2 -b /data/tongji2/InputForRunBT/RepeatMask/${GENOME}_SINE.bed | intersectBed -S -f 0.5 -wa -a - -b /data/tongji2/piRNA/Output/PiRNA_Loci/LiXin_Pachytene.bed | awk 'BEGIN{sum=0} {sum=sum+($4/$5)} END{print sum}' -`
ALL_LTR_PACHY_READS_RVS=`intersectBed -s -f 0.5 -wa -a ${PREFIX}.bed2 -b /data/tongji2/InputForRunBT/RepeatMask/${GENOME}_LTR.bed | intersectBed -S -f 0.5 -wa -a - -b /data/tongji2/piRNA/Output/PiRNA_Loci/LiXin_Pachytene.bed | awk 'BEGIN{sum=0} {sum=sum+($4/$5)} END{print sum}' -`
UNIQ_LINE_PACHY_RVS_READS=`intersectBed -S -f 0.5 -wa -a ${PREFIX}_uniq.bed2 -b /data/tongji2/InputForRunBT/RepeatMask/${GENOME}_LINE.bed | intersectBed -s -f 0.5 -wa -a - -b /data/tongji2/piRNA/Output/PiRNA_Loci/LiXin_Pachytene.bed | awk 'BEGIN{sum=0} {sum=sum+($4/$5)} END{print sum}' -`
UNIQ_SINE_PACHY_RVS_READS=`intersectBed -S -f 0.5 -wa -a ${PREFIX}_uniq.bed2 -b /data/tongji2/InputForRunBT/RepeatMask/${GENOME}_SINE.bed | intersectBed -s -f 0.5 -wa -a - -b /data/tongji2/piRNA/Output/PiRNA_Loci/LiXin_Pachytene.bed | awk 'BEGIN{sum=0} {sum=sum+($4/$5)} END{print sum}' -`
UNIQ_LTR_PACHY_RVS_READS=`intersectBed -S -f 0.5 -wa -a ${PREFIX}_uniq.bed2 -b /data/tongji2/InputForRunBT/RepeatMask/${GENOME}_LTR.bed | intersectBed -s -f 0.5 -wa -a - -b /data/tongji2/piRNA/Output/PiRNA_Loci/LiXin_Pachytene.bed | awk 'BEGIN{sum=0} {sum=sum+($4/$5)} END{print sum}' -`
UNIQ_LINE_PACHY_READS=`intersectBed -s -f 0.5 -wa -a ${PREFIX}_uniq.bed2 -b /data/tongji2/InputForRunBT/RepeatMask/${GENOME}_LINE.bed | intersectBed -s -f 0.5 -wa -a - -b /data/tongji2/piRNA/Output/PiRNA_Loci/LiXin_Pachytene.bed | awk 'BEGIN{sum=0} {sum=sum+($4/$5)} END{print sum}' -`
UNIQ_SINE_PACHY_READS=`intersectBed -s -f 0.5 -wa -a ${PREFIX}_uniq.bed2 -b /data/tongji2/InputForRunBT/RepeatMask/${GENOME}_SINE.bed | intersectBed -s -f 0.5 -wa -a - -b /data/tongji2/piRNA/Output/PiRNA_Loci/LiXin_Pachytene.bed | awk 'BEGIN{sum=0} {sum=sum+($4/$5)} END{print sum}' -`
UNIQ_LTR_PACHY_READS=`intersectBed -s -f 0.5 -wa -a ${PREFIX}_uniq.bed2 -b /data/tongji2/InputForRunBT/RepeatMask/${GENOME}_LTR.bed | intersectBed -s -f 0.5 -wa -a - -b /data/tongji2/piRNA/Output/PiRNA_Loci/LiXin_Pachytene.bed | awk 'BEGIN{sum=0} {sum=sum+($4/$5)} END{print sum}' -`
UNIQ_LINE_PACHY_RVS_READS_RVS=`intersectBed -S -f 0.5 -wa -a ${PREFIX}_uniq.bed2 -b /data/tongji2/InputForRunBT/RepeatMask/${GENOME}_LINE.bed | intersectBed -S -f 0.5 -wa -a - -b /data/tongji2/piRNA/Output/PiRNA_Loci/LiXin_Pachytene.bed | awk 'BEGIN{sum=0} {sum=sum+($4/$5)} END{print sum}' -`
UNIQ_SINE_PACHY_RVS_READS_RVS=`intersectBed -S -f 0.5 -wa -a ${PREFIX}_uniq.bed2 -b /data/tongji2/InputForRunBT/RepeatMask/${GENOME}_SINE.bed | intersectBed -S -f 0.5 -wa -a - -b /data/tongji2/piRNA/Output/PiRNA_Loci/LiXin_Pachytene.bed | awk 'BEGIN{sum=0} {sum=sum+($4/$5)} END{print sum}' -`
UNIQ_LTR_PACHY_RVS_READS_RVS=`intersectBed -S -f 0.5 -wa -a ${PREFIX}_uniq.bed2 -b /data/tongji2/InputForRunBT/RepeatMask/${GENOME}_LTR.bed | intersectBed -S -f 0.5 -wa -a - -b /data/tongji2/piRNA/Output/PiRNA_Loci/LiXin_Pachytene.bed | awk 'BEGIN{sum=0} {sum=sum+($4/$5)} END{print sum}' -`
UNIQ_LINE_PACHY_READS_RVS=`intersectBed -s -f 0.5 -wa -a ${PREFIX}_uniq.bed2 -b /data/tongji2/InputForRunBT/RepeatMask/${GENOME}_LINE.bed | intersectBed -S -f 0.5 -wa -a - -b /data/tongji2/piRNA/Output/PiRNA_Loci/LiXin_Pachytene.bed | awk 'BEGIN{sum=0} {sum=sum+($4/$5)} END{print sum}' -`
UNIQ_SINE_PACHY_READS_RVS=`intersectBed -s -f 0.5 -wa -a ${PREFIX}_uniq.bed2 -b /data/tongji2/InputForRunBT/RepeatMask/${GENOME}_SINE.bed | intersectBed -S -f 0.5 -wa -a - -b /data/tongji2/piRNA/Output/PiRNA_Loci/LiXin_Pachytene.bed | awk 'BEGIN{sum=0} {sum=sum+($4/$5)} END{print sum}' -`
UNIQ_LTR_PACHY_READS_RVS=`intersectBed -s -f 0.5 -wa -a ${PREFIX}_uniq.bed2 -b /data/tongji2/InputForRunBT/RepeatMask/${GENOME}_LTR.bed | intersectBed -S -f 0.5 -wa -a - -b /data/tongji2/piRNA/Output/PiRNA_Loci/LiXin_Pachytene.bed | awk 'BEGIN{sum=0} {sum=sum+($4/$5)} END{print sum}' -`
MULTI_LINE_PACHY_RVS_READS=`intersectBed -S -f 0.5 -wa -a ${PREFIX}_multi.bed2 -b /data/tongji2/InputForRunBT/RepeatMask/${GENOME}_LINE.bed | intersectBed -s -f 0.5 -wa -a - -b /data/tongji2/piRNA/Output/PiRNA_Loci/LiXin_Pachytene.bed | awk 'BEGIN{sum=0} {sum=sum+($4/$5)} END{print sum}' -`
MULTI_SINE_PACHY_RVS_READS=`intersectBed -S -f 0.5 -wa -a ${PREFIX}_multi.bed2 -b /data/tongji2/InputForRunBT/RepeatMask/${GENOME}_SINE.bed | intersectBed -s -f 0.5 -wa -a - -b /data/tongji2/piRNA/Output/PiRNA_Loci/LiXin_Pachytene.bed | awk 'BEGIN{sum=0} {sum=sum+($4/$5)} END{print sum}' -`
MULTI_LTR_PACHY_RVS_READS=`intersectBed -S -f 0.5 -wa -a ${PREFIX}_multi.bed2 -b /data/tongji2/InputForRunBT/RepeatMask/${GENOME}_LTR.bed | intersectBed -s -f 0.5 -wa -a - -b /data/tongji2/piRNA/Output/PiRNA_Loci/LiXin_Pachytene.bed | awk 'BEGIN{sum=0} {sum=sum+($4/$5)} END{print sum}' -`
MULTI_LINE_PACHY_READS=`intersectBed -s -f 0.5 -wa -a ${PREFIX}_multi.bed2 -b /data/tongji2/InputForRunBT/RepeatMask/${GENOME}_LINE.bed | intersectBed -s -f 0.5 -wa -a - -b /data/tongji2/piRNA/Output/PiRNA_Loci/LiXin_Pachytene.bed | awk 'BEGIN{sum=0} {sum=sum+($4/$5)} END{print sum}' -`
MULTI_SINE_PACHY_READS=`intersectBed -s -f 0.5 -wa -a ${PREFIX}_multi.bed2 -b /data/tongji2/InputForRunBT/RepeatMask/${GENOME}_SINE.bed | intersectBed -s -f 0.5 -wa -a - -b /data/tongji2/piRNA/Output/PiRNA_Loci/LiXin_Pachytene.bed | awk 'BEGIN{sum=0} {sum=sum+($4/$5)} END{print sum}' -`
MULTI_LTR_PACHY_READS=`intersectBed -s -f 0.5 -wa -a ${PREFIX}_multi.bed2 -b /data/tongji2/InputForRunBT/RepeatMask/${GENOME}_LTR.bed | intersectBed -s -f 0.5 -wa -a - -b /data/tongji2/piRNA/Output/PiRNA_Loci/LiXin_Pachytene.bed | awk 'BEGIN{sum=0} {sum=sum+($4/$5)} END{print sum}' -`
MULTI_LINE_PACHY_RVS_READS_RVS=`intersectBed -S -f 0.5 -wa -a ${PREFIX}_multi.bed2 -b /data/tongji2/InputForRunBT/RepeatMask/${GENOME}_LINE.bed | intersectBed -S -f 0.5 -wa -a - -b /data/tongji2/piRNA/Output/PiRNA_Loci/LiXin_Pachytene.bed | awk 'BEGIN{sum=0} {sum=sum+($4/$5)} END{print sum}' -`
MULTI_SINE_PACHY_RVS_READS_RVS=`intersectBed -S -f 0.5 -wa -a ${PREFIX}_multi.bed2 -b /data/tongji2/InputForRunBT/RepeatMask/${GENOME}_SINE.bed | intersectBed -S -f 0.5 -wa -a - -b /data/tongji2/piRNA/Output/PiRNA_Loci/LiXin_Pachytene.bed | awk 'BEGIN{sum=0} {sum=sum+($4/$5)} END{print sum}' -`
MULTI_LTR_PACHY_RVS_READS_RVS=`intersectBed -S -f 0.5 -wa -a ${PREFIX}_multi.bed2 -b /data/tongji2/InputForRunBT/RepeatMask/${GENOME}_LTR.bed | intersectBed -S -f 0.5 -wa -a - -b /data/tongji2/piRNA/Output/PiRNA_Loci/LiXin_Pachytene.bed | awk 'BEGIN{sum=0} {sum=sum+($4/$5)} END{print sum}' -`
MULTI_LINE_PACHY_READS_RVS=`intersectBed -s -f 0.5 -wa -a ${PREFIX}_multi.bed2 -b /data/tongji2/InputForRunBT/RepeatMask/${GENOME}_LINE.bed | intersectBed -S -f 0.5 -wa -a - -b /data/tongji2/piRNA/Output/PiRNA_Loci/LiXin_Pachytene.bed | awk 'BEGIN{sum=0} {sum=sum+($4/$5)} END{print sum}' -`
MULTI_SINE_PACHY_READS_RVS=`intersectBed -s -f 0.5 -wa -a ${PREFIX}_multi.bed2 -b /data/tongji2/InputForRunBT/RepeatMask/${GENOME}_SINE.bed | intersectBed -S -f 0.5 -wa -a - -b /data/tongji2/piRNA/Output/PiRNA_Loci/LiXin_Pachytene.bed | awk 'BEGIN{sum=0} {sum=sum+($4/$5)} END{print sum}' -`
MULTI_LTR_PACHY_READS_RVS=`intersectBed -s -f 0.5 -wa -a ${PREFIX}_multi.bed2 -b /data/tongji2/InputForRunBT/RepeatMask/${GENOME}_LTR.bed | intersectBed -S -f 0.5 -wa -a - -b /data/tongji2/piRNA/Output/PiRNA_Loci/LiXin_Pachytene.bed | awk 'BEGIN{sum=0} {sum=sum+($4/$5)} END{print sum}' -`
ALL_LINE_PREPAC_RVS_READS=`intersectBed -S -f 0.5 -wa -a ${PREFIX}.bed2 -b /data/tongji2/InputForRunBT/RepeatMask/${GENOME}_LINE.bed | intersectBed -s -f 0.5 -wa -a - -b /data/tongji2/piRNA/Output/PiRNA_Loci/LiXin_Prepachytene.bed | awk 'BEGIN{sum=0} {sum=sum+($4/$5)} END{print sum}' -`
ALL_SINE_PREPAC_RVS_READS=`intersectBed -S -f 0.5 -wa -a ${PREFIX}.bed2 -b /data/tongji2/InputForRunBT/RepeatMask/${GENOME}_SINE.bed | intersectBed -s -f 0.5 -wa -a - -b /data/tongji2/piRNA/Output/PiRNA_Loci/LiXin_Prepachytene.bed | awk 'BEGIN{sum=0} {sum=sum+($4/$5)} END{print sum}' -`
ALL_LTR_PREPAC_RVS_READS=`intersectBed -S -f 0.5 -wa -a ${PREFIX}.bed2 -b /data/tongji2/InputForRunBT/RepeatMask/${GENOME}_LTR.bed | intersectBed -s -f 0.5 -wa -a - -b /data/tongji2/piRNA/Output/PiRNA_Loci/LiXin_Prepachytene.bed | awk 'BEGIN{sum=0} {sum=sum+($4/$5)} END{print sum}' -`
ALL_LINE_PREPAC_READS=`intersectBed -s -f 0.5 -wa -a ${PREFIX}.bed2 -b /data/tongji2/InputForRunBT/RepeatMask/${GENOME}_LINE.bed | intersectBed -s -f 0.5 -wa -a - -b /data/tongji2/piRNA/Output/PiRNA_Loci/LiXin_Prepachytene.bed | awk 'BEGIN{sum=0} {sum=sum+($4/$5)} END{print sum}' -`
ALL_SINE_PREPAC_READS=`intersectBed -s -f 0.5 -wa -a ${PREFIX}.bed2 -b /data/tongji2/InputForRunBT/RepeatMask/${GENOME}_SINE.bed | intersectBed -s -f 0.5 -wa -a - -b /data/tongji2/piRNA/Output/PiRNA_Loci/LiXin_Prepachytene.bed | awk 'BEGIN{sum=0} {sum=sum+($4/$5)} END{print sum}' -`
ALL_LTR_PREPAC_READS=`intersectBed -s -f 0.5 -wa -a ${PREFIX}.bed2 -b /data/tongji2/InputForRunBT/RepeatMask/${GENOME}_LTR.bed | intersectBed -s -f 0.5 -wa -a - -b /data/tongji2/piRNA/Output/PiRNA_Loci/LiXin_Prepachytene.bed | awk 'BEGIN{sum=0} {sum=sum+($4/$5)} END{print sum}' -`
ALL_LINE_PREPAC_RVS_READS_RVS=`intersectBed -S -f 0.5 -wa -a ${PREFIX}.bed2 -b /data/tongji2/InputForRunBT/RepeatMask/${GENOME}_LINE.bed | intersectBed -S -f 0.5 -wa -a - -b /data/tongji2/piRNA/Output/PiRNA_Loci/LiXin_Prepachytene.bed | awk 'BEGIN{sum=0} {sum=sum+($4/$5)} END{print sum}' -`
ALL_SINE_PREPAC_RVS_READS_RVS=`intersectBed -S -f 0.5 -wa -a ${PREFIX}.bed2 -b /data/tongji2/InputForRunBT/RepeatMask/${GENOME}_SINE.bed | intersectBed -S -f 0.5 -wa -a - -b /data/tongji2/piRNA/Output/PiRNA_Loci/LiXin_Prepachytene.bed | awk 'BEGIN{sum=0} {sum=sum+($4/$5)} END{print sum}' -`
ALL_LTR_PREPAC_RVS_READS_RVS=`intersectBed -S -f 0.5 -wa -a ${PREFIX}.bed2 -b /data/tongji2/InputForRunBT/RepeatMask/${GENOME}_LTR.bed | intersectBed -S -f 0.5 -wa -a - -b /data/tongji2/piRNA/Output/PiRNA_Loci/LiXin_Prepachytene.bed | awk 'BEGIN{sum=0} {sum=sum+($4/$5)} END{print sum}' -`
ALL_LINE_PREPAC_READS_RVS=`intersectBed -s -f 0.5 -wa -a ${PREFIX}.bed2 -b /data/tongji2/InputForRunBT/RepeatMask/${GENOME}_LINE.bed | intersectBed -S -f 0.5 -wa -a - -b /data/tongji2/piRNA/Output/PiRNA_Loci/LiXin_Prepachytene.bed | awk 'BEGIN{sum=0} {sum=sum+($4/$5)} END{print sum}' -`
ALL_SINE_PREPAC_READS_RVS=`intersectBed -s -f 0.5 -wa -a ${PREFIX}.bed2 -b /data/tongji2/InputForRunBT/RepeatMask/${GENOME}_SINE.bed | intersectBed -S -f 0.5 -wa -a - -b /data/tongji2/piRNA/Output/PiRNA_Loci/LiXin_Prepachytene.bed | awk 'BEGIN{sum=0} {sum=sum+($4/$5)} END{print sum}' -`
ALL_LTR_PREPAC_READS_RVS=`intersectBed -s -f 0.5 -wa -a ${PREFIX}.bed2 -b /data/tongji2/InputForRunBT/RepeatMask/${GENOME}_LTR.bed | intersectBed -S -f 0.5 -wa -a - -b /data/tongji2/piRNA/Output/PiRNA_Loci/LiXin_Prepachytene.bed | awk 'BEGIN{sum=0} {sum=sum+($4/$5)} END{print sum}' -`
UNIQ_LINE_PREPAC_RVS_READS=`intersectBed -S -f 0.5 -wa -a ${PREFIX}_uniq.bed2 -b /data/tongji2/InputForRunBT/RepeatMask/${GENOME}_LINE.bed | intersectBed -s -f 0.5 -wa -a - -b /data/tongji2/piRNA/Output/PiRNA_Loci/LiXin_Prepachytene.bed | awk 'BEGIN{sum=0} {sum=sum+($4/$5)} END{print sum}' -`
UNIQ_SINE_PREPAC_RVS_READS=`intersectBed -S -f 0.5 -wa -a ${PREFIX}_uniq.bed2 -b /data/tongji2/InputForRunBT/RepeatMask/${GENOME}_SINE.bed | intersectBed -s -f 0.5 -wa -a - -b /data/tongji2/piRNA/Output/PiRNA_Loci/LiXin_Prepachytene.bed | awk 'BEGIN{sum=0} {sum=sum+($4/$5)} END{print sum}' -`
UNIQ_LTR_PREPAC_RVS_READS=`intersectBed -S -f 0.5 -wa -a ${PREFIX}_uniq.bed2 -b /data/tongji2/InputForRunBT/RepeatMask/${GENOME}_LTR.bed | intersectBed -s -f 0.5 -wa -a - -b /data/tongji2/piRNA/Output/PiRNA_Loci/LiXin_Prepachytene.bed | awk 'BEGIN{sum=0} {sum=sum+($4/$5)} END{print sum}' -`
UNIQ_LINE_PREPAC_READS=`intersectBed -s -f 0.5 -wa -a ${PREFIX}_uniq.bed2 -b /data/tongji2/InputForRunBT/RepeatMask/${GENOME}_LINE.bed | intersectBed -s -f 0.5 -wa -a - -b /data/tongji2/piRNA/Output/PiRNA_Loci/LiXin_Prepachytene.bed | awk 'BEGIN{sum=0} {sum=sum+($4/$5)} END{print sum}' -`
UNIQ_SINE_PREPAC_READS=`intersectBed -s -f 0.5 -wa -a ${PREFIX}_uniq.bed2 -b /data/tongji2/InputForRunBT/RepeatMask/${GENOME}_SINE.bed | intersectBed -s -f 0.5 -wa -a - -b /data/tongji2/piRNA/Output/PiRNA_Loci/LiXin_Prepachytene.bed | awk 'BEGIN{sum=0} {sum=sum+($4/$5)} END{print sum}' -`
UNIQ_LTR_PREPAC_READS=`intersectBed -s -f 0.5 -wa -a ${PREFIX}_uniq.bed2 -b /data/tongji2/InputForRunBT/RepeatMask/${GENOME}_LTR.bed | intersectBed -s -f 0.5 -wa -a - -b /data/tongji2/piRNA/Output/PiRNA_Loci/LiXin_Prepachytene.bed | awk 'BEGIN{sum=0} {sum=sum+($4/$5)} END{print sum}' -`
UNIQ_LINE_PREPAC_RVS_READS_RVS=`intersectBed -S -f 0.5 -wa -a ${PREFIX}_uniq.bed2 -b /data/tongji2/InputForRunBT/RepeatMask/${GENOME}_LINE.bed | intersectBed -S -f 0.5 -wa -a - -b /data/tongji2/piRNA/Output/PiRNA_Loci/LiXin_Prepachytene.bed | awk 'BEGIN{sum=0} {sum=sum+($4/$5)} END{print sum}' -`
UNIQ_SINE_PREPAC_RVS_READS_RVS=`intersectBed -S -f 0.5 -wa -a ${PREFIX}_uniq.bed2 -b /data/tongji2/InputForRunBT/RepeatMask/${GENOME}_SINE.bed | intersectBed -S -f 0.5 -wa -a - -b /data/tongji2/piRNA/Output/PiRNA_Loci/LiXin_Prepachytene.bed | awk 'BEGIN{sum=0} {sum=sum+($4/$5)} END{print sum}' -`
UNIQ_LTR_PREPAC_RVS_READS_RVS=`intersectBed -S -f 0.5 -wa -a ${PREFIX}_uniq.bed2 -b /data/tongji2/InputForRunBT/RepeatMask/${GENOME}_LTR.bed | intersectBed -S -f 0.5 -wa -a - -b /data/tongji2/piRNA/Output/PiRNA_Loci/LiXin_Prepachytene.bed | awk 'BEGIN{sum=0} {sum=sum+($4/$5)} END{print sum}' -`
UNIQ_LINE_PREPAC_READS_RVS=`intersectBed -s -f 0.5 -wa -a ${PREFIX}_uniq.bed2 -b /data/tongji2/InputForRunBT/RepeatMask/${GENOME}_LINE.bed | intersectBed -S -f 0.5 -wa -a - -b /data/tongji2/piRNA/Output/PiRNA_Loci/LiXin_Prepachytene.bed | awk 'BEGIN{sum=0} {sum=sum+($4/$5)} END{print sum}' -`
UNIQ_SINE_PREPAC_READS_RVS=`intersectBed -s -f 0.5 -wa -a ${PREFIX}_uniq.bed2 -b /data/tongji2/InputForRunBT/RepeatMask/${GENOME}_SINE.bed | intersectBed -S -f 0.5 -wa -a - -b /data/tongji2/piRNA/Output/PiRNA_Loci/LiXin_Prepachytene.bed | awk 'BEGIN{sum=0} {sum=sum+($4/$5)} END{print sum}' -`
UNIQ_LTR_PREPAC_READS_RVS=`intersectBed -s -f 0.5 -wa -a ${PREFIX}_uniq.bed2 -b /data/tongji2/InputForRunBT/RepeatMask/${GENOME}_LTR.bed | intersectBed -S -f 0.5 -wa -a - -b /data/tongji2/piRNA/Output/PiRNA_Loci/LiXin_Prepachytene.bed | awk 'BEGIN{sum=0} {sum=sum+($4/$5)} END{print sum}' -`
MULTI_LINE_PREPAC_RVS_READS=`intersectBed -S -f 0.5 -wa -a ${PREFIX}_multi.bed2 -b /data/tongji2/InputForRunBT/RepeatMask/${GENOME}_LINE.bed | intersectBed -s -f 0.5 -wa -a - -b /data/tongji2/piRNA/Output/PiRNA_Loci/LiXin_Prepachytene.bed | awk 'BEGIN{sum=0} {sum=sum+($4/$5)} END{print sum}' -`
MULTI_SINE_PREPAC_RVS_READS=`intersectBed -S -f 0.5 -wa -a ${PREFIX}_multi.bed2 -b /data/tongji2/InputForRunBT/RepeatMask/${GENOME}_SINE.bed | intersectBed -s -f 0.5 -wa -a - -b /data/tongji2/piRNA/Output/PiRNA_Loci/LiXin_Prepachytene.bed | awk 'BEGIN{sum=0} {sum=sum+($4/$5)} END{print sum}' -`
MULTI_LTR_PREPAC_RVS_READS=`intersectBed -S -f 0.5 -wa -a ${PREFIX}_multi.bed2 -b /data/tongji2/InputForRunBT/RepeatMask/${GENOME}_LTR.bed | intersectBed -s -f 0.5 -wa -a - -b /data/tongji2/piRNA/Output/PiRNA_Loci/LiXin_Prepachytene.bed | awk 'BEGIN{sum=0} {sum=sum+($4/$5)} END{print sum}' -`
MULTI_LINE_PREPAC_READS=`intersectBed -s -f 0.5 -wa -a ${PREFIX}_multi.bed2 -b /data/tongji2/InputForRunBT/RepeatMask/${GENOME}_LINE.bed | intersectBed -s -f 0.5 -wa -a - -b /data/tongji2/piRNA/Output/PiRNA_Loci/LiXin_Prepachytene.bed | awk 'BEGIN{sum=0} {sum=sum+($4/$5)} END{print sum}' -`
MULTI_SINE_PREPAC_READS=`intersectBed -s -f 0.5 -wa -a ${PREFIX}_multi.bed2 -b /data/tongji2/InputForRunBT/RepeatMask/${GENOME}_SINE.bed | intersectBed -s -f 0.5 -wa -a - -b /data/tongji2/piRNA/Output/PiRNA_Loci/LiXin_Prepachytene.bed | awk 'BEGIN{sum=0} {sum=sum+($4/$5)} END{print sum}' -`
MULTI_LTR_PREPAC_READS=`intersectBed -s -f 0.5 -wa -a ${PREFIX}_multi.bed2 -b /data/tongji2/InputForRunBT/RepeatMask/${GENOME}_LTR.bed | intersectBed -s -f 0.5 -wa -a - -b /data/tongji2/piRNA/Output/PiRNA_Loci/LiXin_Prepachytene.bed | awk 'BEGIN{sum=0} {sum=sum+($4/$5)} END{print sum}' -`
MULTI_LINE_PREPAC_RVS_READS_RVS=`intersectBed -S -f 0.5 -wa -a ${PREFIX}_multi.bed2 -b /data/tongji2/InputForRunBT/RepeatMask/${GENOME}_LINE.bed | intersectBed -S -f 0.5 -wa -a - -b /data/tongji2/piRNA/Output/PiRNA_Loci/LiXin_Prepachytene.bed | awk 'BEGIN{sum=0} {sum=sum+($4/$5)} END{print sum}' -`
MULTI_SINE_PREPAC_RVS_READS_RVS=`intersectBed -S -f 0.5 -wa -a ${PREFIX}_multi.bed2 -b /data/tongji2/InputForRunBT/RepeatMask/${GENOME}_SINE.bed | intersectBed -S -f 0.5 -wa -a - -b /data/tongji2/piRNA/Output/PiRNA_Loci/LiXin_Prepachytene.bed | awk 'BEGIN{sum=0} {sum=sum+($4/$5)} END{print sum}' -`
MULTI_LTR_PREPAC_RVS_READS_RVS=`intersectBed -S -f 0.5 -wa -a ${PREFIX}_multi.bed2 -b /data/tongji2/InputForRunBT/RepeatMask/${GENOME}_LTR.bed | intersectBed -S -f 0.5 -wa -a - -b /data/tongji2/piRNA/Output/PiRNA_Loci/LiXin_Prepachytene.bed | awk 'BEGIN{sum=0} {sum=sum+($4/$5)} END{print sum}' -`
MULTI_LINE_PREPAC_READS_RVS=`intersectBed -s -f 0.5 -wa -a ${PREFIX}_multi.bed2 -b /data/tongji2/InputForRunBT/RepeatMask/${GENOME}_LINE.bed | intersectBed -S -f 0.5 -wa -a - -b /data/tongji2/piRNA/Output/PiRNA_Loci/LiXin_Prepachytene.bed | awk 'BEGIN{sum=0} {sum=sum+($4/$5)} END{print sum}' -`
MULTI_SINE_PREPAC_READS_RVS=`intersectBed -s -f 0.5 -wa -a ${PREFIX}_multi.bed2 -b /data/tongji2/InputForRunBT/RepeatMask/${GENOME}_SINE.bed | intersectBed -S -f 0.5 -wa -a - -b /data/tongji2/piRNA/Output/PiRNA_Loci/LiXin_Prepachytene.bed | awk 'BEGIN{sum=0} {sum=sum+($4/$5)} END{print sum}' -`
MULTI_LTR_PREPAC_READS_RVS=`intersectBed -s -f 0.5 -wa -a ${PREFIX}_multi.bed2 -b /data/tongji2/InputForRunBT/RepeatMask/${GENOME}_LTR.bed | intersectBed -S -f 0.5 -wa -a - -b /data/tongji2/piRNA/Output/PiRNA_Loci/LiXin_Prepachytene.bed | awk 'BEGIN{sum=0} {sum=sum+($4/$5)} END{print sum}' -`
ALL_LINE_HYBRID_RVS_READS=`intersectBed -S -f 0.5 -wa -a ${PREFIX}.bed2 -b /data/tongji2/InputForRunBT/RepeatMask/${GENOME}_LINE.bed | intersectBed -s -f 0.5 -wa -a - -b /data/tongji2/piRNA/Output/PiRNA_Loci/LiXin_Hybrid.bed | awk 'BEGIN{sum=0} {sum=sum+($4/$5)} END{print sum}' -`
ALL_SINE_HYBRID_RVS_READS=`intersectBed -S -f 0.5 -wa -a ${PREFIX}.bed2 -b /data/tongji2/InputForRunBT/RepeatMask/${GENOME}_SINE.bed | intersectBed -s -f 0.5 -wa -a - -b /data/tongji2/piRNA/Output/PiRNA_Loci/LiXin_Hybrid.bed | awk 'BEGIN{sum=0} {sum=sum+($4/$5)} END{print sum}' -`
ALL_LTR_HYBRID_RVS_READS=`intersectBed -S -f 0.5 -wa -a ${PREFIX}.bed2 -b /data/tongji2/InputForRunBT/RepeatMask/${GENOME}_LTR.bed | intersectBed -s -f 0.5 -wa -a - -b /data/tongji2/piRNA/Output/PiRNA_Loci/LiXin_Hybrid.bed | awk 'BEGIN{sum=0} {sum=sum+($4/$5)} END{print sum}' -`
ALL_LINE_HYBRID_READS=`intersectBed -s -f 0.5 -wa -a ${PREFIX}.bed2 -b /data/tongji2/InputForRunBT/RepeatMask/${GENOME}_LINE.bed | intersectBed -s -f 0.5 -wa -a - -b /data/tongji2/piRNA/Output/PiRNA_Loci/LiXin_Hybrid.bed | awk 'BEGIN{sum=0} {sum=sum+($4/$5)} END{print sum}' -`
ALL_SINE_HYBRID_READS=`intersectBed -s -f 0.5 -wa -a ${PREFIX}.bed2 -b /data/tongji2/InputForRunBT/RepeatMask/${GENOME}_SINE.bed | intersectBed -s -f 0.5 -wa -a - -b /data/tongji2/piRNA/Output/PiRNA_Loci/LiXin_Hybrid.bed | awk 'BEGIN{sum=0} {sum=sum+($4/$5)} END{print sum}' -`
ALL_LTR_HYBRID_READS=`intersectBed -s -f 0.5 -wa -a ${PREFIX}.bed2 -b /data/tongji2/InputForRunBT/RepeatMask/${GENOME}_LTR.bed | intersectBed -s -f 0.5 -wa -a - -b /data/tongji2/piRNA/Output/PiRNA_Loci/LiXin_Hybrid.bed | awk 'BEGIN{sum=0} {sum=sum+($4/$5)} END{print sum}' -`
ALL_LINE_HYBRID_RVS_READS_RVS=`intersectBed -S -f 0.5 -wa -a ${PREFIX}.bed2 -b /data/tongji2/InputForRunBT/RepeatMask/${GENOME}_LINE.bed | intersectBed -S -f 0.5 -wa -a - -b /data/tongji2/piRNA/Output/PiRNA_Loci/LiXin_Hybrid.bed | awk 'BEGIN{sum=0} {sum=sum+($4/$5)} END{print sum}' -`
ALL_SINE_HYBRID_RVS_READS_RVS=`intersectBed -S -f 0.5 -wa -a ${PREFIX}.bed2 -b /data/tongji2/InputForRunBT/RepeatMask/${GENOME}_SINE.bed | intersectBed -S -f 0.5 -wa -a - -b /data/tongji2/piRNA/Output/PiRNA_Loci/LiXin_Hybrid.bed | awk 'BEGIN{sum=0} {sum=sum+($4/$5)} END{print sum}' -`
ALL_LTR_HYBRID_RVS_READS_RVS=`intersectBed -S -f 0.5 -wa -a ${PREFIX}.bed2 -b /data/tongji2/InputForRunBT/RepeatMask/${GENOME}_LTR.bed | intersectBed -S -f 0.5 -wa -a - -b /data/tongji2/piRNA/Output/PiRNA_Loci/LiXin_Hybrid.bed | awk 'BEGIN{sum=0} {sum=sum+($4/$5)} END{print sum}' -`
ALL_LINE_HYBRID_READS_RVS=`intersectBed -s -f 0.5 -wa -a ${PREFIX}.bed2 -b /data/tongji2/InputForRunBT/RepeatMask/${GENOME}_LINE.bed | intersectBed -S -f 0.5 -wa -a - -b /data/tongji2/piRNA/Output/PiRNA_Loci/LiXin_Hybrid.bed | awk 'BEGIN{sum=0} {sum=sum+($4/$5)} END{print sum}' -`
ALL_SINE_HYBRID_READS_RVS=`intersectBed -s -f 0.5 -wa -a ${PREFIX}.bed2 -b /data/tongji2/InputForRunBT/RepeatMask/${GENOME}_SINE.bed | intersectBed -S -f 0.5 -wa -a - -b /data/tongji2/piRNA/Output/PiRNA_Loci/LiXin_Hybrid.bed | awk 'BEGIN{sum=0} {sum=sum+($4/$5)} END{print sum}' -`
ALL_LTR_HYBRID_READS_RVS=`intersectBed -s -f 0.5 -wa -a ${PREFIX}.bed2 -b /data/tongji2/InputForRunBT/RepeatMask/${GENOME}_LTR.bed | intersectBed -S -f 0.5 -wa -a - -b /data/tongji2/piRNA/Output/PiRNA_Loci/LiXin_Hybrid.bed | awk 'BEGIN{sum=0} {sum=sum+($4/$5)} END{print sum}' -`
UNIQ_LINE_HYBRID_RVS_READS=`intersectBed -S -f 0.5 -wa -a ${PREFIX}_uniq.bed2 -b /data/tongji2/InputForRunBT/RepeatMask/${GENOME}_LINE.bed | intersectBed -s -f 0.5 -wa -a - -b /data/tongji2/piRNA/Output/PiRNA_Loci/LiXin_Hybrid.bed | awk 'BEGIN{sum=0} {sum=sum+($4/$5)} END{print sum}' -`
UNIQ_SINE_HYBRID_RVS_READS=`intersectBed -S -f 0.5 -wa -a ${PREFIX}_uniq.bed2 -b /data/tongji2/InputForRunBT/RepeatMask/${GENOME}_SINE.bed | intersectBed -s -f 0.5 -wa -a - -b /data/tongji2/piRNA/Output/PiRNA_Loci/LiXin_Hybrid.bed | awk 'BEGIN{sum=0} {sum=sum+($4/$5)} END{print sum}' -`
UNIQ_LTR_HYBRID_RVS_READS=`intersectBed -S -f 0.5 -wa -a ${PREFIX}_uniq.bed2 -b /data/tongji2/InputForRunBT/RepeatMask/${GENOME}_LTR.bed | intersectBed -s -f 0.5 -wa -a - -b /data/tongji2/piRNA/Output/PiRNA_Loci/LiXin_Hybrid.bed | awk 'BEGIN{sum=0} {sum=sum+($4/$5)} END{print sum}' -`
UNIQ_LINE_HYBRID_READS=`intersectBed -s -f 0.5 -wa -a ${PREFIX}_uniq.bed2 -b /data/tongji2/InputForRunBT/RepeatMask/${GENOME}_LINE.bed | intersectBed -s -f 0.5 -wa -a - -b /data/tongji2/piRNA/Output/PiRNA_Loci/LiXin_Hybrid.bed | awk 'BEGIN{sum=0} {sum=sum+($4/$5)} END{print sum}' -`
UNIQ_SINE_HYBRID_READS=`intersectBed -s -f 0.5 -wa -a ${PREFIX}_uniq.bed2 -b /data/tongji2/InputForRunBT/RepeatMask/${GENOME}_SINE.bed | intersectBed -s -f 0.5 -wa -a - -b /data/tongji2/piRNA/Output/PiRNA_Loci/LiXin_Hybrid.bed | awk 'BEGIN{sum=0} {sum=sum+($4/$5)} END{print sum}' -`
UNIQ_LTR_HYBRID_READS=`intersectBed -s -f 0.5 -wa -a ${PREFIX}_uniq.bed2 -b /data/tongji2/InputForRunBT/RepeatMask/${GENOME}_LTR.bed | intersectBed -s -f 0.5 -wa -a - -b /data/tongji2/piRNA/Output/PiRNA_Loci/LiXin_Hybrid.bed | awk 'BEGIN{sum=0} {sum=sum+($4/$5)} END{print sum}' -`
UNIQ_LINE_HYBRID_RVS_READS_RVS=`intersectBed -S -f 0.5 -wa -a ${PREFIX}_uniq.bed2 -b /data/tongji2/InputForRunBT/RepeatMask/${GENOME}_LINE.bed | intersectBed -S -f 0.5 -wa -a - -b /data/tongji2/piRNA/Output/PiRNA_Loci/LiXin_Hybrid.bed | awk 'BEGIN{sum=0} {sum=sum+($4/$5)} END{print sum}' -`
UNIQ_SINE_HYBRID_RVS_READS_RVS=`intersectBed -S -f 0.5 -wa -a ${PREFIX}_uniq.bed2 -b /data/tongji2/InputForRunBT/RepeatMask/${GENOME}_SINE.bed | intersectBed -S -f 0.5 -wa -a - -b /data/tongji2/piRNA/Output/PiRNA_Loci/LiXin_Hybrid.bed | awk 'BEGIN{sum=0} {sum=sum+($4/$5)} END{print sum}' -`
UNIQ_LTR_HYBRID_RVS_READS_RVS=`intersectBed -S -f 0.5 -wa -a ${PREFIX}_uniq.bed2 -b /data/tongji2/InputForRunBT/RepeatMask/${GENOME}_LTR.bed | intersectBed -S -f 0.5 -wa -a - -b /data/tongji2/piRNA/Output/PiRNA_Loci/LiXin_Hybrid.bed | awk 'BEGIN{sum=0} {sum=sum+($4/$5)} END{print sum}' -`
UNIQ_LINE_HYBRID_READS_RVS=`intersectBed -s -f 0.5 -wa -a ${PREFIX}_uniq.bed2 -b /data/tongji2/InputForRunBT/RepeatMask/${GENOME}_LINE.bed | intersectBed -S -f 0.5 -wa -a - -b /data/tongji2/piRNA/Output/PiRNA_Loci/LiXin_Hybrid.bed | awk 'BEGIN{sum=0} {sum=sum+($4/$5)} END{print sum}' -`
UNIQ_SINE_HYBRID_READS_RVS=`intersectBed -s -f 0.5 -wa -a ${PREFIX}_uniq.bed2 -b /data/tongji2/InputForRunBT/RepeatMask/${GENOME}_SINE.bed | intersectBed -S -f 0.5 -wa -a - -b /data/tongji2/piRNA/Output/PiRNA_Loci/LiXin_Hybrid.bed | awk 'BEGIN{sum=0} {sum=sum+($4/$5)} END{print sum}' -`
UNIQ_LTR_HYBRID_READS_RVS=`intersectBed -s -f 0.5 -wa -a ${PREFIX}_uniq.bed2 -b /data/tongji2/InputForRunBT/RepeatMask/${GENOME}_LTR.bed | intersectBed -S -f 0.5 -wa -a - -b /data/tongji2/piRNA/Output/PiRNA_Loci/LiXin_Hybrid.bed | awk 'BEGIN{sum=0} {sum=sum+($4/$5)} END{print sum}' -`
MULTI_LINE_HYBRID_RVS_READS=`intersectBed -S -f 0.5 -wa -a ${PREFIX}_multi.bed2 -b /data/tongji2/InputForRunBT/RepeatMask/${GENOME}_LINE.bed | intersectBed -s -f 0.5 -wa -a - -b /data/tongji2/piRNA/Output/PiRNA_Loci/LiXin_Hybrid.bed | awk 'BEGIN{sum=0} {sum=sum+($4/$5)} END{print sum}' -`
MULTI_SINE_HYBRID_RVS_READS=`intersectBed -S -f 0.5 -wa -a ${PREFIX}_multi.bed2 -b /data/tongji2/InputForRunBT/RepeatMask/${GENOME}_SINE.bed | intersectBed -s -f 0.5 -wa -a - -b /data/tongji2/piRNA/Output/PiRNA_Loci/LiXin_Hybrid.bed | awk 'BEGIN{sum=0} {sum=sum+($4/$5)} END{print sum}' -`
MULTI_LTR_HYBRID_RVS_READS=`intersectBed -S -f 0.5 -wa -a ${PREFIX}_multi.bed2 -b /data/tongji2/InputForRunBT/RepeatMask/${GENOME}_LTR.bed | intersectBed -s -f 0.5 -wa -a - -b /data/tongji2/piRNA/Output/PiRNA_Loci/LiXin_Hybrid.bed | awk 'BEGIN{sum=0} {sum=sum+($4/$5)} END{print sum}' -`
MULTI_LINE_HYBRID_READS=`intersectBed -s -f 0.5 -wa -a ${PREFIX}_multi.bed2 -b /data/tongji2/InputForRunBT/RepeatMask/${GENOME}_LINE.bed | intersectBed -s -f 0.5 -wa -a - -b /data/tongji2/piRNA/Output/PiRNA_Loci/LiXin_Hybrid.bed | awk 'BEGIN{sum=0} {sum=sum+($4/$5)} END{print sum}' -`
MULTI_SINE_HYBRID_READS=`intersectBed -s -f 0.5 -wa -a ${PREFIX}_multi.bed2 -b /data/tongji2/InputForRunBT/RepeatMask/${GENOME}_SINE.bed | intersectBed -s -f 0.5 -wa -a - -b /data/tongji2/piRNA/Output/PiRNA_Loci/LiXin_Hybrid.bed | awk 'BEGIN{sum=0} {sum=sum+($4/$5)} END{print sum}' -`
MULTI_LTR_HYBRID_READS=`intersectBed -s -f 0.5 -wa -a ${PREFIX}_multi.bed2 -b /data/tongji2/InputForRunBT/RepeatMask/${GENOME}_LTR.bed | intersectBed -s -f 0.5 -wa -a - -b /data/tongji2/piRNA/Output/PiRNA_Loci/LiXin_Hybrid.bed | awk 'BEGIN{sum=0} {sum=sum+($4/$5)} END{print sum}' -`
MULTI_LINE_HYBRID_RVS_READS_RVS=`intersectBed -S -f 0.5 -wa -a ${PREFIX}_multi.bed2 -b /data/tongji2/InputForRunBT/RepeatMask/${GENOME}_LINE.bed | intersectBed -S -f 0.5 -wa -a - -b /data/tongji2/piRNA/Output/PiRNA_Loci/LiXin_Hybrid.bed | awk 'BEGIN{sum=0} {sum=sum+($4/$5)} END{print sum}' -`
MULTI_SINE_HYBRID_RVS_READS_RVS=`intersectBed -S -f 0.5 -wa -a ${PREFIX}_multi.bed2 -b /data/tongji2/InputForRunBT/RepeatMask/${GENOME}_SINE.bed | intersectBed -S -f 0.5 -wa -a - -b /data/tongji2/piRNA/Output/PiRNA_Loci/LiXin_Hybrid.bed | awk 'BEGIN{sum=0} {sum=sum+($4/$5)} END{print sum}' -`
MULTI_LTR_HYBRID_RVS_READS_RVS=`intersectBed -S -f 0.5 -wa -a ${PREFIX}_multi.bed2 -b /data/tongji2/InputForRunBT/RepeatMask/${GENOME}_LTR.bed | intersectBed -S -f 0.5 -wa -a - -b /data/tongji2/piRNA/Output/PiRNA_Loci/LiXin_Hybrid.bed | awk 'BEGIN{sum=0} {sum=sum+($4/$5)} END{print sum}' -`
MULTI_LINE_HYBRID_READS_RVS=`intersectBed -s -f 0.5 -wa -a ${PREFIX}_multi.bed2 -b /data/tongji2/InputForRunBT/RepeatMask/${GENOME}_LINE.bed | intersectBed -S -f 0.5 -wa -a - -b /data/tongji2/piRNA/Output/PiRNA_Loci/LiXin_Hybrid.bed | awk 'BEGIN{sum=0} {sum=sum+($4/$5)} END{print sum}' -`
MULTI_SINE_HYBRID_READS_RVS=`intersectBed -s -f 0.5 -wa -a ${PREFIX}_multi.bed2 -b /data/tongji2/InputForRunBT/RepeatMask/${GENOME}_SINE.bed | intersectBed -S -f 0.5 -wa -a - -b /data/tongji2/piRNA/Output/PiRNA_Loci/LiXin_Hybrid.bed | awk 'BEGIN{sum=0} {sum=sum+($4/$5)} END{print sum}' -`
MULTI_LTR_HYBRID_READS_RVS=`intersectBed -s -f 0.5 -wa -a ${PREFIX}_multi.bed2 -b /data/tongji2/InputForRunBT/RepeatMask/${GENOME}_LTR.bed | intersectBed -S -f 0.5 -wa -a - -b /data/tongji2/piRNA/Output/PiRNA_Loci/LiXin_Hybrid.bed | awk 'BEGIN{sum=0} {sum=sum+($4/$5)} END{print sum}' -`

intersectBed -s -F 0.5 -wo -a /data/tongji2/piRNA/Output/PiRNA_Loci/LiXin_Pachytene.bed -b \
${PREFIX}.bed2 > ./Intersect_Result/All_Pachy.txt
intersectBed -s -F 0.5 -wo -a /data/tongji2/piRNA/Output/PiRNA_Loci/LiXin_Pachytene.bed -b \
${PREFIX}_uniq.bed2 > ./Intersect_Result/Uniq_Pachy.txt
intersectBed -s -F 0.5 -wo -a /data/tongji2/piRNA/Output/PiRNA_Loci/LiXin_Pachytene.bed -b \
${PREFIX}_multi.bed2 > ./Intersect_Result/Multi_Pachy.txt
intersectBed -s -F 0.5 -wo -a /data/tongji2/piRNA/Output/PiRNA_Loci/LiXin_Prepachytene.bed -b \
${PREFIX}.bed2 > ./Intersect_Result/All_Prepachy.txt
intersectBed -s -F 0.5 -wo -a /data/tongji2/piRNA/Output/PiRNA_Loci/LiXin_Prepachytene.bed -b \
${PREFIX}_uniq.bed2 > ./Intersect_Result/Uniq_Prepachy.txt
intersectBed -s -F 0.5 -wo -a /data/tongji2/piRNA/Output/PiRNA_Loci/LiXin_Prepachytene.bed -b \
${PREFIX}_multi.bed2 > ./Intersect_Result/Multi_Prepachy.txt
intersectBed -s -F 0.5 -wo -a /data/tongji2/piRNA/Output/PiRNA_Loci/LiXin_Hybrid.bed -b \
${PREFIX}.bed2 > ./Intersect_Result/All_Hybrid.txt
intersectBed -s -F 0.5 -wo -a /data/tongji2/piRNA/Output/PiRNA_Loci/LiXin_Hybrid.bed -b \
${PREFIX}_uniq.bed2 > ./Intersect_Result/Uniq_Hybrid.txt
intersectBed -s -F 0.5 -wo -a /data/tongji2/piRNA/Output/PiRNA_Loci/LiXin_Hybrid.bed -b \
${PREFIX}_multi.bed2 > ./Intersect_Result/Multi_Hybrid.txt

ALL_LINE_READS=`awk 'BEGIN{sum=0} {sum=sum+($10/$11)} END{print sum}' ./Intersect_Result/All_LINE.txt`
ALL_SINE_READS=`awk 'BEGIN{sum=0} {sum=sum+($10/$11)} END{print sum}' ./Intersect_Result/All_SINE.txt`
ALL_LTR_READS=`awk 'BEGIN{sum=0} {sum=sum+($10/$11)} END{print sum}' ./Intersect_Result/All_LTR.txt`
ALL_DNA_READS=`awk 'BEGIN{sum=0} {sum=sum+($10/$11)} END{print sum}' ./Intersect_Result/All_DNA.txt`
ALL_Satellite_READS=`awk 'BEGIN{sum=0} {sum=sum+($10/$11)} END{print sum}' ./Intersect_Result/All_Satellite.txt`
ALL_Simple_repeat_READS=`awk 'BEGIN{sum=0} {sum=sum+($10/$11)} END{print sum}' ./Intersect_Result/All_Simple_repeat.txt`
ALL_PACHY_READS=`awk 'BEGIN{sum=0} {sum=sum+($10/$11)} END{print sum}' ./Intersect_Result/All_Pachy.txt`
ALL_PREPAC_READS=`awk 'BEGIN{sum=0} {sum=sum+($10/$11)} END{print sum}' ./Intersect_Result/All_Prepachy.txt`
ALL_HYBRID_READS=`awk 'BEGIN{sum=0} {sum=sum+($10/$11)} END{print sum}' ./Intersect_Result/All_Hybrid.txt`
UNIQ_LINE_READS=`awk 'BEGIN{sum=0} {sum=sum+($10/$11)} END{print sum}' ./Intersect_Result/Uniq_LINE.txt`
UNIQ_SINE_READS=`awk 'BEGIN{sum=0} {sum=sum+($10/$11)} END{print sum}' ./Intersect_Result/Uniq_SINE.txt`
UNIQ_LTR_READS=`awk 'BEGIN{sum=0} {sum=sum+($10/$11)} END{print sum}' ./Intersect_Result/Uniq_LTR.txt`
UNIQ_DNA_READS=`awk 'BEGIN{sum=0} {sum=sum+($10/$11)} END{print sum}' ./Intersect_Result/Uniq_DNA.txt`
UNIQ_Satellite_READS=`awk 'BEGIN{sum=0} {sum=sum+($10/$11)} END{print sum}' ./Intersect_Result/Uniq_Satellite.txt`
UNIQ_Simple_repeat_READS=`awk 'BEGIN{sum=0} {sum=sum+($10/$11)} END{print sum}' ./Intersect_Result/Uniq_Simple_repeat.txt`
UNIQ_PACHY_READS=`awk 'BEGIN{sum=0} {sum=sum+($10/$11)} END{print sum}' ./Intersect_Result/Uniq_Pachy.txt`
UNIQ_PREPAC_READS=`awk 'BEGIN{sum=0} {sum=sum+($10/$11)} END{print sum}' ./Intersect_Result/Uniq_Prepachy.txt`
UNIQ_HYBRID_READS=`awk 'BEGIN{sum=0} {sum=sum+($10/$11)} END{print sum}' ./Intersect_Result/Uniq_Hybrid.txt`
MULTI_LINE_READS=`awk 'BEGIN{sum=0} {sum=sum+($10/$11)} END{print sum}' ./Intersect_Result/Multi_LINE.txt`
MULTI_SINE_READS=`awk 'BEGIN{sum=0} {sum=sum+($10/$11)} END{print sum}' ./Intersect_Result/Multi_SINE.txt`
MULTI_LTR_READS=`awk 'BEGIN{sum=0} {sum=sum+($10/$11)} END{print sum}' ./Intersect_Result/Multi_LTR.txt`
MULTI_DNA_READS=`awk 'BEGIN{sum=0} {sum=sum+($10/$11)} END{print sum}' ./Intersect_Result/Multi_DNA.txt`
MULTI_Satellite_READS=`awk 'BEGIN{sum=0} {sum=sum+($10/$11)} END{print sum}' ./Intersect_Result/Multi_Satellite.txt`
MULTI_Simple_repeat_READS=`awk 'BEGIN{sum=0} {sum=sum+($10/$11)} END{print sum}' ./Intersect_Result/Multi_Simple_repeat.txt`
MULTI_PACHY_READS=`awk 'BEGIN{sum=0} {sum=sum+($10/$11)} END{print sum}' ./Intersect_Result/Multi_Pachy.txt`
MULTI_PREPAC_READS=`awk 'BEGIN{sum=0} {sum=sum+($10/$11)} END{print sum}' ./Intersect_Result/Multi_Prepachy.txt`
MULTI_HYBRID_READS=`awk 'BEGIN{sum=0} {sum=sum+($10/$11)} END{print sum}' ./Intersect_Result/Multi_Hybrid.txt`
ALL_LINE_READS_RVS=`awk 'BEGIN{sum=0} {sum=sum+($10/$11)} END{print sum}' ./Intersect_Result/All_LINE_RVS.txt`
ALL_SINE_READS_RVS=`awk 'BEGIN{sum=0} {sum=sum+($10/$11)} END{print sum}' ./Intersect_Result/All_SINE_RVS.txt`
ALL_LTR_READS_RVS=`awk 'BEGIN{sum=0} {sum=sum+($10/$11)} END{print sum}' ./Intersect_Result/All_LTR_RVS.txt`
ALL_DNA_READS_RVS=`awk 'BEGIN{sum=0} {sum=sum+($10/$11)} END{print sum}' ./Intersect_Result/All_DNA_RVS.txt`
ALL_Satellite_READS_RVS=`awk 'BEGIN{sum=0} {sum=sum+($10/$11)} END{print sum}' ./Intersect_Result/All_Satellite_RVS.txt`
ALL_Simple_repeat_READS_RVS=`awk 'BEGIN{sum=0} {sum=sum+($10/$11)} END{print sum}' ./Intersect_Result/All_Simple_repeat_RVS.txt`
UNIQ_LINE_READS_RVS=`awk 'BEGIN{sum=0} {sum=sum+($10/$11)} END{print sum}' ./Intersect_Result/Uniq_LINE_RVS.txt`
UNIQ_SINE_READS_RVS=`awk 'BEGIN{sum=0} {sum=sum+($10/$11)} END{print sum}' ./Intersect_Result/Uniq_SINE_RVS.txt`
UNIQ_LTR_READS_RVS=`awk 'BEGIN{sum=0} {sum=sum+($10/$11)} END{print sum}' ./Intersect_Result/Uniq_LTR_RVS.txt`
UNIQ_DNA_READS_RVS=`awk 'BEGIN{sum=0} {sum=sum+($10/$11)} END{print sum}' ./Intersect_Result/Uniq_DNA_RVS.txt`
UNIQ_Satellite_READS_RVS=`awk 'BEGIN{sum=0} {sum=sum+($10/$11)} END{print sum}' ./Intersect_Result/Uniq_Satellite_RVS.txt`
UNIQ_Simple_repeat_READS_RVS=`awk 'BEGIN{sum=0} {sum=sum+($10/$11)} END{print sum}' ./Intersect_Result/Uniq_Simple_repeat_RVS.txt`
MULTI_LINE_READS_RVS=`awk 'BEGIN{sum=0} {sum=sum+($10/$11)} END{print sum}' ./Intersect_Result/Multi_LINE_RVS.txt`
MULTI_SINE_READS_RVS=`awk 'BEGIN{sum=0} {sum=sum+($10/$11)} END{print sum}' ./Intersect_Result/Multi_SINE_RVS.txt`
MULTI_LTR_READS_RVS=`awk 'BEGIN{sum=0} {sum=sum+($10/$11)} END{print sum}' ./Intersect_Result/Multi_LTR_RVS.txt`
MULTI_DNA_READS_RVS=`awk 'BEGIN{sum=0} {sum=sum+($10/$11)} END{print sum}' ./Intersect_Result/Multi_DNA_RVS.txt`
MULTI_Satellite_READS_RVS=`awk 'BEGIN{sum=0} {sum=sum+($10/$11)} END{print sum}' ./Intersect_Result/Multi_Satellite_RVS.txt`
MULTI_Simple_repeat_READS_RVS=`awk 'BEGIN{sum=0} {sum=sum+($10/$11)} END{print sum}' ./Intersect_Result/Multi_Simple_repeat_RVS.txt`

###output summary
echo -e "\ttotal_reads\ttotal_map_reads\tunique_map_reads\tmultiple_map_reads" > map_summary.txt
echo -e "All_Feature\t${TOTAL_READS}\t${TOTAL_MAP_READS}\t${UNIQ_MAP_READS}\t${MULTI_MAP_READS}" >> map_summary.txt
echo -e "LINE\t${TOTAL_READS}\t${ALL_LINE_READS}\t${UNIQ_LINE_READS}\t${MULTI_LINE_READS}" >> map_summary.txt
echo -e "SINE\t${TOTAL_READS}\t${ALL_SINE_READS}\t${UNIQ_SINE_READS}\t${MULTI_SINE_READS}" >> map_summary.txt
echo -e "LTR\t${TOTAL_READS}\t${ALL_LTR_READS}\t${UNIQ_LTR_READS}\t${MULTI_LTR_READS}" >> map_summary.txt
echo -e "DNA\t${TOTAL_READS}\t${ALL_DNA_READS}\t${UNIQ_DNA_READS}\t${MULTI_DNA_READS}" >> map_summary.txt
echo -e "Satellite\t${TOTAL_READS}\t${ALL_Satellite_READS}\t${UNIQ_Satellite_READS}\t${MULTI_Satellite_READS}" >> map_summary.txt
echo -e "Simple_repeat\t${TOTAL_READS}\t${ALL_Simple_repeat_READS}\t${UNIQ_Simple_repeat_READS}\t${MULTI_Simple_repeat_READS}" >> map_summary.txt
echo -e "LINE_Reverse\t${TOTAL_READS}\t${ALL_LINE_READS_RVS}\t${UNIQ_LINE_READS_RVS}\t${MULTI_LINE_READS_RVS}" >> map_summary.txt
echo -e "SINE_Reverse\t${TOTAL_READS}\t${ALL_SINE_READS_RVS}\t${UNIQ_SINE_READS_RVS}\t${MULTI_SINE_READS_RVS}" >> map_summary.txt
echo -e "LTR_Reverse\t${TOTAL_READS}\t${ALL_LTR_READS_RVS}\t${UNIQ_LTR_READS_RVS}\t${MULTI_LTR_READS_RVS}" >> map_summary.txt
echo -e "DNA_Reverse\t${TOTAL_READS}\t${ALL_DNA_READS_RVS}\t${UNIQ_DNA_READS_RVS}\t${MULTI_DNA_READS_RVS}" >> map_summary.txt
echo -e "Satellite_Reverse\t${TOTAL_READS}\t${ALL_Satellite_READS_RVS}\t${UNIQ_Satellite_READS_RVS}\t${MULTI_Satellite_READS_RVS}" >> map_summary.txt
echo -e "Simple_repeat_Reverse\t${TOTAL_READS}\t${ALL_Simple_repeat_READS_RVS}\t${UNIQ_Simple_repeat_READS_RVS}\t${MULTI_Simple_repeat_READS_RVS}" >> map_summary.txt
echo -e "Pachy\t${TOTAL_READS}\t${ALL_PACHY_READS}\t${UNIQ_PACHY_READS}\t${MULTI_PACHY_READS}" >> map_summary.txt
echo -e "Prepachy\t${TOTAL_READS}\t${ALL_PREPAC_READS}\t${UNIQ_PREPAC_READS}\t${MULTI_PREPAC_READS}" >> map_summary.txt
echo -e "Hybrid\t${TOTAL_READS}\t${ALL_HYBRID_READS}\t${UNIQ_HYBRID_READS}\t${MULTI_HYBRID_READS}" >> map_summary.txt
echo -e "LINE_Reverse_To_Pachy\t${TOTAL_READS}\t${ALL_LINE_PACHY_RVS_READS}\t${UNIQ_LINE_PACHY_RVS_READS}\t${MULTI_LINE_PACHY_RVS_READS}" >> map_summary.txt
echo -e "SINE_Reverse_To_Pachy\t${TOTAL_READS}\t${ALL_SINE_PACHY_RVS_READS}\t${UNIQ_SINE_PACHY_RVS_READS}\t${MULTI_SINE_PACHY_RVS_READS}" >> map_summary.txt
echo -e "LTR_Reverse_To_Pachy\t${TOTAL_READS}\t${ALL_LTR_PACHY_RVS_READS}\t${UNIQ_LTR_PACHY_RVS_READS}\t${MULTI_LTR_PACHY_RVS_READS}" >> map_summary.txt
echo -e "LINE_To_Pachy\t${TOTAL_READS}\t${ALL_LINE_PACHY_READS}\t${UNIQ_LINE_PACHY_READS}\t${MULTI_LINE_PACHY_READS}" >> map_summary.txt
echo -e "SINE_To_Pachy\t${TOTAL_READS}\t${ALL_SINE_PACHY_READS}\t${UNIQ_SINE_PACHY_READS}\t${MULTI_SINE_PACHY_READS}" >> map_summary.txt
echo -e "LTR_To_Pachy\t${TOTAL_READS}\t${ALL_LTR_PACHY_READS}\t${UNIQ_LTR_PACHY_READS}\t${MULTI_LTR_PACHY_READS}" >> map_summary.txt
echo -e "LINE_Reverse_To_Pachy_Reverse\t${TOTAL_READS}\t${ALL_LINE_PACHY_RVS_READS_RVS}\t${UNIQ_LINE_PACHY_RVS_READS_RVS}\t${MULTI_LINE_PACHY_RVS_READS_RVS}" >> map_summary.txt
echo -e "SINE_Reverse_To_Pachy_Reverse\t${TOTAL_READS}\t${ALL_SINE_PACHY_RVS_READS_RVS}\t${UNIQ_SINE_PACHY_RVS_READS_RVS}\t${MULTI_SINE_PACHY_RVS_READS_RVS}" >> map_summary.txt
echo -e "LTR_Reverse_To_Pachy_Reverse\t${TOTAL_READS}\t${ALL_LTR_PACHY_RVS_READS_RVS}\t${UNIQ_LTR_PACHY_RVS_READS_RVS}\t${MULTI_LTR_PACHY_RVS_READS_RVS}" >> map_summary.txt
echo -e "LINE_To_Pachy_Reverse\t${TOTAL_READS}\t${ALL_LINE_PACHY_READS_RVS}\t${UNIQ_LINE_PACHY_READS_RVS}\t${MULTI_LINE_PACHY_READS_RVS}" >> map_summary.txt
echo -e "SINE_To_Pachy_Reverse\t${TOTAL_READS}\t${ALL_SINE_PACHY_READS_RVS}\t${UNIQ_SINE_PACHY_READS_RVS}\t${MULTI_SINE_PACHY_READS_RVS}" >> map_summary.txt
echo -e "LTR_To_Pachy_Reverse\t${TOTAL_READS}\t${ALL_LTR_PACHY_READS_RVS}\t${UNIQ_LTR_PACHY_READS_RVS}\t${MULTI_LTR_PACHY_READS_RVS}" >> map_summary.txt
echo -e "LINE_Reverse_To_Prepachy\t${TOTAL_READS}\t${ALL_LINE_PREPAC_RVS_READS}\t${UNIQ_LINE_PREPAC_RVS_READS}\t${MULTI_LINE_PREPAC_RVS_READS}" >> map_summary.txt
echo -e "SINE_Reverse_To_Prepachy\t${TOTAL_READS}\t${ALL_SINE_PREPAC_RVS_READS}\t${UNIQ_SINE_PREPAC_RVS_READS}\t${MULTI_SINE_PREPAC_RVS_READS}" >> map_summary.txt
echo -e "LTR_Reverse_To_Prepachy\t${TOTAL_READS}\t${ALL_LTR_PREPAC_RVS_READS}\t${UNIQ_LTR_PREPAC_RVS_READS}\t${MULTI_LTR_PREPAC_RVS_READS}" >> map_summary.txt
echo -e "LINE_To_Prepachy\t${TOTAL_READS}\t${ALL_LINE_PREPAC_READS}\t${UNIQ_LINE_PREPAC_READS}\t${MULTI_LINE_PREPAC_READS}" >> map_summary.txt
echo -e "SINE_To_Prepachy\t${TOTAL_READS}\t${ALL_SINE_PREPAC_READS}\t${UNIQ_SINE_PREPAC_READS}\t${MULTI_SINE_PREPAC_READS}" >> map_summary.txt
echo -e "LTR_To_Prepachy\t${TOTAL_READS}\t${ALL_LTR_PREPAC_READS}\t${UNIQ_LTR_PREPAC_READS}\t${MULTI_LTR_PREPAC_READS}" >> map_summary.txt
echo -e "LINE_Reverse_To_Prepachy_Reverse\t${TOTAL_READS}\t${ALL_LINE_PREPAC_RVS_READS_RVS}\t${UNIQ_LINE_PREPAC_RVS_READS_RVS}\t${MULTI_LINE_PREPAC_RVS_READS_RVS}" >> map_summary.txt
echo -e "SINE_Reverse_To_Prepachy_Reverse\t${TOTAL_READS}\t${ALL_SINE_PREPAC_RVS_READS_RVS}\t${UNIQ_SINE_PREPAC_RVS_READS_RVS}\t${MULTI_SINE_PREPAC_RVS_READS_RVS}" >> map_summary.txt
echo -e "LTR_Reverse_To_Prepachy_Reverse\t${TOTAL_READS}\t${ALL_LTR_PREPAC_RVS_READS_RVS}\t${UNIQ_LTR_PREPAC_RVS_READS_RVS}\t${MULTI_LTR_PREPAC_RVS_READS_RVS}" >> map_summary.txt
echo -e "LINE_To_Prepachy_Reverse\t${TOTAL_READS}\t${ALL_LINE_PREPAC_READS_RVS}\t${UNIQ_LINE_PREPAC_READS_RVS}\t${MULTI_LINE_PREPAC_READS_RVS}" >> map_summary.txt
echo -e "SINE_To_Prepachy_Reverse\t${TOTAL_READS}\t${ALL_SINE_PREPAC_READS_RVS}\t${UNIQ_SINE_PREPAC_READS_RVS}\t${MULTI_SINE_PREPAC_READS_RVS}" >> map_summary.txt
echo -e "LTR_To_Prepachy_Reverse\t${TOTAL_READS}\t${ALL_LTR_PREPAC_READS_RVS}\t${UNIQ_LTR_PREPAC_READS_RVS}\t${MULTI_LTR_PREPAC_READS_RVS}" >> map_summary.txt
echo -e "LINE_Reverse_To_Hybrid\t${TOTAL_READS}\t${ALL_LINE_HYBRID_RVS_READS}\t${UNIQ_LINE_HYBRID_RVS_READS}\t${MULTI_LINE_HYBRID_RVS_READS}" >> map_summary.txt
echo -e "SINE_Reverse_To_Hybrid\t${TOTAL_READS}\t${ALL_SINE_HYBRID_RVS_READS}\t${UNIQ_SINE_HYBRID_RVS_READS}\t${MULTI_SINE_HYBRID_RVS_READS}" >> map_summary.txt
echo -e "LTR_Reverse_To_Hybrid\t${TOTAL_READS}\t${ALL_LTR_HYBRID_RVS_READS}\t${UNIQ_LTR_HYBRID_RVS_READS}\t${MULTI_LTR_HYBRID_RVS_READS}" >> map_summary.txt
echo -e "LINE_To_Hybrid\t${TOTAL_READS}\t${ALL_LINE_HYBRID_READS}\t${UNIQ_LINE_HYBRID_READS}\t${MULTI_LINE_HYBRID_READS}" >> map_summary.txt
echo -e "SINE_To_Hybrid\t${TOTAL_READS}\t${ALL_SINE_HYBRID_READS}\t${UNIQ_SINE_HYBRID_READS}\t${MULTI_SINE_HYBRID_READS}" >> map_summary.txt
echo -e "LTR_To_Hybrid\t${TOTAL_READS}\t${ALL_LTR_HYBRID_READS}\t${UNIQ_LTR_HYBRID_READS}\t${MULTI_LTR_HYBRID_READS}" >> map_summary.txt
echo -e "LINE_Reverse_To_Hybrid_Reverse\t${TOTAL_READS}\t${ALL_LINE_HYBRID_RVS_READS_RVS}\t${UNIQ_LINE_HYBRID_RVS_READS_RVS}\t${MULTI_LINE_HYBRID_RVS_READS_RVS}" >> map_summary.txt
echo -e "SINE_Reverse_To_Hybrid_Reverse\t${TOTAL_READS}\t${ALL_SINE_HYBRID_RVS_READS_RVS}\t${UNIQ_SINE_HYBRID_RVS_READS_RVS}\t${MULTI_SINE_HYBRID_RVS_READS_RVS}" >> map_summary.txt
echo -e "LTR_Reverse_To_Hybrid_Reverse\t${TOTAL_READS}\t${ALL_LTR_HYBRID_RVS_READS_RVS}\t${UNIQ_LTR_HYBRID_RVS_READS_RVS}\t${MULTI_LTR_HYBRID_RVS_READS_RVS}" >> map_summary.txt
echo -e "LINE_To_Hybrid_Reverse\t${TOTAL_READS}\t${ALL_LINE_HYBRID_READS_RVS}\t${UNIQ_LINE_HYBRID_READS_RVS}\t${MULTI_LINE_HYBRID_READS_RVS}" >> map_summary.txt
echo -e "SINE_To_Hybrid_Reverse\t${TOTAL_READS}\t${ALL_SINE_HYBRID_READS_RVS}\t${UNIQ_SINE_HYBRID_READS_RVS}\t${MULTI_SINE_HYBRID_READS_RVS}" >> map_summary.txt
echo -e "LTR_To_Hybrid_Reverse\t${TOTAL_READS}\t${ALL_LTR_HYBRID_READS_RVS}\t${UNIQ_LTR_HYBRID_READS_RVS}\t${MULTI_LTR_HYBRID_READS_RVS}" >> map_summary.txt

###plot map summary and distribution
BB_Pip_SRNA_Rscript_Mouse.R map_summary.txt ${PREFIX}_length.txt ${PREFIX}_strand.txt \
${PREFIX}_uniq_length.txt ${PREFIX}_uniq_strand.txt \
${PREFIX}_multi_length.txt ${PREFIX}_multi_strand.txt




