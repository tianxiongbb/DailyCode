#!/bin/bash

##----INTRO-----------##
# Name=BB_DNase_pipeline.sh
# Date=Jan13 ,2017
# Update=Jan13, 2017
# Update information:

########################
# Purpose
# DNase pipeline like ENCODE

#######--Arguments--#######
help_info(){
	echo ""
	echo -e "\033[32m  =========================================================================================================================\033[0m"
	echo "   Author: Tianxiong"
	echo ""
	echo "   Usage:"
	echo "	bash BB_chipseq_pipeline.sh <option>* [-q chip.fq(.gz)] [-g genome]"
	echo ""
	echo "   Optional arguments:"
	echo "  -I input file for macs2 peak calling, can be .fastq or .bam. --default: no input file"
	echo "	-s min score for bam file filter. --default: 30"
	echo "	-p prefix for output file. --default: sample name"
	echo "	-o output directory. --default: ./"
	echo "	-q q-value for macs2 peak calling. --default: 0.05"
	echo "  -S size of genome. human: 2.7e9, mouse: 1.87e9"
	echo "	-c CPU number used for pipeline. --default: 8"
	echo "  -t peak type for macs2 peak calling. can be narrow/broad. --default: narrow"
	echo "  -G genome sequence file. --default: ~/Tianxiong/Annotation/Fasta/genome.fa"
	echo "  -i bwa index. --default: ~/Tianxiong/Annotation/Index/bwa/genome"
	echo "  -C chrom.size file. --default: ~/Tianxiong/Annotation/ChromSize/genome.chrom.size"
	echo ""
	echo -e "\033[31m  !!! Caution: \033[0m"
	echo -e "	The program only can deal with single-end ff-firststrand library. If not, please modify the program."
	echo "	Program needed: macs2 bwa bedops edwstates picards samtools bedtools"
	echo "	Make sure all the program is contained by environment viariable PATH."
	echo "	The pipeline also need index files builded already, if not the pipeline will build it."
	echo ""
	echo "    (￣(工)￣) Enjoy yourself~~~"
	echo -e "\033[32m  =========================================================================================================================\033[0m"
	echo ""
}

if [ $# -lt 2 ];then
	help_info
	exit 1
fi

#############################
# functions #
#############################
function fun_sort_bed()
{
	sort -k1,1 -k2,2n $1 > t_fun_sort && mv t_fun_sort $1
} # sort bed file

#############################
# ARGS reading and checking #
#############################
OUTPATH=./
PREFIX=0
CPU=8
SCORE=30
QVALUE=0.05
PEAK_TYPE=narrow

while getopts "hvq:c:o:g:p:w:s:R:G:d:l:r:t:I:S:" OPTION; do
	case $OPTION in
		h)	help_info && exit 0 ;;
		q)	CHIP_FQ=${OPTARG} ;;
		o)	OUTPATH=${OPTARG} ;;
		c)	CPU=${OPTARG} ;;
		C)	CHROMSIZE=${OPTARG};;
		g)	GENOME=${OPTARG};;
		t)	PEAK_TYPE=${OPTARG};;
		p)	PREFIX=${OPTARG};;
		s)	SCORE=${OPTARG};;
		S)	SIZE_GENOME=${OPTARG};;
		G)	GENOME_FA=${OPTARG};;
		q)	QVALUE=${OPTARG};;
		i)	INDEX_BWA=${OPTARG};;
		I)	INPUT=${OPTARG};;
		*)	usage && exit 1 ;;
	esac
done

if [ ! -f "$CHIP_FQ" ];then
	echo -e "\033[40;31;1mERROR: please use -q to specify the right DNase-seq file\033[0m"
	exit 1
fi

if [ ! -d "$OUTPATH" ];then
	echo -e "\033[40;33;1mWARNING: output path not found, creat one\033[0m"
	mkdir -p ${OUTPATH}
fi

if [ ! -n "$GENOME" ];then
	echo -e "\033[40;31;1mERROR: please use -g to specify genome used, eg. mm10\033[0m"
	exit 1
fi

if [ "$PREFIX" == "0" ];then
	echo -e "\033[40;33;1mWARNING: no prefix name, use sample name as prefix\033[0m"
	PREFIX=${CHIP_FQ%.f*q*}
	PREFIX=${PREFIX##*/}
fi

if [ ! -n "$INDEX_BWA" ];then
	INDEX_BWA=/data/tusers/yutianx/tongji2/Annotation/Index/bwa/${GENOME}
fi

if [ ! -n "$CHROMSIZE" ];then
	CHROMSIZE=/data/tusers/yutianx/tongji2/Annotation/ChromSize/${GENOME}.chrom.size
fi

if [ ! -n "$GENOME_FA" ];then
	GENOME_FA=/data/tusers/yutianx/tongji2/Annotation/Fasta/${GENOME}.fa
fi

if [ -n "${INPUT}" ];then
	if [ ! -f "${INPUT}" ];then
		echo -e "\033[40;31;1mERROR: please use -I to specify the right Input file, or do not use -I\033[0m"
		exit 1
	fi
fi

if [ -n $SIZE_GENOME ];then
	if [ "${GENOME:0:2}" == "mm" ];then
		SIZE_GENOME=1.87e9
	elif [ "${GENOME:0:2}" == "hg" ];then
		SIZE_GENOME=2.7e9
	else
		echo -e "\033[40;31;1mERROR: please use -S to specify the genome size\033[0m"
		exit 1
	fi
fi

###########
# process #
###########


#################
###Preparation###
#################


###make directories
DATE=`date --date="-24 hour"`
echo -e "\033[40;36;1m\033[1m---Preparation......\t"$DATE"\033[0m"
cd ${OUTPATH}
if [ ! -d log_file ];then
	mkdir log_file
fi

if [ ! -d alignments ];then
	mkdir alignments
fi

if [ ! -d figures ];then
	mkdir figures
fi

if [ ! -d bigWig ];then
	mkdir bigWig
fi

if [ ! -d summary ];then
	mkdir summary
fi

if [ ! -d peaks ];then
	mkdir peaks
fi

if [ ! -d macs2 ];then
	mkdir macs2
fi

if [ ! -d bamstat ];then
	mkdir bamstat
fi

###check index
if [ -f ${INDEX_BWA}".amb" -a -f ${INDEX_BWA}".ann" -a -f ${INDEX_BWA}".bwt" -a -f ${INDEX_BWA}".pac" -a -f ${INDEX_BWA}".sa" ];then
	echo -e "\033[40;35;1mIndex: "${INDEX_BWA}"\033[0m"
else
	# no index file, make index
	if [ -f ${GENOME_FA} ];then
		DATE=`date --date="-24 hour"`
		echo -e "\033[40;32mmake index......\t"$DATE"\033[0m"
		set -x
		bwa index -p ${INDEX_BWA} -a bwtsw ${GENOME_FA}
		set +x
	else
		echo -e "\033[40;31;1mERROR: no index file and no genome.fa to build it\033[0m"
		exit 1
	fi
	echo -e "\033[40;35;1mIndex: "${INDEX_BWA}"\033[0m"
fi

###check chrom.size file
if [ -f ${CHROMSIZE} ];then
	echo -e "\033[40;35;1mChromSize: "${CHROMSIZE}"\033[0m"
else
	echo -e "\033[40;31;1mERROR: no chrom.size file, please download it use fetchChromSizes\033[0m"
fi

###check input file
echo -e "\033[40;35;1mInputFile: "${INPUT}"\033[0m"

####################
###Align with bwa###
####################

###align
DATE=`date --date="-24 hour"`
echo -e "\033[40;36;1m\033[1m---Align......---\t"$DATE"\033[0m"
DATE=`date --date="-24 hour"`
echo -e "\033[40;32mbwa alignment......\t"$DATE"\033[0m"
set -x
bwa aln -Y -l 32 -n 0.04 -k 2 -t ${CPU} ${INDEX_BWA} ${CHIP_FQ} > tmp_${PREFIX}.sai \
	2>log_file/${PREFIX}_bwa_aln.log
bwa samse -n 10 ${INDEX_BWA} tmp_${PREFIX}.sai ${CHIP_FQ} > tmp_${PREFIX}.sam \
	2>log_file/${PREFIX}_bwa_samse.log
samtools view -bhS tmp_${PREFIX}.sam > tmp_${PREFIX}.bam
samtools sort -@ ${CPU} -m 4G -f tmp_${PREFIX}.bam alignments/${PREFIX}.bam
rm tmp_${PREFIX}.bam tmp_${PREFIX}.sam tmp_${PREFIX}.sai
samtools index alignments/${PREFIX}.bam
set +x

###collect bam stats
DATE=`date --date="-24 hour"`
echo -e "\033[40;32mcollect bam stats......\t"$DATE"\033[0m"
set -x
samtools flagstat alignments/${PREFIX}.bam > bamstat/${PREFIX}_flagstat.txt
edwBamStats alignments/${PREFIX}.bam bamstat/${PREFIX}_edwBamStats.txt
set +x


#######################
###Filter alignments###
#######################


###mark duplicates with picard
DATE=`date --date="-24 hour"`
echo -e "\033[40;32mrun picard mark duplicates on non-UMI......\t"$DATE"\033[0m"
set -x
time java -Xmx4G -jar /mnt/Storage/home/lijingyi/Tianxiong/bin/MarkDuplicates.jar \
	INPUT=alignments/${PREFIX}.bam OUTPUT=${PREFIX}_marked.bam \
	METRICS_FILE=${PREFIX}_dup_qc.txt ASSUME_SORTED=true \
	VALIDATION_STRINGENCY=SILENT \
	READ_NAME_REGEX='[a-zA-Z0-9]+:[0-9]+:[a-zA-Z0-9]+:[0-9]+:([0-9]+):([0-9]+):([0-9]+).*' \
	> log_file/${PREFIX}_picard.log 2>&1
set +x
mv ${PREFIX}_dup_qc.txt bamstat

###fliter bam on flags and threashold
DATE=`date --date="-24 hour"`
echo -e "\033[40;32mfilter on flag 512 and threadshold score: "$SCORE"......\t"$DATE"\033[0m"
set -x
samtools view -F 512 -q ${SCORE} -b ${PREFIX}_marked.bam > alignments/${PREFIX}_mkdup_filter.bam
set +x
rm ${PREFIX}_marked.bam 

###collect bam stats
DATE=`date --date="-24 hour"`
echo -e "\033[40;32mcollect bam stats......\t"$DATE"\033[0m"
set -x
samtools flagstat alignments/${PREFIX}_mkdup_filter.bam > bamstat/${PREFIX}_mkdup_filter_flagstat.txt
samtools stats alignments/${PREFIX}_mkdup_filter.bam > bamstat/${PREFIX}_mkdup_filter_samstat.txt
grep ^SN bamstat/${PREFIX}_mkdup_filter_samstat.txt | cut -f 2- > bamstat/${PREFIX}_mkdup_filter_samstat_summary.txt
set +x

#################################
###repeat if input.fastq given###
#################################


TAIL_INPUT=${INPUT#*.}
if [ "$TAIL_INPUT" == "fastq" ];then
	DATE=`date --date="-24 hour"`
	echo -e "\033[40;36;1m\033[1m---repeat for input file......---\t"$DATE"\033[0m"
	set -x
	bwa aln -Y -l 32 -n 0.04 -k 2 -t ${CPU} ${INDEX_BWA} ${INPUT} > tmp_${PREFIX}_input.sai \
		2>log_file/${PREFIX}_bwa_aln.log
	bwa samse -n 10 ${INDEX_BWA} tmp_${PREFIX}_input.sai ${INPUT} > tmp_${PREFIX}_input.sam \
		2>log_file/${PREFIX}_input_bwa_samse.log
	samtools view -bhS tmp_${PREFIX}_input.sam > tmp_${PREFIX}_input.bam
	samtools sort -@ ${CPU} -m 4G -f tmp_${PREFIX}_input.bam alignments/${PREFIX}_input.bam
	rm tmp_${PREFIX}_input.bam tmp_${PREFIX}_input.sam tmp_${PREFIX}_input.sai
	samtools index alignments/${PREFIX}_input.bam
	set +x
	set -x
	samtools flagstat alignments/${PREFIX}_input.bam > bamstat/${PREFIX}_input_flagstat.txt
	edwBamStats alignments/${PREFIX}_input.bam bamstat/${PREFIX}_input_edwBamStats.txt
	set +x
	set -x
	time java -Xmx4G -jar /mnt/Storage/home/lijingyi/Tianxiong/bin/MarkDuplicates.jar \
		INPUT=alignments/${PREFIX}_input.bam OUTPUT=${PREFIX}_input_marked.bam \
		METRICS_FILE=${PREFIX}_input_dup_qc.txt ASSUME_SORTED=true \
		VALIDATION_STRINGENCY=SILENT \
		READ_NAME_REGEX='[a-zA-Z0-9]+:[0-9]+:[a-zA-Z0-9]+:[0-9]+:([0-9]+):([0-9]+):([0-9]+).*' \
		> log_file/${PREFIX}_picard.log 2>&1
	set +x
	mv ${PREFIX}_input_dup_qc.txt bamstat
	set -x
	samtools view -F 512 -q ${SCORE} -b ${PREFIX}_input_marked.bam > alignments/${PREFIX}_input_mkdup_filter.bam
	set +x
	rm ${PREFIX}_input_marked.bam 
	set -x
	samtools flagstat alignments/${PREFIX}_input_mkdup_filter.bam > bamstat/${PREFIX}_input_mkdup_filter_flagstat.txt
	samtools stats alignments/${PREFIX}_input_mkdup_filter.bam > bamstat/${PREFIX}_input_mkdup_filter_samstat.txt
	grep ^SN bamstat/${PREFIX}_input_mkdup_filter_samstat.txt | cut -f 2- > bamstat/${PREFIX}_input_mkdup_filter_samstat_summary.txt
	set +x
	INPUT=alignments/${PREFIX}_input_mkdup_filter.bam
fi

#######################
###evaluate bam file###
#######################


:<<!
###filter out chrM
DATE=`date --date="-24 hour"`
echo -e "\033[40;32mfilter out chrM......\t"$DATE"\033[0m"
set -x
edwBamFilter -sponge -chrom=chrM alignments/${PREFIX}_mkdup_filter.bam alignments/${PREFIX}_mkdup_filter_no_chrM.bam
samtools index alignments/${PREFIX}_mkdup_filter_no_chrM.bam
set +x
!


###############
###run macs2###
###############


###run masc2
DATE=`date --date="-24 hour"`
echo -e "\033[40;32mrun macs2......\t"$DATE"\033[0m"
if [ -n $INPUT ];then
	if [ "$PEAK_TYPE" == "broad" ];then
		set -x
		macs2 callpeak -t alignments/${PREFIX}_mkdup_filter.bam -c ${INPUT} \
			-g ${SIZE_GENOME} --keep-dup 1 --outdir macs2 -n ${PREFIX} --SPMR \
			-q ${QVALUE} -B --broad > log_file/${PREFIX}_macs2.log 2>&1
		set +x
	else
		set -x
		macs2 callpeak -t alignments/${PREFIX}_mkdup_filter.bam -c ${INPUT} \
			-g ${SIZE_GENOME} --keep-dup 1 --outdir macs2 -n ${PREFIX} --SPMR \
			-q ${QVALUE} -B > log_file/${PREFIX}_macs2.log 2>&1
		set +x
	fi
else
	if [ "$PEAK_TYPE" == "broad" ];then
		set -x
		macs2 callpeak -t alignments/${PREFIX}_mkdup_filter.bam \
			-g ${SIZE_GENOME} --keep-dup 1 --outdir macs2 -n ${PREFIX} --SPMR \
			-q ${QVALUE} -B --broad > log_file/${PREFIX}_macs2.log 2>&1
		set +x
	else
		set -x
		macs2 callpeak -t alignments/${PREFIX}_mkdup_filter.bam \
			-g ${SIZE_GENOME} --keep-dup 1 --outdir macs2 -n ${PREFIX} --SPMR \
			-q ${QVALUE} -B > log_file/${PREFIX}_macs2.log 2>&1
		set +x
	fi
fi

###convert bedGraph to bigWig
bedGraphToBigWig macs2/${PREFIX}_treat_pileup.bdg ${CHROMSIZE} bigWig/${PREFIX}.bw
if [ "$PEAK_TYPE" == "broad" ];then
	cp macs2/${PREFIX}_peaks.broadPeak peaks/${PREFIX}.bed
else
	cp macs2/${PREFIX}_peaks.narrowPeak peaks/${PREFIX}.bed
fi


###################
###write summary###
###################


###calculate basic information 
echo -e "\033[40;32mcalculate basic information......\t"$DATE"\033[0m"
SIZE_READ=`awk '{if($1=="readSizeMean") {print $2}}' bamstat/${PREFIX}_edwBamStats.txt`
NUM_READ=`awk '{if($1=="readCount") {print $2}}' bamstat/${PREFIX}_edwBamStats.txt`
NUM_MAPPED_READ=`awk '{if($1=="mappedCount") {print $2}}' bamstat/${PREFIX}_edwBamStats.txt`
NUM_UNIQ_MAPPED_READ=`awk '{if($1=="uniqueMappedCount") {print $2}}' bamstat/${PREFIX}_edwBamStats.txt`
NUM_DUP_READ=`awk 'BEGIN{FS="\t"} {if($1 ~ /Library/) {print $5}}' bamstat/${PREFIX}_dup_qc.txt`
NUM_READ_HIGH_QUAL=`awk 'BEGIN{FS="\t"} {if($1=="reads mapped:") {print $2}}' bamstat/${PREFIX}_mkdup_filter_samstat_summary.txt`
NUM_DUP_READ_HIGH_QUAL=`awk 'BEGIN{FS="\t"} {if($1=="reads duplicated:") {print $2}}' bamstat/${PREFIX}_mkdup_filter_samstat_summary.txt`
NUM_PEAK=(`grep "^chr" peaks/${PREFIX}.bed | wc -l`)
NUM_PEAK_5F=(`grep "^chr" peaks/${PREFIX}.bed | awk '$7>=5' | wc -l`)
echo -e "\033[40;32mcalculate read number in chrX and chrY......\t"$DATE"\033[0m"
NUM_READ_CHR_X=(`samtools view alignments/${PREFIX}_mkdup_filter.bam | awk '$3=="chrX"' | wc -l`)
NUM_READ_CHR_Y=(`samtools view alignments/${PREFIX}_mkdup_filter.bam | awk '$3=="chrY"' | wc -l`)

###write information to basic_info.txt
if [ ! -f "summary/basic_info.txt" ];then
	echo -e "\tread_size_mean\tread_num\tmapped_read_num\tuniq_mapped_read_num\tdup_read_num\thigh_qual_read_num\thigh_qual_dup_read_num\tpeak_num\tpeak_num_5f\tchrX_read_num\tchrY_read_num" > summary/basic_info.txt
fi  # initialize basic_info.txt
echo -e "${PREFIX}\t${SIZE_READ}\t${NUM_READ}\t${NUM_MAPPED_READ}\t${NUM_UNIQ_MAPPED_READ}\t${NUM_DUP_READ}\t${NUM_READ_HIGH_QUAL}\t${NUM_DUP_READ_HIGH_QUAL}\t${NUM_PEAK}\t${NUM_PEAK_5F}\t${NUM_READ_CHR_X}\t${NUM_READ_CHR_Y}" >> summary/basic_info.txt

###finished
echo -e "\033[40;32;1mfinished\t"$DATE"\033[0m"

