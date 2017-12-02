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
	echo "	bash BB_DNase_pipeline.sh <option>* [-q dnase.fq(.gz)] [-g genome]"
	echo ""
	echo "   Optional arguments:"
	echo "	-s min score for bam file filter. --default: 30"
	echo "	-p prefix for output file. --default: sample name"
	echo "	-o output directory. --default: ./"
	echo "	-q q-value for hotspot2 peak calling. --default: 0.05"
	echo "	-c CPU number used for pipeline. --default: 8"
	echo "  -G genome sequence file. --default: ~/Tianxiong/Annotation/Fasta/genome.fa"
	echo "  -i bwa index. --default: ~/Tianxiong/Annotation/Index/bwa/genome"
	echo "  -C chrom.size file. --default: ~/Tianxiong/Annotation/ChromSize/genome.chrom.size"
	echo "  -H hotspot source file's root name. --default: ~/Tianxiong/Annotation/Hotspot_CenterSites/genome"
	echo ""
	echo -e "\033[31m  !!! Caution: \033[0m"
	echo -e "	The program only can deal with single-end ff-firststrand library. If not, please modify the program."
	echo "	Program needed: hotspot2 bwa bedops edwstates picards samtools bedtools"
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

while getopts "hvq:c:o:g:p:w:s:R:S:G:d:l:r:" OPTION; do
	case $OPTION in
		h)	help_info && exit 0 ;;
		q)	DNASE_FQ=${OPTARG} ;;
		o)	OUTPATH=${OPTARG} ;;
		c)	CPU=${OPTARG} ;;
		C)	CHROMSIZE=${OPTARG};;
		g)	GENOME=${OPTARG};;
		p)	PREFIX=${OPTARG};;
		s)	SCORE=${OPTARG};;
		G)	GENOME_FA=${OPTARG};;
		q)	QVALUE=${OPTARG};;
		i)	INDEX_BWA=${OPTARG};;
		H)	HOTSPOT_SOURCE=${OPTARG};;
		*)	usage && exit 1 ;;
	esac
done

if [ ! -f "$DNASE_FQ" ];then
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
	PREFIX=${DNASE_FQ%.f*q*}
	PREFIX=${PREFIX##*/}
fi

if [ ! -n "$INDEX_BWA" ];then
	INDEX_BWA=/mnt/Storage/home/lijingyi/Tianxiong/Annotation/Index/bwa/${GENOME}
fi

if [ ! -n "$CHROMSIZE" ];then
	CHROMSIZE=/mnt/Storage/home/lijingyi/Tianxiong/Annotation/ChromSize/${GENOME}.chrom.size
fi

if [ ! -n "$HOTSPOT_SOURCE" ];then
	HOTSPOT_SOURCE=/mnt/Storage/home/lijingyi/Tianxiong/Annotation/Hotspot_CenterSites/${GENOME}
fi

if [ ! -n "$GENOME_FA" ];then
	GENOME_FA=/mnt/Storage/home/lijingyi/Tianxiong/Annotation/Fasta/${GENOME}.fa
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

if [ ! -d hotspots ];then
	mkdir hotspots
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

###check hotspot mappable file
if [ ! -f ${HOTSPOT_SOURCE}".K50.mappable_only_rmblack.bed" ];then
	echo -e "\033[40;31;1mERROR: no mappable file for centersite, please download one\033[0m"
	exit 1
fi

###check hotspot centersite file
if [ -f ${HOTSPOT_SOURCE}".center_sites.starch" ];then
	echo -e "\033[40;35;1mHotspot CenterSite: "${HOTSPOT_SOURCE}".center_sites.starch\033[0m"
else
	# no centersite file, make centersite via extractCenterSite.sh
	if [ ! -f ${HOTSPOT_SOURCE}".chrom_size.bed" ];then
		echo -e "\033[40;33;1mWARNING: no chrom_size file, create one\033[0m"
		if [ -f ${GENOME_FA} ];then
			DATE=`date --date="-24 hour"`
			echo -e "\033[40;32mmake chrom_size file\t"$DATE"\033[0m"
			set -x
			faSize -detailed ${GENOME_FA} | awk 'BEGIN{OFS="\t"} {print $1,0,$2}' | sort-bed - \
				> ${HOTSPOT_SOURCE}.chrom_size.bed
			set +x
		else
			echo -e "\033[40;31;1mERROR: no chrom_size file and no genome.fa to build it\033[0m"
			exit 1
		fi
	fi
	set -x
	extractCenterSites.sh -c ${HOTSPOT_SOURCE}".chrom_size.bed" -M \
		${HOTSPOT_SOURCE}".K50.mappable_only_rmblack.bed" -o \
		${HOTSPOT_SOURCE}".center_sites.starch" > \
		log_file/${PREFIX}_centersites.log 2>&1
	set +x
	echo -e "\033[40;35;1mHotspot CenterSite: "${HOTSPOT_SOURCE}".center_sites.starch\033[0m"
fi

###check chrom.size file
if [ -f ${CHROMSIZE} ];then
	echo -e "\033[40;35;1mChromSize: "${CHROMSIZE}"\033[0m"
else
	echo -e "\033[40;31;1mERROR: no chrom.size file, please download it use fetchChromSizes\033[0m"
fi


####################
###Align with bwa###
####################

###align
DATE=`date --date="-24 hour"`
echo -e "\033[40;36;1m\033[1m---Align......---\t"$DATE"\033[0m"
DATE=`date --date="-24 hour"`
echo -e "\033[40;32mbwa alignment......\t"$DATE"\033[0m"
set -x
bwa aln -Y -l 32 -n 0.04 -k 2 -t ${CPU} ${INDEX_BWA} ${DNASE_FQ} > tmp_${PREFIX}.sai \
	2>log_file/${PREFIX}_bwa_aln.log
bwa samse -n 10 ${INDEX_BWA} tmp_${PREFIX}.sai ${DNASE_FQ} > tmp_${PREFIX}.sam \
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
	READ_NAME_REGEX='[a-zA-Z0-9]+:[0-9]+:[a-zA-Z0-9]+:[0-9]+:([0-9]+):([0-9]+):([0-9]+).*'
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


##################
###run hotspot2###
##################


###run hotspot2
DATE=`date --date="-24 hour"`
echo -e "\033[40;32mrun hotspot2......\t"$DATE"\033[0m"
set -x
hotspot2.sh -c ${HOTSPOT_SOURCE}".chrom_size.bed" -C ${HOTSPOT_SOURCE}".center_sites.starch" \
	-M ${HOTSPOT_SOURCE}.K50.mappable_only_rmblack.bed alignments/${PREFIX}_mkdup_filter.bam \
	${PREFIX}_out/ > log_file/${PREFIX}_hotspot2.log 2>&1
set +x

###convert hotspot peaks to bed and bigBed
DATE=`date --date="-24 hour"`
echo -e "\033[40;32mconvert hotspot to bed and bigBed......\t"$DATE"\033[0m"
unstarch ${PREFIX}_out/*.peaks.narrowpeaks.starch > peaks/${PREFIX}.bed
bedToBigBed -as=/mnt/Storage/home/lijingyi/Tianxiong/bin/narrowPeak.as -type=bed6+4 \
	peaks/${PREFIX}.bed ${CHROMSIZE} peaks/${PREFIX}.bb

###saving other files
mv ${PREFIX}_out/*.density.bw bigWig/${PREFIX}.bw
mv ${PREFIX}_out/*.density.starch bigWig/${PREFIX}.starch
mv ${PREFIX}_out/*.SPOT.txt hotspots/${PREFIX}_SPOT.txt
unstarch ${PREFIX}_out/*.hotspot*.starch > hotspots/${PREFIX}.bed 
mv ${PREFIX}_out/*.allcalls.starch hotspots/${PREFIX}.allcalls.starch

###remove useless files
rm -rf ${PREFIX}_out


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
NUM_HOTSPOT=(`grep "^chr" hotspots/${PREFIX}.bed | wc -l`)
RATIO_SPOT=(`cat hotspots/${PREFIX}_SPOT.txt`)
echo -e "\033[40;32mcalculate read number in chrX and chrY......\t"$DATE"\033[0m"
NUM_READ_CHR_X=(`samtools view alignments/${PREFIX}_mkdup_filter.bam | awk '$3=="chrX"' | wc -l`)
NUM_READ_CHR_Y=(`samtools view alignments/${PREFIX}_mkdup_filter.bam | awk '$3=="chrY"' | wc -l`)

###write information to basic_info.txt
if [ ! -f "summary/basic_info.txt" ];then
	echo -e "\tread_size_mean\tread_num\tmapped_read_num\tuniq_mapped_read_num\tdup_read_num\thigh_qual_read_num\thigh_qual_dup_read_num\tpeak_num\thotspot_num\thotspot_read_num\tchrX_read_num\tchrY_read_num" > summary/basic_info.txt
fi  # initialize basic_info.txt
echo -e "${PREFIX}\t${SIZE_READ}\t${NUM_READ}\t${NUM_MAPPED_READ}\t${NUM_UNIQ_MAPPED_READ}\t${NUM_DUP_READ}\t${NUM_READ_HIGH_QUAL}\t${NUM_DUP_READ_HIGH_QUAL}\t${NUM_PEAK}\t${NUM_HOTSPOT}\t${RATIO_SPOT}\t${NUM_READ_CHR_X}\t${NUM_READ_CHR_Y}" >> summary/basic_info.txt

###finished
echo -e "\033[40;32;1mfinished\t"$DATE"\033[0m"

