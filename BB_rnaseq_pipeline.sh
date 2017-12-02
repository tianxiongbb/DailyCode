#!/bin/bash
#TODO: add gene signal plot

##----INTRO-----------##
# Name=BB_rnaseq_pipeline
VERSION=1.20
# Date=Apr20 ,2016
# Update=Jan 15, 2017
# Update information:
# Update=Oct 27, 2017
# Update information: 1. add CPU for htseq; 2. use BB_bamtobw to make bigWig files for visualizing. 3. add a choice for repbase mapping.

########################
# Purpose
#This file is for doing some pre-analysis for each RNA-seq samples
#1.Mapping: bowtie2 + STAR
#2.Count Feature Signal: htseq-count
#3.Plot Results: to do

#######--Arguments--#######
help_info(){
echo -e "\033[40;33;1m"
cat << EOF
Author: Tianxiong Yu / Bear
usage:
bash BB_rnaseq_pipeline <option>* (-l input_left.fq) (-g genome) [-r input_right.fq]"

Options:

    Input:
	-l left RNA-seq fastq file.
	-r right RNA-seq fastq file. [ set to N or not set for single-end ]
	-g genome used for this pipeline. eg: hg38 or mm10
    Output:
	-p prefix output name. --default: sample name
	-o output directory. --default: ./
    rRNA removing:
	--no-rRNA-removing do not remove rRNA. default: remove rRNA
	--index-rRNA bowtie2 index for removing rRNA.
    mapping:
	--index-genome STAR index for genome mapping.
	-M max allowed multiple mapped reads for STAR map. --default: 100 (-1 mean not allowed)
	-m max mismatch allowed for STAR map. --default: 2
    signal generating and visualization
	-C chrom.size file for bigWig file generate.
	-L library type [reverse|yes|no|miss], see htseq-count help file for more information; you can set to miss if you are not sure. --default: reverse
	--annotation-bed12 annotation file in bed12 format for guessing library type, used with -L miss.
	--no-duplication-removing don't remove duplicates. --default: removed
	-i attribute for htseq feature id. --default: gene_id
	--nonunique [none|all] whether to score reads that are not uniquely aligned or ambiguously assigned to features for htseq. default: none
	--htseq-mode [union|intersection-strict|intersection-nonempty] mode for htseq, see htseq-count document for more detail. default: intersection-nonempty
	-H gtf file for htseq count, htseq will use id for feature name; set to 0 or not set if do not calculate feature signal; if have several gtf for htseq count, use like this: a.gtf,b.gtf,c.gtf. default: 0
	--gene-signal-plot plot gene signal figures. --default: not plot
	--plot-file bed5 file for signal plot. --default: plot genes from gtf files specified by -H
	--transcript-bed7 bed7 file of transcripts for signal plot. --default: plot transcripts from gtf files specified by -H
    transposon analysis:
	--transposon-analysis do transposon analysis. default: do not analysis transposon expression
	--index-transposon index for transposon mapping.
	--size-transposon transposon size file like chromosome size file.
    performance:
	-c CPU used. --default: 8
EOF
echo -e "\033[0m"
}


if [ $# -lt 1 ];then
	help_info
	exit 1
fi

CPU=8
OUTPATH=./
PREFIX=0
LIBRARY=miss
RIGHT=0
MAPNUM=100
MISMATCH=2
ATTRIBUTE=gene_id
HTSEQ=0
NONUNIQUE=none
HTSEQ_MODE=intersection-nonempty
SWITCH_RMDUP=1
SWITCH_RRNA=1
SWITCH_TRANSPOSON=0
SWITCH_PLOT=0
PLOT_FILE=0
TRANSCRIPT_BED7=0

echo0 1 "Command:"
COMMAND="BB_rnaseq_pipeline "$*
echo0 3 "$COMMAND"
echo ""

ARGS=`getopt -a -o l:r:g:p:o:M:m:C:c:L:i:H: -l no-rRNA-removing,index-rRNA:,index-genome:,no-duplication-removing,nonunique:,htseq-mode:,transposon-analysis,index-transposon:,size-transposon:,annotation-bed12:,gene-signal-plot,plot-file:,transcript-bed7: -- "$@"`
[ $? -ne 0 ] && usage
eval set -- "${ARGS}"
while true
do
	case "$1" in
		-h)	help_info && exit 0;;
		-v)	echo0 1 "BB_rnaseq_pipeline"${VERSION} && exit 0;;
		-l)	LEFT=`readlink -f $2` && shift 2;;
		-r)	RIGHT=`readlink -f $2` && shift 2;;
		-g)	GENOME=$2 && shift 2;;
		-p)	PREFIX=$2 && shift 2;;
		-o)	OUTPATH=$2 && shift 2;;
		-M)	MAPNUM=$2 && shift 2;;
		-m)	MISMATCH=$2 && shift 2;;
		-C)	CHROMSIZE=$2 && shift 2;;
		-c)	CPU=$2 && shift 2;;
		-L)	LIBRARY=$2 && shift 2;;
		-i)	ATTRIBUTE=$2 && shift 2;;
		-H)	HTSEQ=$2 && shift 2;;
		--no-rRNA-removing)	SWITCH_RRNA=0 && shift;;
		--index-rRNA)	INDEX_RRNA=$2 && shift 2;;
		--index-genome)	INDEX_STAR=$2 && shift 2;;
		--no-duplication-removing)	SWITCH_RMDUP=0 && shift;;
		--nonunique)	NONUNIQUE=$2 && shift 2;;
		--htseq-mode)	HTSEQ_MODE=$2 && shift 2;;
		--transposon-analysis)	SWITCH_TRANSPOSON=1 && shift;;
		--index-transposon)	INDEX_TRANSPOSON=$2 && shift 2;;
		--size-transposon)	SIZE_TRANSPOSON=$2 && shift 2;;
		--annotation-bed12)	ANNOTATION_BED12=$2 && shift 2;;
		--gene-signal-plot)	SWITCH_PLOT=1 && shift 1;;
		--plot-file)	PLOT_FILE=$2 && shift 2;;
		--transcript-bed7)	TRANSCRIPT_BED7=$2 && shift 2;;
		--)	shift && break;;
		*)	help_info && exit 1;;
	esac
done

######Configure Tools########
function checkTools(){
	if [ `which $1` ];then
		echo0 3 `which $1`
	else
		echo0 0 $1" not found, please install it or add it into your PATH!"
		exit 1
	fi
}

echo0 1 "tools used:"
checkTools STAR
checkTools bowtie2
checkTools samtools
checkTools bedtools
checkTools htseq-count
checkTools ParaFly
checkTools R
checkTools BB_GtfToBed7_For_Sushi.sh
checkTools BB_plot_gene_dis.R
checkTools BB_plot_transposon_dis.R
checkTools BB_GtfToExonForEachGene.py
checkTools BB_GtfToBed7_For_Sushi.sh
if [ ${LIBRARY} == "miss" ];then
	checkTools infer_experiment.py
fi
checkTools ParaFly

######Configure Mode########
echo0 1 "mode and switchs:"
if [ ! -f "$RIGHT" ];then
	echo0 3 "paired-end or single-end: single-end"
else
	echo0 3 "paired-end or single-end: paired-end"
fi

if [ $SWITCH_RRNA -ne 0 ];then
	echo0 3 "remove rRNA: yes"
else
	echo0 3 "remove rRNA: no"
fi

if [ $SWITCH_RMDUP -ne 0 ];then
	echo0 3 "remove duplicates: yes"
else
	echo0 3 "remove duplicates: no"
fi

if [ "$HTSEQ" = "0" ];then
	echo0 3 "calculate signal for features: no"
else
	echo0 3 "calculate signal for features: yes"
fi

if [ $SWITCH_TRANSPOSON -ne 0 ];then
	echo0 3 "analyze transposons: yes"
else
	echo0 3 "analyze transposons: no"
fi

if [ $SWITCH_PLOT -ne 0 ];then
	echo0 3 "plot gene signal: yes"
else
	echo0 3 "plot gene signal: no"
fi

######Configure Parameters########
echo0 1 "configuring parameters......"
if [ ! -f "$LEFT" ];then
	echo0 0 "Error: please specify the correct input fastq file"
	exit 1
fi

if [ "$PREFIX" == "0" ];then
	PREFIX=`basename ${LEFT%.f*q*}`
	echo0 4 "WARNING: no prefix name. set "$PREFIX" as prefix name"
fi

if [ ! -n "$GENOME" ];then
	echo0 0 "Error: please use -g to specify genome used"
	exit 1
fi

if [ ! -d "$OUTPATH" ];then
	echo0 4 "WARNING: output path not found, create one"
	mkdir -p ${OUTPATH}
fi

if [ ! -n "$INDEX_STAR" ];then
	INDEX_STAR=/data/tusers/yutianx/tongji2/Annotation/Index/${GENOME}_star_gtf/
fi

if [ ! -n "$INDEX_RRNA" ];then
	INDEX_RRNA=/data/tusers/yutianx/tongji2/Annotation/Index/${GENOME}_bowtie2/rRNA
fi

if [ ! -n "$INDEX_TRANSPOSON" ];then
	INDEX_TRANSPOSON=/data/tusers/yutianx/tongji2/Annotation/Index/${GENOME}_bowtie2/repbase
fi

if [ ! -n "$CHROMSIZE" ];then
	CHROMSIZE=/data/tusers/yutianx/tongji2/Annotation/ChromSize/${GENOME}.chrom.size.long
fi
if [ ! -f ${CHROMSIZE} ];then
	echo0 0 "Error: no genome chrom.size data in "${CHROMSIZE}", please download via fetchChromSize"
	exit 1
fi

if [ ! -n "${SIZE_TRANSPOSON}" ];then
	SIZE_TRANSPOSON=/data/tusers/yutianx/tongji2/Annotation/ChromSize/${GENOME}.repbase.chrom.size
fi
if [ ! -f ${SIZE_TRANSPOSON} ];then
	echo0 0 "Error: no repbase chrom.size data in "${SIZE_TRANSPOSON}", please download via fetchChromSize"
	exit 1
fi

if [ "$LIBRARY" = "miss" ];then
	if [ ! -n "${ANNOTATION_BED12}" ];then
		ANNOTATION_BED12=/data/tusers/yutianx/tongji2/Annotation/GtfGff/${GENOME}.bed 
	fi
	if [ ! -f ${ANNOTATION_BED12} ];then
		echo0 0 "Error: library type set to miss but no annotation file in "${ANNOTATION_BED12}
		exit 1
	else
		NF=`head -1 ${ANNOTATION_BED12} | awk 'BEGIN{FS="\t"} {print NF}'`
		if [ "$NF" != "12" ]; then
			echo0 0 "Error: annotation file is not bed12 file"
			exit 1
		fi
	fi
fi

if [ ! "${PLOT_FILE}" = "0" ];then
	if [ ! -f ${PLOT_FILE} ];then
		echo0 0 "Error: plot file specified but not find, please input the right plot file"
		exit 1
	fi
fi

HTSEQ_COPY=$HTSEQ
if [ "$HTSEQ" != "0" ];then
	while true
	do
		if [ ! -f ${HTSEQ%%,*} ]; then
			echo0 0 "Error: no gtf file in "${HTSEQ%%,*}	
		fi
		last=$HTSEQ
		HTSEQ=${HTSEQ#*,}
		if [ "$HTSEQ" = "$last" ];then
			break
		fi
	done
fi


###########
# process #
###########


###############
# Preparation #
###############

###make directories
cd ${OUTPATH}

if [ ! -d log_file ];then
	mkdir log_file
fi

if [ ! -d bowtie2 ];then
	mkdir bowtie2
fi

if [ ! -d STAR ];then
	mkdir STAR
fi

if [ ! -d signal ];then
	mkdir signal
fi

if [ ! -d summary ];then
	mkdir summary
fi

if [ ! -d bigWig ];then
	mkdir bigWig
fi

if [ ! -d plot_input ];then
	mkdir plot_input
fi

###check index
echo0 1 "check indexes:"
GENOME_FA="/data/tusers/yutianx/tongji2/Annotation/Fasta/"$GENOME".fa"
if [ -f $INDEX_STAR"chrLength.txt" -a -f $INDEX_STAR"chrNameLength.txt" -a -f $INDEX_STAR"chrName.txt" -a -f $INDEX_STAR"chrStart.txt" -a -f $INDEX_STAR"Genome" -a -f $INDEX_STAR"genomeParameters.txt" -a -f $INDEX_STAR"SA" -a -f $INDEX_STAR"SAindex" ];then
	echo0 3 "genome Index: "$INDEX_STAR
else
	echo0 0 "Error: no STAR index in "${INDEX_STAR}", please build with STAR --runMode genomeGenerate"
	exit 1
fi

if [ ${SWITCH_RRNA} -ne 0 ];then
	if [ -f $INDEX_RRNA".1.bt2" -a -f $INDEX_RRNA".2.bt2" -a -f $INDEX_RRNA".3.bt2" -a -f $INDEX_RRNA".4.bt2" -a -f $INDEX_RRNA".rev.1.bt2" -a -f $INDEX_RRNA".rev.2.bt2" ];then
		echo0 3 "rRNA Index: "$INDEX_RRNA
	else
		echo0 0 "Error: no bowtie2 index in"${INDEX_RRNA}", please build with bowtie2-build"
		exit 1
	fi
fi

if [ -f $INDEX_TRANSPOSON".1.bt2" -a -f $INDEX_TRANSPOSON".2.bt2" -a -f $INDEX_TRANSPOSON".3.bt2" -a -f $INDEX_TRANSPOSON".4.bt2" -a -f $INDEX_TRANSPOSON".rev.1.bt2" -a -f $INDEX_TRANSPOSON".rev.2.bt2" ];then
	echo0 3 "transposon Index: "$INDEX_TRANSPOSON
else
	echo0 0 "Error: no bowtie2 index in"${INDEX_RRNA}", please build with bowtie2-build"
	exit 1
fi

####################################################
###align to rRNA and genome with bowtie2 and STAR###
####################################################


if [ ${RIGHT} == "0" ];then
	###single-end mode
	# bowtie2 remove rRNA
	if [ $SWITCH_RRNA -ne 0 ];then
		echo0 2 "---remove rRNA via bowtie2......"
		set -x
		bowtie2 --very-fast -p ${CPU} -q -x ${INDEX_RRNA} -U ${LEFT} -S ./bowtie2/${PREFIX}_rRNA.sam \
			--un ./bowtie2/${PREFIX}_RNAseq-rRNA.fq --no-unal > log_file/${PREFIX}.rRNA_removing.log 2>&1
		set +x
		LEFT="./bowtie2/"${PREFIX}"_RNAseq-rRNA.fq"
	fi
	# STAR map
	echo0 2 "---align via STAR......"
	set -x
	STAR 	--genomeDir ${INDEX_STAR} \
		--runThreadN ${CPU} \
		--readFilesIn ${LEFT} \
		--outFileNamePrefix ./STAR/${PREFIX} \
		--outFilterMismatchNmax ${MISMATCH} \
		--outFilterMultimapNmax ${MAPNUM} \
		--outSAMattributes All \
		--outFilterIntronMotifs RemoveNoncanonicalUnannotated
	set +x
else
	###pair-end mode
	# bowtie2 remove rRNA
	if [ $SWITCH_RRNA -ne 0 ];then
		echo0 2 "---remove rRNA via bowtie2......"
		set -x
		bowtie2 --very-fast -p ${CPU} -q -x ${INDEX_RRNA} -1 ${LEFT} -2 ${RIGHT} \
			-S ./bowtie2/${PREFIX}_rRNA.sam --un-conc ./bowtie2/${PREFIX}_RNAseq-rRNA.fq \
			--no-unal --no-discordant --no-mixed > log_file/${PREFIX}.rRNA_removing.log 2>&1
		set +x
		LEFT="./bowtie2/"${PREFIX}"_RNAseq-rRNA.1.fq"
		RIGHT="./bowtie2/"${PREFIX}"_RNAseq-rRNA.2.fq"
	fi
	echo -e "\033[40;36;1m---align via STAR......\033[0m"
	set -x
	STAR 	--genomeDir ${INDEX_STAR} \
		--runThreadN ${CPU} \
		--readFilesIn ${LEFT} ${RIGHT} \
		--outFileNamePrefix ./STAR/${PREFIX} \
		--outFilterMismatchNmax ${MISMATCH} \
		--outFilterMultimapNmax ${MAPNUM} \
		--outSAMattributes All \
		--outFilterIntronMotifs RemoveNoncanonicalUnannotated
	set +x
fi


##########################################################
###alignment sort, visualization and signal calculation###
##########################################################

###samtools
echo0 2 "---sam to indexed bam via samtools......"
samtools view -bhS -o ./STAR/${PREFIX}.bam ./STAR/${PREFIX}Aligned.out.sam
samtools sort -@ ${CPU} ./STAR/${PREFIX}.bam ./STAR/${PREFIX}.sort
if [ $SWITCH_RMDUP -ne 0 ];then
	samtools rmdup ./STAR/${PREFIX}.sort.bam ./STAR/${PREFIX}.sort.rmdup.bam > /dev/null 2>&1
	BAM=./STAR/${PREFIX}.sort.rmdup.bam
else
	BAM=./STAR/${PREFIX}.sort.bam
fi
samtools index ${BAM}
if [ -f ${BAM} ];then
	rm -rf ./bowtie2/${PREFIX}.bam ./STAR/${PREFIX}Aligned.out.sam
fi
#rm -rf bowtie2/${PREFIX}_RNAseq-rRNA.*
#rm -rf bowtie2/${PREFIX}_rRNA.sam

###if -L=miss, judge the library type
if [ "$LIBRARY" == "miss" ];then
	infer_experiment.py -i ${BAM} -r \
		${ANNOTATION_BED12} > STAR/${PREFIX}_exp_type.txt
	CHANCE_YES=`head -5 STAR/${PREFIX}_exp_type.txt | tail -1 | awk '{print $7}'`
	CHANCE_REVERSE=`head -6 STAR/${PREFIX}_exp_type.txt | tail -1 | awk '{print $7}'`
	LIBRARY=`awk -v CY=$CHANCE_YES -v CR=$CHANCE_REVERSE 'BEGIN{if(CY>0.6) {print "yes"} else if(CR>0.6) {print "reverse"} else {print "no"}}'`
	echo0 4 "WARNING: library is set to missing, guess is "$LIBRARY"......"
fi

###make bed12 file for unique mapped reads
echo0 2 "---bam to uniq bed via bedtools......"
if [ "${RIGHT}" = "0" ];then
	if [ "${LIBRARY}" = "yes" ];then
		bedtools bamtobed -i ${BAM} -bed12 | awk '$5>=10'> ./STAR/${PREFIX}.uniq.bed12
	else
		bedtools bamtobed -i ${BAM} -bed12 | awk 'BEGIN{FS=OFS="\t"} {if($5>=10){$6=($6=="+"?"-":"+");print $0}}' > ./STAR/${PREFIX}.uniq.bed12
	fi
else
	if [ "${LIBRARY}" = "yes" ];then
		bedtools bamtobed -i ${BAM} -bed12 | awk 'BEGIN{FS=OFS="\t"} {if($5>=10){if(substr($4,length($4))=="2"){if($6=="+"){$6="-"}else{$6="+"}};print $0}}' > ./STAR/${PREFIX}.uniq.bed12
	else
		bedtools bamtobed -i ${BAM} -bed12 | awk 'BEGIN{FS=OFS="\t"} {if($5>=10){if(substr($4,length($4))=="1"){if($6=="+"){$6="-"}else{$6="+"}};print $0}}' > ./STAR/${PREFIX}.uniq.bed12
	fi
fi

###confirm the normalize factor for density
#NUM_UNIQ_READ=`grep "Uniquely mapped reads number" STAR/${PREFIX}Log.final.out | awk 'BEGIN{FS="\t"} {print $2}'`
NUM_UNIQ_READ=(`wc -l ./STAR/${PREFIX}.uniq.bed12`)
FACTOR=`awk -v num=$NUM_UNIQ_READ 'BEGIN{factor=1000000/num;print factor}'`

###make bigWig and normalize
echo0 2 "---make signal files......"
PARA_FILE=${PREFIX}.parafly
if [ "${LIBRARY}" = "no" ];then
	echo "bedtools genomecov -scale ${FACTOR} -split -bg -i ./STAR/${PREFIX}.uniq.bed12 -g ${CHROMSIZE} | sort -k1,1 -k2,2n > ./bigWig/${PREFIX}.uniq.bdg && bedGraphToBigWig ./bigWig/${PREFIX}.uniq.bdg ${CHROMSIZE} ./bigWig/${PREFIX}.uniq.bw" > ${PARA_FILE}
	ParaFly -c ${PARA_FILE} -CPU ${CPU} -failed_cmds ${PARA_FILE}.failed_cmds 1>&2 && rm -rf ${PARA_FILE} ${PARA_FILE}.completed
else
	echo "bedtools genomecov -scale ${FACTOR} -split -bg -strand + -i ./STAR/${PREFIX}.uniq.bed12 -g ${CHROMSIZE} | sort -k1,1 -k2,2n > ./bigWig/${PREFIX}.uniq.watson.bdg && bedGraphToBigWig ./bigWig/${PREFIX}.uniq.watson.bdg ${CHROMSIZE} ./bigWig/${PREFIX}.uniq.watson.bw" > ${PARA_FILE} && echo "bedtools genomecov -scale ${FACTOR} -split -bg -strand - -i ./STAR/${PREFIX}.uniq.bed12 -g ${CHROMSIZE} | sort -k1,1 -k2,2n | awk 'BEGIN{FS=OFS=\"\t\"} {\$4=-\$4;print \$0}' > ./bigWig/${PREFIX}.uniq.crick.bdg && bedGraphToBigWig ./bigWig/${PREFIX}.uniq.crick.bdg ${CHROMSIZE} ./bigWig/${PREFIX}.uniq.crick.bw" >> ${PARA_FILE}
	ParaFly -c ${PARA_FILE} -CPU ${CPU} -failed_cmds ${PARA_FILE}.failed_cmds 1>&2 && rm -rf ${PARA_FILE} ${PARA_FILE}.completed
fi

#if [ "$LIBRARY" == "no" ];then
#	set -x
#	BB_BamToBw ${BAM} bigWig/${PREFIX}_none_nor.bw ${CHROMSIZE}
#	BB_NorBw.sh bigWig/${PREFIX}_none_nor.bw bigWig/${PREFIX}.bw ${FACOTR} ${CHROMSIZE} plus
#	set +x
#	rm bigWig/${PREFIX}_none_nor.bw 
#else
#	set -x
#	BB_SortBamToStrandBw_${TYPE}end.sh ${BAM} bigWig/ ${PREFIX}_none_nor ${CHROMSIZE} ${LIBRARY}
#	BB_NorBw.sh bigWig/${PREFIX}_none_nor_watson.bw bigWig/${PREFIX}_watson.bw ${FACOTR} ${CHROMSIZE} plus
#	BB_NorBw.sh bigWig/${PREFIX}_none_nor_crick.bw bigWig/${PREFIX}_crick.bw ${FACOTR} ${CHROMSIZE} minus 
#	set +x
#	rm bigWig/${PREFIX}_none_nor_watson.bw bigWig/${PREFIX}_none_nor_crick.bw 
#fi

###htseq-count
HTSEQ=$HTSEQ_COPY
NUM=0
if [ "$HTSEQ" != "0" ];then
	while true
	do
		NUM=`expr ${NUM} + 1`
		PARA_FILE=${PREFIX}.parafly
		echo0 2 "---signal calculation via htseq for"${HTSEQ%%,*}"......"
		set -x
		htseq-count -m ${HTSEQ_MODE} -f bam -s ${LIBRARY} -t exon -i ${ATTRIBUTE} --nonunique ${NONUNIQUE} -q ${BAM} ${HTSEQ%%,*} > ./signal/${PREFIX}.${NUM}.sig 2>/dev/null
		set +x
		BB_GtfToExonForEachGene.py ${HTSEQ%%,*} ${PREFIX}_temp.exon.bed
		awk 'BEGIN{FS=OFS="\t"} {a[$4][1]=$5;a[$4][2]+=($3-$2)} END{for(i in a){print i,a[i][1],a[i][2]}}' ${PREFIX}_temp.exon.bed > ${PREFIX}_temp.exonlen
		awk -v factor=${FACTOR} 'BEGIN{FS=OFS="\t"} {if(NR==FNR){a[$1]=$2}else{print $1,$2,a[$1]*1000/$3*factor}}' ./signal/${PREFIX}.${NUM}.sig ${PREFIX}_temp.exonlen > ./signal/${PREFIX}.${NUM}.rpkm 
		if [ ${SWITCH_PLOT} -eq 1 ];then
			echo0 2 "---signal plot via R for"${HTSEQ%%,*}"......"
			if [ "${PLOT_FILE}" == "0" ];then
				awk 'BEGIN{FS=OFS="\t"} {if(!a[$4][0]){a[$4][1]=$5;a[$4][2]=$2;a[$4][3]=$3;a[$4][0]=$1}else{a[$4][3]=$3}} END{for(i in a){print a[i][0],a[i][2],a[i][3],i,a[i][1]}}' ${PREFIX}_temp.exon.bed > ${PREFIX}_gene.bed5
				awk 'BEGIN{FS=OFS="\t"} {if(NR==FNR){if($3>0.01){a[$1]=1}}else{if(a[$4]){print $0}}}' ./signal/${PREFIX}.${NUM}.rpkm ${PREFIX}_gene.bed5 > ${PREFIX}_gene.filter.bed5
			else
				cp ${PLOT_FILE} ${PREFIX}_gene.filter.bed5
			fi
			GENE_NUM_PLOT=(`wc -l ${PREFIX}_gene.filter.bed5`)
			PARA_STEP=`expr ${GENE_NUM_PLOT} / ${CPU}`
			for ((i=0;i<=CPU;i++))
			do
				awk -v step=${PARA_STEP} -v sn=${i} 'NR>(sn)*step && NR<=(sn+1)*step' ${PREFIX}_gene.filter.bed5 > ${PREFIX}_gene.filter.${i}.bed5
			done
			PARA_FILE=${PREFIX}.plot.parafly
			# get trancript bed7 file if not specified
			if [ "${TRANSCRIPT_BED7}" = "0" ];then
				BB_GtfToBed7_For_Sushi.sh ${HTSEQ%%,*} ./plot_input/${PREFIX}_tran.bed7
				TRANSCRIPT_BED7=./plot_input/${PREFIX}_tran.bed7
			fi
			# plot signal
			if [ "${LIBRARY}" = "no" ];then
				for ((i=0;i<=CPU;i++))
				do
					echo "BB_plot_gene_dis.R ${PREFIX}_gene.filter.${i}.bed5 ${TRANSCRIPT_BED7} ./summary/${PREFIX}_gene_dis.${i}.pdf ./bigWig/${PREFIX}.uniq.bdg" >> ${PARA_FILE}
				done
			else
				for ((i=0;i<=CPU;i++))
				do
					echo "BB_plot_gene_dis.R ${PREFIX}_gene.filter.${i}.bed5 ${TRANSCRIPT_BED7} ./summary/${PREFIX}_gene_dis.${i}.pdf ./bigWig/${PREFIX}.uniq.watson.bdg ./bigWig/${PREFIX}.uniq.crick.bdg" >> ${PARA_FILE}
				done
			fi
			ParaFly -c ${PARA_FILE} -CPU ${CPU} -failed_cmds ${PARA_FILE}.failed_cmds 1>&2 && rm -rf ${PARA_FILE} ${PARA_FILE}.completed
			rm ${PREFIX}_gene.*bed5
			mv ${PREFIX}_gene.filter.bed5 ./plot_input/
			pdfjam --papersize '{6in,6in}' --outfile ./summary/${PREFIX}_gene_dis.pdf ./summary/${PREFIX}_gene_dis.*.pdf && rm -rf ./summary/${PREFIX}_gene_dis.*.pdf 
		fi
		rm -rf ${PREFIX}_temp.exon.bed ${PREFIX}_temp.exonlen
		last=$HTSEQ
		HTSEQ=${HTSEQ#*,}
		if [ "$HTSEQ" = "$last" ];then
			break
		fi
	done
fi

#########################
###transposon analysis###
#########################
if [ $SWITCH_TRANSPOSON -eq 0 ];then
	echo0 1 "finished!"
	exit 0
fi

echo0 1 "---start transposon analysis---"
###mapping with bowtie2
echo0 2 "---mapping to transposon consensus sequences via bowtie2......"
if [ "$RIGHT" == "0" ];then
	set -x
	bowtie2 -x ${INDEX_TRANSPOSON} -U ${LEFT} -a -X 800 -p ${CPU} -x ${INDEX_TRANSPOSON} -S ./bowtie2/${PREFIX}.transposon.sam --no-unal > log_file/${PREFIX}.transposon.log 2>&1
	set +x
else
	set -x
	bowtie2 -x ${INDEX_TRANSPOSON} -1 ${LEFT} -2 ${RIGHT} -a -X 800 -p ${CPU} -x ${INDEX_TRANSPOSON} -S ./bowtie2/${PREFIX}.transposon.sam --no-unal --no-discordant --no-mixed > log_file/${PREFIX}.transposon.log 2>&1
	set +x
fi
###samtools
echo0 2 "---sam to indexed bam via samtools......"
samtools view -bhS ./bowtie2/${PREFIX}.transposon.sam > ./bowtie2/${PREFIX}.transposon.bam 
samtools sort -@ ${CPU} ./bowtie2/${PREFIX}.transposon.bam ./bowtie2/${PREFIX}.transposon.sort
if [ ${SWITCH_RMDUP} -eq 0 ];then
	BAM=./bowtie2/${PREFIX}.transposon.sort.bam 
else
	samtools rmdup ./bowtie2/${PREFIX}.transposon.sort.bam ./bowtie2/${PREFIX}.transposon.rmdup.bam  
	BAM=./bowtie2/${PREFIX}.transposon.rmdup.bam 
fi
if [ -f ${BAM} ];then
	rm -rf ./bowtie2/${PREFIX}.transposon.bam ./bowtie2/${PREFIX}.transposon.sam
fi
###make bed12 file for unique mapped reads
echo0 2 "---bam to uniq bed via bedtools......"
if [ "${RIGHT}" = "0" ];then
	if [ "${LIBRARY}" = "yes" ];then
		bedtools bamtobed -i ${BAM} -bed12 | awk '$5>=10'> ./bowtie2/${PREFIX}.transposon.uniq.bed12
	else
		bedtools bamtobed -i ${BAM} -bed12 | awk 'BEGIN{FS=OFS="\t"} {if($5>=10){$6=($6=="+"?"-":"+");print $0}}' > ./bowtie2/${PREFIX}.transposon.uniq.bed12
	fi
else
	if [ "${LIBRARY}" = "yes" ];then
		bedtools bamtobed -i ${BAM} -bed12 | awk 'BEGIN{FS=OFS="\t"} {if($5>=10){if(substr($4,length($4))=="2"){if($6=="+"){$6="-"}else{$6="+"}};print $0}}' > ./bowtie2/${PREFIX}.transposon.uniq.bed12
	else
		bedtools bamtobed -i ${BAM} -bed12 | awk 'BEGIN{FS=OFS="\t"} {if($5>=10){if(substr($4,length($4))=="1"){if($6=="+"){$6="-"}else{$6="+"}};print $0}}' > ./bowtie2/${PREFIX}.transposon.uniq.bed12
	fi
fi
###make signal and bedGraph files
echo0 2 "---plot signal files......"
if [ "$LIBRARY" = "no" ];then
	awk -v factor=${FACOTR} 'BEGIN{FS=OFS="\t"} {if(NR==FNR){a[$1][0]=$2;a[$1][1]=0}else{a[$1][1]++}} END{for(i in a){print i,a[i][1]/factor*1000/a[i][0]}}' ${SIZE_TRANSPOSON} ./bowtie2/${PREFIX}.transposon.uniq.bed12 | sort -k1,1 > ./signal/${PREFIX}.transposon.rpkm
	echo "bedtools genomecov -scale ${FACTOR} -split -bg -i ./bowtie2/${PREFIX}.transposon.uniq.bed12 -g ${SIZE_TRANSPOSON} | sort -k1,1 -k2,2n > ./bigWig/${PREFIX}.transposon.uniq.bdg" > ${PARA_FILE}
	ParaFly -c ${PARA_FILE} -CPU ${CPU} -failed_cmds ${PARA_FILE}.failed_cmds && rm -rf ${PARA_FILE} ${PARA_FILE}.completed
	GENE_NUM_PLOT=(`wc -l ${SIZE_TRANSPOSON}`)
	PARA_STEP=`expr ${GENE_NUM_PLOT} / ${CPU}`
	for ((i=0;i<=CPU;i++))
	do
		awk -v step=${PARA_STEP} -v sn=${i} 'NR>(sn)*step && NR<=(sn+1)*step' ${SIZE_TRANSPOSON} > ${SIZE_TRANSPOSON}.${i}
	done
	PARA_FILE=${PREFIX}.plot.parafly
	for ((i=0;i<=CPU;i++))
	do
		echo "BB_plot_transposon_dis.R ${SIZE_TRANSPOSON}.${i} ./summary/${PREFIX}_transposon_dis.${i}.pdf ./bigWig/${PREFIX}.transposon.uniq.bdg" >> ${PARA_FILE}
	done
	ParaFly -c ${PARA_FILE} -CPU ${CPU} -failed_cmds ${PARA_FILE}.failed_cmds 1>&2 && rm -rf ${PARA_FILE} ${PARA_FILE}.completed
else
	awk -v factor=${FACTOR} 'BEGIN{FS=OFS="\t"} {if(NR==FNR){a[$1][0]=$2;a[$1][1]=0;a[$1][2]=0;a[$1][3]=0}else{a[$1][3]++;if($6=="+"){a[$1][1]++}else{a[$1][2]++}}} END{for(i in a){print i,a[i][1]*factor*1000/a[i][0],a[i][2]*factor*1000/a[i][0],a[i][3]*factor*1000/a[i][0]}}' ${SIZE_TRANSPOSON} ./bowtie2/${PREFIX}.transposon.uniq.bed12 | sort -k1,1 > ./signal/${PREFIX}.transposon.rpkm
	echo "bedtools genomecov -scale ${FACTOR} -split -bg -strand + -i ./bowtie2/${PREFIX}.transposon.uniq.bed12 -g ${SIZE_TRANSPOSON} | sort -k1,1 -k2,2n > ./bigWig/${PREFIX}.transposon.uniq.watson.bdg" > ${PARA_FILE} && echo "bedtools genomecov -scale ${FACTOR} -split -bg -strand - -i ./bowtie2/${PREFIX}.transposon.uniq.bed12 -g ${SIZE_TRANSPOSON} | sort -k1,1 -k2,2n | awk 'BEGIN{FS=OFS=\"\t\"} {\$4=-\$4;print \$0}' > ./bigWig/${PREFIX}.transposon.uniq.crick.bdg" >> ${PARA_FILE}
	ParaFly -c ${PARA_FILE} -CPU ${CPU} -failed_cmds ${PARA_FILE}.failed_cmds && rm -rf ${PARA_FILE} ${PARA_FILE}.completed
	GENE_NUM_PLOT=(`wc -l ${SIZE_TRANSPOSON}`)
	PARA_STEP=`expr ${GENE_NUM_PLOT} / ${CPU}`
	for ((i=0;i<=CPU;i++))
	do
		awk -v step=${PARA_STEP} -v sn=${i} 'NR>(sn)*step && NR<=(sn+1)*step' ${SIZE_TRANSPOSON} > ${SIZE_TRANSPOSON}.${i}
	done
	PARA_FILE=${PREFIX}.plot.parafly
	for ((i=0;i<=CPU;i++))
	do
		echo "BB_plot_transposon_dis.R ${SIZE_TRANSPOSON}.${i} ./summary/${PREFIX}_transposon_dis.${i}.pdf ./bigWig/${PREFIX}.transposon.uniq.watson.bdg ./bigWig/${PREFIX}.transposon.uniq.crick.bdg" >> ${PARA_FILE}
	done
	ParaFly -c ${PARA_FILE} -CPU ${CPU} -failed_cmds ${PARA_FILE}.failed_cmds 1>&2 && rm -rf ${PARA_FILE} ${PARA_FILE}.completed
fi
pdfjam --papersize '{6in,6in}' --outfile ./summary/${PREFIX}_transposon_dis.pdf ./summary/${PREFIX}_transposon_dis.*.pdf && rm -rf ./summary/${PREFIX}_transposon_dis.*.pdf 
rm -rf ${SIZE_TRANSPOSON}.*

###finished
echo0 1 "finished"
