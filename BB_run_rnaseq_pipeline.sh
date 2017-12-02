#!/bin/bash

if [ $# -lt 6 ];then
	echo -e "$0 prefix srr_id num_of_lane single/paired yes/reverse/no/miss genome"
	echo -e "set num_of_lane to 0 if the SRR file had already downloaded and gunzipped"
	exit 1
fi
PREFIX=$1
SRR_ID=$2
NUM_LANE=$3
TYPE=$4
STRAND=$5
GENOME=$6

if [ $NUM_LANE -eq 1 ];then
	TYPE_OTHER_LANE=${TYPE}
elif [ $NUM_LANE -eq 2 ];then
	HEAD=${SRR_ID:0:3}
	TAIL=${SRR_ID:3:7}
	TAIL=`expr $TAIL + 1`
	TYPE_OTHER_LANE=${TYPE}" "${HEAD}${TAIL}
elif [ $NUM_LANE -eq 3 ];then
	HEAD=${SRR_ID:0:3}
	TAIL=${SRR_ID:3:7}
	TAIL=`expr $TAIL + 1`
	TYPE_OTHER_LANE=${TYPE}" "${HEAD}${TAIL}
	TAIL=`expr $TAIL + 1`
	TYPE_OTHER_LANE=${TYPE_OTHER_LANE}" "${HEAD}${TAIL}
elif [ $NUM_LANE -eq 4 ];then
	HEAD=${SRR_ID:0:3}
	TAIL=${SRR_ID:3:7}
	TAIL=`expr $TAIL + 1`
	TYPE_OTHER_LANE=${TYPE}" "${HEAD}${TAIL}
	TAIL=`expr $TAIL + 1`
	TYPE_OTHER_LANE=${TYPE_OTHER_LANE}" "${HEAD}${TAIL}
	TAIL=`expr $TAIL + 1`
	TYPE_OTHER_LANE=${TYPE_OTHER_LANE}" "${HEAD}${TAIL}
elif [ $NUM_LANE -eq 5 ];then
	HEAD=${SRR_ID:0:3}
	TAIL=${SRR_ID:3:7}
	TAIL=`expr $TAIL + 1`
	TYPE_OTHER_LANE=${TYPE}" "${HEAD}${TAIL}
	TAIL=`expr $TAIL + 1`
	TYPE_OTHER_LANE=${TYPE_OTHER_LANE}" "${HEAD}${TAIL}
	TAIL=`expr $TAIL + 1`
	TYPE_OTHER_LANE=${TYPE_OTHER_LANE}" "${HEAD}${TAIL}
	TAIL=`expr $TAIL + 1`
	TYPE_OTHER_LANE=${TYPE_OTHER_LANE}" "${HEAD}${TAIL}
elif [ $NUM_LANE -eq 6 ];then
	HEAD=${SRR_ID:0:3}
	TAIL=${SRR_ID:3:7}
	TAIL=`expr $TAIL + 1`
	TYPE_OTHER_LANE=${TYPE}" "${HEAD}${TAIL}
	TAIL=`expr $TAIL + 1`
	TYPE_OTHER_LANE=${TYPE_OTHER_LANE}" "${HEAD}${TAIL}
	TAIL=`expr $TAIL + 1`
	TYPE_OTHER_LANE=${TYPE_OTHER_LANE}" "${HEAD}${TAIL}
	TAIL=`expr $TAIL + 1`
	TYPE_OTHER_LANE=${TYPE_OTHER_LANE}" "${HEAD}${TAIL}
	TAIL=`expr $TAIL + 1`
	TYPE_OTHER_LANE=${TYPE_OTHER_LANE}" "${HEAD}${TAIL}
else
	echo "Warning: too many lanes"
fi

# download data
cd /data/tusers/yutianx/tongji2/Recycle/RNAseq/Data
if [ $NUM_LANE -gt 0 ];then
	BB_download_GEO.sh ${SRR_ID} ${TYPE_OTHER_LANE}
fi
# run pipeline
if [ "$TYPE" == "single" ];then
	BB_rnaseq_pipeline.sh -l ${SRR_ID}.fastq -g ${GENOME} -o /data/tusers/yutianx/tongji2/Recycle/RNAseq/Output \
		-L ${STRAND} -A gene_name -p ${PREFIX}
elif [ "$TYPE" == "paired" ];then
	BB_rnaseq_pipeline.sh -l ${SRR_ID}_1.fastq -r ${SRR_ID}_2.fastq -g ${GENOME} \
		-o /data/tusers/yutianx/tongji2/Recycle/RNAseq/Output -L ${STRAND} -A gene_name -p ${PREFIX}
else
	echo "wrong type, single/paired"
fi
# check if it has really ran
RESULT_LENGTH=(`wc -l /data/tusers/yutianx/tongji2/Recycle/RNAseq/Output/Htseq/${PREFIX}.rpm`)
if [ $RESULT_LENGTH -lt 2 ];then
	if [ ! -f /data/tusers/yutianx/tongji2/Recycle/RNAseq/error.list ];then
		echo -e "sample_name\tsrr_id\tlane_num\ttype\tstrand\tgenome" > \
			/data/tusers/yutianx/tongji2/Recycle/RNAseq/error.list
	fi
	echo -e "${PREFIX}\t${SRR_ID}\t${NUM_LANE}\t${TYPE}\t${STRAND}\t${GENOME}" \
		>> /data/tusers/yutianx/tongji2/Recycle/RNAseq/error.list
else
	if [ ! -f /data/tusers/yutianx/tongji2/Recycle/RNAseq/finish.list ];then
		echo -e "sample_name\tsrr_id\tlane_num\ttype\tstrand\tgenome" > \
			/data/tusers/yutianx/tongji2/Recycle/RNAseq/finish.list
	fi
	echo -e "${PREFIX}\t${SRR_ID}\t${NUM_LANE}\t${TYPE}\t${STRAND}\t${GENOME}" \
		>> /data/tusers/yutianx/tongji2/Recycle/RNAseq/finish.list
fi
# remove big files
#if [ "$TYPE" == "single" ];then
#	rm /data/tongji2/Recycle/RNAseq/Data/${SRR_ID}*.fastq
#else
#	rm /data/tongji2/Recycle/RNAseq/Data/${SRR_ID}_1.fastq \
#		/data/tongji2/Recycle/RNAseq/Data/${SRR_ID}_2.fastq 
#fi
#rm /data/tongji2/Recycle/RNAseq/Output/STAR/${PREFIX}Log.out \
#	/data/tongji2/Recycle/RNAseq/Output/STAR/${PREFIX}Log.progress.out \
#	/data/tongji2/Recycle/RNAseq/Output/STAR/${PREFIX}*.bam
#
