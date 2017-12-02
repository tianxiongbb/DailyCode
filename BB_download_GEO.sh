#!/bin/bash

# help information
if [ $# -lt 1 ];then
	echo0 1 "usage:"
	echo0 1 "    "$0" SRR/ERR/DRR_ID single/paired prefix.name [SRR/ERR/DRR_ID_IN_DIFF_LANE]"
	exit 1
fi

START_ID=$1
STRAND=$2
PREFIX=$3
# download data
GID=$1
wget -nd -q ftp://ftp-trace.ncbi.nih.gov/sra/sra-instant/reads/ByRun/sra/${GID:0:3}/${GID:0:6}/${GID}/${GID}.sra
# fastq-dump data
if [ "${STRAND}" == "single" ];then
	fastq-dump.2.8.2 ${GID}.sra
else
	fastq-dump.2.8.2 --split-files ${GID}.sra
fi
rm ${GID}.sra
# for each data in different lane
if [ $# -gt 3 ];then
	shift
	GID=$3
	wget -nd -q ftp://ftp-trace.ncbi.nih.gov/sra/sra-instant/reads/ByRun/sra/${GID:0:3}/${GID:0:6}/${GID}/${GID}.sra
	if [ "${STRAND}" == "single" ];then
		fastq-dump.2.8.2 ${GID}.sra
		cat ${GID}.fastq >> ${START_ID}.fastq
	else
		fastq-dump.2.8.2 --split-files ${GID}.sra
		cat ${GID}_1.fastq >> ${START_ID}_1.fastq
		cat ${GID}_2.fastq >> ${START_ID}_2.fastq
	fi
	rm ${GID}.sra
	rm ${GID}*.fastq
fi
# rename data
if [ "$2" == "single" ];then
	mv ${START_ID}.fastq ${PREFIX}.fastq
else
	mv ${START_ID}_1.fastq ${PREFIX}_1.fastq
	mv ${START_ID}_2.fastq ${PREFIX}_2.fastq
fi





