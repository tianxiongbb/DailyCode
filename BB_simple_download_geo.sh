#!/bin/bash

if [ $# -lt 1 ];then
	echo $0" SRR/ERR/SRX_id [sample_name]"
	exit 1
fi

GID=$1
wget -nd -q ftp://ftp-trace.ncbi.nih.gov/sra/sra-instant/reads/ByRun/sra/${GID:0:3}/${GID:0:6}/${GID}/${GID}.sra
if [ $# -gt 1 ];then
	mv ${GID}.sra ${2}.sra
fi
