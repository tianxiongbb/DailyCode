#!/bin/bash

help_info(){
	echo "usage:"
	echo "sh BB_SmallRNAPipeline.sh in.fq out.dir out.prefix species_assemble(eg:mm10)"
	echo ""
}

if [ $# -lt 1 ];then
	help_info
	exit 1
fi

cd /data/tongji2/piPipes/bin/
./piPipes_fastq_to_insert $1 $2/$3.insert
bowtie2 -r -n -x 	



