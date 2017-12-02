#!/bin/bash

help_info(){
	echo "usage:"
	echo "sh BB_PslToGtf.sh in.psl out.gtf"
	echo ""
}

if [ $# -lt 1 ];then
	help_info
	exit 1
fi

pslToBed $1 temp_$1.bed
bedToGenePred temp_$1.bed temp_$1.genePred
genePredToGtf file temp_$1.genePred $2

rm temp_$1.bed temp_$1.genePred