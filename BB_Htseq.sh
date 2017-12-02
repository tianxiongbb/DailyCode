#!/bin/bash
if [ $# -lt 6 ];then
	echo -e "\033[1;31;33m"$0" in.bam gene_id/gene_name reverse/yes/no in.gtf in.exonlen out.prefix\033[0m"
	exit 1
fi

htseq-count -f bam -s $3 -i $2 -t exon -m intersection-nonempty --nonunique all -q $1 $4 > $6.sig
BB_NorHtseqRPM.py $6.sig $6.rpm
awk 'BEGIN{FS=OFS="\t"} {if(NR==FNR){a[$1]=$2}else{if(a[$1]){print $1,$2*1000/a[$1]}}}' $5 $6.rpm > $6.rpkm

