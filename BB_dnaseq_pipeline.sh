#!/bin/bash

if [ $# -lt 3 ];then
	echo0 1 $0" left.fq rigth.fq out.prefix CPU out.dir if_fixed_by_te_cor(0|1)"
	echo0 4 "please modify the configure path in this program"
	exit 1
fi

#configures
TE_FA=/data/tusers/yutianx/tongji2/Annotation/Fasta/phaCin_repbase_denovo.fa
BWA_INDEX=/data/tusers/yutianx/tongji2/Annotation/Index/phaCin0_bwa/genome
BOWTIE2_INDEX=/data/tusers/yutianx/tongji2/Annotation/Index/phaCin0_bowtie2/genome
TE_BED=/data/tusers/yutianx/tongji2/piRNA/Koala_project/RepeatMask/RepearMasker_slow_withrepbase/phaCin_unsw_v4.1.fa.bed
FRAGMENT_SIZE=500

#run
LEFT=`readlink -f $1`
RIGHT=`readlink -f $2`
mkdir $5
cd $5
bwa mem -T 0 -t $4 ${BWA_INDEX} $LEFT $RIGHT > $3.sam
#bowtie2 -x ${BOWTIE2_INDEX} -1 $LEFT -2 $RIGHT --local -S $3.sam --un-conc $3.unpair.uniq.fastq -p 4 > log_$3_bowtie2 2>&1
#bwa aln -n 3 -l 100 -R 10000 -t $4 /data/tusers/yutianx/tongji2/Annotation/Index/dm3_bwa/genome $1 > $3.1.sai
#bwa aln -n 3 -l 100 -R 10000 -t $4 /data/tusers/yutianx/tongji2/Annotation/Index/dm3_bwa/genome $1 > $3.2.sai
#bwa sampe /data/tusers/yutianx/tongji2/Annotation/Index/dm3_bwa/genome $3.1.sai $3.2.sai $1 $2 > $3.sam
samtools view -bhS $3.sam > $3.bam
samtools sort -@ $4 $3.bam $3.sorted
rm $3.sam $3.bam
samtools index $3.sorted.bam
if [ $6 -gt 0 ];then
	awk -v fl=$FRAGMENT_SIZE 'BEGIN{FS=OFS="\t"} {if($1~/^>/){name=substr($1,2);a[name]=""}else{a[name]=a[name]""$1}} END{for(i in a){print ">"i;if(length(a[i])<=2*fl){print a[i]}else{print substr(a[i],1,fl)""substr(a[i],length(a[i])-fl+1)}}}' ${TE_FA} > chopped.te.fa
	TE_FA=chopped.te.fa
fi
/data/tusers/yutianx/tongji2/Software/TEMP/scripts/TEMP_Insertion.sh -i $3.sorted.bam -s /data/tusers/yutianx/tongji2/Software/TEMP/scripts/ -r ${TE_FA} -c $4 -t ${TE_BED} -f $FRAGMENT_SIZE -x 5
awk 'BEGIN{FS=OFS="\t"} {if($1~/^>/){name=substr($1,2);a[name]=0}else{a[name]+=length($1)}} END{for(i in a){print i,a[i]}}' ${TE_FA} > tmp.te.size
awk 'BEGIN{FS=OFS="\t";a["sense"]="+";a["antisense"]="-"} {if(NR>1){print $1,$2,$3,$4,$6,a[$5],$8}}' ${3}.insertion.refined.bp.summary > ${3}.final.bed7 
awk 'BEGIN{FS=OFS="\t"} {if(NR==FNR){a[$4]["all"]+=$7;a[$4][$5]+=$7}else{print $1,a[$1]["all"]?a[$1]["all"]:0,a[$1]["1p1"]?a[$1]["1p1"]:0,a[$1]["2p"]?a[$1]["2p"]:0,a[$1]["singleton"]?a[$1]["singleton"]:0}}' ${3}.final.bed7 tmp.te.size > ${3}.final.result




