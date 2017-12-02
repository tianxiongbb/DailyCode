#!/bin/bash

#cd /data/tongji2/RawData/smallRNA
# nohup cutadapt -a TCGTATGCCGTCTTCTGCTTG -e 0.15 -m 15 -O 4 -o chicken_cutadapt.fastq chicken.fastq &
# nohup cutadapt -a NNNNAGATCGGA -e 0.15 -m 15 -O 7 -o cow_cutadapt_1.fastq cow_1.fastq &
# nohup cutadapt -a NNNNAGATCGGA -e 0.15 -m 15 -O 7 -o cow_cutadapt_2.fastq cow_2.fastq &
# nohup cutadapt -a ATCTCGTATGCCGTCTTCTGCTTG -e 0.15 -m 15 -O 5 -o human_cutadapt.fastq human.fastq &
#nohup cutadapt -a AGATCGGAAGAGCACACGTCTGAACTCC -e 0.15 -m 15 -O 7 -o marmoset_cutadapt.fastq marmoset_adult_smallRNA.fastq &
# nohup cutadapt -a ATCTCGTATGCCGTCTTCTGCTTG -e 0.15 -m 15 -O 5 -o opossum_cutadapt.fastq opossum.fastq &
# nohup cutadapt -a ATCTCGTATGCCGTCTTCTGCTTG -e 0.15 -m 15 -O 5 -o platypus_cutadapt.fastq platypus.fastq &
# nohup cutadapt -a TCGTATGCCGTCTTCTGCTTG -e 0.15 -m 15 -O 2 -o rat_cutadapt.fastq rat.fastq &
# nohup cutadapt -a AGATCGGAAGAGCACACGTCTGAACTCC -e 0.15 -m 15 -O 7 -o rhesus_cutadapt_1.fastq rhesus_1.fastq &
# nohup cutadapt -a AGATCGGAAGAGCACACGTCTGAACTCC -e 0.15 -m 15 -O 7 -o rhesus_cutadapt_2.fastq rhesus_2.fastq &
# nohup cutadapt -a TCGTATGCCGTCTTCTGCTTG -e 0.15 -m 15 -O 4 -o mouse_10.5dpp_cutadapt.fastq mouse_10.5dpp_smallRNA.fastq &
# nohup cutadapt -a TCGTATGCCGTCTTCTGCTTG -e 0.15 -m 15 -O 7 -o mouse_12.5dpp_cutadapt.fastq mouse_12.5dpp_smallRNA.fastq &
# nohup cutadapt -a TCGTATGCCGTCTTCTGCTTG -e 0.15 -m 15 -O 7 -o mouse_14.5dpp_cutadapt.fastq mouse_14.5dpp_smallRNA.fastq &
# nohup cutadapt -a TCGTATGCCGTCTTCTGCTTG -e 0.15 -m 15 -O 7 -o mouse_17.5dpp_cutadapt.fastq mouse_17.5dpp_smallRNA.fastq &
# nohup cutadapt -a TCGTATGCCGTCTTCTGCTTG -e 0.15 -m 15 -O 7 -o mouse_20.5dpp_cutadapt.fastq mouse_20.5dpp_smallRNA.fastq &
# nohup cutadapt -a TCGTATGCCGTCTTCTGCTTG -e 0.15 -m 15 -O 7 -o mouse_6weeks_cutadapt.fastq mouse_6weeks_smallRNA.fastq &

# nohup BB_Pip_SRNA.sh -i /data/tongji2/RawData/smallRNA/opossum_cutadapt.fastq -g monDom5 -c 16 -o /data/tongji2/piRNA/Output/Pip_SRNA/opossum -p opossum > opossum.log &
# nohup BB_Pip_SRNA.sh -i /data/tongji2/RawData/smallRNA/human_cutadapt.fastq -g hg38 -c 16 -o /data/tongji2/piRNA/Output/Pip_SRNA/human -p human > human.log &
# nohup BB_Pip_SRNA.sh -i /data/tongji2/RawData/smallRNA/pig_cutadapt.fastq -g susScr3 -c 8 -o /data/tongji2/piRNA/Output/Pip_SRNA/pig -p pig > pig.log &

# sp=human
# ass=hg38
# nohup BB_IntersectBed2Bed.py /data/tongji2/piRNA/Output/Pip_SRNA/${sp}.piRNA.bed2 /data/tongji2/piRNA/Output/PslMap/20kb_mm10_${ass}.bed /data/tongji2/piRNA/Conservation_Analysis/Intersect/mm10_${ass} plus /data/tongji2/piRNA/Conservation_Analysis/PiRNA_Loci/mm10.list &
# sp=mouse_6weeks
# ass=mm10
# nohup BB_IntersectBed2Bed.py /data/tongji2/piRNA/Output/Pip_SRNA/${sp}.piRNA.bed2 /data/tongji2/piRNA/Output/PslMap/20kb_mm10_${ass}.bed /data/tongji2/piRNA/Conservation_Analysis/Intersect/mm10_${ass} plus /data/tongji2/piRNA/Conservation_Analysis/PiRNA_Loci/mm10.list &
# sp=cow
# ass=bosTau8
# nohup BB_IntersectBed2Bed.py /data/tongji2/piRNA/Output/Pip_SRNA/${sp}.piRNA.bed2 /data/tongji2/piRNA/Output/PslMap/20kb_mm10_${ass}.bed /data/tongji2/piRNA/Conservation_Analysis/Intersect/mm10_${ass} plus /data/tongji2/piRNA/Conservation_Analysis/PiRNA_Loci/mm10.list &
# sp=marmoset
# ass=calJac3
# nohup BB_IntersectBed2Bed.py /data/tongji2/piRNA/Output/Pip_SRNA/${sp}.piRNA.bed2 /data/tongji2/piRNA/Output/PslMap/20kb_mm10_${ass}.bed /data/tongji2/piRNA/Conservation_Analysis/Intersect/mm10_${ass} plus /data/tongji2/piRNA/Conservation_Analysis/PiRNA_Loci/mm10.list &
# sp=chicken
# ass=galGal4
# nohup BB_IntersectBed2Bed.py /data/tongji2/piRNA/Output/Pip_SRNA/${sp}.piRNA.bed2 /data/tongji2/piRNA/Output/PslMap/20kb_mm10_${ass}.bed /data/tongji2/piRNA/Conservation_Analysis/Intersect/mm10_${ass} plus /data/tongji2/piRNA/Conservation_Analysis/PiRNA_Loci/mm10.list &
# sp=platypus
# ass=ornAna1
# nohup BB_IntersectBed2Bed.py /data/tongji2/piRNA/Output/Pip_SRNA/${sp}.piRNA.bed2 /data/tongji2/piRNA/Output/PslMap/20kb_mm10_${ass}.bed /data/tongji2/piRNA/Conservation_Analysis/Intersect/mm10_${ass} plus /data/tongji2/piRNA/Conservation_Analysis/PiRNA_Loci/mm10.list &
# sp=rhesus
# ass=rheMac3
# nohup BB_IntersectBed2Bed.py /data/tongji2/piRNA/Output/Pip_SRNA/${sp}.piRNA.bed2 /data/tongji2/piRNA/Output/PslMap/20kb_mm10_${ass}.bed /data/tongji2/piRNA/Conservation_Analysis/Intersect/mm10_${ass} plus /data/tongji2/piRNA/Conservation_Analysis/PiRNA_Loci/mm10.list &
# sp=rat
# ass=rn6
# nohup BB_IntersectBed2Bed.py /data/tongji2/piRNA/Output/Pip_SRNA/${sp}.piRNA.bed2 /data/tongji2/piRNA/Output/PslMap/20kb_mm10_${ass}.bed /data/tongji2/piRNA/Conservation_Analysis/Intersect/mm10_${ass} plus /data/tongji2/piRNA/Conservation_Analysis/PiRNA_Loci/mm10.list &
# sp=opossum
# ass=monDom5
# nohup BB_IntersectBed2Bed.py /data/tongji2/piRNA/Output/Pip_SRNA/${sp}.piRNA.bed2 /data/tongji2/piRNA/Output/PslMap/20kb_mm10_${ass}.bed /data/tongji2/piRNA/Conservation_Analysis/Intersect/mm10_${ass} plus /data/tongji2/piRNA/Conservation_Analysis/PiRNA_Loci/mm10.list &
# sp=pig
# ass=susScr3
# nohup BB_IntersectBed2Bed.py /data/tongji2/piRNA/Output/Pip_SRNA/${sp}.piRNA.bed2 /data/tongji2/piRNA/Output/PslMap/20kb_mm10_${ass}.bed /data/tongji2/piRNA/Conservation_Analysis/Intersect/mm10_${ass} plus /data/tongji2/piRNA/Conservation_Analysis/PiRNA_Loci/mm10.list &
# sp=rabbit
# ass=oryCun2
# nohup BB_IntersectBed2Bed.py /data/tongji2/piRNA/Output/Pip_SRNA/${sp}.piRNA.bed2 /data/tongji2/piRNA/Output/PslMap/20kb_mm10_${ass}.bed /data/tongji2/piRNA/Conservation_Analysis/Intersect/mm10_${ass} plus /data/tongji2/piRNA/Conservation_Analysis/PiRNA_Loci/mm10.list &

OutPath=/data/tongji2/piRNA/Output/Htseq/AllRNA
InPath=/data/tongji2/piRNA/Output/Transcriptome
AnnoPath=/data/tongji2/InputForRunBT/Reference
#human
nohup htseq-count -f bam -s reverse -t exon -i gene_id -q \
${InPath}/human/human_testes_a/human_testes_a.merge.sort.rmdup.bam \
${AnnoPath}/Homo_sapiens.GRCh38.82_seleno_filtered_final.gtf > \
${OutPath}/human_a_testes_gene.txt 2>human1.log &
nohup htseq-count -f bam -s reverse -t exon -i gene_id -q \
${InPath}/human/human_testes_b/human_testes_b.merge.sort.rmdup.bam \
${AnnoPath}/Homo_sapiens.GRCh38.82_seleno_filtered_final.gtf > \
${OutPath}/human_b_testes_gene.txt 2>human2.log &
#marmoset
nohup htseq-count -f bam -s reverse -t exon -i gene_id -q \
${InPath}/marmoset/marmoset_adult_testes/marmoset_adult_testes.merge.sort.rmdup.bam \
${AnnoPath}/Callithrix_jacchus.C_jacchus3.2.1.83_seleno_filtered_final.gtf > \
${OutPath}/marmoset_testes_gene.txt 2>marmoset1.log &
#opossum
nohup htseq-count -f bam -s yes -t exon -i gene_id -q \
${InPath}/opossum/opossum_testes_a/opossum_testes_a.merge.sort.rmdup.bed \
${AnnoPath}/Monodelphis_domestica.BROADO5.82_seleno_filtered_final.gtf > \
${OutPath}/opossum_a_testes_gene.txt 2>opossum1.log &
nohup htseq-count -f bam -s yes -t exon -i gene_id -q \
${InPath}/opossum/opossum_testes_b/opossum_testes_b.merge.sort.rmdup.bed \
${AnnoPath}/Monodelphis_domestica.BROADO5.82_seleno_filtered_final.gtf > \
${OutPath}/opossum_b_testes_gene.txt 2>opossum1.log &
#platypus
nohup htseq-count -f bam -s yes -t exon -i gene_id -q \
${InPath}/platypus/platypus_testes_a/platypus_testes_a.merge.sort.rmdup.bam \
${AnnoPath}/Ornithorhynchus_anatinus.OANA5.82_seleno_filtered_final.gtf > \
${OutPath}/platypus_a_testes_gene.txt 2>platypus1.log &
nohup htseq-count -f bam -s yes -t exon -i gene_id -q \
${InPath}/platypus/platypus_testes_b/platypus_testes_b.merge.sort.rmdup.bam \
${AnnoPath}/Ornithorhynchus_anatinus.OANA5.82_seleno_filtered_final.gtf > \
${OutPath}/platypus_b_testes_gene.txt 2>platypus1.log &
#rabbit
nohup htseq-count -f bam -s reverse -t exon -i gene_id -q \
${InPath}/rabbit/rabbit_adult_testes/rabbit_adult_testes.merge.sort.rmdup.bam \
${AnnoPath}/Oryctolagus_cuniculus.OryCun2.0.83_seleno_filtered_final.gtf > \
${OutPath}/rabbit_testes_gene.txt 2>rabbit1.log &




