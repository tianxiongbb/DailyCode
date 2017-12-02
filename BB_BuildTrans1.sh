#!/bin/sh

##----INTRO-----------##
# Name=Nov14.BuildTrans
# Date=Nov 14,2015
# Update=Nov 16,2015

########################
# Purpose
#This file is for build new transcriptome by RNA-seq.
#1.blat+psl filter(python)
###process with trinity_out_psl_file, fileter psl and get intron information.
#2.STAR+SAM/BAM+cufflinks
###RNAseq normal mapping and junction mapping, and assemble into transcript.

#######--Arguments--#######
help_info(){
	echo "usage:"
	echo "sh Nov14.BuildTrans.sh <option>* [-1 sample_1.fq] [-2 sample_2.fq] [-a sample.fa] [-l read_len] [-s species] [-n sample_name]"
	echo ""
	echo "Arguments:"
	echo "-1 The first fastq file for paried-end, or the fastq file for single-end"
	echo "-2 The second fastq file if sample is paried-end"
	echo "-a Trinity result fatsa file"
	echo "-s Species used for getting the right annotation files"
	echo "-n Sample name"
	echo "-o Overhang of cufflink arguement, usually read length-1, default: 100"
	echo "-C CPU number used for the pipeline"
	echo "-m minfrags of cufflink arguement, default: 10"
	echo ""
	echo "This file is for process the result after Trinity."
	echo ""
}


if [ $# -lt 1 ];then
	help_info
	exit 1
fi

samp_2_fq="NO"
CPU=1
min_frag=10

while getopts "1:2:a:s:n:l:C:m:" Arg
do
	case $Arg in
		1)	samp_1_fq=$OPTARG;;
		2)	samp_2_fq=$OPTARG;;
		a)	samp_fa=$OPTARG;;
		s)	species=$OPTARG;;
		n)	samp_name=$OPTARG;;
		l)  over_hang=$OPTARG;;
		C)  CPU=$OPTARG;;
		m)  min_frag=$OPTARG;;
		?)	echo "Wrong parameter!!!"
			exit 1;;
	esac
done

# get the right annotations based on species
case $species in
	mouse)		anno="/data/tongji2/annotation/Mus_musculus.GRCm38.82_seleno_filtered_final.gtf"
				geno_fa="/data/tongji2/genome/mm10.fa"
				STAR_inde="/data/tongji2/index/mm10_star_gtf/"
				bowt_inde="/data/tongji2/index/mm10_bowtie2/"
				chro_size="/data/tongji2/chrom.size/mm10.chrom.size";;
	rat)		anno="/data/tongji2/annotation/Rattus_norvegicus.Rnor_6.0.82_seleno_filtered_final.gtf"
				geno_fa="/data/tongji2/genome/rn6.fa"
				STAR_inde="/data/tongji2/index/rn6_star_gtf/"
				bowt_inde="/data/tongji2/index/rn6_bowtie2/"
				chro_size="/data/tongji2/chrom.size/rn6.chrom.size";;
	human)		anno="/data/tongji2/annotation/Homo_sapiens.GRCh38.82_seleno_filtered_final.gtf"
				geno_fa="/data/tongji2/genome/hg38.fa"
				STAR_inde="/data/tongji2/index/hg38_star_gtf/"
				bowt_inde="/data/tongji2/index/hg38_bowtie2/"
				chro_size="/data/tongji2/chrom.size/hg38.chrom.size";;
	cow)		anno="/data/tongji2/annotation/Bos_taurus.UMD3.1.82_seleno_filtered_final.gtf"
				geno_fa="/data/tongji2/genome/bosTau8.fa"
				STAR_inde="/data/tongji2/index/bosTau8_star_gtf/"
				bowt_inde="/data/tongji2/index/bosTau8_bowtie2/"
				chro_size="/data/tongji2/chrom.size/bosTau8.chrom.size";;
	rhesus)		anno="/data/tongji2/annotation/Macaca_mulatta.MMUL_1.82_seleno_filtered_final.gtf"
				geno_fa="/data/tongji2/genome/rheMac3.fa"
				STAR_inde="/data/tongji2/index/rheMac3_star_gtf/"
				bowt_inde="/data/tongji2/index/rheMac3_bowtie2/"
				chro_size="/data/tongji2/chrom.size/rheMac3.chrom.size";;
	chicken)	anno="/data/tongji2/annotation/Gallus_gallus.Galgal4.82_seleno_filtered_final.gtf"
				geno_fa="/data/tongji2/genome/galGal4.fa"
				STAR_inde="/data/tongji2/index/galGal4_star_gtf/"
				bowt_inde="/data/tongji2/index/galGal4_bowtie2/"
				chro_size="/data/tongji2/chrom.size/galGal4.chrom.size";;
	platypus)	anno="/data/tongji2/annotation/Ornithorhynchus_anatinus.OANA5.82_seleno_filtered_final.gtf"
				geno_fa="/data/tongji2/genome/ornAna1.fa"
				STAR_inde="/data/tongji2/index/ornAna1_star_gtf/"
				bowt_inde="/data/tongji2/index/ornAna1_bowtie2/"
				chro_size="/data/tongji2/chrom.size/ornAna1.chrom.size";;
	chimpanzee) anno="/data/tongji2/annotation/Pan_troglodytes.CHIMP2.1.4.82_seleno_filtered_final.gtf"
				geno_fa="/data/tongji2/genome/panTro4.fa"
				STAR_inde="/data/tongji2/index/panTro4_star_gtf/"
				bowt_inde="/data/tongji2/index/panTro4_bowtie2/"
				chro_size="/data/tongji2/chrom.size/panTro4.chrom.size";;
	gorilla)	anno="/data/tongji2/annotation/Gorilla_gorilla.gorGor3.1.82_seleno_filtered_final.gtf"
				geno_fa="/data/tongji2/genome/gorGor3.fa"
				STAR_inde="/data/tongji2/index/gorGor3_star_gtf/"
				bowt_inde="/data/tongji2/index/gorGor3_bowtie2/"
				chro_size="/data/tongji2/chrom.size/gorGor3.chrom.size";;
	opossum)	anno="/data/tongji2/annotation/Monodelphis_domestica.BROADO5.82_seleno_filtered_final.gtf"
				geno_fa="/data/tongji2/genome/monDom5.fa"
				STAR_inde="/data/tongji2/index/monDom5_star_gtf/"
				bowt_inde="/data/tongji2/index/monDom5_bowtie2/"
				chro_size="/data/tongji2/chrom.size/monDom5.chrom.size";;
	?)	echo "No Species!!!"
		exit 1;;
esac

pyth_path="/home/tongji2/piRNA/Code/Tianxiong/python_code/"
shel_path="/home/tongji2/piRNA/Code/Tianxiong/shell_code/"

###########################

echo "Start!"

###--mkdir--###
mkdir $samp_name
cd $samp_name

###--Process--###
#################
#1.blat+psl filter(python)
###process with trinity_out_psl_file, fileter psl and get intron information.
#-blat-#
#from Trinity_output fasta file to psl file
# blat -q=rna -out=psl $geno_fa $samp_fa ./${samp_name}.psl
# #-Python-#
# #filter psl file with %95 fraction mapped and %1 mismatch
# python ${pyth_path}tran_xlstxt_to_txt.py ./${samp_name}.psl ./temp.psl
# mv ./temp.psl ./${samp_name}.psl
# BB_PslFilter.py ./${samp_name}.psl ./temp.psl 0.95 0.01
# mv ./temp.psl ./${samp_name}.psl
# #-Python-#
# #get intron from psl file
# BB_GetIntronFromPsl.py ./${samp_name}.psl ./${samp_name}_trinity_intron.tab
# #merge introns we got and introns annotated
# cat ./${samp_name}_trinity_intron.tab ${STAR_inde}sjdbList.out.tab | sort -u > ${samp_name}_merge_intron.tab
#################
#2.STAR+SAM/BAM+cufflinks
###RNAseq normal mapping and junction mapping, and assemble into transcript.
#-STAR-#
#splicing junction based STAR for all reads from sample.fq
#use 16 CPU, 2 mismatch, 100 max multi-map, , write all sam attributes, write unmapped reads
# mkdir ./temp_index
mkdir ./star_sj
# # build index
# STAR 	--runMode genomeGenerate \
# 		--runThreadN $CPU \
# 		--sjdbFileChrStartEnd ./${samp_name}_merge_intron.tab \
# 		--sjdbOverhang $over_hang \
# 		--genomeFastaFiles $geno_fa \
# 		--genomeDir ./temp_index/

# run STAR
if [ $samp_2_fq = 'NO' ]; then
#sangle-end
	STAR 	--genomeDir /data/tongji2/Pipeline_get_piRNA_transcript/mouse/mouse_testes_star_index/ \
			--runThreadN $CPU \
			--readFilesIn $samp_1_fq \
			--outFileNamePrefix ./star_sj/${samp_name} \
			--outFilterMismatchNmax 2 \
			--outFilterMultimapNmax 100 \
			--outSAMattributes All \
			--outFilterIntronMotifs RemoveNoncanonicalUnannotated
else
#paired-end 
	STAR 	--genomeDir /data/tongji2/Pipeline_get_piRNA_transcript/mouse/mouse_testes_star_index/ \
			--runThreadN $CPU \
			--readFilesIn $samp_1_fq $samp_2_fq \
			--outFileNamePrefix ./star_sj/${samp_name} \
			--outFilterMismatchNmax 2 \
			--outFilterMultimapNmax 100 \
			--outSAMattributes All \
			--outFilterIntronMotifs RemoveNoncanonicalUnannotated
fi

# remove index file for space_free
rm -rf temp_index

#-samtools-#
#transfer sam to bam, sort.
samtools view -bhS -o ./star_sj/${samp_name}.bam ./star_sj/${samp_name}Aligned.out.sam
samtools sort ./star_sj/${samp_name}.bam ./star_sj/${samp_name}.sort
rm -rf ./star_*/*.sam
rm -rf ./star_sj/${samp_name}.bam
#-cufflink-#
#transcriptome assembling
mkdir ./cufflink
cufflinks 	-o ./cufflink/ \
			-g $anno \
			-p $CPU \
			--library-type fr-firststrand \
			-u -j 0.2 \
			--min-frags-per-transfrag $min_frag \
			--overlap-radius 250 \
			./star_sj/${samp_name}.sort.bam > /dev/null 2>&1
#run cuffcompare to check the accuracy of assembled transcript
#the main performance of cuffcompare is in .states file
#and we can get the geneId/transcriptId from .transcripts.gtf.tmap (to be doing)
cuffcompare	-r $anno \
			-o $samp_name \
			./cufflink/transcripts.gtf
#-Python-#
#transfer GTF file to bed file
#python ${pyth_path}gtfToBed_Dec2.py ./cufflink/transcripts.gtf ./${samp_name}.transcripts_assemble.bed

echo "Finish!"

#########################
