#!/usr/bin/env Rscript

###library packages###

###read data###
Args <- commandArgs()
f_map <- Args[6]
f_all_len <- Args[7]
f_all_str <- Args[8]
f_uniq_len <- Args[9]
f_uniq_str <- Args[10]
f_multi_len <- Args[11]
f_multi_str <- Args[12]

map <- read.table(f_map,header=TRUE,row.names=1)
all_len <- read.table(f_all_len,header=FALSE,row.names=1)
all_str <- read.table(f_all_str,header=FALSE,row.names=1)
uniq_len <- read.table(f_uniq_len,header=FALSE,row.names=1)
uniq_str <- read.table(f_uniq_str,header=FALSE,row.names=1)
multi_len <- read.table(f_multi_len,header=FALSE,row.names=1)
multi_str <- read.table(f_multi_str,header=FALSE,row.names=1)

###plot picture###
###open pdf
pdf("Summary.pdf",height=8,width=8)

###mapping rate
l=c(paste(round((map[1,1]-map[1,2])/map[1,1],4)*100,"%",sep=""),
	paste(round(map[1,3]/map[1,1],4)*100,"%",sep=""),
	paste(round(map[1,4]/map[1,1],4)*100,"%",sep=""))

pie(c(map[1,1]-map[1,2],map[1,3],map[1,4]),main="piechart of mapping rate",
	label=l,col=c("#1f78b4","#e31a1c","#33a02c"))
legend("topleft",
	legend=c("Unmapped Reads","Unique Mapped Reads","Multiple Mapped Reads"),
	col=c("#1f78b4","#e31a1c","#33a02c"),pch=15,cex=0.9)

###RepeatMask mapping rate
#all mapping reads
l=c(paste(round((map[2,2]+map[8,2])/sum(map[2:13,2]),4)*100,"%",sep=""),
	paste(round((map[3,2]+map[9,2])/sum(map[2:13,2]),4)*100,"%",sep=""),
	paste(round((map[4,2]+map[10,2])/sum(map[2:13,2]),4)*100,"%",sep=""),
	paste(round((map[5,2]+map[11,2])/sum(map[2:13,2]),4)*100,"%",sep=""),
	paste(round((map[6,2]+map[12,2])/sum(map[2:13,2]),4)*100,"%",sep=""),
	paste(round((map[7,2]+map[13,2])/sum(map[2:13,2]),4)*100,"%",sep=""))
pie(c(map[2:7,2]+map[8:13,2]),main="piechart of RepeatMask Reads Rate of AllReads",
	label=l,col=c("#66c2a5","#fc8d62","#8da0cb","#e78ac3","#a6d854","#ffd92f"))
legend("topleft",
	legend=c("LINE","SINE","LTR","DNA","Satellite","Simple_repeat"),
	col=c("#66c2a5","#fc8d62","#8da0cb","#e78ac3","#a6d854","#ffd92f"),pch=15,cex=0.9)
text(0.8,-1,label=paste(round(sum(map[2:13,2])/map[1,2]*100,2),"% reads in RepeatMask",sep=""),col="red")

l=c(paste(round((map[2,3]+map[8,3])/sum(map[2:13,3]),4)*100,"%",sep=""),
	paste(round((map[3,3]+map[9,3])/sum(map[2:13,3]),4)*100,"%",sep=""),
	paste(round((map[4,3]+map[10,3])/sum(map[2:13,3]),4)*100,"%",sep=""),
	paste(round((map[5,3]+map[11,3])/sum(map[2:13,3]),4)*100,"%",sep=""),
	paste(round((map[6,3]+map[12,3])/sum(map[2:13,3]),4)*100,"%",sep=""),
	paste(round((map[7,3]+map[13,3])/sum(map[2:13,3]),4)*100,"%",sep=""))
pie(c(map[2:7,3]+map[8:13,3]),main="piechart of RepeatMask Reads Rate of UniqMappedReads",
	label=l,col=c("#66c2a5","#fc8d62","#8da0cb","#e78ac3","#a6d854","#ffd92f"))
legend("topleft",
	legend=c("LINE","SINE","LTR","DNA","Satellite","Simple_repeat"),
	col=c("#66c2a5","#fc8d62","#8da0cb","#e78ac3","#a6d854","#ffd92f"),pch=15,cex=0.9)
text(0.8,-1,label=paste(round(sum(map[2:13,3])/map[1,3]*100,2),"% reads in RepeatMask",sep=""),col="red")

l=c(paste(round((map[2,4]+map[8,4])/sum(map[2:13,4]),4)*100,"%",sep=""),
	paste(round((map[3,4]+map[9,4])/sum(map[2:13,4]),4)*100,"%",sep=""),
	paste(round((map[4,4]+map[10,4])/sum(map[2:13,4]),4)*100,"%",sep=""),
	paste(round((map[5,4]+map[11,4])/sum(map[2:13,4]),4)*100,"%",sep=""),
	paste(round((map[6,4]+map[12,4])/sum(map[2:13,4]),4)*100,"%",sep=""),
	paste(round((map[7,4]+map[13,4])/sum(map[2:13,4]),4)*100,"%",sep=""))
pie(c(map[2:7,4]+map[8:13,4]),main="piechart of RepeatMask Reads Rate of MultipleMappedReads",
	label=l,col=c("#66c2a5","#fc8d62","#8da0cb","#e78ac3","#a6d854","#ffd92f"))
legend("topleft",
	legend=c("LINE","SINE","LTR","DNA","Satellite","Simple_repeat"),
	col=c("#66c2a5","#fc8d62","#8da0cb","#e78ac3","#a6d854","#ffd92f"),pch=15,cex=0.9)
text(0.8,-1,label=paste(round(sum(map[2:13,4])/map[1,4]*100,2),"% reads in RepeatMask",sep=""),col="red")

###RepeatMask distribution
all_LINE=c(map[2,2],map[8,2])/sum(map[2,2],map[8,2])*100
all_SINE=c(map[3,2],map[9,2])/sum(map[3,2],map[9,2])*100
all_LTR=c(map[4,2],map[10,2])/sum(map[4,2],map[10,2])*100
uniq_LINE=c(map[2,3],map[8,3])/sum(map[2,3],map[8,3])*100
uniq_SINE=c(map[3,3],map[9,3])/sum(map[3,3],map[9,3])*100
uniq_LTR=c(map[4,3],map[10,3])/sum(map[4,3],map[10,3])*100
multi_LINE=c(map[2,4],map[8,4])/sum(map[2,4],map[8,4])*100
multi_SINE=c(map[3,4],map[9,4])/sum(map[3,4],map[9,4])*100
multi_LTR=c(map[4,4],map[10,4])/sum(map[4,4],map[10,4])*100
x=cbind(all_LINE,all_SINE,all_LTR,uniq_LINE,uniq_SINE,uniq_LTR,multi_LINE,multi_SINE,multi_LTR)

barplot(x,beside=TRUE,legend=c("sense","anti-sense"),args.legend=list(x="topright"),ylim=c(0,100),
	ylab="percent (%)",col=c("#d73027","#4575b4"),main="orientation of piRNA from in RepeatMask",
	cex.names=0.7)

###piLoci mapping rate
l=c(paste(round(map[14,2]/map[1,2],4)*100,"%",sep=""),
	paste(round(map[15,2]/map[1,2],4)*100,"%",sep=""),
	paste(round(map[16,2]/map[1,2],4)*100,"%",sep=""),
	paste(round((map[1,2]-map[14,2]-map[15,2]-map[16,2])/map[1,2],4)*100,"%",sep=""))
pie(c(map[14,2],map[15,2],map[16,2],map[1,2]-map[14,2]-map[15,2]-map[16,2]),
	main="piechart of piRNA Loci Reads Rate of AllMappedReads",
	label=l,col=c("#1f78b4","#e31a1c","#4daf4a","#984ea3"))
legend("topleft",
	legend=c("Pachytene","Prepachytene","Hybrid","Others"),
	col=c("#1f78b4","#e31a1c","#4daf4a","#984ea3"),pch=15,cex=0.9)

l=c(paste(round(map[14,3]/map[1,3],4)*100,"%",sep=""),
	paste(round(map[15,3]/map[1,3],4)*100,"%",sep=""),
	paste(round(map[16,3]/map[1,3],4)*100,"%",sep=""),
	paste(round((map[1,3]-map[14,3]-map[15,3]-map[16,3])/map[1,3],4)*100,"%",sep=""))
pie(c(map[14,3],map[15,3],map[16,3],map[1,3]-map[14,3]-map[15,3]-map[16,3]),
	main="piechart of piRNA Loci Reads Rate of UniqueMappedReads",
	label=l,col=c("#1f78b4","#e31a1c","#4daf4a","#984ea3"))
legend("topleft",
	legend=c("Pachytene","Prepachytene","Hybrid","Others"),
	col=c("#1f78b4","#e31a1c","#4daf4a","#984ea3"),pch=15,cex=0.9)

l=c(paste(round(map[14,4]/map[1,4],4)*100,"%",sep=""),
	paste(round(map[15,4]/map[1,4],4)*100,"%",sep=""),
	paste(round(map[16,4]/map[1,4],4)*100,"%",sep=""),
	paste(round((map[1,4]-map[14,4]-map[15,4]-map[16,4])/map[1,4],4)*100,"%",sep=""))
pie(c(map[14,4],map[15,4],map[16,4],map[1,4]-map[14,4]-map[15,4]-map[16,4]),
	main="piechart of piRNA Loci Reads Rate of MultipleMappedReads",
	label=l,col=c("#1f78b4","#e31a1c","#4daf4a","#984ea3"))
legend("topleft",
	legend=c("Pachytene","Prepachytene","Hybrid","Others"),
	col=c("#1f78b4","#e31a1c","#4daf4a","#984ea3"),pch=15,cex=0.9)

###piLoci mapping rate of RepeatMask Reads
#anti-sense LINE
l=c(paste(round(map[17,2]/map[8,2],4)*100,"%",sep=""),
	paste(round(map[23,2]/map[8,2],4)*100,"%",sep=""),
	paste(round(map[29,2]/map[8,2],4)*100,"%",sep=""),
	paste(round(map[35,2]/map[8,2],4)*100,"%",sep=""),
	paste(round(map[41,2]/map[8,2],4)*100,"%",sep=""),
	paste(round(map[47,2]/map[8,2],4)*100,"%",sep=""),
	paste(round((map[8,2]-map[17,2]-map[23,2]-map[29,2]-map[35,2]-map[41,2]-map[47,2])/map[8,2],4)*100,"%",sep=""))

pie(c(map[17,2],map[23,2],map[29,2],map[35,2],map[41,2],map[47,2],map[8,2]-map[17,2]-map[23,2]-map[29,2]-map[35,2]-map[41,2]-map[47,2]),
	main="piechar of piRNA Loci Reads Rate of Anti-senseLINEAllMappedReads",
	label=l,col=c("#1f78b4","#a6cee3","#33a02c","#b2df8a","#e31a1c","#fb9a99","#984ea3"))
legend("topleft",
	legend=c("Pachytene","Anti-Pachytene","Prepachytene","Anti-Prepachytene","Hybrid","Anti-Hybrid","Others"),
	col=c("#1f78b4","#a6cee3","#33a02c","#b2df8a","#e31a1c","#fb9a99","#984ea3"),
	pch=15,cex=0.8)

l=c(paste(round(map[17,3]/map[8,3],4)*100,"%",sep=""),
    paste(round(map[23,3]/map[8,3],4)*100,"%",sep=""),
    paste(round(map[29,3]/map[8,3],4)*100,"%",sep=""),
    paste(round(map[35,3]/map[8,3],4)*100,"%",sep=""),
    paste(round(map[41,3]/map[8,3],4)*100,"%",sep=""),
    paste(round(map[47,3]/map[8,3],4)*100,"%",sep=""),
    paste(round((map[8,3]-map[17,3]-map[23,3]-map[29,3]-map[35,3]-map[41,3]-map[47,3])/map[8,3],4)*100,"%",sep=""))

pie(c(map[17,3],map[23,3],map[29,3],map[35,3],map[41,3],map[47,3],map[8,3]-map[17,3]-map[23,3]-map[29,3]-map[35,3]-map[41,3]-map[47,3]),
    main="piechar of piRNA Loci Reads Rate of Anti-senseLINEUniqueMappedReads",
    label=l,col=c("#1f78b4","#a6cee3","#33a02c","#b2df8a","#e31a1c","#fb9a99","#984ea3"))
legend("topleft",
	legend=c("Pachytene","Anti-Pachytene","Prepachytene","Anti-Prepachytene","Hybrid","Anti-Hybrid","Others"),
	col=c("#1f78b4","#a6cee3","#33a02c","#b2df8a","#e31a1c","#fb9a99","#984ea3"),
	pch=15,cex=0.8)

l=c(paste(round(map[17,4]/map[8,4],4)*100,"%",sep=""),
    paste(round(map[23,4]/map[8,4],4)*100,"%",sep=""),
    paste(round(map[29,4]/map[8,4],4)*100,"%",sep=""),
    paste(round(map[35,4]/map[8,4],4)*100,"%",sep=""),
    paste(round(map[41,4]/map[8,4],4)*100,"%",sep=""),
    paste(round(map[47,4]/map[8,4],4)*100,"%",sep=""),
    paste(round((map[8,4]-map[17,4]-map[23,4]-map[29,4]-map[35,4]-map[41,4]-map[47,4])/map[8,4],4)*100,"%",sep=""))

pie(c(map[17,4],map[23,4],map[29,4],map[35,4],map[41,4],map[47,4],map[8,4]-map[17,4]-map[23,4]-map[29,4]-map[35,4]-map[41,4]-map[47,4]),
    main="piechar of piRNA Loci Reads Rate of Anti-senseLINEMultipleMappedReads",
    label=l,col=c("#1f78b4","#a6cee3","#33a02c","#b2df8a","#e31a1c","#fb9a99","#984ea3"))
legend("topleft",
	legend=c("Pachytene","Anti-Pachytene","Prepachytene","Anti-Prepachytene","Hybrid","Anti-Hybrid","Others"),
	col=c("#1f78b4","#a6cee3","#33a02c","#b2df8a","#e31a1c","#fb9a99","#984ea3"),
	pch=15,cex=0.8)

#anti-sense SINE
l=c(paste(round(map[18,2]/map[9,2],4)*100,"%",sep=""),
	paste(round(map[24,2]/map[9,2],4)*100,"%",sep=""),
	paste(round(map[30,2]/map[9,2],4)*100,"%",sep=""),
	paste(round(map[36,2]/map[9,2],4)*100,"%",sep=""),
	paste(round(map[42,2]/map[9,2],4)*100,"%",sep=""),
	paste(round(map[48,2]/map[9,2],4)*100,"%",sep=""),
	paste(round((map[9,2]-map[18,2]-map[24,2]-map[30,2]-map[36,2]-map[42,2]-map[48,2])/map[9,2],4)*100,"%",sep=""))

pie(c(map[18,2],map[24,2],map[30,2],map[36,2],map[42,2],map[48,2],map[9,2]-map[18,2]-map[24,2]-map[30,2]-map[36,2]-map[42,2]-map[48,2]),
	main="piechar of piRNA Loci Reads Rate of Anti-senseSINEAllMappedReads",
	label=l,col=c("#1f78b4","#a6cee3","#33a02c","#b2df8a","#e31a1c","#fb9a99","#984ea3"))
legend("topleft",
	legend=c("Pachytene","Anti-Pachytene","Prepachytene","Anti-Prepachytene","Hybrid","Anti-Hybrid","Others"),
	col=c("#1f78b4","#a6cee3","#33a02c","#b2df8a","#e31a1c","#fb9a99","#984ea3"),
	pch=15,cex=0.8)

l=c(paste(round(map[18,3]/map[9,3],4)*100,"%",sep=""),
    paste(round(map[24,3]/map[9,3],4)*100,"%",sep=""),
    paste(round(map[30,3]/map[9,3],4)*100,"%",sep=""),
    paste(round(map[36,3]/map[9,3],4)*100,"%",sep=""),
    paste(round(map[42,3]/map[9,3],4)*100,"%",sep=""),
    paste(round(map[48,3]/map[9,3],4)*100,"%",sep=""),
    paste(round((map[9,3]-map[18,3]-map[24,3]-map[30,3]-map[36,3]-map[42,3]-map[48,3])/map[9,3],4)*100,"%",sep=""))

pie(c(map[18,3],map[24,3],map[30,3],map[36,3],map[42,3],map[48,3],map[9,3]-map[18,3]-map[24,3]-map[30,3]-map[36,3]-map[42,3]-map[48,3]),
    main="piechar of piRNA Loci Reads Rate of Anti-senseSINEUniqMappedReads",
    label=l,col=c("#1f78b4","#a6cee3","#33a02c","#b2df8a","#e31a1c","#fb9a99","#984ea3"))
legend("topleft",
	legend=c("Pachytene","Anti-Pachytene","Prepachytene","Anti-Prepachytene","Hybrid","Anti-Hybrid","Others"),
	col=c("#1f78b4","#a6cee3","#33a02c","#b2df8a","#e31a1c","#fb9a99","#984ea3"),
	pch=15,cex=0.8)

l=c(paste(round(map[18,4]/map[9,4],4)*100,"%",sep=""),
    paste(round(map[24,4]/map[9,4],4)*100,"%",sep=""),
    paste(round(map[30,4]/map[9,4],4)*100,"%",sep=""),
    paste(round(map[36,4]/map[9,4],4)*100,"%",sep=""),
    paste(round(map[42,4]/map[9,4],4)*100,"%",sep=""),
    paste(round(map[48,4]/map[9,4],4)*100,"%",sep=""),
    paste(round((map[9,4]-map[18,4]-map[24,4]-map[30,4]-map[36,4]-map[42,4]-map[48,4])/map[9,4],4)*100,"%",sep=""))

pie(c(map[18,4],map[24,4],map[30,4],map[36,4],map[42,4],map[48,4],map[9,4]-map[18,4]-map[24,4]-map[30,4]-map[36,4]-map[42,4]-map[48,4]),
    main="piechar of piRNA Loci Reads Rate of Anti-senseSINEMultipleMappedReads",
    label=l,col=c("#1f78b4","#a6cee3","#33a02c","#b2df8a","#e31a1c","#fb9a99","#984ea3"))
legend("topleft",
	legend=c("Pachytene","Anti-Pachytene","Prepachytene","Anti-Prepachytene","Hybrid","Anti-Hybrid","Others"),
	col=c("#1f78b4","#a6cee3","#33a02c","#b2df8a","#e31a1c","#fb9a99","#984ea3"),
	pch=15,cex=0.8)

#anti-sense LTR
l=c(paste(round(map[19,2]/map[10,2],4)*100,"%",sep=""),
    paste(round(map[25,2]/map[10,2],4)*100,"%",sep=""),
    paste(round(map[31,2]/map[10,2],4)*100,"%",sep=""),
    paste(round(map[37,2]/map[10,2],4)*100,"%",sep=""),
    paste(round(map[43,2]/map[10,2],4)*100,"%",sep=""),
    paste(round(map[49,2]/map[10,2],4)*100,"%",sep=""),
    paste(round((map[10,2]-map[19,2]-map[25,2]-map[31,2]-map[37,2]-map[43,2]-map[49,2])/map[10,2],4)*100,"%",sep=""))

pie(c(map[19,2],map[25,2],map[31,2],map[37,2],map[43,2],map[49,2],map[10,2]-map[19,2]-map[25,2]-map[31,2]-map[37,2]-map[43,2]-map[49,2]),
    main="piechar of piRNA Loci Reads Rate of Anti-senseLTRAllMappedReads",
    label=l,col=c("#1f78b4","#a6cee3","#33a02c","#b2df8a","#e31a1c","#fb9a99","#984ea3"))
legend("topleft",
	legend=c("Pachytene","Anti-Pachytene","Prepachytene","Anti-Prepachytene","Hybrid","Anti-Hybrid","Others"),
	col=c("#1f78b4","#a6cee3","#33a02c","#b2df8a","#e31a1c","#fb9a99","#984ea3"),
	pch=15,cex=0.8)

l=c(paste(round(map[19,3]/map[10,3],4)*100,"%",sep=""),
    paste(round(map[25,3]/map[10,3],4)*100,"%",sep=""),
    paste(round(map[31,3]/map[10,3],4)*100,"%",sep=""),
    paste(round(map[37,3]/map[10,3],4)*100,"%",sep=""),
    paste(round(map[43,3]/map[10,3],4)*100,"%",sep=""),
    paste(round(map[49,3]/map[10,3],4)*100,"%",sep=""),
    paste(round((map[10,3]-map[19,3]-map[25,3]-map[31,3]-map[37,3]-map[43,3]-map[49,3])/map[10,3],4)*100,"%",sep=""))

pie(c(map[19,3],map[25,3],map[31,3],map[37,3],map[43,3],map[49,3],map[10,3]-map[19,3]-map[25,3]-map[31,3]-map[37,3]-map[43,3]-map[49,3]),
    main="piechar of piRNA Loci Reads Rate of Anti-senseLTRUniqMappedReads",
    label=l,col=c("#1f78b4","#a6cee3","#33a02c","#b2df8a","#e31a1c","#fb9a99","#984ea3"))
legend("topleft",
	legend=c("Pachytene","Anti-Pachytene","Prepachytene","Anti-Prepachytene","Hybrid","Anti-Hybrid","Others"),
	col=c("#1f78b4","#a6cee3","#33a02c","#b2df8a","#e31a1c","#fb9a99","#984ea3"),
	pch=15,cex=0.8)

l=c(paste(round(map[19,4]/map[10,4],4)*100,"%",sep=""),
    paste(round(map[25,4]/map[10,4],4)*100,"%",sep=""),
    paste(round(map[31,4]/map[10,4],4)*100,"%",sep=""),
    paste(round(map[37,4]/map[10,4],4)*100,"%",sep=""),
    paste(round(map[43,4]/map[10,4],4)*100,"%",sep=""),
    paste(round(map[49,4]/map[10,4],4)*100,"%",sep=""),
    paste(round((map[10,4]-map[19,4]-map[25,4]-map[31,4]-map[37,4]-map[43,4]-map[49,4])/map[10,4],4)*100,"%",sep=""))

pie(c(map[19,4],map[25,4],map[31,4],map[37,4],map[43,4],map[49,4],map[10,4]-map[19,4]-map[25,4]-map[31,4]-map[37,4]-map[43,4]-map[49,4]),
    main="piechar of piRNA Loci Reads Rate of Anti-senseLTRMultipleMappedReads",
    label=l,col=c("#1f78b4","#a6cee3","#33a02c","#b2df8a","#e31a1c","#fb9a99","#984ea3"))
legend("topleft",
	legend=c("Pachytene","Anti-Pachytene","Prepachytene","Anti-Prepachytene","Hybrid","Anti-Hybrid","Others"),
	col=c("#1f78b4","#a6cee3","#33a02c","#b2df8a","#e31a1c","#fb9a99","#984ea3"),
	pch=15,cex=0.8)

#sense LINE
l=c(paste(round(map[20,2]/map[2,2],4)*100,"%",sep=""),
    paste(round(map[26,2]/map[2,2],4)*100,"%",sep=""),
    paste(round(map[32,2]/map[2,2],4)*100,"%",sep=""),
    paste(round(map[38,2]/map[2,2],4)*100,"%",sep=""),
    paste(round(map[44,2]/map[2,2],4)*100,"%",sep=""),
    paste(round(map[50,2]/map[2,2],4)*100,"%",sep=""),
    paste(round((map[2,2]-map[20,2]-map[26,2]-map[32,2]-map[38,2]-map[44,2]-map[50,2])/map[2,2],4)*100,"%",sep=""))

pie(c(map[20,2],map[26,2],map[32,2],map[38,2],map[44,2],map[50,2],map[2,2]-map[20,2]-map[26,2]-map[32,2]-map[38,2]-map[44,2]-map[50,2]),
    main="piechar of piRNA Loci Reads Rate of SenseLINEAllMappedReads",
    label=l,col=c("#1f78b4","#a6cee3","#33a02c","#b2df8a","#e31a1c","#fb9a99","#984ea3"))
legend("topleft",
	legend=c("Pachytene","Anti-Pachytene","Prepachytene","Anti-Prepachytene","Hybrid","Anti-Hybrid","Others"),
	col=c("#1f78b4","#a6cee3","#33a02c","#b2df8a","#e31a1c","#fb9a99","#984ea3"),
	pch=15,cex=0.8)

l=c(paste(round(map[20,3]/map[2,3],4)*100,"%",sep=""),
    paste(round(map[26,3]/map[2,3],4)*100,"%",sep=""),
    paste(round(map[32,3]/map[2,3],4)*100,"%",sep=""),
    paste(round(map[38,3]/map[2,3],4)*100,"%",sep=""),
    paste(round(map[44,3]/map[2,3],4)*100,"%",sep=""),
    paste(round(map[50,3]/map[2,3],4)*100,"%",sep=""),
    paste(round((map[2,3]-map[20,3]-map[26,3]-map[32,3]-map[38,3]-map[44,3]-map[50,3])/map[2,3],4)*100,"%",sep=""))

pie(c(map[20,3],map[26,3],map[32,3],map[38,3],map[44,3],map[50,3],map[2,3]-map[20,3]-map[26,3]-map[32,3]-map[38,3]-map[44,3]-map[50,3]),
    main="piechar of piRNA Loci Reads Rate of SenseLINEUniqMappedReads",
    label=l,col=c("#1f78b4","#a6cee3","#33a02c","#b2df8a","#e31a1c","#fb9a99","#984ea3"))
legend("topleft",
	legend=c("Pachytene","Anti-Pachytene","Prepachytene","Anti-Prepachytene","Hybrid","Anti-Hybrid","Others"),
	col=c("#1f78b4","#a6cee3","#33a02c","#b2df8a","#e31a1c","#fb9a99","#984ea3"),
	pch=15,cex=0.8)

l=c(paste(round(map[20,4]/map[2,4],4)*100,"%",sep=""),
    paste(round(map[26,4]/map[2,4],4)*100,"%",sep=""),
    paste(round(map[32,4]/map[2,4],4)*100,"%",sep=""),
    paste(round(map[38,4]/map[2,4],4)*100,"%",sep=""),
    paste(round(map[44,4]/map[2,4],4)*100,"%",sep=""),
    paste(round(map[50,4]/map[2,4],4)*100,"%",sep=""),
    paste(round((map[2,4]-map[20,4]-map[26,4]-map[32,4]-map[38,4]-map[44,4]-map[50,4])/map[2,4],4)*100,"%",sep=""))

pie(c(map[20,4],map[26,4],map[32,4],map[38,4],map[44,4],map[50,4],map[2,4]-map[20,4]-map[26,4]-map[32,4]-map[38,4]-map[44,4]-map[50,4]),
    main="piechar of piRNA Loci Reads Rate of SenseLINEMultipleMappedReads",
    label=l,col=c("#1f78b4","#a6cee3","#33a02c","#b2df8a","#e31a1c","#fb9a99","#984ea3"))
legend("topleft",
	legend=c("Pachytene","Anti-Pachytene","Prepachytene","Anti-Prepachytene","Hybrid","Anti-Hybrid","Others"),
	col=c("#1f78b4","#a6cee3","#33a02c","#b2df8a","#e31a1c","#fb9a99","#984ea3"),
	pch=15,cex=0.8)

#sense SINE
l=c(paste(round(map[21,2]/map[3,2],4)*100,"%",sep=""),
    paste(round(map[27,2]/map[3,2],4)*100,"%",sep=""),
    paste(round(map[33,2]/map[3,2],4)*100,"%",sep=""),
    paste(round(map[39,2]/map[3,2],4)*100,"%",sep=""),
    paste(round(map[45,2]/map[3,2],4)*100,"%",sep=""),
    paste(round(map[51,2]/map[3,2],4)*100,"%",sep=""),
    paste(round((map[3,2]-map[21,2]-map[27,2]-map[33,2]-map[39,2]-map[45,2]-map[51,2])/map[3,2],4)*100,"%",sep=""))

pie(c(map[21,2],map[27,2],map[33,2],map[39,2],map[45,2],map[51,2],map[3,2]-map[21,2]-map[27,2]-map[33,2]-map[39,2]-map[45,2]-map[51,2]),
    main="piechar of piRNA Loci Reads Rate of SenseSINEAllMappedReads",
    label=l,col=c("#1f78b4","#a6cee3","#33a02c","#b2df8a","#e31a1c","#fb9a99","#984ea3"))
legend("topleft",
	legend=c("Pachytene","Anti-Pachytene","Prepachytene","Anti-Prepachytene","Hybrid","Anti-Hybrid","Others"),
	col=c("#1f78b4","#a6cee3","#33a02c","#b2df8a","#e31a1c","#fb9a99","#984ea3"),
	pch=15,cex=0.8)

l=c(paste(round(map[21,3]/map[3,3],4)*100,"%",sep=""),
    paste(round(map[27,3]/map[3,3],4)*100,"%",sep=""),
    paste(round(map[33,3]/map[3,3],4)*100,"%",sep=""),
    paste(round(map[39,3]/map[3,3],4)*100,"%",sep=""),
    paste(round(map[45,3]/map[3,3],4)*100,"%",sep=""),
    paste(round(map[51,3]/map[3,3],4)*100,"%",sep=""),
    paste(round((map[3,3]-map[21,3]-map[27,3]-map[33,3]-map[39,3]-map[45,3]-map[51,3])/map[3,3],4)*100,"%",sep=""))

pie(c(map[21,3],map[27,3],map[33,3],map[39,3],map[45,3],map[51,3],map[3,3]-map[21,3]-map[27,3]-map[33,3]-map[39,3]-map[45,3]-map[51,3]),
    main="piechar of piRNA Loci Reads Rate of SenseSINEUniqMappedReads",
    label=l,col=c("#1f78b4","#a6cee3","#33a02c","#b2df8a","#e31a1c","#fb9a99","#984ea3"))
legend("topleft",
	legend=c("Pachytene","Anti-Pachytene","Prepachytene","Anti-Prepachytene","Hybrid","Anti-Hybrid","Others"),
	col=c("#1f78b4","#a6cee3","#33a02c","#b2df8a","#e31a1c","#fb9a99","#984ea3"),
	pch=15,cex=0.8)

l=c(paste(round(map[21,4]/map[3,4],4)*100,"%",sep=""),
    paste(round(map[27,4]/map[3,4],4)*100,"%",sep=""),
    paste(round(map[33,4]/map[3,4],4)*100,"%",sep=""),
    paste(round(map[39,4]/map[3,4],4)*100,"%",sep=""),
    paste(round(map[45,4]/map[3,4],4)*100,"%",sep=""),
    paste(round(map[51,4]/map[3,4],4)*100,"%",sep=""),
    paste(round((map[3,4]-map[21,4]-map[27,4]-map[33,4]-map[39,4]-map[45,4]-map[51,4])/map[3,4],4)*100,"%",sep=""))

pie(c(map[21,4],map[27,4],map[33,4],map[39,4],map[45,4],map[51,4],map[3,4]-map[21,4]-map[27,4]-map[33,4]-map[39,4]-map[45,4]-map[51,4]),
    main="piechar of piRNA Loci Reads Rate of SenseSINEMultipleMappedReads",
    label=l,col=c("#1f78b4","#a6cee3","#33a02c","#b2df8a","#e31a1c","#fb9a99","#984ea3"))
legend("topleft",
	legend=c("Pachytene","Anti-Pachytene","Prepachytene","Anti-Prepachytene","Hybrid","Anti-Hybrid","Others"),
	col=c("#1f78b4","#a6cee3","#33a02c","#b2df8a","#e31a1c","#fb9a99","#984ea3"),
	pch=15,cex=0.8)

#sense LTR
l=c(paste(round(map[22,2]/map[4,2],4)*100,"%",sep=""),
    paste(round(map[28,2]/map[4,2],4)*100,"%",sep=""),
    paste(round(map[34,2]/map[4,2],4)*100,"%",sep=""),
    paste(round(map[40,2]/map[4,2],4)*100,"%",sep=""),
    paste(round(map[46,2]/map[4,2],4)*100,"%",sep=""),
    paste(round(map[52,2]/map[4,2],4)*100,"%",sep=""),
    paste(round((map[4,2]-map[22,2]-map[28,2]-map[34,2]-map[40,2]-map[46,2]-map[52,2])/map[4,2],4)*100,"%",sep=""))

pie(c(map[22,2],map[28,2],map[34,2],map[40,2],map[46,2],map[52,2],map[4,2]-map[22,2]-map[28,2]-map[34,2]-map[40,2]-map[46,2]-map[52,2]),
    main="piechar of piRNA Loci Reads Rate of SenseLTRAllMappedReads",
    label=l,col=c("#1f78b4","#a6cee3","#33a02c","#b2df8a","#e31a1c","#fb9a99","#984ea3"))
legend("topleft",
	legend=c("Pachytene","Anti-Pachytene","Prepachytene","Anti-Prepachytene","Hybrid","Anti-Hybrid","Others"),
	col=c("#1f78b4","#a6cee3","#33a02c","#b2df8a","#e31a1c","#fb9a99","#984ea3"),
	pch=15,cex=0.8)

l=c(paste(round(map[22,3]/map[4,3],4)*100,"%",sep=""),
    paste(round(map[28,3]/map[4,3],4)*100,"%",sep=""),
    paste(round(map[34,3]/map[4,3],4)*100,"%",sep=""),
    paste(round(map[40,3]/map[4,3],4)*100,"%",sep=""),
    paste(round(map[46,3]/map[4,3],4)*100,"%",sep=""),
    paste(round(map[52,3]/map[4,3],4)*100,"%",sep=""),
    paste(round((map[4,3]-map[22,3]-map[28,3]-map[34,3]-map[40,3]-map[46,3]-map[52,3])/map[4,3],4)*100,"%",sep=""))

pie(c(map[22,3],map[28,3],map[34,3],map[40,3],map[46,3],map[52,3],map[4,3]-map[22,3]-map[28,3]-map[34,3]-map[40,3]-map[46,3]-map[52,3]),
    main="piechar of piRNA Loci Reads Rate of SenseLTRUniqMappedReads",
    label=l,col=c("#1f78b4","#a6cee3","#33a02c","#b2df8a","#e31a1c","#fb9a99","#984ea3"))
legend("topleft",
	legend=c("Pachytene","Anti-Pachytene","Prepachytene","Anti-Prepachytene","Hybrid","Anti-Hybrid","Others"),
	col=c("#1f78b4","#a6cee3","#33a02c","#b2df8a","#e31a1c","#fb9a99","#984ea3"),
	pch=15,cex=0.8)

l=c(paste(round(map[22,4]/map[4,4],4)*100,"%",sep=""),
    paste(round(map[28,4]/map[4,4],4)*100,"%",sep=""),
    paste(round(map[34,4]/map[4,4],4)*100,"%",sep=""),
    paste(round(map[40,4]/map[4,4],4)*100,"%",sep=""),
    paste(round(map[46,4]/map[4,4],4)*100,"%",sep=""),
    paste(round(map[52,4]/map[4,4],4)*100,"%",sep=""),
    paste(round((map[4,4]-map[22,4]-map[28,4]-map[34,4]-map[40,4]-map[46,4]-map[52,4])/map[4,4],4)*100,"%",sep=""))

pie(c(map[22,4],map[28,4],map[34,4],map[40,4],map[46,4],map[52,4],map[4,4]-map[22,4]-map[28,4]-map[34,4]-map[40,4]-map[46,4]-map[52,4]),
    main="piechar of piRNA Loci Reads Rate of SenseLTRMultipleMappedReads",
    label=l,col=c("#1f78b4","#a6cee3","#33a02c","#b2df8a","#e31a1c","#fb9a99","#984ea3"))
legend("topleft",
	legend=c("Pachytene","Anti-Pachytene","Prepachytene","Anti-Prepachytene","Hybrid","Anti-Hybrid","Others"),
	col=c("#1f78b4","#a6cee3","#33a02c","#b2df8a","#e31a1c","#fb9a99","#984ea3"),
	pch=15,cex=0.8)

#length distribution
plot(24:32,as.vector(all_len[,1]/sum(all_len[,1])),type="l",xlab="length (bp)",ylab="percentage (%)",
	main="length distribution of each type of reads",col="#1f78b4",lwd=2,ylim=c(0,0.35))
lines(24:32,as.vector(uniq_len[,1]/sum(uniq_len[,1])),col="#e31a1c",lwd=2)
lines(24:32,as.vector(multi_len[,1]/sum(multi_len[,1])),col="#33a02c",lwd=2)
legend("topleft",legend=c("All Mapped Reads","Unique Mapped Reads","Multiple Mapped Reads"),
	col=c("#1f78b4","#e31a1c","#33a02c"),cex=0.9,lwd=3)

#strand distribution
all=all_str[,1]/sum(all_str[,1])*100
uniq=uniq_str[,1]/sum(uniq_str[,1])*100
multi=multi_str[,1]/sum(multi_str[,1])*100
x=cbind(all,uniq,multi)

barplot(x,beside=TRUE,legend=c("+","-"),args.legend=list(x="topright"),ylim=c(0,80),
	ylab="percent (%)",col=c("#d73027","#4575b4"),main="strand distribution of each mapping category")

###close pdf
dev.off()


