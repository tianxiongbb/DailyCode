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
#all mapping rate
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
uniq_LTR=c(map[4,3],map[10,3])/sum(map[4,3],map[10,4])*100
multi_LINE=c(map[2,4],map[8,4])/sum(map[2,4],map[8,4])*100
multi_SINE=c(map[3,4],map[9,4])/sum(map[3,4],map[9,4])*100
multi_LTR=c(map[4,4],map[10,4])/sum(map[4,4],map[10,4])*100
x=cbind(all_LINE,all_SINE,all_LTR,uniq_LINE,uniq_SINE,uniq_LTR,multi_LINE,multi_SINE,multi_LTR)

barplot(x,beside=TRUE,legend=c("sense","anti-sense"),args.legend=list(x="topright"),ylim=c(0,100),
	ylab="percent (%)",col=c("#d73027","#4575b4"),main="orientation of piRNA from in RepeatMask",
	cex.names=0.7)

###piLoci mapping rate
l=c(paste(round(map[14,2]/map[1,2],4)*100,"%",sep=""),
	paste(round((map[1,2]-map[14,2])/map[1,2],4)*100,"%",sep=""))
pie(c(map[14,2],map[1,2]-map[14,2]),main="piechart of piRNA Loci Reads Rate of AllReads",
	label=l,col=c("#1f78b4","#e31a1c"))
legend("topleft",
	legend=c("Mapped","Unmapped"),col=c("#1f78b4","#e31a1c"),pch=15,cex=0.9)

l=c(paste(round(map[14,3]/map[1,3],4)*100,"%",sep=""),
	paste(round((map[1,3]-map[14,3])/map[1,3],4)*100,"%",sep=""))
pie(c(map[14,3],map[1,3]-map[14,3]),main="piechart of piRNA Loci Reads Rate of UniqMappedReads",
	label=l,col=c("#1f78b4","#e31a1c"))
legend("topleft",
	legend=c("Mapped","Unmapped"),col=c("#1f78b4","#e31a1c"),pch=15,cex=0.9)

l=c(paste(round(map[14,4]/map[1,4],4)*100,"%",sep=""),
	paste(round((map[1,4]-map[14,4])/map[1,4],4)*100,"%",sep=""))
pie(c(map[14,4],map[1,4]-map[14,4]),main="piechart of piRNA Loci Reads Rate of MultipleMappedReads",
	label=l,col=c("#1f78b4","#e31a1c"))
legend("topleft",
	legend=c("Mapped","Unmapped"),col=c("#1f78b4","#e31a1c"),pch=15,cex=0.9)
###piLoci mapping rate of RepeatMask Reads
#anti-sense LINE
a=map[15,2]
b=map[21,2]
c=map[8,2]
l=c(paste(round(a/c,4)*100,"%",sep=""),
	paste(round(b/c,4)*100,"%",sep=""),
	paste(round((c-a-b)/c,4)*100,"%",sep=""))

pie(c(a,b,c-a-b),
	main="piechar of piRNA Loci Reads Rate of Anti-senseLINEAllMappedReads",
	label=l,col=c("#1f78b4","#a6cee3","#984ea3"))
legend("topleft",
	legend=c("PiLoci","Anti-PiLoci","Others"),
	col=c("#1f78b4","#a6cee3","#984ea3"),
	pch=15,cex=0.8)

a=map[15,3]
b=map[21,3]
c=map[8,3]
l=c(paste(round(a/c,4)*100,"%",sep=""),
	paste(round(b/c,4)*100,"%",sep=""),
	paste(round((c-a-b)/c,4)*100,"%",sep=""))

pie(c(a,b,c-a-b),
	main="piechar of piRNA Loci Reads Rate of Anti-senseLINEUniqueMappedReads",
	label=l,col=c("#1f78b4","#a6cee3","#984ea3"))
legend("topleft",
	legend=c("PiLoci","Anti-PiLoci","Others"),
	col=c("#1f78b4","#a6cee3","#984ea3"),
	pch=15,cex=0.8)

a=map[15,4]
b=map[21,4]
c=map[8,4]
l=c(paste(round(a/c,4)*100,"%",sep=""),
	paste(round(b/c,4)*100,"%",sep=""),
	paste(round((c-a-b)/c,4)*100,"%",sep=""))

pie(c(a,b,c-a-b),
	main="piechar of piRNA Loci Reads Rate of Anti-senseLINEMultipleMappedReads",
	label=l,col=c("#1f78b4","#a6cee3","#984ea3"))
legend("topleft",
	legend=c("PiLoci","Anti-PiLoci","Others"),
	col=c("#1f78b4","#a6cee3","#984ea3"),
	pch=15,cex=0.8)

#anti-sense SINE
a=map[16,2]
b=map[22,2]
c=map[9,2]
l=c(paste(round(a/c,4)*100,"%",sep=""),
	paste(round(b/c,4)*100,"%",sep=""),
	paste(round((c-a-b)/c,4)*100,"%",sep=""))
c(a,b,c-a-b)
pie(c(a,b,c-a-b),
	main="piechar of piRNA Loci Reads Rate of Anti-senseSINEAllMappedReads",
	label=l,col=c("#1f78b4","#a6cee3","#984ea3"))
legend("topleft",
	legend=c("PiLoci","Anti-PiLoci","Others"),
	col=c("#1f78b4","#a6cee3","#984ea3"),
	pch=15,cex=0.8)

a=map[16,3]
b=map[22,3]
c=map[9,3]
l=c(paste(round(a/c,4)*100,"%",sep=""),
	paste(round(b/c,4)*100,"%",sep=""),
	paste(round((c-a-b)/c,4)*100,"%",sep=""))

pie(c(a,b,c-a-b),
	main="piechar of piRNA Loci Reads Rate of Anti-senseSINEUniqueMappedReads",
	label=l,col=c("#1f78b4","#a6cee3","#984ea3"))
legend("topleft",
	legend=c("PiLoci","Anti-PiLoci","Others"),
	col=c("#1f78b4","#a6cee3","#984ea3"),
	pch=15,cex=0.8)

a=map[16,4]
b=map[22,4]
c=map[9,4]
l=c(paste(round(a/c,4)*100,"%",sep=""),
	paste(round(b/c,4)*100,"%",sep=""),
	paste(round((c-a-b)/c,4)*100,"%",sep=""))

pie(c(a,b,c-a-b),
	main="piechar of piRNA Loci Reads Rate of Anti-senseSINEMultipleMappedReads",
	label=l,col=c("#1f78b4","#a6cee3","#984ea3"))
legend("topleft",
	legend=c("PiLoci","Anti-PiLoci","Others"),
	col=c("#1f78b4","#a6cee3","#984ea3"),
	pch=15,cex=0.8)

#anti-sense LTR
a=map[17,2]
b=map[23,2]
c=map[10,2]
l=c(paste(round(a/c,4)*100,"%",sep=""),
	paste(round(b/c,4)*100,"%",sep=""),
	paste(round((c-a-b)/c,4)*100,"%",sep=""))

pie(c(a,b,c-a-b),
	main="piechar of piRNA Loci Reads Rate of Anti-senseLTRAllMappedReads",
	label=l,col=c("#1f78b4","#a6cee3","#984ea3"))
legend("topleft",
	legend=c("PiLoci","Anti-PiLoci","Others"),
	col=c("#1f78b4","#a6cee3","#984ea3"),
	pch=15,cex=0.8)

a=map[17,3]
b=map[23,3]
c=map[10,3]
l=c(paste(round(a/c,4)*100,"%",sep=""),
	paste(round(b/c,4)*100,"%",sep=""),
	paste(round((c-a-b)/c,4)*100,"%",sep=""))

pie(c(a,b,c-a-b),
	main="piechar of piRNA Loci Reads Rate of Anti-senseLTRUniqueMappedReads",
	label=l,col=c("#1f78b4","#a6cee3","#984ea3"))
legend("topleft",
	legend=c("PiLoci","Anti-PiLoci","Others"),
	col=c("#1f78b4","#a6cee3","#984ea3"),
	pch=15,cex=0.8)

a=map[17,4]
b=map[23,4]
c=map[10,4]
l=c(paste(round(a/c,4)*100,"%",sep=""),
	paste(round(b/c,4)*100,"%",sep=""),
	paste(round((c-a-b)/c,4)*100,"%",sep=""))

pie(c(a,b,c-a-b),
	main="piechar of piRNA Loci Reads Rate of Anti-senseLTRMultipleMappedReads",
	label=l,col=c("#1f78b4","#a6cee3","#984ea3"))
legend("topleft",
	legend=c("PiLoci","Anti-PiLoci","Others"),
	col=c("#1f78b4","#a6cee3","#984ea3"),
	pch=15,cex=0.8)

#sense LINE
a=map[18,2]
b=map[24,2]
c=map[2,2]
l=c(paste(round(a/c,4)*100,"%",sep=""),
	paste(round(b/c,4)*100,"%",sep=""),
	paste(round((c-a-b)/c,4)*100,"%",sep=""))

pie(c(a,b,c-a-b),
	main="piechar of piRNA Loci Reads Rate of SenseLINEAllMappedReads",
	label=l,col=c("#1f78b4","#a6cee3","#984ea3"))
legend("topleft",
	legend=c("PiLoci","Anti-PiLoci","Others"),
	col=c("#1f78b4","#a6cee3","#984ea3"),
	pch=15,cex=0.8)

a=map[18,3]
b=map[24,3]
c=map[2,3]
l=c(paste(round(a/c,4)*100,"%",sep=""),
	paste(round(b/c,4)*100,"%",sep=""),
	paste(round((c-a-b)/c,4)*100,"%",sep=""))

pie(c(a,b,c-a-b),
	main="piechar of piRNA Loci Reads Rate of SenseLINEUniqueMappedReads",
	label=l,col=c("#1f78b4","#a6cee3","#984ea3"))
legend("topleft",
	legend=c("PiLoci","Anti-PiLoci","Others"),
	col=c("#1f78b4","#a6cee3","#984ea3"),
	pch=15,cex=0.8)

a=map[18,4]
b=map[24,4]
c=map[2,4]
l=c(paste(round(a/c,4)*100,"%",sep=""),
	paste(round(b/c,4)*100,"%",sep=""),
	paste(round((c-a-b)/c,4)*100,"%",sep=""))

pie(c(a,b,c-a-b),
	main="piechar of piRNA Loci Reads Rate of SenseLINEMultipleMappedReads",
	label=l,col=c("#1f78b4","#a6cee3","#984ea3"))
legend("topleft",
	legend=c("PiLoci","Anti-PiLoci","Others"),
	col=c("#1f78b4","#a6cee3","#984ea3"),
	pch=15,cex=0.8)

#sense SINE
a=map[19,2]
b=map[25,2]
c=map[3,2]
l=c(paste(round(a/c,4)*100,"%",sep=""),
	paste(round(b/c,4)*100,"%",sep=""),
	paste(round((c-a-b)/c,4)*100,"%",sep=""))

pie(c(a,b,c-a-b),
	main="piechar of piRNA Loci Reads Rate of SenseSINEAllMappedReads",
	label=l,col=c("#1f78b4","#a6cee3","#984ea3"))
legend("topleft",
	legend=c("PiLoci","Anti-PiLoci","Others"),
	col=c("#1f78b4","#a6cee3","#984ea3"),
	pch=15,cex=0.8)

a=map[19,3]
b=map[25,3]
c=map[3,3]
l=c(paste(round(a/c,4)*100,"%",sep=""),
	paste(round(b/c,4)*100,"%",sep=""),
	paste(round((c-a-b)/c,4)*100,"%",sep=""))

pie(c(a,b,c-a-b),
	main="piechar of piRNA Loci Reads Rate of SenseSINEUniqueMappedReads",
	label=l,col=c("#1f78b4","#a6cee3","#984ea3"))
legend("topleft",
	legend=c("PiLoci","Anti-PiLoci","Others"),
	col=c("#1f78b4","#a6cee3","#984ea3"),
	pch=15,cex=0.8)

a=map[19,4]
b=map[25,4]
c=map[3,4]
l=c(paste(round(a/c,4)*100,"%",sep=""),
	paste(round(b/c,4)*100,"%",sep=""),
	paste(round((c-a-b)/c,4)*100,"%",sep=""))

pie(c(a,b,c-a-b),
	main="piechar of piRNA Loci Reads Rate of SenseSINEMultipleMappedReads",
	label=l,col=c("#1f78b4","#a6cee3","#984ea3"))
legend("topleft",
	legend=c("PiLoci","Anti-PiLoci","Others"),
	col=c("#1f78b4","#a6cee3","#984ea3"),
	pch=15,cex=0.8)

#sense LTR
a=map[20,2]
b=map[26,2]
c=map[4,2]
l=c(paste(round(a/c,4)*100,"%",sep=""),
	paste(round(b/c,4)*100,"%",sep=""),
	paste(round((c-a-b)/c,4)*100,"%",sep=""))

pie(c(a,b,c-a-b),
	main="piechar of piRNA Loci Reads Rate of SenseLTRAllMappedReads",
	label=l,col=c("#1f78b4","#a6cee3","#984ea3"))
legend("topleft",
	legend=c("PiLoci","Anti-PiLoci","Others"),
	col=c("#1f78b4","#a6cee3","#984ea3"),
	pch=15,cex=0.8)

a=map[20,3]
b=map[26,3]
c=map[4,3]
l=c(paste(round(a/c,4)*100,"%",sep=""),
	paste(round(b/c,4)*100,"%",sep=""),
	paste(round((c-a-b)/c,4)*100,"%",sep=""))

pie(c(a,b,c-a-b),
	main="piechar of piRNA Loci Reads Rate of SenseLTRUniqueMappedReads",
	label=l,col=c("#1f78b4","#a6cee3","#984ea3"))
legend("topleft",
	legend=c("PiLoci","Anti-PiLoci","Others"),
	col=c("#1f78b4","#a6cee3","#984ea3"),
	pch=15,cex=0.8)

a=map[20,4]
b=map[26,4]
c=map[4,4]
l=c(paste(round(a/c,4)*100,"%",sep=""),
	paste(round(b/c,4)*100,"%",sep=""),
	paste(round((c-a-b)/c,4)*100,"%",sep=""))

pie(c(a,b,c-a-b),
	main="piechar of piRNA Loci Reads Rate of SenseLTRMultipleMappedReads",
	label=l,col=c("#1f78b4","#a6cee3","#984ea3"))
legend("topleft",
	legend=c("PiLoci","Anti-PiLoci","Others"),
	col=c("#1f78b4","#a6cee3","#984ea3"),
	pch=15,cex=0.8)


#length distribution
plot(24:32,as.vector(all_len[,1]/sum(all_len[,1])),type="l",xlab="length (bp)",ylab="percentage (%)",
	main="length distribution of each type of reads",col="#1f78b4",lwd=2,ylim=0.35)
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


