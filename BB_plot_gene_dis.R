#!/usr/bin/env Rscript

###library packages
library("Sushi")

###arguments
Args=commandArgs()
if(length(Args)<6){
  stop(paste(substr(Args[4],8,100),"gene.bed transcipt.bed7 outFile bedGraph.watson [bedGraph.crick]",sep=" "))
}
gene_cor=read.table(Args[6],header=F,row.names=NULL)
tran_cor=read.table(Args[7],header=F,row.names=NULL)
outFile=Args[8]
bdgFileWatson=read.table(Args[9],header=F,row.names=NULL)

fun_plot1=function(x){
  p <<- p+1
  if(p%%1000==0){
    print(paste(floor(p)," reads processed......",sep=" "))
  }
  chrom=as.character(x[1])
  chromstart=max(1,as.numeric(x[2])-2000)
  chromend=as.numeric(x[3])+2000
  gn=as.character(x[5])
  plotBedgraph(bdgFileWatson,chrom,chromstart,chromend,color="#d53e4f")
  axis(2)
  mtext("Read Depth",2,line=2,cex=1,font=2)
  mtext(gn,3,line=1,cex=1.25,font=2)
  par(mar=c(3,4,0.3,2))
  ti=tran_cor[which(tran_cor[,1]==chrom & tran_cor[,2]<chromend & tran_cor[,3]>chromstart),]
  ti[,1]=as.character(ti[,1])
  ti[,4]=as.character(ti[,4])
  ti[,7]=as.character(ti[,7])
  plotGenes(ti,chrom,chromstart,chromend,labeltext=F,
            types=ti[,7],colorby=ti[,6],colorbyrange=c(-1,1),
            colorbycol=colorRampPalette(c("blue","purple","black","orange","red")))
  labelgenome(chrom,chromstart,chromend,n=5,scale="Kb")
}
fun_plot2=function(x){
  p <<- p+1
  if(p%%1000==0){
    print(paste(floor(p)," reads processed......",sep=" "))
  }
  par(mar=c(0.3,4,3,2))
  chrom=as.character(x[1])
  chromstart=max(1,as.numeric(x[2])-2000)
  chromend=as.numeric(x[3])+2000
  gn=as.character(x[5])
  plotBedgraph(bdgFileWatson,chrom,chromstart,chromend,color="#d53e4f")
  axis(2)
  mtext(gn,3,line=1,cex=1.25,font=2)
  mtext("Read Depth",2,line=2,cex=1,font=2)
  par(mar=c(3,4,0.3,2))
  plotBedgraph(bdgFileCrick,chrom,chromstart,flip=T,chromend,color="#3288bd")
  axis(2)
  mtext("Read Depth",2,line=2,cex=1,font=2)
  par(mar=c(3,4,0.3,2))
  ti=tran_cor[which(tran_cor[,1]==chrom & tran_cor[,2]<chromend & tran_cor[,3]>chromstart),]
  ti[,1]=as.character(ti[,1])
  ti[,4]=as.character(ti[,4])
  ti[,7]=as.character(ti[,7])
  plotGenes(ti,chrom,chromstart,chromend,labeltext=F,
            types=ti[,7],colorby=ti[,6],colorbyrange=c(-1,1),
            colorbycol=colorRampPalette(c("blue","purple","black","orange","red")))
  labelgenome(chrom,chromstart,chromend,n=5,scale="Kb")
}

if(length(Args)<=9){
  pdf(outFile,width=6,height=4.2)
  par(mfrow=c(2,1),mar=c(0.3,4,3,2))
  p=0
  apply(gene_cor,1,fun_plot1)
}else{
  bdgFileCrick=read.table(Args[10],header=F,row.names=NULL)
  bdgFileCrick[,4]=abs(bdgFileCrick[,4])
  pdf(outFile,width=6,height=6)
  par(mfrow=c(3,1))
  p=0
  apply(gene_cor,1,fun_plot2)
}
dev.off()

