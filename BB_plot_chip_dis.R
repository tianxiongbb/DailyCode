#!/usr/bin/env Rscript

###library packages
library("Sushi")

###arguments
Args=commandArgs()
if(length(Args)<6){
  stop(paste(substr(Args[4],8,100),"gene.bed transcipt.bed7 outPrefix bedGraph.watson [bedGraph.crick]",sep=" "))
}
gene_cor=read.table(Args[6],header=F,row.names=NULL)
tran_cor=read.table(Args[7],header=F,row.names=NULL)
outPrefix=Args[8]
bdgFile=read.table(Args[9],header=F,row.names=NULL)

fun_plot=function(x){
  chrom=as.character(x[1])
  chromstart=as.numeric(x[2])
  chromend=as.numeric(x[3])
  leng=chromend-chromstart
  chromstart=max(1,min(chromstart-leng*3,chromstart-30000))
  chromend=max(chromend+leng*3,chromend+30000)
  gn=as.character(x[5])
  pdf(paste(outPrefix,gn,"pdf",sep="."),width=6,height=6)
  layout(matrix(c(1,2,2,3,3),5,1))
  par(mar=c(0.3,4,3,2))
  ti=tran_cor[which(tran_cor[,1]==chrom & tran_cor[,2]<chromend & tran_cor[,3]>chromstart),]
  ti[,1]=as.character(ti[,1])
  ti[,4]=as.character(ti[,4])
  ti[,7]=as.character(ti[,7])
  plotGenes(ti,chrom,chromstart,chromend,labeltext=F,
            types=ti[,7],colorby=ti[,6],colorbyrange=c(-1,1),
            colorbycol=colorRampPalette(c("blue","purple","black","orange","red")))
  mtext("Transcripts",2,line=2,cex=1,font=2)
  mtext(gn,3,line=1,cex=1.25,font=2)
  par(mar=c(3,4,0.3,2))
  plotBedgraph(bdgFile,chrom,chromstart,chromend,colorbycol=SushiColors(5))
  labelgenome(chrom,chromstart,chromend,n=5,scale="Kb")
  axis(2,tcl=0.2)
  mtext("Read Depth",2,line=2,cex=1,font=2)
  chromstart=max(1,as.numeric(x[2])-2000)
  chromend=as.numeric(x[3])+2000
  zoomregion=c(chromstart,chromend)
  zoomsregion(zoomregion,extend=c(0.01,0.13))
  par(mar=c(3,4,0.3,2))
  plotBedgraph(bdgFile,chrom,chromstart,chromend,colorbycol=SushiColors(5))
  labelgenome(chrom,chromstart,chromend,n=5,scale="Kb")
  zoombox()
  mtext("Read Depth",2,line=2,cex=1,font=2)
  axis(2,tcl=0.2)
  dev.off()
}
apply(gene_cor,1,fun_plot)
