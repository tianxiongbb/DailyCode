#!/usr/bin/env Rscript

###library packages
library("Sushi")

###arguments
Args=commandArgs()
if(length(Args)<6){
  stop(paste(substr(Args[4],8,100),"rep.size outFile bedGraph.watson [bedGraph.crick]",sep=" "))
}
rep_size=read.table(Args[6],header=F,row.names=NULL,comment.char="&")
outFile=Args[7]
bdgFileWatson=read.table(Args[8],header=F,row.names=NULL,comment.char="&")

fun_plot1=function(x){
  p <<- p+1
  if(p%%100==0){
    print(paste(floor(p)," reads processed......",sep=" "))
  }
  tran=as.character(x[1])
  plotBedgraph(bdgFileWatson,tran,1,as.numeric(x[2]),color="#d53e4f")
  axis(2)
  mtext("Read Depth",2,line=2,cex=1,font=2)
  mtext(tran,3,line=1,cex=1.25,font=2)
  labelgenome("",1,as.numeric(x[2]),n=5,scale="Kb")
}
fun_plot2=function(x){
  p <<- p+1
  if(p%%100==0){
    print(paste(floor(p)," reads processed......",sep=" "))
  }
  par(mar=c(0.3,4,3,2))
  tran=as.character(x[1])
  plotBedgraph(bdgFileWatson,tran,1,as.numeric(x[2]),color="#d53e4f")
  axis(2)
  mtext(tran,3,line=1,cex=1.25,font=2)
  mtext("Read Depth",2,line=2,cex=1,font=2)
  par(mar=c(3,4,0.3,2))
  plotBedgraph(bdgFileCrick,tran,1,as.numeric(x[2]),flip=T,color="#3288bd")
  axis(2)
  mtext("Read Depth",2,line=2,cex=1,font=2)
  labelgenome("",1,as.numeric(x[2]),n=5,scale="Kb")
}

if(length(Args)<=8){
  pdf(outFile,width=6,height=2.5)
  par(mar=c(3,4,3,2))
  p=0
  apply(rep_size,1,fun_plot1)
}else{
  bdgFileCrick=read.table(Args[9],header=F,row.names=NULL)
  bdgFileCrick[,4]=abs(bdgFileCrick[,4])
  pdf(outFile,width=6,height=4.2)
  par(mfrow=c(2,1))
  p=0
  apply(rep_size,1,fun_plot2)
}
dev.off()

