#!/usr/bin/R

args <- commandArgs()
inputf <- args[6]
outputf <- args[7]

dat <- read.table(inputf,header=F,stringsAsFactors=F)
dat <- dat[dat[,2]==1,]
re <- dat[,-c(1:8)]
write.table(re,file=outputf,col.names=F,row.names=F,sep="\t",quote=F)
