setwd("D:/research/gaw19/G19")

source("hwu.core.R")

Gexpr<-read.table("G.tab",header=T,row.names=1)
Snp<-read.table("S.tab",header=T,row.names=1)
Xcov<-read.table("X.tab",header=T,row.names=1)
Y<-read.table("Y.tab",header=T,row.names=1)

geno<-collapse.burden(Snp)
K<-weight.gaussian(as.matrix(Gexpr))

rst<-HWU(Y[,1],K,geno,as.matrix(Xcov))




