#!/usr/bin/env Rscript
kgp<-read.table(snakemake@input[[1]], header=T, stringsAsFactors=F)
 admx<-read.table(snakemake@input[[2]], header=F, stringsAsFactors=F)
 fam<-read.table(snakemake@input[[3]], header=F, stringsAsFactors=F)
 pop<-read.table(snakemake@input[[4]], header=F, stringsAsFactors=F)
 admx<-cbind(fam, pop, admx)
 admx<-admx[,-(3:6)]
colnames(admx)<-c("FID", "IID", "POP", "EUR", "EAS", "AMR", "SAS", "AFR")
covid19<-admx[admx[,"FID"]!=0,]
covid19<-as.matrix(covid19)
 covid19[covid19[,"EUR"]>=0.8,"POP"]<-"EUR"
 covid19[covid19[,"EAS"]>=0.8,"POP"]<-"EAS"
 covid19[covid19[,"SAS"]>=0.8,"POP"]<-"SAS"
 covid19[covid19[,"AFR"]>=0.8,"POP"]<-"AFR"
 covid19[covid19[,"AMR"]>=0.8,"POP"]<-"AMR"
write.table(covid19, col.names=T, row.names=F, quote=F, file=snakemake@output[[1]])