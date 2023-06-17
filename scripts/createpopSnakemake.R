#!/usr/bin/env Rscript

fam<-read.table(snakemake@input[[1]], header=F, stringsAsFactors=F) #fam file of merged covid 1000GP
kgp<-read.table(snakemake@input[[2]], header=T, stringsAsFactors=F) #kgp.superpop file
 pop<-cbind(fam[,1:2], "-")
pop<-as.matrix(pop)
 rownames(pop)<-pop[,2]
 pop[kgp[,2],3]<-kgp[,"POP"]
 pop<-pop[fam[,2],]

write.table(pop[,3], col.names=F, row.names=F, quote=F, file=snakemake@output[[1]]) #"merged_covid19_kgp_pruned.pop"