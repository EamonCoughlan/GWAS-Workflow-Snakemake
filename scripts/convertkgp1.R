#!/usr/bin/env Rscript
bim<-read.table(snakemake@input[[1]], stringsAsFactors=F)
bim<-cbind(bim, paste(bim[,1], bim[,4], bim[,6], bim[,5], sep=":"))
kgp<-read.table(snakemake@input[[2]], stringsAsFactors=F)
write.table(bim[!bim[,7]%in%kgp[,2],2], col.names=F, row.names=F, quote=F, file=snakemake@output[[1]])