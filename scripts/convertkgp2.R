#!/usr/bin/env Rscript
bim<-read.table(snakemake@input[[1]], stringsAsFactors=F)
bim<-cbind(bim, paste(bim[,1], bim[,4], bim[,6], bim[,5], sep=":"))
kgp<-read.table(snakemake@input[[2]], stringsAsFactors=F)

write.table(bim[bim[,7]%in%kgp[,2],c(2,7)], col.names=F, row.names=F, quote=F, file=snakemake@output[[1]])
write.table(bim[bim[,7]%in%kgp[,2],7], col.names=F, row.names=F, quote=F, file=snakemake@output[[2]])
bim[,2]<-bim[,7]
write.table(bim[duplicated(bim[,2]),2], col.names=F, row.names=F, quote=F, file=snakemake@output[[3]])
write.table(bim[,1:6], col.names=F, row.names=F, quote=F, file=snakemake@output[[4]])