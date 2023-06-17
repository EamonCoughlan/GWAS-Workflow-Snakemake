#!/usr/bin/env Rscript

ids<-read.table(snakemake@input[[1]], header=T, stringsAsFactors=FALSE)

id2 <- ids[c("FID2", "ID2")] #id2 <- ids[c("FID2", "ID2")] for 2 column version; id2 <- ids["FID2"] for 1 column version

write.table(id2, file = snakemake@output[[1]], row.names = FALSE, col.names = FALSE, quote = FALSE)