#Import pcl file
pcloadings <- read.table(snakemake@input[[1]], header = TRUE, stringsAsFactors = FALSE)

#Merge columns Pc-1 -2 -3 into one column and sort low-high
newcolumn <- pcloadings$pc1_loading + pcloadings$pc2_loading + pcloadings$pc3_loading
sorted <- sort(newcolumn)
#Find value x at Nth percentile (eg. index length*0.67) N = threshold value from config
index <- (length(sorted) * as.numeric(snakemake@params[['PCLthreshold']]))
value <- sorted[index]

#From original table, filter out any SNP with value > x in any of Pc-1 -2 -3

nopc1 <- pcloadings[pcloadings$pc1_loading < value,]
nopc2 <- nopc1[nopc1$pc2_loading < value,]
nopc3 <- nopc2[nopc2$pc3_loading < value,]

lowpcSNPs <- nopc3$SNP

#Write SNPs to a new file to be extracted via plink
write.table(lowpcSNPs, file = snakemake@output[[1]], sep = '\n', quote = FALSE, row.names = FALSE, col.names = FALSE)

