import pandas as pd

#read in covars and subset it according to config
covars = pd.read_table(snakemake.input[0],header=0,index_col=None)

#read in pca values and generate column names
pca = pd.read_table(snakemake.input[1],header=None,index_col=None, sep='\s+')
list = []
for i in pca.columns:
    list.append(i-1) #reindex columns to count 1 to n with space for FID and IID columns (replacing -1 and 0 index)
    list[i] = 'PC' + str(list[i]) #name with 'PC' prefix for clarity
list[0] = 'FID'
list[1] = 'IID'
pca.columns = list

#merge covar and pca values
addedpca = pd.merge(covars, pca, how='inner', on=['IID', 'FID'])

addedpca.to_csv(snakemake.output[0], sep='\t', index=False)