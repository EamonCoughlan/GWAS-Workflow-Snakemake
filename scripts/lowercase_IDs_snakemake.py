import pandas as pd

famfile = pd.read_table(snakemake.input[0],header=None,index_col=None, sep='\s+') #input = fam file
ids = famfile.iloc[:, [0, 1]]
ids[2]=ids[0]
ids[3]=ids[1].str.lower()
ids.to_csv(snakemake.output[0], sep='\t', header=False, index=False) #output = remap file for plink