import pandas as pd
import ast

covar_columns = ast.literal_eval(snakemake.config['covar_columns']) #read columns param as list

#read in ancestry file and subset the IDs and ancestry assignment, assign config ID
ancestry = pd.read_table(snakemake.input[0],header=0,index_col=None, sep='\s+')
ancestry = ancestry[['FID','IID','POP']]
ancestry.rename(columns = {'IID':snakemake.config['iid_column']}, inplace = True) #iid from config

#read in covars and subset it according to config
covars = pd.read_csv(snakemake.input[1],header=0,index_col=None)
covars = covars[covar_columns]
covars[snakemake.config['iid_column']] = covars[snakemake.config['iid_column']].str.lower() #force lowercase all IDs to ensure match with other files

#merge with other data
allmerge = pd.merge(ancestry, covars, how='inner', on=[snakemake.config['iid_column']])

#write merge to new file
allmerge.to_csv(snakemake.output[0], sep=',')