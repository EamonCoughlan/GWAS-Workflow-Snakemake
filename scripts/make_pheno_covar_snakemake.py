import pandas as pd
import ast

#read in data file and config variables, prepare empty variables
data = pd.read_table(snakemake.input[0], header=0, index_col=0, sep=',', na_values=' ')
groupslist = ast.literal_eval(snakemake.config['grouping_values']) #config input for group variables; eval as list instead of string
covars = ['FID'] + ast.literal_eval(snakemake.config['covar_columns']) #second list = config input for covars
covars.remove(snakemake.config['case_control']) #remove pheno column from the covars list. It's requested in config but not wanted here.
phenocodes = str.split(snakemake.config['pheno_code'], sep=' ') #string from config phenotype codes
NaNvalIDs = pd.DataFrame() #empty dataframe to store IDs containing NaN values for mean imputation
if snakemake.config['mi_columns'] != "NULL":
    mi_columns = ast.literal_eval(snakemake.config['mi_columns']) #eval as list
else:
    mi_columns = snakemake.config['mi_columns']

#cleaning data
#replace NaN values with mean where mean-impute specified
if mi_columns != "NULL":
    for i in mi_columns: #config input for mean impute columns
        data[i] = pd.to_numeric(data[i], errors='coerce')
    #some way to list number of NA values for Markdown?
        NaNvals = data[data[i].isnull()]
        NaNvalIDs = pd.concat([NaNvalIDs, NaNvals[snakemake.config['iid_column']]])
        name = 'data/' + snakemake.wildcards.sample + '/' + i + '.NaNvals'
        NaNvals.to_csv(name, sep='\t', index=False)
        data[i].fillna(data[i].mean(), inplace=True)
        #data[i] = data[i].round(0) #if need/want to round values to integer?
        #print(data[i])
    name = 'data/' + snakemake.wildcards.sample + '/NaNvalsReplaced.IDs'
    NaNvalIDs.to_csv(name, index=False, header=False)

#remove missing values where no mean-impute desired
#not actually necessary since means are imputed first
#no_mi = data.drop(['AGE', 'EXTRA'], axis=1, inplace=True) #config variable for MI
#print(no_mi)

#collect stats on how many individuals removed
NaNremoved = data[data.isnull().any(axis=1)]
RemovedIDs = NaNremoved[snakemake.config['iid_column']]
name = 'data/' + snakemake.wildcards.sample + '/NaNvalsRemoved.IDs'
RemovedIDs.to_csv(name, index=False, header=False)
stats = {'n mean Imputed': [len(NaNvalIDs)], 'n Removed': [len(RemovedIDs)]}
ungrouped = len(data[snakemake.config['iid_column']]) #to subtract groups from for final count

#recode phenotypes to 0 1
data[snakemake.config['case_control']] = data[snakemake.config['case_control']].astype(str) #encodes phenotypes as strings for matching with config
data[snakemake.config['case_control']] = data[snakemake.config['case_control']].map({phenocodes[0]: 0, phenocodes[1]: 1}) #replaces phenotypes with 0 1 (numeric)
print(data)

def WriteFiles(group_df, group):
    phenofile = group_df[['FID',snakemake.config['iid_column'],snakemake.config['case_control']]]
    phenofile.rename(columns={snakemake.config['iid_column']: 'IID'}, inplace=True) #sets column name back to IID for regenie
    name = 'data/' + snakemake.wildcards.sample + '.' + group + '.pheno'
    phenofile.to_csv(name, sep='\t', index=False)
    idfile = phenofile[['FID', 'IID']]
    name = 'data/' + snakemake.wildcards.sample + '.' + group + '.idlist'
    idfile.to_csv(name, sep='\t', index=False)
    covarfile = group_df[covars]
    covarfile.rename(columns={snakemake.config['iid_column']: 'IID'}, inplace=True) #as above
    name = 'data/' + snakemake.wildcards.sample + '.' + group + '.covar'
    covarfile.to_csv(name, sep='\t', index=False)
    stats[group] = len(idfile['IID']) #record how many in this group for Markdown
    global ungrouped
    ungrouped -= len(idfile['IID']) #remove this number from total count
#subsetting by group
for group in groupslist:
    subset = data[data[snakemake.config['grouping_column']]==group]
    WriteFiles(subset, group)

stats['Other_or_Unassigned'] = ungrouped
stats = pd.DataFrame(stats)
stats.to_csv(snakemake.output[3], sep=',', index=False) #name via snakemake output
