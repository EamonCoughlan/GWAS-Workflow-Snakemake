#Make pandas dataframe from king '.kin0' file and extract ID columns:
import pandas as pd

file = pd.read_table(snakemake.input[0], sep='\t', header = 0)

df = pd.DataFrame(file)
ids = df[['ID1', 'ID2']]

#Konrad's code
# 'You need a pandas dataframe with each row a pair of related individuals encoded by two columns (ID1, ID2)
# which are the ids of individuals in the pair. You call related_pivots with that dataframe as only argument
# and it returns a list of individual IDs such that if these IDs are excluded then there are no pairs in the data.'
def reduce(rel):
    s = pd.concat([rel.groupby('ID1').size(), rel.groupby('ID2').size()],axis=1,keys=[0,1]).fillna(0.0).sum(axis=1)
    who = s.idxmax()
    return who,s[who],rel[~(rel[['ID1','ID2']] == who).any(axis=1)].copy()

def related_pivots(rel):
    r = rel.copy()
    who = []
    while r.shape[0] > 0:
        w,c,r = reduce(r)
        who += [w]
        if c == 1.0:
            who += list(r.ID1)
            break
    return who

#Apply functions to the IDs and writes with corresponding FID1 to an output file
relatedids = related_pivots(ids) #finds related individuals to be removed using function
relateddf = df.query('ID1 in @relatedids') #extracts only these individuals from original df
output = relateddf[['FID1', 'ID1']].drop_duplicates(subset=['ID1']) #takes FID1 and ID1 columns, removing duplicates
output.to_csv(snakemake.output[0], sep='\t', index = False, header = False) #writes to PLINK-readable file