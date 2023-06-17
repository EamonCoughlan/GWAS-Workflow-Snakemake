#import the bim file as a dataframe
import pandas as pd

file = pd.read_table(snakemake.input[0], sep='\t', header = None)

df = pd.DataFrame(file)

#iterate through rows and make a list of SNP IDs (column 1) where SNP strand is ambiguous (complementary bases)
ambiguous = []
for index, row in df.iterrows():
    if (row[4] == 'A' and row[5] == 'T')\
            or (row[4] == 'T' and row[5] == 'A')\
            or (row[4] == 'C' and row[5] == 'G')\
            or (row[4] == 'G' and row[5] == 'C'):
                ambiguous.append(row[1])

#write the list of ambiguous SNP IDs to a text document to exclude using PLINK:
with open(snakemake.output[0], 'w') as outfile:
    for SNP in ambiguous:
        outfile.write(SNP + '\n')