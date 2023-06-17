import pandas as pd
import os

#import chr{i}.bim into pandas df
bim = pd.read_table(snakemake.input[0], header=None, index_col=None, sep='\s+')

#get genpos (col 4) at every 500000 snps
#make each chr+genpos into a range (/batch) for regenie
batches = []
i = 0
while i <=len(bim):
    start = bim[3][i]
    i += 500000
    if i <= len(bim):
        end = bim[3][i]
        range = str(bim[0][0]) + ':' + str(start) + '-' + str(end)
        batches.append(range)
    else:
        end = bim[3][len(bim)-1]
        range = str(bim[0][0]) + ':' + str(start) + '-' + str(end)
        batches.append(range)
#print(batches) - check output is accurate


#Write each 'batch' to a file stored under the appropriate chromosome and named with the batch region
outdir = 'data/' + str(snakemake.wildcards['sample']) + '/chr' + str(snakemake.wildcards['i']) + '_batches/'
makecommand = 'mkdir ' + outdir
os.system(makecommand)
for batch in batches:
    outfile = outdir + batch + '.batch'
    with open(outfile, 'a') as file:
        file.write(batch+'\n')