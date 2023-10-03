SAMPLES = glob_wildcards('data/{condition}.ped').condition
configfile: 'files/gwas_config.yaml'
#comment out this print statement if you are creating an image of the dag or it won't work
#print('Config is ', config)

#TARGET RULES FOR DIVISION OF PIPELINE
#This runs when you run 'snakemake' with no arguments - as long as it is kept as the first rule in the file.
rule all:
    input:
        expand('results/{sample}-results_Markdown.html', sample=SAMPLES),
        expand('results/{sample}.BasicQCSummary.html', sample=SAMPLES),
        expand('results/{sample}.PCA-QCSummary.html', sample=SAMPLES)

#Run 'snakemake basicqc' to perform only basicQC (up to output file PostBasicQC.bed / .bim / .fam)
rule basicqc:
    input:
        expand('data/{sample}.PostBasicQC.bed', sample=SAMPLES),
        expand('results/{sample}.BasicQCSummary.html', sample=SAMPLES)

#Run 'snakemake pca_check' to perform everything up to PCA and plotting
rule pca_check:
    input:
        expand('results/{sample}.PCA-QCSummary.html', sample=SAMPLES)

#Run 'snakemake ancestry_check' to conduct ancestry analysis with Admixture software
rule ancestry_check:
    input:
        expand('data/{sample}.estimatedancestry.txt', sample=SAMPLES)

#Run 'snakemake topmed_imputation' to submit data to Topmed for imputation and automatically download (warning; large files)
rule topmed_imputation:
    input:
        expand('data/{sample}_imputation', sample = SAMPLES)

#R Markdown rules (For visual summaries of QC)
rule RMarkdown_BasicQC:
    input:
        'data/{sample}.frq',
        'data/{sample}.samplesout.log',
        'data/{sample}.SNPsout.log',
        'data/{sample}.samplesout2.frq',
        'data/{sample}.samplesout2.log',
        'data/{sample}.hwefiltered.log',
        'data/{sample}.duplicate_samples.txt'
    output:
        'results/{sample}.BasicQCSummary.html'
    params:
        hwevalue = config['hwe'],
        mind1value = config['mind1'],
        genovalue = config['geno1'],
        mafvalue = config['maf1'],
        mind2value = config['mind2']
    script:
        'scripts/BasicQCMarkdownSnakemake.Rmd'

rule RMarkdown_PCA:
    input:
        'data/{sample}.PostBasicQC.bim',
        'data/{sample}.highLDexcluded.bim',
        'data/{sample}.ambiguousSNPs.txt',
        'data/{sample}.prekinship.log',
        'data/{sample}.prekinshipmaf.frq',
        'data/{sample}.related_d3.txt',
        'data/{sample}.unrelated.fam',
        'data/{sample}.prune.in',
        'data/{sample}.prune.out',
        'data/{sample}.lowPCSNPs',
        'data/{sample}.related_d3_round2.txt',
        'data/{sample}.highPCexcluded.fam',
        'data/{sample}_snps_toflip.txt',
        'data/{sample}.SNPsCommonKGP',
        'data/{sample}.SNPsdup',
        'data/{sample}-1000GP.prune.in',
        'data/{sample}-1000GP.prune.out',
        'data/{sample}.unrelated_merge.eigenvec',
        'data/{sample}.mergedpruned.pop'
    output:
        'results/{sample}.PCA-QCSummary.html'
    params:
        maf2value = config['maf2'],
        geno2value = config['geno2'],
        ipvalue = config['indep_pairwise'],
        PCLthreshold = config['pcl_threshold'],
        ipvalue2 = config['indep_pairwise2'],
        PCnumber = config['pc_number']
    script:
        'scripts/PCAMarkdownSnakemake.Rmd'

#BASIC-QC

#generate binary filetypes from ped / map & generate pre-QC frequency table of minor alleles
#also putting into allele ACGT format at this point, if it isn't already; this is needed for later ambiguous SNP detection
rule make_binary_file:
    input:
        'data/{sample}.ped',
        'data/{sample}.map'
    output:
        'data/{sample}.bed',
        'data/{sample}.bim',
        'data/{sample}.fam',
        temp('data/{sample}.frq')
    shell:
        'plink --noweb --make-bed --file data/{wildcards.sample} --alleleACGT --freq --out data/{wildcards.sample}'

#generate ID-remapping file for normalising all IDs to lowercase:
rule make_id_remap:
    input:
        'data/{sample}.fam'
    output:
        temp('data/{sample}.id-remap')
    script:
        'scripts/lowercase_IDs_snakemake.py'

rule lowercase_IDs:
    input:
        'data/{sample}.bed',
        'data/{sample}.bim',
        'data/{sample}.fam',
        'data/{sample}.id-remap'
    output:
        temp('data/{sample}.idnorm.bed'),
        temp('data/{sample}.idnorm.bim'),
        temp('data/{sample}.idnorm.fam')
    shell:
        'plink --bfile data/{wildcards.sample} --update-ids data/{wildcards.sample}.id-remap --make-bed --out data/{wildcards.sample}.idnorm'

#remove samples with high number of missing SNPs
rule remove_highmissing_samples:
    input:
        'data/{sample}.idnorm.bed',
        'data/{sample}.idnorm.bim',
        'data/{sample}.idnorm.fam'
    output:
        temp('data/{sample}.samplesout.bed'),
        temp('data/{sample}.samplesout.bim'),
        temp('data/{sample}.samplesout.fam'),
        temp('data/{sample}.samplesout.log'), #stats for markdown
    params:
        mindvalue = config['mind1']
    shell:
        'plink --bfile data/{wildcards.sample}.idnorm --mind {params.mindvalue} --make-bed --out data/{wildcards.sample}.samplesout'

#remove SNPs with high rate of missingness across samples
rule remove_highmissing_SNPs:
    input:
        'data/{sample}.samplesout.bed',
        'data/{sample}.samplesout.bim',
        'data/{sample}.samplesout.fam'
    output:
        temp('data/{sample}.SNPsout.bed'),
        temp('data/{sample}.SNPsout.bim'),
        temp('data/{sample}.SNPsout.fam'),
        temp('data/{sample}.SNPsout.log'),
    params:
        genovalue = config['geno1']
    shell:
        'plink --bfile data/{wildcards.sample}.samplesout --geno {params.genovalue} --make-bed --out data/{wildcards.sample}.SNPsout'

#remove SNPs with low minor allele frequency (only feed into to PCA not GWAS?)
rule minor_allele_frequency_cutoff:
    input:
        'data/{sample}.SNPsout.bed',
        'data/{sample}.SNPsout.bim',
        'data/{sample}.SNPsout.fam'
    output:
        temp('data/{sample}.MAFcutoff.bed'),
        temp('data/{sample}.MAFcutoff.bim'),
        temp('data/{sample}.MAFcutoff.fam'),
    params:
        mafvalue = config['maf1']        
    shell:
        'plink --bfile data/{wildcards.sample}.SNPsout --maf {params.mafvalue} --make-bed --out data/{wildcards.sample}.MAFcutoff'

#remove samples with high number of missing SNPs (mind 0.03) - and take minor allele frequency (freq) for R Markdown to measure effect of previous step.
rule remove_highmissing_samples_2:
    input:
        'data/{sample}.MAFcutoff.bed',
        'data/{sample}.MAFcutoff.bim',
        'data/{sample}.MAFcutoff.fam'
    output:
        temp('data/{sample}.samplesout2.bed'),
        temp('data/{sample}.samplesout2.bim'),
        temp('data/{sample}.samplesout2.fam'),
        temp('data/{sample}.samplesout2.log'),
        temp('data/{sample}.samplesout2.frq')
    params:
        mind2value = config['mind2']
    shell:
        r'''
        plink --bfile data/{wildcards.sample}.MAFcutoff --freq --out data/{wildcards.sample}.samplesout2
        plink --bfile data/{wildcards.sample}.MAFcutoff --mind {params.mind2value} --make-bed --out data/{wildcards.sample}.samplesout2
        '''

#filter for HWE
rule filter_hwe:
    input:
        'data/{sample}.samplesout2.bed',
        'data/{sample}.samplesout2.bim',
        'data/{sample}.samplesout2.fam'
    output:
        temp('data/{sample}.hwefiltered.bed'),
        temp('data/{sample}.hwefiltered.bim'),
        temp('data/{sample}.hwefiltered.fam'),
        temp('data/{sample}.hwefiltered.log')
    params:
        hwevalue = config['hwe']
    shell:
        'plink --bfile data/{wildcards.sample}.samplesout2 --hwe {params.hwevalue} --make-bed --out data/{wildcards.sample}.hwefiltered'

#remove duplicates (king-dependent - possibly not required depending on data source)
rule find_duplicates:
    input:
        'data/{sample}.hwefiltered.bed',
        'data/{sample}.hwefiltered.bim',
        'data/{sample}.hwefiltered.fam'
    output:
        'data/{sample}.con'
    shell:
        'king -b data/{wildcards.sample}.hwefiltered.bed --duplicate --prefix data/{wildcards.sample}'

rule isolate_duplicates:
    input:
        'data/{sample}.con'
    output:
        'data/{sample}.duplicate_samples.txt'
    script:
        'scripts/remove_duplicates.R'

#Remove duplicates from the text file
#Was using FID and --remove-fam; adapted to 2 column.
#NB --remove will require a 2 column file with FID + ID, --remove fam only needs 1 column.
rule remove_duplicates:
    input:
        'data/{sample}.hwefiltered.bed',
        'data/{sample}.hwefiltered.bim',
        'data/{sample}.hwefiltered.fam',
        'data/{sample}.duplicate_samples.txt'
    output:
        'data/{sample}.PostBasicQC.bed',
        'data/{sample}.PostBasicQC.bim',
        'data/{sample}.PostBasicQC.fam',
    shell:
        'plink --bfile data/{wildcards.sample}.hwefiltered --remove data/{wildcards.sample}.duplicate_samples.txt --make-bed --out data/{wildcards.sample}.PostBasicQC'

#SAMPLE QC (PCA)

#Exclude SNPs in high LD regions <- note we will come back to file before this later (PostBasicQC) for analysis!
rule exclude_highLD_regions:
    input:
        'data/{sample}.PostBasicQC.bed',
        'data/{sample}.PostBasicQC.bim',
        'data/{sample}.PostBasicQC.fam',
        'files/snps.tofilter'
    output:
        temp('data/{sample}.highLDexcluded.bed'),
        temp('data/{sample}.highLDexcluded.bim'),
        temp('data/{sample}.highLDexcluded.fam'),
    shell:
        'plink --bfile data/{wildcards.sample}.PostBasicQC --exclude range files/snps.tofilter --make-bed --out data/{wildcards.sample}.highLDexcluded'

#find ambiguous SNPs and write to file:
rule find_amb_SNPs:
    input:
        'data/{sample}.highLDexcluded.bim'
    output:
        'data/{sample}.ambiguousSNPs.txt'
    script:
        'scripts/FindAmbiguousSNPs_Snakemake.py'

#remove ambiguous SNPs
rule remove_amb_SNPS:
    input:
        'data/{sample}.highLDexcluded.bed',
        'data/{sample}.highLDexcluded.bim',
        'data/{sample}.highLDexcluded.fam',
        'data/{sample}.ambiguousSNPs.txt'
    output:
        temp('data/{sample}.ambSNPsremoved.bed'),
        temp('data/{sample}.ambSNPsremoved.bim'),
        temp('data/{sample}.ambSNPsremoved.fam')    
    shell:
        'plink --bfile data/{wildcards.sample}.highLDexcluded --exclude data/{wildcards.sample}.ambiguousSNPs.txt --make-bed --out data/{wildcards.sample}.ambSNPsremoved' 

#Filter for Minor Allele Frequency and SNPs with high missingness again
rule pre_kinship_filters:
    input:
        'data/{sample}.ambSNPsremoved.bed',
        'data/{sample}.ambSNPsremoved.bim',
        'data/{sample}.ambSNPsremoved.fam'
    output:
        temp('data/{sample}.prekinship.bed'),
        temp('data/{sample}.prekinship.bim'),
        temp('data/{sample}.prekinship.fam'),
        temp('data/{sample}.prekinship.log'),
    params:
        mafvalue = config['maf2'],
        genovalue = config['geno2']
    shell:
        'plink --bfile data/{wildcards.sample}.ambSNPsremoved --maf {params.mafvalue} --geno {params.genovalue} --make-bed --out data/{wildcards.sample}.prekinship'

#Generate MAF stats for Markdown
rule maf_graph_prekinship:
    input:
        'data/{sample}.prekinship.bed',
        'data/{sample}.prekinship.bim',
        'data/{sample}.prekinship.fam'
    output:
        temp('data/{sample}.prekinshipmaf.frq')
    shell:
        'plink --bfile data/{wildcards.sample}.prekinship --freq --out data/{wildcards.sample}.prekinshipmaf'

#Kinship testing (King)
rule kinship_testing:
    input:
        'data/{sample}.prekinship.bed',
        'data/{sample}.prekinship.bim',
        'data/{sample}.prekinship.fam'
    output:
        temp('data/{sample}.kin0')
    shell:
        'king -b data/{wildcards.sample}.prekinship.bed --kinship --degree 3 --prefix data/{wildcards.sample}'

#calculate largest independent set of unrelated individuals using Konrad's script
rule calculate_largest_independent_unrelated_group:
    input:
        'data/{sample}.kin0',
        'data/{sample}.prekinship.fam'
    output:
        temp('data/{sample}.related_d3.txt')
    script:
        'scripts/FindRelatedIndividualsSnakemake.py'

#remove related individuals 
rule remove_related_individuals:
    input:
        'data/{sample}.prekinship.bed',
        'data/{sample}.prekinship.bim',
        'data/{sample}.prekinship.fam',
        'data/{sample}.related_d3.txt'
    output:
        temp('data/{sample}.unrelated.bed'),
        temp('data/{sample}.unrelated.bim'),
        temp('data/{sample}.unrelated.fam')
    shell:       
        'plink --bfile data/{wildcards.sample}.prekinship --remove data/{wildcards.sample}.related_d3.txt --make-bed --out data/{wildcards.sample}.unrelated'

#pruning for SNPs in LD - 1000 80 0.1 (specs controlled via config)
rule prune_LD_SNPs:
    input:
        'data/{sample}.unrelated.bed',
        'data/{sample}.unrelated.bim',
        'data/{sample}.unrelated.fam',
    output:
        temp('data/{sample}.prune.out'),   
        temp('data/{sample}.prune.in')
    params:
        ipvalue = config['indep_pairwise']
    shell:
        'plink --bfile data/{wildcards.sample}.unrelated --indep-pairwise {params.ipvalue} --out data/{wildcards.sample}'

#remove pruned LD SNPs
rule remove_pruned:
    input:
        'data/{sample}.unrelated.bed',
        'data/{sample}.unrelated.bim',
        'data/{sample}.unrelated.fam',
        'data/{sample}.prune.out',
        'data/{sample}.prune.in'
    output:
        temp('data/{sample}.unrelatedpruned.bed'),
        temp('data/{sample}.unrelatedpruned.bim'),
        temp('data/{sample}.unrelatedpruned.fam')
    shell:
        'plink --bfile data/{wildcards.sample}.unrelated --extract data/{wildcards.sample}.prune.in --make-bed --out data/{wildcards.sample}.unrelatedpruned' 
        
#calculate pca 
rule pca:
    input:
        'data/{sample}.unrelatedpruned.bed',
        'data/{sample}.unrelatedpruned.bim',
        'data/{sample}.unrelatedpruned.fam'
    output:
        'data/{sample}.unrelatedpruned.eigenvec',
        'data/{sample}.unrelatedpruned.eigenval'
    shell:
        'plink --bfile data/{wildcards.sample}.unrelatedpruned --pca 10 --out data/{wildcards.sample}.unrelatedpruned'

#calculate pca loadings
rule pca_loadings:
    input:
        'data/{sample}.unrelatedpruned.bed',
        'data/{sample}.unrelatedpruned.bim',
        'data/{sample}.unrelatedpruned.fam',
        'data/{sample}.unrelatedpruned.eigenvec',
        'data/{sample}.unrelatedpruned.eigenval'
    output:
        temp('data/{sample}.pcl')
    shell:
        'gcta64 --bfile data/{wildcards.sample}.unrelatedpruned --pc-loading data/{wildcards.sample}.unrelatedpruned --out data/{wildcards.sample}'

#This section 'often manual' - can adjust the percentile to the config file to allow user to tweak settings
#Find & Remove SNPs with high PC loading PC1-3 (by extracting only the other SNPs)
rule find_highPC_SNPs:
    input:
        'data/{sample}.pcl'
    output:
        'data/{sample}.lowPCSNPs'
    params:
        PCLthreshold = config['pcl_threshold']
    script:
        'scripts/PCfilterSnakemake.R'

rule exclude_highPC_SNPs:
    input:
        'data/{sample}.prekinship.bed',
        'data/{sample}.prekinship.bim',
        'data/{sample}.prekinship.fam',
        'data/{sample}.lowPCSNPs'
    output:
        temp('data/{sample}.highPCexcluded.bed'),
        temp('data/{sample}.highPCexcluded.bim'),
        temp('data/{sample}.highPCexcluded.fam')
    shell:
        'plink --bfile data/{wildcards.sample}.prekinship --extract data/{wildcards.sample}.lowPCSNPs --make-bed --out data/{wildcards.sample}.highPCexcluded'

#Kinship testing 2 (for final list of unrelated individuals)
rule kinship_testing_round2:
    input:
        'data/{sample}.highPCexcluded.bed',
        'data/{sample}.highPCexcluded.bim',
        'data/{sample}.highPCexcluded.fam'
    output:
        temp('data/{sample}.round2.kin0')
    shell:
        'king -b data/{wildcards.sample}.highPCexcluded.bed --kinship --degree 3 --prefix data/{wildcards.sample}.round2'

#calculate largest independent set of unrelated individuals again
rule calculate_largest_independent_unrelated_group_round2:
    input:
        'data/{sample}.round2.kin0',
        'data/{sample}.highPCexcluded.fam'
    output:
        'data/{sample}.related_d3_round2.txt'
    script:
        'scripts/FindRelatedIndividualsSnakemake.py'

#DOWNLOAD AND PROCESS 1000GP DATA - avoid these steps (if you've already got the 1000GP plink files)
#by having a folder 1000GP/ with files chr1-22.bim / bed / fam in the folder with snakefile.

#Define Chromosomes as wildcard
rule get_chr:
    input: 
        expand('1000GP/CCDG_14151_B01_GRM_WGS_2020-08-05_chr{i}.filtered.shapeit2-duohmm-phased.bed', i = range(1,23)), #(it doesn't include the final number in the range)
        expand('1000GP/CCDG_14151_B01_GRM_WGS_2020-08-05_chr{i}.filtered.shapeit2-duohmm-phased.bim', i = range(1,23)), #bim and fam added for the chromosome merger step
        expand('1000GP/CCDG_14151_B01_GRM_WGS_2020-08-05_chr{i}.filtered.shapeit2-duohmm-phased.fam', i = range(1,23)) 

#Download 1000GP files as vcf.gz - these are the GRCh38 build files
rule download_1000GP_vcf:
    output: 
        '1000GP/CCDG_14151_B01_GRM_WGS_2020-08-05_chr{chr}.filtered.shapeit2-duohmm-phased.vcf.gz'
    shell:
        'wget -P 1000GP/ ftp://ftp.1000genomes.ebi.ac.uk/vol1/ftp/data_collections/1000G_2504_high_coverage/working/20201028_3202_phased//CCDG_14151_B01_GRM_WGS_2020-08-05_chr{wildcards.chr}.filtered.shapeit2-duohmm-phased.vcf.gz ftp://ftp.1000genomes.ebi.ac.uk/vol1/ftp/data_collections/1000G_2504_high_coverage/working/20201028_3202_phased//CCDG_14151_B01_GRM_WGS_2020-08-05_chr{wildcards.chr}.filtered.shapeit2-duohmm-phased.vcf.gz.tbi'

#Download 1000GP ped file
rule download_1000GP_ped:
    output:
        '1000GP/20130606_g1k.ped'
    shell:
        'wget -P 1000GP/ ftp://ftp.1000genomes.ebi.ac.uk/vol1/ftp/technical/working/20130606_sample_info/20130606_g1k.ped'

#Download the GRCh38 reference genome
rule download_GRCh38_reference:
    output: 
        '1000GP/Homo_sapiens.GRCh38.dna.chromosome.{chr}.fa'
    shell: #the sed command converts chromosomes from 1,2,3 etc. to chr1, chr2, chr3 for consistency with the vcf file
        r'''
        wget -P 1000GP/ ftp://ftp.ensembl.org/pub/release-88/fasta/homo_sapiens/dna/Homo_sapiens.GRCh38.dna.chromosome.{wildcards.chr}.fa.gz
        gunzip {output}.gz
        sed -i -e 's/^>/>chr/' {output}
        '''

#Convert the 1000 Genomes files to BCF
rule convert_to_BCF:
    input:
        fa = '1000GP/Homo_sapiens.GRCh38.dna.chromosome.{chr}.fa',
        vcf = '1000GP/CCDG_14151_B01_GRM_WGS_2020-08-05_chr{chr}.filtered.shapeit2-duohmm-phased.vcf.gz'
    output:
        temp('1000GP/CCDG_14151_B01_GRM_WGS_2020-08-05_chr{chr}.filtered.shapeit2-duohmm-phased.bcf')
    shell:
        r'''
        bcftools norm -m-any --check-ref w -f {input.fa} {input.vcf} | \
        bcftools annotate -x ID -I +'%CHROM:%POS:%REF:%ALT' | \
        bcftools norm -Ob --rm-dup both > {output} ; bcftools index {output}
        '''

#Convert the BCF files to PLINK format
rule convert_to_PLINK:
    input:
        '1000GP/CCDG_14151_B01_GRM_WGS_2020-08-05_chr{chr}.filtered.shapeit2-duohmm-phased.bcf'
    output:
        temp('1000GP/CCDG_14151_B01_GRM_WGS_2020-08-05_chr{chr}.filtered.shapeit2-duohmm-phased.bed'),
        temp('1000GP/CCDG_14151_B01_GRM_WGS_2020-08-05_chr{chr}.filtered.shapeit2-duohmm-phased.bim'),
        temp('1000GP/CCDG_14151_B01_GRM_WGS_2020-08-05_chr{chr}.filtered.shapeit2-duohmm-phased.fam')
    shell:
        'plink --bcf {input} --keep-allele-order --vcf-idspace-to _ --const-fid --allow-extra-chr 0 --split-x b38 no-fail --make-bed --out 1000GP/CCDG_14151_B01_GRM_WGS_2020-08-05_chr{wildcards.chr}.filtered.shapeit2-duohmm-phased'

#Merge all chromosomes to single file
rule list_chromosome_files:
    input:
        expand('1000GP/CCDG_14151_B01_GRM_WGS_2020-08-05_chr{i}.filtered.shapeit2-duohmm-phased.bed', i = range(1,23)),
        expand('1000GP/CCDG_14151_B01_GRM_WGS_2020-08-05_chr{i}.filtered.shapeit2-duohmm-phased.bim', i = range(1,23)),
        expand('1000GP/CCDG_14151_B01_GRM_WGS_2020-08-05_chr{i}.filtered.shapeit2-duohmm-phased.fam', i = range(1,23))
    output:
        temp('1000GP/chromosome_list')
    shell:
        r'''
        for chr in {{1..22}}; do
            echo 1000GP/CCDG_14151_B01_GRM_WGS_2020-08-05_chr${{chr}}.filtered.shapeit2-duohmm-phased.bed 1000GP/CCDG_14151_B01_GRM_WGS_2020-08-05_chr${{chr}}.filtered.shapeit2-duohmm-phased.bim 1000GP/CCDG_14151_B01_GRM_WGS_2020-08-05_chr${{chr}}.filtered.shapeit2-duohmm-phased.fam >> 1000GP/chromosome_list;
        done
        '''

rule merge_chromosomes: #needs quite a bit of memory, probably requires cluster.
    input:
        '1000GP/chromosome_list',
        expand('1000GP/CCDG_14151_B01_GRM_WGS_2020-08-05_chr{i}.filtered.shapeit2-duohmm-phased.bed', i = range(1,23)),
        expand('1000GP/CCDG_14151_B01_GRM_WGS_2020-08-05_chr{i}.filtered.shapeit2-duohmm-phased.bim', i = range(1,23)),
        expand('1000GP/CCDG_14151_B01_GRM_WGS_2020-08-05_chr{i}.filtered.shapeit2-duohmm-phased.fam', i = range(1,23))
    output:
        '1000GP/chr1-22.bed',
        '1000GP/chr1-22.bim',
        '1000GP/chr1-22.fam'
    shell:
        'plink --merge-list 1000GP/chromosome_list --make-bed --out 1000GP/chr1-22'

#remove the 698 related individuals from 1000GP dataset
rule remove_related_1000GP:
    input:
        '1000GP/chr1-22.bed',
        '1000GP/chr1-22.bim',
        '1000GP/chr1-22.fam',
        'files/1000GPrelated.txt'
    output:
        temp('1000GP/chr1-22unrelated.bed'),
        temp('1000GP/chr1-22unrelated.bim'),
        temp('1000GP/chr1-22unrelated.fam')
    shell:
        'plink --bfile 1000GP/chr1-22 --remove files/1000GPrelated.txt --make-bed --out 1000GP/chr1-22unrelated'

rule fix_chr_codes:
    input:
        '1000GP/chr1-22unrelated.bed',
        '1000GP/chr1-22unrelated.bim',
        '1000GP/chr1-22unrelated.fam'
    output:
        '1000GP/chr1-22-tomerge.bim',
        '1000GP/chr1-22-tomerge.bed',
        '1000GP/chr1-22-tomerge.fam'
    shell:
        r'''
        sed 's/chr//' 1000GP/chr1-22unrelated.bim > 1000GP/chr1-22-tomerge.bim
        cp 1000GP/chr1-22unrelated.bed 1000GP/chr1-22-tomerge.bed
        cp 1000GP/chr1-22unrelated.fam 1000GP/chr1-22-tomerge.fam
        '''

#Merge back with own dataset

#Change IDs to 1000GP ids and identify SNPs requiring strand-flipping
rule convertkgp1:
    input:
        'data/{sample}.prekinship.bim',
        '1000GP/chr1-22-tomerge.bim'
    output:
        temp('data/{sample}_snps_toflip.txt')
    script:
        'scripts/convertkgp1.R'

#flip SNPs
rule flip_SNPs:
    input:
        'data/{sample}.prekinship.bed',
        'data/{sample}.prekinship.bim',
        'data/{sample}.prekinship.fam',
        'data/{sample}_snps_toflip.txt'
    output:
        temp('data/{sample}.flip.bed'),
        temp('data/{sample}.flip.bim'),
        temp('data/{sample}.flip.fam'),
    shell:
        'plink --bfile data/{wildcards.sample}.prekinship --flip data/{wildcards.sample}_snps_toflip.txt --make-bed --out data/{wildcards.sample}.flip'

rule convertkgp2:
    input:
        'data/{sample}.flip.bim',
        '1000GP/chr1-22-tomerge.bim'
    output:
        temp('data/{sample}.Mappingtokgp'),
        temp('data/{sample}.SNPsCommonKGP'),
        temp('data/{sample}.SNPsdup'),
        temp('data/{sample}.KGP_mapping.bim')
    script:
        'scripts/convertkgp2.R'

#create new bed/bim/fam for study data with 1000GP mapping & snps, excluding duplicates
#This is a rare case where we can use the {input.xyz} because we are using the .bed/bim/fam as independent inputs rather than having plink assume them from a prefix
rule map_to_1000GP:
    input:
        bed = 'data/{sample}.flip.bed',
        bim = 'data/{sample}.KGP_mapping.bim',
        fam = 'data/{sample}.flip.fam',
        common = 'data/{sample}.SNPsCommonKGP',
        dup = 'data/{sample}.SNPsdup'
    output:
        temp('data/{sample}.1000GPSNPs.bed'),
        temp('data/{sample}.1000GPSNPs.bim'),
        temp('data/{sample}.1000GPSNPs.fam')
    shell:
        'plink --bed {input.bed} --bim {input.bim} --fam {input.fam} --extract {input.common} --exclude {input.dup} --make-bed --out data/{wildcards.sample}.1000GPSNPs'

#filter 1000GP data for SNPs common to study data, removing duplicates
rule filter_1000GP:
    input:
        '1000GP/chr1-22-tomerge.bed',
        '1000GP/chr1-22-tomerge.bim',
        '1000GP/chr1-22-tomerge.fam',
        'data/{sample}.SNPsCommonKGP',
        'data/{sample}.SNPsdup'
    output:
        temp('1000GP/{sample}.1000GPCommonSNPs.bed'),
        temp('1000GP/{sample}.1000GPCommonSNPs.bim'),
        temp('1000GP/{sample}.1000GPCommonSNPs.fam'),
    shell:
        'plink --bfile 1000GP/chr1-22-tomerge --extract data/{wildcards.sample}.SNPsCommonKGP --exclude data/{wildcards.sample}.SNPsdup --out 1000GP/{wildcards.sample}.1000GPCommonSNPs --make-bed'

rule merge_with_1000GP:
    input:
        'data/{sample}.1000GPSNPs.bed',
        'data/{sample}.1000GPSNPs.bim',
        'data/{sample}.1000GPSNPs.fam',
        '1000GP/{sample}.1000GPCommonSNPs.bed',
        '1000GP/{sample}.1000GPCommonSNPs.bim',
        '1000GP/{sample}.1000GPCommonSNPs.fam'
    output:
        'data/{sample}.mergedwith1000GP.bed',
        'data/{sample}.mergedwith1000GP.bim',
        'data/{sample}.mergedwith1000GP.fam'
    shell:
        'plink --bfile data/{wildcards.sample}.1000GPSNPs --bmerge 1000GP/{wildcards.sample}.1000GPCommonSNPs --out data/{wildcards.sample}.mergedwith1000GP'

#pruning for SNPs in LD - 1000 50 0.05
rule prune_LD_SNPs_2:
    input:
        'data/{sample}.mergedwith1000GP.bed',
        'data/{sample}.mergedwith1000GP.bim',
        'data/{sample}.mergedwith1000GP.fam'
    output:
        temp('data/{sample}-1000GP.prune.out'),   
        temp('data/{sample}-1000GP.prune.in')
    params:
        ipvalue2 = config['indep_pairwise2']
    shell:
        'plink --bfile data/{wildcards.sample}.mergedwith1000GP --indep-pairwise {params.ipvalue2} --out data/{wildcards.sample}-1000GP'

#remove pruned LD SNPs
rule remove_pruned_round2:
    input:
        'data/{sample}.mergedwith1000GP.bed',
        'data/{sample}.mergedwith1000GP.bim',
        'data/{sample}.mergedwith1000GP.fam',
        'data/{sample}-1000GP.prune.out',   
        'data/{sample}-1000GP.prune.in'
    output:
        temp('data/{sample}.mergedpruned.bed'),
        temp('data/{sample}.mergedpruned.bim'),
        temp('data/{sample}.mergedpruned.fam')
    shell:
        'plink --bfile data/{wildcards.sample}.mergedwith1000GP --extract data/{wildcards.sample}-1000GP.prune.in --make-bed --out data/{wildcards.sample}.mergedpruned' 

#remove related individuals from the pre-pca files
rule remove_related_individuals_round2:
    input:
        'data/{sample}.mergedpruned.bed',
        'data/{sample}.mergedpruned.bim',
        'data/{sample}.mergedpruned.fam',
        'data/{sample}.related_d3_round2.txt'
    output:
        temp('data/{sample}.unrelated_merge.bed'),
        temp('data/{sample}.unrelated_merge.bim'),
        temp('data/{sample}.unrelated_merge.fam')
    shell:       
        'plink --bfile data/{wildcards.sample}.mergedpruned --remove data/{wildcards.sample}.related_d3_round2.txt --make-bed --out data/{wildcards.sample}.unrelated_merge'

#calculate pca 
rule pca2:
    input:
        'data/{sample}.unrelated_merge.bed',
        'data/{sample}.unrelated_merge.bim',
        'data/{sample}.unrelated_merge.fam'
    output:
        'data/{sample}.unrelated_merge.eigenvec',
        'data/{sample}.unrelated_merge.eigenval'
    params:
        PCnumber = config['pc_number']
    shell:
        'plink --bfile data/{wildcards.sample}.unrelated_merge --pca {params.PCnumber} --out data/{wildcards.sample}.unrelated_merge'

#ANCESTRY ANALYSIS

#obtain pop file for superpopulations in 1000 genomes
rule get_kgpop:
    output:
        '1000GP/integrated_call_samples_v3.20130502.ALL.panel'
    shell:
        'wget -P 1000GP/ ftp://ftp.1000genomes.ebi.ac.uk/vol1/ftp/release/20130502/integrated_call_samples_v3.20130502.ALL.panel'

rule create_superpop_file:
    input:
        '1000GP/integrated_call_samples_v3.20130502.ALL.panel'
    output:
        '1000GP/kgp.superpop'
    shell:
        r'''awk 'BEGIN {{print "FID\tIID\tPOP"}} NR>1 {{print "0\t"$1"\t"$3}}' {input} > {output}'''

rule create_pop:
    input:
        'data/{sample}.mergedpruned.fam',
        '1000GP/kgp.superpop'
    output:
        temp('data/{sample}.mergedpruned.pop')
    script:
        'scripts/createpopSnakemake.R'

rule admixture:
    input:
        bed = 'data/{sample}.mergedpruned.bed',
        bim = 'data/{sample}.mergedpruned.bim',
        fam = 'data/{sample}.mergedpruned.fam',
        popfile = 'data/{sample}.mergedpruned.pop'
    output:
        temp('{sample}.mergedpruned.5.Q'),
        temp('{sample}.mergedpruned.5.P')
    shell:
        'admixture {input.bed} 5 --supervised -C 0.00001' # 5 (here and in output files) is the number of populations

rule estimate_ancestry:
    input:
        '1000GP/kgp.superpop',
        '{sample}.mergedpruned.5.Q',
        'data/{sample}.mergedpruned.fam',
        'data/{sample}.mergedpruned.pop'
    output:
        'data/{sample}.estimatedancestry.txt'
    script:
        'scripts/estimateancestrySnakemake.R'

#placeholder - ancestry analysis (more complicated for global data)

#IMPUTATION STEPS

#add instance for accessing Topmed API (once per user/machine)
rule add_API_instance:
    params:
        APIkey = config['API']
    shell:
        r'''
        echo -e 'https://imputation.biodatacatalyst.nhlbi.nih.gov/\n{params.APIkey}' | ./imputationbot add-instance
        '''

#Pre-Imputation QC - Filter for flipped SNPs and others that can't be read by TopMed, and apply maf/hwe
#Note that the lists of SNPs for flip/removal were based on TopMed failure with our data. They may vary with other datasets.
#In which case make these files empty, attempt imputation, and generate new lists based on TopMed's failure report.
rule imputation_QC:
    input:
       'data/{sample}.PostBasicQC.bed',
       'data/{sample}.PostBasicQC.bim',
       'data/{sample}.PostBasicQC.fam',
       'files/flipthesesnps.txt',
       'files/removethesesnps.txt'
    output:
       temp('data/{sample}.preimputation.bed'),
       temp('data/{sample}.preimputation.bim'),
       temp('data/{sample}.preimputation.fam')
    shell:
        'plink --bfile data/{wildcards.sample}.PostBasicQC --maf 0.01 --hwe 0.006 --flip files/flipthesesnps.txt --exclude files/removethesesnps.txt --make-bed --out data/{wildcards.sample}.preimputation'

#create variable for chromosome number for use in the convert_to_VCF rule
rule get_chr_again:
    input:
        expand('data/{sample}vcfs/{i}.vcf.gz', sample = SAMPLES, i = range(1,23))

#split by chromosome and convert to vcf for TopMed submission
rule convert_to_VCF:
    input:
       'data/{sample}.preimputation.bed',
       'data/{sample}.preimputation.bim',
       'data/{sample}.preimputation.fam'
    output: 
        'data/{sample}vcfs/{i}.vcf.gz'
    shell:
        r'''
        plink --bfile data/{wildcards.sample}.preimputation --recode vcf --chr {wildcards.i} --output-chr chr26 --out data/{wildcards.sample}vcfs/{wildcards.i}
        bgzip data/{wildcards.sample}vcfs/{wildcards.i}.vcf
        '''

#submit for imputation with imputation bot
checkpoint submit_for_imputation:
    input:
        expand('data/{sample}vcfs/{i}.vcf.gz', sample = SAMPLES, i = range(1,23))
    output:
        directory('data/{sample}_imputation')
    params:
        password = config['imputation_password']
    shell:
        r'''
        ./imputationbot impute --files data/{wildcards.sample}vcfs/*.vcf.gz --refpanel topmed-r2 --population all --build hg38 --autoDownload --name {wildcards.sample} --password {params.password}
        mv job-*{wildcards.sample} data/{wildcards.sample}_imputation
        '''

def get_imputed_files(wildcards):
    checkpoint_output = checkpoints.submit_for_imputation.get(**wildcards).output[0]
    return expand("data/{sample}/imputationcleaned/chr{i}.cleaned.vcf.gz", sample = SAMPLES, i = glob_wildcards(os.path.join(checkpoint_output, "/local/chr{i}.dose.vcf")).i, allow_missing=True)

#Post-imputation QC
#remove low imputation score SNPs (R2 < 0.4) and monoallelic SNPs with bcftools
rule remove_monoallelic_and_low_score:
    input:
        get_imputed_files
    output:
        temp('data/{sample}/imputationcleaned/chr{i}.cleaned.vcf.gz')
    shell:
        r'''bcftools view -m2 -e '( R2 < 0.4)' data/{wildcards.sample}_imputation/local/chr{wildcards.i}.dose.vcf.gz -o {output}'''

#convert to bgen format {remove maf later - remove monoallelic snps using bcftools}
rule convert_to_pgen:
    input:
        'data/{sample}/imputationcleaned/chr{i}.cleaned.vcf.gz'
    output:
        'data/{sample}/pgens/chr{i}.pgen',
        'data/{sample}/pgens/chr{i}.pvar',
        'data/{sample}/pgens/chr{i}.psam'
    shell:
        'plink2 --vcf {input} --id-delim --make-pgen --out data/{wildcards.sample}/pgens/chr{wildcards.i}'

#ANALYSIS

#Prepare data (covariates etc.) for input to regenie
#Merge ancestry information with pc values and other information (age, sex, phenotype, etc.)
rule merge_data_files:
    input:
        'data/{sample}.estimatedancestry.txt',
        'data/{sample}.info'
    output:
        'data/{sample}.alldata.csv'
    script:
        'scripts/buildmerge_snakemake.py'

#define groups from config (& convert from string to list)
import ast
grouptitles = config['grouping_values']
GROUPS = ast.literal_eval(grouptitles)

rule define_groups:
    input:
        expand('data/{sample}.{group}.idlist', sample = SAMPLES, group = GROUPS)

#clean data, group by variable (eg. pop), extract phenotype info to .pheno file per group, extract covariates to .covar file per group.
rule prep_for_analysis:
    input:
        'data/{sample}.alldata.csv'
    output:
        'data/{sample}.{group}.idlist',
        'data/{sample}.{group}.pheno',
        'data/{sample}.{group}.covar',
        'data/{sample}.{group}.stats'
    script:
        'scripts/make_pheno_covar_snakemake.py'

#split plink files by group
rule split_plink_by_group:
    input:
       'data/{sample}.PostBasicQC.bed',
       'data/{sample}.PostBasicQC.bim',
       'data/{sample}.PostBasicQC.fam',
       'data/{sample}.{group}.idlist'
    output:
        'data/{sample}/{group}.bed',
        'data/{sample}/{group}.bim',
        'data/{sample}/{group}.fam'
    wildcard_constraints:
        group = '\w+'
    shell:
        r'''
        plink --bfile data/{wildcards.sample}.PostBasicQC --keep data/{wildcards.sample}.{wildcards.group}.idlist --make-bed --out data/{wildcards.sample}/{wildcards.group}
        plink --bfile data/{wildcards.sample}/{wildcards.group} --maf 0.01 --make-bed --out data/{wildcards.sample}/{wildcards.group}
        '''
#Calculate PCs by specific group
rule pca_by_group:
    input:
        'data/{sample}/{group}.bed',
        'data/{sample}/{group}.bim',
        'data/{sample}/{group}.fam'
    output:
        'data/{sample}/{group}.eigenvec',
        'data/{sample}/{group}.eigenval'
    wildcard_constraints:
        group = '\w+'
    shell:
        'plink --bfile data/{wildcards.sample}/{wildcards.group} --pca 10 --out data/{wildcards.sample}/{wildcards.group}'

#Add the new PCs to the covar file
rule add_pcs_to_covar:
    input:
        'data/{sample}.{group}.covar',
        'data/{sample}/{group}.eigenvec'
    output:
        'data/{sample}.{group}.covarpcs'
    script:
        'scripts/add_pc_covar_snakemake.py'

#define categorical covariates
if config['cat_columns'] == 'NULL':
    catcovarstring = ''
else:
    catcovars = ast.literal_eval(config['cat_columns'])
    catcovarstring = '--catCovarList '
for covar in catcovars:
    catcovarstring += covar + ', '


#association analysis in regenie
rule regenie_step_1:
    conda: 'files/regenie_env.yaml'
    input:
        bed = 'data/{sample}/{group}.bed',
        bim = 'data/{sample}/{group}.bim',
        fam = 'data/{sample}/{group}.fam',
        pheno = 'data/{sample}.{group}.pheno',
        covar = 'data/{sample}.{group}.covarpcs'
    output:
        'data/{sample}/{group}.regenie_1.loco', #temp these eventually once debugged - also check that it outputs regenie_1 (index) not regenie_pheno
        'data/{sample}/{group}.regenie_pred.list'
    shell:
        'regenie --step 1 --bed data/{wildcards.sample}/{wildcards.group} --phenoFile data/{wildcards.sample}.{wildcards.group}.pheno --covarFile data/{wildcards.sample}.{wildcards.group}.covarpcs {catcovarstring}--bsize 100 --bt --out data/{wildcards.sample}/{wildcards.group}.regenie'

#convert imputation pgens to plink to get the bims for input to batch script
rule pgen_to_bim:
    input:
        'data/{sample}/pgens/chr{i}.pvar'
    output:
        'data/{sample}/bims/chr{i}.bim',
    shell:
        'plink2 --pvar data/{wildcards.sample}/pgens/chr{wildcards.i}.pvar --make-just-bim --out data/{wildcards.sample}/bims/chr{wildcards.i}'

#define pheno (for regenie outputs):
pheno = config['case_control']    

#Generate batches from bim file and run each batch on regenie via script; output to directory per chromosome
checkpoint generate_batches:
    input:
        'data/{sample}/bims/chr{i}.bim'
    output:
        directory('data/{sample}/chr{i}_batches/')
    script:
        'scripts/generate_batches.py'

rule regenie_on_batches:
    conda: 'files/regenie_env.yaml'
    input:
        'data/{sample}/chr{i}_batches/{x}.batch',
        'data/{sample}/pgens/chr{i}.pgen', #does the pgen need to be divided up by group as well?
        pred = 'data/{sample}/{group}.regenie_pred.list',
        phenofile = 'data/{sample}.{group}.pheno',
        covar = 'data/{sample}.{group}.covarpcs',
    params:
        cat = catcovarstring.strip(', ')
    output:
        'data/{sample}/{group}.results/chr{i}/{x}_{pheno}.regenie'
    shell:
        'regenie --step 2 --pgen data/{wildcards.sample}/pgens/chr{wildcards.i} --phenoFile {input.phenofile} --covarFile {input.covar} {catcovarstring}--bt --firth --approx --pThresh 0.01 --pred {input.pred} --bsize 400 --out data/{wildcards.sample}/{wildcards.group}.results/chr{wildcards.i}/{wildcards.x}'
        #if need to swap to bgen later- 'regenie --step 2 --bgen data/{wildcards.sample}/bgens/chr{wildcards.i}.bgen --ref-first --sample data/{wildcards.sample}/bgens/chr{wildcards.i}.sample --phenoFile {input.phenofile} --covarFile {input.covar} {catcovarstring}--bt --firth --approx --pThresh 0.01 --pred {input.pred} --bsize 400 --out data/{wildcards.sample}/{wildcards.group}.results/chr{wildcards.i}/{wildcards.batch}'

#function to make wildcard {x} (batch) out of the batch files from the output of the checkpoint rule above
def get_batches(wildcards):
    checkpoint_output = checkpoints.generate_batches.get(**wildcards).output[0]
    return expand("data/{sample}/{group}.results/chr{{i}}/{x}_{pheno}.regenie", sample = SAMPLES, group = GROUPS, i = range(1,23), pheno = pheno, x = glob_wildcards(os.path.join(checkpoint_output, "{x}.batch")).x, allow_missing=True)

rule list_results_files:
    input:
        get_batches
    output:
        'data/{sample}/{group}.results/chr{i}/results_list'
    shell:
        r'''for batch in {input}; do
                echo ${{batch}} >> {output}
            done'''

rule collate_results_lists:
    input:
        expand('data/{sample}/{group}.results/chr{i}/results_list', sample = SAMPLES, group = GROUPS, i = range(1,23))
    output:
        'data/{sample}/{group}.results/results_list'
    shell:
        r'''for file in {input}; do
                cat ${{file}} >> {output}
            done'''


#plot results from regenie as Manhattan plots (all chromosomes / bins merged per sample)
rule plot_results:
    input:
        'data/{sample}/{group}.results/results_list'
    output:
        'data/{sample}/{group}.results/manhattan.png'
    script:
        'scripts/plot_manhattan_snakemake.py'

#replace 'regenie analysis' with this once pipeline confirmed working; it should still run appropriately
rule get_results:
    input:
        expand('results/{sample}-results_Markdown.html', sample=SAMPLES)

#Create list of groups for collating results per sample
rule collate_groups:
    input:
        png = expand('data/{sample}/{group}.results/manhattan.png', sample = SAMPLES, group = GROUPS),
        stats = expand('data/{sample}.{group}.stats', sample = SAMPLES, group = GROUPS)
    output:
        manhattans = 'data/{sample}/groups_list',
        stats = 'data/{sample}/stats'
    shell:
        r'''
        for group in {{GROUPS}}; do
            echo {input.png} >> {output.manhattans}
            echo {input.stats} >> {output.stats}
        done
        '''

#target rule for analysis & generation of Markdown doc
rule analysis:
    input:
        'data/{sample}/stats',
        'data/{sample}/groups_list'
    output:
        'results/{sample}-results_Markdown.html'
    script:
        'scripts/AnalysisMarkdown_Snakemake.Rmd'
    