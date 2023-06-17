# GWAS-workflow-Snakemake

Snakemake automation of GWAS workflows (QC, ancestry, imputation, association analysis)

**Caveat Emptor:** This is a work in progress with design and instructions to run on Baillie Lab data using the University of Edinburgh computing cluster (Eddie). It comes (for now) with no user guarantee or implementation instructions for other systems.

**Prerequisites:**
- TopMed account
- Anaconda installed
- around 12GB RAM for most steps, but large amounts required a few specific steps (potentially more for very large sample groups):
            - the merge of chromosomes from the 1000GP files (avoided by using pre-prepared files) needs 80-100GB RAM.
            - the combination of all results to generate a manhattan plot needs 80-100GB RAM.
            - some other steps, such as the ancestry/admixture calculations and merger of 1000GP and test datasets, probably require 48-64GB RAM. 
- Enough local (to the snakefile/data) storage space for 1000GP genomes and imputed genomes from your data (eg. around 500GB)
- Adjust config file (files/gwas_config.yaml) to specify labelling conventions, phenotype codes, covariates, etc. for your own dataset.
            - please read the comments on config; some values require specific input (eg. 'NULL' as input if value should be empty)

To run on Eddie (needs at least 64GB, regardless of sample size, unless using the pre-prepared 1000GP files):

1. **have data (ped/map files in Hg38 build) in a folder 'data' in the same directory as snakefile** and related files/folders from GitHub repo.
      - the workflow should run on every ped/map in the folder, so you can run chromosomes separately or run multiple GWAS, in theory.
2. 'module load anaconda/5.3.1' in your qlogin session
3. 'source activate workflow'* (if it's the first time you're using, you have to build 'workflow' env from the workflow.yaml, instructions below)

3.5. If building environment; 'conda env create -n nameofenviroment -f files/workflow-env.yaml' . May take some time to solve/install everything.

4. Snakemake runs:

      -```snakemake basicqc -c*``` (* = number of cores to use, eg. 2) - does the basic qc and creates a markdown .html file tabulating/graphing filtering steps.
      
      -```snakemake pca_check -c*``` - does the pca qc and kinship analysis, creates a markdown .html file of filtering steps and PCs
      
      -```snakemake ancestry_check -c*``` - does ancestry analysis and creates text file
      
      4.5 - ```snakemake add_API_instance -c*``` - do this once per user (eg. on your Eddie login account) to set up the API for TopMed. Needs API key in config.
      
      -```snakemake topmed_imputation -c*``` - submits to TopMed for imputation and automatically downloads results
      
      -```snakemake get_results -c* --use-conda --conda-frontend conda``` - performs association analysis in regenie. --use-conda because regenie needs its                                                                              own environment, handled by snakemake. --conda-frontend conda needed for running
                                                                             on Eddie.
            
      **this step requires an additional .info file in the 'data' folder with the same name as the map/ped (ie. data/{filename}.info), containing phenotype                       and covariate / other info for the samples (eg. Age, Sex, etc.). Must be in CSV format, with missing values either empty or whitespaced.**
    
These steps are generally consecutive - if you run one (eg. ancestry_check) without running the earlier ones, it will still do necessary steps from the earlier
part of the pipeline, but it won't produce the Markdown files.

QC/PCA Workflow runs in around 1 hour on 2 32GB cores on Eddie with a 5000-person dataset, assuming 1000GP data already downloaded/processed (I have a folder containing these files to avoid download). The workflow will skip these download steps if the end-files are present, unless --force-all flag is used.

Ancestry step takes around 1-2 hours extra.

Topmed submission steps can be flaky due to server issues. If the workflow fails/exits, check the job status on the Topmed site. If the job completes successfully, you can manually download the files and re-enter the workflow pipeline by organising them as follows from the working directory containing the snakefile:

      data > {yoursampleID}_imputation > local > chr1-22 files (dose.vcf.gz, info.gz,etc.)

This should be achievable by downloading from Topmed and copying the unzipped files to a folder as described.

**Feeling brave?**
```snakemake -c* --use-conda --conda-frontend conda``` (ie. no target rule specified) should run everything. Topmed submission/download step can take 7-10 hours (or more with larger datasets), so recommended as an overnight run. Post-imputation steps also take several hours due to large files. The final regenie step is optimised for running on many small cores in parallel, but could take many days if run on a single core.

**Imputationbot Documentation**
The Topmed Imputationbot has been repackaged here for ease of use. Documentation / original source for the imputationbot software can be found here https://github.com/lukfor/imputationbot or here https://imputationbot.readthedocs.io/en/latest/
