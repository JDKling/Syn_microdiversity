# Syn_microdiversity

This pipeline was designed to run on USCs HPC. Therefore some steps may have to be optimized for your specific computational environment.





### Before running the pipeline
Really you will need two things to get this running: 1) Anaconda and 2) Snakemake. Here you can find instructions for downloading anaconda and snakemake. HYPERLINK THESE LATER

Many HPCs using slurm do not let you access the internet within a cluster. This mostly fine, it does mean that steps such as downloading raw data or databases needs to be done in the head node. I've tried to make this easy to do with a simple batch script.

I would use a screen for this because some of the databases are quite large.
```{bash}
 SCRIPTS/DL_data.sh
```
