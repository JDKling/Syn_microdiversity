# Syn_microdiversity

This pipeline was designed to run on USC's HPC. Therefore some steps may have to be optimized for your specific computational environment.

### Before running the pipeline
Really you will need two things to get this running: 1) Anaconda and 2) Snakemake. Here you can find instructions for downloading [anaconda](https://docs.anaconda.com/anaconda/install/linux/) and [snakemake](https://snakemake.readthedocs.io/en/stable/getting_started/installation.html).

Many HPCs using slurm do not let you access the internet within a cluster. This mostly fine, it does mean that steps such as downloading raw data or databases needs to be done in the head node. I've tried to make this easy to do with a simple batch script.

I would use a screen for this because some of the databases take a while to download.
```{bash}
 bash SCRIPTS/DL_data.sh
```

In addition to making sure all of the reads for the samples listed in the sra_accessions.txt file are downloaded into DATA/RAW/ it will also install any conda environments that will need to be run outside of the snakemake framework (e.g. the epigenomics pipeline here).

It will also make a directory for storing the error logs for each step.

### Running snakemake
Once all of the raw data is present and each database has been downloaded, this pipeline is ready to fly. Most of the functionality (read qc, digital normalization, genome assembly, SNP profiling, and community analysis) is run inside of the snakemake workflow management system. This tool is really great for managing large datasets where running steps in parallel is important. It also integrates well into larger computational frameworks (like academic HPC's).

With snakemake installed, you can run majority of the pipeline like this.
```{bash}
conda activate snakemake
mkdir -p LOGS/SLURM_LOG
snakemake -s snakefile --cluster-config cluster.yml --use-conda --max-jobs-per-second 10 --cluster 'sbatch --output=LOGS/SLURM_LOG/%j.out --error=LOGS/SLURM_LOG/%j.out -t {cluster.time} --mem={cluster.mem} -c {cluster.cpus}' -j 10 --rerun-incomplete
conda deactivate
```

### Epigenomics
INSTRUCTIONS COMING SOON 
