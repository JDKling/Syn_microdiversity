#!/usr/bin/env bash
mkdir -p DBs

# download sequence reads
## create and environment that will let us download the raw seq data we are interested in
if conda env list | grep -q 'dl_data'; then
  echo "dl_data environment is installed, checking to see if we have our raw data"
else
  echo "Setting up the dl_data environment"
  conda create -y -f CONDA-ENVS/env0.yml -n dl_data
  conda activate dl_data
  conda install -y -c conda-forge -c bioconda -c defaults -c astrobiomike entrez-direct sra-tools bit
fi

if ls RAW/*fq.gz
mkdir RAW
cd RAW
SAMN=$(grep 'short_read' sra_accessions.txt | cut -f2)
STRAIN=$(grep 'short_read' sra_accessions.txt | cut -f1)
#ONT_SAMN=
#ONT_STRAIN=
REF_ACC=$(cut -f2 refseq_accessions.txt)
REF_STRAIN=$(cut -f1 refseq_accessions.txt)

## Download short reads
count=0
for acc in ${SAMN[@]}
do
  if [ ! -f "${STRAIN[$count]}_1.fq.gz" ]
  then
      echo "$0: File '${STRAIN[$count]}' not found."
      echo 'now downloading' ${STRAIN[$count]}
      esearch -db sra -query $acc  | efetch --format runinfo | cut -d ',' -f 1 | grep SRR | xargs -n 1 -P 12 fastq-dump --split-files --gzip
      mv SRR*_1.fastq.gz ${STRAIN[$count]}_1.fq.gz
      mv SRR*_2.fastq.gz ${STRAIN[$count]}_2.fq.gz
  fi
count=$(expr $count + 1)
done

## Download ONT reads

### FILL IN LATER!!!!

##Download reference genome
esearch -db sra -query $ref_samn  | efetch --format runinfo | cut -d ',' -f 1 | grep SRR | xargs -n 1 -P 12 fastq-dump --split-files --gzip
cat $REF_ACC > tmp_acc.txt
bit-dl-ncbi-assemblies -w tmp_acc.txt -f fasta
mv GCF*.fa > ${REF_STRAIN}_ref_genome.fa
rm tmp_acc.txt

## Download db's for quast
quast-download-silva
quast-download-busco

## exit
cd ..
conda deactivate

# setup anvio
#curl https://merenlab.org/files/anvio-conda-environments/anvio-environment-7-LINUX.yaml \
#     --output anvio-environment-7.yaml
if conda env list | grep -q 'anvio-7'; then
  echo "anvio-7 environment is installed, let's look at some pangenomes!"
else
  echo "Looks like Anvi'o isn't setup yet. Hold on, this is going to take a while..."
  conda env create -f CONDA-ENVS/env3.yaml -n anvio-7
  conda activate anvio-7
  conda install -y -c bioconda hmmer=3.2.1 muscle

  ## on the hpc and my machine, some python utils are not installing well
  ## using the fix here: https://www.gitmemory.com/issue/merenlab/anvio/1643/761546004
  conda install -c conda-forge python-levenshtein  -y
  conda install -c bioconda pysam  -y
  conda install -c anaconda psutil==5.4.3 -y


  curl -L https://github.com/merenlab/anvio/releases/download/v7/anvio-7.tar.gz --output anvio-7.tar.gz
  pip install anvio-7.tar.gz
  pip install mistune==0.8.4

  ## setup annotation dbs
  cd ../DBs
  anvi-setup-ncbi-cogs -T 4 --cog-data-dir './COG-DATA-DIR'
  anvi-setup-kegg-kofams --kegg-data-dir './KEGG-DATA-DIR'
  anvi-setup-pfams --pfam-data-dir './PFAM_DATA_DIR'
  anvi-setup-scg-taxonomy -T 4 --scgs-taxonomy-data-dir .

  cd ..
  conda deactivate
fi

# Phlann dbs
## checking for neccessary conda install
if conda env list | grep -q 'phlann'; then
  echo "Phlann is setup, now checking if the databases have been downloaded."
else
  echo "Phlann is not setup, installing from conda."
  conda create -y -f CONDA-ENVS/env6.yml -n phlann python=3.7
fi

## Download the humann databases
conda activate phlann
cd DBs/
DIR="/.humann_db/chocophlan"
if [ ! -d "$DIR" ]; then
  echo "$0: File '${DIR}' not found. Downloading... this is gonna take a sec..."
  mkdir -p humann_db
  humann_databases --download chocophlan full humann_db --update-config yes
  humann_databases --download uniref uniref90_diamond humann_db --update-config yes
  humann_databases --download utility_mapping full humann_db --update-config yes
else
  echo "Humann DBs are installed. Checking Metaphlan now."
fi
cd ..
conda deactivate

## Download metaphlan database
file='/.metaphlan_db/mpa_v30_CHOCOPhlAn_201901.fna.bz2'
if [ ! -d "$FILE" ]; then
  echo "$0: File '${FILE}' not found. Downloading... this is gonna take a sec..."
  mkdir -p metaphlan_db
  metaphlan --install --bowtie2db metaphlan_db
else
  echo "Metaphlan DBs installed. Moving on..."
fi
