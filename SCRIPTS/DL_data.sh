#!/usr/bin/env bash
mkdir -p DBs
mkdir -p LOGS
mkdir -p DATA

# Setup log files
exec 3>&1 4>&2
trap 'exec 2>&4 1>&3' 0 1 2 3
exec 1>LOGS/DL_data_log.out 2>&1

# create and environment that will let us download the raw seq data we are interested in
if conda env list | grep -q 'env0'; then
  echo "dl_data environment is installed, checking to see if we have our raw data"
else
  echo "Setting up the dl_data environment"
  conda create -y -f CONDA-ENVS/env0.yml
fi

# download sequence reads
if [ ! -d DATA/RAW ]; then
  conda activate dl_data
  echo "RAW directory does not exist"
  mkdir RAW
  cd RAW
  SAMN=$(grep 'short_read' ../sra_accessions.txt | cut -f2)
  STRAIN=$(grep 'short_read' ../sra_accessions.txt | cut -f1)

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
  conda deactivate
else
  echo "RAW directory exists"
fi

# Download ONT reads
if [ ! -d DATA/RAW/JTC_unzip_fast5 ]; then
  # Download fast5 files (this takes a while, definitely open a screen)
  fileid="19T7YU0aWuZyN6WggQdslsP1VmzYL8BAW"
  filename="DATA/ONT/raw_ONT_fast5.tar.gz"
  curl -c ./cookie -s -L "https://drive.google.com/uc?export=download&id=${fileid}" > /dev/null
  curl -Lb ./cookie "https://drive.google.com/uc?export=download&confirm=`awk '/download/ {print $NF}' ./cookie`&id=${fileid}" -o ${filename}
  tar -xzvf raw_ONT_fast5.tar.gz
  mkdir DATA/RAW/JTC_unzip_fast5
  unzip "JTC_Kling_fast5/*zip" DATA/RAW/JTC_unzip_fast5/
  rm raw_ONT_fast5.tar.gz
else
  echo "Looks like you've already downloaded the ONT fast5 file."
fi


#Download reference genome
if [ ! -d DATA/REF ]; then
  conda activate dl_data
  mkdir DATA/REF

  # get accession and strain name
  REF_ACC=$(grep 'reference' ../refseq_accessions.txt | cut -f3)
  REF_STRAIN=$(grep 'reference' ../refseq_accessions.txt | cut -f1)
  cat $REF_ACC > tmp_acc.txt

  # Download genome
  bit-dl-ncbi-assemblies -w tmp_acc.txt -f fasta
  mv G*.fa > DATA/REF/${REF_STRAIN}_ref_genome.fa
  rm tmp_acc.txt
else
  echo "Reference genome is downloaded and can be found in: DATA/REF_GENOME/ "
fi

# setup anvio
if conda env list | grep -q 'env3'; then
  echo "anvio-7 environment is installed, let's look at some pangenomes!"
else
  echo "Looks like Anvi'o isn't setup yet. Hold on, this is going to take a while..."
  conda env create -f CONDA-ENVS/env3.yaml -n anvio-7
  conda activate anvio-7
  # need to install the actual anvio package from their website
  cd ../DBs
  curl -L https://github.com/merenlab/anvio/releases/download/v7/anvio-7.tar.gz --output anvio-7.tar.gz
  pip install anvio-7.tar.gz

  # setup annotation dbs
  anvi-setup-ncbi-cogs -T 4 --cog-data-dir './COG-DATA-DIR'
  anvi-setup-kegg-kofams --kegg-data-dir './KEGG-DATA-DIR'
  anvi-setup-pfams --pfam-data-dir './PFAM_DATA_DIR'
  anvi-setup-scg-taxonomy -T 4 --scgs-taxonomy-data-dir .
  cd ..
  conda deactivate
fi

# Phlann dbs
if conda env list | grep -q 'env4'; then
  echo "Phlann is setup."
else
  echo "Phlann is not setup, installing from conda."
  conda create -y -f CONDA-ENVS/env4.yml python=3.7
  # Setup phlan data bases
  conda activate env4
  cd DBs
  # humann
  echo "Downloading humann databases... this is gonna take a sec..."
  mkdir humann_db
  humann_databases --download chocophlan full humann_db --update-config yes
  humann_databases --download uniref uniref90_diamond humann_db --update-config yes
  humann_databases --download utility_mapping full humann_db --update-config yes

  # metaphlann
  mkdir metaphlan_db
  metaphlan --install --index mpa_v30_CHOCOPhlAn_201901 --bowtie2db metaphlan_db
  # need to install metaphlan dbs outside of SNAKEMAKE
  ## /scratch/joshuakl/syn_final/.snakemake/conda/0eea8dbb/bin/metaphlan --install --index mpa_v30_CHOCOPhlAn_201901 --bowtie2db metaphlan_db
  cd ..
  conda deactivate
fi

# setup envs for epigenomics
if conda env list | grep -q 'env5'; then
  echo "All programs for gene calling and annotation is setup."
else
  echo "Installing packages for gene calling and annotation"
  conda create -y -f CONDA-ENVS/env5.yml
  # get kofamscan DBs
  wget ftp://ftp.genome.jp/pub/db/kofam/ko_list.gz
  gunzip -c ko_list.gz > DBs/
  rm ko_list.gz
  wget ftp://ftp.genome.jp/pub/db/kofam/profiles.tar.gz
  tar xzf profiles.tar.gz -C DBs/
  rm profiles.tar.gz
fi

# nanopolish
if conda env list | grep -q 'env6'; then
  echo "All programs for calling methylation sites is setup."
else
  echo "Installing packages for calling methylation sites"
  conda create -y -f CONDA-ENVS/env6.yml
fi

# r
if conda env list | grep -q 'env7'; then
  echo "R and all packages that you will need are setup."
else
  echo "Installing R and all of the packages you will need"
  conda create -y -f CONDA-ENVS/env7.yml
fi
