#!/usr/bin/env bash

# salloc --cpus-per-task=16 --time=3:00:00 --mem=64GB

# Setup log files
exec 3>&1 4>&2
trap 'exec 2>&4 1>&3' 0 1 2 3
exec 1>LOGS/epigenomics_log.out 2>&1

# Define reference strain
REF_STRAIN=$(grep 'reference' ../refseq_accessions.txt | cut -f1)

# Make directory for the output
mkdir -p DATA/EPIGENOMICS

# Getting coding sequences and from reference genome
conda activate env5
# Gene calls
prodigal -i DATA/REF/${REF_STRAIN}.fa -o DATA/EPIGENOMICS/${REF_STRAIN}.gff -a DATA/EPIGENOMICS/${REF_STRAIN}.faa -f gff
# Annotation with kegg
exec_annotation -o DATA/EPIGENOMICS/${REF_STRAIN}_annotations.tsv --profile=profiles/ --ko-list=ko_list --cpu=16 -f detail-tsv --report-unannotated REF/${REF_STRAIN}.faa

bit-filter-KOFamScan-results -i DATA/REF/${REF_STRAIN}_annotations.tsv -o DATA/REF/filt_${REF_STRAIN}_annotations.tsv
conda deactivate


# Switch conda environments for methylation calls
conda deactivate
conda activate env6

# Index
nanopolish index -d DATA/RAW/JTC_unzip_fast5/ DATA/RAW/LA127_ONT.fq.gz
nanopolish index -d DATA/RAW/JTC_unzip_fast5/ DATA/RAW/LA31_ONT.fq.gz

# Map reads to the index
minimap2 -a -x map-ont DATA/REF/${REF_STRAIN}.fa DATA/RAW/LA127_ONT.fq.gz | samtools sort -T tmp -o DATA/EPIGENOMICS/LA127_ONT_forNP.sorted.bam
minimap2 -a -x map-ont DATA/REF/${REF_STRAIN}.fa DATA/RAW/LA31_ONT.fq.gz | samtools sort -T tmp -o DATA/EPIGENOMICS/LA31_ONT_forNP.sorted.bam

# Make bam
samtools index DATA/EPIGENOMICS/LA127_ONT_forNP.sorted.bam
samtools index DATA/EPIGENOMICS/LA31_ONT_forNP.sorted.bam

# make methylation calls
nanopolish call-methylation -t 8 -r DATA/RAW/LA127_ONT.fq.gz -b DATA/EPIGENOMICS/LA127_ONT_forNP.sorted.bam -g DATA/REF/${REF_STRAIN}.fa > DATA/EPIGENOMICS/LA127_methylation_calls.tsv
nanopolish call-methylation -t 8 -r DATA/RAW/LA31_ONT.fq.gz -b DATA/EPIGENOMICS/LA31_ONT_forNP.sorted.bam -g DATA/REF/${REF_STRAIN}.fa > DATA/EPIGENOMICS/LA31_methylation_calls.tsv
conda deactivate

# subset for stranded
# nanopolish output here contains strand specific information; however, if conglomerates them when calculating frequency. Splitting here.
for tsv in DATA/EPIGENOMICS/*_methylation_calls.tsv
do
  conda activate env6
  # split into strands
  BASE=$(basename $tsv _methylation_calls.tsv)
  head $tsv -n1 > tmp_headers1
  cp tmp_headers1 tmp_headers2
  cat $tsv | awk '$2 ~ /\+/' >> tmp_headers1
  cat $tsv | awk '$2 ~ /\-/' >> tmp_headers2
  mv tmp_headers1 ${BASE}_plus_methylation.tsv
  mv tmp_headers2 ${BASE}_minus_methylation.tsv

  # Call meth frequency
  python calculate_methylation_frequency.py ${BASE}_plus_methylation.tsv > ${BASE}_plus_methylation_frequency.tsv
  python calculate_methylation_frequency.py ${BASE}_minus_methylation.tsv > ${BASE}_minus_methylation_frequency.tsv
  conda deactivate


  # format for parsing
  conda activate env7
  Rscript Scripts/assign_coding_methylation.R ${BASE} ${BASE}_plus_methylation_frequency.tsv ${BASE}_minus_methylation_frequency.tsv DATA/EPIGENOMICS/${REF_STRAIN}.gff
  # output {prefix}_nano_cds.tsv
  conda deactivate

  # subsetting the output to get enriched ko's
  ## Originally, I tried .09 as my cutoff but did not see any major differences, now being less stringent with a 70% cutoff
  awk '$8 > 0.7' ${BASE}_nano_cds.tsv > tmp
  tr -d '"' < tmp > ${BASE}_top_methylation_sites.tsv
  rm tmp
done
