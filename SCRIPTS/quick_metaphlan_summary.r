#!/usr/bin/env Rscript
# last updated: 081021
args = commandArgs(trailingOnly=TRUE)

message("arguments:PATH/TO/merged_abundance_table_species.txt \n")

# Import data
dat <- read.table(args[1],header = T,row.names = 1)
tmp <- dat

# make summary table
tmp_mean <- apply(tmp,1,FUN = mean)
tmp_sd <- apply(tmp,1,FUN = sd)
tmp_max <- apply(tmp,1,FUN=max)
tmp_min <- apply(tmp,1,FUN=min)

tmp_out <- data.frame(mean=tmp_mean,sd=tmp_sd,min=tmp_min,max=tmp_max)

# make output table
write.csv(tmp_out,'DATA/METAPHLAN/merged_abundance_summary.csv')


