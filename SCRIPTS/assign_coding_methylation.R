#!/usr/bin/env Rscript
args = commandArgs(trailingOnly=TRUE)

message("arguments: output_prefix plus_strand minus_strand .gff \n")

library(stringr)
if("stringr" %in% rownames(installed.packages()) == FALSE){print('stringr not installed. Please install in conda env8 and try again.')}

if (length(args)==0) {
  stop("At least one argument must be supplied (input file).n", call.=FALSE)
} #else if (length(args)==1) {
  # default output file
  #args[2] = "out.txt"
#}

# inputs
output_pre <- args[1]
nano_plus <- read.table(args[2],header = T)
nano_minus <- read.table(args[3],header = T)
gff_tmp <- read.table(args[4])[,c(1,4,5,7,9)]

# formatting gff file
colnames(gff_tmp) <- c('contig','cds_start','cds_end','strand','ID')
ids <- str_extract(str_extract(gff_tmp$ID,"[^;]+"),"_[:digit:]+")
gff_tmp$gene_ID <- paste(gff_tmp$contig,ids,sep = '')

assign_CDS <- function(nanopolish_out,gff_file){
  unlist(lapply(sapply(1:nrow(nanopolish_out),function(i){
    start <- nanopolish_out[i,]$start
    stop <- nanopolish_out[i,]$end
    gff_file[gff_file$cds_start < start & gff_file$cds_end > stop,]$gene_ID[1]
  }),function(y) if(identical(y, character(0))) NA_character_ else y))
}

# 84 over for plus and 84 over for minus
# for right now, deciding that's an acceptable margin of error so going to just take the first hit.
nano_plus$gene_ID <- assign_CDS(nano_plus,gff_tmp)
nano_minus$gene_ID <- assign_CDS(nano_minus,gff_tmp)

# add strand information and combine
nano_plus$strand <- c(rep('+',nrow(nano_plus)))
nano_minus$strand <- c(rep('-',nrow(nano_minus)))
nano_all <- rbind(nano_plus,nano_minus)
nano_all <- nano_all[order(nano_all$methylated_frequency,decreasing = TRUE),]

# setting output
write.table(nano_all,paste(output_pre,'nano_cds.tsv',sep='_'))
