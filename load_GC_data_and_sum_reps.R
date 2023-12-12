# Name: load GC data and sum reps
# Author: EY
# Date: 12/9/2022
# Version:4.1.2
# Description: will load in the parsed gene count data
# set cwd
library(DESeq2)

# load the DESeq data set
# we read in the data using the DESeqDataSetFromHTSeqCount bc our dataset is in the same format as HTSeq data
dds_set<-DESeqDataSetFromHTSeqCount(sampleTable=data_table,directory='/scratch/ely67071/sunflower_data/gene_count_data/',design=~0)

# save the data set 
saveRDS(dds_set,file="/scratch/ely67071/sunflower_data/dev_deseq_dataset.Rdata")
