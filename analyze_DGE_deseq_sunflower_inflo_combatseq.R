# Name: analyze DGE deseq sunflower inflo ruvseq
# Author: EY (based off of code written by E. Dittmar)
# Date: 08/29/2023
# Version:4.2.1
# Description: Will analyze the output from the DESeq DGE wit ruvseq for sunflower inflo stages
# need Functions.R written by ED



setwd('/home/ely67071/sunflower_inflo_dev_analysis/')

library(dplyr)
library(ggplot2)
library(UpSetR)
library(Glimma)
source("Functions.R")


# now analyze 
# read in the data
DEData_pairwise_cs<-ImportCSVs('deseq_results/combatseq/pairwise/',0.05)
# filter out significant results
mydataSig_pairwise_cs<-lapply(DEData_pairwise_cs,SigDEdf,PvaluesCol=7,CritP=0.05)

# see which genes overlap
SigOverlap_pairwise_cs<-GeneSets(mydataSig_pairwise_cs$result_10D_v_20D[1], mydataSig_pairwise_cs$result_20D_v_30D[1],mydataSig_pairwise_cs$result_30D_v_35D[1])
names(SigOverlap_pairwise_cs)
lapply(SigOverlap_pairwise_cs,function(x) {length(x$Gene)})

SigOverlapGraph_pairwise_cs<-lapply(mydataSig_pairwise_cs, function(x) {x$Gene})

png("plots/sequential_pairwise_upset_cs_larger.png", res=215, width = 1800, height=1000)
upset(fromList(SigOverlapGraph_pairwise_cs),order.by="freq",nsets=13,nintersects=20, text.scale = 1.5)
dev.off()


# now analyze from LRT
# read in the data
DEData_pairwise_lrt<-ImportCSVs('deseq_results/combatseq/lrt/',0.05)
# filter out significant results
mydataSig_pairwise_lrt<-lapply(DEData_pairwise_lrt,SigDEdf,PvaluesCol=6,CritP=0.05)

# see which genes overlap
SigOverlap_pairwise_lrt<-GeneSets(mydataSig_pairwise_lrt$result_10D_v_20D[1], mydataSig_pairwise_lrt$result_20D_v_30D[1],mydataSig_pairwise_lrt$result_30D_v_35D[1])
names(SigOverlap_pairwise_lrt)
lapply(SigOverlap_pairwise_lrt,function(x) {length(x$Gene)})

SigOverlapGraph_pairwise_lrt<-lapply(mydataSig_pairwise_lrt, function(x) {x$Gene})

png("plots/sequential_pairwise_upset_cs_lrt.png", res=215, width = 1200, height=1000)
upset(fromList(SigOverlapGraph_pairwise_lrt),order.by="freq",nsets=13,nintersects=20)
dev.off()