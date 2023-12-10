# Name: run_DGE_sunflower_inflo_combatseq.R
# Author: EY
# Date: 08/29/2023
# Version:4.2.1
# Description: Will run the differential gene expression on the infloresence stages accounting for variation with combatseq

library(dplyr)
library(DESeq2)
library(Glimma)
library(sva)
library(edgeR)
source("Functions.R")



# read in and process data

setwd('/home/ely67071/sunflower_inflo_dev_analysis/')

# read in the data matrix
summed_counts<-readRDS("/scratch/ely67071/sunflower_inflo_dev_data/collapsed_replicates_deseq_dataset.Rdata")

samples=c("10D_REP1_ATTACTCG", "20D_REP2_TCCGGAGA" ,"30D_REP2_CGCTCATT", "35D_REP1_GAGATTCC", 
          "HA_10D_2_ACCTTGGC", "HA_10D_3_ATATCTCG", "HA_20D_2_GCGCTCTA", 
          "HA_20D_3_AACAGGTT", "HA_30D_2_GGTGAACC", "HA_30D_3_CAACAATG", "HA_35D_2_TGGTGGCA", "HA_35D_3_AGGCAGAG")

dev_stage<-sub(".*([0-9]{2,2}D).*", "\\1",samples)
dev_stage <- as.factor(dev_stage)

metadata<-data.frame(samples, dev_stage)

# create the factors of interest
metadata$dev_stage<-factor(metadata$dev_stage)

# create the model
summed_counts<-DESeqDataSetFromMatrix(counts(summed_counts),colData = metadata, design=~0+dev_stage)

# pre-filter for reads where at least 3 samples have a count of 1 or higher
keep<-rowSums(counts(summed_counts)>=1)>=3
length(which(keep==1))
summed_counts_filt<-summed_counts[keep,]

# get a dataframe of the counts
count_matrix <- as.matrix(counts(summed_counts_filt))

batch <- c(1,1,1,1,2,3,2,3,2,3,2,3)
group <- c(1,2,3,4,1,1,2,2,3,3,4,4)

adjusted_counts <- ComBat_seq(count_matrix, batch=batch, group=group)

write.csv(as.data.frame((adjusted_counts)), file='deseq_results/adjusted_counts_combatseq.csv')

# plot MDS of adjusted counts
glimmaMDS(adjusted_counts, group=metadata)

# plot PCA of adjusted counts
png("plots/pca_combseq.png", res=215, width = 1200, height=1000)
par(mar=c(5.1, 4.1, 4.1, 8.1), xpd=TRUE)
plotPCA(adjusted_counts,labels=TRUE, col=dev_stage)
legend("topright",inset=c(-0.4,0),legend=unique(dev_stage), fill=dev_stage)
dev.off()

# create the model
adjusted_counts_deseq<-DESeqDataSetFromMatrix((adjusted_counts),colData = metadata, design=~0+dev_stage)


# run the DGE analysis
DESeq_dataset_results_combatseq<-DESeq(adjusted_counts_deseq,parallel=TRUE)
saveRDS(DESeq_dataset_results_combatseq, file='deseq_results/combatseq/pairwise/deseq_dataset_results_combatseq.RData')


result_10D_v_20D_combatseq<-results(DESeq_dataset_results_combatseq,contrast=c("dev_stage","10D","20D"),alpha=0.05,parallel=TRUE)
result_20D_v_30D_combatseq<-results(DESeq_dataset_results_combatseq,contrast=c("dev_stage","20D","30D"),alpha=0.05,parallel=TRUE)
result_30D_v_35D_combatseq<-results(DESeq_dataset_results_combatseq,contrast=c("dev_stage","30D","35D"),alpha=0.05,parallel=TRUE)


write.csv(as.data.frame(result_10D_v_20D_combatseq), file='deseq_results/combatseq/pairwise/result_10D_v_20D.csv')
write.csv(as.data.frame(result_20D_v_30D_combatseq), file='deseq_results/combatseq/pairwise/result_20D_v_30D.csv')
write.csv(as.data.frame(result_30D_v_35D_combatseq), file='deseq_results/combatseq/pairwise/result_30D_v_35D.csv')




# try in edgeR


design<-model.matrix(~0+dev_stage, data=metadata)
y<- DGEList(counts=(adjusted_counts), group=dev_stage)
y<-calcNormFactors(y,method="TMM")

# estimate dispersion
y <-estimateGLMCommonDisp(y, design)
y <-estimateGLMTagwiseDisp(y,design)
plotBCV(y)

# fit the model
fit<-glmQLFit(y,design)

# use the output from the line below to design the contrasts 
colnames(fit)


res_10v20<-glmQLFTest(fit, contrast=c(-1,1,0,0))
summary(decideTests(res_10v20))
plotMD(res_10v20)
top_10v20<-topTags(res_10v20, n=Inf)
write.csv(as.data.frame(top_10v20), file='deseq_results/combatseq/lrt/result_10D_v_20D.csv')



res_20v30<-glmQLFTest(fit, contrast=c(0,-1,1,0))
summary(decideTests(res_20v30))
plotMD(res_20v30)
top_20v30<-topTags(res_20v30, n=Inf)
write.csv(as.data.frame(top_20v30), file='deseq_results/combatseq/lrt/result_20D_v_30D.csv')




res_30v35<-glmQLFTest(fit, contrast=c(0,0,-1,1))
summary(decideTests(res_30v35))
plotMD(res_30v35)
top_30v35<-topTags(res_30v35, n=Inf)
write.csv(as.data.frame(top_30v35), file='deseq_results/combatseq/lrt/result_30D_v_35D.csv')



