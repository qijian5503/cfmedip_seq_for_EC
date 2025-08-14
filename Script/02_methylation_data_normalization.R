rm(list=ls())
gc()
options(stringsAsFactors = F)

library(limma)
library(edgeR)
library(data.table)
library(tidyverse)

load("medip_esca_peak_samp146_exp.RData")
load("medip_esca_samp146_phen.RData")
load("medip_esca_peak_samp146_valid_samp.RData")
peak=peak[!grepl("chrX|chrY",rownames(peak)),]
peak=peak[!grepl("chrUn|hap|random",rownames(peak)),]

table(phen_esca$Group)
train_test=phen_esca[!phen_esca$Sample %in% valid_sample,]
train_test_samp=train_test$Sample

train_test=peak[,train_test_samp]
keep <- rowSums(train_test > 10) >= ncol(train_test)*0.1

peak_nam=rownames(train_test)[keep]
peak_nam=str_split(peak_nam,"_",simplify = T)
write.table(peak_nam,"peak_keep_name.bed",quote = F,row.names = F,col.names = F,sep = "\t")

train_test=train_test[keep,]
f75 <- rep_len(1,ncol(train_test))
for (j in seq_len(ncol(train_test))) f75[j] <- quantile(train_test[,j], probs=0.75)
lib.size=colSums(train_test)
f75=f75 / lib.size
refColumn <- colnames(train_test)[which.min(abs(f75-mean(f75)))]
valid=peak[keep,c(valid_sample,refColumn)]


d=DGEList(counts = train_test) %>% calcNormFactors(. ,method = "TMM" ,refColumn = refColumn) %>% voom()
norm=as.data.frame(d$E)
save(norm,file = "medip_esca_peak_samp146_traintest_norm_voom_10_01_noxy.RData")

d=DGEList(counts = valid) %>% calcNormFactors(. ,method = "TMM" ,refColumn = refColumn) %>% voom()
norm_v=as.data.frame(d$E)
norm_v[,refColumn]=NULL
save(norm_v,file = "medip_esca_peak_samp146_norm_voom_10_01_valid_noxy.RData")

