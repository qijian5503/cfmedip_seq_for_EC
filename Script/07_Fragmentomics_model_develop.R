rm(list=ls())
gc()

load("medip_esca_samp146_phen.RData")
load("medip_esca_peak_samp146_valid_samp.RData")
load("medip_esca_peak_samp146_subsamp.RData")

load("end_motif_all_de_res_samp146.RData")
all_de=all_de[order(all_de$p),]
de_edm=all_de[1:50,1]
load("end_motif_all_freq.RData")
rownames(a)=a$Sample
a$Sample=NULL
exp_edm=a[,de_edm]

load("frag_prop_all_de_res_samp146.RData")
all_de=all_de[order(all_de$p),]
de_frag=as.character(all_de[1:50,1])
load("frag_prop_all_exp_samp146.RData")
exp_frag=as.data.frame(t(all_exp[de_frag,]))

load("bpm_end_motif_all_de_res_samp146.RData")
all_de=all_de[order(all_de$p),]
de_bpm=as.character(all_de[1:50,1])
load("bpm_end_motif_all_freq.RData")
exp_bpm=bpm_motif_freq_all[,c(de_bpm,"Group")]

comb_data=cbind(exp_edm,exp_frag,exp_bpm)

library(glmnet)
de_exp=comb_data[!rownames(comb_data) %in% valid_sample,]
de_exp$Group=factor(de_exp$Group,levels = c("Normal","Tumor"))
cv_fit = cv.glmnet(as.matrix(de_exp[,!colnames(de_exp) %in% "Group"]), de_exp$Group, family = "binomial",type.measure = "class")
coef.min = coef(cv_fit, s = "lambda.min")
lasso_gene=coef.min@Dimnames[[1]][coef.min@i+1][-1]  ## 筛选得到的基因
lasso_gene
# save(cv_fit,file = "lasso_cv_fit_frag.RData")

lasso_gene=gsub("-","_",lasso_gene)
colnames(comb_data)=gsub("-","_",colnames(comb_data))
rownames(comb_data)=phen_esca[match(rownames(comb_data),phen_esca$Sample2),"Sample"]

ctrl <- trainControl(method = "repeatedcv",number = 10,repeats = 5,summaryFunction = twoClassSummary,classProbs = TRUE)
limma_modlist=list()
for(i in 1:100){
    train_samp=subsample[[i]]$train_samp
    train_samp=train_samp[train_samp %in% rownames(comb_data)]
    test_samp=subsample[[i]]$test_samp
    test_samp=test_samp[test_samp %in% rownames(comb_data)]
    train_data=comb_data[rownames(comb_data) %in% train_samp,c(lasso_gene,"Group")]
    test_data=comb_data[rownames(comb_data) %in% test_samp,c(lasso_gene,"Group")]
    
    rf_mod=randomForest(Group~. ,train_data,important=TRUE,proximity=TRUE,ntree=100,mtry=10)
    test_pre=predict(rf_mod,test_data,type="prob")%>%data.frame
    test_pre$group=as.factor(phen_esca[rownames(test_pre),"Group"])
    
    valid_pre=predict(rf_mod,comb_data[valid_sample,],type="prob")%>%data.frame
    valid_pre$group=as.factor(phen_esca[rownames(valid_pre),"Group"])
    
    limma_pre=rbind(data.frame(type="test",sample=rownames(test_pre),group=test_pre$group,pre=test_pre$Tumor),
                    data.frame(type="valid",sample=rownames(valid_pre),group=valid_pre$group,pre=valid_pre$Tumor))
    limma_modlist[[i]]= list(limma_pre=limma_pre,rf_mod=rf_mod)
}
# save(limma_modlist,file="medip_esca_peak_samp146_frag_comb_data_rf_modlist.RData") # 0.9987 0.9938

