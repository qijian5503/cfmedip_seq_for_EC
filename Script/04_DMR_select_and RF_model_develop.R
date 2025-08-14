rm(list=ls())
gc()
options(stringsAsFactors = F)
load("medip_esca_peak_samp146_subsamp.RData")
load("medip_esca_peak_samp146_traintest_norm_voom_10_01_noxy.RData")
load("medip_esca_samp146_phen.RData")
load("medip_esca_peak_samp146_norm_voom_10_01_valid_noxy.RData")
load("medip_esca_peak_samp146_norm_voom_10_01_limma_PSM_de_noxy_all.RData")

all_de2=data.frame()
for (r in 1:100) {
    res=limma_de_all[[r]]
    res=res[abs(res$logFC)> 3 & res$adj.P.Val<0.0001,]
    message(nrow(res))
    res$regule=ifelse(res$logFC>0,"Hypermethylation","Hypomethylation")
    res$peak=rownames(res)
    if(nrow(res)>0){all_de2=rbind(all_de2,data.frame(round=r,res))}
}
rownames(all_de2)=NULL

imp_de_gene=as.data.frame(table(all_de2$peak,all_de2$regule))
imp_de_gene=reshape2::dcast(imp_de_gene,Var1~Var2,value.var = "Freq")
de_gene=as.character(imp_de_gene[imp_de_gene$Hypermethylation>90 | imp_de_gene$Hypomethylation>90,"Var1"])

phen_esca$Group=factor(phen_esca$Group,levels = c("Normal","Tumor"))
de_exp=as.data.frame(t(norm[de_gene,]))
de_exp$group=phen_esca[rownames(de_exp),"Group"]
valid_data=as.data.frame(t(norm_v[de_gene,]))
valid_data$group=as.factor(phen_esca[rownames(valid_data),"Group"])

library(glmnet)
cv_fit = cv.glmnet(as.matrix(de_exp[,de_gene]), de_exp$Group, family = "binomial",type.measure = "class")
coef.min = coef(cv_fit, s = "lambda.min")
lasso_gene=coef.min@Dimnames[[1]][coef.min@i+1][-1]  ## 筛选得到的基因
# save(cv_fit,file = "lasso_gene_methy.RData")


ctrl <- trainControl(method = "repeatedcv",number = 10,repeats = 5,summaryFunction = twoClassSummary,classProbs = TRUE)
limma_modlist=list()
for(i in 1:100){
    train_samp=subsample[[i]]$train_samp
    train_samp=train_samp[train_samp %in% rownames(de_exp)]
    test_samp=subsample[[i]]$test_samp
    test_samp=test_samp[test_samp %in% rownames(de_exp)]
    train_data=de_exp[rownames(de_exp) %in% train_samp,c(lasso_gene,"group")]
    test_data=de_exp[rownames(de_exp) %in% test_samp,c(lasso_gene,"group")]
    
    rf_mod=randomForest(group~. ,train_data,important=TRUE,proximity=TRUE,ntree=100,mtry=10)
    test_pre=predict(rf_mod,test_data,type="prob")%>%data.frame
    test_pre$group=as.factor(phen_esca[rownames(test_pre),"Group"])
    
    valid_pre=predict(rf_mod,valid_data,type="prob")%>%data.frame
    valid_pre$group=as.factor(phen_esca[rownames(valid_pre),"Group"])
    
    limma_pre=rbind(data.frame(type="test",sample=rownames(test_pre),group=test_pre$group,pre=test_pre$Tumor),
                    data.frame(type="valid",sample=rownames(valid_pre),group=valid_pre$group,pre=valid_pre$Tumor))
    limma_modlist[[i]]= list(limma_pre=limma_pre,rf_mod=rf_mod)
}
save(limma_modlist,file="medip_esca_peak_samp146_norm_voom_10_01_rf_modlist.RData")


