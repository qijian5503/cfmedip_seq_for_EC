rm(list=ls())
gc()
options(stringsAsFactors = F)

load("medip_esca_peak_samp146_subsamp.RData")
load("medip_esca_peak_samp146_traintest_norm_voom_10_01_noxy.RData")
load("medip_esca_samp146_phen.RData")

phen_esca=na.omit(phen_esca[,c(2,3,5)])
phen_esca$Gender=as.numeric(factor(phen_esca$Gender,levels = c("F","M")))-1
phen_esca$Group=factor(phen_esca$Group,levels = c("Normal","Tumor"))

limma_de_all=list()
for(i in 1:100){
    train_samp=subsample[[i]]$train_samp
    train_samp=train_samp[train_samp %in% rownames(phen_esca)]
    coldata=phen_esca[train_samp,]
    psm <- matchit(Group~Age+Gender,data=coldata,method="nearest",distance = "logit",ratio = 2)
    coldata <- match.data(psm)
    design=model.matrix(~0+coldata$subclass+coldata$Group)
    colnames(design)=gsub("coldata[$]","",colnames(design))
    res=lmFit(norm[,rownames(coldata)],design = design,) %>% eBayes(.,trend = T) %>% topTable(., number =Inf,adjust.method = "BH",coef="GroupTumor") %>% .[order(.$adj.P.Val),]  %>% na.omit(.)
    limma_de_all[[i]]=res[abs(res$logFC)>1 & res$adj.P.Val<0.01,]
    message(nrow(limma_de_all[[i]]))
}
save(limma_de_all,file="medip_esca_peak_samp146_norm_voom_10_01_limma_PSM_de_noxy_all.RData")


