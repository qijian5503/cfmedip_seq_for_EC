## lasso               ##########################
rm(list=ls())
gc()

load("lasso_gene_methy.RData")
plot(cv_fit)


library(pROC)
library(rms)
load("medip_esca_peak_samp146_norm_voom_10_01_rf_modlist.RData")
load("medip_esca_peak_samp146_subsamp.RData")
all_test_pre=data.frame()
all_valid_pre=data.frame()
for (i in 1:100) {
    limma_pre=limma_modlist[[i]]$limma_pre
    test_samp=subsample[[i]]$test_samp
    all_test_pre=rbind(all_test_pre,limma_pre[limma_pre$sample %in% test_samp,])
    all_valid_pre=rbind(all_valid_pre,limma_pre[!limma_pre$sample %in% test_samp,])
}

g=pROC::ggroc(list(`Test`=roc(as.factor(all_test_pre$group), as.numeric(all_test_pre$pre)),
                   `Validation`=roc(as.factor(all_valid_pre$group), as.numeric(all_valid_pre$pre))))
auc1=roc(as.factor(all_test_pre$group), as.numeric(all_test_pre$pre))
auc2=roc(as.factor(all_valid_pre$group), as.numeric(all_valid_pre$pre))
p2=g+ggsci::scale_color_lancet()+
    geom_segment(aes(x=0,xend=1,y=1,yend=0),color="darkgrey",linetype=4)+
    theme_classic() + theme(legend.title = element_blank(),legend.position = "")+
    annotate("text", y = 0.6, x=0.5, label = paste0("Test set AUC = ",round(auc1$auc[[1]],4),"\n95% CI ( ",round(ci(auc1)[1],4)," - ",round(ci(auc1)[3],4)," )"),colour ="#00468BFF")+
    annotate("text", y = 0.45, x=0.5, label = paste0("Validation set AUC = ",round(auc2$auc[[1]],4),"\n95% CI ( ",round(ci(auc2)[1],4)," - ",round(ci(auc2)[3],4)," )"),colour ="#ED0000FF")
p2




load("medip_esca_peak_samp146_norm_voom_10_01_rf_modlist.RData")
all_pre=data.frame()
for (i in 1:100) {
    all_pre=rbind(all_pre,limma_modlist[[i]]$limma_pre)
}
all_pre=do.call(rbind,lapply(unique(all_pre$sample),function(x){
    data.frame(sample=x,risk=median(all_pre[all_pre$sample %in% x,"pre"]))
}))

load("medip_esca_samp146_phen.RData")
table(phen_esca$Stage)
phen_esca$Stage=as.character(phen_esca$Stage)
phen_esca$Risk_Score=all_pre[match(phen_esca$Sample,all_pre$sample),"risk"]
phen_esca$Stage=ifelse(phen_esca$Stage %in% c("StageI","StageII"),"Early",
                       ifelse(phen_esca$Stage %in% c("StageIII"),"Later",phen_esca$Stage))
median(phen_esca$Age,na.rm = T)
phen_esca$Age=ifelse(phen_esca$Age<52,"<52",ifelse(phen_esca$Age>=52,">=52",NA))
median(phen_esca$Concentration,na.rm = T)
phen_esca$Concentration=ifelse(phen_esca$Concentration<4.21,"<4.21",ifelse(phen_esca$Concentration>=4.21,">=4.21",NA))

colnames(phen_esca)
dat=reshape2::melt(phen_esca[,c(1:5,11)],c("Sample","Risk_Score","Group"))
dat$sub_group=paste0(dat$variable," ",dat$value)
dat=dat[!is.na(dat$value),]

dat3=reshape2::melt(phen_esca[,c(1,5,6,11)],c("Sample","Risk_Score","Group"))
dat3=dat3[!is.na(dat3$value),]
dat3=dat3[!dat3$value %in% "Normal",]
dat3$sub_group=paste0(dat3$variable," ",dat3$value)

p3=ggplot(dat, aes(sub_group, Risk_Score,fill=sub_group)) + 
    geom_boxplot()+geom_jitter(width = 0.1, alpha = 0.5, color = 'black')+
    scale_fill_manual(values = c("#0571B0","#CA0020","#00468BFF","#ED0000FF","#42B540FF","#0099B4FF","#925E9FFF","#FDAF91FF"))+
    stat_compare_means(aes(group=sub_group),label.y = 1.1,label = "p")+
    facet_grid(.~Group)+theme_bw()+ylab("Predicted risk score")+
    theme(axis.title.x = element_blank(),axis.text.x = element_blank(),panel.grid = element_blank(),legend.title = element_blank())
p3

p4=ggplot(dat3, aes(sub_group, Risk_Score,fill=sub_group)) + 
    geom_boxplot()+geom_jitter(width = 0.1, alpha = 0.5, color = 'black')+
    scale_fill_manual(values = c("#0571B0","#CA0020","#00468BFF","#ED0000FF","#42B540FF","#0099B4FF","#925E9FFF","#FDAF91FF"))+
    stat_compare_means(aes(group=sub_group),label.y = 1.1,label = "p")+
    theme_bw()+ylab("Predicted risk score")+
    theme(axis.title.x = element_blank(),axis.text.x = element_blank(),panel.grid = element_blank(),legend.title = element_blank())
p4




load("medip_esca_peak_samp146_norm_voom_10_01_rf_modlist.RData")
all_pre=data.frame()
for (i in 1:100) {
    all_pre=rbind(all_pre,limma_modlist[[i]]$limma_pre)
}
all_pre=do.call(rbind,lapply(unique(all_pre$sample),function(x){
    data.frame(sample=x,risk=median(all_pre[all_pre$sample %in% x,"pre"]))
}))

load("medip_esca_samp146_phen.RData")
table(phen_esca$Stage)
phen_esca$Stage=as.character(phen_esca$Stage)
phen_esca$Risk_Score=all_pre[match(phen_esca$Sample,all_pre$sample),"risk"]
phen_esca$Stage=ifelse(phen_esca$Stage %in% c("StageI","StageII"),"Early",
                       ifelse(phen_esca$Stage %in% c("StageIII"),"Later",phen_esca$Stage))

for (i in c("Early","Later")) {
    dat=phen_esca[phen_esca$Stage %in% c("Normal",i),]
    assign(paste0(i,"_auc"),roc(as.factor(dat$Group), as.numeric(dat$Risk_Score)))
    assign(paste0(i,"_auc1"),roc(as.factor(dat$Group), as.numeric(dat$Risk_Score))$auc[[1]])
}
g1=pROC::ggroc(list(`Early`=Early_auc,`Later`=Later_auc))

roc.test(Later_auc,Early_auc)$p.value
p5=g1+ggsci::scale_color_lancet() + labs(x = "1 - Specificity", y = "Sensitivity") + 
    theme_classic() + theme(legend.title = element_blank())+
    geom_segment(aes(x=0,xend=1,y=1,yend=0),color="darkgrey",linetype=4)+
    annotate("text", y = 0.6, x=0.5, label = paste0("p=",format(roc.test(Later_auc,Early_auc)$p.value,scientific=T,digits=4)),colour ="black")+
    annotate("text", y = 0.5, x=0.5, label = paste0("Later AUC = ",round(Later_auc1,4)," \nEarly AUC = ",round(Early_auc1,4)),colour ="#00468B99")
p5




phen_esca$Gender=factor(phen_esca$Gender,levels = c("F","M"))
g2=pROC::ggroc(list(Age=roc(as.factor(phen_esca$Group), as.numeric(phen_esca$Age)),
                    Gender=roc(as.factor(phen_esca$Group), as.numeric(phen_esca$Gender)),
                    Concentration=roc(as.factor(phen_esca$Group), as.numeric(phen_esca$Concentration)),
                    Risk=roc(as.factor(phen_esca$Group), as.numeric(phen_esca$Risk_Score))))
auc1=roc(as.factor(phen_esca$Group), as.numeric(phen_esca$Age))$auc[[1]]
auc2=roc(as.factor(phen_esca$Group), as.numeric(phen_esca$Gender))$auc[[1]]
auc3=roc(as.factor(phen_esca$Group), as.numeric(phen_esca$Concentration))$auc[[1]]
auc4=roc(as.factor(phen_esca$Group), as.numeric(phen_esca$Risk_Score))$auc[[1]]


p6=g2+ggsci::scale_color_lancet() + labs(x = "1 - Specificity", y = "Sensitivity") + 
    theme_classic() + theme(legend.title = element_blank())+
    geom_segment(aes(x=0,xend=1,y=1,yend=0),color="darkgrey",linetype=4)+
    annotate("text", y = 0.5, x=0.5, label = paste0("Age AUC = ",round(auc1,4),"\nGender AUC = ",round(auc2,4),"\nConcentration AUC = ",round(auc3,4),"\nRisk AUC = ",round(auc4,4)),colour ="#00468B99")
p6
