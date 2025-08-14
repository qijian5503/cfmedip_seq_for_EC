## edm pheatmap        ###################################################################################
rm(list=ls())
gc()

library(pheatmap)
library(gplots)
library(colorRamps)
library(RColorBrewer)

load("end_motif_all_de_res_samp146.RData")
all_de=all_de[order(all_de$p),]
de_gene_edm=all_de[1:50,1]
load("end_motif_all_freq.RData")
rownames(a)=a$Sample
a$Sample=NULL
exp_edm=a[,de_gene_edm]

load("bpm_end_motif_all_de_res_samp146.RData")
all_de=all_de[order(all_de$p),]
de_gene_bpm=as.character(all_de[1:50,1])
load("bpm_end_motif_all_freq.RData")
exp_bpm=bpm_motif_freq_all[,c(de_gene_bpm,"Group")]


## bpm
nmf_exp=exp_bpm
nmf_exp$Group=NULL
nmf_exp=as.data.frame(t(nmf_exp))

load("medip_esca_samp146_phen.RData")
phen=phen_esca[phen_esca$Sample2 %in% colnames(nmf_exp),]
rownames(phen)=phen$Sample2
colnames(phen)
phen=phen[,c(5,6)]
phen=phen[order(phen$Group),]
phen$Stage=NULL
ann_colors = list(
    Group=c(Tumor="#CA0020", Normal="#0571B0"),
    BPM=c(`A`="#00468BFF",`T`="#FDAF91FF",`C`="#0099B4FF",`G`="#925E9FFF")
)

split_row=rownames(nmf_exp)
split_row=gsub("-6-motif_r","",split_row)
split_row=data.frame(edm=split_row,BPM=substr(split_row,1,1))
rownames(split_row)=split_row$edm
split_row$edm=NULL

p1=ComplexHeatmap::pheatmap(as.matrix(nmf_exp[,rownames(phen)]), scale = "row",cluster_cols = F,
                            show_rownames = F,show_colnames = F,annotation_colors = ann_colors,
                            col = colorRampPalette(c("navy","white","firebrick3"))(10),
                            annotation_col = phen,annotation_row = split_row,
                            row_split = split_row$BPM,annotation_names_row = T,
                            annotation_names_col = T ,column_title = NULL,row_title = NULL)
p1


## edm 
nmf_exp=exp_edm
nmf_exp$Group=NULL
nmf_exp=as.data.frame(t(nmf_exp))

ann_colors2 = list(
    Group=c(Tumor="#CA0020", Normal="#0571B0"),
    EDM=c(`A`="#00468BFF",`T`="#FDAF91FF",`C`="#0099B4FF",`G`="#925E9FFF")
)

split_row2=rownames(nmf_exp)
split_row2=gsub("-motif_r","",split_row2)
split_row2=data.frame(edm=split_row2,EDM=substr(split_row2,1,1))
rownames(split_row2)=split_row2$edm
split_row2$edm=NULL

p2=ComplexHeatmap::pheatmap(as.matrix(nmf_exp[,rownames(phen)]), scale = "row",cluster_cols = F,
                            show_rownames = F,show_colnames = F,annotation_colors = ann_colors2,
                            col = colorRampPalette(c("navy","white","firebrick3"))(10),
                            annotation_col = phen,annotation_row = split_row2,
                            row_split = split_row2$EDM,annotation_names_row = T,
                            annotation_names_col = T ,column_title = NULL,row_title = NULL)
p2


library(ggplot2)
library(ggpubr)
p3=ggplot(exp_bpm, aes(Group, `GCTCAC-6-motif_r`)) + 
    geom_boxplot(aes(fill=Group))+
    stat_compare_means(hide.ns = T)+
    xlab("")+ ylab("Percentage of GCTCAC")+geom_jitter(width = 0.3, alpha = 0.7, color = 'black')+
    scale_fill_manual(values = c("#0571B0","#CA0020","#00468BFF","#42B540FF","#0099B4FF","#925E9FFF","#FDAF91FF"))+theme(legend.position = "")
p3


Normal_ESCA=roc(as.factor(exp_bpm$Group), as.numeric(exp_bpm$`GCTCAC-6-motif_r`))
g=pROC::ggroc(list(`Normal VS Tumor`=Normal_ESCA))
p4=g+ggsci::scale_color_lancet()+geom_segment(aes(x=0,xend=1,y=1,yend=0),color="darkgrey",linetype=4)+theme_classic() +
    annotate("text", y = 0.4, x=0.5, label = paste0("motif GCTCAC \nAUC = ",round(Normal_ESCA$auc,4)),colour ="#00468BFF")+
    theme(legend.position = "")
p4




library(ggplot2)
library(ggpubr)
exp_edm$Group=phen_esca[match(rownames(exp_edm),phen_esca$Sample2),"Group"]
exp_edm$Group=factor(exp_edm$Group,levels = c("Normal","Tumor"))
p5=ggplot(exp_edm, aes(Group, `CGATCT-motif_r`)) + 
    geom_boxplot(aes(fill=Group))+
    stat_compare_means(hide.ns = T)+
    xlab("")+ ylab("Percentage of CGATCT")+geom_jitter(width = 0.3, alpha = 0.7, color = 'black')+
    scale_fill_manual(values = c("#0571B0","#CA0020","#00468BFF","#42B540FF","#0099B4FF","#925E9FFF","#FDAF91FF"))+theme(legend.position = "")
p5


Normal_ESCA=roc(as.factor(exp_edm$Group), as.numeric(exp_edm$`CGATCT-motif_r`))
g=pROC::ggroc(list(`Normal VS Tumor`=Normal_ESCA))
p6=g+ggsci::scale_color_lancet()+geom_segment(aes(x=0,xend=1,y=1,yend=0),color="darkgrey",linetype=4)+theme_classic() +
    annotate("text", y = 0.4, x=0.5, label = paste0("motif CGATCT \nAUC = ",round(Normal_ESCA$auc,4)),colour ="#00468BFF")+
    theme(legend.position = "")
p6


library(pROC)
library(rms)
load("medip_esca_peak_samp146_frag_comb_data_rf_modlist.RData")
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
round(ci(auc1)[1],4)

library(ggsci)
pal_lancet(palette = c("lanonc"), alpha = 0.6)(9)
c("#00468B99","#ED000099","#42B54099","#0099B499","#925E9F99","#FDAF9199","#AD002A99","#ADB6B699","#1B191999") 

p7=g+ggsci::scale_color_lancet()+
    geom_segment(aes(x=0,xend=1,y=1,yend=0),color="darkgrey",linetype=4)+
    theme_classic() + theme(legend.title = element_blank(),legend.position = "")+
    annotate("text", y = 0.75, x=0.5, label = paste0("Fragmentomics model"),colour ="black")+
    annotate("text", y = 0.6, x=0.5, label = paste0("Test set AUC = ",round(auc1$auc[[1]],4),"\n95% CI ( ",round(ci(auc1)[1],4)," - ",round(ci(auc1)[3],4)," )"),colour ="#00468B99")+
    annotate("text", y = 0.45, x=0.5, label = paste0("Validation set AUC = ",round(auc2$auc[[1]],4),"\n95% CI ( ",round(ci(auc2)[1],4)," - ",round(ci(auc2)[3],4)," )"),colour ="#ED000099")
p7



load("medip_esca_peak_samp136_valid_samp.RData")
load("medip_esca_peak_samp146_comb_data_rf_modlist.RData")
all_pre=data.frame()
for (i in 1:100) {
    all_pre=rbind(all_pre,data.frame(group2="Comb",limma_modlist[[i]]$limma_pre))
}

load("medip_esca_peak_samp146_frag_comb_data_rf_modlist.RData")
for (i in 1:100) {
    all_pre=rbind(all_pre,data.frame(group2="Frag",limma_modlist[[i]]$limma_pre))
}

load("medip_esca_peak_samp146_norm_voom_10_01_rf_modlist.RData")
for (i in 1:100) {
    all_pre=rbind(all_pre,data.frame(group2="Methy",limma_modlist[[i]]$limma_pre))
}

g=pROC::ggroc(list(`Methy`=roc(as.factor(all_pre[all_pre$group2 %in% "Methy","group"]), as.numeric(all_pre[all_pre$group2 %in% "Methy","pre"])),
                   `Frag`=roc(as.factor(all_pre[all_pre$group2 %in% "Frag","group"]), as.numeric(all_pre[all_pre$group2 %in% "Frag","pre"])),
                   `Comb`=roc(as.factor(all_pre[all_pre$group2 %in% "Comb","group"]), as.numeric(all_pre[all_pre$group2 %in% "Comb","pre"]))))
auc1=roc(as.factor(all_pre[all_pre$group2 %in% "Methy","group"]), as.numeric(all_pre[all_pre$group2 %in% "Methy","pre"]))
auc2=roc(as.factor(all_pre[all_pre$group2 %in% "Frag","group"]), as.numeric(all_pre[all_pre$group2 %in% "Frag","pre"]))
auc3=roc(as.factor(all_pre[all_pre$group2 %in% "Comb","group"]), as.numeric(all_pre[all_pre$group2 %in% "Comb","pre"]))

c("#00468B99","#ED000099","#42B54099","#0099B499","#925E9F99","#FDAF9199","#AD002A99","#ADB6B699","#1B191999") 
p8=g+ggsci::scale_color_lancet()+
    geom_segment(aes(x=0,xend=1,y=1,yend=0),color="darkgrey",linetype=4)+
    theme_classic() + theme(legend.title = element_blank(),legend.position = "")+
    annotate("text", y = 0.75, x=0.5, label = paste0("Combine Model AUC = ",round(auc3$auc[[1]],4),"\n95% CI ( ",round(ci(auc3)[1],4)," - ",round(ci(auc3)[3],4)," )"),colour ="#42B54099")+
    annotate("text", y = 0.6, x=0.5, label = paste0("Fragmentomics Model AUC = ",round(auc2$auc[[1]],4),"\n95% CI ( ",round(ci(auc2)[1],4)," - ",round(ci(auc2)[3],4)," )"),colour ="#ED000099")+
    annotate("text", y = 0.45, x=0.5, label = paste0("Methylation Model AUC = ",round(auc1$auc[[1]],4),"\n95% CI ( ",round(ci(auc1)[1],4)," - ",round(ci(auc1)[3],4)," )"),colour ="#00468B99")
p8



