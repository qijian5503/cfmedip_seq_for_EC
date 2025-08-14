rm(list=ls())
gc()

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


library(pheatmap)
library(gplots)
library(colorRamps)
library(RColorBrewer)
load("medip_esca_peak_samp146_traintest_norm_voom_10_01_noxy.RData")
load("medip_esca_peak_samp146_norm_voom_10_01_valid_noxy.RData")
norm=cbind(norm,norm_v)

nmf_exp=norm[de_gene,]
phen=phen_esca[phen_esca$Sample %in% colnames(nmf_exp),]
colnames(phen)
phen=phen[,c(5,6,3,2)]
phen=phen[order(phen$Group),]

phen$Group=factor(phen$Group,levels = c("Normal","Tumor"))
phen$Stage=factor(phen$Stage,levels = c("Normal","StageI","StageII","StageIII"))
median(phen$Age,na.rm = T)
phen$Age=ifelse(phen$Age>=52,">=52",ifelse(phen$Age<52,"<52",NA))
phen$Age=factor(phen$Age,levels = c("<52",">=52"))

ann_colors = list(
    Group=c(Tumor="#CA0020", Normal="#0571B0"),
    change=c(Hyper="#CA0020", Hypo="#0571B0"),
    Gender=c(`F`="#FB9A99",`M`="#A65628"),
    Age=c(`>=52`="#4DAF4A",`<52`="#984EA3"),
    Stage=c(`Normal`="#0571B0",`StageI`="#C7EAE5",`StageII`="#01665E",`StageIII`="#003C30")
)
phen1=data.frame(row.names = rownames(phen),Group=phen$Group)

imp_de_gene$change=ifelse(imp_de_gene$Hypermethylation>0,"Hyper",
                          ifelse(imp_de_gene$Hypomethylation>0,"Hypo",NA))
colnames(imp_de_gene)[1]="pos"
split_row=unique(imp_de_gene[de_gene,c("change","pos")])
rownames(split_row)=split_row$pos
split_row$pos=NULL

p1=ComplexHeatmap::pheatmap(as.matrix(nmf_exp[,rownames(phen)]), scale = "row",cluster_cols = F,
                            show_rownames = F,show_colnames = F,annotation_colors = ann_colors,
                            col = colorRampPalette(c("navy","white","firebrick3"))(10),
                            annotation_col = phen,annotation_row = split_row,
                            row_split = split_row$change,annotation_names_row = T,
                            annotation_names_col = T ,column_title = NULL,row_title = NULL)
p1


load("medip_esca_peak_samp146_peakAnno.RData")
de_gene=peak_ana[peak_ana$pos %in% de_gene,]
de_gene$regule=imp_de_gene[match(de_gene$pos,imp_de_gene$Var1),"change"]
colnames(de_gene)

de_gene2=unique(de_gene[,c("pos","regule","annotation1")])
b= de_gene2[,c("regule","annotation1")] %>% dplyr::group_by(regule,annotation1) %>% 
    dplyr::count() %>% group_by(regule) %>% mutate(Freq=n/sum(n)*100)
b$n_lab=ifelse(b$Freq<5,"",b$n)
b$Freq_lab=ifelse(b$Freq<5,"",paste(round(b$Freq,2),"%"))
b$group=gsub("methylation","",b$regule)
p_value2=chisq.test(table(de_gene2$regule,de_gene2$annotation1))$p.value
p2=ggplot(b,aes(x=regule,y=Freq,fill=annotation1))+geom_col()+geom_text(aes(label=Freq_lab),position=position_stack(vjust=0.5))+theme_classic()+
    scale_fill_manual(values=c("#E64B35FF","#4DBBD5FF","#00A087FF","#3C5488FF","#F39B7FFF","#8491B4FF","#91D1C2FF","#DC0000FF","#7E6148FF","#B09C85FF"))+
    labs(subtitle = paste0("p = ",format(p_value2,scientific=T,digits=3)),title="", x="", y="Frenquency")+theme(legend.title = element_blank(),axis.text.x = element_text(angle = 45, hjust = 0.1, vjust = 0))
p2


de_gene2=unique(de_gene[,c("pos","regule","annot.type.x")])
d= de_gene2[,c("regule","annot.type.x")] %>% dplyr::group_by(regule,annot.type.x) %>% 
    dplyr::count() %>% group_by(regule) %>% mutate(Freq=n/sum(n)*100)
d$n_lab=ifelse(d$Freq<5,"",d$n)
d$Freq_lab=ifelse(d$Freq<5,"",paste(round(d$Freq,2),"%"))
p_value4=chisq.test(table(de_gene2$regule,de_gene2$annot.type.x))$p.value
p3=ggplot(d,aes(x=regule,y=Freq,fill=annot.type.x))+geom_col()+geom_text(aes(label=Freq_lab),position=position_stack(vjust=0.5))+theme_classic()+
    scale_fill_manual(values=c("#E64B35FF","#4DBBD5FF","#00A087FF","#3C5488FF","#F39B7FFF","#8491B4FF","#91D1C2FF","#DC0000FF","#7E6148FF","#B09C85FF"))+
    labs(subtitle = paste0("p = ",format(p_value4,scientific=T,digits=3)),title="", x="", y="Frenquency")+theme(legend.title = element_blank(),axis.text.x = element_text(angle = 45, hjust = 0.1, vjust = 0))
p3



down_gene=unique(de_gene[de_gene$regule %in% "Hypo","SYMBOL"])
up_gene=unique(de_gene[de_gene$regule %in% "Hyper","SYMBOL"])

library(clusterProfiler)
library(org.Hs.eg.db)
ensem_down=bitr(down_gene,fromType = "SYMBOL",toType = "ENTREZID",OrgDb = "org.Hs.eg.db")
go_down=enrichGO(ensem_down$ENTREZID,OrgDb = org.Hs.eg.db,ont = "ALL",pAdjustMethod = "none",pvalueCutoff = 0.1,qvalueCutoff = 0.2,keyType = "ENTREZID")
kegg_down=enrichKEGG(ensem_down$ENTREZID, organism = 'hsa', keyType = 'kegg',pAdjustMethod = "none", pvalueCutoff = 0.1,qvalueCutoff = 0.2)
dotplot(kegg_down, showCategory=15) 
kegg_down = setReadable(kegg_down, OrgDb = org.Hs.eg.db, keyType="ENTREZID") ## 显示基因symbol

ensem_up=bitr(up_gene,fromType = "SYMBOL",toType = "ENTREZID",OrgDb = "org.Hs.eg.db")
go_up=enrichGO(ensem_up$ENTREZID,OrgDb = org.Hs.eg.db,ont = "ALL",pAdjustMethod = "none",pvalueCutoff = 0.1,qvalueCutoff = 0.2,keyType = "ENTREZID")
kegg_up=enrichKEGG(ensem_up$ENTREZID, organism = 'hsa', keyType = 'kegg',pAdjustMethod = "none", pvalueCutoff = 0.1,qvalueCutoff = 0.2)
dotplot(kegg_up, showCategory=15) 
kegg_up = setReadable(kegg_up, OrgDb = org.Hs.eg.db, keyType="ENTREZID") ## 显示基因symbol

enrich_res=list(go_down=go_down,kegg_down=kegg_down,go_up=go_up,kegg_up=kegg_up)
save(enrich_res,file = "enrich_up_down_res.RData")


load("enrich_up_down_res.RData")
go_up=enrich_res$go_up
kegg_up=enrich_res$kegg_up
go_down=enrich_res$go_down
kegg_down=enrich_res$kegg_down
p4=dotplot(go_up, split="ONTOLOGY", showCategory=5)+facet_grid(ONTOLOGY~., scale="free")
p4
p5=dotplot(go_down, split="ONTOLOGY", showCategory=5)+facet_grid(ONTOLOGY~., scale="free")
p5

