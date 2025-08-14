## fragment length     ########################
rm(list=ls())
gc()

load("medip_esca_samp146_phen.RData")
load("length_all_samp146.RData")

length_all_samp$Group=factor(length_all_samp$Group,levels = c("Normal","Tumor"))
p1=ggplot(length_all_samp[length_all_samp$V2<500,], aes(V2, freq,fill=V3,color=Group))+geom_line()+xlab("CfDNA length")+ylab("Freq")
p1

a=reshape2::dcast(length_all_samp,V2~Group,value.var = "V1",fill = 0,fun.aggregate = sum)
a$Tumor=a$Tumor/sum(a$Tumor)
a$Normal=a$Normal/sum(a$Normal)
a=reshape2::melt(a,id.var="V2")
colnames(a)=c("CfDNA length","Group","Freq")
p2=ggplot(a[a$`CfDNA length`<100 & !a$`Group` %in% "ALL",], aes(`CfDNA length`, Freq,fill=`Group`,color=`Group`))+geom_line()+scale_color_manual(values=c("#00468BFF","#CA0020","#42B540FF","#0099B4FF","#925E9FFF","#FDAF91FF"))
p2




b=reshape2::dcast(length_all_samp,V2~V3,value.var = "freq",fill = 0,fun.aggregate = sum)
dat=data.frame(Group=phen_esca[match(colnames(b)[2:ncol(b)],phen_esca$Sample2),"Group"],
               prop=colSums(b[b$V2>=20 & b$V2<61,2:ncol(b)]))
dat$Group=factor(dat$Group,levels = c("Normal","Tumor"))
library(ggplot2)
library(ggpubr)
p3=ggplot(dat, aes(Group, prop)) + 
    geom_boxplot(aes(fill=Group))+
    stat_compare_means(hide.ns = T)+
    xlab("")+ ylab("Percentage of fragment length below 60")+geom_jitter(width = 0.3, alpha = 0.7, color = 'black')+
    scale_fill_manual(values = c("#0571B0","#CA0020","#00468BFF","#42B540FF","#0099B4FF","#925E9FFF","#FDAF91FF"))+theme(legend.position = "")
p3


Normal_ESCA=roc(as.factor(dat$Group), as.numeric(dat$prop))
g=pROC::ggroc(list(`Normal VS Tumor`=Normal_ESCA))
p4=g+ggsci::scale_color_lancet()+geom_segment(aes(x=0,xend=1,y=1,yend=0),color="darkgrey",linetype=4)+theme_classic() +
    annotate("text", y = 0.4, x=0.5, label = paste0("Normal VS Tumor AUC = ",round(Normal_ESCA$auc,4)),colour ="#00468BFF")+
    theme(legend.position = "")
p4


