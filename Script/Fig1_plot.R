rm(list=ls())
options(stringsAsFactors = F)

load("medip_esca_samp146_phen.RData")
colnames(phen_esca)=c("Sample","Age","Gender","cfDNA Concentration(ng/ml)","Group","Stage","Raw_Reads","Raw_Data","Freq_bed","Sample2")
table(phen_esca$Stage)
phen_esca$Stage=ifelse(phen_esca$Group %in% "Normal","Normal",phen_esca$Stage)

library(ggpubr)
phen_esca$Group=factor(phen_esca$Group,levels = c("Normal","Tumor"))
phen_esca$Stage=factor(phen_esca$Stage,levels = c("Normal","StageI","StageII","StageIII"))
median(phen_esca$Age,na.rm = T)
phen_esca$Age=ifelse(phen_esca$Age>=52,">=52",ifelse(phen_esca$Age<52,"<52",NA))
phen_esca$Age=factor(phen_esca$Age,levels = c("<52",">=52"))
phen_esca$Raw_Reads=phen_esca$Raw_Reads/1000000
colnames(phen_esca)[7]=c("Reads Number(Millions)")
colnames(phen_esca)[2]="Age(years)"
# save(phen_esca,file = "medip_esca_samp146_phen2.RData")


load("medip_esca_samp146_phen2.RData")

p2=ggplot(phen_esca[!is.na(phen_esca$Group),], aes(Group, `cfDNA Concentration(ng/ml)`)) + 
    geom_boxplot(aes(fill=Group))+xlab(label = "")+
    scale_fill_manual(values=c("#0571B0","#CA0020"))+
    stat_compare_means(aes(group=Group),hide.ns = T,label = "p")+labs(title = "")+ylab("cfDNA Concentration(ng/ml)")+
    geom_jitter(width = 0.3, alpha = 0.7, color = 'black')+
    theme(axis.text.x = element_text(angle = 90, hjust = 0.5, vjust = 0.5))
p2

p3=ggplot(phen_esca[!is.na(phen_esca$Gender) & !is.na(phen_esca$Group),], aes(Group, `cfDNA Concentration(ng/ml)`)) + 
    geom_boxplot(aes(fill=Gender))+xlab(label = "")+
    stat_compare_means(aes(group=Gender),hide.ns = T,label = "p")+labs(title = "")+ylab("cfDNA Concentration(ng/ml)")+
    geom_jitter(width = 0.3, alpha = 0.7, color = 'black')+
    theme(axis.text.x = element_text(angle = 90, hjust = 0.5, vjust = 0.5))+
    scale_fill_manual(values=c("#0571B0","#CA0020","#ffad60"))
p3
p4=ggplot(phen_esca[!is.na(phen_esca$`Age(years)`) & !is.na(phen_esca$Group),], aes(Group, `cfDNA Concentration(ng/ml)`)) + 
    geom_boxplot(aes(fill=`Age(years)`))+xlab(label = "")+
    stat_compare_means(aes(group=`Age(years)`),hide.ns = T,label = "p")+labs(title = "")+ylab("cfDNA Concentration(ng/ml)")+
    geom_jitter(width = 0.3, alpha = 0.7, color = 'black')+
    theme(axis.text.x = element_text(angle = 90, hjust = 0.5, vjust = 0.5))+
    scale_fill_manual(values=c("#0571B0","#CA0020","#ffad60"))
p4
table(phen_esca$Stage)
my_comparisons=list(c("StageII","StageI"),c("StageIII","StageI"),c("StageIII","StageII"),c("StageI","Normal"),c("StageII","Normal"),c("StageIII","Normal"))
p5=ggplot(phen_esca[!is.na(phen_esca$Stage) ,], aes(Stage, `cfDNA Concentration(ng/ml)`)) + 
    geom_boxplot(aes(fill=Stage))+xlab(label = "")+
    stat_compare_means(comparisons = my_comparisons,aes(group=Stage),hide.ns = T)+labs(title = "")+ylab("cfDNA Concentration(ng/ml)")+
    geom_jitter(width = 0.3, alpha = 0.7, color = 'black')+
    theme(axis.text.x = element_text(angle = 90, hjust = 0.5, vjust = 0.5),legend.position = "")+
    scale_fill_manual(values = c("#0571B0","#CA0020","#00468BFF","#42B540FF","#0099B4FF","#925E9FFF","#FDAF91FF"))
p5


