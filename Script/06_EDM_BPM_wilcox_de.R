## end motif (EDM)     ###########################################
rm(list=ls())
gc()

load("medip_esca_samp146_phen.RData")
load("end_motif_all_freq.RData")
rownames(a)=a$Sample
a$Sample=NULL
a=a[,!grepl("^_",colnames(a))]
motif=colnames(a)[2:ncol(a)]

load("medip_esca_samp146_phen.RData")
load("medip_esca_peak_samp146_subsamp.RData")
table(rownames(a) %in% phen_esca$Sample2)
rownames(a)=phen_esca[match(rownames(a),phen_esca$Sample2),"Sample"]
a=a[subsample[[1]]$train_samp,]
all_de=data.frame()
for (i in colnames(a)[2:ncol(a)]) {
    all_de=rbind(all_de,data.frame(prop=i,p=wilcox.test(as.numeric(a[,i])~a$Group)$p.value))
}
save(all_de,file = "end_motif_all_de_res_samp146.RData")


## breakpoint motif    ###########################################
load("medip_esca_samp146_phen.RData")
load("medip_esca_peak_samp146_subsamp.RData")
load("bpm_end_motif_all_freq.RData")
a=bpm_motif_freq_all
rownames(a)=phen_esca[match(rownames(a),phen_esca$Sample2),"Sample"]
a=a[subsample[[1]]$train_samp,]
a$Group=factor(a$Group,levels = c("Normal","Tumor"))
a=a[,!grepl("^_",colnames(a))]
all_de=data.frame()
for (i in colnames(a)[1:(ncol(a)-1)]) {
    all_de=rbind(all_de,data.frame(prop=i,p=wilcox.test(as.numeric(a[,i])~a$Group)$p.value))
}
save(all_de,file = "bpm_end_motif_all_de_res_samp146.RData")



