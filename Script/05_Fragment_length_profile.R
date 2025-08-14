rm(list=ls())
gc()

find_peaks <- function (x, m = 3){
    shape <- diff(sign(diff(x, na.pad = FALSE)))
    pks <- sapply(which(shape < 0), FUN = function(i){
        z <- i - m + 1
        z <- ifelse(z > 0, z, 1)
        w <- i + m + 1
        w <- ifelse(w < length(x), w, length(x))
        if(all(x[c(z : i, (i + 2) : w)] <= x[i + 1])) return(i + 1) else return(numeric(0))
    })
    pks <- unlist(pks)
    pks
}

load("medip_esca_samp146_phen.RData")
load("length_all_samp146.RData")
table(length_all_samp$V3 %in% phen_esca$Sample2)
length_all_samp=length_all_samp[length_all_samp$V3 %in% phen_esca$Sample2,]
a=reshape2::dcast(length_all_samp,V2~Group,value.var = "V1",fill = 0,fun.aggregate = median)
a$ALL=a$Normal+a$Tumor
a$Normal=a$Normal/sum(a$Normal)
a$Tumor=a$Tumor/sum(a$Tumor)
a$ALL=a$ALL/sum(a$ALL)

peak=unique(c(20,find_peaks(a$ALL)+19,153,find_peaks(-a$ALL)+19))
peak=peak[peak<500]

a=reshape2::dcast(length_all_samp,V2~V3,value.var = "freq",fill = 0,fun.aggregate = sum)
rownames(a)=a$V2

all_exp=data.frame()
for (s in peak) {
    for (e in peak) {
        if(s+5 < e){
            prop=colSums(a[a$V2>=s & a$V2<e,])
            dat=data.frame(sample=names(prop),prop=as.numeric(prop),length=paste(s,e,sep = "_"))
            all_exp=rbind(all_exp,dat)
        }
    }
}
all_exp=reshape2::dcast(all_exp,length~sample,fill = 0,value.var = "prop",fun.aggregate = sum)
rownames(all_exp)=paste0("length_",all_exp$length)
all_exp$length=NULL
all_exp=all_exp[,colnames(all_exp) %in% phen_esca$Sample2]
colnames(all_exp)=phen_esca[match(colnames(all_exp),phen_esca$Sample2),"Sample"]
# save(all_exp,file = "frag_prop_all_exp_samp146.RData")


all_exp=as.data.frame(t(all_exp[,colnames(all_exp) %in% subsample[[1]]$train_samp]))
all_exp$Group=phen_esca[match(rownames(all_exp),phen_esca$Sample),"Group"]
all_exp$Group=factor(all_exp$Group,levels = c("Normal","Tumor"))
all_de=data.frame()
for (i in colnames(all_exp)[1:(ncol(all_exp)-1)]) {
    all_de=rbind(all_de,data.frame(prop=i,p=wilcox.test(as.numeric(all_exp[,i])~all_exp$Group)$p.value))
}
# save(all_de,file = "frag_prop_all_de_res_samp146.RData")
