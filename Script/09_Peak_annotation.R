rm(list=ls())
options(stringsAsFactors = F)

load("medip_esca_peak_samp146_exp.RData")
peak_ana=data.frame(pos=rownames(peak),max=rowSums(peak),count=rowSums(peak>0))
peak_ana=peak_ana[!grepl("chrUn|random|hap",peak_ana$pos),]

require(ChIPseeker)
require(GenomicRanges)
require(TxDb.Hsapiens.UCSC.hg19.knownGene)
require(org.Hs.eg.db)
require(tidyverse)
txdb=TxDb.Hsapiens.UCSC.hg19.knownGene

pos=do.call(rbind,lapply(as.character(peak_ana$pos), function(x){unlist(strsplit(x,"_"))})) %>% as.data.frame()
pos=na.omit(data.frame(seqnames=pos$V1,start=as.numeric(pos$V2),end=as.numeric(pos$V3)))
peakAnno = GRanges(seqnames=Rle(pos[,1]),ranges=IRanges(pos[,2], pos[,3]), strand=rep(c("*"), nrow(pos))) %>% annotatePeak(. ,TxDb=txdb, annoDb="org.Hs.eg.db") %>% as.data.frame()
peakAnno$annotation1=str_split(peakAnno$annotation,'[(]',simplify = T)[,1]
table(peakAnno$annotation1)

peakAnno$pos=paste(peakAnno$seqnames,peakAnno$start,peakAnno$end,sep = "_")
peak_ana=merge(peak_ana,peakAnno[,c(6,14:16,18:19)],by="pos")

hgnc_complete_set <- read.delim("E:/database/database/genome_annotation/hgnc_complete_set.txt")
peak_ana$type=hgnc_complete_set[match(peak_ana$SYMBOL,hgnc_complete_set$symbol),"locus_group"]

gpl = GRanges(seqnames=Rle(pos$seqnames),ranges=IRanges(pos$start, pos$end), strand=rep(c("*"), nrow(pos))) 
annots = c( 'hg19_basicgenes', 'hg19_genes_intergenic','hg19_genes_intronexonboundaries')
annots = build_annotations(genome = 'hg19', annotations = annots)
annots = annotate_regions(regions = gpl,annotations = annots,ignore.strand = TRUE,quiet = FALSE)
print(annots)
annots$annot$type
save(annots,file="hg19_basicgenes_anno_all_peak.RData")

annots1=c('hg19_cpgs','hg19_cpgs_islands','hg19_cpgs_shores','hg19_cpgs_shelves','hg19_cpgs_inter')
annots1 = build_annotations(genome = 'hg19', annotations = annots1)
annots1 = annotate_regions(regions = gpl,annotations = annots1,ignore.strand = TRUE,quiet = FALSE)
print(annots1)
annots1$annot$type
save(annots1,file="hg19_cpgs_anno_all_peak.RData")

## repeat_ucsc_table_hg19 
extraCols = c(repName = 'character', repClass = 'character', repFamily = 'character')
dm_regions = read_regions(con = "repeat_ucsc_table_hg19.bed", genome = 'hg19', extraCols = extraCols, format = 'bed')
repeat_annotated = annotate_regions(regions = gpl,annotations = dm_regions,ignore.strand = TRUE,quiet = FALSE)
save(repeat_annotated,file="repeat_annotated_all_peak.RData")

## enhancer_ucsc_table_hg19 
enhancer_regions = read_regions(con = "enhancer_ucsc_table_hg19.bed", genome = 'hg19', format = 'bed')
enhancer_annotated = annotate_regions(regions = gpl,annotations = enhancer_regions,ignore.strand = TRUE,quiet = FALSE)
save(enhancer_annotated,file="enhancer_annotated_all_peak.RData")

## encode_hg19 
extraCols2 = c(class = 'character')
encode_regions = read_regions(con = "hg19_encode.bed", genome = 'hg19', extraCols = extraCols2 ,format = 'bed')
encode_annotated = annotate_regions(regions = gpl,annotations = encode_regions,ignore.strand = TRUE,quiet = FALSE)
save(encode_annotated,file="encode_annotated_all_peak.RData")


load("hg19_basicgenes_anno_all_peak.RData")
annots=as.data.frame(annots)
annots=annots[,c(1:3,11:15)]
annots$pos=paste(annots$seqnames,annots$start,annots$end,sep = "_")

load("encode_annotated_all_peak.RData")
encode_annotated=as.data.frame(encode_annotated)
encode_annotated=encode_annotated[,c(1:3,11)]
encode_annotated$pos=paste(encode_annotated$seqnames,encode_annotated$start,encode_annotated$end,sep = "_")

load("enhancer_annotated_all_peak.RData")
enhancer_annotated=as.data.frame(enhancer_annotated)
enhancer_annotated=enhancer_annotated[,c(1:3)]
enhancer_annotated$enhancer_type="enhancer"
enhancer_annotated$pos=paste(enhancer_annotated$seqnames,enhancer_annotated$start,enhancer_annotated$end,sep = "_")


load("repeat_annotated_all_peak.RData")
repeat_annotated=as.data.frame(repeat_annotated)
repeat_annotated=repeat_annotated[,c(1:3,11:13)]
repeat_annotated$pos=paste(repeat_annotated$seqnames,repeat_annotated$start,repeat_annotated$end,sep = "_")


load("hg19_cpgs_anno_all_peak.RData")
annots1=as.data.frame(annots1)
annots1=annots1[,c(1,2,3,15)]
annots1$pos=paste(annots1$seqnames,annots1$start,annots1$end,sep = "_")

load("medip_esca_peak_samp146_peakAnno.RData")
peak_ana=merge(peak_ana,annots1[,4:5],by="pos",all=T)
peak_ana=merge(peak_ana,enhancer_annotated[,4:5],by="pos",all=T)
peak_ana=merge(peak_ana,encode_annotated[,4:5],by="pos",all=T)
peak_ana=merge(peak_ana,repeat_annotated[,4:7],by="pos",all=T)
peak_ana=merge(peak_ana,annots[,4:9],by="pos",all=T)
colnames(peak_ana)

peak_ana$annot.tx_id=NULL
peak_ana$annot.gene_id=NULL
peak_ana$annot.id=NULL
# save(peak_ana,file="medip_esca_peak_samp146_peakAnno.RData")

