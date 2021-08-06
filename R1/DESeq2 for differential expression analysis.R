# BiocManager::install("DESeq2")
library(DESeq2)

gse15_anno<-read.csv("GSE15_anno.txt",sep = "\t",header = F)
colnames(gse15_anno)<-c("sample","type")
dat<-read.csv("GSE151371_raw_gene_counts_de-ID.csv",header = T,sep = ",",row.names = 1)
 bb<-apply(dat,1,function(x){
   sum(x>0)
 })
dat<-dat[bb>0,]
row.names(gse15_anno)<-gse15_anno$sample
########################hc_sci
anno_hc_sci<-gse15_anno[gse15_anno$type%in%c("Control","Spinal_Cord_Injury"),]
anno_hc_sci$type<-factor(anno_hc_sci$type,levels = c(unique(anno_hc_sci$type)))
dat_hc_sci<-dat[,row.names(anno_hc_sci)]
table(row.names(anno_hc_sci)==colnames(dat_hc_sci))
dd<-apply(dat_hc_sci,1,function(x){sum(x>0)})
dat_hc_sci<-dat_hc_sci[dd>5,]
dds_hc_sci <- DESeqDataSetFromMatrix(dat_hc_sci, anno_hc_sci, design= ~ type)
dds_hc_sci <- DESeq(dds_hc_sci)
rld <- rlogTransformation(dds_hc_sci)
norm_exp<-assay(rld)
# save(norm_exp,anno_hc_sci,file = "hc_sci_anno_exp.rdata")
res_hc_sci = results(dds_hc_sci, contrast=c("type","Control","Spinal_Cord_Injury"))
res_hc_sci = res_hc_sci[order(res_hc_sci$pvalue),]
head(res_hc_sci)
summary(res_hc_sci)
diff_gene_hc_sci<-as.data.frame(res_hc_sci)
type<-rep("Normal",dim(diff_gene_hc_sci)[1])
type[diff_gene_hc_sci$log2FoldChange>1&diff_gene_hc_sci$padj<0.05]<-"Down"
type[diff_gene_hc_sci$log2FoldChange< -1&diff_gene_hc_sci$padj<0.05]<-"Up"
diff_gene_hc_sci$type<-type
save(diff_gene_hc_sci,file = "hc_sci_diff_gene.rdata")

#######################################Draw a volcano map
library(ggplot2)
library(ggrepel)
library(ggthemes)
library(RColorBrewer)
color<-brewer.pal(n = 11, "RdYlGn")
plot_mode <- "advanced"
logFCcut <- log2(1.5) 
logFCcut2 <- 2 
logFCcut3 <- 4
pvalCut <- 0.05
pvalCut2 <- 0.001 #for advanced mode
pvalCut3 <- 0.00001 #for advanced mode
adjPcut <- 0.05
x<-as.data.frame(res_hc_sci)
x$label<- rownames(x)
colnames(x)[c(2,5,6)]<-c("logFC","P.Value","adj.P.Val")
x$logFC<-x$logFC*(-1)
x<-x[!is.na(x$logFC),]
x<-x[!is.na(x$P.Value),]
x<-x[!is.na(x$adj.P.Val),]
re.adj <- x[x$adj.P.Val <adjPcut & (x$logFC > logFCcut | x$logFC < -logFCcut),] 
re.p = x[x$P.Val <pvalCut & (x$logFC > logFCcut | x$logFC < -logFCcut),]

if (plot_mode == "classic"){
  x[,6] <- ifelse((x$P.Value < pvalCut & x$logFC > logFCcut), "red", ifelse((x$P.Value < pvalCut & x$logFC < -logFCcut), "blue","grey30"))
  size <- ifelse((x$P.Value < pvalCut & abs(x$logFC) > logFCcut), 4, 2)
  
} else if (plot_mode == "advanced") {
  n1 <- length(x[,1])
  cols <- rep("grey30",n1)
  names(cols)<- rownames(x)
  cols[x$P.Value < pvalCut & x$logFC >0]<- color[4]
  cols[x$P.Value < pvalCut & x$logFC >logFCcut]<- color[3]
  cols[x$P.Value < pvalCut2 & x$logFC > logFCcut2]<- color[2]
  
  cols[x$P.Value < pvalCut & x$logFC < 0]<- color[8]
  cols[x$P.Value < pvalCut & x$logFC < -logFCcut]<- color[9]
  cols[x$P.Value < pvalCut2 & x$logFC < -logFCcut2]<- color[10]
  #cols[names(cols)==genelist]<- "red"
  color_transparent <- adjustcolor(cols, alpha.f = 0.5)  # set color transparence
  x[,6] <- color_transparent
  n1 <- length(x[,1])
  size <- rep(1,n1)
  size[x$P.Value < pvalCut & x$logFC > logFCcut]<- 1
  size[x$P.Value < pvalCut2 & x$logFC > logFCcut2]<- 2
  size[x$P.Value < pvalCut3 & x$logFC > logFCcut3]<- 5
  size[x$P.Value < pvalCut & x$logFC < -logFCcut]<- 1
  size[x$P.Value < pvalCut2 & x$logFC < -logFCcut2]<- 2
  size[x$P.Value < pvalCut3 & x$logFC < -logFCcut3]<- 5
  
} else {
  stop("Unsupport mode")
}
ggplot(data=x, aes(x=logFC, y=-log10(P.Value), label = label)) +
  geom_point(alpha = 0.6, size=size, colour=x[,6]) +
  scale_color_manual(values = c("lightgrey", "navy", "red"))

###################################Draw a heat map

library(reshape2)
library(ggsignif)
library(ggsci)
library(RColorBrewer)
library(ggplot2)
library(pheatmap)
color<-brewer.pal(n = 11, "RdYlGn")
load("selectgenes.rdata")

m6a_exp<-as.data.frame(norm_exp[row.names(selectgenes),])
table(row.names(anno_hc_sci)==colnames(m6a_exp))
colnames(m6a_exp)<-anno_hc_sci$type

dat_plot<-melt(as.matrix(m6a_exp))
colnames(dat_plot)[1:2]<-c("symbol","Type")
ggplot(dat_plot, aes(x=Type, y=value,color=Type))+
  facet_wrap(~symbol,ncol = 3)+
  stat_boxplot(geom = 'errorbar',width=0.15)+
  geom_boxplot()+
  # stat_summary(fun.data = calc_stat, geom="boxplot")+
  geom_jitter(width = 0.2, alpha = 0.5, size=1)+
  geom_signif(comparisons = list(c("Control","Spinal_Cord_Injury")),step_increase=0.05,map_signif_level = F,test = wilcox.test,color="black")+
  # theme_classic()+
  scale_color_manual(values = color[c(2,10)])+
  labs(x="Sample type",y= "The expression of genes")+
  theme(legend.position = "none")+
  theme_bw()


anno_col<-data.frame(Sample_type=as.character(anno_hc_sci$type))
color<-RColorBrewer::brewer.pal(n = 11, name = 'Spectral')
row.names(anno_col)<-colnames(m6a_exp)
anno_color<-list(Sample_type=color[c(1,9)])
names(anno_color[[1]])<-unique(as.character(anno_hc_sci$type))
m6a_exp<-m6a_exp[,order(m6a_exp[1,],decreasing = T)]
pheatmap(m6a_exp,
         scale = "row",
         cluster_rows = T,
         cluster_cols = T,
         show_colnames = FALSE,
         show_rownames = T,
         color = c(colorRampPalette(c("#3C7FAC", "white"))(50),colorRampPalette(c("white","#931F23"))(50)),
         annotation_col = anno_col,
         # annotation_row = anno_row,
         annotation_colors = anno_color,
         border_color = NA,
         cutree_cols = 3
)