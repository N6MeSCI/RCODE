library(ggplot2)
library(fgsea)
library(enrichplot)
library(clusterProfiler)
library(GSVA)
load("lncRNA_exp_anno_dat.rdata")
library(pheatmap)
library(ConsensusClusterPlus)
library(reshape2)
sci_exp<-norm_exp[,as.character(anno_hc_sci$sample[anno_hc_sci$type%in%"Spinal_Cord_Injury"])]
zero_a<-apply(sci_exp,1,function(x){
  sum(x>0)
})
sci_exp<-sci_exp[zero_a>0.2*28,]
# 
source("utils.R")
pathways <- gmtPathways("c2.cp.kegg.v7.4.symbols.gmt")


sig_path<-pathways[grep("_SIGNALING",names(pathways))]
sig_gene<-unique(unlist(sig_path))
gsva_result<-gsva(as.matrix(sci_exp), sig_path,
     min.sz=10, max.sz=500, verbose=TRUE)
results2 = ConsensusClusterPlus(as.matrix(gsva_result),maxK=10,reps=1000,pItem=0.8,pFeature=1,
                               title = "bb",clusterAlg="km",distance="pearson",seed=1000,plot="pdf")
consen_class<-results2[[4]][["consensusClass"]]
sci_anno<-data.frame(sample=names(consen_class),type=paste("class",consen_class,sep = ""))
sci_anno<-sci_anno[order(sci_anno$type),]
gsva_result<-gsva_result[,sci_anno$sample]

library(pheatmap)
table(colnames(gsva_result)==sci_anno$sample)

anno_col<-data.frame(Sample_type=as.character(sci_anno$type))
color<-RColorBrewer::brewer.pal(n = 11, name = 'Spectral')
row.names(anno_col)<-sci_anno$sample
anno_color<-list(Sample_type=color[c(1,3,6,9)])
names(anno_color[[1]])<-unique(as.character(sci_anno$type))

pheatmap(gsva_result,
         # scale = "row",
         cluster_rows = T,
         cluster_cols = F,
         show_colnames = FALSE,
         show_rownames = T,
         color = c(colorRampPalette(c("#3C7FAC", "white"))(50),colorRampPalette(c("white","#931F23"))(50)),
         annotation_col = anno_col,
         # annotation_row = anno_row,
         annotation_colors = anno_color,
         border_color = NA,
         cutree_cols = 4
)

site<-c(grep("MAPK",row.names(gsva_result)),
        grep("NOTCH",row.names(gsva_result)),
        grep("MTOR",row.names(gsva_result)),
        grep("WNT",row.names(gsva_result)))

colnames(gsva_result)<-sci_anno$type
all_p<-list()
for (i in 1:dim(gsva_result)[1]) {
  dat<-melt(as.matrix(gsva_result[i,]))
  fit <- aov(value ~ Var1, data=dat)
  pvalue<-summary(fit)[[1]][1,5]
  pvalue<-signif(pvalue,3)
  p<-ggplot(dat, aes(x=Var1, y=value,color=Var1))+
      # facet_wrap(~Var1,ncol = 5)+
      stat_boxplot(geom = 'errorbar',width=0.15)+
      geom_boxplot()+
      # stat_summary(fun.data = calc_stat, geom="boxplot")+
      geom_jitter(width = 0.2, alpha = 1, size=1)+
      ggtitle(row.names(gsva_result)[i]) +
      theme_classic() +
      theme(legend.position = "none",
            axis.text.x = element_text(angle = 0, hjust = 1)) +
      labs(x="", y="The activity of pathway")+
      annotate("text",x=2,y=max(dat$value),label = paste("one-way anova, p =",pvalue))+
      scale_color_manual(values = color[c(2,5,8,10)])
    all_p<-c(all_p,list(p))
  
}
library(cowplot)
plot_grid(all_p[[site[1]]],all_p[[site[2]]],all_p[[site[3]]],all_p[[site[4]]])
##########################pathway->gene
sci_exp<-sci_exp[,sci_anno$sample]
four_path<-unique(unlist(sig_path[site]))
four_path_gene<-diff_gene_hc_sci[row.names(diff_gene_hc_sci)%in%four_path,]
four_path_gene<-four_path_gene[four_path_gene$pvalue<0.05,]
four_path_exp2<-sci_exp[row.names(sci_exp)%in%row.names(four_path_gene),]
table(colnames(four_path_exp2)==sci_anno$sample)
colnames(four_path_exp2)<-sci_anno$type
four_path_exp2<-melt(as.matrix(four_path_exp2))
four_path_exp2<-split(four_path_exp2,four_path_exp2$Var1)

all_p<-vector()
for (j in 1:length(four_path_exp2)) {
  dat2<-four_path_exp2[[j]]
  fit <- aov(value ~ Var2, data=dat2)
  pvalue<-summary(fit)[[1]][1,5]
  all_p<-c(all_p,pvalue)
}
names(all_p)<-names(four_path_exp2)
all_p<-all_p[all_p<0.05]

four_path_exp<-sci_exp[row.names(sci_exp)%in%names(all_p),]
load("selectgenes.rdata")
m6a_exp<-as.data.frame(sci_exp[row.names(selectgenes),])
table(colnames(four_path_exp)==colnames(m6a_exp))

person_cor<-apply(m6a_exp,1,function(x){
  aa=apply(four_path_exp,1,function(y){
    dd<-cor.test(x,y)
    list(dd[["p.value"]],dd[["estimate"]][["cor"]])
  })
  aa[unlist(lapply(aa,function(x){
    x[[1]]<0.01
  }))]
})
save(person_cor,file = "person_cor_pathway.rdata")
aa<-person_cor
all_cor_dat<-vector()
for (i in 1:length(aa)) {
  dat<-aa[[i]]
  if(length(dat)!=0){
    cor_dat<-data.frame(inter_gene=names(dat),cor=unlist(lapply(dat,function(x){x[[2]]})),row.names = NULL)
    cor_dat$select_gene<-names(aa)[i]
    all_cor_dat<-rbind(all_cor_dat,cor_dat)
  }
}
all_cor_dat<-all_cor_dat[abs(all_cor_dat$cor)>0.7 & all_cor_dat$cor<1,]
all_cor_dat<-all_cor_dat[all_cor_dat$inter_gene!=all_cor_dat$select_gene,]

######################Drug Target Data Research
gene_drug<-read.csv("drug_gene.txt",header = T,sep = "\t")
all_dat_drug<-vector()
for (i in 1:dim(all_cor_dat)[1]) {
  aa<-grep(all_cor_dat$inter_gene[i],gene_drug$V2)
  if(length(aa)!=0){
    gene_drug2<-gene_drug[aa,]
    colnames(gene_drug2)<-c("drug","target_gene")
    new_dat<-data.frame(m6a_gene=rep(all_cor_dat[i,3],dim(gene_drug2)[1]),
                        target=rep(all_cor_dat[i,1],dim(gene_drug2)[1]),
                        cor=rep(all_cor_dat[i,2],dim(gene_drug2)[1]))
    new_dat<-cbind(new_dat,gene_drug2)
    all_dat_drug<-rbind(all_dat_drug,new_dat)
  }
}

inter_dat1<-all_dat_drug[,c(1,2)]
inter_dat2<-all_dat_drug[,c(2,4)]
colnames(inter_dat1)<-c("a","b")
colnames(inter_dat2)<-c("a","b")
inter_dat<-rbind(inter_dat1,inter_dat2)
inter_dat<-inter_dat[!duplicated(inter_dat),]
inter_anno<-data.frame(symbol=c(unique(all_dat_drug$m6a_gene),unique(all_dat_drug$target),unique(all_dat_drug$drug)),
                       type=c(rep("m6a_gene",length(unique(all_dat_drug$m6a_gene))),
                              rep("target",length(unique(all_dat_drug$target))),
                              rep("drug",length(unique(all_dat_drug$drug)))))
write.table(inter_dat,"inter_drug_dat.txt",sep = "\t",row.names = F,quote = F)
write.table(inter_anno,"inter_drug_anno.txt",sep = "\t",row.names = F,quote = F)
