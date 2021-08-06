library(ggplot2)

tf_dat<-read.csv("Homo_sapiens_TF.txt",sep = "\t",header = T)
load("m6a_cor_dat.rdata")
k=4

tf<-all_cor_dat[all_cor_dat$inter_gene%in%tf_dat$Symbol,]
tf<-split(tf,tf$select_gene)
tf_fto<-tf[[k]]
fto_cor_dat<-all_cor_dat[all_cor_dat$select_gene%in%names(tf)[k],]
tf_fto$type<-"positive"
tf_fto$type[tf_fto$cor<0]<-"negative"

po_tf<-as.vector(tf_fto$cor[tf_fto$cor>0])
names(po_tf)<-tf_fto$inter_gene[tf_fto$cor>0]
po_tf<-sort(po_tf,decreasing = T)
po_tf_dat<-data.frame(symbol=names(po_tf),num=as.vector(po_tf),type=rep("a",length(po_tf)))
po_tf_dat$symbol<-factor(po_tf_dat$symbol,levels = rev(names(po_tf)))
p1<-ggplot(po_tf_dat,aes(x=symbol,y=num,fill=type))+
  geom_bar(stat  = "identity")+
  coord_flip()+
  scale_fill_manual(values = "#F8AA8E")+
  theme_bw()+
  labs(x="",y="Correlation")+
  theme(panel.grid.major = element_blank(),   
        panel.grid.minor = element_blank(),
        panel.border = element_blank(),
        axis.line = element_line(colour = "black",size=1)
  )+guides(fill=FALSE)

ne_tf<-as.vector(tf_fto$cor[tf_fto$cor<0])
names(ne_tf)<-tf_fto$inter_gene[tf_fto$cor<0]
ne_tf<-sort(ne_tf,decreasing = T)
ne_tf_dat<-data.frame(symbol=names(ne_tf),num=as.vector(ne_tf),type=rep("a",length(ne_tf)))
ne_tf_dat$symbol<-factor(ne_tf_dat$symbol,levels = names(ne_tf))
p2<-ggplot(ne_tf_dat,aes(x=symbol,y=num,fill=type))+
  geom_bar(stat  = "identity")+
  coord_flip()+
  scale_fill_manual(values = "#9A9BCC")+
  theme_bw()+
  labs(x="",y="Correlation")+
  theme(panel.grid.major = element_blank(),   
        panel.grid.minor = element_blank(),
        panel.border = element_blank(),
        axis.line = element_line(colour = "black",size=1)
  )+guides(fill=FALSE)+
  scale_x_discrete(position = "top")
plot_grid(p1,p2,nrow = 1)
#######################
load("TF_TG_interaction.rdata")
po_inter<-tf_tg[tf_tg$TF.Name%in%names(po_tf),]
po_inter2<-po_inter[po_inter$TG.Name%in%fto_cor_dat$inter_gene[fto_cor_dat$cor>0],]
po_inter2<-po_inter2[!po_inter2$TG.Name%in%names(po_tf),]
po_inter_anno<-data.frame(symbol=c(unique(po_inter2$TF.Name),unique(po_inter2$TG.Name)),
                          type=c(rep("po_TF",length(unique(po_inter2$TF.Name))),rep("po_TG",length(unique(po_inter2$TG.Name)))),
                          cor=as.vector(fto_cor_dat$cor[match(c(unique(po_inter2$TF.Name),unique(po_inter2$TG.Name)),fto_cor_dat$inter_gene)]))
# write.table(po_inter2,"po_inter.txt",sep = "\t",quote = F,row.names = F)
# write.table(po_inter_anno,"po_inter_anno.txt",sep = "\t",quote = F,row.names = F)

ne_inter<-tf_tg[tf_tg$TF.Name%in%names(ne_tf),]
ne_inter2<-ne_inter[ne_inter$TG.Name%in%fto_cor_dat$inter_gene[fto_cor_dat$cor<0],]
ne_inter2<-ne_inter2[!ne_inter2$TG.Name%in%names(ne_tf),]
ne_inter_anno<-data.frame(symbol=c(unique(ne_inter2$TF.Name),unique(ne_inter2$TG.Name)),
                          type=c(rep("ne_TF",length(unique(ne_inter2$TF.Name))),rep("ne_TG",length(unique(ne_inter2$TG.Name)))),
                          cor=abs(as.vector(fto_cor_dat$cor[match(c(unique(ne_inter2$TF.Name),unique(ne_inter2$TG.Name)),fto_cor_dat$inter_gene)])))
# write.table(ne_inter2,"ne_inter.txt",sep = "\t",quote = F,row.names = F)
# write.table(ne_inter_anno,"ne_inter_anno.txt",sep = "\t",quote = F,row.names = F)
all_inter<-rbind(po_inter2,ne_inter2)
all_inter_anno<-rbind(po_inter_anno,ne_inter_anno)
write.table(all_inter,paste(names(tf)[k],"all_inter.txt",sep = "_"),sep = "\t",quote = F,row.names = F)
write.table(all_inter_anno,paste(names(tf)[k],"all_inter_anno.txt",sep = "_"),sep = "\t",quote = F,row.names = F)
##################绘制boxplot
library(reshape2)
library(RColorBrewer)
library(ggsignif)
color<-brewer.pal(n = 11, "RdYlGn")
tf_exp<-norm_exp[all_inter_anno$symbol,]
table(colnames(tf_exp)==anno_hc_sci$sample)
colnames(tf_exp)<-anno_hc_sci$type
tf_exp<-as.data.frame(tf_exp)
aa<-unique(all_inter$TF.Name)
tf_exp<-tf_exp[aa,]

all_p<-list()
for (i in 1:dim(tf_exp)[1]) {
  plot_dat<-melt(as.matrix(tf_exp[i,]))
  p<-ggplot(plot_dat, aes(x=Var2, y=value,color=Var2))+
    stat_boxplot(geom = 'errorbar',width=0.15)+
    geom_boxplot()+
    geom_jitter(width = 0.2, alpha = 0.5, size=1)+
    geom_signif(comparisons = list(c("Control","Spinal_Cord_Injury")),step_increase=0.05,map_signif_level = F,test = wilcox.test,color="black")+
    theme_bw()+
    scale_color_manual(values = color[c(2,10)])+
    labs(x="Sample type",y= paste("The expression of ",row.names(tf_exp)[i],sep = ""))+
    theme(legend.position = "none")
  all_p<-c(all_p,list(p))
}
names(all_p)<-row.names(tf_exp)
save(all_p,file = paste(names(tf)[k],"gene_boxplot.rdata",sep = "_"))
library(cowplot)
plot_grid(all_p[[1]],all_p[[2]],all_p[[3]],all_p[[4]],ncol = 2)
###########################GO
color<-brewer.pal(n = 11, "Spectral")
library(clusterProfiler)
library(org.Hs.eg.db)
library(GOplot)
library(enrichplot)
inter_anno<-all_inter_anno
go_cluster<-enrichGO(inter_anno$symbol,OrgDb = org.Hs.eg.db,
                     keyType = "SYMBOL",
                     ont = "BP",
                     pvalueCutoff = 0.05,
                     qvalueCutoff = 0.05,
                     pAdjustMethod = "BH")
go_dat<-go_cluster@result[1:10,]
go_dat$Description<-factor(go_dat$Description,levels = go_dat$Description[10:1])
p1<-ggplot(go_dat,aes(x=Description,y=Count,fill= -log10(pvalue)))+
  geom_bar(stat="identity")+
  coord_flip()+
  theme_bw()+
  scale_fill_gradientn(colours =rev(color))

eg <- bitr(inter_anno$symbol,fromType="SYMBOL", toType="ENTREZID", OrgDb="org.Hs.eg.db")
kegg_ne <- enrichKEGG(eg$ENTREZID, organism = 'hsa', keyType = 'ncbi-geneid', pvalueCutoff = 0.05,pAdjustMethod = 'BH', minGSSize = 10,maxGSSize = 500,qvalueCutoff = 0.2,use_internal_data = FALSE)
b<-kegg_ne@result
e<-b[1:10,]
e$Description<-factor(e$Description,levels = e$Description[10:1])
p2<-ggplot(e,aes(x=Description,y=Count,fill= -log10(pvalue)))+
  geom_bar(stat="identity")+
  coord_flip()+
  theme_bw()+
  scale_fill_gradientn(colours =rev(color))
library(cowplot)
plot_grid(p1,p2,ncol = 1)
