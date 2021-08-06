library(ggplot2)
load("hc_sci_anno_exp.rdata")
load("m6a_cor_dat.rdata")
lanno<- as.data.frame(rtracklayer::import("gencode.v29.long_noncoding_RNAs.gtf"))
lncrna_name<-unique(lanno$gene_name)
m6a_cor_lncRNA<-all_cor_dat[all_cor_dat$inter_gene%in%lncrna_name,]
length(unique(m6a_cor_lncRNA$inter_gene))

cerna_dat<-read.csv("lncRNA_ceRNA.txt",header = T,sep = "\t")
cor_cerna<-cerna_dat[cerna_dat$geneName%in%unique(m6a_cor_lncRNA$inter_gene),]
cor_lnc<-unique(cor_cerna$geneName)

lnc_exp<-norm_exp[cor_lnc,]
person_cor<-apply(lnc_exp,1,function(x){
  aa=apply(norm_exp,1,function(y){
    dd<-cor.test(x,y)
    list(dd[["p.value"]],dd[["estimate"]][["cor"]])
  })
  aa[unlist(lapply(aa,function(x){
    x[[1]]<0.01
  }))]
})
save(person_cor,file = "person_cor_lncRNA.rdata")

aa<-person_cor
all_cor_lnc<-vector()
for (i in 1:length(aa)) {
  dat<-aa[[i]]
  if(length(dat)!=0){
    cor_dat<-data.frame(inter_gene=names(dat),cor=unlist(lapply(dat,function(x){x[[2]]})),row.names = NULL)
    cor_dat$select_gene<-names(aa)[i]
    cor_dat<-cor_dat[cor_dat$inter_gene%in%cerna_dat$ceRNAname[cerna_dat$geneName%in%names(aa)[i]],]
    cor_dat$mirna<-cerna_dat$hitMiRNAFamily[cerna_dat$geneName%in%names(aa)[i] & cerna_dat$ceRNAname%in%unique(cor_dat$inter_gene)]
    all_cor_lnc<-rbind(all_cor_lnc,cor_dat)
  }
}
all_cor_lnc<-all_cor_lnc[all_cor_lnc$cor>0.7 & all_cor_lnc$cor<1,]
all_cor_lnc<-all_cor_lnc[all_cor_lnc$inter_gene!=all_cor_lnc$select_gene,]
save(all_cor_lnc,file = "lnc_cor_dat.rdata")
all_cor_lnc2<-all_cor_lnc[,c(3,1,2,4)]
inter_dat<-all_cor_lnc2[,1:2]
inter_anno<-data.frame(symbol=c(unique(all_cor_lnc2$select_gene),unique(all_cor_lnc2$inter_gene)),
                       type=c(rep("lncRNA",length(unique(all_cor_lnc2$select_gene))),
                              rep("mRNA",length(unique(all_cor_lnc2$inter_gene)))))
write.table(inter_dat,"inter_ceRNA_dat.txt",sep = "\t",row.names = F,quote = F)
write.table(inter_anno,"inter_ceRNA_anno.txt",sep = "\t",row.names = F,quote = F)

stat<-data.frame(type=c("mRNA","lncRNA"),
                 num=c(length(unique(all_cor_lnc2$inter_gene)),length(unique(all_cor_lnc2$select_gene))))
ggplot(stat,aes(x=type,y=num,fill=type))+
  geom_bar(stat = "identity", width = 0.7, position = "identity")+
  scale_fill_manual(values = color[c(2,8)])+
  theme_bw()+
  labs(x="",y="The number of genes")+
  theme(legend.position = "none")

################################functional enrichment analysis
library(clusterProfiler)
library(org.Hs.eg.db)
library(GOplot)
library(enrichplot)
library(fgsea)
library(stringr)
library(ggsci)
library(RColorBrewer)
color<-brewer.pal(n = 11, "Spectral")
dat<-all_cor_lnc2
ego_BP <- enrichGO(gene = unique(as.character(dat$inter_gene)),
                   keyType = "SYMBOL",
                   OrgDb= org.Hs.eg.db,
                   ont = "ALL",
                   pAdjustMethod = "BH",
                   minGSSize = 1,
                   pvalueCutoff = 0.05,
                   qvalueCutoff = 0.05)
ego <- as.data.frame(ego_BP@result)
colnames(ego)
ego <- ego[,c(2,3,9,7)] 

ego$geneID <- str_replace_all(ego$geneID,"/",",") 
names(ego)=c("ID","Term","Genes","adj_pval")
ego$Category = "BP"
head(ego)
genes = data.frame(ID=as.character(dat$Target.Gene),
                   logFC=rnorm(length(as.character(dat$Target.Gene))))
head(genes)

color<-RColorBrewer::brewer.pal(n = 11, name = "Spectral")
circ <- circle_dat(ego,genes)
chord <- chord_dat(data=circ, genes=genes,process = ego$Term) # 
GOChord(chord, space = 0.02, gene.order = 'logFC', gene.space = 0.25, gene.size = 3,lfc.min= -5,lfc.max=5,process.label=2,ribbon.col=color[c(1,3,5,7,9)],lfc.col = c("#4DFDFD","black","#FAD103"),
        border.size = 0.05)
#############################KEGG
entr_id<-bitr(unique(as.character(dat$inter_gene)),fromType = "SYMBOL",toType = "ENTREZID",OrgDb = "org.Hs.eg.db")
ego_kegg <- enrichKEGG(gene = as.character(entr_id$ENTREZID),
                       keyType = "kegg",
                       organism = "hsa",
                       pAdjustMethod = "BH",
                       minGSSize = 1,
                       pvalueCutoff = 0.05,
                       qvalueCutoff = 0.05)


ego <- as.data.frame(ego_kegg@result)
ego<-ego[ego$pvalue<0.05,]
colnames(ego)
aa<-c(grep("pathway",ego$Description))
kegg_dot<-ego[aa[1:20],]
ggplot(kegg_dot,aes(x=pvalue,y=Description))+
  geom_point(aes(size=Count,color=pvalue))+  
  scale_color_gradientn(colours = color)+
  labs(color=expression(Pvalue),size="Count",  
       x="Pvalue",y="Pathway",title="KEGG Pathway enrichment")+
  theme_bw()

ego <- ego[aa[1:10],c(1,2,8,6)] 
ego$geneID <- str_replace_all(ego$geneID,"/",",") 
names(ego)=c("ID","Term","Genes","adj_pval")
ego$Category = "KEGG"
for (i in 1:dim(ego)[1]) {
  entr<-ego$Genes[i]
  entr<-unlist(strsplit(entr,split = ","))
  sym_name<-entr_id$SYMBOL[entr_id$ENTREZID%in%entr]
  ego$Genes[i]<-paste0(c(sym_name),collapse = ",")
}

head(ego)
genes = data.frame(ID=unique(as.character(dat$inter_gene)),
                   logFC=rnorm(length(unique(as.character(dat$inter_gene)))))
head(genes)

circ <- circle_dat(ego,genes)
chord <- chord_dat(data=circ, genes=genes,process = ego$Term) # 

GOChord(chord, space = 0.02, gene.order = 'logFC', gene.space = 0.25, gene.size = 4,lfc.min= -5,lfc.max=5,process.label=2,ribbon.col=color[2:11],lfc.col = c("#4DFDFD","black","#FAD103"),
        border.size = 0.05)
