library(clusterProfiler)
library(org.Hs.eg.db)
library(GOplot)
library(enrichplot)
library(fgsea)
library(stringr)
library(ggsci)
color<-brewer.pal(n = 11, "RdYlGn")
load("selectgenes.rdata")
load("hc_sci_anno_exp.rdata")
m6a_exp<-as.data.frame(norm_exp[row.names(selectgenes),])
table(row.names(anno_hc_sci)==colnames(m6a_exp))


person_cor<-apply(m6a_exp,1,function(x){
  aa=apply(norm_exp,1,function(y){
    dd<-cor.test(x,y)
    list(dd[["p.value"]],dd[["estimate"]][["cor"]])
  })
  aa[unlist(lapply(aa,function(x){
    x[[1]]<0.01
  }))]
})

save(person_cor,file = "person_cor.rdata")

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
save(all_cor_dat,file = "m6a_cor_dat.rdata")

all_cor_dat2<-all_cor_dat[,c(3,1,2)]

inter_dat<-all_cor_dat2[,1:2]
inter_anno<-data.frame(symbol=c(unique(all_cor_dat2$select_gene),unique(all_cor_dat2$inter_gene)),
                       type=c(rep("select",length(unique(all_cor_dat2$select_gene))),
                              rep("interact",length(unique(all_cor_dat2$inter_gene)))))
write.table(inter_dat,"inter_dat.txt",sep = "\t",row.names = F,quote = F)
write.table(inter_anno,"inter_anno.txt",sep = "\t",row.names = F,quote = F)
###########################GO
dat<-inter_anno
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
