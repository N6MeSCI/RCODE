library(Boruta)
par(mfrow=c(2,1))
load("hc_sci_anno_exp.rdata")
load("m6a_cor_dat.rdata")
lanno<- as.data.frame(rtracklayer::import("gencode.v29.long_noncoding_RNAs.gtf"))
lncrna_name<-unique(lanno$gene_name)
m6a_cor_lncRNA<-all_cor_dat[all_cor_dat$inter_gene%in%lncrna_name,]
m6a_lnc<-unique(m6a_cor_lncRNA$inter_gene)

lnc_exp<-norm_exp[m6a_lnc,]
lnc_exp2<-as.data.frame(t(lnc_exp))
table(row.names(lnc_exp2)==row.names(anno_hc_sci))
lnc_exp2$class<-anno_hc_sci$type
save(lnc_exp,anno_hc_sci,file = "lncRNA_exp_anno_dat.rdata")
###############Boruta
gse10_boruta<-Boruta(class~.,data=lnc_exp2,doTrace = 2)
print(gse10_boruta)
plot(gse10_boruta,xlab = "", xaxt = "n")
lz<-lapply(1:ncol(gse10_boruta$ImpHistory),function(i)
  
  gse10_boruta$ImpHistory[is.finite(gse10_boruta$ImpHistory[,i]),i])
names(lz) <- colnames(gse10_boruta$ImpHistory)
Labels <- sort(sapply(lz,median))
axis(side = 1,las=2,labels = names(Labels),
     at = 1:ncol(gse10_boruta$ImpHistory), cex.axis = 0.5)
final_gse10_boruta <- TentativeRoughFix(gse10_boruta)
print(final_gse10_boruta)
plot(final_gse10_boruta,xlab = "", xaxt = "n")
axis(side = 1,las=2,labels = names(Labels),
     at = 1:ncol(final_gse10_boruta$ImpHistory), cex.axis = 0.5)

confirm_result<-getSelectedAttributes(final_gse10_boruta, withTentative = F)
save(confirm_result,file = "confirm_result.rdata")
