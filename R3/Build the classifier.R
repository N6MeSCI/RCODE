library(pROC)
library(e1071)
library(RColorBrewer)
library(pheatmap)
par(mfrow=c(1,2))
col<-brewer.pal(11,"Spectral")
load("lncRNA_exp_anno_dat.rdata")
load("confirm_result.rdata")
confirm_result<-gsub("`","",confirm_result)

m6a_exp<-lnc_exp[confirm_result,]
anno_col<-data.frame(Sample_type=as.character(anno_hc_sci$type))
row.names(anno_col)<-colnames(m6a_exp)
anno_color<-list(Sample_type=col[c(1,9)])
names(anno_color[[1]])<-unique(as.character(anno_hc_sci$type))
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

feature_exp<-as.data.frame(t(lnc_exp[confirm_result,]))
table(row.names(feature_exp)==row.names(anno_hc_sci))
feature_exp$class<-anno_hc_sci$type

train_sub = sample(nrow(feature_exp),5/10*nrow(feature_exp))
save(train_sub,file = "train_site.rdata")
train_data = feature_exp[train_sub,]
test_data = feature_exp[-train_sub,]
#######################SVM
ciber10_svm<- svm(as.numeric(class)-1 ~ ., 
                  data = train_data,
                  type = 'eps',kernel = 'radial')
summary(ciber10_svm)
ciber10_pre_svm <- predict(ciber10_svm,newdata = test_data)
svm_roc <- roc(test_data$class,ciber10_pre_svm)
plot(svm_roc, print.auc=TRUE, auc.polygon=TRUE, grid=c(0.1, 0.2),grid.col=c(col[2], col[9]), max.auc.polygon=TRUE,auc.polygon.col=col[5], print.thres=TRUE,main='SVM ROC_curve')
######################Decision tree
library("rpart")
library("rpart.plot")
library(survival)
model_tree <- rpart(Surv(as.numeric(class))~.,data = train_data, method = 'exp')
rpart.plot(model_tree)
test_tree = predict(model_tree, test_data) 
# table(test_tree,test_data$class)
test_roc <- roc(test_data$class,test_tree)
plot(test_roc, print.auc=TRUE, auc.polygon=TRUE, grid=c(0.1, 0.2),grid.col=c(col[2], col[9]), max.auc.polygon=TRUE,auc.polygon.col=col[5], print.thres=TRUE,main='Decision tree ROC_curve')
