# Set work path, read data and data statistics
path<-'/Users/jiaming/Desktop/BIO211_pre/Lab4/gene_expression_data_analysis/datas'
setwd(paste(path))
rm(list=ls()) # 清空env
a=read.csv('LIHC_counts_lab4.csv')
dim(a)
row.names(a)<-a[,1]
a<-a[,-1]
group_list<-ifelse(as.numeric(substr(colnames(a),14,15))<10,'tumor','normal')
# Separate samples into two groups
# Check the hints on TCGA Barcode for how to differentiate tumor/normal samples
group_list<-factor(group_list,levels=c("normal","tumor"))
table(group_list)
View(a)
# Check the data matrix

# Select significant genes exhibiting differential expression
## Method 1: edgeR
library(edgeR)
dge<-DGEList(counts=a, group=group_list)
dge$samples$lib.size<-colSums(dge$counts)
design<-model.matrix(~0+group_list)
rownames(design)<-colnames(dge)
colnames(design)<-levels(group_list)
dge<-estimateGLMCommonDisp(dge,design)
dge<-estimateGLMTrendedDisp(dge,design)
dge<-estimateGLMTagwiseDisp(dge,design)
fit<-glmFit(dge,design)
fit2<-glmLRT(fit,contrast=c(-1,1))
DEG2=topTags(fit2, n=nrow(a))
DEG2=as.data.frame(DEG2)
logFC_cutoff2<-with(DEG2, mean(abs(logFC))+2*sd(abs(logFC)))
DEG2$change=as.factor(ifelse(DEG2$PValue<0.05 & abs(DEG2$logFC) > logFC_cutoff2, ifelse(DEG2$logFC > logFC_cutoff2,"UP","DOWN"),"NOT"))
head(DEG2)
table(DEG2$change)
## Method 2: limma-voom
library(limma)
design<-model.matrix(~0+group_list)
colnames(design)=levels(group_list)
rownames(design)=colnames(a)
dge<-DGEList(counts=a)
dge<-calcNormFactors(dge)
logCPM<-cpm(dge,log=TRUE,prior.count=3)
v<-voom(dge,design,normalize="quantile")
fit<-lmFit(v,design)
constra=paste(rev(levels(group_list)),collapse="-")
cont.matrix<-makeContrasts(contrasts=constra,levels=design)
fit3=contrasts.fit(fit,cont.matrix)
fit3=eBayes(fit3)
DEG3=topTable(fit3,coef=constra,n=Inf)
DEG3=na.omit(DEG3)
logFC_cutoff3<-with(DEG3,mean(abs(logFC)+2*sd(abs(logFC))))
DEG3$change=as.factor(
  ifelse(DEG3$P.Value<0.05 & abs(DEG3$logFC) > logFC_cutoff3,         
         ifelse(DEG3$logFC > logFC_cutoff3,"UP","DOWN"),"NOT")
  )
head(DEG3)
table(DEG3$change)
# Count the number of up- and down-regulated genes

## Save your results
edgeR_DEG<-DEG2
limma_voom_DEG<-DEG3
save(edgeR_DEG,limma_voom_DEG,group_list,file='DEG.Rdata')

##################   Data Visualization by Volcano Plot   ######################
# Data preparation
library(stringr)
rownames(DEG2)<-str_sub(rownames(DEG2),start=1,end=15)
DEG2$ENSEMBL<-rownames(DEG2)
diff_gene<-DEG2[DEG2$change!='NOT',]
write.table(diff_gene,
            file=path,
            sep='\t',quote=F)
# Construct a volcano plot
# install.packages('ggplot2')
library(cowplot)
library(patchwork)
library(ggplotify)
library(ggplot2)
loc_up<-intersect(which(DEG2$FDR<0.05),which(DEG2$logFC>=1))
loc_down<-intersect(which(DEG2$FDR<0.05),which(DEG2$logFC<(-1)))
significant<-rep('normal',times=nrow(DEG2))
significant[loc_up]<-'up'
significant[loc_down]<-'down'
significant<-factor(significant,levels=c('up','down','normal'))
p<-qplot(x=DEG2$logFC,y=-log10(DEG2$FDR),xlab='log2(FC)',ylab='-log10(FDR)',         + size=I(0.7),colour=significant)
p<-p+scale_color_manual(values=c('up'='red','normal'='grey','down'='blue'))
## Add cut-off lines & export image
xline=c(-log2(2),log2(2))
p<-p+geom_vline(xintercept = xline,lty=2,size=I(0.2),color='grey11')
yline=-log(0.05,10)
p<-p+geom_hline(yintercept = yline,lty=2,size=I(0.2),color='grey11')
p<-p+theme_bw()+theme(panel.background = element_rect(colour = 'black',                                                      + size=1,fill='white'),panel.grid=element_blank())
pdf('deg2_volcano.pdf')
print(p)
dev.off()



