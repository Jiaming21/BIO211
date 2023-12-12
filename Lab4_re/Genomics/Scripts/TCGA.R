chooseCRANmirror()

##########################  Prepare R Packages  ################################
install.packages('stringr')
install.packages('ggplotify')
install.packages('patchwork')
install.packages('cowplot')

# install.packages("BiocManager")
BiocManager::install(c("edgeR","limma")) 

################  Set work path, read data and data statistics   ###############
path <- '/Users/jiaming/Desktop/BIO211_pre/Lab4_re/Genomics/Datas'
setwd(path)

rm(list=ls())
a=read.csv('LIHC_counts_lab4.csv')
dim(a) # 60488 425 ｜ 行是每个基因表达量，列是样本名（是否是normal还是tumor样本）

row.names(a) <- a[,1] # 选出第一列作为行名
a <- a[,-1] # 去除第一列
group_list <- ifelse(as.numeric(substr(colnames(a),14,15)) < 10,'tumor','normal') # 分为两组
# 1-9 normal; 10-19 normal; 20-29 control

group_list <- factor(group_list, levels=c("normal","tumor"))
table(group_list)

View(a) # Check the data matrix

###################       选择出差异表达的显著基因      ########################
## Method 1: edgeR
library(edgeR)
dge<-DGEList(counts=a, group=group_list)

# 查看DGEList的属性
attributes(dge)
attributes(dge$samples)
dge$samples$lib.size<-colSums(dge$counts) # 为samples加上列的总计数信息
# dge$samples$lib.size有425个数，对应425个样本
design<-model.matrix(~0+group_list) # 425行，两列分别对应normal与tumor
# 相当于就是给每组是tumor还是normal做编码
# 创建没有截距的设计矩阵 
# 其中每一列对应于样本分组（group_list）的水平。
# 这个设计矩阵将用于线性模型，以描述基因表达在不同分组之间的变化。
# 在差异表达分析中，通过比较不同分组的系数，可以识别基因的差异表达。
rownames(design)<-colnames(dge) # 将TCGA编号给矩阵的行名
colnames(design)<-levels(group_list) # 将factor的"normal"和"tumor"给design的列名
dge<-estimateGLMCommonDisp(dge,design) # 估算基因表达数据的离散度（dispersion）
dge<-estimateGLMTrendedDisp(dge,design) # 估算每个基因的趋势离散度（trended dispersion）
dge<-estimateGLMTagwiseDisp(dge,design) # 估算基因表达数据的标记依赖离散度（tagwise dispersion）

fit<-glmFit(dge,design) # 拟合一个广义线性模型（GLM）到基因表达数据
# fit 对象通常包含以下信息：
# 估计的模型系数：表示不同样本组之间的平均表达水平差异。
# 估计的标准差：表示每个基因表达值的方差或离散度。
# 模型拟合的统计信息：例如对数似然值等。
# fit 对象通常会用于进行后续的差异表达分析，如似然比检验、负二项分布回归等。

fit2<-glmLRT(fit,contrast=c(-1,1)) # 进行广义线性模型（GLM）的似然比检验（Likelihood Ratio Test，LRT）
# contrast=c(-1,1) 设置是指定似然比检验中的对比（contrast）系数
DEG2=topTags(fit2, n=nrow(a))
# fit2 是先前通过 glmLRT 进行似然比检验得到的结果对象，其中包含了每个基因的统计信息。
# topTags 是 edgeR 包中的函数，用于提取在检验中具有最高统计显著性的基因。
# nrow(a) 给定了你想要提取的基因的数量，通常是数据集中所有的基因的数量。
# 这里使用 nrow(a) 是为了提取所有的基因。
DEG2=as.data.frame(DEG2) # 只用logFC属性
# 结果对象 DEG2 包含了似然比检验中具有最高统计显著性的基因的信息，
# 包括基因的名称
# 对应的 p-value
# 调整过的 p-value（如果进行了多重检验校正的话）
# 估计的对数折算比（logFC）
logFC_cutoff2<-with(DEG2, mean(abs(logFC))+2*sd(abs(logFC)))
# 计算一个对数折算比（logFC）的阈值
# 用于确定哪些基因在似然比检验中被认为是显著差异表达的
DEG2$change=as.factor(
  ifelse(DEG2$PValue<0.05 & abs(DEG2$logFC)>logFC_cutoff2,
  ifelse(DEG2$logFC>logFC_cutoff2,"UP","DOWN"),"NOT")
  )
# 如果满足DEG2$PValue<0.05 & abs(DEG2$logFC)>logFC_cutoff2 进行后续判断是否是"UP"还是
# "DOWN",不满足直接"NOT"

# 为 DEG2 数据框添加一个名为 change 的列，该列用于表示每个基因的差异表达情况
head(DEG2)

## Method 2: limma-voom
library(limma)
design<-model.matrix(~0+group_list) # 搞一个基准为0的矩阵，对tumor还是normal进行编码，一共425个样本（是normal/tumor）
colnames(design)=levels(group_list) # 列名注释为normal还是tumor
rownames(design)=colnames(a) # 给矩阵的行以样本名
dge<-DGEList(counts=a) # 创建一个DGEList给dge
dge<-calcNormFactors(dge) 
# 计算规范化因子（normalization factors）的函数
# 在差异表达分析中，规范化是为了消除不同样本之间的技术偏差和批次效应，
# 以便更准确地比较基因表达水平
logCPM<-cpm(dge,log=TRUE,prior.count=3) # 计算每百万对数
# log counts per million，logCPM
# dge 是一个 DGEList 对象，包含原始的基因计数数据以及实验设计和样本信息
# cpm 是 edgeR 包中的函数，用于计算每百万对数计数。
# log = TRUE 参数指定返回的值是对数变换后的，即 logCPM
# prior.count = 3 参数指定加入一个小的平滑值，以避免在计算对数时遇到零值
# 这是为了稳定对数变换，尤其是在计算对数时避免零计数导致无穷大或负无穷大的问题
v<-voom(dge,design,normalize="quantile")
# 使用 voom 函数进行基因表达数据的转换，将数据转换为适用于线性模型的形式
fit<-lmFit(v,design)
# 使用 lmFit 函数对 voom 转换后的数据进行线性模型的拟合
constra=paste(rev(levels(group_list)),collapse="-")
cont.matrix<-makeContrasts(contrasts=constra,levels=design)
# 创建一个对比矩阵（contrast matrix）用于在 limma 包中进行线性模型的分析
fit3=contrasts.fit(fit,cont.matrix)
# fit3 是包含了应用对比矩阵后的线性模型的结果对象，用于进行后续的统计分析
fit3=eBayes(fit3) # 对已经应用了对比矩阵的线性模型进行贝叶斯估计
DEG3=topTable(fit3,coef=constra,n=Inf) # 从贝叶斯估计后的线性模型结果中提取差异表达基因（DEG）
DEG3=na.omit(DEG3) # 使用 na.omit 函数去除结果对象 DEG3 中包含缺失值（NA，Not Available）的行
logFC_cutoff3<-with(DEG3,mean(abs(logFC)+2*sd(abs(logFC))))
DEG3$change=as.factor(
  ifelse(DEG3$P.Value<0.05 & abs(DEG3$logFC) > logFC_cutoff3,
  ifelse(DEG3$logFC > logFC_cutoff3,"UP","DOWN"),"NOT")
  )
head(DEG3)
# Count the number of up- and down-regulated genes
table(DEG3$change)
# Save your results
edgeR_DEG<-DEG2
limma_voom_DEG<-DEG3
save(edgeR_DEG,limma_voom_DEG,group_list,file='DEG.Rdata')

####################  Data Visualization by Volcano Plot   #####################

# P 值的负对数（通常以10为底）绘制在 y 轴上（通常以 10 为基数）
# P 值低（高度显著）的数据点就会出现在图的顶部
# 出现在图的顶端。x 轴是两种情况下折叠变化的对数

# x 轴是两种情况下折叠变化的对数
# 使用折减变化的对数是为了使两个方向的变化看起来与中心等距
# 以这种方式绘制点会在图中产生两个感兴趣的区域
# 这些点位于图的顶端，远离左侧或右侧。
# 这些点显示出较大的折叠变化幅度（因此位于中心的左侧或右侧）以及较高的统计显著性（因此位于顶部）

# 1. Data preparation
library(stringr)
rownames(DEG2)<-str_sub(rownames(DEG2),start=1,end=15)
DEG2$ENSEMBL<-rownames(DEG2)
diff_gene<-DEG2[DEG2$change!='NOT',]
write.csv(diff_gene,
            file='/Users/jiaming/Desktop/BIO211_pre/Lab4_re/Genomics/Datas/LIHC_diff_gene.csv',
            sep='\t',
            quote=F)

# 2.  Construct a volcano plot
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
p<-qplot(x=DEG2$logFC,
         y=-log10(DEG2$FDR),
         xlab='log2(FC)',
         ylab='-log10(FDR)',
         size=I(0.7),colour=significant)
p<-p+scale_color_manual(values=c('up'='red','normal'='grey','down'='blue'))
## Add cut-off lines & export image
xline=c(-log2(2),log2(2))
p<-p+geom_vline(xintercept = xline,lty=2,size=I(0.2),color='grey11')
yline=-log(0.05,10)
p<-p+geom_hline(yintercept = yline,lty=2,size=I(0.2),color='grey11')
p<-p+theme_bw()+theme(panel.background = element_rect(colour = 'black',size=1,fill='white'),panel.grid=element_blank())
pdf('deg2_volcano.pdf')
print(p)
dev.off()
############    Data Organization by Hierarchical Clustering    ################

# 在对基因表达数据进行聚类分析时，有各种各样的算法数据进行聚类分析的算法多种多样。
# 分层聚类是最常用的方法之一。
# 我们将在本节中对其进行简要研究。
# 每当绘制热图时，heatmap() 方法实际上都会对表达式数据应用默认的分层聚类方法。
# 使用欧氏距离和完全关联。
# 这并不是 BioConductor 中聚类分析的理想方案，因为已有几种复杂的方法（例如：agnes(), diana(), hopach()），
# 但为了举例和简单起见，我们将坚持使用热图实现。
# 让我们通过基因聚类和样本聚类来研究热图的聚类功能。

# install.packages('pheatmap')
cg1=rownames(edgeR_DEG)[edgeR_DEG$change!='NOT']
cg2=rownames(limma_voom_DEG)[limma_voom_DEG$change!='NOT']
library(pheatmap)
library(RColorBrewer)
color<-colorRampPalette(c('#436EEE','white','#EE0000'))(100)

## edgeR: Hierarchical clustering produced significant genes
mat1=a[cg1,]
n1=t(scale(t(mat1)))
n1[n1>1]=1
n1[n1<-1]=-1
ac=data.frame(group=group_list)
rownames(ac)=colnames(mat1)
ht1<-pheatmap(n1,show_rownames = F,show_colnames = F,
                cluster_rows = F,cluster_cols = T,
                annotation_col = ac, color = color)
print(ht1)
# dev.off()

## limma-voom: Hierarchical clustering produced significant genes
mat2=a[cg2,]
n2=t(scale(t(mat2)))
n2[n2>1]=1
n2[n2<-1]=-1
ac=data.frame(group=group_list)
rownames(ac)=colnames(mat2)
ht2=pheatmap(n2,show_rownames = F, 
             show_colnames = F, 
             cluster_rows = F, 
             cluster_cols = T, 
             annotation_col = ac, 
             color = color)
print(ht2)
dev.off()

###################     KEGG pathway enrichment results.    ####################
# Read enrichment data
pathway_enrich<-read.table("enrichment_results_wg_kegg.txt",header=T,sep="\t",stringsAsFactors=F)
# Draw a scatter plot
library(ggplot2)
colnames(pathway_enrich)
p<-ggplot(pathway_enrich,aes(x=enrichmentRatio,y=description))
p+geom_point()
# Map the size of the point to overlap and map the color to FDR

p+geom_point(aes(size=overlap,color=FDR))+
  scale_color_gradient(low="red",high="green")+
  labs(title="Statistics of Pathway Enrichment",x="Enrichment Ratio",y="",color="gvalue",size="Gene-number")+
  theme_bw()

ggsave('kegg.png',width=7,height=7)

dev.off()

########################     GO enrichment results    ##########################     
# Read enrichment data
go_enrich<-read.table("enrichment_results_wg_go.txt",header=T,sep="\t",stringsAsFactors=F)
# Draw a scatter plot
library(ggplot2)
colnames(go_enrich)
p<-ggplot(go_enrich,aes(x=enrichmentRatio,y=description))
p+geom_point()
# Map the size of the point to overlap and map the color to FDR

p+geom_point(aes(size=overlap,color=FDR))+
  scale_color_gradient(low="red",high="green")+
  labs(title="Statistics of GO Enrichment",x="Enrichment Ratio",y="",color="gvalue",size="Gene-number")+
  theme_bw()

ggsave('go.png',width=9,height=7)

dev.off()







