library(BiocManager)
library(clusterProfiler)
library(org.Dr.eg.db)
datas <- read.csv('/Users/Jiaming/Desktop/BIO211_pre/ZebraFish/GO_KEGG/Datas/bioDBnet.csv')
genes <- datas$Ensembl.Gene.ID

# GO富集分析
ego1 <- enrichGO(gene         = genes,
                OrgDb        = org.Dr.eg.db,
                keyType      = "ENSEMBL", # "SYMBOL"表示基因的标准名称或符号
                ont          = "BP", # 选择BP, CC或MF
                pAdjustMethod = "BH", # 使用Benjamini-Hochberg方法调整p值
                qvalueCutoff = 0.05, # 设置显著性阈值为0.05
                readable     = TRUE) # 将结果转换为可读的基因名
ego2 <- enrichGO(gene         = genes,
                 OrgDb        = org.Dr.eg.db,
                 keyType      = "ENSEMBL", # "SYMBOL"表示基因的标准名称或符号
                 ont          = "CC", # 选择BP, CC或MF
                 pAdjustMethod = "BH", # 使用Benjamini-Hochberg方法调整p值
                 qvalueCutoff = 0.05, # 设置显著性阈值为0.05
                 readable     = TRUE) # 将结果转换为可读的基因名
ego3 <- enrichGO(gene         = genes,
                 OrgDb        = org.Dr.eg.db,
                 keyType      = "ENSEMBL", # "SYMBOL"表示基因的标准名称或符号
                 ont          = "MF", # 选择BP, CC或MF
                 pAdjustMethod = "BH", # 使用Benjamini-Hochberg方法调整p值
                 qvalueCutoff = 0.05, # 设置显著性阈值为0.05
                 readable     = TRUE) # 将结果转换为可读的基因名

GOplotIn_BP <- ego1[1:8,c(1,2,6,8)] # 提取GO富集BP的前29行,提取ID,Description,p.adjust,GeneID四列
GOplotIn_CC <- ego2[1:3,c(1,2,6,8)] # 提取GO富集CC的前29行,提取ID,Description,p.adjust,GeneID四列
GOplotIn_MF <- ego3[1:8,c(1,2,6,8)] # 提取GO富集MF的前29行,提取ID,Description,p.adjust,GeneID四列

GOplotIn_BP$geneID <-str_replace_all(GOplotIn_BP$geneID,'/',',') # 把GeneID列中的’/’替换成‘,’
GOplotIn_CC$geneID <-str_replace_all(GOplotIn_CC$geneID,'/',',')
GOplotIn_MF$geneID <-str_replace_all(GOplotIn_MF$geneID,'/',',')

names(GOplotIn_BP)<-c('ID','Term','adj_pval','Genes') # 修改列名,后面弦图绘制的时候需要这样的格式
names(GOplotIn_CC)<-c('ID','Term','adj_pval','Genes')
names(GOplotIn_MF)<-c('ID','Term','adj_pval','Genes')

GOplotIn_BP$Category = "BP" # 分类信息
GOplotIn_CC$Category = "CC"
GOplotIn_MF$Category = "MF"



GOplotIn_BP <- data.frame(
  ID = c("GO:0008150", "GO:0009987", "GO:0008152"),
  Term = c("biological_process", "cellular process", "metabolic process"),
  adj_pval = c(0.01, 0.05, 0.02),
  Genes = c("Gene1,Gene2", "Gene3,Gene4", "Gene5,Gene6")
)

genes <- data.frame(
  ID = c("Gene1", "Gene2", "Gene3", "Gene4", "Gene5", "Gene6"),
  logFC = c(1.5, -2.0, 0.5, -1.2, 1.8, -2.2)
)



# genes <- data.frame(ID = names(gene), logFC = gene)
class(genes)
circ_BP <- GOplot::circle_dat(GOplotIn_BP, genes) # GOplot导入数据格式整理
circ_CC <- GOplot::circle_dat(GOplotIn_CC, genes) 
circ_MF <- GOplot::circle_dat(GOplotIn_MF, genes) 

chord_BP <- chord_dat(data = circ_BP, genes = genes) #生成含有选定基因的数据框
chord_CC <- chord_dat(data = circ_CC, genes = genes) 
chord_MF <- chord_dat(data = circ_MF, genes = genes) 
s
GOChord(data = chord_BP,#弦图
        title = 'GO-Biological Process',space = 0.01,#GO Term间距
        limit = c(1,1),gene.order = 'logFC',gene.space = 0.25,gene.size = 5,
        lfc.col = c('red','white','blue'), #上下调基因颜色
        process.label = 10) #GO Term字体大小
GOChord(data = chord_CC,title = 'GO-Cellular Component',space = 0.01,
        limit = c(1,1),gene.order = 'logFC',gene.space = 0.25,gene.size = 5,
        lfc.col = c('red','white','blue'), 
        process.label = 10) 
GOChord(data = chord_MF,title = 'GO-Mollecular Function',space = 0.01,
        limit = c(1,1),gene.order = 'logFC',gene.space = 0.25,gene.size = 5,
        lfc.col = c('red','white','blue'), 
        process.label = 10)


# 假设您有一些基因表达数据，这里是一个示例向量
trans_id <- read.csv('/Users/Jiaming/Desktop/BIO211_pre/ZebraFish/GO_KEGG/Datas/trans_id.csv')
PvD_trans <- read.csv('/Users/Jiaming/Desktop/BIO211_pre/ZebraFish/GO_KEGG/Datas/PvD_trans.csv')

PvD_trans$Expression <- PvD_trans$Proximal - PvD_trans$Distal
# 加载库
library(dplyr)
# 合并数据
merged_data <- merge(trans_id, PvD_trans, by.x = "x", by.y = "id")
head(merged_data)
# 构建gene向量
merged_data_entrez <- bitr(merged_data$x[-22], fromType = 'ENSEMBL', toType = 'ENTREZID', OrgDb = 'org.Dr.eg.db')
gene <- setNames(merged_data$Expression[-22], merged_data$symbol[-22]) # 替换为您的基因表达数据
# 查看结果
head(gene)




genedata<-data.frame(ID=info$gene_symbol,logFC=info$log2FoldChange)
write.table(GO$ONTOLOGY, file = "/Users/ZYP/Downloads/KEGG_GO/GO_ONTOLOGYs.txt", #将所有GO富集到的基因集所对应的类型写入本地文件从而得到BP/CC/MF各自的起始位置如我的数据里是1，2103，2410
            append = FALSE, quote = TRUE, sep = " ",
            eol = "\n", na = "NA", dec = ".", row.names = TRUE,
            col.names = TRUE, qmethod = c("escape", "double"),
            fileEncoding = "")

GOplotIn_BP<-GO[1:10,c(2,3,7,9)] #提取GO富集BP的前10行,提取ID,Description,p.adjust,GeneID四列
GOplotIn_CC<-GO[2103:2112,c(2,3,7,9)]#提取GO富集CC的前10行,提取ID,Description,p.adjust,GeneID四列
GOplotIn_MF<-GO[2410:2419,c(2,3,7,9)]#提取GO富集MF的前10行,提取ID,Description,p.adjust,GeneID四列
GOplotIn_BP$geneID <-str_replace_all(GOplotIn_BP$geneID,'/',',') #把GeneID列中的’/’替换成‘,’
GOplotIn_CC$geneID <-str_replace_all(GOplotIn_CC$geneID,'/',',')
GOplotIn_MF$geneID <-str_replace_all(GOplotIn_MF$geneID,'/',',')
names(GOplotIn_BP)<-c('ID','Term','adj_pval','Genes')#修改列名,后面弦图绘制的时候需要这样的格式
names(GOplotIn_CC)<-c('ID','Term','adj_pval','Genes')
names(GOplotIn_MF)<-c('ID','Term','adj_pval','Genes')
GOplotIn_BP$Category = "BP"#分类信息
GOplotIn_CC$Category = "CC"
GOplotIn_MF$Category = "MF"
circ_BP<-GOplot::circle_dat(GOplotIn_BP,genedata) #GOplot导入数据格式整理
circ_CC<-GOplot::circle_dat(GOplotIn_CC,genedata) 
circ_MF<-GOplot::circle_dat(GOplotIn_MF,genedata) 
chord_BP<-chord_dat(data = circ_BP,genes = genedata) #生成含有选定基因的数据框
chord_CC<-chord_dat(data = circ_CC,genes = genedata) 
chord_MF<-chord_dat(data = circ_MF,genes = genedata) 
GOChord(data = chord_BP,#弦图
        title = 'GO-Biological Process',space = 0.01,#GO Term间距
        limit = c(1,1),gene.order = 'logFC',gene.space = 0.25,gene.size = 5,
        lfc.col = c('red','white','blue'), #上下调基因颜色
        process.label = 10) #GO Term字体大小
GOChord(data = chord_CC,title = 'GO-Cellular Component',space = 0.01,
        limit = c(1,1),gene.order = 'logFC',gene.space = 0.25,gene.size = 5,
        lfc.col = c('red','white','blue'), 
        process.label = 10) 
GOChord(data = chord_MF,title = 'GO-Mollecular Function',space = 0.01,
        limit = c(1,1),gene.order = 'logFC',gene.space = 0.25,gene.size = 5,
        lfc.col = c('red','white','blue'), 
        process.label = 10)



