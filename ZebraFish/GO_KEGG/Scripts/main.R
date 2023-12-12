# 安装和加载必要的包
BiocManager::install("clusterProfiler")
BiocManager::install("org.Dr.eg.db") # 斑马鱼的注释包

library(BiocManager)
library(clusterProfiler)
library(org.Dr.eg.db)

# 准备基因列表
datas <- read.csv('/Users/Jiaming/Desktop/BIO211_pre/ZebraFish/GO_KEGG/Datas/bioDBnet.csv')
genes <- datas$Ensembl.Gene.ID # 用您的基因ID替换这里的示例

# GO富集分析
ego <- enrichGO(gene         = genes,
                OrgDb        = org.Dr.eg.db,
                keyType      = "ENSEMBL", # "SYMBOL"表示基因的标准名称或符号
                ont          = "CC", # 选择BP, CC或MF
                pAdjustMethod = "BH", # 使用Benjamini-Hochberg方法调整p值
                qvalueCutoff = 0.05, # 设置显著性阈值为0.05
                readable     = TRUE) # 将结果转换为可读的基因名

# KEGG富集分析
genes_2 <- bitr(genes, fromType = 'ENSEMBL', toType = 'ENTREZID', OrgDb = 'org.Dr.eg.db')
kegg <- enrichKEGG(gene = genes_2$ENTREZID, # KEGG富集分析
                   organism = "dre",
                   pAdjustMethod = "BH",
                   qvalueCutoff = 0.05)

# 结果可视化
library(ggplot2)
dotplot(ego) + ggtitle("GO Enrichment Analysis") # 加载ggplot2包，用于可视化
dotplot(kegg) + ggtitle("KEGG Enrichment Analysis") # 绘制KEGG富集分析的点图

# GO富集分析的柱状图
barplot(ego, showCategory=20) + ggtitle("GO Enrichment Analysis") # 显示前20个富集的类别
# KEGG富集分析的柱状图
barplot(kegg, showCategory=20) + ggtitle("KEGG Enrichment Analysis") # 显示前20个富集的类别

# pathview
# 安装和加载必要的包
library(pathview)

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
gene <- setNames(merged_data$Expression[-22], merged_data_entrez$ENTREZID) # 替换为您的基因表达数据
# 查看结果
head(gene)



library(KEGGREST)
zebrafish_pathways <- keggList("pathway", "dre") # 查询所有斑马鱼的KEGG通路
pathways_with_keyword <- zebrafish_pathways[grep("Retinol metabolism", zebrafish_pathways)]
# 搜索包含“Retinol metabolism”这一关键词的通路
pathway_id <- "00830" # 替换为查询到的通路ID

# 生成通路图
setwd('/Users/Jiaming/Desktop/BIO211_pre/ZebraFish/GO_KEGG/Datas/Pathways')
pathview(gene.data = gene,
         pathway.id = pathway_id,
         species = "dre", # 设置物种代码，"hsa"代表人类
         keytype = "ncbi", # 设置基因ID类型，这里使用ENTREZ ID
         limit = list(gene = max(abs(gene_data)), cpd = 1))

# Find ENSEMBL IDs with positive expression
gene_plus_ids <- merged_data_entrez$ENSEMBL[merged_data$Expression > 0]
# Find ENSEMBL IDs with negative expression
gene_minus_ids <- merged_data_entrez$ENSEMBL[merged_data$Expression < 0]
# 通过ENSEMBL ID获取正表达基因的符号
gene_plus_symbols <- merged_data$symbol[-22][merged_data_entrez$ENSEMBL %in% gene_plus_ids]
# 通过ENSEMBL ID获取负表达基因的符号
gene_minus_symbols <- merged_data$symbol[-22][merged_data_entrez$ENSEMBL %in% gene_minus_ids]
