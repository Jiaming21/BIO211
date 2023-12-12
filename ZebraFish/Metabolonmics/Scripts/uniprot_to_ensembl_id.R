meta <- read.csv("/Users/Jiaming/Desktop/BIO211_pre/ZebraFish/Metabolonmics/Datas/compound_id.csv")
library(biomaRt)
# 创建biomaRt对象，选择要查询的数据库和数据集
mart <- useMart("ensembl", dataset = "drerio_gene_ensembl")
# 定义输入的Uniprot ID列表
uniprot_ids <- c("P12345", "Q98765", "O87654")
# 查询UniProt ID对应的Ensembl ID
ensembl_ids <- getBM(attributes = c("uniprot_id","ensembl_gene_id"),
                     filters = "uniprotswissprot",
                     values = uniprot_ids,
                     mart = mart)
# 提取Ensembl ID列表
ensembl_id_list <- ensembl_ids$ensembl_gene_id
# 进行GO富集分析
ego <- enrichGO(gene = ensembl_id_list, OrgDb = org.Dr.eg.db, keyType = "ENSEMBL", ont = "BP")
print(ego)
# 进行KEGG富集分析
kegg <- enrichKEGG(gene = ensembl_id_list, organism = "drerio", keyType = "kegg")
print(kegg)

