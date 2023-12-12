library(pheatmap)
library(readxl)

proteins <- read_excel('/Users/jiaming/Desktop/BIO211_pre/ZebraFish/Proteomics/Datas/proteins.xlsx')
# colnames(proteins)
# [1] "Uniprot"               "Gene Name"            
# [3] "Fasta header"          "Sequence coverage [%]"
# [5] "Distal#1_Ave"          "Distal#2_Ave"         
# [7] "Distal#3_Ave"          "Proximal#1_Ave"       
# [9] "Proximal#2_Ave"        "Proximal#3_Ave"   

proteins$Distal <- rowMeans(proteins[,c(5,6,7)])
proteins$Proximal <- rowMeans(proteins[,c(8,9,10)])

write.csv(proteins, '/Users/jiaming/Desktop/BIO211_pre/ZebraFish/Proteomics/Datas/proteins.csv')

# 改变第二列的Gene Name为symbol

# 查找重复的基因
colnames(proteins)[2] <- "symbol"
duplicated_symbols <- proteins$symbol[duplicated(proteins$symbol)]
proteins <- proteins[!proteins$symbol %in% duplicated_symbols, ]
df <- data.frame(Proximal = proteins$Proximal,
                 Distal = proteins$Distal)
rownames(df) <- proteins$symbol

# 无标准化

# 对在远端和近端进行标准化以更好显示差异
# 聚类前
pheatmap(df, scale = "row", cluster_rows = FALSE)
# 聚类后
pheatmap(df, scale = "row", cluster_rows = TRUE)









