library(pheatmap)

trans <- read.csv('/Users/jiaming/Desktop/BIO211_pre/ZebraFish/Transcriptomics/Datas/PvD_trans.csv', header = T, stringsAsFactors = F)

# 查找重复的基因
duplicated_symbols <- trans$symbol[duplicated(trans$symbol)]
trans <- trans[!trans$symbol %in% duplicated_symbols, ]
df <- data.frame(Proximal = trans$Proximal,
                 Distal = trans$Distal)
rownames(df) <- trans$symbol

# 无标准化

# 对在远端和近端进行标准化以更好显示差异
# 聚类前
pheatmap(df, scale = "row", cluster_rows = FALSE)
# 聚类后
pheatmap(df, scale = "row", cluster_rows = TRUE)









