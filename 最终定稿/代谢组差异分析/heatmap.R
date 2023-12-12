library(pheatmap)

trans <- read.csv('/Users/jiaming/Desktop/代谢组差异分析/meta_for_heatmap.csv', header = T, stringsAsFactors = F)

# 查找重复的基因
df <- data.frame(Proximal = trans$P,
                 Distal = trans$D)
rownames(df) <- trans$meta

# 无标准化

# 对在远端和近端进行标准化以更好显示差异
# 聚类前
pheatmap(df, scale = "row", cluster_rows = FALSE)
# 聚类后
pheatmap(df, scale = "row", cluster_rows = TRUE)









