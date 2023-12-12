# getwd(): 获取工作目录
trans <- read.csv('/Users/jiaming/Desktop/BIO211_pre/ZebraFish/Transcriptomics/Datas/PvD_trans.csv', header = T, stringsAsFactors = F)

df <- data.frame(symbol = trans$symbol,
                 logFC = trans$log2FoldChange,
                 P.Value = trans$pval)

# Distal vs Proximal
df$group <- ifelse(df$logFC > log2(3) & df$P.Value <= 0.05, "Distal",
                   ifelse(df$logFC <= -log2(3) & df$P.Value <= 0.05, "Proximal", "Normal"))

# Labeled
df$pvalue_log10 <- -log10(df$P.Value)
df1 <- df[df$pvalue_log10 >= 100,]

library(ggplot2)
library(ggrepel)

ggplot(df,aes(x=logFC, y=-log10(P.Value)))+
  geom_point(aes(color=group))+
  scale_color_manual(values=c("dodgerblue", "gray", "firebrick"))+
  geom_label_repel(data=df1,aes(x=logFC, y=-log(P.Value),label=symbol))+
  theme_bw()+
  xlim(-30,30)+
  ylim(-5,500)





