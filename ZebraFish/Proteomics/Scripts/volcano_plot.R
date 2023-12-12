# getwd(): 获取工作目录
proteins <- read.csv('/Users/jiaming/Desktop/BIO211_pre/ZebraFish/Proteomics/Datas/proteins.csv')

df <- data.frame(symbol = proteins$symbol,
                 DF = proteins$Difference,
                 minus_log10_p_value = proteins$minus_LOG10_p_value,
                 p_value = proteins$p_value)

# Distal vs Proximal
df$group <- ifelse(df$DF > log2(3) & df$minus_log10_p_value >= -log10(0.05), "Distal",
                   ifelse(df$DF <= -log2(3) & df$minus_log10_p_value >= -log10(0.05), "Proximal", "Normal"))

# df$group <- ifelse(df$DF > log2(3) & df$minus_log10_p_value <= -log10(0.05), "Distal",
#                    ifelse(df$DF <= -log2(3) & df$minus_log10_p_value <= -log10(0.05), "Proximal", "Normal"))
# 
# df$group <- ifelse(df$p_value <= 0.05, "Distal",
#                    ifelse(df$p_value <= -log10(0.05), "Proximal", "Normal"))

# Labeled
df1 <- df[df$minus_log10_p_value >= -log10(0.05),]

library(ggplot2)
library(ggrepel)

ggplot(df,aes(x=DF, y=minus_log10_p_value))+
  geom_point(aes(color=group))+
  scale_color_manual(values=c("dodgerblue", "gray", "firebrick"))+
  geom_label_repel(data=df1,aes(x=DF, y=minus_log10_p_value, label=symbol))+
  theme_bw()+
  xlim(-10,10)+
  ylim(0,5)








