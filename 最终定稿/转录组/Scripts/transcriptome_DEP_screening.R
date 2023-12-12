# Task 1.2
# 辨认出远端和近段富集的蛋白质
# Read data
path <- '/Users/jiaming/Desktop/BIO211_pre/Datas'
setwd(paste(path))
rm(list = ls()) # 清空environment
protein <- read.csv('PvD_combine.csv', header = T, sep = ',', stringsAsFactors = F, 
                    quote = "", row.names = 1)
LFQ <- protein[, -1] # 去除基因那一列
P_D <- 2^LFQ

# DEP analysis
pvalue <- apply(P_D[, c(1:30)], 1, function(x) {
  a <- factor(c(rep("Distal", 15), rep("Proximal", 15)))
  t.test(x~a, var.equal = T, alternative = "two.sided", conf.level = 0.95)
})

result_t.test <- data.frame(gene_symbol=protein$Gene.Name, 
                            Pvalue = as.numeric(unlist(lapply(pvalue, function(x) x$p.value))),
                            log2FC = log2(as.numeric(unlist(lapply(pvalue, 
                                                                   function(x) x$estimate[1]/x$estimate[2]))))
)

rownames(result_t.test) <- rownames(protein)
t.test_fdr <- p.adjust(result_t.test$Pvalue, method = 'fdr')
result_t.test_fdr <- cbind(result_t.test, t.test_fdr)
dep <- result_t.test_fdr[result_t.test_fdr$t.test_fdr < 0.05 & abs(result_t.test_fdr$log2FC) > 1, ]
dep1 <- dep[order(abs(dep$log2FC), decreasing=T),] 
write.table(dep1, file = 'dep_fdr_raw.txt', sep = '\t',quote = F, row.names = F)

#######################      远端和近端差异化表达       ##################### 
# volcano plot: quickly identify changes in large data sets
# Volcano plot
# install.packages('ggplot2')
# install.packages('cowplot')
# install.packages('patchwork')
# install.packages('ggplotify')
# install.packages("ggrepel")
library(cowplot)
library(patchwork)
library(ggplotify)
library(ggplot2)
library(ggrepel)

logFC_cutoff <- 1
result_t.test_fdr <- result_t.test_fdr[which (result_t.test_fdr$gene_symbol != ""), ] 
# omit rows with unknown ID 












anno <- read.csv('volcano_anno.csv', header = T, sep = ',', stringsAsFactors = F, 
                 quote = "", row.names = 1)
result_t.test_fdr$change = as.factor(
  ifelse (result_t.test_fdr$t.test_fdr < 0.05 & abs(result_t.test_fdr$log2FC) > logFC_cutoff,
          ifelse (result_t.test_fdr$log2FC > logFC_cutoff,"UP","DOWN"),"NOT")
)

up_down <- table(result_t.test_fdr$change)
loc_up <- which(result_t.test_fdr$change == 'UP')
loc_down <- which(result_t.test_fdr$change == 'DOWN')
significant <- rep('normal', times=nrow(result_t.test_fdr))
significant[loc_up] <- 'up'
significant[loc_down] <- 'down'
significant <- factor(significant, levels = c('up', 'down', 'normal'))
result_t.test_fdr_new = cbind(result_t.test_fdr, anno)
p <- qplot(x = result_t.test_fdr_new$log2FC, 
           y = -log10(result_t.test_fdr_new$t.test_fdr), 
           xlab = 'log2(FC)', ylab = '-log10(FDR)', 
           size=I(2), alpha = I(1/3), colour=significant)
p <- p+scale_color_manual(values = c('up'='red', 'normal' = 'gray', 'down'='blue'))

# Add cutoff line
xline = c(-logFC_cutoff,logFC_cutoff)
p <- p+geom_vline(xintercept = xline, lty = 2, size = I(0.2), color = 'grey11')
yline = -log(0.05,10)
p <- p+geom_hline(yintercept = yline, lty = 2, size = I(0.2), color = 'grey11')

# Label gene names of interest
p <- p+geom_text_repel(aes(label = result_t.test_fdr_new$label), 
                       point.padding = unit(0.25, "lines"), 
                       arrow = arrow(length = unit(0.01, "npc")), 
                       nudge_y = 0.1)
p <- p+theme_bw() + theme(panel.background = element_rect(colour = 'black', 
                                                          size = 1, fill = 'white'), 
                          panel.grid = element_blank(), 
                          axis.title.x = element_text(size = 15) , 
                          axis.title.y = element_text(size = 15), 
                          axis.text.x = element_text(size = 15), 
                          axis.text.y = element_text(size = 15), 
                          legend.title = element_text(size = 12), 
                          legend.text = element_text(size = 10))

ggsave('volcano_dep.pdf')
