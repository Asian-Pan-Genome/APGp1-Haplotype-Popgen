library(ggplot2)
library(gggenes)
library(ggsci)
test <- read.table("hap.gggene.bed", header = TRUE, sep = "\t")
maxend <- max(test$end)
bound <- ceiling(maxend / 10000) * 10000
plot <-ggplot(test, aes(xmin = start, xmax = end, y = molecule, fill = gene, forward = orientation)) +
geom_gene_arrow() +
#geom_blank(data = dummies) +
facet_wrap(~ molecule, scales = "free", ncol = 1) +
scale_x_continuous(limits = c(0, bound), breaks=c(0,bound)) +
scale_color_aaas()+
theme_genes()
plot
