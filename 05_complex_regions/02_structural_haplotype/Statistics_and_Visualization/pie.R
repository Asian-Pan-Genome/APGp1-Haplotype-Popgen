library(ggplot2)
library(ggpubr)

# Create Data
data <-  read.table("pop.stat", header = TRUE, sep = "\t", comment.char = "")
pop <- colnames(data)[!colnames(data) %in% "HAP"]
# Basic piechart

plots <- lapply(1:length(pop), function(i) {
popname <- pop[i]
ggplot(data, aes_string(x = 1, y = popname, fill="HAP")) +
  geom_bar(stat="identity", width=1) +
  coord_polar("y", start=0) +
  theme_void() + ylab("") + xlab("") + ggtitle(popname)

})

p <-  ggplot(data, aes(x = HAP, y = EAS, color = HAP)) +
    geom_point(size = 3) +
    labs(color = "Experimental Groups") +
    theme(legend.position = "bottom")

legend <- get_legend(p)

ggarrange(plotlist=plots, legend.grob = legend, ncol = length(pop), nrow = 1)
