library(HaploBlocker)
maptbl <- read.table("case.map", sep="\t", header=FALSE)
map <- maptbl$V4 - min(maptbl$V4)
numbers <- scan("order", what = numeric())
order <- as.character(numbers)
blockl <- block_calculation("case.ped", target_coverage = 0.999, min_share = 0.99, bp_map = map, window_size = 8, overlap_remove = TRUE)

pdf("no_overlap.pdf", width = 8, height = 8)
fine_order <- plot_block(blockl, type="bp", export_order = TRUE, orientation="mid") #orientation="mid")#,  import_order=order$V2)
#fine_order <- plot_block(blockl, type="bp", import_order = order) #orientation="mid")#,  import_order=order$V2)
cat(unlist(fine_order), sep = "\n")
dev.off()

blockl <- block_calculation("./case.ped", target_coverage = 0.999, min_share = 0.99, bp_map = map, window_size = 8, overlap_remove = FALSE)

pdf("overlap.pdf", width = 8, height = 8)
fine_order <- plot_block(blockl, type="bp", export_order = TRUE, orientation="mid") #orientation="mid")#,  import_order=order$V2)

dev.off()
