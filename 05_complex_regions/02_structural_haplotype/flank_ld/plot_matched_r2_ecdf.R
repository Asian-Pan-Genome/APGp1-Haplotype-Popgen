#!/usr/bin/env Rscript
# Compare target and equally weighted matched-background ECDFs (base R only).

usage <- function() {
  cat("Usage: Rscript plot_matched_r2_ecdf.R INPUT_DIR OUTPUT_PREFIX [--backgrounds-per-target 10,20,50,all] [--ks-permutations 1000] [--seed 1] [--include-incomplete]\n", file=stderr())
}
args <- commandArgs(trailingOnly=TRUE)
if (length(args) < 2L) { usage(); quit(status=2L) }
input_dir <- args[1]
output_prefix <- args[2]
per_target <- "20"
seed <- 1L
ks_permutations <- 1000L
include_incomplete <- FALSE
i <- 3L
while (i <= length(args)) {
  key <- args[i]
  if (key == "--include-incomplete") {
    include_incomplete <- TRUE
    i <- i + 1L
  } else if (key %in% c("--backgrounds-per-target", "--ks-permutations", "--seed")) {
    if (i == length(args)) stop("Missing value for ", key)
    value <- args[i + 1L]
    if (key == "--backgrounds-per-target") per_target <- value
    if (key == "--ks-permutations") ks_permutations <- as.integer(value)
    if (key == "--seed") seed <- as.integer(value)
    i <- i + 2L
  } else if (key %in% c("-h", "--help")) {
    usage(); quit(status=0L)
  } else stop("Unknown option: ", key)
}
if (is.na(seed)) stop("--seed must be an integer")
if (is.na(ks_permutations) || ks_permutations < 1L) stop("--ks-permutations must be positive")
sample_labels <- unique(trimws(strsplit(per_target, ",", fixed=TRUE)[[1L]]))
if (!length(sample_labels) || any(!nzchar(sample_labels))) stop("Empty background sample size")
sample_sizes <- suppressWarnings(as.integer(sample_labels))
invalid_size <- sample_labels != "all" & (is.na(sample_sizes) | sample_sizes < 1L)
if (any(invalid_size)) stop("--backgrounds-per-target requires positive integers or all")
if (!dir.exists(input_dir)) stop("Input directory does not exist: ", input_dir)
files <- list.files(input_dir, pattern="_flanking_r2[.]summary[.]tsv$", full.names=TRUE)
if (!length(files)) stop("No *_flanking_r2.summary.tsv files found")
read_result_tsv <- function(path, nrows=-1L) {
  table <- read.delim(path, nrows=nrows, check.names=FALSE,
                      na.strings=c("NA", ""), blank.lines.skip=TRUE)
  if (length(names(table))) names(table)[1] <- sub("^[^A-Za-z_]+", "", names(table)[1])
  table
}
set.seed(seed)
targets <- numeric()
background_groups <- list()
used <- 0L
skipped <- 0L
for (summary_path in files) {
  background_path <- sub("[.]summary[.]tsv$", ".backgrounds.tsv", summary_path)
  if (!file.exists(background_path)) { skipped <- skipped + 1L; next }
  summary <- read_result_tsv(summary_path, nrows=1L)
  if (nrow(summary) != 1L || !"status" %in% names(summary)) stop("Invalid summary: ", summary_path)
  if (!include_incomplete && as.character(summary$status[1]) != "ok") {
    skipped <- skipped + 1L; next
  }
  if (!"median_r2" %in% names(summary)) stop("Missing median_r2 in: ", summary_path)
  target <- suppressWarnings(as.numeric(summary$median_r2[1]))
  if (!is.finite(target)) { skipped <- skipped + 1L; next }
  bg <- read_result_tsv(background_path)
  if (!"median_r2" %in% names(bg)) stop("Missing median_r2 in: ", background_path)
  keep <- is.finite(suppressWarnings(as.numeric(bg$median_r2)))
  if (!include_incomplete && "full_flanks" %in% names(bg)) {
    keep <- keep & tolower(as.character(bg$full_flanks)) == "true"
  }
  values <- as.numeric(bg$median_r2[keep])
  if (!length(values)) { skipped <- skipped + 1L; next }
  targets <- c(targets, target)
  background_groups[[length(background_groups) + 1L]] <- values
  used <- used + 1L
}
if (!used) stop("No usable target/background pairs")
dir.create(dirname(output_prefix), recursive=TRUE, showWarnings=FALSE)

# Shuffle once per target: size 10 is a subset of size 20, and so on.
shuffled_groups <- lapply(background_groups, function(values) values[sample.int(length(values))])
select_groups <- function(label) {
  if (label == "all") return(shuffled_groups)
  n <- as.integer(label)
  lapply(shuffled_groups, function(values) values[seq_len(min(n, length(values)))])
}
flatten_groups <- function(groups) {
  values <- unlist(groups, use.names=FALSE)
  weights <- unlist(lapply(groups, function(group) rep(1 / (used * length(group)), length(group))),
                    use.names=FALSE)
  list(values=values, weights=weights)
}
weighted_ecdf <- function(points, values, weights) {
  order_values <- order(values)
  sorted_values <- values[order_values]
  cumulative <- cumsum(weights[order_values])
  index <- findInterval(points, sorted_values)
  ifelse(index == 0L, 0, cumulative[pmax(index, 1L)])
}
all_values <- unlist(background_groups, use.names=FALSE)
plot_max <- max(1e-6, as.numeric(quantile(c(targets, all_values), 0.995, names=FALSE)))
bandwidth_data <- c(targets, unlist(lapply(shuffled_groups, function(values)
  values[seq_len(min(20L, length(values)))]), use.names=FALSE))
common_bandwidth <- max(1e-6, bw.nrd0(bandwidth_data))
rm(all_values, bandwidth_data)

x <- seq(0, 1, length.out=1001L)
target_ecdf <- ecdf(targets)(x)
ecdf_table <- data.frame(r2=x, target_ecdf=target_ecdf)
density_x <- seq(0, plot_max, length.out=512L)
target_density <- density(targets, weights=rep(1 / used, used), bw=common_bandwidth,
                          from=0, to=plot_max, n=512L)$y
density_table <- data.frame(r2=density_x, target_density=target_density)
colors <- grDevices::rainbow(length(sample_labels), s=0.75, v=0.75)
ks_rows <- vector("list", length(sample_labels))

for (group_index in seq_along(sample_labels)) {
  label <- sample_labels[group_index]
  groups <- select_groups(label)
  pooled <- flatten_groups(groups)
  values <- pooled$values
  weights <- pooled$weights
  ecdf_column <- paste0("background_ecdf_n_", label)
  density_column <- paste0("background_density_n_", label)
  ecdf_table[[ecdf_column]] <- weighted_ecdf(x, values, weights)
  density_table[[density_column]] <- density(values, weights=weights,
                                             bw=common_bandwidth, from=0,
                                             to=plot_max, n=512L)$y

  # Exact observed D on the union of target and sampled background values.
  support <- sort(unique(c(targets, values)))
  background_at_support <- weighted_ecdf(support, values, weights)
  target_at_support <- findInterval(support, sort(targets)) / used
  weighted_D <- max(abs(target_at_support - background_at_support))
  null_D <- numeric(ks_permutations)
  for (iteration in seq_len(ks_permutations)) {
    pseudo_targets <- vapply(groups, function(group)
      group[sample.int(length(group), 1L)], numeric(1L))
    pseudo_at_support <- findInterval(support, sort(pseudo_targets)) / used
    null_D[iteration] <- max(abs(pseudo_at_support - background_at_support))
  }
  matched_p <- (1 + sum(null_D >= weighted_D)) / (ks_permutations + 1)
  nominal <- suppressWarnings(ks.test(targets, values, exact=FALSE))
  ks_rows[[group_index]] <- data.frame(
    backgrounds_per_target=label, n_targets=used, n_background_values=length(values),
    ks_D_weighted=weighted_D, matched_empirical_p=matched_p,
    nominal_unweighted_ks_D=as.numeric(nominal$statistic),
    nominal_unweighted_ks_p=nominal$p.value, ks_permutations=ks_permutations)
}

write.table(ecdf_table, paste0(output_prefix, ".ecdf.tsv"),
            sep="\t", row.names=FALSE, quote=FALSE)
write.table(density_table, paste0(output_prefix, ".distribution.tsv"),
            sep="\t", row.names=FALSE, quote=FALSE)
write.table(do.call(rbind, ks_rows), paste0(output_prefix, ".ks.tsv"),
            sep="\t", row.names=FALSE, quote=FALSE)

pdf(paste0(output_prefix, ".ecdf.pdf"), width=7, height=5)
plot(x, target_ecdf, type="l", lwd=2, col="#000000", xlim=c(0, 1), ylim=c(0, 1),
     xlab="Flanking median r²", ylab="ECDF",
     main=sprintf("Matched flanking LD: %d target regions", used))
for (group_index in seq_along(sample_labels)) {
  lines(x, ecdf_table[[paste0("background_ecdf_n_", sample_labels[group_index])]],
        lwd=2, col=colors[group_index])
}
legend("bottomright", c(sprintf("Targets (n=%d)", used),
                        paste0("Background / target: ", sample_labels)),
       col=c("#000000", colors), lwd=2, bty="n", cex=0.8)
dev.off()

pdf(paste0(output_prefix, ".distribution.pdf"), width=7, height=5)
plot(density_x, target_density, type="l", lwd=2, col="#000000",
     xlim=c(0, plot_max), ylim=c(0, max(as.matrix(density_table[-1L])) * 1.08),
     xlab="Flanking median r²", ylab="Density",
     main=sprintf("Median r² density: %d target regions", used),
     sub="Display range ends at pooled 99.5th percentile")
for (group_index in seq_along(sample_labels)) {
  lines(density_x,
        density_table[[paste0("background_density_n_", sample_labels[group_index])]],
        lwd=2, col=colors[group_index])
}
legend("topright", c("Targets", paste0("Background / target: ", sample_labels)),
       col=c("#000000", colors), lwd=2, bty="n", cex=0.8)
dev.off()
cat(sprintf("Wrote ECDF, density, and KS outputs at %s; targets=%d, sample sizes=%s, skipped=%d\n",
            output_prefix, used, paste(sample_labels, collapse=","), skipped))
