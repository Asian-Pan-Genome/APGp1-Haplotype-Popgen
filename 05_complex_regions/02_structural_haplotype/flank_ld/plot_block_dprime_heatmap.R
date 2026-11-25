#!/usr/bin/env Rscript

suppressPackageStartupMessages({
  if (!requireNamespace("ggplot2", quietly = TRUE)) {
    stop(
      "Package 'ggplot2' is required. Install it with: ",
      "install.packages('ggplot2')",
      call. = FALSE
    )
  }
})

usage <- function(status = 0L) {
  cat(paste0(
    "Draw a block-D heatmap from fixed_window_block_ld.py output.\n\n",
    "Usage:\n",
    "  Rscript plot_block_dprime_heatmap.R --input FILE --output FILE [options]\n\n",
    "Required:\n",
    "  --input FILE          TSV output from fixed_window_block_ld.py\n",
    "  --output FILE         Output image (.pdf, .png, .svg, ...); format follows suffix\n\n",
    "Options:\n",
    "  --value COLUMN        Value used for fill [multiallelic_dprime]\n",
    "                        Common choice: dprime_excess\n",
    "  --region-id ID        Plot only this region_id [all regions]\n",
    "  --layout MODE         auto, triangle, or cross [auto]\n",
    "  --width NUMBER        Figure width in inches [8]\n",
    "  --height NUMBER       Figure height in inches [7]\n",
    "  --dpi INTEGER         Raster resolution [300]\n",
    "  --no-labels           Do not print values inside cells\n",
    "  --include-non-ok      Plot values even when status is not 'ok'\n",
    "  -h, --help            Show this help\n\n",
    "By default, status != 'ok' and absent block pairs are grey.\n",
    "Input made with --pairing all-blocks is drawn as a triangular LD matrix.\n"
  ))
  quit(save = "no", status = status)
}

parse_args <- function(args) {
  opts <- list(
    input = NULL,
    output = NULL,
    value = "multiallelic_dprime",
    region_id = NULL,
    layout = "auto",
    width = 8,
    height = 7,
    dpi = 300L,
    labels = TRUE,
    include_non_ok = FALSE
  )

  i <- 1L
  while (i <= length(args)) {
    arg <- args[[i]]
    if (arg %in% c("-h", "--help")) usage(0L)
    if (arg == "--no-labels") {
      opts$labels <- FALSE
      i <- i + 1L
      next
    }
    if (arg == "--include-non-ok") {
      opts$include_non_ok <- TRUE
      i <- i + 1L
      next
    }
    if (!arg %in% c(
      "--input", "--output", "--value", "--region-id", "--layout",
      "--width", "--height", "--dpi"
    )) {
      stop("Unknown option: ", arg, call. = FALSE)
    }
    if (i == length(args)) stop("Missing value after ", arg, call. = FALSE)
    value <- args[[i + 1L]]
    if (arg == "--input") opts$input <- value
    if (arg == "--output") opts$output <- value
    if (arg == "--value") opts$value <- value
    if (arg == "--region-id") opts$region_id <- value
    if (arg == "--layout") opts$layout <- value
    if (arg == "--width") opts$width <- as.numeric(value)
    if (arg == "--height") opts$height <- as.numeric(value)
    if (arg == "--dpi") opts$dpi <- as.integer(value)
    i <- i + 2L
  }

  if (is.null(opts$input) || is.null(opts$output)) usage(1L)
  if (!file.exists(opts$input)) {
    stop("Input file does not exist: ", opts$input, call. = FALSE)
  }
  if (!opts$layout %in% c("auto", "triangle", "cross")) {
    stop("--layout must be auto, triangle, or cross", call. = FALSE)
  }
  if (!is.finite(opts$width) || opts$width <= 0 ||
      !is.finite(opts$height) || opts$height <= 0 ||
      is.na(opts$dpi) || opts$dpi <= 0) {
    stop("--width, --height, and --dpi must be positive numbers", call. = FALSE)
  }
  opts
}

opts <- parse_args(commandArgs(trailingOnly = TRUE))

dat <- read.delim(
  opts$input,
  header = TRUE,
  sep = "\t",
  quote = "",
  comment.char = "",
  stringsAsFactors = FALSE,
  check.names = FALSE,
  na.strings = c("NA", "NaN", "")
)

required <- c("region_id", opts$value)
missing_columns <- setdiff(required, names(dat))
if (length(missing_columns) > 0L) {
  stop(
    "Input is missing required column(s): ",
    paste(missing_columns, collapse = ", "),
    call. = FALSE
  )
}

if (!is.null(opts$region_id)) {
  dat <- dat[dat$region_id == opts$region_id, , drop = FALSE]
  if (nrow(dat) == 0L) {
    stop("No rows found for region_id: ", opts$region_id, call. = FALSE)
  }
}
if (nrow(dat) == 0L) stop("Input contains no data rows", call. = FALSE)

dat$region_id <- as.character(dat$region_id)
dat$plot_value <- suppressWarnings(as.numeric(dat[[opts$value]]))

if (!opts$include_non_ok && "status" %in% names(dat)) {
  dat$plot_value[is.na(dat$status) | dat$status != "ok"] <- NA_real_
}

region_ids <- unique(dat$region_id)
triangle_columns <- c(
  "block1_side", "block1_label", "block1_genomic_index",
  "block2_side", "block2_label", "block2_genomic_index"
)
has_triangle_columns <- all(triangle_columns %in% names(dat))
layout <- opts$layout
if (layout == "auto") {
  has_within_side_pairs <- has_triangle_columns && any(
    dat$block1_side == dat$block2_side,
    na.rm = TRUE
  )
  layout <- if (has_within_side_pairs) "triangle" else "cross"
}

common_theme <- ggplot2::theme_minimal(base_size = 11) +
  ggplot2::theme(
    panel.grid = ggplot2::element_blank(),
    axis.ticks = ggplot2::element_line(color = "grey50"),
    strip.text = ggplot2::element_text(face = "bold"),
    plot.title.position = "plot",
    legend.position = "right"
  )

if (layout == "triangle") {
  missing_triangle_columns <- setdiff(triangle_columns, names(dat))
  if (length(missing_triangle_columns) > 0L) {
    stop(
      "Triangular layout requires column(s): ",
      paste(missing_triangle_columns, collapse = ", "),
      call. = FALSE
    )
  }
  dat$block1_genomic_index <- suppressWarnings(as.integer(dat$block1_genomic_index))
  dat$block2_genomic_index <- suppressWarnings(as.integer(dat$block2_genomic_index))
  bad_index <- is.na(dat$block1_genomic_index) |
    is.na(dat$block2_genomic_index) |
    dat$block1_genomic_index < 1L |
    dat$block2_genomic_index < 1L
  if (any(bad_index)) {
    stop("Genomic block indices must be positive integers", call. = FALSE)
  }

  dat$x_index <- pmin(dat$block1_genomic_index, dat$block2_genomic_index)
  dat$y_index <- pmax(dat$block1_genomic_index, dat$block2_genomic_index)
  key <- paste(dat$region_id, dat$x_index, dat$y_index, sep = "\r")
  if (anyDuplicated(key)) {
    stop("Input contains duplicate unordered block pairs", call. = FALSE)
  }

  block_map <- unique(rbind(
    data.frame(
      block_index = dat$block1_genomic_index,
      block_label = dat$block1_label,
      block_side = dat$block1_side,
      stringsAsFactors = FALSE
    ),
    data.frame(
      block_index = dat$block2_genomic_index,
      block_label = dat$block2_label,
      block_side = dat$block2_side,
      stringsAsFactors = FALSE
    )
  ))
  mapping_key <- paste(block_map$block_index, block_map$block_label, sep = "\r")
  block_map <- block_map[!duplicated(mapping_key), , drop = FALSE]
  labels_per_index <- table(block_map$block_index)
  if (any(labels_per_index != 1L)) {
    stop("A genomic block index maps to more than one block label", call. = FALSE)
  }
  block_map <- block_map[order(block_map$block_index), , drop = FALSE]
  block_count <- max(block_map$block_index)
  block_labels <- block_map$block_label[match(seq_len(block_count), block_map$block_index)]

  grids <- lapply(region_ids, function(region_id) {
    grid <- expand.grid(
      region_id = region_id,
      x_index = seq_len(block_count),
      y_index = seq_len(block_count),
      stringsAsFactors = FALSE
    )
    grid[grid$x_index <= grid$y_index, , drop = FALSE]
  })
  grid <- do.call(rbind, grids)
  plot_dat <- merge(
    grid,
    dat[, c("region_id", "x_index", "y_index", "plot_value")],
    by = c("region_id", "x_index", "y_index"),
    all.x = TRUE,
    sort = FALSE
  )
  plot_dat$region_id <- factor(plot_dat$region_id, levels = region_ids)
  left_count <- length(unique(block_map$block_index[block_map$block_side == "left"]))

  p <- ggplot2::ggplot(
    plot_dat,
    ggplot2::aes(x = x_index, y = y_index, fill = plot_value)
  ) +
    ggplot2::geom_tile(color = "white", linewidth = 0.3) +
    ggplot2::scale_x_continuous(
      breaks = seq_len(block_count),
      labels = block_labels,
      position = "top",
      expand = ggplot2::expansion(mult = 0)
    ) +
    ggplot2::scale_y_reverse(
      breaks = seq_len(block_count),
      labels = block_labels,
      expand = ggplot2::expansion(mult = 0)
    ) +
    ggplot2::coord_fixed() +
    ggplot2::labs(
      title = paste("Triangular flanking block LD:", opts$value),
      subtitle = "Genomic order is Lmax ... L1 | R1 ... Rmax; grey = failed QC or unavailable",
      x = NULL,
      y = NULL,
      fill = opts$value
    ) +
    common_theme +
    ggplot2::theme(
      axis.text.x = ggplot2::element_text(angle = 45, hjust = 0)
    )
  if (left_count > 0L && left_count < block_count) {
    p <- p +
      ggplot2::geom_vline(xintercept = left_count + 0.5, linewidth = 0.45) +
      ggplot2::geom_hline(yintercept = left_count + 0.5, linewidth = 0.45)
  }
} else {
  cross_columns <- c("left_block_rank", "right_block_rank")
  missing_cross_columns <- setdiff(cross_columns, names(dat))
  if (length(missing_cross_columns) > 0L) {
    stop(
      "Cross layout requires column(s): ",
      paste(missing_cross_columns, collapse = ", "),
      call. = FALSE
    )
  }
  dat$left_block_rank <- suppressWarnings(as.integer(dat$left_block_rank))
  dat$right_block_rank <- suppressWarnings(as.integer(dat$right_block_rank))
  bad_rank <- is.na(dat$left_block_rank) | is.na(dat$right_block_rank) |
    dat$left_block_rank < 1L | dat$right_block_rank < 1L
  if (any(bad_rank)) stop("Block ranks must be positive integers", call. = FALSE)

  key <- paste(dat$region_id, dat$left_block_rank, dat$right_block_rank, sep = "\r")
  if (anyDuplicated(key)) {
    stop("Input contains duplicate left/right block pairs", call. = FALSE)
  }
  grids <- lapply(region_ids, function(region_id) {
    region_data <- dat[dat$region_id == region_id, , drop = FALSE]
    max_left <- max(region_data$left_block_rank)
    max_right <- max(region_data$right_block_rank)
    expected_pairs <- max_left * max_right
    if (nrow(region_data) < expected_pairs) {
      warning(
        sprintf(
          paste0(
            "Region '%s' has %d of %d possible L-R pairs; ",
            "use --pairing all upstream for a complete cross heatmap."
          ),
          region_id, nrow(region_data), expected_pairs
        ),
        call. = FALSE
      )
    }
    expand.grid(
      region_id = region_id,
      left_block_rank = seq_len(max_left),
      right_block_rank = seq_len(max_right),
      stringsAsFactors = FALSE
    )
  })
  grid <- do.call(rbind, grids)
  plot_dat <- merge(
    grid,
    dat[, c("region_id", "left_block_rank", "right_block_rank", "plot_value")],
    by = c("region_id", "left_block_rank", "right_block_rank"),
    all.x = TRUE,
    sort = FALSE
  )
  plot_dat$region_id <- factor(plot_dat$region_id, levels = region_ids)
  all_left_ranks <- sort(unique(plot_dat$left_block_rank))
  all_right_ranks <- sort(unique(plot_dat$right_block_rank))

  p <- ggplot2::ggplot(
    plot_dat,
    ggplot2::aes(x = right_block_rank, y = left_block_rank, fill = plot_value)
  ) +
    ggplot2::geom_tile(color = "white", linewidth = 0.35) +
    ggplot2::scale_x_continuous(
      breaks = all_right_ranks,
      position = "top",
      expand = ggplot2::expansion(mult = 0)
    ) +
    ggplot2::scale_y_reverse(
      breaks = all_left_ranks,
      expand = ggplot2::expansion(mult = 0)
    ) +
    ggplot2::coord_fixed() +
    ggplot2::labs(
      title = paste("Left-by-right flanking block LD:", opts$value),
      subtitle = "Rank 1 is nearest to the target region; grey = absent or failed QC",
      x = "Right block rank (near to far)",
      y = "Left block rank (near to far)",
      fill = opts$value
    ) +
    common_theme
}

label_dat <- plot_dat[!is.na(plot_dat$plot_value), , drop = FALSE]
label_dat$cell_label <- formatC(label_dat$plot_value, format = "f", digits = 2)

finite_values <- plot_dat$plot_value[is.finite(plot_dat$plot_value)]
is_diverging <- identical(opts$value, "dprime_excess") ||
  (length(finite_values) > 0L && any(finite_values < 0))

if (is_diverging) {
  color_limit <- if (length(finite_values) > 0L) max(abs(finite_values)) else 1
  if (!is.finite(color_limit) || color_limit == 0) color_limit <- 1
  p <- p + ggplot2::scale_fill_gradient2(
    low = "#2166AC",
    mid = "white",
    high = "#B2182B",
    midpoint = 0,
    limits = c(-color_limit, color_limit),
    na.value = "grey85",
    oob = scales::squish
  )
} else {
  fixed_01 <- opts$value %in% c(
    "multiallelic_dprime", "dprime_permutation_mean", "dprime_empirical_p"
  )
  fill_limits <- if (fixed_01) c(0, 1) else NULL
  p <- p + ggplot2::scale_fill_viridis_c(
    option = "C",
    direction = -1,
    limits = fill_limits,
    na.value = "grey85",
    oob = scales::squish
  )
}

if (opts$labels && nrow(label_dat) > 0L) {
  p <- p + ggplot2::geom_text(
    data = label_dat,
    ggplot2::aes(label = cell_label),
    size = 2.7,
    color = "black",
    na.rm = TRUE
  )
}

if (length(region_ids) > 1L) {
  p <- p + ggplot2::facet_wrap(ggplot2::vars(region_id))
} else {
  p <- p + ggplot2::labs(caption = paste("Region:", region_ids[[1L]]))
}

output_dir <- dirname(normalizePath(opts$output, mustWork = FALSE))
if (!dir.exists(output_dir)) {
  dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)
}

ggplot2::ggsave(
  filename = opts$output,
  plot = p,
  width = opts$width,
  height = opts$height,
  units = "in",
  dpi = opts$dpi,
  bg = "white"
)

message(sprintf(
  "Saved heatmap for %d region(s) to %s",
  length(region_ids), opts$output
))
