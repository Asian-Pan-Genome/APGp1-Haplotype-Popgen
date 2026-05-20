#!/usr/bin/env Rscript
# Step 4: SV PheWAS Analysis
# ============================================
# Run phenome-wide association scans for each extracted SV using the PheWAS package.
# Model: logistic regression
# Formula: PheCode ~ SV_Genotype + Age + Sex
#
# Usage:
#   Rscript 04.sv_phewas.R \
#     --geno_cov 03.results/03.phewas_data/phewas_covariates_genotypes.tsv \
#     --icd10 03.results/03.phewas_data/filtered_icd10.csv \
#     --outdir 03.results/04.phewas

library(optparse)
library(PheWAS)
library(dplyr)
library(tidyr)

option_list = list(
  make_option(c("--geno_cov"), type="character", default=NULL, help="Path to phewas_covariates_genotypes.tsv", metavar="character"),
  make_option(c("--icd10"), type="character", default=NULL, help="Path to filtered_icd10.csv", metavar="character"),
  make_option(c("--outdir"), type="character", default=NULL, help="Output directory", metavar="character"),
  make_option(c("--sv_index"), type="integer", default=NULL, help="Run only for the Nth SV (for array jobs)", metavar="integer"),
  make_option(c("--cores"), type="integer", default=1, help="Number of cores to use in phewas()", metavar="integer")
)

opt_parser = OptionParser(option_list=option_list)
opt = parse_args(opt_parser)

if (is.null(opt$geno_cov) || is.null(opt$icd10) || is.null(opt$outdir)) {
  print_help(opt_parser)
  stop("Missing required arguments: --geno_cov, --icd10 and --outdir", call.=FALSE)
}

dir.create(opt$outdir, showWarnings = FALSE, recursive = TRUE)

cat("[", format(Sys.time(), "%H:%M:%S"), "] Loading genotypes and covariates...\n")
geno_cov <- read.table(opt$geno_cov, header=TRUE, sep="\t", stringsAsFactors=FALSE)

cat("[", format(Sys.time(), "%H:%M:%S"), "] Loading and processing ICD10 data...\n")
icd10_data <- read.csv(opt$icd10)

# Build the id.sex data frame to filter sex-specific false positives (UKB: 1=male, 0=female -> M/F)
id_sex <- data.frame(
  id = geno_cov$id,
  sex = ifelse(geno_cov$Sex == 1, "M", "F")
)

# Ensure icd10_data has the four columns required by the package (id, vocabulary_id, code, count)
if(!"count" %in% names(icd10_data)) {
  icd10_data$count <- 1
}
icd10_data <- icd10_data[, c("id", "vocabulary_id", "code", "count")]

# Convert ICD10 records to a PheCode matrix
cat("[", format(Sys.time(), "%H:%M:%S"), "] Creating Phenotype Matrix (min.code.count=2)...\n")
phenotypes <- createPhenotypes(
  id.vocab.code.index = icd10_data,
  min.code.count = 2,
  add.phecode.exclusions = TRUE,
  translate = TRUE,
  vocabulary.map = phecode_map_icd10,
  id.sex = id_sex
)

# Determine all SV column names (excluding id, Age, Sex and PC1-10)
exclude_cols <- c("id", "Age", "Sex", paste0("PC", 1:10))
sv_cols <- setdiff(colnames(geno_cov), exclude_cols)

# If sv_index is provided, process only that specific SV
if (!is.null(opt$sv_index)) {
  if (opt$sv_index < 1 || opt$sv_index > length(sv_cols)) {
    stop("sv_index is out of range.", call.=FALSE)
  }
  sv_cols <- sv_cols[opt$sv_index]
}

cat("  Number of SVs to process:", length(sv_cols), "\n")

all_results <- list()

for (i in seq_along(sv_cols)) {
  sv <- sv_cols[i]
  cat(sprintf("\n[%s] === Testing SV: %s ===\n", format(Sys.time(), "%H:%M:%S"), sv))
  
  # Prepare predictors for the current SV, retaining only required columns to reduce memory use
  current_cov <- geno_cov[, c("id", sv, "Age", "Sex", paste0("PC", 1:10))]
  # Ensure the SV column is numeric (additive coding)
  current_cov[[sv]] <- as.numeric(current_cov[[sv]])
  
  # Merge phenotypes and covariates
  current_data <- inner_join(phenotypes, current_cov, by="id")
  
  # Run PheWAS
  results <- tryCatch({
    phewas(phenotypes = names(phenotypes)[-1], 
           genotypes = sv, 
           data = current_data,
           covariates = c("Age", "Sex", paste0("PC", 1:10)),
           significance.threshold = c("bonferroni", "fdr"),
           cores = opt$cores) # Limit the number of cores
  }, error = function(e) {
    cat("  Error during execution:", e$message, "\n")
    return(NULL)
  })
  
  if (!is.null(results)) {
    # Add human-readable phenotype labels using the PheWAS dictionary mapping
    results <- addPhecodeInfo(results)
    
    # Add the SV identifier column
    results$SV_ID <- sv
    
    # Sort by p-value
    results <- results[order(results$p),]
    
    all_results[[sv]] <- results
    
    # Count nominally significant results
    sig_count <- sum(results$p < 0.05, na.rm = TRUE) 
    cat("  Completed; found", sig_count, "potential associations (uncorrected P < 0.05).\n")
  }
}

cat("\n[", format(Sys.time(), "%H:%M:%S"), "] Combining and saving results...\n")
if (length(all_results) > 0) {
  final_results <- bind_rows(all_results)
  
  # Apply p-value-based multiple-testing correction within each SV
  # The PheWAS package reports FDR/Bonferroni thresholds; here we additionally compute exact p.adjust values
  final_results <- final_results %>%
    group_by(SV_ID) %>%
    mutate(
      FDR_P = p.adjust(p, method = "fdr"),
      Bonferroni_P = p.adjust(p, method = "bonferroni")
    ) %>%
    ungroup() %>%
    arrange(p)
  
  # If running in single-SV mode for array jobs, write an indexed output file
  if (!is.null(opt$sv_index)) {
    out_path <- file.path(opt$outdir, paste0("sv_phewas_results_part_", opt$sv_index, ".tsv"))
  } else {
    out_path <- file.path(opt$outdir, "sv_phewas_results_all.tsv")
  }
  
  write.table(final_results, out_path, sep="\t", quote=FALSE, row.names=FALSE)
  cat("  → Saved:", out_path, "\n")
  
} else {
  cat("  No results were generated.\n")
}

cat("[", format(Sys.time(), "%H:%M:%S"), "] PheWAS completed!\n")

