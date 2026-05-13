# 读取数据
df <- read.table("136.autosome.merge.stat.tsv",
                 header = TRUE,
                 sep = "\t",
                 stringsAsFactors = FALSE)

############################
# function of P-value
############################

format_p <- function(p) {
  sprintf("%.3E", p)
}

############################
# 1. ngs_ser vs assem_ser
############################

cat("\n===== ngs_ser vs assem_ser =====\n")

diff_ser <- df$ngs_ser - df$assem_ser

# test of normal distribution
shapiro_ser <- shapiro.test(diff_ser)

cat("Shapiro-Wilk P-value:",
    format_p(shapiro_ser$p.value), "\n")

# choose test
if (shapiro_ser$p.value > 0.05) {

  cat("Use paired t-test\n")

  test_ser <- t.test(df$ngs_ser,
                     df$assem_ser,
                     paired = TRUE)

} else {

  cat("Use Wilcoxon signed-rank test\n")

  test_ser <- wilcox.test(df$ngs_ser,
                          df$assem_ser,
                          paired = TRUE)
}

cat("Test P-value:",
    format_p(test_ser$p.value), "\n")


############################
# 2. ngs_fer vs assem_fer
############################

cat("\n===== ngs_fer vs assem_fer =====\n")

diff_fer <- df$ngs_fer - df$assem_fer

# test of normal distribution
shapiro_fer <- shapiro.test(diff_fer)

cat("Shapiro-Wilk P-value:",
    format_p(shapiro_fer$p.value), "\n")

# choose test
if (shapiro_fer$p.value > 0.05) {

  cat("Use paired t-test\n")

  test_fer <- t.test(df$ngs_fer,
                     df$assem_fer,
                     paired = TRUE)

} else {

  cat("Use Wilcoxon signed-rank test\n")

  test_fer <- wilcox.test(df$ngs_fer,
                          df$assem_fer,
                          paired = TRUE)
}

cat("Test P-value:",
    format_p(test_fer$p.value), "\n")
