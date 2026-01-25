#!/usr/bin/env Rscript

suppressPackageStartupMessages({
  library(tidyverse)
  library(data.table)
  library(pheatmap)
})

args <- commandArgs(trailingOnly = TRUE)
project_dir <- if (length(args) > 0) normalizePath(args[1]) else normalizePath(".")

metrics_long <- file.path(project_dir, "results", "metrics_long.tsv")
runtime_tsv  <- file.path(project_dir, "results", "runtime.tsv")
bench_dir    <- file.path(project_dir, "results", "benchmarks")
plot_dir     <- file.path(project_dir, "results", "plots")

if (!file.exists(metrics_long)) stop("Missing results/metrics_long.tsv. Run hap.py and metrics export first.")
if (!file.exists(runtime_tsv))  stop("Missing results/runtime.tsv. Run runtime export first.")

dir.create(plot_dir, showWarnings = FALSE, recursive = TRUE)

caller_levels <- c("gatk", "deepvariant", "strelka2", "freebayes")

# -------------------------
# 1) Load summary metrics (SNP/INDEL, ALL)
# -------------------------
m <- fread(metrics_long) %>%
  as_tibble() %>%
  mutate(
    Type = factor(Type, levels = c("SNP", "INDEL")),
    caller = factor(caller, levels = caller_levels)
  )

# -------------------------
# 2) Barplots: Precision/Recall/F1 (separate SNP vs INDEL)
# -------------------------
m_long <- m %>%
  select(
    caller,
    Type,
    recall = `METRIC.Recall`,
    precision = `METRIC.Precision`,
    f1 = `METRIC.F1_Score`
  ) %>%
  pivot_longer(
    cols = c(recall, precision, f1),
    names_to = "metric",
    values_to = "value"
  )

p_bar <- ggplot(m_long, aes(x = caller, y = value, fill = metric)) +
  geom_col(position = position_dodge(width = 0.8), width = 0.7) +
  facet_wrap(~Type, nrow = 1) +
  scale_y_continuous(limits = c(0, 1)) +
  labs(
    title = "Benchmark (ALL): Precision / Recall / F1 theo SNP vs INDEL",
    x = "Variant caller",
    y = "Gia tri"
  ) +
  theme_bw() +
  theme(axis.text.x = element_text(angle = 20, hjust = 1))

ggsave(file.path(plot_dir, "01_bar_precision_recall_f1.png"), p_bar, width = 12, height = 4, dpi = 200)

# -------------------------
# 3) Heatmap: F1 only (SNP/INDEL)
# -------------------------
hm <- m %>%
  select(caller, Type, f1 = `METRIC.F1_Score`) %>%
  pivot_wider(names_from = Type, values_from = f1) %>%
  arrange(caller)

hm_mat <- hm %>%
  select(-caller) %>%
  as.matrix()

rownames(hm_mat) <- hm$caller

png(file.path(plot_dir, "02_heatmap_f1.png"), width = 900, height = 450)
pheatmap(
  hm_mat,
  cluster_rows = FALSE,
  cluster_cols = FALSE,
  main = "F1 heatmap (ALL) - SNP vs INDEL"
)
dev.off()

# -------------------------
# 4) Runtime plot (wall + max RSS)
# -------------------------
rt <- fread(runtime_tsv) %>%
  as_tibble() %>%
  mutate(caller = factor(caller, levels = caller_levels))

p_time <- ggplot(rt, aes(x = caller, y = wall_s, fill = step)) +
  geom_col(position = position_dodge(width = 0.8), width = 0.7) +
  labs(title = "Runtime (wall seconds) theo buoc", x = "Caller", y = "Wall time (s)") +
  theme_bw()

ggsave(file.path(plot_dir, "03_runtime_wall.png"), p_time, width = 10, height = 4, dpi = 200)

p_mem <- ggplot(rt, aes(x = caller, y = max_rss_kb / 1024, fill = step)) +
  geom_col(position = position_dodge(width = 0.8), width = 0.7) +
  labs(title = "Max RSS theo buoc", x = "Caller", y = "Max RSS (MB)") +
  theme_bw()

ggsave(file.path(plot_dir, "04_runtime_maxrss.png"), p_mem, width = 10, height = 4, dpi = 200)

# -------------------------
# 5) PR curves from hap.py roc outputs
#    roc.*.csv.gz contains METRIC.Recall & METRIC.Precision over QUAL thresholds
# -------------------------
roc_files <- list.files(bench_dir, pattern = "\\.roc\\..*\\.csv\\.gz$", recursive = TRUE, full.names = TRUE)

read_roc <- function(f) {
  x <- suppressMessages(readr::read_csv(f, show_col_types = FALSE))
  caller <- basename(dirname(dirname(f)))
  x %>% mutate(caller = caller, file = basename(f))
}

if (length(roc_files) > 0) {
  roc_df <- map_dfr(roc_files, read_roc) %>%
    filter(Filter %in% c("PASS", "ALL")) %>%
    mutate(caller = factor(caller, levels = caller_levels))

  roc_pick <- roc_df %>%
    filter(str_detect(file, "Locations")) %>%
    filter((Type == "SNP" & Filter == "PASS") | (Type == "INDEL" & Filter == "PASS")) %>%
    filter(is.finite(`METRIC.Recall`), is.finite(`METRIC.Precision`))

  plot_pr <- function(dat, vtype) {
    ggplot(dat %>% filter(Type == vtype), aes(x = `METRIC.Recall`, y = `METRIC.Precision`, color = caller)) +
      geom_line(linewidth = 1) +
      labs(
        title = paste0("Precision-Recall curve (hap.py, PASS): ", vtype),
        x = "Recall",
        y = "Precision"
      ) +
      theme_bw()
  }

  p_pr_snp <- plot_pr(roc_pick, "SNP")
  p_pr_indel <- plot_pr(roc_pick, "INDEL")

  ggsave(file.path(plot_dir, "05_pr_curve_snp.png"), p_pr_snp, width = 7, height = 5, dpi = 200)
  ggsave(file.path(plot_dir, "06_pr_curve_indel.png"), p_pr_indel, width = 7, height = 5, dpi = 200)
}

# -------------------------
# 6) Per-variant QUAL distributions (TP vs FP) from exported per_variant.tsv.gz
# -------------------------
per_files <- list.files(bench_dir, pattern = "per_variant\\.tsv\\.gz$", recursive = TRUE, full.names = TRUE)

if (length(per_files) > 0) {
  per_df <- map_dfr(per_files, function(f) {
    d <- fread(f) %>% as_tibble()
    colnames(d) <- c("caller", "label", "Type", "QUAL", "CHROM", "POS")
    d$caller <- basename(dirname(dirname(f)))
    d
  }) %>%
    mutate(
      caller = factor(caller, levels = caller_levels),
      label = factor(label, levels = c(1, 0), labels = c("TP", "FP")),
      Type = factor(Type, levels = c("SNP", "INDEL")),
      QUAL = as.numeric(QUAL)
    )

  p_qual <- ggplot(per_df %>% filter(is.finite(QUAL)), aes(x = caller, y = QUAL, fill = label)) +
    geom_boxplot(outlier.size = 0.4) +
    facet_wrap(~Type, scales = "free_y", nrow = 1) +
    labs(
      title = "QUAL distribution (TP vs FP) theo caller (tu hap.py classification VCF)",
      x = "Caller",
      y = "QUAL"
    ) +
    theme_bw() +
    theme(axis.text.x = element_text(angle = 20, hjust = 1))

  ggsave(file.path(plot_dir, "07_qual_tp_fp_boxplot.png"), p_qual, width = 12, height = 4, dpi = 200)
}

cat("[DONE] Plots written to results/plots/\n")
