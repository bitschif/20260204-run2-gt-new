#!/usr/bin/env Rscript

suppressPackageStartupMessages({
  library(tidyverse)
  library(data.table)
  library(pheatmap)
  library(jsonlite)
})

args <- commandArgs(trailingOnly = TRUE)
project_dir <- if (length(args) > 0) normalizePath(args[1]) else normalizePath(".")

metrics_long <- file.path(project_dir, "results", "metrics_long.tsv")
runtime_tsv  <- file.path(project_dir, "results", "runtime.tsv")
runtime_csv  <- file.path(project_dir, "logs", "runtime.csv")
bench_dir    <- file.path(project_dir, "results", "benchmarks")
plot_dir     <- file.path(project_dir, "results", "plots")
metrics_combined <- file.path(project_dir, "results", "metrics_combined.tsv")

if (!file.exists(metrics_long)) {
  metrics_json <- list.files(
    bench_dir,
    pattern = "\\.metrics\\.json\\.gz$",
    recursive = TRUE,
    full.names = TRUE
  )
  summary_csv <- list.files(
    bench_dir,
    pattern = "\\.summary\\.csv$",
    recursive = TRUE,
    full.names = TRUE
  )
  metric_files <- list.files(bench_dir, pattern = "_metrics\\.tsv$", recursive = TRUE, full.names = TRUE)

  if (length(metrics_json) == 0 && length(summary_csv) == 0 && length(metric_files) == 0) {
    stop("Missing results/metrics_long.tsv and no hap.py metrics files found under results/benchmarks.")
  }

  get_caller <- function(path) {
    parent <- dirname(path)
    if (basename(parent) == "happy") {
      return(basename(dirname(parent)))
    }
    basename(parent)
  }

  read_metrics_json <- function(path) {
    raw <- jsonlite::fromJSON(gzfile(path), flatten = TRUE)
    metrics <- raw
    if ("metrics" %in% names(raw)) {
      metrics <- raw$metrics
    } else if ("summary" %in% names(raw)) {
      metrics <- raw$summary
    }
    metrics_tbl <- as_tibble(metrics)
    if (!"Type" %in% names(metrics_tbl) && "type" %in% names(metrics_tbl)) {
      metrics_tbl <- metrics_tbl %>% rename(Type = type)
    }
    caller <- get_caller(path)
    metrics_tbl %>% mutate(caller = caller)
  }

  read_summary_csv <- function(path) {
    df <- suppressMessages(readr::read_csv(path, show_col_types = FALSE))
    caller <- get_caller(path)
    df %>% mutate(caller = caller)
  }

  metrics_raw <- list()
  if (length(metrics_json) > 0) {
    metrics_raw <- append(metrics_raw, list(map_dfr(metrics_json, read_metrics_json)))
  }
  if (length(summary_csv) > 0) {
    metrics_raw <- append(metrics_raw, list(map_dfr(summary_csv, read_summary_csv)))
  }
  if (length(metric_files) > 0) {
    metrics_raw <- append(metrics_raw, list(map_dfr(metric_files, ~fread(.x) %>% as_tibble())))
  }

  metrics_long_data <- bind_rows(metrics_raw)

  if (nrow(metrics_long_data) == 0) {
    stop("No metrics data found after loading available hap.py metrics.")
  }

  if (!is.null(metrics_long_data) && nrow(metrics_long_data) > 0) {
    fwrite(metrics_long_data, metrics_combined, sep = "\t")
  }

  metrics_long_data <- metrics_long_data %>%
    rename(
      Type = any_of(c("VariantType", "Type")),
      `METRIC.Precision` = any_of(c("Precision", "METRIC.Precision")),
      `METRIC.Recall` = any_of(c("Recall", "METRIC.Recall")),
      `METRIC.F1_Score` = any_of(c("F1", "METRIC.F1_Score", "F1_Score")),
      `METRIC.ROC_AUC` = any_of(c("ROC", "ROC_AUC", "METRIC.ROC_AUC", "METRIC.RocAuc"))
    ) %>%
    select(
      caller,
      any_of("Type"),
      any_of(c("METRIC.Precision", "METRIC.Recall", "METRIC.F1_Score", "METRIC.ROC_AUC"))
    )

  if (!"Type" %in% names(metrics_long_data)) {
    metrics_long_data$Type <- NA_character_
  }

  for (col in c("METRIC.Precision", "METRIC.Recall", "METRIC.F1_Score", "METRIC.ROC_AUC")) {
    if (!col %in% names(metrics_long_data)) {
      metrics_long_data[[col]] <- NA_character_
    }
  }

  coerce_atomic <- function(x) {
    if (is.list(x)) {
      vapply(x, function(value) {
        if (length(value) == 0) {
          NA_character_
        } else {
          as.character(value[[1]])
        }
      }, character(1))
    } else {
      x
    }
  }

  mean_or_na <- function(x) {
    if (all(is.na(x))) {
      NA_real_
    } else {
      mean(x, na.rm = TRUE)
    }
  }

  metrics_long_data <- metrics_long_data %>%
    mutate(across(everything(), coerce_atomic)) %>%
    mutate(
      `METRIC.Precision` = as.numeric(`METRIC.Precision`),
      `METRIC.Recall` = as.numeric(`METRIC.Recall`),
      `METRIC.F1_Score` = as.numeric(`METRIC.F1_Score`),
      `METRIC.ROC_AUC` = as.numeric(`METRIC.ROC_AUC`)
    ) %>%
    mutate(
      caller = tolower(caller),
      Type = toupper(Type)
    ) %>%
    group_by(caller, Type) %>%
    summarise(
      `METRIC.Precision` = mean_or_na(`METRIC.Precision`),
      `METRIC.Recall` = mean_or_na(`METRIC.Recall`),
      `METRIC.F1_Score` = mean_or_na(`METRIC.F1_Score`),
      `METRIC.ROC_AUC` = mean_or_na(`METRIC.ROC_AUC`),
      .groups = "drop"
    )

  fwrite(metrics_long_data, metrics_long, sep = "\t")
}

has_runtime_csv <- file.exists(runtime_csv)
has_runtime_tsv <- file.exists(runtime_tsv)
if (!has_runtime_csv && !has_runtime_tsv) {
  warning("Missing logs/runtime.csv and results/runtime.tsv. Skipping runtime plots.")
}

dir.create(plot_dir, showWarnings = FALSE, recursive = TRUE)

caller_levels <- c("gatk", "deepvariant", "strelka2", "freebayes")

# -------------------------
# 1) Load summary metrics (SNP/INDEL, ALL)
# -------------------------
metric_cols <- c("METRIC.Recall", "METRIC.Precision", "METRIC.F1_Score", "METRIC.ROC_AUC")
m_raw <- fread(metrics_long) %>%
  as_tibble()

for (col in metric_cols) {
  if (!col %in% names(m_raw)) {
    m_raw[[col]] <- NA_real_
  }
}

mean_or_na <- function(x) {
  x <- as.numeric(x)
  if (all(is.na(x))) {
    NA_real_
  } else {
    mean(x, na.rm = TRUE)
  }
}

m <- m_raw %>%
  mutate(
    Type = factor(Type, levels = c("SNP", "INDEL")),
    caller = factor(caller, levels = caller_levels)
  ) %>%
  group_by(caller, Type) %>%
  summarise(
    across(all_of(metric_cols), mean_or_na),
    .groups = "drop"
  ) %>%
  filter(if_any(all_of(metric_cols), is.finite))

# -------------------------
# 2) Barplots: Precision/Recall/F1 (separate SNP vs INDEL)
# -------------------------
m_long <- m %>%
  select(
    caller,
    Type,
    recall = `METRIC.Recall`,
    precision = `METRIC.Precision`,
    f1 = `METRIC.F1_Score`,
    roc = `METRIC.ROC_AUC`
  ) %>%
  pivot_longer(
    cols = c(recall, precision, f1, roc),
    names_to = "metric",
    values_to = "value"
  ) %>%
  filter(is.finite(value))

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

ggsave(file.path(plot_dir, "01_bar_precision_recall_f1_roc.png"), p_bar, width = 12, height = 4, dpi = 200)

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
if (has_runtime_csv) {
  rt <- suppressMessages(readr::read_csv(runtime_csv, col_names = c("caller", "time_s"), show_col_types = FALSE)) %>%
    mutate(caller = factor(caller, levels = caller_levels))

  p_time <- ggplot(rt %>% filter(is.finite(time_s)), aes(x = caller, y = time_s)) +
    geom_col(width = 0.7, fill = "#4C78A8") +
    labs(title = "Runtime (seconds)", x = "Caller", y = "Time (s)") +
    theme_bw()

  ggsave(file.path(plot_dir, "03_runtime_time.png"), p_time, width = 10, height = 4, dpi = 200)
} else if (has_runtime_tsv) {
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
}

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
