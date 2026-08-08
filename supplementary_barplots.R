# ==============================================================================
# Figure 4 - Top Genes by Variant Count, per Cancer Type
# ==============================================================================
#
# Description:
#   For each cancer type, plots a side-by-side (dodged) bar chart of the top
#   N genes by variant count, split by germline confirmation status
#   (Yes/No), with each bar's share of the panel's total (%) and raw count
#   labeled above it. One PDF is produced per cancer type.
#
# ==============================================================================

# ---- 1. Setup ----------------------------------------------------------------

library(dplyr)
library(tidyr)
library(readxl)
library(ggplot2)

# ---- Paths (EDIT THESE before running) -----------------------------------
data_dir   <- "data"
output_dir <- "output"

dir.create(output_dir, showWarnings = FALSE, recursive = TRUE)

f_gene_data <- file.path(data_dir, "input_file.xlsx")

# Number of top genes (by variant count) to show per cancer type panel
n_top_genes <- 10


# ==============================================================================
# 2. Load data and keep only variants that completed germline testing
# ==============================================================================

gene_data <- read_excel(f_gene_data)


gene_data <- gene_data %>% filter(GERM_TEST == "Yes, test performed")

gene_data$NEOPLASM <- toupper(gene_data$NEOPLASM)

message("Cancer types present in the input data:")
print(unique(gene_data$NEOPLASM))

all_genes    <- unique(gene_data$GENE)
cancer_types <- unique(gene_data$NEOPLASM)


# ==============================================================================
# 3. Build one plot per cancer type
# ==============================================================================

for (cancer_type in cancer_types) {

  message("Processing: ", cancer_type)

  variants_this_type <- gene_data %>%
    filter(NEOPLASM == cancer_type)

  # Variant count per gene, split by germline confirmation status
  counts_by_gene <- variants_this_type %>%
    group_by(GENE, GERMLINE) %>%
    summarise(count = n(), .groups = "drop") %>%
    mutate(GERMLINE = ifelse(is.na(GERMLINE), "No", GERMLINE))

  # Add explicit zero-count rows for panel genes with no variants at all in
  # this cancer type, so the top-N ranking considers the full gene panel
  missing_genes <- setdiff(all_genes, counts_by_gene$GENE)
  if (length(missing_genes) > 0) {
    counts_by_gene <- bind_rows(
      counts_by_gene,
      data.frame(GENE = missing_genes, GERMLINE = "No", count = 0)
    )
  }

  # Total variants per gene, used to rank and select the top N genes
  gene_totals <- counts_by_gene %>%
    group_by(GENE) %>%
    summarise(total_count = sum(count), .groups = "drop") %>%
    arrange(desc(total_count))

  # Total variants across ALL genes in this cancer type (shown in the panel
  # title annotation) - computed before restricting to the top N genes
  total_variants_all_genes <- sum(gene_totals$total_count)

  gene_totals <- gene_totals %>% slice_head(n = n_top_genes)

  counts_by_gene <- counts_by_gene %>%
    filter(GENE %in% gene_totals$GENE) %>%
    left_join(gene_totals, by = "GENE")

  # y-axis headroom, scaled to leave room for the percentage/count labels
  y_axis_limit <- if (max(gene_totals$total_count) > 20) {
    max(gene_totals$total_count) + 15
  } else {
    max(gene_totals$total_count) + 2.5
  }

  # Ensure every one of the top N genes has both a "Yes" and "No" row (even
  # if one status has zero variants), and compute each bar's percentage
  # share of the top-N-gene total - NOTE this is relative to the top N genes
  # shown, not the cancer type's overall total (that's the separate
  # "N = ..." panel annotation below, based on total_variants_all_genes).
  top_n_total <- sum(gene_totals$total_count)

  counts_by_gene <- counts_by_gene %>%
    tidyr::complete(GENE, GERMLINE = c("Yes", "No"),
                     fill = list(count = 0, total_count = 0)) %>%
    mutate(
      percentage_positive = ifelse(GERMLINE == "Yes", (count / top_n_total) * 100, 0),
      percentage_negative = ifelse(GERMLINE == "No", (count / top_n_total) * 100, 0),
      GENE     = factor(GENE, levels = gene_totals$GENE),
      GERMLINE = factor(GERMLINE, levels = c("Yes", "No"))
    )

  # ---- Plot: side-by-side (dodged) bars, Yes vs No ----
  p <- ggplot(counts_by_gene, aes(x = GENE, y = count, fill = GERMLINE)) +
    geom_bar(stat = "identity", position = position_dodge(width = 0.9), color = "black") +
    # Labels for bars with at least one variant
    geom_text(
      data = counts_by_gene %>% filter(GERMLINE == "Yes" & count > 0),
      aes(x = GENE, y = count, color = "brown3",
          label = paste0(round(percentage_positive, 2), "% - (", count, ")")),
      position = position_dodge(width = 0.9), vjust = -1, hjust = -0.1,
      size = 7, angle = 90, fontface = "bold"
    ) +
    geom_text(
      data = counts_by_gene %>% filter(GERMLINE == "No" & count > 0),
      aes(x = GENE, y = count,
          label = paste0(round(percentage_negative, 2), "% - (", count, ")")),
      position = position_dodge(width = 0.9), vjust = 2, hjust = -0.1,
      size = 7, angle = 90
    ) +
    # Labels for genes with zero variants of a given status
    geom_text(
      data = counts_by_gene %>% filter(GERMLINE == "Yes" & count == 0),
      aes(x = GENE, y = count, color = "brown3", label = "0.00% - (0)"),
      position = position_dodge(width = 0.9), vjust = -1, hjust = -0.1,
      size = 7, angle = 90, fontface = "bold"
    ) +
    geom_text(
      data = counts_by_gene %>% filter(GERMLINE == "No" & count == 0),
      aes(x = GENE, y = count, label = "0.00% - (0)"),
      position = position_dodge(width = 0.9), vjust = 2, hjust = -0.1,
      size = 7, angle = 90
    ) +
    # Panel title: cancer type + total variant count across ALL genes (not
    # just the top N shown)
    annotate(
      "text",
      x = 0.5 + length(gene_totals$GENE) / 2,
      y = y_axis_limit - (y_axis_limit * 0.05),
      label = paste0(cancer_type, "\nN = ", total_variants_all_genes),
      size = 6.5, fontface = "bold", hjust = 0.5
    ) +
    scale_y_continuous(limits = c(0, y_axis_limit)) +
    scale_fill_manual(values = c("Yes" = "brown3", "No" = "lightsteelblue3")) +
    scale_color_manual(values = c("brown3" = "brown3"), guide = "none") +
    labs(x = "Gene", y = "Number of Variants", fill = "Test result") +
    theme_minimal() +
    theme(
      axis.text.x    = element_text(angle = 90, vjust = 0.5, hjust = 0.5, size = 18, color = "black"),
      axis.text.y    = element_text(size = 16, color = "black"),
      axis.title.x   = element_blank(),
      axis.title.y   = element_text(size = 18),
      text           = element_text(size = 18),
      panel.grid     = element_blank(),
      panel.border   = element_rect(color = "black", fill = NA, linewidth = 0.5),
      plot.title     = element_text(hjust = 0.5, face = "bold", size = 20),
      legend.position = "none"
    )

  pdf_path <- file.path(output_dir, paste0(tolower(cancer_type), "_neoplasm.pdf"))
  pdf(file = pdf_path, width = 12, height = 8)
  print(p)
  dev.off()
  message("Saved: ", pdf_path)
}
