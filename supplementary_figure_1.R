# ==============================================================================
# Figure 2 - Confirmed Germline Variants by Gene
# ==============================================================================
#
# Description:
#   Stacked bar chart showing the number of confirmed germline variants per
#   gene, stacked by variant classification (missense, nonsense, frameshift,
#   splice site, etc.), with the total variant count and percentage of the
#   cohort's overall variant burden labeled above each bar. Genes are ordered
#   by descending total variant count.

#
# Output file (written to `output_dir`):
#   - figure2_germline_variants_by_gene.pdf
#
# Required packages: dplyr, ggplot2
#
# ==============================================================================

# ---- 1. Setup ----------------------------------------------------------------

library(dplyr)
library(ggplot2)

# ---- Paths (EDIT THESE before running) -----------------------------------
data_dir   <- "data"
output_dir <- "output_dir"

dir.create(output_dir, showWarnings = FALSE, recursive = TRUE)

f_mutation_data <- file.path(data_dir, "d_mut_germ.txt")
f_pdf_out       <- file.path(output_dir, "figure2_germline_variants_by_gene.pdf")


# ==============================================================================
# 2. Load confirmed germline variants
# ==============================================================================

mutation_data <- read.delim(f_mutation_data, header = TRUE, sep = "\t", stringsAsFactors = FALSE)

confirmed_germline_variants <- mutation_data %>% filter(GERMLINE == "Yes")
message("Confirmed germline variants: ", nrow(confirmed_germline_variants))


# ==============================================================================
# 3. Count variants per gene and variant classification
# ==============================================================================

counts_by_gene <- confirmed_germline_variants %>%
     group_by(Hugo_Symbol, Variant_Classification) %>%
     summarise(count = n(), .groups = "drop")

# Total variants per gene - used both to order genes on the x-axis (most to
# least variants) and as the basis for the percentage labels in Section 4
gene_totals <- counts_by_gene %>%
     group_by(Hugo_Symbol) %>%
     summarise(total_count = sum(count)) %>%
     arrange(desc(total_count))

counts_by_gene <- counts_by_gene %>%
     left_join(gene_totals, by = "Hugo_Symbol") %>%
     mutate(Hugo_Symbol = factor(Hugo_Symbol, levels = gene_totals$Hugo_Symbol))


# ==============================================================================
# 4. Compute per-gene percentage labels
# ==============================================================================

total_variants <- sum(counts_by_gene$count)

gene_percentages <- gene_totals %>%
     mutate(
          percentage = round((total_count / total_variants) * 100, 1),
          label      = paste0(percentage, "% (", total_count, ")")
     )

# Standardize missing variant classifications for a consistent legend entry
message("Variant classifications before standardization:")
print(unique(counts_by_gene$Variant_Classification))

counts_by_gene <- counts_by_gene %>%
     mutate(Variant_Classification = ifelse(is.na(Variant_Classification), "Unknown", Variant_Classification))

message("Variant classifications after standardization:")
print(unique(counts_by_gene$Variant_Classification))


# ==============================================================================
# 5. Build and save the plot
# ==============================================================================

p <- ggplot(counts_by_gene, aes(x = Hugo_Symbol, y = count, fill = Variant_Classification)) +
     geom_bar(stat = "identity", position = "stack", color = "black") +
     # Percentage + count label above each gene's total bar
     geom_text(
          data = gene_percentages,
          aes(x = Hugo_Symbol, y = total_count, label = label),
          inherit.aes = FALSE, vjust = 0.5, hjust = -0.1, size = 6, fontface = "bold", angle = 90
     ) +
     # NOTE: fixed at 0-125 to fit the original dataset - adjust if any gene's
     # total variant count now exceeds this, or bars will be clipped.
     scale_y_continuous(limits = c(0, 125)) +
     scale_fill_manual(values = c(
          "Missense_Mutation" = "darkcyan",
          "Nonsense_Mutation"  = "tomato",
          "Frame_Shift_Del"    = "turquoise",
          "Frame_Shift_Ins"    = "gold",
          "Splice_Site"        = "orange",
          "In_Frame_Del"       = "seagreen3",
          "Splice_Region"      = "deepskyblue2",
          "In_Frame_Ins"       = "mediumpurple1",
          "Unknown"            = "red"
     )) +
     labs(x = "Gene", y = "Number of Variants", fill = "Variant Type") +
     theme_minimal() +
     theme(
          axis.text.x        = element_text(angle = 90, vjust = 0.5, hjust = 1, size = 15),
          text                = element_text(size = 15),
          panel.grid.major.x  = element_blank(),
          panel.grid.minor.x  = element_blank(),
          axis.title.x        = element_text(size = 14),
          axis.title.y        = element_text(size = 14)
     ) +
     ggtitle(paste("Total No. of germline variants:", total_variants))

pdf(f_pdf_out, height = 10, width = 15)
print(p)
dev.off()
message("PDF saved to: ", f_pdf_out)
