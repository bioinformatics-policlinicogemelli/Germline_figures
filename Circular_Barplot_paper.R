# ==============================================================================
# Circular Barplot - Germline Variant Frequency by Gene and Cancer Type
# ==============================================================================
#
# Description:
#   For each selected cancer type, plots the percentage of germline variants
#   attributable to each gene as a circular (polar) bar chart. Bars are
#   colored by gene; green/red arcs mark variants classified as ON/OFF tumor
#   status. Output is saved as both a vector PDF and a transparent SVG (for
#   easy editing in PowerPoint/Illustrator).
#
# Input file:
#   - germline_mutations_raw.xlsx   Same source file used by Oncoplot.R
#
#   NOTE: file name above is a placeholder. Either rename your local file to
#   match, or edit the path in the "Paths" section below.
#
# Output files (written to `output_dir`):
#   - circular_barplot_paper.pdf
#   - circular_barplot_paper.svg
#
# Required packages: dplyr, readxl, ggplot2, svglite
#
# ==============================================================================

# ---- 1. Setup ----------------------------------------------------------------

library(dplyr)
library(readxl)
library(ggplot2)
library(svglite)

# ---- Paths (EDIT THESE before running) ----------------------------------
data_dir   <- "/Users/simonerossi/Desktop/GERMLINE/REVISIONI"
output_dir <- "/Users/simonerossi/Desktop/GERMLINE/REVISIONI"

dir.create(output_dir, showWarnings = FALSE, recursive = TRUE)

f_raw     <- file.path(data_dir, "d_mut_germ_on_off_tumor.xlsx")
f_pdf_out <- file.path(output_dir, "circular_barplot_paper.pdf")
f_svg_out <- file.path(output_dir, "circular_barplot_paper.svg")


# ==============================================================================
# 2. Load data and select germline variants
# ==============================================================================

variants_tested <- read_excel(f_raw)
variants_tested$NEOPLASM <- toupper(variants_tested$NEOPLASM)

message("Cancer types present in the input data:")
print(unique(variants_tested$NEOPLASM))

# Cancer types to include in the plot
cancer_types_to_plot <- c("OVARY", "ENDOMETRIUM", "LUNG", "COLORECTUM", "PANCREAS")

# Keep only variants with a positive germline test
germline_variants <- variants_tested[variants_tested$GERMLINE == "Yes", ]
message("Germline variant rows retained: ", nrow(germline_variants))


# ==============================================================================
# 3. Gene color palette
# ==============================================================================

gene_names <- c("ATM", "BARD1", "BRCA1", "BRCA2", "BRIP1", "CHEK2", "MLH1", "MSH2", "MSH6", "MUTYH",
                "PALB2", "PMS2", "POLE", "RAD51C", "RAD51D", "RET", "SDHA", "SDHB", "PTEN", "FH")

gene_colors <- c("#D1B0D9", "#8DD3C7", "#9BC8E6", "#377EB8", "darkorchid", "#FD8D3C", "forestgreen", "#FFE2A9",
                 "#F4A582", "#A6D96A", "maroon", "#DE77AE", "darkblue", "antiquewhite2", "red", "pink",
                 "orange", "royalblue", "darkgray", "darkgoldenrod")

gene_color_lookup <- setNames(gene_colors, gene_names)

message("Gene -> color mapping:")
print(data.frame(gene = gene_names, color = gene_colors))

# Attach each variant's gene color as its own column
germline_variants$gene_color <- gene_color_lookup[germline_variants$Hugo_Symbol]


# ==============================================================================
# 4. Filter to cancer types of interest and compute per-gene percentages
# ==============================================================================

germline_variants_filtered <- germline_variants %>%
     filter(NEOPLASM %in% cancer_types_to_plot)

# Count variants per gene, split by ON/OFF tumor status
gene_counts <- germline_variants_filtered %>%
     group_by(NEOPLASM, Hugo_Symbol, TYPE) %>%
     summarise(Count = n(), .groups = "drop")

# Percentage of each gene's variants relative to the cancer type's total
totals <- gene_counts %>%
     group_by(NEOPLASM) %>%
     summarise(Total = sum(Count))

gene_counts <- gene_counts %>%
     left_join(totals, by = "NEOPLASM") %>%
     mutate(Percentage = round((Count / Total) * 100, 2))

# Re-attach gene colors (a direct lookup, since color depends only on gene
# name) and build the bar label
gene_counts <- gene_counts %>%
     mutate(
          gene_color = gene_color_lookup[Hugo_Symbol],
          Category   = paste0(Hugo_Symbol, " - ", Percentage, "%")
     )


# ==============================================================================
# 5. Assemble the plotting data frame
# ==============================================================================
# A "NEUTRAL" spacer bar is added twice per cancer type, purely to create a
# small visual gap between sectors in the polar plot.

cancer_types_present <- unique(gene_counts$NEOPLASM)

spacer_bars <- tibble(
     NEOPLASM    = cancer_types_present,
     Hugo_Symbol = NA_character_,
     TYPE        = "NEUTRAL",
     Count       = NA_integer_,
     Total       = NA_integer_,
     Percentage  = 0,
     gene_color  = NA_character_,
     Category    = ""
)
spacer_bars <- spacer_bars[rep(seq_len(nrow(spacer_bars)), each = 2), ]

# Order: cancer type -> variant status (ON, OFF, spacer) -> gene
df <- bind_rows(gene_counts, spacer_bars) %>%
     mutate(
          NEOPLASM    = factor(NEOPLASM, levels = cancer_types_to_plot),
          TYPE        = factor(TYPE, levels = c("ON", "OFF", "NEUTRAL")),
          Hugo_Symbol = factor(Hugo_Symbol, levels = gene_names)
     ) %>%
     arrange(NEOPLASM, TYPE, Hugo_Symbol) %>%
     mutate(ID = row_number())

number_of_bars <- nrow(df)
max_val_sqrt   <- max(sqrt(df$Percentage), na.rm = TRUE)


# ==============================================================================
# 6. Compute arcs and label positions
# ==============================================================================

# ON/OFF status arcs, drawn just outside the plot center
arc_data <- df %>%
     filter(Percentage > 0) %>%
     group_by(NEOPLASM, TYPE) %>%
     summarize(start = min(as.numeric(ID)) - 0.45, end = max(as.numeric(ID)) + 0.45, .groups = "drop") %>%
     rowwise() %>%
     do(data.frame(NEOPLASM = .$NEOPLASM, TYPE = .$TYPE, x = seq(.$start, .$end, length.out = 50)))


neoplasm_label_data <- df %>%
     mutate(id_num = as.numeric(ID)) %>%
     group_by(NEOPLASM) %>%
     summarize(id_avg = mean(id_num), .groups = "drop") %>%
     mutate(
          angle_deg  = (id_avg - 0.5) / number_of_bars * 360,
          text_angle = angle_deg,
          # Flip upside-down labels the right way up
          text_angle = ifelse(angle_deg > 90 & angle_deg < 270, text_angle + 180, text_angle)
     )

# Radial percentage labels, one per bar
label_data <- df %>%
     mutate(
          id_num = as.numeric(ID),
          angle  = 90 - 360 * (id_num - 0.5) / number_of_bars,
          hjust  = ifelse(angle < -90, 1, 0),
          angle  = ifelse(angle < -90, angle + 180, angle)
     )


# ==============================================================================
# 7. Build and save the plot
# ==============================================================================

p <- ggplot(df, aes(x = ID, y = sqrt(Percentage), fill = Hugo_Symbol)) +
     # Gene bars
     geom_col(width = 0.9, color = "black", linewidth = 0.05) +
     coord_polar(start = 0) +
     # ON/OFF status arcs
     geom_path(
          data = arc_data,
          aes(x = x, y = -max_val_sqrt * 0.06, color = TYPE, group = interaction(NEOPLASM, TYPE)),
          linewidth = 0.8, inherit.aes = FALSE
     ) +
     # Percentage labels
     geom_text(
          data = label_data %>% filter(Percentage > 0),
          aes(x = ID, y = sqrt(Percentage) + max_val_sqrt * 0.05, label = paste0(Percentage, "%"),
              angle = angle, hjust = hjust),
          color = "black", size = 3.5, fontface = "bold", inherit.aes = FALSE, show.legend = FALSE
     ) +
     scale_fill_manual(
          values = gene_color_lookup, breaks = gene_names,
          name = "Genes", guide = guide_legend(ncol = 1, order = 1)
     ) +
     scale_color_manual(
          values = c("ON" = "#228B22", "OFF" = "#CD0000", "NEUTRAL" = "transparent"),
          labels = c("ON" = "ON Tumor", "OFF" = "OFF Tumor"),
          breaks = c("ON", "OFF"), name = "Variant Status",
          guide = guide_legend(override.aes = list(linewidth = 4), order = 2)
     ) +
     geom_hline(yintercept = 0, color = "gray50", linewidth = 0.4) +
     ylim(c(-max_val_sqrt * 0.8, max_val_sqrt * 1.3)) +
     theme_void() +
     theme(
          plot.margin     = margin(20, 20, 20, 20),
          legend.position = "right",
          legend.title    = element_text(face = "bold", size = 12),
          legend.text     = element_text(size = 11, face = "italic")
     )

# Vector PDF (useDingbats = FALSE avoids symbol rendering issues in some readers)
pdf(file = f_pdf_out, width = 15, height = 12, useDingbats = FALSE, onefile = TRUE)
print(p)
dev.off()
message("PDF saved to: ", f_pdf_out)

# Transparent-background SVG (editable in PowerPoint/Illustrator)
ggsave(
     filename = f_svg_out,
     plot = p,
     width = 15, height = 12, units = "in",
     device = "svg", bg = "transparent"
)
message("SVG saved to: ", f_svg_out)
