# ==============================================================================
# Supplementary Figure - Oncoplot of Germline and Second-Hit Variants
# ==============================================================================
#
# Description:
#   Prepares the germline variant dataset (merging VAF data, second-hit
#   somatic variant requests, and germline-confirmation labels) and generates
#   an oncoplot (ComplexHeatmap::oncoPrint) summarizing germline variants,
#   second-hit events, VAF status, MUTYH biallelic variants, and secondary
#   findings across the patient cohort, stratified by cancer type.
#
# ==============================================================================

# ---- 1. Setup ----------------------------------------------------------------

library(readxl)
library(writexl)
library(openxlsx)
library(dplyr)
library(tidyr)
library(purrr)
library(stringr)
library(ComplexHeatmap)
library(grid)
library(RColorBrewer)

# Close any open graphics device before starting
if (!is.null(dev.list())) dev.off()

# ---- Paths (EDIT THESE before running) ----------------------------------
# Point data_dir / output_dir at your local folders. All three input files
# listed above are expected inside data_dir.

data_dir   <- "data"
output_dir <- "output"

dir.create(output_dir, showWarnings = FALSE, recursive = TRUE)

# Input / intermediate / output file paths, defined once and reused
# throughout so write and read steps always stay in sync.
f_raw             <- file.path(data_dir, "germline_mutations_raw.xlsx")
f_vaf_excluded    <- file.path(data_dir, "vaf_excluded_second_hit_cases.xlsx")
f_second_hit_req  <- file.path(data_dir, "second_hit_requests.xlsx")

f_merged_vaf      <- file.path(output_dir, "dataset_merged_vaf.xlsx")
f_with_second_hit <- file.path(output_dir, "dataset_with_second_hit.xlsx")
f_final_annotated <- file.path(output_dir, "dataset_final_annotated.xlsx")
f_oncoplot_pdf     <- file.path(output_dir, "oncoplot_supplementary_figure1.pdf")


# ==============================================================================
# 2. Merge main dataset with VAF values for second-hit-excluded cases
# ==============================================================================

# --- 2.1 Load and standardize the main variant dataset ---
main_data <- read_excel(f_raw) %>%
  mutate(
    PATIENT_ID = as.character(PATIENT_ID),
    PATIENT_ID = str_trim(PATIENT_ID),
    PATIENT_ID = str_pad(PATIENT_ID, width = 8, side = "left", pad = "0")
  )

message("Patient ID length after zero-padding (all values should be 8):")
print(table(nchar(main_data$PATIENT_ID)))

# --- 2.2 Load and standardize the VAF-excluded second-hit cases ---
vaf_excluded_clean <- read_excel(f_vaf_excluded) %>%
  rename(
    PATIENT_ID  = `Internal ID patient (MRN)`,
    Hugo_Symbol = Gene,
    TIER        = Tier,
    NOM_C       = `Nomenclature c.`,
    NOM_P       = `Nomenclature p.`,
    VAF         = `Variant allele frequency (VAF)`,
    TYPE        = `ON/OFF`
  ) %>%
  mutate(PATIENT_ID = as.character(PATIENT_ID)) %>%
  # Keep only the first record per variant (drop duplicate second-hit entries)
  distinct(PATIENT_ID, Hugo_Symbol, TIER, NOM_C, NOM_P, VAF, TYPE, .keep_all = TRUE) %>%
  select(PATIENT_ID, Hugo_Symbol, TIER, NOM_C, NOM_P, VAF, TYPE, `VAF > 60 %`)

# --- 2.3 Standardize VAF formatting across both datasets ---
# Handles values using a comma as decimal separator, a trailing "%" sign,
# and proportions (0-1) that need converting to percentages (0-100).
clean_vaf <- function(x) {
  x %>%
    as.character() %>%
    str_replace(",", ".") %>%
    str_remove("%") %>%
    str_trim() %>%
    as.numeric() %>%
    { ifelse(. <= 1, . * 100, .) } %>%
    round(2)
}

main_data <- main_data %>%
  mutate(VAF = clean_vaf(VAF))

vaf_excluded_clean <- vaf_excluded_clean %>%
  mutate(VAF = clean_vaf(VAF))

# --- 2.4 Merge (expected to be one-to-one / many-to-one, no row duplication) ---
merged_data <- main_data %>%
  left_join(
    vaf_excluded_clean,
    by = c("PATIENT_ID", "Hugo_Symbol", "TIER", "NOM_C", "NOM_P", "VAF", "TYPE")
  )

# Sanity check: row count should be unchanged after the merge
message("Rows before merge: ", nrow(main_data))
message("Rows after merge:  ", nrow(merged_data))

write_xlsx(merged_data, f_merged_vaf)


# ==============================================================================
# 3. Integrate second-hit variant requests
# ==============================================================================

merged_data <- read_excel(f_merged_vaf)

second_hit_clean <- read_excel(f_second_hit_req) %>%
  rename(
    PATIENT_ID = `Internal ID patient (MRN)`,
    GENE       = Gene,
    TIER       = Tier,
    NEOPLASM   = Neoplasm,
    NOM_C      = `Nomenclature c.`,
    NOM_P      = `Nomenclature p.`
  ) %>%
  mutate(PATIENT_ID = as.character(PATIENT_ID)) %>%
  select(PATIENT_ID, NEOPLASM, GENE, TIER, NOM_C, NOM_P, VAF, `Variant Identification`)

# Keep only the rows explicitly flagged as "Second Hit"
second_hit_only <- second_hit_clean[second_hit_clean$`Variant Identification` == "Second Hit", ]

# Add the tracking column to the main dataset if not already present
if (!"Variant Identification" %in% colnames(merged_data)) {
  merged_data$`Variant Identification` <- NA
}

second_hit_only <- second_hit_only %>%
  mutate(`Variant Identification` = "Second Hit")

# Align columns between the two datasets before binding
missing_cols <- setdiff(colnames(merged_data), colnames(second_hit_only))
for (col_name in missing_cols) {
  second_hit_only[[col_name]] <- NA
}
second_hit_only <- second_hit_only %>%
  select(all_of(colnames(merged_data)))

# Combine the main dataset with the newly added second-hit rows
combined_data <- bind_rows(merged_data, second_hit_only)

wb <- createWorkbook()
addWorksheet(wb, "Combined")
writeData(wb, "Combined", combined_data)

blue_style     <- createStyle(fontColour = "#0000FF")
new_rows_start <- nrow(merged_data) + 2
new_rows_end   <- nrow(combined_data) + 1

addStyle(
  wb, sheet = "Combined", style = blue_style,
  rows = new_rows_start:new_rows_end,
  cols = 1:ncol(combined_data),
  gridExpand = TRUE
)

saveWorkbook(wb, f_with_second_hit, overwrite = TRUE)


# ==============================================================================
# 4. Flag germline variants confirmed by a matching second-hit event
# ==============================================================================

full_data <- read.xlsx(f_with_second_hit)

second_hit_rows <- full_data %>%
  filter(Variant.Identification == "Second Hit")

# Find germline variants (unlabeled rows) that share patient/gene/neoplasm
# with a confirmed second-hit event
germline_confirmed_matches <- full_data %>%
  filter(
    is.na(Variant.Identification),  # not yet labeled
    GERMLINE == "Yes"                # true germline variants only
  ) %>%
  inner_join(
    second_hit_rows %>% select(PATIENT_ID, GENE, NEOPLASM),
    by = c("PATIENT_ID", "GENE", "NEOPLASM")
  )

message("Germline / second-hit matches found for confirmation:")
print(table(paste(germline_confirmed_matches$PATIENT_ID,
                   germline_confirmed_matches$GENE,
                   germline_confirmed_matches$NEOPLASM)))

full_data <- full_data %>%
  mutate(
    Variant.Identification = ifelse(
      is.na(Variant.Identification) &
        GERMLINE == "Yes" &
        paste(PATIENT_ID, GENE, NEOPLASM) %in%
          paste(germline_confirmed_matches$PATIENT_ID,
                germline_confirmed_matches$GENE,
                germline_confirmed_matches$NEOPLASM),
      "Germline Confirmed",
      Variant.Identification
    )
  )

wb <- createWorkbook()
addWorksheet(wb, "Updated")
writeData(wb, "Updated", full_data)

blue_style  <- createStyle(fontColour = "#0000FF")
green_style <- createStyle(fontColour = "#008000")

blue_rows  <- which(full_data$Variant.Identification == "Second Hit") + 1
green_rows <- which(full_data$Variant.Identification == "Germline Confirmed") + 1

addStyle(wb, "Updated", blue_style,  rows = blue_rows,  cols = 1:ncol(full_data), gridExpand = TRUE)
addStyle(wb, "Updated", green_style, rows = green_rows, cols = 1:ncol(full_data), gridExpand = TRUE)

saveWorkbook(wb, f_final_annotated, overwrite = TRUE)


# ==============================================================================
# 5. Build the oncoplot alteration matrix
# ==============================================================================

annotated_data <- read_excel(f_final_annotated)


has_secondary_findings <- any(grepl("(?i)^secondary[ ._]findings$", names(annotated_data), perl = TRUE))

if (has_secondary_findings) {
  names(annotated_data) <- sub("(?i)^secondary[ ._]findings$", "SECONDARY_FINDINGS",
                                names(annotated_data), perl = TRUE)
} else {
  message("NOTE: no 'Secondary Findings' column found in the input data - ",
          "the Secondary Findings row will be omitted from the oncoplot. ",
          "Check whether this information needs to be merged in from a ",
          "separate source file.")
  annotated_data$SECONDARY_FINDINGS <- NA_character_
}

# --- 5.1 Keep germline and second-hit variants only ---
annotated_data <- annotated_data %>%
  filter(GERMLINE == "Yes" | Variant.Identification == "Second Hit") %>%
  mutate(PATIENT_NEOP = paste(PATIENT_ID, NEOPLASM, sep = "_"))

# --- 5.2 Identify "double" variants: germline + second hit in the same gene ---
double_variants <- annotated_data %>%
  group_by(PATIENT_ID, GENE) %>%
  filter(any(GERMLINE == "Yes") & any(Variant.Identification == "Second Hit")) %>%
  ungroup()

# --- 5.3 Remaining ("pure") germline variants, excluding the doubles above ---
additional_germline <- annotated_data %>%
  filter(GERMLINE == "Yes") %>%
  filter(!paste(PATIENT_ID, GENE) %in% paste(double_variants$PATIENT_ID, double_variants$GENE))

# --- 5.4 Combine into the final variant set used for plotting ---
combined_variants <- bind_rows(double_variants, additional_germline) %>%
  select(PATIENT_ID, GENE, PATIENT_NEOP, NEOPLASM, VAF, Variant.Identification,
         GERMLINE, `VAF.>.60.%`, SECONDARY_FINDINGS)

# --- 5.5 Assign an alteration type per (gene, patient/neoplasm) pair ---
mat_standard <- combined_variants %>%
  select(GENE, PATIENT_NEOP, Variant.Identification, GERMLINE, `VAF.>.60.%`) %>%
  mutate(types = pmap(list(GERMLINE, Variant.Identification, `VAF.>.60.%`), function(g, v, vf) {
    res <- c()
    if (!is.na(g) && g == "Yes") res <- c(res, "GERMLINE")
    if (!is.na(v) && v == "Second Hit") res <- c(res, "Second Hit")
    if (!is.na(vf) && vf == "MUTYH") {
      res <- c(res, "MUTYH Biallelic Variant")
    } else {
      if (!is.na(vf) && vf == "Yes") res <- c(res, "VAF > 60%")
      if (!is.na(vf) && vf == "No")  res <- c(res, "VAF < 60%")
    }
    return(res)
  })) %>%
  unnest(types) %>%
  select(GENE, PATIENT_NEOP, types) %>%
  distinct()

# --- 5.6 Secondary findings shown as a single dedicated row ---
mat_secondary_findings <- combined_variants %>%
  filter(SECONDARY_FINDINGS == "Yes") %>%
  mutate(GENE = "Secondary Findings", types = "Secondary Finding") %>%
  select(GENE, PATIENT_NEOP, types) %>%
  distinct()

# --- 5.7 Order genes by mutation frequency (Secondary Findings row last) ---
gene_freq_order <- mat_standard %>%
  group_by(GENE) %>%
  summarise(n = n()) %>%
  arrange(desc(n)) %>%
  pull(GENE)

mat_total_df <- bind_rows(
  mat_standard %>%
    pivot_wider(names_from = PATIENT_NEOP, values_from = types,
                values_fn = function(x) paste(unique(x), collapse = ";")),
  mat_secondary_findings %>%
    pivot_wider(names_from = PATIENT_NEOP, values_from = types)
) %>%
  filter(!is.na(GENE)) %>%
  mutate(GENE = factor(GENE, levels = c(gene_freq_order, "Secondary Findings"))) %>%
  arrange(GENE) %>%
  replace(is.na(.), "") %>%
  as.data.frame()

mat_data <- as.matrix(mat_total_df[, -1])
rownames(mat_data) <- mat_total_df$GENE


# ==============================================================================
# 6. Prepare sample metadata and annotations
# ==============================================================================

# --- 6.1 Standardize cancer type labels for display ---
cancer_type_map <- c(
  "cholangiocarcinoma" = "Cholangiocarcinoma",
  "Colorectum"          = "Colorectal Cancer",
  "Endometrium"         = "Endometrial Cancer",
  "GIST"                = "Gastrointestinal Stromal Tumor",
  "Lung"                = "Lung Cancer",
  "Melanoma"            = "Melanoma",
  "Ovary"               = "Ovarian Cancer",
  "Pancreas"            = "Pancreatic Cancer",
  "Prostate"            = "Prostate Cancer",
  "Thyroid"             = "Thyroid Cancer"
)

metadata_sync <- combined_variants %>%
  select(PATIENT_NEOP, NEOPLASM) %>%
  distinct() %>%
  mutate(NEOPLASM = ifelse(NEOPLASM %in% names(cancer_type_map),
                            cancer_type_map[NEOPLASM], NEOPLASM))

# Align metadata row order with the matrix columns
metadata_sync <- metadata_sync[match(colnames(mat_data), metadata_sync$PATIENT_NEOP), ]

# --- 6.2 Color palette for alteration types ---
alteration_colors <- c(
  "GERMLINE"                = "#377EB8",
  "Second Hit"               = "#FF7F00",
  "VAF > 60%"                = "#4DAF4A",
  "VAF < 60%"                = "darkgray",
  "MUTYH Biallelic Variant"  = "#f5f5f5",
  "Secondary Finding"        = "red"
)

# --- 6.3 Color palette for the cancer type annotation ---
cancer_types_unique <- sort(unique(metadata_sync$NEOPLASM))
n_cancer_types <- length(cancer_types_unique)
cancer_type_colors <- setNames(
  colorRampPalette(brewer.pal(min(n_cancer_types, 12), "Paired"))(n_cancer_types),
  cancer_types_unique
)
message("Cancer type color mapping:")
print(cancer_type_colors)

bottom_annotation <- HeatmapAnnotation(
  Neoplasia = metadata_sync$NEOPLASM,
  col = list(Neoplasia = cancer_type_colors),
  show_legend = TRUE,
  show_annotation_name = FALSE,
  annotation_legend_param = list(
    Neoplasia = list(
      title = "Cancer Type",
      title_gp = gpar(fontsize = 16, fontface = "bold"),
      labels_gp = gpar(fontsize = 14),
      grid_height = unit(8, "mm"),
      grid_width = unit(8, "mm")
    )
  )
)

# --- 6.4 Top annotation: per-sample max VAF and multi-tumor patient markers ---


vaf_per_sample <- combined_variants %>%
  filter(!Variant.Identification %in% c("Second Hit")) %>%
  filter(!is.na(VAF)) %>%
  group_by(PATIENT_NEOP) %>%
  summarize(vaf_val = max(VAF, na.rm = TRUE)) %>%
  ungroup()

vaf_ordered <- vaf_per_sample$vaf_val[match(colnames(mat_data), vaf_per_sample$PATIENT_NEOP)]

# Identify patients with more than one neoplasm (synchronous tumors)
patient_ids_per_column <- sub("_.*", "", colnames(mat_data))

multi_neoplasm_patients <- combined_variants %>%
  group_by(PATIENT_ID) %>%
  summarize(n = n_distinct(NEOPLASM)) %>%
  filter(n > 1) %>%
  pull(PATIENT_ID)

# Distinct color per multi-neoplasm patient, for visual linking across columns
set.seed(123)  # for reproducible color assignment
multi_neoplasm_colors <- setNames(
  colorRampPalette(brewer.pal(8, "Dark2"))(length(multi_neoplasm_patients)),
  multi_neoplasm_patients
)

column_point_colors <- ifelse(
  patient_ids_per_column %in% multi_neoplasm_patients,
  multi_neoplasm_colors[patient_ids_per_column],
  NA
)

top_annotation <- HeatmapAnnotation(
  # Row 1: colored point marking samples from multi-neoplasm patients
  MultiTumor = anno_points(
    ifelse(!is.na(column_point_colors), 1, NA),
    pch = 16,
    gp = gpar(col = column_point_colors),
    size = unit(4, "mm"),
    ylim = c(0, 2),
    axis = FALSE
  ),
  # Row 2: VAF barplot
  VAF = anno_barplot(
    vaf_ordered,
    height = unit(2, "cm"),
    gp = gpar(fill = "black", col = NA),
    axis_param = list(side = "left", at = c(0, 20, 40, 60, 80, 100), gp = gpar(fontsize = 12))
  ),
  annotation_label = c("Synchronous Tumors", "VAF (%)"),
  show_annotation_name = TRUE
)


# ==============================================================================
# 7. Custom alteration glyphs (split-cell rendering)
# ==============================================================================


alter_fun <- list(
  background = function(x, y, w, h) {
    grid.rect(x, y, w, h, gp = gpar(fill = "#f5f5f5", col = "white"))
  },
  # Upper-left triangle: germline variant
  GERMLINE = function(x, y, w, h) {
    grid.polygon(
      unit.c(x - w * 0.5, x - w * 0.5, x + w * 0.5),
      unit.c(y - h * 0.5, y + h * 0.5, y + h * 0.5),
      gp = gpar(fill = alteration_colors["GERMLINE"], col = NA)
    )
  },
  # Lower-right triangle: second hit
  `Second Hit` = function(x, y, w, h) {
    grid.polygon(
      unit.c(x - w * 0.5, x + w * 0.5, x + w * 0.5),
      unit.c(y - h * 0.5, y - h * 0.5, y + h * 0.5),
      gp = gpar(fill = alteration_colors["Second Hit"], col = NA)
    )
  },
  # Lower-right triangle: VAF > 60%
  `VAF > 60%` = function(x, y, w, h) {
    grid.polygon(
      unit.c(x - w * 0.5, x + w * 0.5, x + w * 0.5),
      unit.c(y - h * 0.5, y - h * 0.5, y + h * 0.5),
      gp = gpar(fill = alteration_colors["VAF > 60%"], col = NA)
    )
  },
  # Lower-right triangle: VAF < 60%
  `VAF < 60%` = function(x, y, w, h) {
    grid.polygon(
      unit.c(x - w * 0.5, x + w * 0.5, x + w * 0.5),
      unit.c(y - h * 0.5, y - h * 0.5, y + h * 0.5),
      gp = gpar(fill = alteration_colors["VAF < 60%"], col = NA)
    )
  },
  # Lower-right triangle + asterisk: MUTYH biallelic variant
  `MUTYH Biallelic Variant` = function(x, y, w, h) {
    grid.polygon(
      unit.c(x - w * 0.5, x + w * 0.5, x + w * 0.5),
      unit.c(y - h * 0.5, y - h * 0.5, y + h * 0.5),
      gp = gpar(fill = alteration_colors["MUTYH Biallelic Variant"], col = NA)
    )
    grid.text("*", x = x + w * 0.1, y = y - h * 0.4,
              gp = gpar(col = "black", fontsize = 20, fontface = "bold"))
  },
  # Secondary finding row: asterisk only
  `Secondary Finding` = function(x, y, w, h) {
    grid.text("*", x = x, y = y - h * 0.2,
              gp = gpar(col = "red", fontsize = 20, fontface = "bold"))
  }
)


# ==============================================================================
# 8. Generate and save the oncoplot
# ==============================================================================

pdf(f_oncoplot_pdf, width = 26, height = 14)

p <- oncoPrint(
  mat_data,
  alter_fun = alter_fun,
  col = alteration_colors,
  top_annotation = top_annotation,
  bottom_annotation = bottom_annotation,
  column_split = metadata_sync$NEOPLASM,
  column_title = NULL,
  show_pct = FALSE,
  row_order = 1:nrow(mat_data),
  show_row_names = TRUE,
  row_names_side = "left",
  row_names_gp = gpar(
    fontface = "italic",
    fontsize = 14,
    # Hide placeholder row labels (e.g. spacer rows), if present
    col = ifelse(rownames(mat_data) %in% c("buffer", "Secondary Findings"), "transparent", "black")
  ),
  right_annotation = NULL,
  column_gap = unit(1, "mm"),
  alter_fun_is_vectorized = FALSE,
  remove_empty_columns = TRUE,
  remove_empty_rows = TRUE,
  heatmap_legend_param = list(
    title = "Alterations & VAF Status",
    at = c("GERMLINE", "Second Hit", "MUTYH Biallelic Variant", "VAF > 60%", "VAF < 60%", "Secondary Finding"),
    labels = c("Germline Variant", "Second Hit (Somatic)", "MUTYH Biallelic Variants", "VAF > 60%", "VAF < 60%", "Secondary Finding"),
    grid_height = unit(8, "mm"),
    grid_width = unit(8, "mm"),
    labels_gp = gpar(fontsize = 14),
    title_gp = gpar(fontsize = 16, fontface = "bold")
  )
)

# Large left padding leaves room for gene name labels
draw(p, merge_legend = TRUE, heatmap_legend_side = "right",
     padding = unit(c(10, 10, 10, 10), "mm"))

# Dashed red reference line at VAF = 60%, drawn on each cancer-type panel
n_panels <- length(unique(metadata_sync$NEOPLASM))
for (s in 1:n_panels) {
  decorate_annotation("VAF", slice = s, {
    grid.lines(
      x = unit(c(0, 1), "npc"),
      y = unit(c(60, 60), "native"),
      gp = gpar(col = "red", lty = "dashed", lwd = 1.5)
    )
  })
}

dev.off()


# ==============================================================================
# 9. Summary: patients with multiple neoplasms
# ==============================================================================

multi_neoplasm_summary <- combined_variants %>%
  group_by(PATIENT_ID) %>%
  summarize(
    n_neoplasms = n_distinct(NEOPLASM),
    neoplasms   = paste(sort(unique(NEOPLASM)), collapse = " & ")
  ) %>%
  filter(n_neoplasms > 1) %>%
  arrange(desc(n_neoplasms))

message("\n--- Patients with multiple (synchronous) neoplasms ---")
print(as.data.frame(multi_neoplasm_summary))
message("-------------------------------------------------------")
