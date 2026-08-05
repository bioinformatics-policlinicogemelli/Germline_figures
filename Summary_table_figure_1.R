# ==============================================================================
# Germline Testing Summary Tables
# Table 1: Cases / Tested / Positive by Cancer Type
# Table 2: Confirmed-Variant Percentage by Gene x Cancer Type
# ==============================================================================
#
# Description:
#   Produces the two underlying tables behind the manuscript's germline
#   testing summary figure:
#
#   Table 1 (one row per cancer type, plus an "Any Cancer" total row):
#     - Total cases referred for genetic counseling (ESMO-eligible)
#     - Total cases that actually completed germline testing
#     - Positive test results (>=1 confirmed germline variant)
#     - Positive rate (%)
#
#   Table 2 (cancer type x gene matrix, plus an "Any Cancer" row):
#     - For each gene in the panel, the % of tested variants confirmed
#       germline-positive, broken down by cancer type.
#
#   NOTE: the published figure additionally overlays a horizontal bar chart
#   next to Table 1 and merges Table 1 + Table 2 into a single sheet. That
#   final layout step is done manually in Excel and is not reproduced here -
#   this script only produces the underlying numbers
#
#
# Input files (place in `data_dir`):
#   - clinical_sample_germline.csv   Per-patient/sample clinical file
#                                     (cBioPortal-style "data_clinical_sample"
#                                     export), used for Table 1
#   - mutation_data_germline.txt     Per-variant mutation file (tab-delimited),
#                                     used for Table 2
#
#   NOTE: file names above are placeholders - edit the paths in the "Paths"
#   section below to match your local files.
#
# Output files (written to `output_dir`):
#   - table1_cases_tested_positive.txt
#   - table2_gene_by_cancer_type.txt
#
# Required packages: readxl, dplyr, tidyr
#
# ==============================================================================

# ---- 1. Setup -----------------------------------------------------------

library(readxl)
library(dplyr)
library(tidyr)

# ---- Paths (EDIT THESE before running) -----------------------------------
data_dir   <- "data"
output_dir <- "output"

dir.create(output_dir, showWarnings = FALSE, recursive = TRUE)

f_clinical_sample <- file.path(data_dir, "clinical_sample_germline.csv")
f_mutation_data    <- file.path(data_dir, "mutation_data_germline.txt")

f_table1_out <- file.path(output_dir, "table1_cases_tested_positive.txt")
f_table2_out <- file.path(output_dir, "table2_gene_by_cancer_type.txt")

# Cancer type display order, shared by both tables
cancer_type_order <- c(
     "Ovary", "Endometrium", "Lung", "Colorectum", "Pancreas",
     "Cholangiocarcinoma", "Prostate", "Melanoma", "Thyroid", "GIST",
     "Breast", "Any Cancer"
)

# Gene panel included in Table 2. The original script referenced an
# `esmo_genes` list that was never actually defined - reconstructed here
# from the validated published table (26-gene ESMO hereditary cancer panel).
esmo_genes <- c(
     "ATM", "BAP1", "BARD1", "BRCA1", "BRCA2", "BRIP1", "CHEK2", "DICER1",
     "FH", "FLCN", "MLH1", "MSH2", "MSH6", "MUTYH", "NF1", "PALB2", "PMS2",
     "POLE", "PTCH1", "PTEN", "RAD51C", "RAD51D", "RET", "SDHA", "SDHB", "SDHC"
)


# ==============================================================================
# 2. Table 1 - Cases, tested, positive, and positive rate by cancer type
# ==============================================================================

# Per-patient/sample clinical file. PATIENT_ID is zero-padded to 8 digits,
# consistent with the other scripts in this repo.
clinical_sample <- read.csv(f_clinical_sample, skip = 4, stringsAsFactors = FALSE)
clinical_sample$PATIENT_ID <- sprintf("%08s", as.character(clinical_sample$PATIENT_ID))
clinical_sample$PATIENT_ID <- gsub(" ", "0", clinical_sample$PATIENT_ID)

message("Patient ID length after zero-padding (should be all 8):")
print(table(nchar(clinical_sample$PATIENT_ID)))
message("Unique patients in clinical sample file: ", length(unique(clinical_sample$PATIENT_ID)))

table1 <- clinical_sample %>%
     group_by(NEOPLASM) %>%
     summarise(
          # This file is pre-filtered to ESMO-eligible patients, so every row
          # already satisfies GENETIC_CONS == "Yes" - this is effectively a count
          # of patients per cancer type.
          `Total cases (ESMO)` = sum(GENETIC_CONS == "Yes", na.rm = TRUE),
          `Total cases tested` = sum(GERMLINE_TEST == "Yes, test performed", na.rm = TRUE),
          `Positive tests`     = sum(GERMLINE == "Yes", na.rm = TRUE)
     ) %>%
     mutate(`Positive rate (%)` = round(`Positive tests` / `Total cases tested` * 100, 1)) %>%
     arrange(desc(`Total cases (ESMO)`))

# Cancer types with zero tested patients (e.g. Breast) produce a 0/0 divide,
# which R evaluates to NaN rather than 0 - this is expected and gets
# corrected to a displayed "0" manually during the Excel finishing step.
any_cancer_row <- table1 %>%
     summarise(
          NEOPLASM             = "Any Cancer",
          `Total cases (ESMO)` = sum(`Total cases (ESMO)`),
          `Total cases tested` = sum(`Total cases tested`),
          `Positive tests`     = sum(`Positive tests`)
     ) %>%
     mutate(`Positive rate (%)` = round(`Positive tests` / `Total cases tested` * 100, 1))

table1 <- bind_rows(table1, any_cancer_row)

write.table(table1, file = f_table1_out, sep = "\t", row.names = FALSE, quote = FALSE)
message("Table 1 saved to: ", f_table1_out)


# ==============================================================================
# 3. Table 2 - Confirmed-variant percentage by gene and cancer type
# ==============================================================================

# Per-variant mutation file (one row per tested variant)
mutation_data <- read.delim(f_mutation_data, header = TRUE, sep = "\t", stringsAsFactors = FALSE)

gene_summary <- mutation_data %>%
     group_by(NEOPLASM, GENE) %>%
     summarise(
          total     = n(),
          confirmed = sum(GERMLINE == "Yes", na.rm = TRUE),
          .groups = "drop"
     ) %>%
     mutate(`Confirmed %` = ifelse(total == 0, NA, round(confirmed / total * 100, 1)))

any_cancer_genes <- gene_summary %>%
     group_by(GENE) %>%
     summarise(
          total         = sum(total),
          confirmed     = sum(confirmed),
          `Confirmed %` = round(confirmed / total * 100, 1),
          .groups = "drop"
     ) %>%
     mutate(NEOPLASM = "Any Cancer") %>%
     select(NEOPLASM, GENE, `Confirmed %`)

gene_summary <- bind_rows(
     gene_summary %>% select(NEOPLASM, GENE, `Confirmed %`),
     any_cancer_genes
)

# Pivot to wide format: one row per cancer type, one column per gene
table2 <- gene_summary %>%
     pivot_wider(names_from = GENE, values_from = `Confirmed %`, values_fill = NA)

# Keep only genes in the ESMO panel, alphabetically ordered
table2 <- table2 %>%
     select(NEOPLASM, sort(intersect(colnames(.), esmo_genes)))

# Apply the same cancer type ordering as Table 1, round to whole percentages
table2$NEOPLASM <- factor(table2$NEOPLASM, levels = cancer_type_order)
table2 <- table2 %>%
     arrange(NEOPLASM) %>%
     mutate(across(where(is.numeric), ~ round(.)))

write.table(table2, file = f_table2_out, sep = "\t", row.names = FALSE, quote = FALSE)
message("Table 2 saved to: ", f_table2_out)