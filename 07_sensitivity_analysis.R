# 07_sensitivity_analysis.R
# Sensitivity analyses using African-ancestry data
# Addresses reviewer concerns about TCGA population representativeness

library(tidyverse)
library(data.table)
library(here)
library(corrplot)
library(ggrepel)
library(patchwork)
library(ggsci)

# Load configuration and set dplyr precedence
library(conflicted)
conflict_prefer("select", "dplyr")
conflict_prefer("filter", "dplyr")

cat("\n")
cat(strrep("=", 80), "\n")
cat("SENSITIVITY ANALYSIS USING AFRICAN-ANCESTRY DATA\n")
cat(strrep("=", 80), "\n\n")

# Load existing data

panel_70 <- read_csv(here("results", "tables", "panel_70_genes.csv"),
                     show_col_types = FALSE)

integrated_scores <- read_csv(here("data", "processed", "integrated_gene_scores.csv"),
                              show_col_types = FALSE)

pancancer_all <- read_csv(here("data", "processed", "pancancer_mutation_scores.csv"),
                          show_col_types = FALSE)

tcga_mutations <- readRDS(here("data", "processed", "tcga_mutations_top15_africa.rds"))

# Source helper for TCGA study mapping
source(here("scripts", "02_database_mining.R"), local = FALSE)
# Alternatively, define the function inline:
get_tcga_studies_for_top15_local <- function() {
  tribble(
    ~study_id, ~cancer_type, ~globocan_rank,
    "brca_tcga_pan_can_atlas_2018", "Breast", 1,
    "prad_tcga_pan_can_atlas_2018", "Prostate", 2,
    "cesc_tcga_pan_can_atlas_2018", "Cervical", 3,
    "lihc_tcga_pan_can_atlas_2018", "Liver", 4,
    "coadread_tcga_pan_can_atlas_2018", "Colorectal", 5,
    "luad_tcga_pan_can_atlas_2018", "Lung", 6,
    "lusc_tcga_pan_can_atlas_2018", "Lung", 6,
    "ov_tcga_pan_can_atlas_2018", "Ovarian", 7,
    "blca_tcga_pan_can_atlas_2018", "Bladder", 9,
    "stad_tcga_pan_can_atlas_2018", "Stomach", 10,
    "esca_tcga_pan_can_atlas_2018", "Esophageal", 11,
    "ucec_tcga_pan_can_atlas_2018", "Corpus uteri", 12,
    "paad_tcga_pan_can_atlas_2018", "Pancreatic", 14
  )
}



# SECTION 1: SCORING COMPONENT CORRELATION ----


cat("SECTION 1: Scoring component correlation analysis\n")
cat(strrep("-", 60), "\n\n")

# Extract scoring components for correlation analysis
score_components <- integrated_scores %>%
  filter(!is.na(combined_score), !is.na(african_score)) %>%
  select(
    `Pan-Cancer Score` = combined_score,
    `African-Weighted Score` = african_score,
    `Literature Score` = literature_score,
    `Expert Curation` = curated_score
  )

cor_matrix <- cor(score_components, use = "pairwise.complete.obs", method = "spearman")

cat("Spearman correlation matrix between scoring components:\n")
print(round(cor_matrix, 3))
cat("\n")

# Save correlation heatmap
pdf(here("results", "figures", "fig_supp_scoring_correlation.pdf"),
    width = 8, height = 7)

corrplot(cor_matrix,
         method = "color",
         type = "upper",
         addCoef.col = "black",
         tl.col = "black",
         tl.srt = 45,
         col = colorRampPalette(c("#4DBBD5", "white", "#E64B35"))(200),
         title = "Spearman Correlation Between IPS Scoring Components",
         mar = c(0, 0, 2, 0),
         number.cex = 1.2,
         cl.cex = 0.8)

dev.off()
cat("Saved: fig_supp_scoring_correlation.pdf\n\n")



# SECTION 2: TCGA BLACK/AFRICAN AMERICAN SUBSET ----


cat("SECTION 2: TCGA Black/African American subset analysis\n")
cat(strrep("-", 60), "\n\n")

tcga_dir <- here("data", "raw", "tcga_maf")
tcga_studies <- get_tcga_studies_for_top15_local()

# Read clinical data and extract Black patients
black_samples <- map_dfr(1:nrow(tcga_studies), function(i) {
  study_id <- tcga_studies$study_id[i]
  cancer_type <- tcga_studies$cancer_type[i]
  
  # Try multiple possible clinical file locations
  clinical_files <- c(
    file.path(tcga_dir, study_id, "data_clinical_patient.txt"),
    file.path(tcga_dir, study_id, "data_clinical_sample.txt")
  )
  
  clinical_file <- clinical_files[file.exists(clinical_files)][1]
  
  if (!is.na(clinical_file)) {
    # Read file, skipping metadata comment lines starting with #
    lines <- readLines(clinical_file, n = 10)
    skip_n <- sum(grepl("^#", lines))
    
    clinical <- tryCatch({
      fread(clinical_file, skip = skip_n, header = TRUE)
    }, error = function(e) {
      cat("  Warning: Could not read", clinical_file, "\n")
      return(NULL)
    })
    
    if (is.null(clinical)) return(tibble())
    
    # Find race column
    race_col <- grep("RACE|race", names(clinical), value = TRUE, ignore.case = TRUE)[1]
    patient_col <- grep("PATIENT_ID|SAMPLE_ID", names(clinical), value = TRUE, ignore.case = TRUE)[1]
    
    if (!is.na(race_col) && !is.na(patient_col)) {
      black_pts <- clinical %>%
        as_tibble() %>%
        filter(grepl("BLACK|AFRICAN", toupper(.data[[race_col]]))) %>%
        pull(.data[[patient_col]])
      
      if (length(black_pts) > 0) {
        cat("  ", study_id, ":", length(black_pts), "Black patients\n")
        return(tibble(
          study_id = study_id,
          cancer_type = cancer_type,
          patient_id = black_pts
        ))
      }
    }
  }
  return(tibble())
})

cat("\nTotal Black/African American patients identified:", nrow(black_samples), "\n")
cat("By cancer type:\n")
print(black_samples %>% count(cancer_type, sort = TRUE))
cat("\n")

# Filter mutations to Black patients
# TCGA barcodes: TCGA-XX-XXXX-01A-... Patient IDs: TCGA-XX-XXXX
tcga_mutations_black <- tcga_mutations %>%
  mutate(patient_id = str_extract(sample_id, "TCGA-[A-Z0-9]+-[A-Z0-9]+")) %>%
  inner_join(black_samples %>% select(patient_id), by = "patient_id")

cat("Mutations in Black patients:", nrow(tcga_mutations_black), "\n")
cat("Unique Black samples:", n_distinct(tcga_mutations_black$sample_id), "\n")
cat("Unique genes mutated:", n_distinct(tcga_mutations_black$gene), "\n\n")

# Calculate pan-cancer frequencies for Black patients
total_black_samples <- n_distinct(tcga_mutations_black$sample_id)

pancancer_black <- tcga_mutations_black %>%
  group_by(gene) %>%
  summarise(
    total_mutated_black = n_distinct(sample_id),
    .groups = "drop"
  ) %>%
  mutate(
    total_black_samples = total_black_samples,
    pancancer_freq_black = total_mutated_black / total_black_samples
  ) %>%
  arrange(desc(pancancer_freq_black))

# Compare with full cohort
comparison_tcga <- pancancer_black %>%
  select(gene, pancancer_freq_black) %>%
  inner_join(
    pancancer_all %>% select(gene, pancancer_frequency),
    by = "gene"
  ) %>%
  mutate(
    freq_ratio = ifelse(pancancer_frequency > 0,
                        pancancer_freq_black / pancancer_frequency, NA),
    enriched_in_black = freq_ratio > 1.2,
    depleted_in_black = freq_ratio < 0.8
  ) %>%
  arrange(desc(freq_ratio))

# Correlation test
cor_test_tcga <- cor.test(comparison_tcga$pancancer_frequency,
                          comparison_tcga$pancancer_freq_black,
                          method = "spearman")

cat("Spearman correlation (All TCGA vs Black-only frequencies):\n")
cat("  rho =", round(cor_test_tcga$estimate, 3),
    "  p =", format(cor_test_tcga$p.value, digits = 3), "\n\n")

# Filter to panel-70 genes
panel_genes_comparison <- comparison_tcga %>%
  filter(gene %in% panel_70$gene) %>%
  mutate(
    enrichment_status = case_when(
      enriched_in_black ~ "Higher in Black",
      depleted_in_black ~ "Lower in Black",
      TRUE ~ "Similar"
    )
  )

cat("Panel-70 genes with Black frequency data:", nrow(panel_genes_comparison), "\n")
cat("  Higher in Black patients:", sum(panel_genes_comparison$enriched_in_black, na.rm = TRUE), "\n")
cat("  Lower in Black patients:", sum(panel_genes_comparison$depleted_in_black, na.rm = TRUE), "\n")
cat("  Similar:", sum(!panel_genes_comparison$enriched_in_black &
                        !panel_genes_comparison$depleted_in_black, na.rm = TRUE), "\n\n")

# Save
write_csv(comparison_tcga,
          here("results", "tables", "sensitivity_tcga_black_vs_all.csv"))
write_csv(panel_genes_comparison,
          here("results", "tables", "sensitivity_panel70_black_frequencies.csv"))

# Correlation plot
p_tcga_corr <- ggplot(panel_genes_comparison,
                      aes(x = pancancer_frequency, y = pancancer_freq_black)) +
  geom_abline(slope = 1, intercept = 0, linetype = "dashed", color = "gray50") +
  geom_point(aes(color = enrichment_status), size = 3, alpha = 0.8) +
  geom_text_repel(aes(label = gene), size = 2.5, max.overlaps = 25,
                  segment.color = "grey60", segment.size = 0.3) +
  scale_color_manual(values = c("Higher in Black" = "#E64B35",
                                "Lower in Black" = "#4DBBD5",
                                "Similar" = "grey50"),
                     name = "Frequency\nComparison") +
  labs(
    title = "Mutation Frequencies: All TCGA vs Black/African American Subset",
    subtitle = paste("Panel-70 genes | Spearman rho =",
                     round(cor_test_tcga$estimate, 3)),
    x = "Pan-cancer frequency (All TCGA)",
    y = "Pan-cancer frequency (Black patients only)"
  ) +
  theme_minimal(base_size = 11) +
  coord_fixed(ratio = 1)

ggsave(here("results", "figures", "fig_supp_tcga_black_correlation.pdf"),
       p_tcga_corr, width = 10, height = 10, dpi = 300)

cat("Saved: fig_supp_tcga_black_correlation.pdf\n\n")



# SECTION 3: GENIE NHB DATA FROM WEN ET AL. ----


cat("SECTION 3: GENIE NHB data from Wen et al. (PMID: 40996301)\n")
cat(strrep("-", 60), "\n\n")

# Data extracted from Supplementary Table S4, Wen et al.
# Cancer Epidemiol Biomarkers Prev 2025;34(12), GENIE v13
# Values are mutation frequencies (%) for NHB and NHW patients
# Cancer types mapped to our 12 priority cancers:
#   NSCLC -> Lung, Breast -> Breast, Colorectal -> Colorectal,
#   Ovarian -> Ovarian, Prostate -> Prostate, Endometrial -> Corpus uteri,
#   Pancreatic -> Pancreatic, Bladder -> Bladder,
#   Esophagogastric -> Esophageal/Gastric, Hepatobiliary -> Liver

# NHB frequencies (%)
genie_nhb <- tribble(
  ~gene,      ~NSCLC, ~Breast, ~Colorectal, ~Ovarian, ~Prostate, ~Endometrial,
  ~Pancreatic, ~Bladder, ~Esophagogastric, ~Hepatobiliary,
  "TP53",        54,    56,    74,    70,    24,    73,    74,    51,    70,    44,
  "KRAS",        23,     1,    56,    12,     1,    10,    81,     4,     5,    10,
  "PIK3CA",       5,    26,    22,    12,     3,    31,     3,    16,     8,     4,
  "PTEN",         5,     5,     3,     3,     8,    39,     1,     3,     1,     2,
  "APC",          2,     2,    71,     3,     2,     4,     1,     4,     4,     6,
  "NF1",          7,     3,     4,     3,     5,     5,     3,     5,     2,     3,
  "ARID1A",       7,     5,     5,     5,     4,    18,     7,    14,    11,    11,
  "ATM",          7,     5,     7,     3,     5,     5,     5,     7,     4,     5,
  "BRCA2",        4,     4,     5,     5,     5,     5,     3,    10,     5,     3,
  "BRCA1",        2,     3,     2,     5,     1,     2,     2,     3,     1,     2,
  "FBXW7",        2,     2,    12,     1,     1,    17,     0,     5,     2,     1,
  "KMT2D",        5,     4,     9,     3,     3,    17,     2,    12,     4,     4,
  "ERBB2",        4,     5,     4,     3,     3,     4,     2,     9,     5,     4,
  "RET",          2,     1,     1,     1,     1,     2,     1,     2,     2,     1,
  "KIT",          2,     1,     2,     1,     1,     1,     1,     2,     2,     2,
  "FGFR3",        1,     1,     1,     2,     1,     1,     0,    18,     0,     0,
  "MET",          4,     2,     2,     2,     2,     3,     1,     3,     2,     2,
  "ALK",          4,     2,     4,     1,     2,     3,     1,     4,     5,     1,
  "BRAF",         4,     1,     9,     3,     1,     6,     1,     1,     1,     2,
  "CTNNB1",       1,     1,     1,     1,     1,     4,     1,     2,     1,    16,
  "EGFR",        24,     1,     1,     1,     1,     1,     1,     2,     1,     1,
  "POLE",         2,     1,     3,     2,     1,     7,     0,     3,     1,     1,
  "CREBBP",       4,     4,     4,     2,     2,     5,     1,    14,     5,     1,
  "ERBB4",        4,     1,     2,     1,     1,     2,     2,     2,     3,     4,
  "NRAS",         1,     0,     5,     1,     0,     1,     0,     2,     0,     1,
  "NOTCH1",       5,     3,     3,     4,     2,     4,     3,     8,     5,     3,
  "MSH6",         1,     1,     5,     1,     1,     8,     1,     4,     1,     1,
  "MSH2",         1,     1,     4,     1,     1,     5,     1,     3,     1,     1,
  "MLH1",         1,     1,     3,     1,     1,     5,     1,     2,     1,     1,
  "PALB2",        1,     1,     2,     1,     2,     2,     1,     3,     1,     1,
  "CHEK2",        1,     2,     2,     1,     2,     2,     1,     3,     1,     1,
  "FGFR2",        2,     2,     2,     2,     1,     5,     1,     3,     2,     1,
  "ESR1",         1,     1,     0,     0,     0,     0,     0,     0,     0,     0,
  "SMARCA4",      8,     2,     3,     2,     1,     2,     2,     7,     5,     2
)

# NHW frequencies (%)
genie_nhw <- tribble(
  ~gene,      ~NSCLC, ~Breast, ~Colorectal, ~Ovarian, ~Prostate, ~Endometrial,
  ~Pancreatic, ~Bladder, ~Esophagogastric, ~Hepatobiliary,
  "TP53",        47,    39,    75,    77,    26,    42,    65,    50,    73,    31,
  "KRAS",        31,     1,    44,     8,     1,    20,    80,     7,     5,    14,
  "PIK3CA",       5,    34,    17,     8,     4,    44,     2,    19,     7,     4,
  "PTEN",         4,     6,     3,     3,     8,    44,     1,     3,     1,     1,
  "APC",          3,     3,    74,     3,     2,     4,     1,     6,     4,     6,
  "NF1",          8,     3,     4,     3,     3,     4,     3,     6,     3,     3,
  "ARID1A",       6,     6,     6,     8,     2,    38,     8,    24,    13,    14,
  "ATM",          7,     6,     7,     3,     6,     5,     5,     8,     5,     4,
  "BRCA2",        4,     4,     4,     6,     5,     4,     3,     6,     4,     3,
  "BRCA1",        3,     3,     1,     8,     1,     2,     2,     5,     1,     1,
  "FBXW7",        1,     1,    11,     1,     0,    13,     1,     5,     3,     1,
  "KMT2D",        5,     4,     8,     3,     2,    14,     2,    12,     4,     3,
  "ERBB2",        4,     6,     5,     4,     5,     5,     2,    10,     6,     4,
  "RET",          2,     1,     1,     1,     1,     2,     1,     2,     2,     1,
  "KIT",          2,     1,     2,     1,     1,     1,     1,     2,     2,     2,
  "FGFR3",        1,     1,     1,     1,     1,     2,     1,    22,     1,     1,
  "MET",          4,     2,     2,     2,     2,     3,     1,     3,     2,     2,
  "ALK",          4,     2,     3,     1,     1,     3,     1,     3,     3,     2,
  "BRAF",         5,     1,    11,     4,     1,     7,     1,     1,     2,     2,
  "CTNNB1",       1,     1,     1,     1,     1,     4,     1,     2,     1,    16,
  "EGFR",        20,     2,     1,     1,     1,     2,     1,     2,     2,     1,
  "POLE",         2,     1,     3,     2,     1,     8,     0,     3,     1,     1,
  "CREBBP",       4,     4,     4,     2,     2,     5,     1,    14,     5,     1,
  "ERBB4",        4,     1,     4,     1,     1,     3,     2,     3,     6,     2,
  "NRAS",         1,     0,     4,     1,     0,     2,     0,     1,     0,     2,
  "NOTCH1",       4,     3,     3,     4,     2,     5,     2,     5,     7,     2,
  "MSH6",         1,     1,     5,     1,     1,     8,     1,     4,     1,     1,
  "MSH2",         1,     1,     4,     1,     1,     5,     1,     3,     1,     1,
  "MLH1",         1,     1,     3,     1,     1,     5,     1,     2,     1,     1,
  "PALB2",        1,     1,     2,     1,     2,     2,     1,     3,     1,     1,
  "CHEK2",        1,     2,     2,     1,     2,     2,     1,     3,     1,     1,
  "FGFR2",        2,     2,     2,     2,     1,     5,     1,     3,     2,     1,
  "ESR1",         1,     1,     0,     0,     0,     0,     0,     0,     0,     0,
  "SMARCA4",      7,     2,     3,     2,     1,     5,     2,     7,     6,     2
)

# Convert to long format
genie_nhb_long <- genie_nhb %>%
  pivot_longer(-gene, names_to = "cancer_type", values_to = "nhb_freq_pct") %>%
  mutate(nhb_freq = nhb_freq_pct / 100)

genie_nhw_long <- genie_nhw %>%
  pivot_longer(-gene, names_to = "cancer_type", values_to = "nhw_freq_pct") %>%
  mutate(nhw_freq = nhw_freq_pct / 100)

# Join NHB and NHW
nhb_vs_nhw <- genie_nhb_long %>%
  inner_join(genie_nhw_long, by = c("gene", "cancer_type")) %>%
  mutate(
    freq_diff_pct = nhb_freq_pct - nhw_freq_pct,
    freq_diff = nhb_freq - nhw_freq,
    fold_change = ifelse(nhw_freq > 0, nhb_freq / nhw_freq, NA)
  )

# Summary by gene (mean across cancer types)
nhb_nhw_summary <- nhb_vs_nhw %>%
  group_by(gene) %>%
  summarise(
    mean_nhb_pct = mean(nhb_freq_pct, na.rm = TRUE),
    mean_nhw_pct = mean(nhw_freq_pct, na.rm = TRUE),
    mean_diff_pct = mean(freq_diff_pct, na.rm = TRUE),
    n_cancer_types = n(),
    n_higher_nhb = sum(freq_diff_pct > 1),
    n_higher_nhw = sum(freq_diff_pct < -1),
    .groups = "drop"
  ) %>%
  arrange(desc(abs(mean_diff_pct)))

cat("NHB vs NHW frequency comparison (mean across cancer types):\n\n")

cat("Genes with higher mean NHB frequency (>1% difference):\n")
print(nhb_nhw_summary %>% filter(mean_diff_pct > 1) %>%
        select(gene, mean_nhb_pct, mean_nhw_pct, mean_diff_pct) %>%
        arrange(desc(mean_diff_pct)))

cat("\nGenes with higher mean NHW frequency (>1% difference):\n")
print(nhb_nhw_summary %>% filter(mean_diff_pct < -1) %>%
        select(gene, mean_nhb_pct, mean_nhw_pct, mean_diff_pct) %>%
        arrange(mean_diff_pct))

# Save
write_csv(nhb_vs_nhw, here("results", "tables", "sensitivity_genie_nhb_vs_nhw.csv"))
write_csv(nhb_nhw_summary, here("results", "tables", "sensitivity_genie_nhb_nhw_summary.csv"))

# Calculate mean NHB frequency per gene for re-ranking
genie_nhb_mean <- genie_nhb_long %>%
  group_by(gene) %>%
  summarise(
    mean_nhb_freq = mean(nhb_freq, na.rm = TRUE),
    max_nhb_freq = max(nhb_freq, na.rm = TRUE),
    .groups = "drop"
  )

# Heatmap: NHB - NHW frequency difference
nhb_nhw_for_heat <- nhb_vs_nhw %>%
  filter(gene %in% panel_70$gene) %>%
  select(gene, cancer_type, freq_diff_pct)

# Order genes by mean absolute difference
gene_order_heat <- nhb_nhw_for_heat %>%
  group_by(gene) %>%
  summarise(mean_abs_diff = mean(abs(freq_diff_pct))) %>%
  arrange(desc(mean_abs_diff)) %>%
  pull(gene)

nhb_nhw_for_heat <- nhb_nhw_for_heat %>%
  mutate(gene = factor(gene, levels = rev(gene_order_heat)))

p_nhb_heatmap <- ggplot(nhb_nhw_for_heat, aes(x = cancer_type, y = gene, fill = freq_diff_pct)) +
  geom_tile(color = "white", size = 0.3) +
  scale_fill_gradient2(
    low = "#4DBBD5", mid = "white", high = "#E64B35",
    midpoint = 0,
    name = "NHB - NHW\nFreq Diff (%)",
    limits = c(-25, 35)
  ) +
  labs(
    title = "Mutation Frequency Differences: NHB vs NHW (GENIE v13)",
    subtitle = "Positive (red) = higher in NHB; Negative (blue) = higher in NHW",
    x = "Cancer Type", y = "Gene"
  ) +
  theme_minimal(base_size = 10) +
  theme(
    axis.text.x = element_text(angle = 45, hjust = 1, size = 9),
    axis.text.y = element_text(size = 7),
    legend.position = "right"
  )

ggsave(here("results", "figures", "fig_supp_nhb_nhw_heatmap.pdf"),
       p_nhb_heatmap, width = 12, height = 14, dpi = 300)

cat("\nSaved: fig_supp_nhb_nhw_heatmap.pdf\n\n")


# SECTION 4: GENE RE-RANKING WITH AFRICAN-ANCESTRY FREQUENCIES ----

cat("SECTION 4: Gene re-ranking with African-ancestry frequencies\n")
cat(strrep("-", 60), "\n\n")

# Merge Black frequencies into integrated scores
rescored <- integrated_scores %>%
  left_join(
    pancancer_black %>% select(gene, pancancer_freq_black),
    by = "gene"
  ) %>%
  left_join(
    genie_nhb_mean %>% select(gene, mean_nhb_freq),
    by = "gene"
  ) %>%
  mutate(
    # Use average of TCGA-Black and GENIE-NHB where both available
    african_ancestry_freq = case_when(
      !is.na(pancancer_freq_black) & !is.na(mean_nhb_freq) ~
        (pancancer_freq_black + mean_nhb_freq) / 2,
      !is.na(pancancer_freq_black) ~ pancancer_freq_black,
      !is.na(mean_nhb_freq) ~ mean_nhb_freq,
      TRUE ~ NA_real_
    ),
    
    # Recalculate African score using actual ancestry data
    african_score_revised = percent_rank(coalesce(african_ancestry_freq,
                                                  pancancer_frequency)),
    
    # Recalculate IPS
    priority_score_revised = (
      combined_score * 0.3 +
        african_score_revised * 0.3 +
        literature_score * 0.2 +
        curated_score * 0.2
    )
  ) %>%
  arrange(desc(priority_score_revised))

# Compare rankings for panel-70 genes
rank_comparison <- rescored %>%
  mutate(
    original_rank = rank(-priority_score, ties.method = "first"),
    revised_rank = rank(-priority_score_revised, ties.method = "first")
  ) %>%
  filter(gene %in% panel_70$gene) %>%
  mutate(rank_change = original_rank - revised_rank) %>%
  select(gene, original_rank, revised_rank, rank_change,
         priority_score, priority_score_revised) %>%
  arrange(revised_rank)

# Rank correlation
cor_ranks <- cor.test(rank_comparison$original_rank,
                      rank_comparison$revised_rank,
                      method = "spearman")

cat("Gene ranking stability (Panel-70 genes):\n")
cat("  Spearman rho:", round(cor_ranks$estimate, 3), "\n")
cat("  p-value:", format(cor_ranks$p.value, digits = 3), "\n")
cat("  Genes improving >5 ranks:", sum(rank_comparison$rank_change > 5), "\n")
cat("  Genes declining >5 ranks:", sum(rank_comparison$rank_change < -5), "\n")
cat("  Genes stable (within 5 ranks):",
    sum(abs(rank_comparison$rank_change) <= 5), "\n\n")

write_csv(rank_comparison,
          here("results", "tables", "sensitivity_rank_comparison.csv"))

# Top movers
cat("Genes with largest rank improvements (higher in African-ancestry):\n")
print(rank_comparison %>% arrange(desc(rank_change)) %>% head(10) %>%
        select(gene, original_rank, revised_rank, rank_change))

cat("\nGenes with largest rank declines (lower in African-ancestry):\n")
print(rank_comparison %>% arrange(rank_change) %>% head(10) %>%
        select(gene, original_rank, revised_rank, rank_change))



# SECTION 5: PANEL COVERAGE IN BLACK PATIENTS ----


cat("\n\nSECTION 5: Panel coverage in Black patients only\n")
cat(strrep("-", 60), "\n\n")

if (nrow(tcga_mutations_black) > 0) {
  
  panel_gene_list <- panel_70$gene
  
  # Calculate coverage by cancer type
  coverage_black <- tcga_mutations_black %>%
    mutate(in_panel = gene %in% panel_gene_list) %>%
    group_by(study_id, cancer_type) %>%
    summarise(
      total_samples = n_distinct(sample_id),
      samples_with_panel_mutation = n_distinct(sample_id[in_panel]),
      sample_coverage_pct = (samples_with_panel_mutation / total_samples) * 100,
      .groups = "drop"
    ) %>%
    filter(total_samples >= 10) %>%
    arrange(desc(sample_coverage_pct))
  
  cat("Panel-70 coverage in Black patients (cancer types with >=10 patients):\n")
  print(coverage_black)
  cat("\nMean sample coverage (Black patients):",
      round(mean(coverage_black$sample_coverage_pct), 1), "%\n\n")
  
  write_csv(coverage_black,
            here("results", "tables", "sensitivity_coverage_black_patients.csv"))
} else {
  cat("Insufficient Black patient data for coverage analysis.\n\n")
}



# SECTION 6: WEIGHTING SCHEME SENSITIVITY ----


cat("SECTION 6: Weighting scheme sensitivity analysis\n")
cat(strrep("-", 60), "\n\n")

weighting_schemes <- list(
  original           = c(0.30, 0.30, 0.20, 0.20),
  equal              = c(0.25, 0.25, 0.25, 0.25),
  african_dominant   = c(0.20, 0.50, 0.15, 0.15),
  frequency_dominant = c(0.50, 0.20, 0.15, 0.15),
  no_literature      = c(0.35, 0.35, 0.00, 0.30),
  no_curation        = c(0.35, 0.35, 0.30, 0.00)
)

panel_70_genes <- panel_70$gene

sensitivity_results <- map_dfr(names(weighting_schemes), function(scheme_name) {
  w <- weighting_schemes[[scheme_name]]
  
  rescored_scheme <- integrated_scores %>%
    mutate(
      new_priority = (combined_score * w[1]) +
        (african_score * w[2]) +
        (literature_score * w[3]) +
        (curated_score * w[4])
    ) %>%
    arrange(desc(new_priority))
  
  top_70 <- rescored_scheme %>% slice_head(n = 70) %>% pull(gene)
  
  # Jaccard similarity with original
  intersection_n <- length(base::intersect(top_70, panel_70_genes))
  union_n <- length(base::union(top_70, panel_70_genes))
  jaccard <- intersection_n / union_n
  
  tibble(
    scheme = scheme_name,
    weights = paste(w, collapse = "/"),
    genes_retained = intersection_n,
    genes_new = 70 - intersection_n,
    jaccard_similarity = round(jaccard, 3)
  )
})

cat("Weighting Sensitivity Analysis Results:\n\n")
print(sensitivity_results)
cat("\n")

write_csv(sensitivity_results,
          here("results", "tables", "sensitivity_weighting_schemes.csv"))



# SECTION 7: PANEL TIER COMPARISON TABLE ----


cat("SECTION 7: Quantitative comparison across panel tiers\n")
cat(strrep("-", 60), "\n\n")

panel_30 <- read_csv(here("results", "tables", "panel_30_genes.csv"),
                     show_col_types = FALSE)
panel_130 <- read_csv(here("results", "tables", "panel_130_genes.csv"),
                      show_col_types = FALSE)

# Calculate comparison metrics
tier_comparison <- tibble(
  Panel = c("30-gene", "70-gene", "130-gene"),
  Size = c(30, 70, 130),
  Mean_Priority_Score = c(
    mean(integrated_scores$priority_score[integrated_scores$gene %in% panel_30$gene], na.rm = TRUE),
    mean(integrated_scores$priority_score[integrated_scores$gene %in% panel_70$gene], na.rm = TRUE),
    mean(integrated_scores$priority_score[integrated_scores$gene %in% panel_130$gene], na.rm = TRUE)
  ),
  SD_Priority_Score = c(
    sd(integrated_scores$priority_score[integrated_scores$gene %in% panel_30$gene], na.rm = TRUE),
    sd(integrated_scores$priority_score[integrated_scores$gene %in% panel_70$gene], na.rm = TRUE),
    sd(integrated_scores$priority_score[integrated_scores$gene %in% panel_130$gene], na.rm = TRUE)
  ),
  Mean_African_Score = c(
    mean(integrated_scores$african_score[integrated_scores$gene %in% panel_30$gene], na.rm = TRUE),
    mean(integrated_scores$african_score[integrated_scores$gene %in% panel_70$gene], na.rm = TRUE),
    mean(integrated_scores$african_score[integrated_scores$gene %in% panel_130$gene], na.rm = TRUE)
  ),
  SD_African_Score = c(
    sd(integrated_scores$african_score[integrated_scores$gene %in% panel_30$gene], na.rm = TRUE),
    sd(integrated_scores$african_score[integrated_scores$gene %in% panel_70$gene], na.rm = TRUE),
    sd(integrated_scores$african_score[integrated_scores$gene %in% panel_130$gene], na.rm = TRUE)
  ),
  Actionable_Genes = c(
    sum(panel_30$is_actionable, na.rm = TRUE),
    sum(panel_70$is_actionable, na.rm = TRUE),
    sum(panel_130$is_actionable, na.rm = TRUE)
  ),
  Pct_Actionable = round(Actionable_Genes / Size * 100, 1),
  Mean_Pancancer_Freq = c(
    mean(panel_30$pancancer_frequency, na.rm = TRUE),
    mean(panel_70$pancancer_frequency, na.rm = TRUE),
    mean(panel_130$pancancer_frequency, na.rm = TRUE)
  ),
  SD_Pancancer_Freq = c(
    sd(panel_30$pancancer_frequency, na.rm = TRUE),
    sd(panel_70$pancancer_frequency, na.rm = TRUE),
    sd(panel_130$pancancer_frequency, na.rm = TRUE)
  ),
  Estimated_Cost = c("$300", "$500", "$800")
)

cat("Panel Tier Comparison:\n")
print(tier_comparison)
cat("\n")

write_csv(tier_comparison,
          here("results", "tables", "sensitivity_panel_tier_comparison.csv"))

