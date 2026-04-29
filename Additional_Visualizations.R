library(dplyr)
library(ggplot2)
library(tidyr)

### Supplementary Table ----

Gene_Scores <- comprehensive_gene_database_top15 %>%
  select(gene, clinical_relevance, tier) %>%
  left_join(categorized_scores %>% 
              select(gene, mean_mutation_freq, priority_score), 
            by = "gene")

### Remove genes without priority scores
Gene_Scores <- Gene_Scores %>%
  filter(!is.na(priority_score))

### Remove duplicate gene rows
Gene_Scores <- Gene_Scores %>%
  distinct(gene, .keep_all = TRUE)

Gene_Scores <- Gene_Scores %>%
  select(gene, mean_mutation_freq, priority_score, clinical_relevance, tier)

write.csv(Gene_Scores, "Gene_Scores.csv", row.names = FALSE)




#### Figure 1b ----

tcga_mutation_frequencies <- tcga_mutation_frequencies %>%
  arrange(desc(total_samples)) %>%
  mutate(cancer_type = factor(cancer_type, levels = unique(cancer_type)))

### Filter unique total samples
tcga_frequencies <- tcga_mutation_frequencies %>%
  distinct(cancer_type, .keep_all = TRUE)


### Make the Plot 
ggplot(tcga_frequencies) +
  aes(x = cancer_type, y = total_samples) +
  geom_col(fill = "blue") +
  geom_text(aes(label = total_samples), 
            vjust = -0.5, 
            size = 3) +
  scale_y_continuous(expand = expansion(mult = c(0, 0.1))) +  # Adds space for labels
  theme_minimal() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
  labs(x = "Cancer Type", y = "Total Samples")

#### Figure 1c ----
## FIlter genes to top 20 and arrange by literature count
literature_frequency <- gene_literature_frequency %>%
  arrange(desc(literature_count)) %>%
  slice_head(n = 20) %>%
  mutate(genes = factor(genes, levels = unique(genes)))


## Make the plot
ggplot(literature_frequency) +
  aes(x = genes, y = literature_count) +
  geom_col(fill = "blue") +
  geom_text(aes(label = literature_count), 
            vjust = -0.5, 
            size = 3) +
  scale_y_continuous(expand = expansion(mult = c(0, 0.1))) +  # Adds space for labels
  theme_minimal() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
  labs(x = "Genes", y = "Literature Count")

#### Figure 2b ----
# Filter genes to top 20 and arrange by priority score
priority_score <- integrated_gene_scores %>%
  arrange(desc(priority_score)) %>%
  slice_head(n = 30) %>%
  mutate(genes = factor(gene, levels = unique(gene)))

## Make the plot
ggplot(priority_score) +
  aes(reorder(gene, -priority_score), y = priority_score) +
  geom_col(fill = "blue") +
  #geom_text(aes(label = priority_score), 
   #         vjust = -0.5, 
   #         size = 3) +
  scale_y_continuous(expand = expansion(mult = c(0, 0.1))) +  # Adds space for labels
  theme_minimal() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
  labs(x = "Genes", y = "Gene Priority Score")

#### Figure 2c ----
# Filter genes to top 20 and arrange by clinical_relevance
clinical_score <- comprehensive_gene_database_top15 %>%
  arrange(desc(clinical_relevance)) %>%
  slice_head(n = 30) %>%
  mutate(genes = factor(gene, levels = unique(gene)))

## Make the plot
ggplot(clinical_score) +
  aes(reorder(gene, -clinical_relevance), y = clinical_relevance) +
  geom_col(fill = "blue") +
  #geom_text(aes(label = priority_score), 
  #         vjust = -0.5, 
  #         size = 3) +
  scale_y_continuous(expand = expansion(mult = c(0, 0.1))) +  # Adds space for labels
  theme_minimal() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
  labs(x = "Genes", y = "Gene Clinical Relevance Score")

#### Figure 3a ----

# Find the 80% coverage inflection point
inflection_80 <- panel_size_elbow_analysis %>%
  filter(mean_coverage > 81) %>%
  arrange(panel_size) %>%
  filter(panel_size %% 10 == 0) %>%
  slice_head(n = 1)

if (nrow(inflection_80) == 0) {
  inflection_80 <- panel_size_elbow_analysis %>%
    filter(mean_coverage > 81) %>%
    arrange(panel_size) %>%
    slice_head(n = 1)
}

# Calculate elbow point from marginal benefit (as in your original code)
elbow_analysis_with_marginal <- panel_size_elbow_analysis %>%
  arrange(panel_size) %>%
  mutate(
    marginal_coverage = c(NA, diff(mean_coverage)),
    marginal_coverage_per_gene = marginal_coverage / 10  # Assuming increments of 10
  )

# Find elbow point using your threshold method
elbow_threshold <- max(elbow_analysis_with_marginal$marginal_coverage_per_gene, na.rm = TRUE) * 0.2
elbow_point_data <- elbow_analysis_with_marginal %>%
  filter(marginal_coverage_per_gene < elbow_threshold) %>%
  slice_head(n = 1)

if (nrow(elbow_point_data) == 0) {
  elbow_point_data <- panel_size_elbow_analysis %>%
    filter(panel_size == 80) %>%  # Default as in your code
    mutate(marginal_coverage = NA, marginal_coverage_per_gene = NA)
}

# Determine y-axis maximum
y_max <- max(panel_size_elbow_analysis$mean_coverage, na.rm = TRUE)
y_max <- ceiling(y_max / 20) * 20

# Create the plot
ggplot(panel_size_elbow_analysis) +
  aes(x = panel_size, y = mean_coverage) +
  geom_line(colour = "#112446") +
  geom_point(color = "#112446") +
  
  # Highlight 80% inflection point (in blue)
  geom_point(data = inflection_80,
             color = "blue", size = 4, shape = 21, fill = "white", stroke = 1.5) +
  geom_vline(data = inflection_80,
             aes(xintercept = panel_size),
             linetype = "dashed", color = "blue", alpha = 0.5) +
  
  # Highlight elbow point (in red)
  geom_point(data = elbow_point_data,
             color = "red", size = 4, shape = 24, fill = "white", stroke = 1.5) +
  geom_vline(data = elbow_point_data,
             aes(xintercept = panel_size),
             linetype = "dashed", color = "red", alpha = 0.5) +
  
  # Add horizontal line at 80%
  geom_hline(yintercept = 80, linetype = "dashed", color = "gray", alpha = 0.7) +
  
  # Add annotations
  geom_text(data = inflection_80,
            aes(x = panel_size, y = mean_coverage,
                label = paste0("80% Coverage\n", panel_size, " genes")),
            hjust = -0.1, vjust = -0.5, size = 3, color = "blue") +
  
  geom_text(data = elbow_point_data,
            aes(x = panel_size, y = mean_coverage,
                label = paste0("Elbow Point\n", panel_size, " genes")),
            hjust = -0.1, vjust = 1.5, size = 3, color = "red") +
  
  # Y-axis with intervals of 20
  scale_y_continuous(
    limits = c(0, y_max),
    breaks = seq(0, y_max, by = 20),
    labels = function(x) paste0(x, "%")
  ) +
  
  theme_minimal() +
  labs(title = "Elbow Analysis: Coverage by Panel Size",
       subtitle = paste("Blue: First panel size >80% coverage",
                        "Red: Elbow point (marginal benefit <20% of max)"),
       x = "Panel Size",
       y = "Mean Coverage (%)",
       caption = paste("Elbow threshold:", round(elbow_threshold, 3), 
                       "% per gene increase"))

#### Figure 3b ----

# Create combined plot with coverage and marginal benefit
coverage_plot <- ggplot(panel_size_elbow_analysis) +
  aes(x = panel_size, y = mean_coverage) +
  geom_line(colour = "#112446") +
  geom_point(color = "#112446") +
  geom_point(data = inflection_80, color = "blue", size = 4) +
  geom_point(data = elbow_point_data, color = "red", size = 4) +
  geom_hline(yintercept = 80, linetype = "dashed", color = "gray") +
  scale_y_continuous(
    limits = c(0, y_max),
    breaks = seq(0, y_max, by = 20),
    labels = function(x) paste0(x, "%"),
    name = "Coverage (%)"
  ) +
  labs(title = "Coverage Elbow Analysis",
       x = "Panel Size") +
  theme_minimal()

# Create marginal benefit plot
marginal_plot <- ggplot(elbow_analysis_with_marginal) +
  aes(x = panel_size, y = marginal_coverage_per_gene) +
  geom_line(color = "darkgreen") +
  geom_point(color = "darkgreen") +
  geom_hline(yintercept = elbow_threshold, 
             linetype = "dashed", color = "red") +
  geom_vline(data = elbow_point_data,
             aes(xintercept = panel_size),
             linetype = "dashed", color = "red", alpha = 0.5) +
  labs(title = "Marginal Benefit per Additional 10 Genes",
       x = "Panel Size",
       y = "Coverage Increase per Gene (%)") +
  theme_minimal()

# Combine plots
library(patchwork)
coverage_plot / marginal_plot +
  plot_annotation(title = "Complete Elbow Analysis",
                  subtitle = "Top: Coverage curve | Bottom: Marginal benefit with elbow threshold")

#### Figure 3c ----
# Calculate the percentage and keep both values
panel_size_elbow_analysis <- panel_size_elbow_analysis %>%
  mutate(
    n_actionable_pct = (n_actionable / max(n_actionable, na.rm = TRUE)) * 100,
    max_actionable = max(n_actionable, na.rm = TRUE)
  )

# Find inflection point
inflection_point <- panel_size_elbow_analysis %>%
  filter(n_actionable_pct > 80) %>%
  arrange(panel_size) %>%
  filter(panel_size %% 10 == 0) %>%
  slice_head(n = 1)

if (nrow(inflection_point) == 0) {
  inflection_point <- panel_size_elbow_analysis %>%
    filter(n_actionable_pct > 80) %>%
    arrange(panel_size) %>%
    slice_head(n = 1)
}

# Get the maximum value for y-axis scaling
max_raw <- max(panel_size_elbow_analysis$n_actionable, na.rm = TRUE)
y_max_pct <- max(panel_size_elbow_analysis$n_actionable_pct, na.rm = TRUE)
y_max_pct <- ceiling(y_max_pct / 20) * 20

ggplot(panel_size_elbow_analysis) +
  aes(x = panel_size, y = n_actionable_pct) +
  geom_line(colour = "#112446") +
  geom_point(color = "#112446") +
  geom_point(data = inflection_point,
             color = "red", size = 4, shape = 21, fill = "white", stroke = 1.5) +
  geom_hline(yintercept = 80, linetype = "dashed", color = "gray", alpha = 0.7) +
  geom_vline(data = inflection_point,
             aes(xintercept = panel_size),
             linetype = "dashed", color = "red", alpha = 0.7) +
  # Annotation showing both percentage and raw value
  geom_text(data = inflection_point,
            aes(x = panel_size, y = n_actionable_pct,
                label = paste0("Panel: ", panel_size, 
                               "\n", round(n_actionable_pct, 1), "% of max",
                               "\n(", n_actionable, " / ", max_actionable, ")")),
            hjust = -0.1, vjust = -0.5, size = 3) +
  scale_y_continuous(
    limits = c(0, y_max_pct),
    breaks = seq(0, y_max_pct, by = 20),
    labels = function(x) paste0(x, "%"),
    # Add secondary axis showing raw values
    sec.axis = sec_axis(
      ~ . * max_raw / 100,
      name = "Number of Actionable Variants",
      breaks = scales::pretty_breaks(n = 6)
    )
  ) +
  theme_minimal() +
  labs(title = "Elbow Analysis: Actionable Variants by Panel Size",
       subtitle = paste("Maximum actionable variants:", max_raw),
       x = "Panel Size",
       y = "Actionable Variants (% of maximum)")


#### Figure 3d ----

# Sort data and find top 2 highest efficiency ratios
panel_size_cost_effectiveness <- panel_size_cost_effectiveness %>%
  arrange(desc(efficiency_ratio))

top_2 <- panel_size_cost_effectiveness %>%
  slice_max(efficiency_ratio, n = 2)

ggplot(panel_size_cost_effectiveness) +
  aes(x = panel_size, y = efficiency_ratio) +
  geom_line(colour = "#112446") +
  geom_point(color = "#112446") +
  # Highlight top 2 points in red
  geom_point(data = top_2, 
             color = "red", 
             size = 3) +
  # Add labels for top 2 points
  geom_text(data = top_2,
            aes(label = paste0("Size: ", panel_size, "\nRatio: ", 
                               round(efficiency_ratio, 2))),
            hjust = -0.1, 
            vjust = 0.5,
            size = 3,
            color = "red") +
  # Start y-axis from 0
  scale_y_continuous(limits = c(0, NA), 
                     expand = expansion(mult = c(0, 0.1))) +
  theme_minimal() +
  labs(title = "Cost Efficiency Ratio by Panel Size",
       subtitle = "Red points highlight top 2 highest efficiency ratios",
       x = "Panel Size",
       y = "Cost Efficiency Ratio")

#### Figure 6b ----
# Boxplot by gene
ggplot(simulation_results) +
  aes(x = reorder(gene, mutation_rate, median), y = mutation_rate) +
  geom_boxplot(fill = "steelblue", alpha = 0.7) +
  coord_flip() +  # Flip for better gene name readability
  theme_minimal() +
  labs(title = "Mutation Rate Distribution Across Simulations",
       x = "Gene",
       y = "Mutation Rate")

# Violin plot for better distribution visualization
ggplot(simulation_results) +
  aes(x = reorder(gene, mutation_rate, median), y = mutation_rate) +
  geom_violin(fill = "steelblue", alpha = 0.7) +
  geom_boxplot(width = 0.1, fill = "white", alpha = 0.7) +
  coord_flip() +
  theme_minimal() +
  labs(title = "Mutation Rate Distribution by Gene",
       subtitle = "Violin shows distribution, box shows quartiles")

### Top genes heatmap
# Calculate average mutation rate per gene
gene_avg <- simulation_results %>%
  group_by(gene) %>%
  summarise(
    avg_rate = mean(mutation_rate),
    avg_mutated = mean(n_mutated)
  ) %>%
  arrange(desc(avg_rate)) %>%
  slice_head(n = 20)  # Top 20 genes

# Filter for top genes
top_genes_data <- simulation_results %>%
  filter(gene %in% gene_avg$gene) %>%
  mutate(gene = factor(gene, levels = gene_avg$gene))

# Heatmap
ggplot(top_genes_data) +
  aes(x = simulation, y = gene, fill = mutation_rate) +
  geom_tile() +
  scale_fill_viridis_c(option = "plasma", name = "Mutation Rate") +
  theme_minimal() +
  labs(title = "Mutation Rate Heatmap (Top 20 Genes)",
       x = "Simulation",
       y = "Gene")

## Faceted Histograms
# Top 12 genes faceted
top_12_genes <- simulation_results %>%
  group_by(gene) %>%
  summarise(avg_rate = mean(mutation_rate)) %>%
  arrange(desc(avg_rate)) %>%
  slice_head(n = 12) %>%
  pull(gene)

simulation_results %>%
  filter(gene %in% top_12_genes) %>%
  ggplot() +
  aes(x = mutation_rate) +
  geom_histogram(bins = 30, fill = "steelblue", alpha = 0.7) +
  facet_wrap(~ gene, scales = "free_y", ncol = 3) +
  theme_minimal() +
  labs(title = "Mutation Rate Distribution for Top 12 Genes",
       x = "Mutation Rate",
       y = "Frequency")

## Simulated Trends line plot
# Plot mutation rate trends across simulations for top genes
top_5_genes <- simulation_results %>%
  group_by(gene) %>%
  summarise(avg_rate = mean(mutation_rate)) %>%
  arrange(desc(avg_rate)) %>%
  slice_head(n = 5) %>%
  pull(gene)

simulation_results %>%
  filter(gene %in% top_5_genes) %>%
  ggplot() +
  aes(x = simulation, y = mutation_rate, color = gene, group = gene) +
  geom_line(alpha = 0.7) +
  geom_point(size = 0.5) +
  theme_minimal() +
  labs(title = "Mutation Rate Trends Across Simulations",
       x = "Simulation",
       y = "Mutation Rate",
       color = "Gene")

### Supplementary table Article and Gene List ----
# Create Supplementary Table 6: Complete article list with PMIDs and genes

library(tidyverse)
library(here)

# Load data

articles <- read_csv(
  here("data", "processed", "pubmed_african_cancer_articles.csv"),
  show_col_types = FALSE
)

genes_lit <- read_csv(
  here("data", "processed", "genes_in_literature.csv"),
  show_col_types = FALSE
)

cat("Articles loaded:", nrow(articles), "\n")
cat("Gene-article rows loaded:", nrow(genes_lit), "\n\n")

# Collapse genes per PMID into one row

# Each PMID can have multiple gene rows — collapse to comma-separated
genes_collapsed <- genes_lit %>%
  filter(!is.na(genes), genes != "") %>%
  group_by(pmid) %>%
  summarise(
    genes_detected = paste(sort(unique(genes)), collapse = ", "),
    n_genes = n_distinct(genes),
    .groups = "drop"
  )

cat("Unique PMIDs with genes:", nrow(genes_collapsed), "\n")

# Extract first author from authors string

# The authors column uses "; " as separator
# First author format: "Surname Initials" or "Surname I; ..."
# We want: "Surname et al" (or just the name if single author)

articles_processed <- articles %>%
  mutate(
    # Extract the first author (everything before the first ";")
    first_author_raw = str_trim(str_extract(authors, "^[^;]+")),
    
    # Extract surname: take everything before the last space 
    # (which would be initials like "C" or "DW")
    # For names like "Christowitz C", surname = "Christowitz"
    # For names like "van der Merwe N", we need the part before the last token
    first_author_surname = str_trim(
      str_replace(first_author_raw, "\\s+[A-Z]+$", "")
    ),
    
    # Count number of authors (number of ";" separators + 1)
    n_authors = str_count(authors, ";") + 1L,
    
    # Format: "Surname et al" if >1 author, just "Surname" if solo
    first_author = if_else(
      n_authors > 1,
      paste0(first_author_surname, " et al"),
      first_author_surname
    )
  )

# Merge articles with collapsed genes

supp_table_6 <- articles_processed %>%
  left_join(genes_collapsed, by = "pmid") %>%
  mutate(
    # Articles without detected genes get "None detected"
    genes_detected = if_else(is.na(genes_detected), "None detected", genes_detected),
    n_genes = if_else(is.na(n_genes), 0L, n_genes)
  ) %>%
  select(
    PMID = pmid,
    First_Author = first_author,
    Journal = journal,
    Cancer_Type = cancer_type,
    Genes = genes_detected,
    N_Genes = n_genes
  ) %>%
  arrange(Cancer_Type, desc(N_Genes), PMID)

## arrange in descending order of PMID (most recent first) for better readability in the table
supp_table_6 <- supp_table_6 %>%
  arrange(desc(PMID))

# Summary statistics

cat("\n")
cat(strrep("=", 60), "\n")
cat("SUPPLEMENTARY TABLE 6 SUMMARY\n")
cat(strrep("=", 60), "\n\n")

cat("Total articles:", nrow(supp_table_6), "\n")
cat("Articles with ≥1 gene detected:", sum(supp_table_6$N_Genes > 0), "\n")
cat("Articles with no genes detected:", sum(supp_table_6$N_Genes == 0), "\n\n")

cat("By cancer type:\n")
print(
  supp_table_6 %>%
    group_by(Cancer_Type) %>%
    summarise(
      N_Articles = n(),
      Articles_With_Genes = sum(N_Genes > 0),
      Mean_Genes_Per_Article = round(mean(N_Genes), 1),
      .groups = "drop"
    ) %>%
    arrange(desc(N_Articles))
)

cat("\nTop 10 most gene-rich articles:\n")
print(
  supp_table_6 %>%
    arrange(desc(N_Genes)) %>%
    select(PMID, First_Author, Cancer_Type, N_Genes, Genes) %>%
    head(10)
)

# Save

write_csv(
  supp_table_6,
  here("results", "tables", "supplementary_table_6_articles_with_genes.csv")
)

cat("\nSaved: supplementary_table_6_articles_with_genes.csv\n")
cat("Location:", here("results", "tables"), "\n")

# Preview first few rows

cat("\nPreview (first 5 rows):\n")
print(supp_table_6 %>% head(5), width = 120)

### Black Patient Sample Coverage Analysis ----
# Generate Figure 4b - Coverage analysis restricted to Black/African American patients in TCGA
# This script uses locally downloaded TCGA clinical files (data_clinical_patient.txt)
# and filters mutations to self-reported Black individuals.

library(tidyverse)
library(ggplot2)
library(patchwork)
library(here)

# Set theme
theme_set(theme_minimal(base_size = 12))

# 1. Load data

cat("Loading TCGA mutation data...\n")
tcga_mutations <- readRDS(here("data", "processed", "tcga_mutations_top15_africa.rds"))

# Load panel gene lists (generated by 04_panel_optimization.R)
cat("Loading panel gene lists...\n")
panel_files <- list.files(here("results", "tables"), pattern = "^panel_[0-9]+_genes\\.csv$", full.names = TRUE)
if (length(panel_files) == 0) {
  stop("No panel gene files found in results/tables/. Run 04_panel_optimization.R first.")
}

# Read all panels and extract sizes
panels <- list()
for (f in panel_files) {
  size <- as.integer(str_extract(basename(f), "[0-9]+"))
  panels[[as.character(size)]] <- read_csv(f, show_col_types = FALSE)
}
cat("Loaded panels with sizes:", paste(names(panels), collapse = ", "), "\n")

# 2. Extract Black patient IDs from local clinical files

# TCGA raw data directory
tcga_raw_dir <- here("data", "raw", "tcga_maf")

# Find all study directories (those with _tcga_pan_can_atlas_2018 suffix)
study_dirs <- list.dirs(tcga_raw_dir, recursive = FALSE, full.names = TRUE)
study_dirs <- study_dirs[grepl("_tcga_pan_can_atlas_2018$", study_dirs)]

cat("\nProcessing clinical files from", length(study_dirs), "studies...\n")

black_patient_ids <- c()

for (study_dir in study_dirs) {
  # Locate clinical patient file
  clinical_file <- file.path(study_dir, "data_clinical_patient.txt")
  
  if (!file.exists(clinical_file)) {
    cat("  Skipping", basename(study_dir), "- no data_clinical_patient.txt\n")
    next
  }
  
  # Read clinical data
  # The file is tab-delimited with comments starting with '#'
  clinical <- tryCatch({
    read_tsv(clinical_file, comment = "#", show_col_types = FALSE)
  }, error = function(e) {
    cat("  Error reading", basename(study_dir), ":", e$message, "\n")
    return(NULL)
  })
  
  if (is.null(clinical)) next
  
  # Identify race column (could be "race", "patient.race", "RACE", etc.)
  race_col <- NULL
  possible_race <- c("race", "patient.race", "RACE", "RACE_1")
  for (col in possible_race) {
    if (col %in% colnames(clinical)) {
      race_col <- col
      break
    }
  }
  
  if (is.null(race_col)) {
    cat("  Skipping", basename(study_dir), "- no race column found\n")
    next
  }
  
  # Identify patient barcode column
  barcode_col <- NULL
  possible_barcode <- c("patient.bcr_patient_barcode", "bcr_patient_barcode", "patient_id", "PATIENT_ID")
  for (col in possible_barcode) {
    if (col %in% colnames(clinical)) {
      barcode_col <- col
      break
    }
  }
  
  if (is.null(barcode_col)) {
    cat("  Skipping", basename(study_dir), "- no patient barcode column found\n")
    next
  }
  
  # Filter Black/African American patients
  # Common values: "BLACK OR AFRICAN AMERICAN", "black", "BLACK", etc.
  black_patients <- clinical %>%
    filter(tolower(.data[[race_col]]) %in% c("black or african american", "black", "african american")) %>%
    pull(.data[[barcode_col]])
  
  if (length(black_patients) > 0) {
    black_patient_ids <- c(black_patient_ids, black_patients)
    cat("  ", basename(study_dir), ":", length(black_patients), "Black patients\n")
  } else {
    cat("  ", basename(study_dir), ": no Black patients found\n")
  }
}

# Remove possible duplicates (same patient may appear in multiple studies? unlikely but safe)
black_patient_ids <- unique(black_patient_ids)

cat("\nTotal unique Black/African American patients found:", length(black_patient_ids), "\n")

if (length(black_patient_ids) == 0) {
  stop("No Black patient IDs extracted. Check clinical file formats and race column names.")
}


# 3. Filter mutations to Black patients only


cat("\nFiltering mutations to Black patients...\n")

# TCGA sample barcodes: first 12 characters are patient barcode (TCGA-XX-XXXX)
tcga_mutations <- tcga_mutations %>%
  mutate(patient_barcode = str_sub(sample_id, 1, 12))

black_mutations <- tcga_mutations %>%
  filter(patient_barcode %in% black_patient_ids)

cat("Original mutations:", nrow(tcga_mutations), "\n")
cat("Black-only mutations:", nrow(black_mutations), "\n")
cat("Number of unique Black samples:", n_distinct(black_mutations$sample_id), "\n")

if (nrow(black_mutations) == 0) {
  stop("No mutations remain after filtering to Black patients. Check patient ID matching (first 12 chars).")
}


# 4. Recalculate panel coverage for Black patients


calculate_panel_coverage_black <- function(panel_genes, mutations_data) {
  
  panel_gene_list <- panel_genes$gene
  
  # Identify column names
  gene_col <- if ("gene" %in% names(mutations_data)) "gene" else "Hugo_Symbol"
  sample_col <- if ("sample_id" %in% names(mutations_data)) "sample_id" else "Tumor_Sample_Barcode"
  project_col <- if ("study_id" %in% names(mutations_data)) "study_id" else "cancer_type"
  globocan_col <- if ("globocan_rank" %in% names(mutations_data)) "globocan_rank" else NULL
  
  # Calculate coverage
  mutations_data <- mutations_data %>%
    mutate(in_panel = .data[[gene_col]] %in% panel_gene_list)
  
  # By project with GLOBOCAN ranking
  if (!is.null(globocan_col)) {
    coverage_stats <- mutations_data %>%
      group_by(
        project = .data[[project_col]],
        globocan_rank = .data[[globocan_col]]
      ) %>%
      summarise(
        total_mutations = n(),
        panel_mutations = sum(in_panel),
        coverage_percent = (panel_mutations / total_mutations) * 100,
        total_samples = n_distinct(.data[[sample_col]]),
        samples_with_panel_mutation = n_distinct(.data[[sample_col]][in_panel]),
        sample_coverage_percent = (samples_with_panel_mutation / total_samples) * 100,
        .groups = "drop"
      ) %>%
      arrange(globocan_rank)
  } else {
    coverage_stats <- mutations_data %>%
      group_by(project = .data[[project_col]]) %>%
      summarise(
        total_mutations = n(),
        panel_mutations = sum(in_panel),
        coverage_percent = (panel_mutations / total_mutations) * 100,
        total_samples = n_distinct(.data[[sample_col]]),
        samples_with_panel_mutation = n_distinct(.data[[sample_col]][in_panel]),
        sample_coverage_percent = (samples_with_panel_mutation / total_samples) * 100,
        .groups = "drop"
      )
  }
  
  return(coverage_stats)
}

# Calculate coverage for each panel size
coverage_black_list <- list()
for (size in names(panels)) {
  cat("Calculating coverage for", size, "-gene panel on Black patients...\n")
  coverage_black_list[[size]] <- calculate_panel_coverage_black(panels[[size]], black_mutations) %>%
    mutate(panel_size = as.integer(size))
}

coverage_black <- bind_rows(coverage_black_list)


# 5. Create plots similar to Figure 4 (coverage analysis)


# Re-annotate cancer types for pretty labels
coverage_black <- coverage_black %>%
  mutate(
    cancer_label = case_when(
      grepl("brca", project, ignore.case = TRUE)      ~ "Breast (BRCA)",
      grepl("prad", project, ignore.case = TRUE)      ~ "Prostate (PRAD)",
      grepl("cesc", project, ignore.case = TRUE)      ~ "Cervical (CESC)",
      grepl("lihc", project, ignore.case = TRUE)      ~ "Liver (LIHC)",
      grepl("coadread", project, ignore.case = TRUE)  ~ "Colorectal (COAD/READ)",
      grepl("luad", project, ignore.case = TRUE)      ~ "Lung Adeno (LUAD)",
      grepl("lusc", project, ignore.case = TRUE)      ~ "Lung Squam (LUSC)",
      grepl("ov_", project, ignore.case = TRUE)       ~ "Ovarian (OV)",
      grepl("blca", project, ignore.case = TRUE)      ~ "Bladder (BLCA)",
      grepl("stad", project, ignore.case = TRUE)      ~ "Stomach (STAD)",
      grepl("esca", project, ignore.case = TRUE)      ~ "Esophageal (ESCA)",
      grepl("ucec", project, ignore.case = TRUE)      ~ "Endometrial (UCEC)",
      grepl("paad", project, ignore.case = TRUE)      ~ "Pancreatic (PAAD)",
      TRUE ~ project
    ),
    panel_size = as.integer(panel_size)
  )

# A) Mutation coverage by cancer type
p4a_black <- coverage_black %>%
  ggplot(aes(x = reorder(cancer_label, coverage_percent), y = coverage_percent, 
             fill = base::as.factor(panel_size))) +
  geom_col(position = "dodge") +
  coord_flip() +
  scale_fill_nejm(name = "Panel Size") +
  labs(
    title = "A) Mutation Coverage by Cancer Type (Black Patients)",
    x = "Cancer Type",
    y = "Mutation Coverage (%)"
  ) +
  theme(legend.position = "bottom")

# B) Sample coverage by cancer type
p4b_black <- coverage_black %>%
  ggplot(aes(x = reorder(cancer_label, sample_coverage_percent), 
             y = sample_coverage_percent,
             fill = base::as.factor(panel_size))) +
  geom_col(position = "dodge") +
  coord_flip() +
  scale_fill_nejm(name = "Panel Size") +
  labs(
    title = "B) Sample Coverage by Cancer Type (Black Patients)",
    x = "Cancer Type",
    y = "Samples with Panel Mutation (%)"
  ) +
  theme(legend.position = "bottom")

# C) Coverage vs panel size (average across cancers)
p4c_black <- coverage_black %>%
  group_by(panel_size) %>%
  summarise(
    mean_mutation_cov = mean(coverage_percent, na.rm = TRUE),
    se_mutation_cov = sd(coverage_percent, na.rm = TRUE) / sqrt(n()),
    mean_sample_cov = mean(sample_coverage_percent, na.rm = TRUE),
    se_sample_cov = sd(sample_coverage_percent, na.rm = TRUE) / sqrt(n())
  ) %>%
  ggplot(aes(x = panel_size)) +
  geom_line(aes(y = mean_mutation_cov, color = "Mutation Coverage"), size = 1) +
  geom_point(aes(y = mean_mutation_cov, color = "Mutation Coverage"), size = 3) +
  geom_errorbar(aes(ymin = mean_mutation_cov - se_mutation_cov,
                    ymax = mean_mutation_cov + se_mutation_cov,
                    color = "Mutation Coverage"),
                width = 5) +
  geom_line(aes(y = mean_sample_cov, color = "Sample Coverage"), size = 1) +
  geom_point(aes(y = mean_sample_cov, color = "Sample Coverage"), size = 3) +
  geom_errorbar(aes(ymin = mean_sample_cov - se_sample_cov,
                    ymax = mean_sample_cov + se_sample_cov,
                    color = "Sample Coverage"),
                width = 5) +
  scale_color_manual(values = c("#E64B35", "#4DBBD5"), name = "") +
  labs(
    title = "C) Average Coverage vs Panel Size (Black Patients)",
    x = "Panel Size (number of genes)",
    y = "Coverage (%)"
  ) +
  theme(legend.position = "bottom")

# Combine panels into one figure
p4_black <- (p4a_black / p4b_black) | p4c_black
p4_black <- p4_black +
  plot_annotation(
    title = "Figure 4b: Panel Coverage Analysis in Black/African American Patients",
    subtitle = "Coverage restricted to self-reported Black individuals in TCGA",
    theme = theme(plot.title = element_text(size = 16, face = "bold"))
  )


# 6. Save the figure


output_dir <- here("results", "figures")
if (!dir.exists(output_dir)) dir.create(output_dir, recursive = TRUE)

ggsave(
  here(output_dir, "fig4b_coverage_analysis_in_black.pdf"),
  p4_black, width = 15, height = 10, dpi = 300
)

cat("\nFigure saved to:", here(output_dir, "fig4b_coverage_analysis_in_black.pdf"), "\n")

# Save the coverage data for reference
write_csv(coverage_black, here("results", "tables", "coverage_stats_black_patients.csv"))

cat("\nDone.\n")
