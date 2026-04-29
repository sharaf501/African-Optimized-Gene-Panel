# 05_visualization.R
# Create comprehensive visualizations

library(tidyverse)
library(ggplot2)
library(plotly)
library(pheatmap)
library(ComplexHeatmap)
library(circlize)
library(UpSetR)
library(scales)
library(patchwork)
library(ggsci)
library(ggrepel)
library(here)

# Load data
gene_scores <- read_csv(here("data", "processed", "categorized_gene_scores.csv"))
panel_70 <- read_csv(here("results", "tables", "panel_70_genes.csv"))
panel_130 <- read_csv(here("results", "tables", "panel_130_genes.csv"))
panel_30 <- read_csv(here("results", "tables", "panel_30_genes.csv"))
tcga_freq <- read_csv(here("data", "processed", "tcga_mutation_frequencies.csv"))
coverage_stats <- read_csv(here("results", "tables", "panel_coverage_stats.csv"))
cost_benefit <- read_csv(here("results", "tables", "cost_benefit_analysis.csv"))

# Set theme
theme_set(theme_minimal(base_size = 12))

# Color palettes
cancer_colors <- pal_nejm()(8)
tier_colors <- c("#E64B35", "#4DBBD5", "#00A087", "#3C5488", "#F39B7F")


# Figure 1: Gene prioritization overview ----


# Get top 20 African-prioritized genes (ranked by priority score)
top20_african_pri <- gene_scores %>%
  filter(!is.na(priority_score)) %>%
  mutate(
    score_difference = african_score - combined_score,
    in_panel_70 = gene %in% panel_70$gene,
    priority_type = case_when(
      score_difference > 0.05 ~ "African-prioritized",
      score_difference < -0.05 ~ "Pan-cancer-prioritized",
      TRUE ~ "Balanced"
    )
  ) %>%
  filter(priority_type == "African-prioritized") %>%
  arrange(desc(priority_score)) %>%
  slice_head(n = 20) %>%
  mutate(group = "African-Prioritized", group_rank = row_number())

# Get top 20 Pan-cancer-prioritized genes (ranked by priority score)
top20_pancancer_pri <- gene_scores %>%
  filter(!is.na(priority_score)) %>%
  mutate(
    score_difference = african_score - combined_score,
    in_panel_70 = gene %in% panel_70$gene,
    priority_type = case_when(
      score_difference > 0.05 ~ "African-prioritized",
      score_difference < -0.05 ~ "Pan-cancer-prioritized",
      TRUE ~ "Balanced"
    )
  ) %>%
  filter(priority_type == "Pan-cancer-prioritized") %>%
  arrange(desc(priority_score)) %>%
  slice_head(n = 20) %>%
  mutate(group = "Pan-Cancer-Prioritized", group_rank = row_number())

# Combine both groups
top40_combined <- bind_rows(top20_african_pri, top20_pancancer_pri) %>%
  mutate(
    group = factor(group, levels = c("African-Prioritized", "Pan-Cancer-Prioritized"))
  )

# Create scatter plot
p_top20_scatter <- ggplot(top40_combined,
                          aes(x = combined_score,
                              y = african_score,
                              color = group)) +
  # Add diagonal reference line (y = x)
  geom_abline(
    intercept = 0,
    slope = 1,
    linetype = "dashed",
    color = "gray30",
    size = 1,
    alpha = 0.7
  ) +
  # Plot points
  geom_point(
    aes(size = literature_count,
        shape = in_panel_70),
    alpha = 0.8
  ) +
  # Add gene labels
  geom_text_repel(
    aes(label = gene),
    size = 3.5,
    fontface = "bold",
    max.overlaps = Inf,
    box.padding = 0.5,
    point.padding = 0.3,
    segment.color = "grey50",
    segment.size = 0.3,
    segment.alpha = 0.6,
    force = 2,
    force_pull = 0.5,
    max.time = 3,
    max.iter = 100000,
    seed = 123
  ) +
  # Color scale
  scale_color_manual(
    values = c(
      "African-Prioritized" = "#E64B35",
      "Pan-Cancer-Prioritized" = "#4DBBD5"
    ),
    name = "Prioritization Group"
  ) +
  # Size scale for literature support
  scale_size_continuous(
    range = c(3, 12),
    name = "Literature\nCount",
    breaks = pretty_breaks(n = 4)
  ) +
  # Shape for panel inclusion
  scale_shape_manual(
    values = c(16, 18),  # Circle and diamond
    name = "In 70-gene\npanel",
    labels = c("No", "Yes")
  ) +
  # Set axis limits to 0-1
  scale_x_continuous(limits = c(0, 1), breaks = seq(0, 1, 0.2)) +
  scale_y_continuous(limits = c(0, 1), breaks = seq(0, 1, 0.2)) +
  # Labels
  labs(
    title = "Top 20 African-Prioritized vs Top 20 Pan-Cancer-Prioritized Genes",
    subtitle = "Each group shows top 20 genes by priority score\nGenes above diagonal (y=x) have higher African score; below have higher pan-cancer score",
    x = "Pan-Cancer Score (mutation frequency × breadth)",
    y = "African-Weighted Score (GLOBOCAN 2022 weighted)",
    caption = "Point size = literature support | Shape = panel inclusion\nDashed line represents equal scores (y = x)"
  ) +
  # Theme
  theme_minimal(base_size = 12) +
  theme(
    legend.position = "right",
    legend.box = "vertical",
    legend.key.size = unit(0.7, "cm"),
    plot.title = element_text(size = 15, face = "bold", hjust = 0),
    plot.subtitle = element_text(size = 10, color = "gray40", hjust = 0, lineheight = 1.2),
    plot.caption = element_text(size = 9, color = "gray50", hjust = 0, lineheight = 1.2),
    panel.grid.minor = element_blank(),
    panel.grid.major = element_line(color = "gray90", size = 0.3),
    panel.background = element_rect(fill = "white", color = NA),
    plot.background = element_rect(fill = "white", color = NA),
    plot.margin = margin(20, 20, 20, 20)
  ) +
  # Equal aspect ratio with fixed 0-1 limits
  coord_fixed(ratio = 1, xlim = c(0, 1), ylim = c(0, 1))

# Save scatter plot
ggsave(
  here("results", "figures", "fig1_top20_gene_prioritization_both.pdf"),
  p_top20_scatter,
  width = 14,
  height = 12,
  dpi = 300,
  bg = "white"
)

# Print summary statistics
cat("\n=== Top 20 African-Prioritized Genes (by priority score) ===\n")
top20_african_pri %>%
  select(group_rank, gene, priority_score, score_difference, african_score, combined_score, in_panel_70, literature_count) %>%
  print()

cat("\n\n=== Top 20 Pan-Cancer-Prioritized Genes (by priority score) ===\n")
top20_pancancer_pri %>%
  select(group_rank, gene, priority_score, score_difference, african_score, combined_score, in_panel_70, literature_count) %>%
  print()

cat("\n=== Score Statistics ===\n")
cat("\nAfrican-Prioritized Genes:\n")
top20_african_pri %>%
  summarise(
    n = n(),
    mean_priority = mean(priority_score),
    mean_african = mean(african_score),
    mean_pancancer = mean(combined_score),
    mean_difference = mean(score_difference),
    in_panel = sum(in_panel_70)
  ) %>%
  print()

cat("\nPan-Cancer-Prioritized Genes:\n")
top20_pancancer_pri %>%
  summarise(
    n = n(),
    mean_priority = mean(priority_score),
    mean_african = mean(african_score),
    mean_pancancer = mean(combined_score),
    mean_difference = mean(score_difference),
    in_panel = sum(in_panel_70)
  ) %>%
  print()

# Save tables
write_csv(top40_combined, here("results", "tables", "top20_each_prioritization_group.csv"))


# Figure 2: Panel composition ----


# A) Gene category distribution
p2a <- panel_70 %>%
  count(gene_category) %>%
  mutate(percentage = n / sum(n) * 100) %>%
  ggplot(aes(x = reorder(gene_category, n), y = n, fill = gene_category)) +
  geom_col() +
  geom_text(aes(label = paste0(n, "\n(", round(percentage, 1), "%)")),
            hjust = -0.1, size = 3) +
  coord_flip() +
  scale_fill_nejm() +
  labs(
    title = "A) Gene Categories in 70-gene Panel",
    x = "",
    y = "Number of Genes"
  ) +
  theme(legend.position = "none") +
  ylim(0, max(table(panel_70$gene_category)) * 1.2)

# B) Tier distribution
p2b <- panel_70 %>%
  count(tier) %>%
  mutate(percentage = n / sum(n) * 100) %>%
  ggplot(aes(x = reorder(tier, n), y = n, fill = tier)) +
  geom_col() +
  geom_text(aes(label = paste0(n, "\n(", round(percentage, 1), "%)")),
            hjust = -0.1, size = 3) +
  coord_flip() +
  scale_fill_manual(values = tier_colors) +
  labs(
    title = "B) Priority Tiers in 70-gene Panel",
    x = "",
    y = "Number of Genes"
  ) +
  theme(legend.position = "none") +
  ylim(0, max(table(panel_70$tier)) * 1.2)

# C) Actionability
p2c <- panel_70 %>%
  mutate(actionable = ifelse(is_actionable, "Actionable", "Non-Actionable")) %>%
  count(actionable) %>%
  ggplot(aes(x = "", y = n, fill = actionable)) +
  geom_col(width = 1) +
  geom_text(aes(label = paste0(n, "\n(", round(n/sum(n)*100, 1), "%)")),
            position = position_stack(vjust = 0.5)) +
  coord_polar("y") +
  scale_fill_manual(values = c("#E64B35", "#4DBBD5")) +
  labs(
    title = "C) Actionability",
    fill = ""
  ) +
  theme_void() +
  theme(plot.title = element_text(hjust = 0.5))

# Combine panels
p2 <- (p2a | p2b | p2c) +
  plot_annotation(
    title = "Figure 2: 70-gene Panel Composition",
    theme = theme(plot.title = element_text(size = 16, face = "bold"))
  )

ggsave(
  here("results", "figures", "fig2_panel_composition.pdf"),
  p2, width = 15, height = 10, dpi = 300
)


# Figure 3: Mutation frequency heatmap ----


# Prepare matrix for heatmap
# Function to prepare heatmap matrix using base R

prepare_heatmap_matrix_base <- function(freq_df, panel_genes) {
  
  cat("Preparing heatmap matrix (base R)...\n")
  
  # Filter
  freq_filtered <- freq_df[freq_df$gene %in% panel_genes$gene, ]
  
  cat("  Genes:", length(unique(freq_filtered$gene)), "\n")
  cat("  Cancer types:", length(unique(freq_filtered$cancer_type)), "\n")
  
  # Get unique genes and cancer types
  genes <- unique(freq_filtered$gene)
  cancer_types <- unique(freq_filtered$cancer_type)
  
  # Create empty matrix
  freq_matrix <- matrix(0, nrow = length(genes), ncol = length(cancer_types))
  rownames(freq_matrix) <- genes
  colnames(freq_matrix) <- cancer_types
  
  # Fill matrix
  for (i in 1:nrow(freq_filtered)) {
    gene <- freq_filtered$gene[i]
    cancer <- freq_filtered$cancer_type[i]
    freq <- freq_filtered$mutation_frequency[i]
    
    if (!is.na(freq)) {
      freq_matrix[gene, cancer] <- freq
    }
  }
  
  cat("  Matrix:", nrow(freq_matrix), "x", ncol(freq_matrix), "\n\n")
  
  return(freq_matrix)
}

# Use base R version (most reliable)
freq_matrix <- prepare_heatmap_matrix_base(tcga_freq, panel_70)

# Order by total frequency
gene_order <- rowSums(freq_matrix) %>% sort(decreasing = TRUE) %>% names()
freq_matrix <- freq_matrix[gene_order, ]

# Create annotation
gene_annotation <- panel_70 %>%
  filter(gene %in% rownames(freq_matrix)) %>%
  select(gene, gene_category, is_actionable, tier) %>%
  column_to_rownames("gene")

# Create heatmap

# Palette - equential Reds (darker version, avoiding white)
col_fun <- colorRamp2(
  c(0, 0.1, 0.3, 0.5, 1),
  c("#fef0d9", "#fdcc8a", "#fc8d59", "#e34a33", "#b30000")
)

# Save as PDF
pdf(
  here("results", "figures", "fig3_mutation_heatmap.pdf"),
  width = 10, 
  height = 15
)

ht <- Heatmap(
  freq_matrix,
  name = "Mutation\nFrequency",
  col = col_fun,
  
  # Row parameters
  cluster_rows = FALSE,
  show_row_names = TRUE,
  row_names_side = "left",
  row_names_gp = gpar(fontsize = 8),
  
  # Column parameters
  cluster_columns = FALSE,
  show_column_names = TRUE,
  column_names_side = "bottom",
  column_names_rot = 45,
  
  # Annotations
  left_annotation = rowAnnotation(
    Category = gene_annotation$gene_category,
    Actionable = gene_annotation$is_actionable,
    Tier = gene_annotation$tier,
    col = list(
      Category = setNames(cancer_colors[1:length(unique(gene_annotation$gene_category))],
                          unique(gene_annotation$gene_category)),
      Actionable = c("TRUE" = "#E64B35", "FALSE" = "grey80")
    ),
    annotation_legend_param = list(
      Category = list(title = "Gene Category"),
      Actionable = list(title = "Actionable")
    )
  ),
  
  # Cell size
  width = ncol(freq_matrix) * unit(8, "mm"),
  height = nrow(freq_matrix) * unit(4, "mm"),
  
  # Title
  column_title = "Mutation Frequencies Across Cancer Types (70-gene Panel)",
  column_title_gp = gpar(fontsize = 14, fontface = "bold")
)

ComplexHeatmap::draw(ht)

dev.off()

cat("✓ Heatmap saved as PDF!\n")


# Figure 4: Panel coverage analysis ----
# Re-annotate cancer types:
coverage_stats <- coverage_stats %>%
  mutate(
    cancer_label = case_when(
      grepl("brca", project)      ~ "Breast (BRCA)",
      grepl("prad", project)      ~ "Prostate (PRAD)",
      grepl("cesc", project)      ~ "Cervical (CESC)",
      grepl("lihc", project)      ~ "Liver (LIHC)",
      grepl("coadread", project)  ~ "Colorectal (COAD/READ)",
      grepl("luad", project)      ~ "Lung Adeno (LUAD)",
      grepl("lusc", project)      ~ "Lung Squam (LUSC)",
      grepl("ov_", project)       ~ "Ovarian (OV)",
      grepl("blca", project)      ~ "Bladder (BLCA)",
      grepl("stad", project)      ~ "Stomach (STAD)",
      grepl("esca", project)      ~ "Esophageal (ESCA)",
      grepl("ucec", project)      ~ "Endometrial (UCEC)",
      grepl("paad", project)      ~ "Pancreatic (PAAD)",
      TRUE ~ project
    )
  )

# A) Mutation coverage by cancer type
p4a <- coverage_stats %>%
  ggplot(aes(x = reorder(cancer_label, coverage_percent), y = coverage_percent, 
             fill = base::as.factor(panel_size))) +
  geom_col(position = "dodge") +
  coord_flip() +
  scale_fill_nejm(name = "Panel Size") +
  labs(
    title = "A) Mutation Coverage by Cancer Type",
    x = "Cancer Type",
    y = "Mutation Coverage (%)"
  ) +
  theme(legend.position = "bottom")

# B) Sample coverage by cancer type
p4b <- coverage_stats %>%
  ggplot(aes(x = reorder(cancer_label, sample_coverage_percent), 
             y = sample_coverage_percent,
             fill = base::as.factor(panel_size))) +
  geom_col(position = "dodge") +
  coord_flip() +
  scale_fill_nejm(name = "Panel Size") +
  labs(
    title = "B) Sample Coverage by Cancer Type",
    x = "Cancer Type",
    y = "Samples with Panel Mutation (%)"
  ) +
  theme(legend.position = "bottom")

# C) Coverage vs panel size
p4c <- coverage_stats %>%
  group_by(panel_size) %>%
  summarise(
    mean_mutation_cov = mean(coverage_percent),
    se_mutation_cov = sd(coverage_percent) / sqrt(n()),
    mean_sample_cov = mean(sample_coverage_percent),
    se_sample_cov = sd(sample_coverage_percent) / sqrt(n())
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
    title = "C) Average Coverage vs Panel Size",
    x = "Panel Size (number of genes)",
    y = "Coverage (%)"
  ) +
  theme(legend.position = "bottom")

# Combine
p4 <- (p4a / p4b) | p4c
p4 <- p4 +
  plot_annotation(
    title = "Figure 4: Panel Coverage Analysis",
    theme = theme(plot.title = element_text(size = 16, face = "bold"))
  )

print(p4)

ggsave(
  here("results", "figures", "fig4_coverage_analysis.pdf"),
  p4, width = 15, height = 10, dpi = 300
)


# Figure 5: Cost-benefit analysis ----


# First, let's check the data range to set appropriate axis breaks
cost_range <- range(cost_benefit$cost_per_actionable, na.rm = TRUE)
benefit_range <- range(cost_benefit$benefit_score, na.rm = TRUE)
coverage_range <- range(cost_benefit$mean_sample_coverage, na.rm = TRUE)

# Calculate min and max size for the scatter plot
min_size <- min(cost_benefit$n_actionable, na.rm = TRUE)
max_size <- max(cost_benefit$n_actionable, na.rm = TRUE)


# A) Cost per sample - Y-axis starts at 0
p5a <- cost_benefit %>%
  ggplot(aes(x = panel_size, y = cost_per_actionable)) +
  geom_line(color = "#E64B35", size = 1.5) +
  geom_point(color = "#E64B35", size = 4) +
  geom_text(aes(label = paste0("$", round(cost_per_actionable, 1))),
            vjust = -1, size = 3.5, check_overlap = TRUE) +
  scale_x_continuous(
    breaks = seq(min(cost_benefit$panel_size), 
                 max(cost_benefit$panel_size), 
                 length.out = 5) %>% round(0)
  ) +
  scale_y_continuous(
    limits = c(0, max(cost_range) * 1.1),  # Starts at 0
    breaks = scales::pretty_breaks(n = 6),
    expand = expansion(mult = c(0, 0.05))  # Remove padding at bottom
  ) +
  labs(
    title = "A) Cost per actionable target",
    x = "Panel Size (genes)",
    y = "Cost per actionable target (USD)"
  ) +
  theme_minimal(base_size = 12) +
  theme(
    panel.grid.major = element_line(color = "grey90"),
    panel.grid.minor = element_blank()
  )

# B) Coverage per dollar - Y-axis starts at 0
p5b <- cost_benefit %>%
  ggplot(aes(x = panel_size, y = coverage_per_dollar)) +
  geom_line(color = "#4DBBD5", size = 1.5) +
  geom_point(color = "#4DBBD5", size = 4) +
  geom_text(aes(label = round(cost_benefit$coverage_per_dollar, 2)),
            vjust = -1, size = 3.5, check_overlap = TRUE) +
  scale_x_continuous(
    breaks = seq(min(cost_benefit$panel_size), 
                 max(cost_benefit$panel_size), 
                 length.out = 5) %>% round(0)
  ) +
  scale_y_continuous(
    limits = c(0, max(cost_benefit$coverage_per_dollar) * 1.1),  # Starts at 0
    breaks = scales::pretty_breaks(n = 6),
    expand = expansion(mult = c(0, 0.05))  # Remove padding at bottom
  ) +
  labs(
    title = "B) Coverage efficiency",
    x = "Panel Size (genes)",
    y = "Sample coverage per dollar (%/USD)"
  ) +
  theme_minimal(base_size = 12) +
  theme(
    panel.grid.major = element_line(color = "grey90"),
    panel.grid.minor = element_blank()
  )

# C) Benefit score - Y-axis starts at 0
p5c <- cost_benefit %>%
  ggplot(aes(x = panel_size, y = benefit_score)) +
  geom_line(color = "#00A087", size = 1.5) +
  geom_point(color = "#00A087", size = 4) +
  geom_text(aes(label = round(benefit_score, 3)),
            vjust = -1, size = 3.5, check_overlap = TRUE) +
  scale_x_continuous(
    breaks = seq(min(cost_benefit$panel_size), 
                 max(cost_benefit$panel_size), 
                 length.out = 5) %>% round(0)
  ) +
  scale_y_continuous(
    limits = c(0, max(benefit_range) * 1.1),  # Starts at 0
    breaks = scales::pretty_breaks(n = 6),
    expand = expansion(mult = c(0, 0.05))  # Remove padding at bottom
  ) +
  labs(
    title = "C) Overall benefit score",
    x = "Panel Size (genes)",
    y = "Benefit score"
  ) +
  theme_minimal(base_size = 12) +
  theme(
    panel.grid.major = element_line(color = "grey90"),
    panel.grid.minor = element_blank()
  )

# D) Cost vs coverage scatter - Y-axis starts at 0, X-axis starts at 0
p5d <- cost_benefit %>%
  ggplot(aes(x = cost_per_actionable, y = mean_sample_coverage)) +
  geom_point(aes(size = n_actionable, color = base::as.factor(panel_size)),
             alpha = 0.8) +
  geom_text_repel(
    aes(label = paste0(panel_size, " genes")), 
    size = 3.5,
    box.padding = 0.5,
    point.padding = 0.3,
    max.overlaps = Inf
  ) +
  scale_size_continuous(
    range = c(3, 12),
    breaks = pretty(c(min_size, max_size), n = 4),
    name = "Actionable\nGenes"
  ) +
  scale_y_continuous(
    limits = c(0, 100),  # Starts at 0
    breaks = seq(0, 100, 20),
    labels = function(x) paste0(x, "%"),
    expand = expansion(mult = c(0, 0.05))  # Remove padding at bottom
  ) +
  scale_x_continuous(
    limits = c(0, max(cost_range) * 1.1),  # Starts at 0
    labels = scales::dollar_format(),
    expand = expansion(mult = c(0, 0.05))  # Remove padding at left
  ) +
  scale_color_brewer(
    palette = "Set1",
    name = "Panel Size"
  ) +
  labs(
    title = "D) Cost vs coverage trade-off",
    x = "Cost per actionable target (USD)",
    y = "Mean sample coverage (%)"
  ) +
  theme_minimal(base_size = 12) +
  theme(
    panel.grid.major = element_line(color = "grey90"),
    panel.grid.minor = element_blank(),
    legend.position = "right",
    legend.box = "vertical"
  )

# Combine plots with better spacing
p5 <- (p5a | p5b | p5c) / p5d +
  plot_layout(heights = c(1, 1.5))

# Add annotation
p5 <- p5 +
  plot_annotation(
    title = "Figure 5: Cost-benefit analysis of gene panel sizing",
    subtitle = "Analysis of optimal panel size based on cost efficiency and clinical benefit",
    caption = "Data source: Panel optimization analysis | Note: All costs are per sample",
    theme = theme(
      plot.title = element_text(size = 16, face = "bold", hjust = 0.5),
      plot.subtitle = element_text(size = 12, hjust = 0.5, color = "grey40"),
      plot.caption = element_text(size = 10, color = "grey60", hjust = 1)
    )
  )

# Save with appropriate dimensions
ggsave(
  here("results", "figures", "fig5_cost_benefit.pdf"),
  p5, 
  width = 20, 
  height = 15, 
  dpi = 300
)



# Figure 6: Panel Comparison ----


cat("Creating Figure 6: Panel comparison...\n")

# Define existing panels properly
existing_panels <- list(
  FoundationOne = c(
    "ABL1", "ACVR1B", "AKT1", "AKT2", "AKT3", "ALK", "ALOX12B", "AMER1", "APC",
    "AR", "ARAF", "ARFRP1", "ARID1A", "ASXL1", "ATM", "ATR", "ATRX", "AURKA",
    "AURKB", "AXIN1", "AXL", "BAP1", "BARD1", "BCL2", "BCL2L1", "BCL2L2", "BCL6",
    "BCOR", "BCORL1", "BRAF", "BRCA1", "BRCA2", "BRD4", "BRIP1", "BTG1", "BTG2",
    "BTK", "C11ORF30", "CALR", "CARD11", "CASP8", "CBFB", "CBL", "CCND1", "CCND2",
    "CCND3", "CCNE1", "CD22", "CD274", "CD70", "CD79A", "CD79B", "CDC73", "CDH1",
    "CDK12", "CDK4", "CDK6", "CDK8", "CDKN1A", "CDKN1B", "CDKN2A", "CDKN2B", "CDKN2C",
    "CEBPA", "CHEK1", "CHEK2", "CIC", "CREBBP", "CRKL", "CSF1R", "CSF3R", "CTCF",
    "CTNNA1", "CTNNB1", "CUL3", "CUL4A", "CXCR4", "CYP17A1", "DAXX", "DDR1", "DDR2",
    "DIS3", "DNMT3A", "DOT1L", "EED", "EGFR", "EP300", "EPHA3", "EPHB1", "EPHB4",
    "ERBB2", "ERBB3", "ERBB4", "ERCC4", "ERG", "ERRFI1", "ESR1", "EZH2", "FAM46C",
    "FANCA", "FANCC", "FANCG", "FANCL", "FAS", "FBXW7", "FGF10", "FGF12", "FGF14",
    "FGF19", "FGF23", "FGF3", "FGF4", "FGF6", "FGFR1", "FGFR2", "FGFR3", "FGFR4",
    "FH", "FLCN", "FLT1", "FLT3", "FOXL2", "FUBP1", "GATA3", "GATA4", "GATA6",
    "GNA11", "GNA13", "GNAQ", "GNAS", "HRAS", "IDH1", "IDH2", "JAK1", "JAK2", 
    "JAK3", "KDR", "KEAP1", "KIT", "KRAS", "MAP2K1", "MAP2K2", "MAP3K1", "MCL1",
    "MDM2", "MDM4", "MET", "MLH1", "MSH2", "MSH6", "MTOR", "MYC", "MYCL", "MYCN",
    "NBN", "NF1", "NF2", "NFE2L2", "NOTCH1", "NOTCH2", "NOTCH3", "NPM1", "NRAS",
    "NTRK1", "NTRK2", "NTRK3", "PALB2", "PDGFRA", "PDGFRB", "PIK3CA", "PIK3R1",
    "PMS2", "PTEN", "PTPN11", "RAD51", "RAD51B", "RAD51C", "RAD51D", "RAF1",
    "RB1", "RET", "ROS1", "SDHA", "SDHB", "SDHC", "SDHD", "SETD2", "SF3B1",
    "SMAD2", "SMAD4", "SMARCA4", "SMARCB1", "SMO", "STK11", "TP53", "TSC1", 
    "TSC2", "VHL"
  ),
  MSK_IMPACT = c(
    "ABL1", "AKT1", "AKT2", "AKT3", "ALK", "ALOX12B", "AMER1", "APC", "AR",
    "ARAF", "ARID1A", "ARID1B", "ARID2", "ASXL1", "ASXL2", "ATM", "ATR", "ATRX",
    "AURKA", "AURKB", "AXIN1", "AXIN2", "AXL", "B2M", "BAP1", "BARD1", "BCL2",
    "BCL2L1", "BCL2L11", "BCL6", "BCOR", "BRAF", "BRCA1", "BRCA2", "BRD4",
    "BRIP1", "BTK", "CALR", "CARD11", "CASP8", "CBFB", "CBL", "CCND1", "CCND2",
    "CCND3", "CCNE1", "CD274", "CD79A", "CD79B", "CDC73", "CDH1", "CDK4", "CDK6",
    "CDK8", "CDKN1A", "CDKN1B", "CDKN2A", "CDKN2B", "CDKN2C", "CEBPA", "CHEK1",
    "CHEK2", "CIC", "CREBBP", "CRKL", "CSF1R", "CSF3R", "CTCF", "CTLA4", 
    "CTNNB1", "CUL3", "CXCR4", "DAXX", "DDR2", "DICER1", "DIS3", "DNMT3A",
    "DOT1L", "EED", "EGFR", "EP300", "EPHA3", "EPHA5", "EPHA7", "EPHB1",
    "ERBB2", "ERBB3", "ERBB4", "ERCC2", "ERCC3", "ERCC4", "ERCC5", "ERF",
    "ERG", "ERRFI1", "ESR1", "ETV1", "ETV6", "EZH1", "EZH2", "FANCA", "FANCC",
    "FBXW7", "FGF19", "FGF3", "FGF4", "FGFR1", "FGFR2", "FGFR3", "FGFR4",
    "FH", "FLCN", "FLT1", "FLT3", "FLT4", "FOXA1", "FOXL2", "FOXO1", "FOXP1",
    "FUBP1", "GATA1", "GATA2", "GATA3", "GNA11", "GNAQ", "GNAS", "GNB1",
    "H3F3A", "H3F3B", "H3F3C", "HGF", "HNF1A", "HRAS", "ID3", "IDH1", "IDH2",
    "IGF1R", "IKBKE", "IKZF1", "IL7R", "INPP4A", "INPP4B", "IRF4", "IRS1",
    "IRS2", "JAK1", "JAK2", "JAK3", "JUN", "KDM5A", "KDM5C", "KDM6A", "KDR",
    "KEAP1", "KIT", "KMT2A", "KMT2C", "KMT2D", "KRAS", "LATS1", "LATS2",
    "MAP2K1", "MAP2K2", "MAP2K4", "MAP3K1", "MAP3K13", "MAPK1", "MAPK3",
    "MAX", "MCL1", "MDM2", "MDM4", "MED12", "MEF2B", "MEN1", "MET", "MGA",
    "MITF", "MLH1", "MPL", "MSH2", "MSH3", "MSH6", "MST1R", "MTAP", "MTOR",
    "MUTYH", "MYC", "MYCL", "MYCN", "MYD88", "NBN", "NF1", "NF2", "NFE2L2",
    "NFKBIA", "NKX2-1", "NOTCH1", "NOTCH2", "NOTCH3", "NOTCH4", "NPM1", "NRAS",
    "NTRK1", "NTRK2", "NTRK3", "PALB2", "PARP1", "PAX5", "PBRM1", "PDCD1",
    "PDCD1LG2", "PDGFRA", "PDGFRB", "PIK3CA", "PIK3CB", "PIK3R1", "PIM1",
    "PMS2", "POLD1", "POLE", "PPARG", "PPM1D", "PPP2R1A", "PRDM1", "PRKCI",
    "PTEN", "PTPN11", "RAC1", "RAD21", "RAD50", "RAD51", "RAD51B", "RAD51C",
    "RAD51D", "RAF1", "RARA", "RB1", "RBM10", "REL", "RET", "RICTOR", "RNF43",
    "ROS1", "RPTOR", "RUNX1", "SDHA", "SDHB", "SDHC", "SDHD", "SETD2", "SF3B1",
    "SMAD2", "SMAD3", "SMAD4", "SMARCA4", "SMARCB1", "SMO", "SOCS1", "SOX2",
    "SOX9", "SPEN", "SPOP", "SRC", "SRSF2", "STAG2", "STAT3", "STK11", "SUFU",
    "TET2", "TGFBR2", "TP53", "TSC1", "TSC2", "U2AF1", "VEGFA", "VHL", "WT1",
    "XPO1", "ZRSR2"
  )
)


cat("Creating Figure 6: Panel comparison for all panel sizes...\n\n")

# Ensure existing_panels has proper names
if (is.null(names(existing_panels))) {
  names(existing_panels) <- c("FoundationOne", "MSK_IMPACT")
}


# Create overlap data for each panel size

# 70-gene panel overlap
overlap_data_70 <- tibble(
  gene = base::union(base::union(panel_70$gene, existing_panels$FoundationOne), 
               existing_panels$MSK_IMPACT)
) %>%
  mutate(
    African_panel_70 = as.integer(gene %in% panel_70$gene),
    FoundationOne = as.integer(gene %in% existing_panels$FoundationOne),
    MSK_IMPACT = as.integer(gene %in% existing_panels$MSK_IMPACT)
  ) %>%
  as.data.frame()

cat("70-gene panel overlap:\n")
cat("  Total unique genes:", nrow(overlap_data_70), "\n")
cat("  In African Panel (50):", sum(overlap_data_70$African_panel_70), "\n")
cat("  In FoundationOne:", sum(overlap_data_70$FoundationOne), "\n")
cat("  In MSK_IMPACT:", sum(overlap_data_70$MSK_IMPACT), "\n")
cat("  In all three:", sum(overlap_data_70$African_panel_70 == 1 & 
                             overlap_data_70$FoundationOne == 1 & 
                             overlap_data_70$MSK_IMPACT == 1), "\n\n")

# 130-gene panel overlap
overlap_data_130 <- tibble(
  gene = base::union(base::union(panel_130$gene, existing_panels$FoundationOne), 
               existing_panels$MSK_IMPACT)
) %>%
  mutate(
    African_panel_130 = as.integer(gene %in% panel_130$gene),
    FoundationOne = as.integer(gene %in% existing_panels$FoundationOne),
    MSK_IMPACT = as.integer(gene %in% existing_panels$MSK_IMPACT)
  ) %>%
  as.data.frame()

cat("130-gene panel overlap:\n")
cat("  Total unique genes:", nrow(overlap_data_130), "\n")
cat("  In African Panel (100):", sum(overlap_data_130$African_panel_130), "\n")
cat("  In FoundationOne:", sum(overlap_data_130$FoundationOne), "\n")
cat("  In MSK_IMPACT:", sum(overlap_data_130$MSK_IMPACT), "\n")
cat("  In all three:", sum(overlap_data_130$African_panel_130 == 1 & 
                             overlap_data_130$FoundationOne == 1 & 
                             overlap_data_130$MSK_IMPACT == 1), "\n\n")

# 30-gene panel overlap
overlap_data_30 <- tibble(
  gene = base::union(base::union(panel_30$gene, existing_panels$FoundationOne), 
               existing_panels$MSK_IMPACT)
) %>%
  mutate(
    African_panel_30 = as.integer(gene %in% panel_30$gene),
    FoundationOne = as.integer(gene %in% existing_panels$FoundationOne),
    MSK_IMPACT = as.integer(gene %in% existing_panels$MSK_IMPACT)
  ) %>%
  as.data.frame()

cat("30-gene panel overlap:\n")
cat("  Total unique genes:", nrow(overlap_data_30), "\n")
cat("  In African Panel (30):", sum(overlap_data_30$African_panel_30), "\n")
cat("  In FoundationOne:", sum(overlap_data_30$FoundationOne), "\n")
cat("  In MSK_IMPACT:", sum(overlap_data_30$MSK_IMPACT), "\n")
cat("  In all three:", sum(overlap_data_30$African_panel_30 == 1 & 
                             overlap_data_30$FoundationOne == 1 & 
                             overlap_data_30$MSK_IMPACT == 1), "\n\n")

# Figure 6A: 70-gene panel comparison
cat("Creating Figure 6A (70-gene panel)...\n")

pdf(
  here("results", "figures", "fig6a_panel_comparison_70genes.pdf"),
  width = 10, height = 6
)

upset(
  overlap_data_70,
  sets = c("African_panel_70", "FoundationOne", "MSK_IMPACT"),
  sets.x.label = "Genes per Panel",
  mainbar.y.label = "Gene Intersections",
  keep.order = TRUE,
  order.by = "freq",
  main.bar.color = "#E64B35",
  sets.bar.color = "#4DBBD5",
  text.scale = c(1.5, 1.3, 1.3, 1.2, 1.5, 1.3),
  point.size = 3.5,
  line.size = 1,
  mb.ratio = c(0.6, 0.4)
)

dev.off()

cat("✓ Figure 6A saved\n")

# Figure 6B: 100-gene panel comparison
cat("Creating Figure 6B (130-gene panel)...\n")

pdf(
  here("results", "figures", "fig6b_panel_comparison_130genes.pdf"),
  width = 10, height = 6
)

upset(
  overlap_data_130,
  sets = c("African_panel_130", "FoundationOne", "MSK_IMPACT"),
  sets.x.label = "Genes per Panel",
  mainbar.y.label = "Gene Intersections",
  keep.order = TRUE,
  order.by = "freq",
  main.bar.color = "#E64B35",
  sets.bar.color = "#4DBBD5",
  text.scale = c(1.5, 1.3, 1.3, 1.2, 1.5, 1.3),
  point.size = 3.5,
  line.size = 1,
  mb.ratio = c(0.6, 0.4)
)

dev.off()

cat("✓ Figure 6B saved\n")

# Figure 6C: 30-gene panel comparison
cat("Creating Figure 6C (30-gene panel)...\n")

pdf(
  here("results", "figures", "fig6c_panel_comparison_30genes.pdf"),
  width = 10, height = 6
)

upset(
  overlap_data_30,
  sets = c("African_panel_30", "FoundationOne", "MSK_IMPACT"),
  sets.x.label = "Genes per Panel",
  mainbar.y.label = "Gene Intersections",
  keep.order = TRUE,
  order.by = "freq",
  main.bar.color = "#E64B35",
  sets.bar.color = "#4DBBD5",
  text.scale = c(1.5, 1.3, 1.3, 1.2, 1.5, 1.3),
  point.size = 3.5,
  line.size = 1,
  mb.ratio = c(0.6, 0.4)
)

dev.off()

cat("✓ Figure 6C saved\n\n")


# Create combined comparison figure

cat("Creating combined comparison figure...\n")

# Calculate overlap statistics
overlap_stats <- tibble(
  Panel_Size = c(70, 130, 30),
  African_Genes = c(
    sum(overlap_data_70$African_panel_70),
    sum(overlap_data_130$African_panel_130),
    sum(overlap_data_30$African_panel_30)
  ),
  In_FoundationOne = c(
    sum(overlap_data_70$African_panel_70 == 1 & overlap_data_70$FoundationOne == 1),
    sum(overlap_data_130$African_panel_130 == 1 & overlap_data_130$FoundationOne == 1),
    sum(overlap_data_30$African_panel_30 == 1 & overlap_data_30$FoundationOne == 1)
  ),
  In_MSK_IMPACT = c(
    sum(overlap_data_70$African_panel_70 == 1 & overlap_data_70$MSK_IMPACT == 1),
    sum(overlap_data_130$African_panel_130 == 1 & overlap_data_130$MSK_IMPACT == 1),
    sum(overlap_data_30$African_panel_30 == 1 & overlap_data_30$MSK_IMPACT == 1)
  ),
  In_Both_Panels = c(
    sum(overlap_data_70$African_panel_70 == 1 & 
          overlap_data_70$FoundationOne == 1 & 
          overlap_data_70$MSK_IMPACT == 1),
    sum(overlap_data_130$African_panel_130 == 1 & 
          overlap_data_130$FoundationOne == 1 & 
          overlap_data_130$MSK_IMPACT == 1),
    sum(overlap_data_30$African_panel_30 == 1 & 
          overlap_data_30$FoundationOne == 1 & 
          overlap_data_30$MSK_IMPACT == 1)
  ),
  Unique_to_African = c(
    sum(overlap_data_70$African_panel_70 == 1 & 
          overlap_data_70$FoundationOne == 0 & 
          overlap_data_70$MSK_IMPACT == 0),
    sum(overlap_data_130$African_panel_130 == 1 & 
          overlap_data_130$FoundationOne == 0 & 
          overlap_data_130$MSK_IMPACT == 0),
    sum(overlap_data_30$African_panel_30 == 1 & 
          overlap_data_30$FoundationOne == 0 & 
          overlap_data_30$MSK_IMPACT == 0)
  )
) %>%
  mutate(
    Overlap_Percent = (In_Both_Panels / African_Genes) * 100,
    Unique_Percent = (Unique_to_African / African_Genes) * 100
  )

# Save statistics
write_csv(overlap_stats, here("results", "tables", "panel_overlap_statistics.csv"))

# Create summary plot
p_overlap_summary <- overlap_stats %>%
  tidyr::pivot_longer(
    cols = c(In_FoundationOne, In_MSK_IMPACT, In_Both_Panels, Unique_to_African),
    names_to = "Category",
    values_to = "Count"
  ) %>%
  mutate(
    Category = factor(Category, levels = c(
      "Unique_to_African", "In_FoundationOne", "In_MSK_IMPACT", "In_Both_Panels"
    )),
    Category_Label = case_when(
      Category == "Unique_to_African" ~ "Unique to African Panel",
      Category == "In_FoundationOne" ~ "Overlap with FoundationOne",
      Category == "In_MSK_IMPACT" ~ "Overlap with MSK-IMPACT",
      Category == "In_Both_Panels" ~ "In All Three Panels"
    )
  ) %>%
  ggplot(aes(x = base::as.factor(Panel_Size), y = Count, fill = Category_Label)) +
  geom_col(position = "dodge") +
  geom_text(aes(label = Count), position = position_dodge(width = 0.9), 
            vjust = -0.5, size = 3.5) +
  scale_fill_manual(
    values = c("#F39B7F", "#4DBBD5", "#00A087", "#E64B35"),
    name = ""
  ) +
  labs(
    title = "Panel Overlap with Commercial Panels",
    x = "African Panel Size (genes)",
    y = "Number of Genes"
  ) +
  theme_minimal(base_size = 12) +
  theme(legend.position = "bottom")

ggsave(
  here("results", "figures", "fig6d_overlap_summary.pdf"),
  p_overlap_summary, width = 10, height = 6, device = "pdf"
)

cat("✓ Figure 6D (summary) saved\n")


# Create comparison table

cat("\nCreating panel comparison table...\n")

comparison_table <- overlap_stats %>%
  dplyr::select(
    `Panel Size` = Panel_Size,
    `Total Genes` = African_Genes,
    `Overlap with FM` = In_FoundationOne,
    `Overlap with MSK` = In_MSK_IMPACT,
    `In All Three` = In_Both_Panels,
    `Unique to African` = Unique_to_African,
    `Unique %` = Unique_Percent
  ) %>%
  mutate(across(where(is.numeric), ~round(., 1)))

write_csv(comparison_table, here("results", "tables", "panel_comparison_table.csv"))

cat("\nPanel Comparison Summary:\n")
print(comparison_table)
cat("\n")


# Summary visualization comparing all three


cat("Creating summary comparison figure...\n")

# Overlap percentage plot
p_overlap_pct <- overlap_stats %>%
  tidyr::pivot_longer(
    cols = c(Overlap_Percent, Unique_Percent),
    names_to = "Type",
    values_to = "Percentage"
  ) %>%
  mutate(
    Type = ifelse(Type == "Overlap_Percent", 
                  "Overlap with Commercial Panels", 
                  "Unique to African Panel")
  ) %>%
  ggplot(aes(x = base::as.factor(Panel_Size), y = Percentage, fill = Type)) +
  geom_col(position = "stack", width = 0.7) +
  geom_text(aes(label = paste0(round(Percentage, 1), "%")), 
            position = position_stack(vjust = 0.5),
            color = "white", fontface = "bold", size = 4) +
  scale_fill_manual(values = c("#4DBBD5", "#E64B35"), name = "") +
  labs(
    title = "Overlap with Commercial Panels by Panel Size",
    x = "African Panel Size (genes)",
    y = "Percentage of Genes (%)"
  ) +
  theme_minimal(base_size = 12) +
  theme(legend.position = "bottom")

ggsave(
  here("results", "figures", "fig6e_overlap_percentages.pdf"),
  p_overlap_pct, width = 8, height = 6, device = "pdf"
)

# Absolute numbers plot
p_overlap_abs <- overlap_stats %>%
  dplyr::select(Panel_Size, In_FoundationOne, In_MSK_IMPACT, 
                In_Both_Panels, Unique_to_African) %>%
  tidyr::pivot_longer(
    cols = -Panel_Size,
    names_to = "Category",
    values_to = "Count"
  ) %>%
  mutate(
    Category = factor(Category, 
                      levels = c("Unique_to_African", "In_FoundationOne", 
                                 "In_MSK_IMPACT", "In_Both_Panels"),
                      labels = c("Unique to African", "In FoundationOne Only",
                                 "In MSK-IMPACT Only", "In All Three"))
  ) %>%
  ggplot(aes(x = base::as.factor(Panel_Size), y = Count, fill = Category)) +
  geom_col(position = "dodge", width = 0.8) +
  geom_text(aes(label = Count), 
            position = position_dodge(width = 0.8),
            vjust = -0.5, size = 3.5) +
  scale_fill_manual(
    values = c("#F39B7F", "#4DBBD5", "#00A087", "#E64B35"),
    name = ""
  ) +
  labs(
    title = "Gene Distribution Across Panels",
    x = "African Panel Size (genes)",
    y = "Number of Genes"
  ) +
  theme_minimal(base_size = 12) +
  theme(legend.position = "bottom")

ggsave(
  here("results", "figures", "fig6f_overlap_counts.pdf"),
  p_overlap_abs, width = 10, height = 6, device = "pdf"
)

cat("✓ Summary figures saved\n\n")



# Summary table


summary_table <- tibble(
  Panel = c("African 70-gene", "African 130-gene", "FoundationOne", "MSK-IMPACT"),
  `Number of Genes` = c(
    nrow(panel_70),
    nrow(panel_130),
    length(existing_panels$FoundationOne),
    length(existing_panels$MSK_IMPACT)
  ),
  `Estimated Cost` = c("$260-320", "$410-520", "$500-700", "$800-1000"),
  `African-Specific Genes` = c(
    sum(!panel_70$gene %in% c(existing_panels$FoundationOne, existing_panels$MSK_IMPACT)),
    sum(!panel_130$gene %in% c(existing_panels$FoundationOne, existing_panels$MSK_IMPACT)),
    0,
    0
  ),
  `Actionable Genes` = c(
    sum(panel_70$is_actionable, na.rm = TRUE),
    sum(panel_130$is_actionable, na.rm = TRUE),
    35,
    42
  )
)

write_csv(summary_table, here("results", "tables", "panel_comparison_summary.csv"))

cat("Figures saved in:", here("results", "figures"), "\n")
