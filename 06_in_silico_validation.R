# 06_in_silico_validation.R
# In silico performance analyses and power calculations

library(tidyverse)
library(pwr)
library(here)


# Function: Simulate panel performance ----


simulate_panel_performance <- function(
    panel_genes,
    mutation_rates,
    n_samples = 500,
    n_simulations = 10000
) {
  
  results <- map_dfr(1:n_simulations, function(sim) {
    
    # Simulate mutations
    mutations <- map_dfr(panel_genes$gene, function(current_gene) {
      
      gene_mut_rate <- mutation_rates %>%
        filter(gene == current_gene) %>%  # Use current_gene to avoid name conflict
        pull(mean_mutation_freq)
      
      if (length(gene_mut_rate) == 0) gene_mut_rate <- 0.05
      
      n_mutated <- rbinom(1, n_samples, gene_mut_rate)
      
      tibble(
        simulation = sim,
        gene = current_gene,
        n_mutated = n_mutated,
        mutation_rate = n_mutated / n_samples
      )
    })
    
    mutations
  })
  
  return(results)
}



# Run simulation
panel_70 <- read_csv(here("results", "tables", "panel_70_genes.csv"))
mutation_rates <- read_csv(here("data", "processed", "pancancer_mutation_scores.csv"))

simulation_results <- simulate_panel_performance(
  panel_70,
  mutation_rates,
  n_samples = 500,
  n_simulations = 10000
)

write_csv(simulation_results, here("results", "tables", "simulation_results.csv"))


# Function: Power analysis ----


calculate_power_analysis <- function(
    effect_sizes = c(0.01, 0.05, 0.10, 0.20),
    sample_sizes = seq(100, 1000, 100),
    alpha = 0.05,
    n_tests = 50
) {
  
  # Bonferroni correction
  adjusted_alpha <- alpha / n_tests
  
  power_results <- expand.grid(
    effect_size = effect_sizes,
    sample_size = sample_sizes
  ) %>%
    as_tibble() %>%
    mutate(
      power = map2_dbl(effect_size, sample_size, function(es, n) {
        pwr.p.test(
          h = ES.h(es, 0),
          n = n,
          sig.level = adjusted_alpha,
          alternative = "greater"
        )$power
      })
    )
  
  return(power_results)
}

power_results <- calculate_power_analysis()
write_csv(power_results, here("results", "tables", "power_analysis.csv"))

# Visualize
p_power <- ggplot(power_results, 
                  aes(x = sample_size, y = power, 
                      color = base::as.factor(effect_size))) +
  geom_line(size = 1) +
  geom_point(size = 2) +
  geom_hline(yintercept = 0.8, linetype = "dashed", color = "red") +
  scale_color_nejm(name = "Effect Size\n(MAF)") +
  labs(
    title = "Statistical Power Analysis for Variant Detection",
    x = "Sample Size",
    y = "Power",
    caption = "Red line indicates 80% power threshold"
  ) +
  theme_minimal(base_size = 12)

print (p_power)

ggsave(
  here("results", "figures", "fig7_power_analysis.pdf"),
  p_power, width = 10, height = 6, dpi = 300
)

