# 01_data_acquisition.R
# Mine PubMed for African cancer genomics literature
# Updated with GLOBOCAN 2022 data (version 1.1) - 08.02.2024

library(tidyverse)
library(rentrez)
library(europepmc)
library(here)

# Load configuration
config <- readRDS(here("data", "config.rds"))

# Set NCBI API key (increases rate limit)
# NCBI API key can be obtained from https://www.ncbi.nlm.nih.gov/account/settings/
# Will already be set up in config.rds via the setup script

if (!is.null(config$ncbi_api_key) && config$ncbi_api_key != "YOUR_NCBI_API_KEY") {
  set_entrez_key(config$ncbi_api_key)
}

# Helper function for separator lines
print_separator <- function(char = "=", width = 80) {
  cat(strrep(char, width), "\n")
}

# Function: Search PubMed ----

search_pubmed_african_cancer <- function(
    cancer_type = NULL,
    max_results = 1000,
    year_start = 2000
) {
  
  # Build search query
  base_query <- '("cancer"[Title/Abstract] OR "neoplasm"[Title/Abstract] OR "tumor"[Title/Abstract])'
  
  genomic_terms <- '("genomic"[Title/Abstract] OR "mutation"[Title/Abstract] OR 
                     "sequencing"[Title/Abstract] OR "exome"[Title/Abstract] OR 
                     "variant"[Title/Abstract] OR "SNP"[Title/Abstract])'
  
  african_terms <- '("Africa"[Title/Abstract] OR "African"[Title/Abstract] OR 
                     "Sub-Saharan"[Title/Abstract] OR "Nigeria"[Title/Abstract] OR 
                     "Kenya"[Title/Abstract] OR "South Africa"[Title/Abstract] OR
                     "Ghana"[Title/Abstract] OR "Ethiopia"[Title/Abstract] OR
                     "Uganda"[Title/Abstract] OR "Tanzania"[Title/Abstract] OR
                     "Black"[Title/Abstract] OR "African ancestry"[Title/Abstract])'
  
  query <- paste(base_query, "AND", genomic_terms, "AND", african_terms)
  
  if (!is.null(cancer_type)) {
    query <- paste(query, "AND", cancer_type, "[Title/Abstract]")
  }
  
  # Add date filter
  query <- paste0(query, ' AND ', year_start, ':3000[PDAT]')
  
  # Search PubMed
  cat("Searching PubMed with query:\n", query, "\n\n")
  
  search_results <- tryCatch({
    entrez_search(
      db = "pubmed",
      term = query,
      retmax = max_results,
      use_history = TRUE
    )
  }, error = function(e) {
    cat("Error searching PubMed:", e$message, "\n")
    return(NULL)
  })
  
  if (is.null(search_results) || search_results$count == 0) {
    cat("No results found\n")
    return(NULL)
  }
  
  cat("Found", search_results$count, "articles\n")
  
  # Fetch article details in batches
  batch_size <- 100
  num_batches <- ceiling(length(search_results$ids) / batch_size)
  
  all_summaries <- list()
  
  for (i in 1:num_batches) {
    cat("Fetching batch", i, "of", num_batches, "\n")
    
    start_idx <- (i - 1) * batch_size + 1
    end_idx <- min(i * batch_size, length(search_results$ids))
    batch_ids <- search_results$ids[start_idx:end_idx]
    
    summaries <- tryCatch({
      entrez_summary(db = "pubmed", id = batch_ids)
    }, error = function(e) {
      cat("Error fetching batch:", e$message, "\n")
      return(NULL)
    })
    
    if (!is.null(summaries)) {
      all_summaries <- c(all_summaries, summaries)
    }
    
    Sys.sleep(0.5)  # Be nice to NCBI servers !!!!
  }
  
  if (length(all_summaries) == 0) {
    return(NULL)
  }
  
  # Extract relevant information
  articles_df <- map_dfr(all_summaries, function(article) {
    tryCatch({
      tibble(
        pmid = article$uid %||% NA_character_,
        title = article$title %||% NA_character_,
        authors = paste(article$authors$name, collapse = "; ") %||% NA_character_,
        journal = article$source %||% NA_character_,
        pub_date = article$pubdate %||% NA_character_,
        doi = article$elocationid %||% 
          article$articleids$value[article$articleids$idtype == "doi"] %||% 
          NA_character_
      )
    }, error = function(e) {
      tibble(
        pmid = NA_character_,
        title = NA_character_,
        authors = NA_character_,
        journal = NA_character_,
        pub_date = NA_character_,
        doi = NA_character_
      )
    })
  }) %>%
    filter(!is.na(pmid))
  
  return(articles_df)
}

# Function: Extract genes from abstracts ----

extract_genes_from_abstracts <- function(pmids, batch_size = 100) {
  
  num_batches <- ceiling(length(pmids) / batch_size)
  all_genes <- list()
  
  # Common cancer genes pattern
  gene_pattern <- paste0(
    "\\b(",
    paste(c(
      "TP53", "BRCA1", "BRCA2", "KRAS", "NRAS", "BRAF", "EGFR", "PIK3CA", 
      "PTEN", "APC", "ATM", "CHEK2", "PALB2", "MLH1", "MSH2", "MSH6", "PMS2",
      "ERBB2", "MET", "ALK", "ROS1", "CDKN2A", "CTNNB1", "VHL", "STK11",
      "SMAD4", "FBXW7", "FGFR1", "FGFR2", "FGFR3", "NOTCH1", "IDH1", "IDH2",
      "AR", "ESR1", "GATA3", "MAP3K1", "CDH1", "SPOP", "TERT", "AXIN1",
      "ID3", "TCF3", "NFE2L2", "EP300", "MYC", "MYCN", "ARID1A", "SMARCA4",
      "RET", "RB1", "NF1", "HRAS", "JAK2", "KIT", "PDGFRA", "FLT3"
    ), collapse = "|"),
    ")\\b"
  )
  
  for (i in 1:num_batches) {
    cat("Processing batch", i, "of", num_batches, "\n")
    
    start_idx <- (i - 1) * batch_size + 1
    end_idx <- min(i * batch_size, length(pmids))
    batch_ids <- pmids[start_idx:end_idx]
    
    # Fetch abstracts as XML so we can parse per-article
    abstract_xml <- tryCatch({
      entrez_fetch(
        db = "pubmed",
        id = batch_ids,
        rettype = "xml",
        retmode = "xml"
      )
    }, error = function(e) {
      cat("Error fetching batch:", e$message, "\n")
      return(NULL)
    })
    
    if (is.null(abstract_xml)) {
      next
    }
    
    # Parse the XML
    xml_doc <- xml2::read_xml(abstract_xml)
    articles <- xml2::xml_find_all(xml_doc, ".//PubmedArticle")
    
    batch_results <- map_dfr(articles, function(article) {
      # Extract PMID for this specific article
      article_pmid <- xml2::xml_text(
        xml2::xml_find_first(article, ".//PMID")
      )
      
      # Extract abstract text for this specific article
      abstract_nodes <- xml2::xml_find_all(
        article, ".//AbstractText"
      )
      abstract_text <- paste(
        xml2::xml_text(abstract_nodes), collapse = " "
      )
      
      # Also get the title
      title_node <- xml2::xml_find_first(article, ".//ArticleTitle")
      title_text <- xml2::xml_text(title_node)
      
      # Combine title and abstract for gene searching
      full_text <- paste(title_text, abstract_text)
      
      # Extract genes from THIS article only
      genes_found <- unique(toupper(
        unlist(str_extract_all(full_text, regex(gene_pattern, ignore_case = TRUE)))
      ))
      
      if (length(genes_found) == 0) {
        return(tibble(pmid = as.numeric(article_pmid), genes = NA_character_))
      }
      
      tibble(
        pmid = as.numeric(article_pmid),
        genes = genes_found
      )
    })
    
    all_genes[[i]] <- batch_results
    Sys.sleep(0.5)
  }
  
  genes_df <- bind_rows(all_genes) %>%
    filter(!is.na(genes), genes != "") %>%
    distinct()
  
  return(genes_df)
}


# Execute searches for TOP 15 cancer types from GLOBOCAN 2022 ----


# Top 15 cancer types ranked by ASR Incidence from GLOBOCAN 2022
cancer_types <- c(
  "breast cancer",              # 40.49
  "prostate cancer",            # 30.28
  "cervical cancer",            # 26.42
  "liver cancer",               # 8.52
  "colorectal cancer",          # 8.24
  "lung cancer",                # 6.28
  "ovarian cancer",             # 5.28
  "non-Hodgkin lymphoma",       # 4.99
  "bladder cancer",             # 4.66
  "stomach cancer",             # 4.02
  "oesophageal cancer",         # 3.64
  "corpus uteri cancer",        # 3.48
  "leukaemia",                  # 3.10
  "pancreatic cancer",          # 2.35
  "Kaposi sarcoma"              # 2.25
)

cat("Starting PubMed searches for African cancer genomics...\n")
cat("Using GLOBOCAN 2022 (version 1.1) - Top 15 cancer types by ASR Incidence\n\n")

all_articles <- map_dfr(cancer_types, function(cancer) {
  cat("\n")
  print_separator()
  cat("Searching for:", cancer, "\n")
  print_separator()
  cat("\n")
  
  results <- search_pubmed_african_cancer(
    cancer_type = cancer,
    max_results = 500,
    year_start = 2010
  )
  
  if (!is.null(results)) {
    results$cancer_type <- cancer
  }
  
  Sys.sleep(1)
  return(results)
}) %>%
  distinct(pmid, .keep_all = TRUE)

# Save results
write_csv(all_articles, here("data", "processed", "pubmed_african_cancer_articles.csv"))

cat("\nTotal unique articles found:", nrow(all_articles), "\n")


# Extract genes mentioned in articles ----


if (nrow(all_articles) > 0) {
  cat("\nExtracting gene mentions from abstracts...\n")
  
  genes_in_literature <- extract_genes_from_abstracts(all_articles$pmid)
  
  write_csv(genes_in_literature, here("data", "processed", "genes_in_literature.csv"))
  
  # Summarize gene frequencies
  gene_summary <- genes_in_literature %>%
    count(genes, sort = TRUE) %>%
    rename(literature_count = n)
  
  write_csv(gene_summary, here("data", "processed", "gene_literature_frequency.csv"))
  
  cat("\nGenes extracted from", nrow(genes_in_literature), "abstracts\n")
  cat("Unique genes found:", nrow(gene_summary), "\n")
  
  # Display top genes
  cat("\nTop 20 genes mentioned in literature:\n")
  print(gene_summary %>% head(20))
  
} else {
  cat("\nNo articles found. Creating empty gene list...\n")
  
  # Create empty files so downstream scripts don't fail
  tibble(
    pmid = character(),
    genes = character()
  ) %>%
    write_csv(here("data", "processed", "genes_in_literature.csv"))
  
  tibble(
    genes = character(),
    literature_count = integer()
  ) %>%
    write_csv(here("data", "processed", "gene_literature_frequency.csv"))
}


# Load GLOBOCAN 2022 data from provided file ----


load_globocan_2022_data <- function() {
  
  cat("\nLoading GLOBOCAN 2022 data for Africa...\n")
  cat("Source: https://gco.iarc.who.int\n")
  cat("Dataset: GLOBOCAN 2022 (version 1.1) - 08.02.2024\n\n")
  
  # Top 15 cancer types by ASR Incidence with actual GLOBOCAN 2022 data
  cancer_incidence_africa <- tibble(
    rank = 1:15,
    cancer_type = c(
      "Breast",
      "Prostate", 
      "Cervix uteri",
      "Liver and intrahepatic bile ducts",
      "Colorectum",
      "lung",
      "Ovary",
      "Non-Hodgkin lymphoma",
      "Bladder",
      "Stomach",
      "Oesophagus",
      "Corpus uteri",
      "Leukaemia",
      "Pancreas",
      "Kaposi sarcoma"
    ),
    cancer_code = c(20, 27, 23, 11, 41, 15, 25, 34, 30, 7, 6, 24, 36, 13, 19),
    asr_incidence = c(40.49, 30.28, 26.42, 8.52, 8.24, 6.28, 5.28, 4.99, 4.66, 
                      4.02, 3.64, 3.48, 3.10, 2.35, 2.25),
    crude_rate_incidence = c(28.21, 14.66, 17.87, 5.25, 5.01, 3.54, 3.66, 3.59, 
                             2.63, 2.37, 2.13, 2.12, 2.35, 1.35, 1.89),
    cumulative_risk_incidence = c(4.30, 3.67, 2.91, 0.97, 0.94, 0.75, 0.58, 0.51, 
                                  0.55, 0.46, 0.43, 0.44, 0.30, 0.28, 0.20),
    total_cases = c(198342, 102977, 125615, 73759, 70362, 49765, 25736, 50449, 
                    37012, 33317, 29937, 14891, 32971, 18980, 26524),
    asr_mortality = c(19.16, 17.32, 17.58, 8.19, 5.55, 5.77, 3.97, 3.26, 2.74, 
                      3.50, 3.48, 1.13, 2.54, 2.23, 1.16),
    crude_rate_mortality = c(12.97, 7.93, 11.46, 5.00, 3.28, 3.23, 2.56, 2.19, 
                             1.47, 2.04, 2.01, 0.66, 1.74, 1.26, 0.99),
    cumulative_risk_mortality = c(2.05, 1.71, 1.99, 0.93, 0.59, 0.68, 0.46, 0.33, 
                                  0.28, 0.40, 0.41, 0.14, 0.26, 0.26, 0.10),
    total_deaths = c(91173, 55711, 80558, 70231, 46061, 45402, 18010, 30739, 
                     20716, 28700, 28251, 4662, 24487, 17758, 13908),
    mortality_to_incidence_ratio = round(total_deaths / total_cases, 3),
    priority_score = c(10, 10, 10, 9, 8, 7, 7, 6, 6, 6, 6, 5, 5, 5, 5),
    notes = c(
      "Leading cancer in women; highest incidence",
      "Leading cancer in men; second highest incidence",
      "Third highest burden; high mortality",
      "HBV/HCV endemic; high mortality rate",
      "Rising incidence; significant burden",
      "Increasing rates; respiratory cancers",
      "Gynecological cancer; moderate burden",
      "Hematological malignancy",
      "Urological cancer; moderate incidence",
      "H. pylori endemic; digestive system",
      "Endemic in East Africa; high mortality",
      "Gynecological cancer; lower mortality",
      "Hematological malignancy; blood cancers",
      "Highly lethal; digestive system",
      "HIV-associated; infectious etiology"
    ),
    search_term = c(
      "breast cancer",
      "prostate cancer",
      "cervical cancer",
      "liver cancer",
      "colorectal cancer",
      "lung cancer",
      "ovarian cancer",
      "non-Hodgkin lymphoma",
      "bladder cancer",
      "stomach cancer",
      "oesophageal cancer",
      "corpus uteri cancer",
      "leukaemia",
      "pancreatic cancer",
      "Kaposi sarcoma"
    )
  ) %>%
    mutate(
      data_source = "GLOBOCAN 2022 (v1.1)",
      extraction_date = "2024-02-08",
      region = "Africa"
    )
  
  return(cancer_incidence_africa)
}

cancer_incidence <- load_globocan_2022_data()

# Save the data
write_csv(cancer_incidence, here("data", "processed", "cancer_incidence_africa_globocan2022.csv"))

cat("\nCancer incidence data loaded and saved\n")
cat("Data source: GLOBOCAN 2022 (version 1.1) - 08.02.2024\n")
cat("Region: Africa\n")
cat("Cancer types included: Top 15 by ASR Incidence\n\n")

# Display summary
cat("Summary of Top 15 Cancers in Africa:\n")
print_separator("-")
print(cancer_incidence %>% 
        select(rank, cancer_type, asr_incidence, asr_mortality, total_cases, total_deaths) %>%
        mutate(across(where(is.numeric), ~round(., 2))))
