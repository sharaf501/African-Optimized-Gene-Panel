# African-targeted-Cancer-Panel

# Introduction

## Background

Africa bears 15% of the global cancer burden but has minimal representation in cancer genomics databases. This disparity creates challenges for precision oncology implementation in African healthcare systems.

## Objectives

1. Design a cost-effective targeted gene panel
2. Maximize clinical actionability for African populations
3. Include pan-cancer and African-specific driver genes
4. Ensure applicability across resource-limited settings

## These script presents the design and validation of an African-specific targeted cancer gene panel optimized for resource-limited settings. 

## Notes:
- Pan-Cancer Score: The average of (a) frequency score, defined as the gene's pan-cancer mutation frequency divided by the maximum observed frequency, and (b) breadth score, defined as the proportion of the 12 cancer types in which the gene harbors non-silent mutations.
- African Relevance Score (30%): Mutation frequencies across TCGA cancer types were weighted by GLOBOCAN 2022 African priority scores (Table 1). The weighted frequencies were summed per gene and converted to percentile ranks to yield a 0–1 scale.
- Literature Support Score (20%): The normalized frequency of gene mentions across 1,349 African cancer genomics articles retrieved through per-abstract PubMed mining, with each gene's article count divided by the maximum observed count.
- Expert Curation Score (20%): A binary score reflecting presence in the combined Foundation Medicine CDx and MSK-IMPACT gene list (1 if present; 0 if absent).

## Limitations

1. **Data Scarcity**: Limited African genomic data necessitated TCGA and AACR Project GENIE extrapolation
2. **Ethnic Diversity**: Africa's genetic diversity may require regional adaptations
3. **Validation Needed**: Clinical validation in African cohorts required
4. **Actionability Gap**: Limited access to targeted therapies in some regions

# Conclusions

This study presents a rigorously designed, evidence-based cancer gene panel optimized for African populations. The 70-Gene panel offers:
✓ Comprehensive coverage of clinically actionable mutations
✓ African cancer-relevant gene selection
✓ Cost-effective implementation
✓ Pan-cancer applicability
✓ Pharmacogenomic variant inclusion

**Recommendation**: Implementation of the 70-Gene panel for clinical use in African cancer centers, with the 130-gene panel for research applications.
