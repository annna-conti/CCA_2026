# ============================================================
# CCA vs WATER INDICATOR ANALYSIS
# ============================================================

library(phyloseq)
library(vegan)
library(indicspecies)
library(permute)
library(dplyr)
library(tidyr)
library(ggplot2)
library(ggrepel)

# ------------------------------------------------------------
# Load cleaned biological-sample phyloseq object
# ------------------------------------------------------------

ps_analysis <- readRDS(
  "../annabelle_16S_analysis_no_controls_phyloseq.rds"
)

ps_analysis
phyloseq::nsamples(ps_analysis)
phyloseq::ntaxa(ps_analysis)


# ------------------------------------------------------------
# Extract abundance matrix/recreate relative abundance matrix
# ------------------------------------------------------------

beta_mat <- as(
  phyloseq::otu_table(ps_analysis),
  "matrix"
)

# vegan expects samples = rows and taxa = columns
if (phyloseq::taxa_are_rows(ps_analysis)) {
  beta_mat <- t(beta_mat)
}

# Convert each sample to relative abundance
beta_rel <- sweep(
  beta_mat,
  1,
  rowSums(beta_mat),
  "/"
)

dim(beta_rel)


# ------------------------------------------------------------
# Extract metadata
# ------------------------------------------------------------

meta_df <- as.data.frame(
  phyloseq::sample_data(ps_analysis)
)

meta_df$SampleID <- rownames(meta_df)

table(meta_df$Environment)


# ------------------------------------------------------------
# Extract bacterial taxonomy
# ------------------------------------------------------------

tax_df <- as.data.frame(
  as(
    phyloseq::tax_table(ps_analysis),
    "matrix"
  )
)

tax_df$Taxon <- rownames(tax_df)


# ============================================================
# CCA vs WATER INDICATOR ANALYSIS
# ============================================================

# Create grouping variable
cca_water_group <- factor(
  meta_df$Environment,
  levels = c("CCA", "Water")
)

names(cca_water_group) <- meta_df$SampleID

# Align group order with abundance matrix
cca_water_group <- cca_water_group[
  rownames(beta_rel)
]

table(cca_water_group)


# ------------------------------------------------------------
# Prevalence filter
# Keep taxa present in at least 3 biological samples
# ------------------------------------------------------------

taxon_prevalence <- colSums(
  beta_rel > 0
)

beta_cca_water_ind <- beta_rel[
  ,
  taxon_prevalence >= 3,
  drop = FALSE
]

dim(beta_cca_water_ind)



# ------------------------------------------------------------
# Indicator value analysis
# ------------------------------------------------------------

set.seed(123)

indval_cca_water <- indicspecies::multipatt(
  beta_cca_water_ind,
  cca_water_group,
  func = "IndVal.g",
  duleg = TRUE,
  control = permute::how(
    nperm = 999
  )
)

summary(
  indval_cca_water,
  indvalcomp = TRUE,
  alpha = 0.05
)


# ============================================================
# Extract significant CCA indicator taxa
# ============================================================

indval_cca_water_raw <- indval_cca_water$sign

indval_cca_water_raw$Taxon <- rownames(
  indval_cca_water_raw
)

# Attach taxonomy
indval_cca_water_tax <- dplyr::left_join(
  indval_cca_water_raw,
  tax_df,
  by = "Taxon"
)

# Pull only significant CCA indicators
CCA_indicators <- indval_cca_water_tax %>%
  dplyr::filter(
    s.CCA == 1,
    p.value < 0.05
  ) %>%
  dplyr::select(
    Taxon,
    phylum,
    class,
    order,
    family,
    genus,
    species,
    stat,
    p.value
  ) %>%
  dplyr::arrange(
    dplyr::desc(stat)
  )

nrow(CCA_indicators)

CCA_indicators


# ============================================================
# Add abundance and prevalence information
# ============================================================

cca_ids <- meta_df$SampleID[
  meta_df$Environment == "CCA"
]

water_ids <- meta_df$SampleID[
  meta_df$Environment == "Water"
]

CCA_indicator_summary <- CCA_indicators %>%
  dplyr::rowwise() %>%
  dplyr::mutate(
    
    Mean_CCA_Abundance =
      mean(beta_rel[cca_ids, Taxon]),
    
    Mean_Water_Abundance =
      mean(beta_rel[water_ids, Taxon]),
    
    CCA_Prevalence =
      sum(beta_rel[cca_ids, Taxon] > 0),
    
    Water_Prevalence =
      sum(beta_rel[water_ids, Taxon] > 0)
    
  ) %>%
  dplyr::ungroup() %>%
  dplyr::mutate(
    
    Mean_CCA_Percent =
      Mean_CCA_Abundance * 100,
    
    Mean_Water_Percent =
      Mean_Water_Abundance * 100
    
  ) %>%
  dplyr::arrange(
    dplyr::desc(stat)
  )


CCA_indicator_summary %>%
  dplyr::select(
    Taxon,
    genus,
    species,
    stat,
    p.value,
    CCA_Prevalence,
    Water_Prevalence,
    Mean_CCA_Percent,
    Mean_Water_Percent
  ) %>%
  print(n = 30)


#narrowing it down from 90 taxa --> rank them by mean abundance on CCA
# ============================================================
# Prioritize strong, abundant CCA indicator taxa
# ============================================================

CCA_indicators_priority <- CCA_indicator_summary %>%
  dplyr::filter(
    stat >= 0.95
  ) %>%
  dplyr::arrange(
    dplyr::desc(Mean_CCA_Percent)
  ) %>%
  dplyr::slice_head(
    n = 15
  )

CCA_indicators_priority %>%
  dplyr::select(
    Taxon,
    genus,
    species,
    stat,
    p.value,
    CCA_Prevalence,
    Water_Prevalence,
    Mean_CCA_Percent,
    Mean_Water_Percent
  )
