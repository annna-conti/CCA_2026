#prepare the dataset for their workflow 
#the paper retained bacterial features that reached >0.01% relative abundance in at least one sample, then rarefied the filtered dataset to the depth of the smallest retained sample 
library(phyloseq)
library(vegan)
library(dplyr)
library(ggplot2)

# Start from cleaned, decontaminated biological dataset
ps_turnlund <- ps_analysis

# Keep only CCA samples
ps_turnlund_CCA <- phyloseq::prune_samples(
  phyloseq::sample_data(ps_turnlund)$Environment == "CCA",
  ps_turnlund
)

nsamples(ps_turnlund_CCA)

# Convert temporarily to relative abundance
ps_turnlund_rel <- phyloseq::transform_sample_counts(
  ps_turnlund_CCA,
  function(x) x / sum(x)
)

rel_mat <- as(
  phyloseq::otu_table(ps_turnlund_rel),
  "matrix"
)

# Keep taxa that exceed 0.01% relative abundance
# in at least one CCA sample
if (phyloseq::taxa_are_rows(ps_turnlund_rel)) {
  
  keep_taxa <- apply(
    rel_mat,
    1,
    function(x) any(x > 0.0001)
  )
  
} else {
  
  keep_taxa <- apply(
    rel_mat,
    2,
    function(x) any(x > 0.0001)
  )
}

ps_turnlund_CCA_filt <- phyloseq::prune_taxa(
  keep_taxa,
  ps_turnlund_CCA
)

ntaxa(ps_turnlund_CCA_filt)
#retained 693 taxa after doing this filtering 

# ------------------------------------------------------------
# 2. Prepare rounded copy for alpha-diversity rarefaction
# ------------------------------------------------------------

ps_turnlund_alpha <- ps_turnlund_CCA_filt

otu_alpha <- as(
  phyloseq::otu_table(ps_turnlund_alpha),
  "matrix"
)

# Round EMU estimated counts ONLY for rarefaction
otu_alpha <- round(otu_alpha)

phyloseq::otu_table(ps_turnlund_alpha) <- phyloseq::otu_table(
  otu_alpha,
  taxa_are_rows = phyloseq::taxa_are_rows(ps_turnlund_alpha)
)

# Remove any taxa that became zero everywhere after rounding
ps_turnlund_alpha <- phyloseq::prune_taxa(
  phyloseq::taxa_sums(ps_turnlund_alpha) > 0,
  ps_turnlund_alpha
)

#calculating sequencing depth of each sample
alpha_depths <- phyloseq::sample_sums(
  ps_turnlund_alpha
)

sort(alpha_depths)

#get rarefaction depths to use 
min_alpha_depth <- min(alpha_depths)

min_alpha_depth
#47,033 reads is the rarefaction depth I will use to rarefy all 14 CCA samples for alpha diversity calcs


#because EMU outputs fractional estimated counts, a rounded copy of the filtered 
#dataset was created only for rarefaction based alpha diversity calcs

#rarefy all samples to min depth 
set.seed(123)

ps_turnlund_alpha_rare <- phyloseq::rarefy_even_depth(
  ps_turnlund_alpha,
  sample.size = min_alpha_depth,
  rngseed = 123,
  replace = FALSE,
  verbose = TRUE
)
#original unrounded data will be used for nMDS, bray-curtis, permanova, etc
phyloseq::sample_sums(
  ps_turnlund_alpha_rare
)
phyloseq::nsamples(
  ps_turnlund_alpha_rare
)


#next step is to calculate observed richness and then evenness = Shannon / log(observed richness)
alpha_mat_turnlund <- as(
  phyloseq::otu_table(ps_turnlund_alpha_rare),
  "matrix"
)

if (phyloseq::taxa_are_rows(ps_turnlund_alpha_rare)) {
  alpha_mat_turnlund <- t(alpha_mat_turnlund)
}

alpha_turnlund <- data.frame(
  SampleID = rownames(alpha_mat_turnlund),
  
  Observed = vegan::specnumber(
    alpha_mat_turnlund
  ),
  
  Shannon = vegan::diversity(
    alpha_mat_turnlund,
    index = "shannon"
  )
)

alpha_turnlund$Evenness <- alpha_turnlund$Shannon /
  log(alpha_turnlund$Observed)

alpha_turnlund

#attach CCA metadata to alpha diversity table 
alpha_meta_turnlund <- data.frame(
  phyloseq::sample_data(ps_turnlund_alpha_rare),
  stringsAsFactors = FALSE,
  check.names = FALSE
)

alpha_meta_turnlund$SampleID <- rownames(alpha_meta_turnlund)

alpha_turnlund <- dplyr::left_join(
  alpha_turnlund,
  alpha_meta_turnlund,
  by = "SampleID"
)

alpha_turnlund %>%
  dplyr::select(
    SampleID,
    CCA_Morphology,
    Host_Genus,
    Host_Species,
    Observed,
    Shannon,
    Evenness
  )


#kruskal-wallis to see differences in observed richness and evenness
table(alpha_turnlund$CCA_Morphology)

#observed richness
kruskal_richness <- kruskal.test(
  Observed ~ CCA_Morphology,
  data = alpha_turnlund
)

kruskal_richness

#evenness
kruskal_evenness <- kruskal.test(
  Evenness ~ CCA_Morphology,
  data = alpha_turnlund
)

kruskal_evenness

#making the figures to show this - Plot richness and evenness by CCA morphology
library(tidyr)

alpha_plot_df <- alpha_turnlund %>%
  dplyr::select(
    SampleID,
    CCA_Morphology,
    Observed,
    Evenness
  ) %>%
  tidyr::pivot_longer(
    cols = c(Observed, Evenness),
    names_to = "Metric",
    values_to = "Value"
  )

alpha_plot_df$Metric <- factor(
  alpha_plot_df$Metric,
  levels = c("Observed", "Evenness"),
  labels = c("Observed richness", "Evenness")
)

#actual plot
alpha_turnlund_plot <- ggplot(
  alpha_plot_df,
  aes(
    x = CCA_Morphology,
    y = Value,
    color = CCA_Morphology
  )
) +
  
  geom_boxplot(
    aes(fill = CCA_Morphology),
    alpha = 0.25,
    outlier.shape = NA,
    width = 0.6
  ) +
  
  geom_jitter(
    width = 0.10,
    size = 3,
    alpha = 0.9
  ) +
  
  facet_wrap(
    ~ Metric,
    scales = "free_y",
    nrow = 1
  ) +
  
  scale_color_manual(
    values = viridisLite::viridis(
      2,
      option = "viridis",
      begin = 0.15,
      end = 0.85
    )
  ) +
  
  scale_fill_manual(
    values = viridisLite::viridis(
      2,
      option = "viridis",
      begin = 0.15,
      end = 0.85
    )
  ) +
  
  labs(
    x = "CCA morphology",
    y = NULL,
    color = "CCA morphology",
    fill = "CCA morphology"
  ) +
  
  theme_classic() +
  
  theme(
    legend.position = "none",
    strip.background = element_blank(),
    strip.text = element_text(
      face = "bold",
      size = 12
    ),
    axis.title.x = element_text(
      size = 11
    )
  )

alpha_turnlund_plot

# ------------------------------------------------------------
# Bray-Curtis NMDS on non-rarefied, log-transformed counts
# ------------------------------------------------------------

# Extract the filtered, NON-RAREFIED CCA count table
nmds_mat_turnlund <- as(
  phyloseq::otu_table(ps_turnlund_CCA_filt),
  "matrix"
)

# Make rows = samples and columns = taxa
if (phyloseq::taxa_are_rows(ps_turnlund_CCA_filt)) {
  nmds_mat_turnlund <- t(nmds_mat_turnlund)
}

dim(nmds_mat_turnlund)

# Log-transform counts
nmds_log_turnlund <- log1p(nmds_mat_turnlund)

set.seed(123)

nmds_turnlund <- vegan::metaMDS(
  nmds_log_turnlund,
  distance = "bray",
  k = 2,
  trymax = 200,
  autotransform = FALSE,
  trace = FALSE
)

nmds_turnlund
#stress = 0.0926

# ------------------------------------------------------------
#Extract NMDS scores and attach metadata
# ------------------------------------------------------------

nmds_scores_turnlund <- as.data.frame(
  vegan::scores(
    nmds_turnlund,
    display = "sites"
  )
)

nmds_scores_turnlund$SampleID <- rownames(
  nmds_scores_turnlund
)

nmds_meta_turnlund <- data.frame(
  phyloseq::sample_data(ps_turnlund_CCA_filt),
  stringsAsFactors = FALSE,
  check.names = FALSE
)

nmds_meta_turnlund$SampleID <- rownames(
  nmds_meta_turnlund
)

nmds_scores_turnlund <- dplyr::left_join(
  nmds_scores_turnlund,
  nmds_meta_turnlund,
  by = "SampleID"
)

nmds_scores_turnlund %>%
  dplyr::select(
    SampleID,
    NMDS1,
    NMDS2,
    CCA_Morphology,
    Host_Genus,
    Host_Species
  )

##nmds colored by morophology##

morph_colors <- viridisLite::viridis(
  n = 2,
  option = "viridis",
  begin = 0.15,
  end = 0.85
)

names(morph_colors) <- c(
  "branching",
  "encrusting"
)

nmds_morph_plot <- ggplot(
  nmds_scores_turnlund,
  aes(
    x = NMDS1,
    y = NMDS2,
    color = CCA_Morphology
  )
) +
  
  geom_point(
    size = 4,
    alpha = 0.9
  ) +
  
  geom_text(
    aes(label = SampleID),
    vjust = -0.8,
    size = 3.2,
    show.legend = FALSE
  ) +
  
  scale_color_manual(
    values = morph_colors
  ) +
  
  labs(
    title = "CCA bacterial community structure",
    subtitle = "Bray-Curtis NMDS of log-transformed bacterial abundances",
    x = "NMDS1",
    y = "NMDS2",
    color = "CCA morphology"
  ) +
  
  theme_classic() +
  
  theme(
    plot.title = element_text(
      face = "bold"
    ),
    legend.position = "right"
  )

nmds_morph_plot

# ------------------------------------------------------------
# 10. NMDS plot colored by CCA host identity
# ------------------------------------------------------------

nmds_scores_turnlund <- nmds_scores_turnlund %>%
  dplyr::mutate(
    Host_Label = dplyr::case_when(
      
      !is.na(Host_Species) &
        Host_Species != "" &
        Host_Species != "NF" ~ Host_Species,
      
      !is.na(Host_Genus) &
        Host_Genus != "" &
        Host_Genus != "NF" ~ Host_Genus,
      
      TRUE ~ "Unidentified"
    )
  )

table(nmds_scores_turnlund$Host_Label)


host_levels <- unique(nmds_scores_turnlund$Host_Label)

host_colors <- viridisLite::viridis(
  n = length(host_levels),
  option = "viridis",
  begin = 0.05,
  end = 0.95
)

names(host_colors) <- host_levels

nmds_host_plot <- ggplot(
  nmds_scores_turnlund,
  aes(
    x = NMDS1,
    y = NMDS2,
    color = Host_Label,
    shape = CCA_Morphology
  )
) +
  
  geom_point(
    size = 4,
    alpha = 0.9
  ) +
  
  geom_text(
    aes(label = SampleID),
    vjust = -0.8,
    size = 3.2,
    show.legend = FALSE
  ) +
  
  scale_color_manual(
    values = host_colors
  ) +
  
  labs(
    title = "CCA bacterial community structure by host identity",
    subtitle = "Bray-Curtis NMDS of log-transformed bacterial abundances",
    x = "NMDS1",
    y = "NMDS2",
    color = "CCA host",
    shape = "CCA morphology"
  ) +
  
  theme_classic() +
  
  theme(
    plot.title = element_text(
      face = "bold"
    ),
    legend.position = "right"
  )

nmds_host_plot


# ------------------------------------------------------------
# Beta dispersion by CCA morphology
# ------------------------------------------------------------

# Bray-Curtis distance matrix from log-transformed counts
bray_turnlund <- vegan::vegdist(
  nmds_log_turnlund,
  method = "bray"
)

# Make sure morphology is aligned with the matrix rows
morph_group_turnlund <- nmds_scores_turnlund$CCA_Morphology

names(morph_group_turnlund) <- nmds_scores_turnlund$SampleID

morph_group_turnlund <- morph_group_turnlund[
  rownames(nmds_log_turnlund)
]

# Test dispersion
disp_morph_turnlund <- vegan::betadisper(
  bray_turnlund,
  group = morph_group_turnlund,
  bias.adjust = TRUE
)

disp_morph_turnlund

#test whether the mean distance to centroid differs
anova(disp_morph_turnlund)


# ------------------------------------------------------------
# PERMANOVA by CCA morphology
# ------------------------------------------------------------

permanova_meta_turnlund <- nmds_scores_turnlund %>%
  dplyr::select(
    SampleID,
    CCA_Morphology
  )

# Put metadata in exactly the same order as the community matrix
permanova_meta_turnlund <- permanova_meta_turnlund[
  match(
    rownames(nmds_log_turnlund),
    permanova_meta_turnlund$SampleID
  ),
]

rownames(permanova_meta_turnlund) <-
  permanova_meta_turnlund$SampleID

# Check alignment
all(
  rownames(permanova_meta_turnlund) ==
    rownames(nmds_log_turnlund)
)

set.seed(123)

permanova_morph_turnlund <- vegan::adonis2(
  bray_turnlund ~ CCA_Morphology,
  data = permanova_meta_turnlund,
  permutations = 999
)

permanova_morph_turnlund


# ------------------------------------------------------------
# Pairwise beta diversity using Whittaker's index
# ------------------------------------------------------------

whittaker_turnlund <- vegan::betadiver(
  nmds_log_turnlund,
  method = "w"
)

whittaker_turnlund




# ------------------------------------------------------------
# Define core100 and core75 CCA bacterial communities
# ------------------------------------------------------------

core_mat_turnlund <- as(
  phyloseq::otu_table(ps_turnlund_CCA_filt),
  "matrix"
)

# Make rows = samples and columns = taxa
if (phyloseq::taxa_are_rows(ps_turnlund_CCA_filt)) {
  core_mat_turnlund <- t(core_mat_turnlund)
}

# Number of CCA samples in which each taxon is present
core_prevalence <- colSums(core_mat_turnlund > 0)

# Total number of CCA samples
n_CCA <- nrow(core_mat_turnlund)

n_CCA

# Present in ALL 14 CCA samples
core100_ids <- names(
  core_prevalence[core_prevalence == n_CCA]
)

# Present in at least 75% of CCA samples
# 75% of 14 = 10.5, so this requires at least 11 samples
core75_ids <- names(
  core_prevalence[core_prevalence >= ceiling(0.75 * n_CCA)]
)

length(core100_ids)
length(core75_ids)

#adding in a water enrichment requirement - which core taxa are enriched relative to seawater
# ------------------------------------------------------------
# Compare core100 CCA taxa with seawater
# ------------------------------------------------------------

# Convert all biological samples to relative abundance
ps_rel_turnlund <- phyloseq::transform_sample_counts(
  ps_analysis,
  function(x) x / sum(x)
)

rel_all_turnlund <- as(
  phyloseq::otu_table(ps_rel_turnlund),
  "matrix"
)

# Make rows = samples, columns = taxa
if (phyloseq::taxa_are_rows(ps_rel_turnlund)) {
  rel_all_turnlund <- t(rel_all_turnlund)
}

# Metadata
meta_all_turnlund <- data.frame(
  phyloseq::sample_data(ps_rel_turnlund),
  stringsAsFactors = FALSE,
  check.names = FALSE
)

meta_all_turnlund$SampleID <- rownames(meta_all_turnlund)

table(meta_all_turnlund$Environment)


core100_rel <- rel_all_turnlund[
  ,
  core100_ids,
  drop = FALSE
]


core100_long <- as.data.frame(core100_rel) %>%
  tibble::rownames_to_column("SampleID") %>%
  tidyr::pivot_longer(
    cols = -SampleID,
    names_to = "Taxon",
    values_to = "RelativeAbundance"
  ) %>%
  dplyr::left_join(
    meta_all_turnlund %>%
      dplyr::select(
        SampleID,
        Environment
      ),
    by = "SampleID"
  )

core100_water_tests <- core100_long %>%
  dplyr::group_by(Taxon) %>%
  dplyr::summarise(
    
    Mean_CCA = mean(
      RelativeAbundance[Environment == "CCA"]
    ),
    
    Mean_Water = mean(
      RelativeAbundance[Environment == "Water"]
    ),
    
    Median_CCA = median(
      RelativeAbundance[Environment == "CCA"]
    ),
    
    Median_Water = median(
      RelativeAbundance[Environment == "Water"]
    ),
    
    p_value = wilcox.test(
      RelativeAbundance[Environment == "CCA"],
      RelativeAbundance[Environment == "Water"],
      exact = FALSE
    )$p.value,
    
    .groups = "drop"
  ) %>%
  
  dplyr::mutate(
    Mean_CCA_Percent = Mean_CCA * 100,
    Mean_Water_Percent = Mean_Water * 100,
    
    CCA_enriched = Mean_CCA > Mean_Water,
    
    Significant = p_value < 0.05
  )

core100_water_tests %>%
  dplyr::filter(
    CCA_enriched,
    Significant
  ) %>%
  nrow()
#47 taxa


core100_CCA_enriched <- core100_water_tests %>%
  dplyr::filter(
    CCA_enriched,
    Significant
  ) %>%
  dplyr::arrange(p_value)

core100_CCA_enriched


# ------------------------------------------------------------
# Attach taxonomy to CCA-enriched core100 taxa
# ------------------------------------------------------------

tax_turnlund <- as.data.frame(
  phyloseq::tax_table(ps_analysis),
  stringsAsFactors = FALSE
)

tax_turnlund$Taxon <- rownames(tax_turnlund)

core100_CCA_enriched_tax <- core100_CCA_enriched %>%
  dplyr::left_join(
    tax_turnlund,
    by = "Taxon"
  )

core100_CCA_enriched_tax <- core100_CCA_enriched_tax %>%
  dplyr::mutate(
    
    Fold_Enrichment = dplyr::if_else(
      Mean_Water > 0,
      Mean_CCA / Mean_Water,
      Inf
    )
    
  ) %>%
  dplyr::arrange(
    dplyr::desc(Mean_CCA)
  )

core100_CCA_enriched_tax %>%
  dplyr::select(
    Taxon,
    family,
    genus,
    species,
    Mean_CCA_Percent,
    Mean_Water_Percent,
    Fold_Enrichment,
    p_value
  )



# ------------------------------------------------------------
# Compare core75 CCA taxa with seawater
# ------------------------------------------------------------

core75_rel <- rel_all_turnlund[
  ,
  core75_ids,
  drop = FALSE
]

core75_long <- as.data.frame(core75_rel) %>%
  tibble::rownames_to_column("SampleID") %>%
  tidyr::pivot_longer(
    cols = -SampleID,
    names_to = "Taxon",
    values_to = "RelativeAbundance"
  ) %>%
  dplyr::left_join(
    meta_all_turnlund %>%
      dplyr::select(
        SampleID,
        Environment
      ),
    by = "SampleID"
  )

core75_water_tests <- core75_long %>%
  dplyr::group_by(Taxon) %>%
  dplyr::summarise(
    
    Mean_CCA = mean(
      RelativeAbundance[Environment == "CCA"]
    ),
    
    Mean_Water = mean(
      RelativeAbundance[Environment == "Water"]
    ),
    
    Median_CCA = median(
      RelativeAbundance[Environment == "CCA"]
    ),
    
    Median_Water = median(
      RelativeAbundance[Environment == "Water"]
    ),
    
    p_value = wilcox.test(
      RelativeAbundance[Environment == "CCA"],
      RelativeAbundance[Environment == "Water"],
      exact = FALSE
    )$p.value,
    
    .groups = "drop"
  ) %>%
  
  dplyr::mutate(
    Mean_CCA_Percent = Mean_CCA * 100,
    Mean_Water_Percent = Mean_Water * 100,
    
    CCA_enriched = Mean_CCA > Mean_Water,
    
    Significant = p_value < 0.05
  )

core75_CCA_enriched <- core75_water_tests %>%
  dplyr::filter(
    CCA_enriched,
    Significant
  ) %>%
  dplyr::arrange(p_value)

nrow(core75_CCA_enriched)


# ------------------------------------------------------------
# Identify core75-enriched taxa not already in core100
# ------------------------------------------------------------

core75_only_enriched <- core75_CCA_enriched %>%
  dplyr::filter(
    !Taxon %in% core100_CCA_enriched$Taxon
  )

nrow(core75_only_enriched)

core75_only_enriched_tax <- core75_only_enriched %>%
  dplyr::left_join(
    tax_turnlund,
    by = "Taxon"
  ) %>%
  dplyr::mutate(
    Fold_Enrichment = dplyr::if_else(
      Mean_Water > 0,
      Mean_CCA / Mean_Water,
      Inf
    )
  ) %>%
  dplyr::arrange(
    dplyr::desc(Mean_CCA)
  )

core75_only_enriched_tax %>%
  dplyr::select(
    Taxon,
    family,
    genus,
    species,
    Mean_CCA_Percent,
    Mean_Water_Percent,
    Fold_Enrichment,
    p_value
  )


# ------------------------------------------------------------
# Core100 Bray-Curtis NMDS
# ------------------------------------------------------------

# Keep only the 133 core100 taxa
ps_core100_turnlund <- phyloseq::prune_taxa(
  core100_ids,
  ps_turnlund_CCA_filt
)

phyloseq::ntaxa(ps_core100_turnlund)
phyloseq::nsamples(ps_core100_turnlund)


core100_mat_turnlund <- as(
  phyloseq::otu_table(ps_core100_turnlund),
  "matrix"
)

# Make rows = samples, columns = taxa
if (phyloseq::taxa_are_rows(ps_core100_turnlund)) {
  core100_mat_turnlund <- t(core100_mat_turnlund)
}

core100_log_turnlund <- log1p(
  core100_mat_turnlund
)

#bray-curtis nmds
set.seed(123)

nmds_core100_turnlund <- vegan::metaMDS(
  core100_log_turnlund,
  distance = "bray",
  k = 2,
  trymax = 200,
  autotransform = FALSE,
  trace = FALSE
)

nmds_core100_turnlund


#making the plot 
# ------------------------------------------------------------
# Plot core100 NMDS by morphology
# ------------------------------------------------------------

core100_scores <- as.data.frame(
  vegan::scores(
    nmds_core100_turnlund,
    display = "sites"
  )
)

core100_scores$SampleID <- rownames(core100_scores)

core100_scores <- dplyr::left_join(
  core100_scores,
  nmds_meta_turnlund,
  by = "SampleID"
)

core100_morph_plot <- ggplot(
  core100_scores,
  aes(
    x = NMDS1,
    y = NMDS2,
    color = CCA_Morphology
  )
) +
  
  geom_point(
    size = 4,
    alpha = 0.9
  ) +
  
  geom_text(
    aes(label = SampleID),
    vjust = -0.8,
    size = 3.2,
    show.legend = FALSE
  ) +
  
  scale_color_manual(
    values = morph_colors
  ) +
  
  labs(
    title = "Core100 CCA bacterial community structure",
    subtitle = "Bray-Curtis NMDS of log-transformed core100 bacterial abundances",
    x = "NMDS1",
    y = "NMDS2",
    color = "CCA morphology"
  ) +
  
  theme_classic() +
  
  theme(
    plot.title = element_text(face = "bold"),
    legend.position = "right"
  )

core100_morph_plot

#plot by host
core100_scores <- core100_scores %>%
  dplyr::mutate(
    Host_Label = dplyr::case_when(
      !is.na(Host_Species) &
        Host_Species != "" &
        Host_Species != "NF" ~ Host_Species,
      
      !is.na(Host_Genus) &
        Host_Genus != "" &
        Host_Genus != "NF" ~ Host_Genus,
      
      TRUE ~ "Unidentified"
    )
  )

host_levels_core100 <- unique(
  core100_scores$Host_Label
)

host_colors_core100 <- viridisLite::viridis(
  n = length(host_levels_core100),
  option = "viridis",
  begin = 0.05,
  end = 0.95
)

names(host_colors_core100) <- host_levels_core100

core100_host_plot <- ggplot(
  core100_scores,
  aes(
    x = NMDS1,
    y = NMDS2,
    color = Host_Label,
    shape = CCA_Morphology
  )
) +
  
  geom_point(
    size = 4,
    alpha = 0.9
  ) +
  
  geom_text(
    aes(label = SampleID),
    vjust = -0.8,
    size = 3.2,
    show.legend = FALSE
  ) +
  
  scale_color_manual(
    values = host_colors_core100
  ) +
  
  labs(
    title = "Core100 CCA bacterial community structure by host identity",
    subtitle = "Bray-Curtis NMDS of log-transformed core100 bacterial abundances",
    x = "NMDS1",
    y = "NMDS2",
    color = "CCA host",
    shape = "CCA morphology"
  ) +
  
  theme_classic()

core100_host_plot



# ------------------------------------------------------------
# Core75 Bray-Curtis NMDS
# ------------------------------------------------------------

ps_core75_turnlund <- phyloseq::prune_taxa(
  core75_ids,
  ps_turnlund_CCA_filt
)

phyloseq::ntaxa(ps_core75_turnlund)
phyloseq::nsamples(ps_core75_turnlund)


core75_mat_turnlund <- as(
  phyloseq::otu_table(ps_core75_turnlund),
  "matrix"
)

if (phyloseq::taxa_are_rows(ps_core75_turnlund)) {
  core75_mat_turnlund <- t(core75_mat_turnlund)
}


core75_log_turnlund <- log1p(
  core75_mat_turnlund
)

set.seed(123)

nmds_core75_turnlund <- vegan::metaMDS(
  core75_log_turnlund,
  distance = "bray",
  k = 2,
  trymax = 200,
  autotransform = FALSE,
  trace = FALSE
)

nmds_core75_turnlund
#good stress value

# ------------------------------------------------------------
# Plot core75 NMDS by morphology
# ------------------------------------------------------------

core75_scores <- as.data.frame(
  vegan::scores(
    nmds_core75_turnlund,
    display = "sites"
  )
)

core75_scores$SampleID <- rownames(core75_scores)

core75_scores <- dplyr::left_join(
  core75_scores,
  nmds_meta_turnlund,
  by = "SampleID"
)

core75_morph_plot <- ggplot(
  core75_scores,
  aes(
    x = NMDS1,
    y = NMDS2,
    color = CCA_Morphology
  )
) +
  
  geom_point(
    size = 4,
    alpha = 0.9
  ) +
  
  geom_text(
    aes(label = SampleID),
    vjust = -0.8,
    size = 3.2,
    show.legend = FALSE
  ) +
  
  scale_color_manual(
    values = morph_colors
  ) +
  
  labs(
    title = "Core75 CCA bacterial community structure",
    subtitle = "Bray-Curtis NMDS of log-transformed core75 bacterial abundances",
    x = "NMDS1",
    y = "NMDS2",
    color = "CCA morphology"
  ) +
  
  theme_classic() +
  
  theme(
    plot.title = element_text(face = "bold"),
    legend.position = "right"
  )

core75_morph_plot


#host version
core75_scores <- core75_scores %>%
  dplyr::mutate(
    Host_Label = dplyr::case_when(
      !is.na(Host_Species) &
        Host_Species != "" &
        Host_Species != "NF" ~ Host_Species,
      
      !is.na(Host_Genus) &
        Host_Genus != "" &
        Host_Genus != "NF" ~ Host_Genus,
      
      TRUE ~ "Unidentified"
    )
  )

host_levels_core75 <- unique(
  core75_scores$Host_Label
)

host_colors_core75 <- viridisLite::viridis(
  n = length(host_levels_core75),
  option = "viridis",
  begin = 0.05,
  end = 0.95
)

names(host_colors_core75) <- host_levels_core75

core75_host_plot <- ggplot(
  core75_scores,
  aes(
    x = NMDS1,
    y = NMDS2,
    color = Host_Label,
    shape = CCA_Morphology
  )
) +
  
  geom_point(
    size = 4,
    alpha = 0.9
  ) +
  
  geom_text(
    aes(label = SampleID),
    vjust = -0.8,
    size = 3.2,
    show.legend = FALSE
  ) +
  
  scale_color_manual(
    values = host_colors_core75
  ) +
  
  labs(
    title = "Core75 CCA bacterial community structure by host identity",
    subtitle = "Bray-Curtis NMDS of log-transformed core75 bacterial abundances",
    x = "NMDS1",
    y = "NMDS2",
    color = "CCA host",
    shape = "CCA morphology"
  ) +
  
  theme_classic()

core75_host_plot


# ------------------------------------------------------------
# Prepare CCA-enriched core100 taxa for visualization
# ------------------------------------------------------------

core100_enriched_ids <- core100_CCA_enriched_tax$Taxon

core100_enriched_mat <- rel_all_turnlund[
  rownames(rel_all_turnlund) %in%
    rownames(meta_all_turnlund)[
      meta_all_turnlund$Environment == "CCA"
    ],
  core100_enriched_ids,
  drop = FALSE
]

dim(core100_enriched_mat)

core100_labels <- core100_CCA_enriched_tax %>%
  dplyr::mutate(
    Taxon_Label = dplyr::case_when(
      
      !is.na(genus) &
        genus != "" &
        genus != "NF" ~ genus,
      
      !is.na(family) &
        family != "" &
        family != "NF" ~ paste0(family, " (family)"),
      
      TRUE ~ Taxon
    )
  ) %>%
  dplyr::select(
    Taxon,
    Taxon_Label
  )

core100_labels %>%
  dplyr::count(Taxon_Label) %>%
  dplyr::filter(n > 1)
#this shows there are 3 duplicated genus labels 
#for heatmap can combine those EMU taxa at genus level 
#if five distinct EMU taxa are all classified as Litoreibacter, the heatmap shows their 
#total relative abundance as Litoreibacter rather than five separate rows with the same name

# ------------------------------------------------------------
# Collapse enriched core100 taxa to interpretable labels
# ------------------------------------------------------------

core100_label_map <- core100_CCA_enriched_tax %>%
  dplyr::mutate(
    Taxon_Label = dplyr::case_when(
      
      !is.na(genus) &
        genus != "" &
        genus != "NF" ~ genus,
      
      !is.na(family) &
        family != "" &
        family != "NF" ~ paste0(family, " (family)"),
      
      TRUE ~ Taxon
    )
  ) %>%
  dplyr::select(
    Taxon,
    Taxon_Label
  )

core100_enriched_all <- rel_all_turnlund[
  ,
  core100_enriched_ids,
  drop = FALSE
]

core100_heatmap_long <- as.data.frame(
  core100_enriched_all
) %>%
  tibble::rownames_to_column("SampleID") %>%
  
  tidyr::pivot_longer(
    cols = -SampleID,
    names_to = "Taxon",
    values_to = "RelativeAbundance"
  ) %>%
  
  dplyr::left_join(
    core100_label_map,
    by = "Taxon"
  ) %>%
  
  dplyr::group_by(
    SampleID,
    Taxon_Label
  ) %>%
  
  dplyr::summarise(
    RelativeAbundance = sum(RelativeAbundance),
    .groups = "drop"
  )

dplyr::n_distinct(
  core100_heatmap_long$Taxon_Label
)


# ------------------------------------------------------------
# 25. Attach metadata and order samples/taxa for heatmap
# ------------------------------------------------------------

heatmap_meta <- meta_all_turnlund %>%
  dplyr::select(
    SampleID,
    Environment,
    CCA_Morphology,
    Host_Genus,
    Host_Species
  ) %>%
  dplyr::mutate(
    Sample_Order_Group = dplyr::case_when(
      Environment == "CCA" & CCA_Morphology == "encrusting" ~ 1,
      Environment == "CCA" & CCA_Morphology == "branching" ~ 2,
      Environment == "Water" ~ 3,
      TRUE ~ 4
    )
  )

taxon_order_df <- core100_heatmap_long %>%
  dplyr::left_join(
    heatmap_meta,
    by = "SampleID"
  ) %>%
  dplyr::filter(Environment == "CCA") %>%
  dplyr::group_by(Taxon_Label) %>%
  dplyr::summarise(
    Mean_CCA_Abundance = mean(RelativeAbundance),
    .groups = "drop"
  ) %>%
  dplyr::arrange(dplyr::desc(Mean_CCA_Abundance))

sample_order <- heatmap_meta %>%
  dplyr::arrange(
    Sample_Order_Group,
    Host_Genus,
    Host_Species,
    SampleID
  ) %>%
  dplyr::pull(SampleID)

core100_heatmap_plot_df <- core100_heatmap_long %>%
  dplyr::left_join(
    heatmap_meta,
    by = "SampleID"
  ) %>%
  dplyr::mutate(
    Taxon_Label = factor(
      Taxon_Label,
      levels = rev(taxon_order_df$Taxon_Label)
    ),
    SampleID = factor(
      SampleID,
      levels = sample_order
    ),
    
    # convert to percent for plotting
    RelativePercent = RelativeAbundance * 100,
    
    # log-transform for visualization
    Log10Percent = log10(RelativePercent + 0.01)
  )



core100_heatmap_plot <- ggplot(
  core100_heatmap_plot_df,
  aes(
    x = SampleID,
    y = Taxon_Label,
    fill = Log10Percent
  )
) +
  geom_tile(color = "white", linewidth = 0.2) +
  scale_x_discrete(labels = sample_label_vector) +
  scale_fill_viridis_c(
    option = "viridis",
    name = "log10(% abundance + 0.01)"
  ) +
  labs(
    title = "CCA-enriched core100 bacterial taxa across CCA and seawater samples",
    subtitle = "40 genus/family-level labels derived from 47 enriched core100 EMU taxa",
    x = "Sample",
    y = "Bacterial taxon"
  ) +
  theme_classic() +
  theme(
    axis.text.x = element_text(
      angle = 60,
      hjust = 1,
      vjust = 1,
      size = 9
    ),
    axis.text.y = element_text(size = 9),
    plot.title = element_text(face = "bold")
  )

core100_heatmap_plot


#next want to make mean relative abundance on CCA and CCA:water enrichment
# ------------------------------------------------------------
# Rank CCA-enriched core100 bacterial labels
# ------------------------------------------------------------

core100_rank_df <- core100_heatmap_long %>%
  dplyr::left_join(
    heatmap_meta %>%
      dplyr::select(
        SampleID,
        Environment
      ),
    by = "SampleID"
  ) %>%
  dplyr::group_by(
    Taxon_Label,
    Environment
  ) %>%
  dplyr::summarise(
    MeanRelativeAbundance = mean(RelativeAbundance),
    .groups = "drop"
  ) %>%
  tidyr::pivot_wider(
    names_from = Environment,
    values_from = MeanRelativeAbundance,
    values_fill = 0
  ) %>%
  dplyr::mutate(
    Mean_CCA_Percent = CCA * 100,
    Mean_Water_Percent = Water * 100
  )

# Very small pseudocount used ONLY for the enrichment ratio
pseudocount <- 0.000001

core100_rank_df <- core100_rank_df %>%
  dplyr::mutate(
    Fold_Enrichment = 
      (CCA + pseudocount) /
      (Water + pseudocount),
    
    Log2_Fold_Enrichment =
      log2(Fold_Enrichment)
  )

core100_rank_df <- core100_rank_df %>%
  dplyr::arrange(
    dplyr::desc(Mean_CCA_Percent)
  )

core100_rank_df


core100_rank_df %>%
  dplyr::select(
    Taxon_Label,
    Mean_CCA_Percent,
    Mean_Water_Percent,
    Fold_Enrichment,
    Log2_Fold_Enrichment
  ) %>%
  dplyr::mutate(
    dplyr::across(
      c(
        Mean_CCA_Percent,
        Mean_Water_Percent,
        Fold_Enrichment,
        Log2_Fold_Enrichment
      ),
      ~ round(.x, 3)
    )
  )


#abundance vs enrichment 
candidate_plot <- ggplot(
  core100_rank_df,
  aes(
    x = Mean_CCA_Percent,
    y = Log2_Fold_Enrichment
  )
) +
  
  geom_point(
    size = 3,
    alpha = 0.8
  ) +
  
  ggrepel::geom_text_repel(
    aes(label = Taxon_Label),
    size = 3,
    max.overlaps = 20
  ) +
  
  labs(
    title = "Candidate CCA-associated core bacteria",
    subtitle = "Core100 taxa ranked by abundance on CCA and enrichment relative to seawater",
    x = "Mean relative abundance on CCA (%)",
    y = "log2 CCA:water enrichment"
  ) +
  
  theme_classic()

candidate_plot


#simpler graph
# ------------------------------------------------------------
# Plot CCA abundance vs water abundance
# ------------------------------------------------------------

candidate_plot_abundance <- ggplot(
  core100_rank_df,
  aes(
    x = Mean_CCA_Percent,
    y = Mean_Water_Percent
  )
) +
  
  geom_abline(
    slope = 1,
    intercept = 0,
    linetype = "dashed",
    linewidth = 0.7
  ) +
  
  geom_point(
    size = 3,
    alpha = 0.8
  ) +
  
  ggrepel::geom_text_repel(
    aes(label = Taxon_Label),
    size = 3,
    max.overlaps = 20
  ) +
  
  labs(
    title = "CCA-enriched core bacterial taxa",
    subtitle = "Core100 taxa compared by mean relative abundance in CCA and seawater",
    x = "Mean relative abundance in CCA (%)",
    y = "Mean relative abundance in seawater (%)"
  ) +
  
  theme_classic()

candidate_plot_abundance
#most should fall below the 1:1 dashed line 

# ------------------------------------------------------------
# Core100 relative-abundance stacked bar plot
# ------------------------------------------------------------

ps_core100_rel <- phyloseq::transform_sample_counts(
  ps_core100_turnlund,
  function(x) x / sum(x)
)

ps_core100_family <- phyloseq::tax_glom(
  ps_core100_rel,
  taxrank = "family",
  NArm = FALSE
)

core100_family_long <- phyloseq::psmelt(
  ps_core100_family
)

names(core100_family_long)


core100_family_long <- core100_family_long %>%
  dplyr::mutate(
    family = dplyr::case_when(
      is.na(family) ~ "Unclassified",
      family == "" ~ "Unclassified",
      family == "NF" ~ "Unclassified",
      TRUE ~ family
    )
  )


top30_core100_families <- core100_family_long %>%
  dplyr::group_by(family) %>%
  dplyr::summarise(
    MeanAbundance = mean(Abundance),
    .groups = "drop"
  ) %>%
  dplyr::arrange(
    dplyr::desc(MeanAbundance)
  ) %>%
  dplyr::slice_head(n = 30) %>%
  dplyr::pull(family)


core100_family_plot_df <- core100_family_long %>%
  dplyr::mutate(
    Family_Plot = dplyr::if_else(
      family %in% top30_core100_families,
      family,
      "Other"
    )
  ) %>%
  dplyr::group_by(
    Sample,
    Family_Plot
  ) %>%
  dplyr::summarise(
    Abundance = sum(Abundance),
    .groups = "drop"
  )

core100_family_plot_df %>%
  dplyr::group_by(Sample) %>%
  dplyr::summarise(
    Total = sum(Abundance)
  )


#relative abundance stacked plot
# ------------------------------------------------------------
# Plot core100 relative abundance by family
# ------------------------------------------------------------

# Attach sample metadata
core100_sample_meta <- nmds_meta_turnlund %>%
  dplyr::select(
    SampleID,
    CCA_Morphology,
    Host_Genus,
    Host_Species
  ) %>%
  dplyr::mutate(
    Host_Label = dplyr::case_when(
      !is.na(Host_Species) &
        Host_Species != "" &
        Host_Species != "NF" ~ Host_Species,
      
      !is.na(Host_Genus) &
        Host_Genus != "" &
        Host_Genus != "NF" ~ Host_Genus,
      
      TRUE ~ "Unidentified"
    )
  )


core100_family_plot_df <- core100_family_plot_df %>%
  dplyr::left_join(
    core100_sample_meta,
    by = c("Sample" = "SampleID")
  )

core100_sample_order <- core100_sample_meta %>%
  dplyr::arrange(
    Host_Label,
    SampleID
  ) %>%
  dplyr::pull(SampleID)

core100_family_plot_df$Sample <- factor(
  core100_family_plot_df$Sample,
  levels = core100_sample_order
)

core100_barplot <- ggplot(
  core100_family_plot_df,
  aes(
    x = Sample,
    y = Abundance * 100,
    fill = Family_Plot
  )
) +
  
  geom_col(
    width = 0.9
  ) +
  
  facet_grid(
    ~ Host_Label,
    scales = "free_x",
    space = "free_x"
  ) +
  
  scale_fill_viridis_d(
    option = "viridis",
    begin = 0,
    end = 1
  ) +
  
  labs(
    title = "Core100 CCA bacterial community composition",
    subtitle = "Relative abundance of the 30 most abundant bacterial families",
    x = "CCA sample",
    y = "Relative abundance (%)",
    fill = "Bacterial family"
  ) +
  
  theme_classic() +
  
  theme(
    axis.text.x = element_text(
      angle = 60,
      hjust = 1,
      vjust = 1,
      size = 8
    ),
    
    strip.background = element_rect(
      fill = "grey95",
      color = NA
    ),
    
    strip.text = element_text(
      face = "bold"
    ),
    
    legend.position = "right",
    
    plot.title = element_text(
      face = "bold"
    )
  )

core100_barplot


#create the same plot with a different color palette


# Convert core100 community to relative abundance
ps_core100_rel <- phyloseq::transform_sample_counts(
  ps_core100_turnlund,
  function(x) x / sum(x)
)

# Collapse taxa to bacterial family
ps_core100_family <- phyloseq::tax_glom(
  ps_core100_rel,
  taxrank = "family",
  NArm = FALSE
)

# Convert to long format
core100_family_long <- phyloseq::psmelt(
  ps_core100_family
)

# Clean missing family labels
core100_family_long <- core100_family_long %>%
  dplyr::mutate(
    family = dplyr::case_when(
      is.na(family) ~ "Unclassified",
      family == "" ~ "Unclassified",
      family == "NF" ~ "Unclassified",
      TRUE ~ family
    )
  )

# Identify the 20 most abundant bacterial families
top20_core100_families <- core100_family_long %>%
  dplyr::group_by(family) %>%
  dplyr::summarise(
    MeanAbundance = mean(Abundance),
    .groups = "drop"
  ) %>%
  dplyr::arrange(
    dplyr::desc(MeanAbundance)
  ) %>%
  dplyr::slice_head(n = 20) %>%
  dplyr::pull(family)

# Collapse all remaining families into "Other"
core100_family_plot_df <- core100_family_long %>%
  dplyr::mutate(
    Family_Plot = dplyr::if_else(
      family %in% top20_core100_families,
      family,
      "Other"
    )
  ) %>%
  dplyr::group_by(
    Sample,
    Family_Plot
  ) %>%
  dplyr::summarise(
    Abundance = sum(Abundance),
    .groups = "drop"
  )

# ------------------------------------------------------------
# Attach CCA metadata
# ------------------------------------------------------------

core100_sample_meta <- nmds_meta_turnlund %>%
  dplyr::select(
    SampleID,
    CCA_Morphology,
    Host_Genus,
    Host_Species
  ) %>%
  dplyr::mutate(
    Host_Label = dplyr::case_when(
      
      !is.na(Host_Species) &
        Host_Species != "" &
        Host_Species != "NF" ~ Host_Species,
      
      !is.na(Host_Genus) &
        Host_Genus != "" &
        Host_Genus != "NF" ~ Host_Genus,
      
      TRUE ~ "Unidentified"
    )
  )

# Join metadata to plotting dataframe
core100_family_plot_df <- core100_family_plot_df %>%
  dplyr::left_join(
    core100_sample_meta,
    by = c("Sample" = "SampleID")
  )

# ------------------------------------------------------------
# Order samples by host identity
# ------------------------------------------------------------

core100_sample_order <- core100_sample_meta %>%
  dplyr::arrange(
    Host_Label,
    SampleID
  ) %>%
  dplyr::pull(SampleID)

core100_family_plot_df$Sample <- factor(
  core100_family_plot_df$Sample,
  levels = core100_sample_order
)

# ------------------------------------------------------------
# Set family order and Babel categorical colors
# ------------------------------------------------------------
install.packages("biopalette")
library(phyloseq)
library(dplyr)
library(ggplot2)
library(biopalette)
family_levels_core100 <- c(
  top20_core100_families,
  "Other"
)

core100_family_plot_df$Family_Plot <- factor(
  core100_family_plot_df$Family_Plot,
  levels = family_levels_core100
)

babel_colors_core100 <- biopalette::get_palette(
  "babel",
  n = 21
)

names(babel_colors_core100) <- family_levels_core100

# ------------------------------------------------------------
# Check that each sample still sums to 100%
# ------------------------------------------------------------

core100_family_plot_df %>%
  dplyr::group_by(Sample) %>%
  dplyr::summarise(
    Total = sum(Abundance)
  )

# ------------------------------------------------------------
# Plot
# ------------------------------------------------------------

core100_barplot <- ggplot(
  core100_family_plot_df,
  aes(
    x = Sample,
    y = Abundance * 100,
    fill = Family_Plot
  )
) +
  
  geom_col(
    width = 0.9
  ) +
  
  facet_grid(
    ~ Host_Label,
    scales = "free_x",
    space = "free_x"
  ) +
  
  scale_fill_manual(
    values = babel_colors_core100,
    drop = FALSE
  ) +
  
  labs(
    title = "Core100 CCA bacterial community composition",
    subtitle = "Relative abundance of the 20 most abundant bacterial families",
    x = "CCA sample",
    y = "Relative abundance (%)",
    fill = "Bacterial family"
  ) +
  
  theme_classic() +
  
  theme(
    axis.text.x = element_text(
      angle = 60,
      hjust = 1,
      vjust = 1,
      size = 8
    ),
    
    strip.background = element_rect(
      fill = "grey95",
      color = NA
    ),
    
    strip.text = element_text(
      face = "bold"
    ),
    
    legend.position = "right",
    
    plot.title = element_text(
      face = "bold"
    )
  )

core100_barplot
