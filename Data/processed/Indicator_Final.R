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



# ============================================================
# CCA INDICATOR TAXA ABUNDANCE FIGURE
# Grouped by morphology and CCA host identity
# ============================================================

library(dplyr)
library(tidyr)
library(ggplot2)
library(stringr)
library(MoMAColors)


# ------------------------------------------------------------
# 1. Make sure metadata is a plain dataframe
# ------------------------------------------------------------

meta_df <- data.frame(
  meta_df,
  check.names = FALSE
)

meta_df$SampleID <- as.character(
  meta_df$SampleID
)


# ------------------------------------------------------------
# 2. Select prioritized CCA indicator taxa
# ------------------------------------------------------------

indicator_ids <- CCA_indicators_priority$Taxon


# ------------------------------------------------------------
# 3. Create readable bacterial taxon labels
# ------------------------------------------------------------

indicator_tax_labels <- CCA_indicators_priority %>%
  dplyr::mutate(
    
    TaxonLabel = dplyr::case_when(
      
      # Genus + species
      !is.na(genus) &
        genus != "" &
        genus != "NF" &
        !is.na(species) &
        species != "" &
        species != "NF" ~
        paste(genus, species),
      
      # Genus only
      !is.na(genus) &
        genus != "" &
        genus != "NF" ~
        genus,
      
      # Family only
      !is.na(family) &
        family != "" &
        family != "NF" ~
        paste0(family, " (family)"),
      
      # Fall back to EMU taxon ID
      TRUE ~ Taxon
    )
  ) %>%
  
  dplyr::select(
    Taxon,
    TaxonLabel
  )


# ------------------------------------------------------------
# 4. Pull indicator taxa from relative abundance matrix
# ------------------------------------------------------------

indicator_mat <- beta_rel[
  ,
  indicator_ids,
  drop = FALSE
]

indicator_df <- as.data.frame(
  indicator_mat
)

indicator_df$SampleID <- rownames(
  indicator_df
)


# ------------------------------------------------------------
# 5. Create sample, morphology, and host labels
# ------------------------------------------------------------

meta_plot <- meta_df %>%
  dplyr::mutate(
    
    # barcode01 -> Sample 1
    SampleLabel = paste0(
      "Sample ",
      as.integer(
        stringr::str_extract(
          SampleID,
          "[0-9]+"
        )
      )
    ),
    
    # Morphology
    Morphology_Label = dplyr::case_when(
      
      Environment == "Water" ~
        "Water",
      
      tolower(as.character(CCA_Morphology)) ==
        "encrusting" ~
        "Encrusting",
      
      tolower(as.character(CCA_Morphology)) ==
        "branching" ~
        "Branching",
      
      TRUE ~
        "Unknown"
    ),
    
    # CCA host identity
    Host_Label = dplyr::case_when(
      
      Environment == "Water" ~
        "r",
      
      # Genus + species
      !is.na(Host_Genus) &
        Host_Genus != "" &
        Host_Genus != "NF" &
        !is.na(Host_Species) &
        Host_Species != "" &
        Host_Species != "NF" ~
        paste(
          Host_Genus,
          Host_Species
        ),
      
      # Genus only
      !is.na(Host_Genus) &
        Host_Genus != "" &
        Host_Genus != "NF" ~
        paste0(
          Host_Genus,
          " sp."
        ),
      
      # No host identification
      TRUE ~
        "Unidentified"
    )
  )


# ------------------------------------------------------------
# 6. Convert bacterial abundances to long format
# and attach sample metadata
# ------------------------------------------------------------

indicator_long <- indicator_df %>%
  
  tidyr::pivot_longer(
    cols = -SampleID,
    names_to = "Taxon",
    values_to = "RelativeAbundance"
  ) %>%
  
  dplyr::left_join(
    indicator_tax_labels,
    by = "Taxon"
  ) %>%
  
  dplyr::left_join(
    meta_plot %>%
      dplyr::select(
        SampleID,
        SampleLabel,
        Environment,
        Morphology_Label,
        Host_Label
      ),
    by = "SampleID"
  )


# ------------------------------------------------------------
# 7. Set morphology order
# ------------------------------------------------------------

indicator_long$Morphology_Label <- factor(
  indicator_long$Morphology_Label,
  levels = c(
    "Encrusting",
    "Branching",
    "Water",
    "Unknown"
  )
)


# ------------------------------------------------------------
# 8. Order CCA host groups
# ------------------------------------------------------------

host_order <- meta_plot %>%
  
  dplyr::distinct(
    Host_Label,
    Morphology_Label
  ) %>%
  
  dplyr::mutate(
    Morphology_Label = factor(
      Morphology_Label,
      levels = c(
        "Encrusting",
        "Branching",
        "Water",
        "Unknown"
      )
    )
  ) %>%
  
  dplyr::arrange(
    Morphology_Label,
    Host_Label
  ) %>%
  
  dplyr::pull(
    Host_Label
  ) %>%
  
  unique()


indicator_long$Host_Label <- factor(
  indicator_long$Host_Label,
  levels = host_order
)


# ------------------------------------------------------------
# 9. Order individual samples
# ------------------------------------------------------------

sample_order <- meta_plot %>%
  
  dplyr::mutate(
    Morphology_Label = factor(
      Morphology_Label,
      levels = c(
        "Encrusting",
        "Branching",
        "Water",
        "Unknown"
      )
    )
  ) %>%
  
  dplyr::arrange(
    Morphology_Label,
    Host_Label,
    SampleID
  ) %>%
  
  dplyr::pull(
    SampleLabel
  )


indicator_long$SampleLabel <- factor(
  indicator_long$SampleLabel,
  levels = sample_order
)


# ------------------------------------------------------------
# 10. Original MoMAColors "Warhol" qualitative palette
# ------------------------------------------------------------

indicator_levels <- unique(
  indicator_tax_labels$TaxonLabel
)

indicator_colors <- MoMAColors::moma.colors(
  "Warhol",
  n = length(indicator_levels),
  type = "discrete"
)

names(indicator_colors) <- indicator_levels

# ------------------------------------------------------------
# 11. ONLY palette modification:
# Swap Portibacter and Litoreibacter janthinus
# ------------------------------------------------------------

portibacter_color <-
  indicator_colors["Portibacter"]

litoreibacter_janthinus_color <-
  indicator_colors["Litoreibacter janthinus"]


indicator_colors["Portibacter"] <-
  litoreibacter_janthinus_color

indicator_colors["Litoreibacter janthinus"] <-
  portibacter_color

# ------------------------------------------------------------
# 12. Make indicator abundance plot
# ------------------------------------------------------------
install.packages("ggh4x")
library(ggh4x)


CCA_indicator_plot <- ggplot2::ggplot(
  indicator_long,
  ggplot2::aes(
    x = SampleLabel,
    y = RelativeAbundance,
    fill = TaxonLabel
  )
) +
  
  ggplot2::geom_col(
    width = 0.85
  ) +
  
  # Morphology + host identity
  ggh4x::facet_nested(
    cols = ggplot2::vars(
      Morphology_Label,
      Host_Label
    ),
    scales = "free_x",
    space = "free_x",
    
    # Draw nesting lines so species boundaries are obvious
    nest_line = ggplot2::element_line(
      colour = "grey30",
      linewidth = 0.7
    ),
    
    # Also draw the hierarchy line when a morphology
    # contains only one host category
    solo_line = TRUE
  ) +
  
  # Warhol categorical palette
  ggplot2::scale_fill_manual(
    values = indicator_colors
  ) +
  
  # Relative abundance shown as percent
  ggplot2::scale_y_continuous(
    labels = scales::percent,
    expand = ggplot2::expansion(
      mult = c(
        0,
        0.05
      )
    )
  ) +
  
  ggplot2::labs(
    title = "Relative abundance of CCA indicator taxa",
    subtitle = "Strong CCA-associated bacterial taxa grouped by host identity and morphology",
    x = "Sample",
    y = "Relative abundance",
    fill = "CCA indicator taxon"
  ) +
  
  ggplot2::theme_classic() +
  
  ggplot2::theme(
    
    # Individual sample labels
    axis.text.x = ggplot2::element_text(
      angle = 55,
      hjust = 1,
      vjust = 1,
      size = 8
    ),
    
    # Facet strip formatting
    strip.background = ggplot2::element_rect(
      fill = "grey95",
      color = "black",
      linewidth = 0.7
    ),
    
    strip.text.x = ggplot2::element_text(
      size = 9,
      margin = ggplot2::margin(
        t = 5,
        r = 6,
        b = 5,
        l = 6
      )
    ),
    
    # Borders around host groups
    panel.border = ggplot2::element_rect(
      color = "grey30",
      fill = NA,
      linewidth = 0.6
    ),
    
    # Space between host-species panels
    panel.spacing.x = grid::unit(
      0.18,
      "cm"
    ),
    
    legend.position = "right",
    
    plot.title = ggplot2::element_text(
      face = "bold"
    )
  )


# Display figure
CCA_indicator_plot



###nMDS which uses the full bacterial community to put the points in ordination space bu the vectors are based on the CCA indicators 
### vectors are significantly associated with the main community gradients
# ============================================================
# 12. NMDS + VECTORS FOR PRIORITIZED CCA INDICATOR TAXA
# ============================================================

library(vegan)
library(dplyr)
library(ggplot2)
library(ggrepel)

# ------------------------------------------------------------
# Run Bray-Curtis NMDS on the FULL bacterial community
# ------------------------------------------------------------

set.seed(123)

nmds_fit <- vegan::metaMDS(
  beta_rel,
  distance = "bray",
  k = 2,
  trymax = 200,
  autotransform = FALSE,
  trace = FALSE
)

# Check stress
nmds_fit$stress
#stress = 0.072

# ------------------------------------------------------------
# Extract sample coordinates
# ------------------------------------------------------------

nmds_scores <- as.data.frame(
  vegan::scores(
    nmds_fit,
    display = "sites"
  )
)

nmds_scores$SampleID <- rownames(
  nmds_scores
)

# Attach metadata
# ------------------------------------------------------------
# 2. Extract sample coordinates + plotting metadata
# ------------------------------------------------------------

nmds_scores <- as.data.frame(
  vegan::scores(
    nmds_fit,
    display = "sites"
  )
)

nmds_scores$SampleID <- rownames(
  nmds_scores
)

nmds_scores <- nmds_scores %>%
  dplyr::left_join(
    meta_plot %>%
      dplyr::select(
        SampleID,
        Environment,
        CCA_Morphology
      ),
    by = "SampleID"
  ) %>%
  
  dplyr::mutate(
    
    # barcode01 -> Sample 1
    SampleLabel = paste0(
      "Sample ",
      as.integer(
        stringr::str_extract(
          SampleID,
          "[0-9]+"
        )
      )
    ),
    
    # Create plotting groups
    PlotGroup = dplyr::case_when(
      Environment == "Water" ~ "Water",
      CCA_Morphology == "encrusting" ~ "Encrusting",
      CCA_Morphology == "branching" ~ "Branching",
      TRUE ~ "Unknown"
    )
  )

# ------------------------------------------------------------
# Extract vector coordinates, R2, and p-values
# ------------------------------------------------------------

indicator_vectors <- as.data.frame(
  vegan::scores(
    indicator_envfit,
    display = "vectors"
  )
)

indicator_vectors$Taxon <- rownames(
  indicator_vectors
)

indicator_vectors$R2 <-
  indicator_envfit$vectors$r

indicator_vectors$p_value <-
  indicator_envfit$vectors$pvals


# ------------------------------------------------------------
# Attach readable bacterial names
# ------------------------------------------------------------

indicator_vectors <- indicator_vectors %>%
  dplyr::left_join(
    indicator_tax_labels,
    by = "Taxon"
  ) %>%
  dplyr::arrange(
    dplyr::desc(R2)
  )

# View results
indicator_vectors %>%
  dplyr::select(
    Taxon,
    TaxonLabel,
    R2,
    p_value
  )


# ------------------------------------------------------------
# Keep only vectors significantly associated
# with the NMDS configuration
# ------------------------------------------------------------

indicator_vectors_sig <- indicator_vectors %>%
  dplyr::filter(
    p_value < 0.05
  ) %>%
  dplyr::arrange(
    dplyr::desc(R2)
  )

indicator_vectors_sig %>%
  dplyr::select(
    Taxon,
    TaxonLabel,
    R2,
    p_value
  )


# ------------------------------------------------------------
# Scale vectors so they fit cleanly on the NMDS
# ------------------------------------------------------------

nmds_x_range <- diff(
  range(nmds_scores$NMDS1)
)

nmds_y_range <- diff(
  range(nmds_scores$NMDS2)
)

vector_x_max <- max(
  abs(indicator_vectors_sig$NMDS1)
)

vector_y_max <- max(
  abs(indicator_vectors_sig$NMDS2)
)

arrow_multiplier <- min(
  0.40 * nmds_x_range / vector_x_max,
  0.40 * nmds_y_range / vector_y_max
)

indicator_vectors_sig <- indicator_vectors_sig %>%
  dplyr::mutate(
    NMDS1_plot =
      NMDS1 * arrow_multiplier,
    
    NMDS2_plot =
      NMDS2 * arrow_multiplier
  )


# ------------------------------------------------------------
# 8. Create Encrusting / Branching / Water colors
# ------------------------------------------------------------

group_colors <- viridisLite::viridis(
  n = 3,
  option = "viridis",
  begin = 0.10,
  end = 0.90
)

names(group_colors) <- c(
  "Encrusting",
  "Branching",
  "Water"
)


# ------------------------------------------------------------
# 9. Plot NMDS with significant indicator vectors
# ------------------------------------------------------------

CCA_nmds_vectors_plot <- ggplot2::ggplot(
  nmds_scores,
  ggplot2::aes(
    x = NMDS1,
    y = NMDS2,
    color = PlotGroup
  )
) +
  
  ggplot2::geom_point(
    size = 4,
    alpha = 0.9
  ) +
  
  ggrepel::geom_text_repel(
    ggplot2::aes(
      label = SampleLabel
    ),
    size = 3,
    show.legend = FALSE
  ) +
  
  ggplot2::geom_segment(
    data = indicator_vectors_sig,
    ggplot2::aes(
      x = 0,
      y = 0,
      xend = NMDS1_plot,
      yend = NMDS2_plot
    ),
    inherit.aes = FALSE,
    color = "black",
    linewidth = 0.7,
    arrow = grid::arrow(
      length = grid::unit(
        0.18,
        "cm"
      )
    )
  ) +
  
  ggrepel::geom_text_repel(
    data = indicator_vectors_sig,
    ggplot2::aes(
      x = NMDS1_plot,
      y = NMDS2_plot,
      label = TaxonLabel
    ),
    inherit.aes = FALSE,
    color = "black",
    size = 3.2,
    fontface = "italic"
  ) +
  
  ggplot2::scale_color_manual(
    values = group_colors
  ) +
  
  ggplot2::labs(
    title = "CCA and seawater bacterial community structure",
    subtitle = "Bray-Curtis NMDS with significant CCA indicator taxa vectors",
    x = "NMDS1",
    y = "NMDS2",
    color = "Sample type"
  ) +
  
  ggplot2::theme_classic()

CCA_nmds_vectors_plot
