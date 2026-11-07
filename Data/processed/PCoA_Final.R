# ============================================================
# CCA VS SEAWATER BRAY-CURTIS PCoA
# Final presentation version
# ============================================================

library(phyloseq)
library(vegan)
library(ape)
library(dplyr)
library(ggplot2)
library(ggrepel)
library(stringr)
library(MoMAColors)


# ------------------------------------------------------------
# 1. Load cleaned biological-sample phyloseq object
# ------------------------------------------------------------

ps_analysis <- readRDS(
  "../annabelle_16S_analysis_no_controls_phyloseq.rds"
)


# ------------------------------------------------------------
# 2. Extract microbial abundance matrix
# ------------------------------------------------------------

beta_mat <- as(
  phyloseq::otu_table(ps_analysis),
  "matrix"
)

# vegan requires samples as rows and taxa as columns
if (phyloseq::taxa_are_rows(ps_analysis)) {
  beta_mat <- t(beta_mat)
}


# ------------------------------------------------------------
# 3. Convert to relative abundance
# ------------------------------------------------------------

beta_rel <- sweep(
  beta_mat,
  1,
  rowSums(beta_mat),
  "/"
)


# ------------------------------------------------------------
# 4. Extract metadata
# ------------------------------------------------------------

beta_meta <- data.frame(
  phyloseq::sample_data(ps_analysis),
  check.names = FALSE
)

# Match metadata order to abundance matrix
beta_meta <- beta_meta[
  rownames(beta_rel),
  ,
  drop = FALSE
]


# ------------------------------------------------------------
# 5. Calculate Bray-Curtis distances
# ------------------------------------------------------------

bray_env <- vegan::vegdist(
  beta_rel,
  method = "bray"
)

# ------------------------------------------------------------
# PERMANOVA: CCA vs Water
# ------------------------------------------------------------

set.seed(123)

permanova_env <- vegan::adonis2(
  bray_env ~ Environment,
  data = beta_meta,
  permutations = 999
)

permanova_env

# Extract p-value for plotting
permanova_p <- permanova_env$`Pr(>F)`[1]

# Make presentation-friendly label
p_label <- ifelse(
  permanova_p < 0.001,
  "PERMANOVA p < 0.001",
  paste0(
    "PERMANOVA p = ",
    format(
      permanova_p,
      digits = 3
    )
  )
)


# ------------------------------------------------------------
# 6. Run PCoA
# ------------------------------------------------------------

pcoa_env <- ape::pcoa(
  bray_env
)


# ------------------------------------------------------------
# 7. Extract PCoA coordinates
# ------------------------------------------------------------

pcoa_env_scores <- data.frame(
  SampleID = rownames(
    pcoa_env$vectors
  ),
  
  PCoA1 = pcoa_env$vectors[, 1],
  
  PCoA2 = pcoa_env$vectors[, 2]
)


# ------------------------------------------------------------
# 8. Add metadata
# ------------------------------------------------------------

meta_env_join <- beta_meta

meta_env_join$SampleID <- rownames(
  meta_env_join
)

pcoa_env_scores <- dplyr::left_join(
  pcoa_env_scores,
  meta_env_join,
  by = "SampleID"
)


# ------------------------------------------------------------
# 9. Create presentation sample labels
# barcode01 -> Sample 1
# ------------------------------------------------------------

pcoa_env_scores <- pcoa_env_scores %>%
  dplyr::mutate(
    
    SampleLabel = paste0(
      "Sample ",
      as.integer(
        stringr::str_extract(
          SampleID,
          "[0-9]+"
        )
      )
    )
  )


# ------------------------------------------------------------
# 10. Calculate percent variation explained
# ------------------------------------------------------------

var_exp_env <-
  pcoa_env$values$Relative_eig * 100

var_exp_env[1:2]


# ------------------------------------------------------------
# 11. MoMAColors Warhol palette
# ------------------------------------------------------------

env_colors <- MoMAColors::moma.colors(
  "Warhol",
  n = 2,
  type = "discrete"
)

names(env_colors) <- c(
  "CCA",
  "Water"
)


# ------------------------------------------------------------
# 12. Final PCoA figure
# ------------------------------------------------------------

CCA_water_PCoA <- ggplot2::ggplot(
  pcoa_env_scores,
  ggplot2::aes(
    x = PCoA1,
    y = PCoA2,
    color = Environment
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
    max.overlaps = Inf,
    show.legend = FALSE
  ) +
  
  ggplot2::annotate(
    "text",
    x = Inf,
    y = Inf,
    label = p_label,
    hjust = 1.1,
    vjust = 1.5,
    size = 4
  ) +
  
  ggplot2::scale_color_manual(
    values = env_colors
  ) +
  
  ggplot2::labs(
    
    title =
      "Bray-Curtis PCoA of CCA and Seawater Bacterial Communities",
    
    x = paste0(
      "PCoA1 (",
      round(
        var_exp_env[1],
        1
      ),
      "%)"
    ),
    
    y = paste0(
      "PCoA2 (",
      round(
        var_exp_env[2],
        1
      ),
      "%)"
    ),
    
    color = "Sample type"
  ) +
  
  ggplot2::theme_classic() +
  
  ggplot2::theme(
    plot.title = ggplot2::element_text(
      face = "bold"
    ),
    legend.position = "right"
  )


# Display figure
CCA_water_PCoA
