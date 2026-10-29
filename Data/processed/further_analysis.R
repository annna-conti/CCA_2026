# ============================================================
# MAKING THE RAREFACTION CURVES
# ============================================================
qc_before_after


ps_analysis <- phyloseq::prune_samples(
  phyloseq::sample_data(ps_decontam)$ControlType == "none",
  ps_decontam
)

phyloseq::nsamples(ps_analysis)
phyloseq::sample_names(ps_analysis)

saveRDS(
  ps_analysis,
  "../annabelle_16S_analysis_no_controls_phyloseq.rds"
)

library(phyloseq)
#needed to add some stuff to metadata
sample_data(ps_analysis)$Environment <- ifelse(
  sample_data(ps_analysis)$SampleType == "CCA",
  "CCA",
  "Water"
)
phyloseq::sample_data(ps_analysis)$Concentrated <- ifelse(
  phyloseq::sample_data(ps_analysis)$SampleOG == "6B" |
    phyloseq::sample_data(ps_analysis)$Environment == "Water",
  "Yes",
  "No"
)
phyloseq::sample_data(ps_analysis)$Concentrated <- factor(
  phyloseq::sample_data(ps_analysis)$Concentrated,
  levels = c("No", "Yes")
)
phyloseq::sample_data(ps_analysis)[
  ,
  c(
    "SampleOG",
    "SampleType",
    "Environment",
    "Concentrated"
  )
]

saveRDS(
  ps_analysis,
  "../annabelle_16S_analysis_no_controls_phyloseq.rds"
)

####need to make a copy of this for rarefaction only because it needs integer counts and EMU gives fractional estimates, so will make a rounded copy only for visualizing rarefaction
#####
ps_rarefy <- ps_analysis

otu_round <- as(
  phyloseq::otu_table(ps_rarefy),
  "matrix"
)

# Round EMU estimated counts to nearest whole number
otu_round <- round(otu_round)

# Put rounded matrix back into phyloseq
phyloseq::otu_table(ps_rarefy) <- phyloseq::otu_table(
  otu_round,
  taxa_are_rows = phyloseq::taxa_are_rows(ps_rarefy)
)

# Remove any taxa that became zero everywhere after rounding
ps_rarefy <- phyloseq::prune_taxa(
  phyloseq::taxa_sums(ps_rarefy) > 0,
  ps_rarefy
)

ps_rarefy
phyloseq::nsamples(ps_rarefy)
phyloseq::ntaxa(ps_rarefy)

#extract matrix in vegan format
rare_mat <- as(
  phyloseq::otu_table(ps_rarefy),
  "matrix"
)

if (phyloseq::taxa_are_rows(ps_rarefy)) {
  rare_mat <- t(rare_mat)
}

dim(rare_mat)

head(rownames(rare_mat))

#check sample depth and observed richness 
rare_qc <- data.frame(
  SampleID = rownames(rare_mat),
  EstimatedReads = rowSums(rare_mat),
  ObservedTaxa = vegan::specnumber(rare_mat)
)

rare_qc

#make rarefaction curves
library(vegan)

rarecurve(
  rare_mat,
  step = 1000,
  label = TRUE,
  xlab = "Estimated reads",
  ylab = "Observed taxa"
)

rarecurve(
  rare_mat,
  step = 100,
  label = TRUE,
  xlab = "Estimated reads",
  ylab = "Observed taxa"
)
#step controls how densely vegan calculates along each curve - will use step 100 because it gives smoother curves 
#sequencing depth appears sufficient to characterize the observed taxonomic richness of all 17 samples
#a sample can have fewer reads but still adequately characterize less-rich communities if the curve plateaus
#some samples plateau at very different observed richness values, which I will quantify with alpha diversity

plateau_check <- data.frame(
  SampleID = rownames(rare_mat),
  TotalReads = rowSums(rare_mat),
  ObservedTaxa = vegan::specnumber(rare_mat)
)
#creates a summary table --> totalreads is rounded EMU estimated sequencing depth, observedtaxa is the number of taxa with nonzero abundance
#vegan::specnumber() is counting how many taxa are present 

plateau_check$Depth90 <- floor(
  plateau_check$TotalReads * 0.90
)
#calculates 90% of each sample's sequencing depth --> to compare richness at 90% of the reads with richness at 100% of the reads

plateau_check$Richness90 <- sapply(
  seq_len(nrow(rare_mat)),
  function(i) {
    vegan::rarefy(
      rare_mat[i, , drop = FALSE],
      sample = plateau_check$Depth90[i]
    )
  }
)
#"if i had sequenced only 90% as deeply, how many taxa would I expect to observe 

plateau_check$GainLast10Percent <-
  plateau_check$ObservedTaxa -
  plateau_check$Richness90
#calculates how many additional taxa were gained during the final 10% of sequencing 

plateau_check$PercentGainLast10Percent <-
  100 *
  plateau_check$GainLast10Percent /
  plateau_check$ObservedTaxa
#expresses that gain as a percentage os the sample's final richness 

plateau_check
plateau_check[
  order(
    plateau_check$PercentGainLast10Percent,
    decreasing = TRUE
  ),
]
#plateau check: compares expected richness at 90% sequencing depth with observed richness at full depth
#small percent gain in the last 10% indicates that the rarefaction curve is approaching a plateau and sequencing depth is likely adequate 

####ALPHA DIVERSITY####
# Extract original CLEAN, unrounded EMU abundances
alpha_mat <- as(
  phyloseq::otu_table(ps_analysis),
  "matrix"
)
#extracts the cleaned microbial abundance table from phyloseq object and turns it into normal numeric matrix


if (phyloseq::taxa_are_rows(ps_analysis)) {
  alpha_mat <- t(alpha_mat)
}
#makes sure matrix is in the orientation that vegan expects --> vegan calculates diversity across each row

# Align metadata
alpha_meta <- as.data.frame(
  phyloseq::sample_data(ps_analysis)
)

alpha_meta <- alpha_meta[
  rownames(alpha_mat),
  ,
  drop = FALSE
]

alpha_df <- data.frame(
  SampleID = rownames(alpha_mat),
  
  Observed = vegan::specnumber(alpha_mat),
  
  Shannon = vegan::diversity(
    alpha_mat,
    index = "shannon"
  ),
  
  Simpson = vegan::diversity(
    alpha_mat,
    index = "simpson"
  )
)
#calculating three alpha diversity metrics 
#observed:number of taxa detected in a sample
#Shannon:combines richness and evenness --> increases when a sample has more taxa and when those taxa are more evenly distributed
#--> high Shannon: many taxa + relatively even abundances 
#--> low Shannon: fewer taxa/strong dominance by a few taxa 
#Simpson: rare taxa influence simpsion less strongly than they influence richness 
#--> how dominated is the community by common taxa


alpha_df <- cbind(
  alpha_df,
  alpha_meta[alpha_df$SampleID, , drop = FALSE]
)

alpha_df

library(tidyr)
library(ggplot2)

library(dplyr)

# Diversity metrics only
alpha_metrics <- data.frame(
  SampleID = rownames(alpha_mat),
  Observed = vegan::specnumber(alpha_mat),
  Shannon = vegan::diversity(
    alpha_mat,
    index = "shannon"
  ),
  Simpson = vegan::diversity(
    alpha_mat,
    index = "simpson"
  )
)

# Metadata as a normal dataframe
alpha_meta <- as.data.frame(
  phyloseq::sample_data(ps_analysis)
)

# Make sure SampleID exists as a regular column
alpha_meta$SampleID <- rownames(alpha_meta)

# Join by SampleID
alpha_df <- alpha_metrics %>%
  left_join(
    alpha_meta,
    by = "SampleID"
  )

alpha_CCA <- alpha_df %>%
  filter(Environment == "CCA")

table(alpha_CCA$CCA_Morphology)

#cbind attaches columns side by side and can accidentally create duplicate column names
#left_join... matches rows by sample ID and is safer for combining diversity results with metadata

#CCA morphology comparison
library(dplyr)
library(tidyr)
library(ggplot2)

alpha_CCA <- alpha_df %>%
  filter(Environment == "CCA")

alpha_CCA_long <- alpha_CCA %>%
  select(
    SampleID,
    CCA_Morphology,
    Observed,
    Shannon,
    Simpson
  ) %>%
  pivot_longer(
    cols = c(Observed, Shannon, Simpson),
    names_to = "Metric",
    values_to = "Value"
  )

#plot encrusting vs branching
ggplot(
  alpha_CCA_long,
  aes(
    x = CCA_Morphology,
    y = Value
  )
) +
  geom_boxplot(
    outlier.shape = NA,
    width = 0.5
  ) +
  geom_jitter(
    width = 0.10,
    size = 3,
    alpha = 0.8
  ) +
  facet_wrap(
    ~ Metric,
    scales = "free_y"
  ) +
  labs(
    x = "CCA morphology",
    y = NULL,
    title = "Alpha diversity of CCA-associated bacterial communities"
  ) +
  theme_classic()

library(ggrepel)

ggplot(
  alpha_CCA_long,
  aes(
    x = CCA_Morphology,
    y = Value
  )
) +
  geom_boxplot(
    outlier.shape = NA,
    width = 0.5
  ) +
  geom_jitter(
    width = 0.10,
    size = 3,
    alpha = 0.8
  ) +
  geom_text_repel(
    aes(label = SampleID),
    size = 3,
    max.overlaps = Inf,
    show.legend = FALSE
  ) +
  facet_wrap(
    ~ Metric,
    scales = "free_y"
  ) +
  labs(
    x = "CCA morphology",
    y = NULL,
    title = "Alpha diversity of CCA-associated bacterial communities"
  ) +
  theme_classic()

#CCA vs water 
alpha_env_long <- alpha_df %>%
  select(
    SampleID,
    Environment,
    Observed,
    Shannon,
    Simpson
  ) %>%
  pivot_longer(
    cols = c(Observed, Shannon, Simpson),
    names_to = "Metric",
    values_to = "Value"
  )

ggplot(
  alpha_env_long,
  aes(
    x = Environment,
    y = Value
  )
) +
  geom_boxplot(
    outlier.shape = NA,
    width = 0.5
  ) +
  geom_jitter(
    width = 0.10,
    size = 3,
    alpha = 0.8
  ) +
  facet_wrap(
    ~ Metric,
    scales = "free_y"
  ) +
  labs(
    x = NULL,
    y = NULL,
    title = "Alpha diversity of CCA and water bacterial communities"
  ) +
  theme_classic()

#wilcoxin rank-sum test for each metric to get p values 
alpha_tests <- data.frame(
  Metric = c("Observed", "Shannon", "Simpson"),
  p_value = c(
    wilcox.test(
      Observed ~ CCA_Morphology,
      data = alpha_CCA,
      exact = FALSE
    )$p.value,
    
    wilcox.test(
      Shannon ~ CCA_Morphology,
      data = alpha_CCA,
      exact = FALSE
    )$p.value,
    
    wilcox.test(
      Simpson ~ CCA_Morphology,
      data = alpha_CCA,
      exact = FALSE
    )$p.value
  )
)

alpha_tests
#test has really low power with such an unbalanced deisgn with only two branching samples, so need to interpret these p values cautiously


###BETA DIVERSITY###
library(phyloseq)
library(vegan)
library(ape)
library(dplyr)
library(ggplot2)
library(ggrepel)

beta_mat <- as(
  phyloseq::otu_table(ps_analysis),
  "matrix"
)

# vegan wants samples as rows, taxa as columns
if (phyloseq::taxa_are_rows(ps_analysis)) {
  beta_mat <- t(beta_mat)
}

#extract the cleaned microbial abundance table and orient it as samples x taxa
#unlike rarefaction, can use original unrounded EMU abundance estimates 

#convert to relative abundance to make Bray-Curtis describe differences in community composition instead of differences in total sequencing depth
beta_rel <- sweep(
  beta_mat,
  1,
  rowSums(beta_mat),
  "/"
)

beta_meta <- as.data.frame(
  phyloseq::sample_data(ps_analysis)
)

beta_meta <- beta_meta[
  rownames(beta_rel),
  ,
  drop = FALSE
]

CCA_samples <- rownames(beta_meta)[
  beta_meta$Environment == "CCA"
]

beta_CCA <- beta_rel[
  CCA_samples,
  ,
  drop = FALSE
]

meta_CCA <- beta_meta[
  CCA_samples,
  ,
  drop = FALSE
]

table(meta_CCA$CCA_Morphology)

beta_mat <- as(
  phyloseq::otu_table(ps_analysis),
  "matrix"
)

if (phyloseq::taxa_are_rows(ps_analysis)) {
  beta_mat <- t(beta_mat)
}

beta_rel <- sweep(
  beta_mat,
  1,
  rowSums(beta_mat),
  "/"
)

bray_CCA <- vegan::vegdist(
  beta_CCA,
  method = "bray"
)

pcoa_CCA <- ape::pcoa(bray_CCA)

pcoa_scores <- data.frame(
  SampleID = rownames(pcoa_CCA$vectors),
  PCoA1 = pcoa_CCA$vectors[, 1],
  PCoA2 = pcoa_CCA$vectors[, 2]
)

meta_CCA_for_join <- meta_CCA
meta_CCA_for_join$SampleID <- rownames(meta_CCA_for_join)

pcoa_scores <- dplyr::left_join(
  pcoa_scores,
  meta_CCA_for_join,
  by = "SampleID"
)
var_exp <- pcoa_CCA$values$Relative_eig * 100

var_exp[1:2]

#make the plot
ggplot2::ggplot(
  pcoa_scores,
  ggplot2::aes(
    x = PCoA1,
    y = PCoA2,
    color = CCA_Morphology
  )
) +
  ggplot2::geom_point(
    size = 4,
    alpha = 0.9
  ) +
  ggrepel::geom_text_repel(
    ggplot2::aes(label = SampleID),
    size = 3,
    show.legend = FALSE
  ) +
  ggplot2::labs(
    title = "Bray-Curtis PCoA of CCA-associated bacterial communities",
    x = paste0(
      "PCoA1 (",
      round(var_exp[1], 1),
      "%)"
    ),
    y = paste0(
      "PCoA2 (",
      round(var_exp[2], 1),
      "%)"
    ),
    color = "CCA morphology"
  ) +
  ggplot2::theme_classic()

#morphology PERMANOVA
meta_CCA_df <- as(meta_CCA, "data.frame")

identical(
  rownames(meta_CCA_df),
  labels(bray_CCA)
)
str(meta_CCA_df$Quantification)


permanova_morphology <- vegan::adonis2(
  bray_CCA ~ CCA_Morphology,
  data = meta_CCA_df,
  permutations = 999
)

permanova_morphology

#DNA concentration
permanova_DNA <- vegan::adonis2(
  bray_CCA ~ Quantification,
  data = meta_CCA_df,
  permutations = 999
)

permanova_DNA

#both together 
permanova_morph_DNA <- vegan::adonis2(
  bray_CCA ~ CCA_Morphology + Quantification,
  data = meta_CCA_df,
  permutations = 999,
  by = "margin"
)

permanova_morph_DNA

#beta dispersion
disp_CCA <- vegan::betadisper(
  bray_CCA,
  group = meta_CCA_df$CCA_Morphology,
  bias.adjust = TRUE
)

vegan::permutest(
  disp_CCA,
  permutations = 999
)

#morphology only: no sig difference between branching and encrusting 
#DNA conc: does not significantly explain community composition 
#morphology controlling for DNA concentration: morph still not significant 
#concentration controlling for morphology: conc not significant 
#beta dispersion: no evidence of unequal dispersion between morphologies 

################################################################################
#################CCA vs SEAWATER BRAY-CURTIS PCoA + PERMANOVA###################
################################################################################

bray_env <- vegan::vegdist(
  beta_rel,
  method = "bray"
)

meta_env_df <- as.data.frame(beta_meta)
identical(
  rownames(meta_env_df),
  labels(bray_env)
)

table(meta_env_df$Environment)

phyloseq::sample_data(ps_analysis)$Host_Genus[
  phyloseq::sample_data(ps_analysis)$Host_Genus == "Lithopyllum"
] <- "Lithophyllum"

meta_env_df <- as(
  phyloseq::sample_data(ps_analysis),
  "data.frame"
)

class(meta_env_df)
identical(
  rownames(meta_env_df),
  labels(bray_env)
)
table(meta_env_df$Environment)

#CCA vs water PERMANOVA
permanova_environment <- vegan::adonis2(
  bray_env ~ Environment,
  data = meta_env_df,
  permutations = 999
)

permanova_environment

#beta dispersion
disp_env <- vegan::betadisper(
  bray_env,
  group = meta_env_df$Environment,
  bias.adjust = TRUE
)

vegan::permutest(
  disp_env,
  permutations = 999
)

#plot CCA vs seawater PCoA
pcoa_env <- ape::pcoa(bray_env)

pcoa_env_scores <- data.frame(
  SampleID = rownames(pcoa_env$vectors),
  PCoA1 = pcoa_env$vectors[, 1],
  PCoA2 = pcoa_env$vectors[, 2]
)
#samples that plot cloes together have more similar bacterial communities; samples farther apart are more compositionally different

meta_env_join <- meta_env_df

meta_env_join$SampleID <- rownames(meta_env_join)

pcoa_env_scores <- dplyr::left_join(
  pcoa_env_scores,
  meta_env_join,
  by = "SampleID"
)

pcoa_env_scores[
  ,
  c("SampleID", "Environment")
]

var_exp_env <- pcoa_env$values$Relative_eig * 100

var_exp_env[1:2]
#this tells you how much of the total Bray-Curtis variation is represented by PCoA1 and PCoA2

ggplot2::ggplot(
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
    ggplot2::aes(label = SampleID),
    size = 3,
    max.overlaps = Inf,
    show.legend = FALSE
  ) +
  ggplot2::labs(
    title = "Bray-Curtis PCoA of CCA and Seawater Bacterial Communities",
    x = paste0(
      "PCoA1 (",
      round(var_exp_env[1], 1),
      "%)"
    ),
    y = paste0(
      "PCoA2 (",
      round(var_exp_env[2], 1),
      "%)"
    ),
    color = "Environment"
  ) +
  ggplot2::theme_classic()

#verifying story 
disp_env$group.distances
#this gives average distance to the group centroid for CCA and water

aggregate(
  disp_env$distances,
  by = list(Environment = disp_env$group),
  FUN = mean
)
#group with the larger value here has greater multivariate dispersion --> CCA has greater dispersion

boxplot(
  disp_env,
  xlab = "Environment",
  ylab = "Distance to group centroid",
  main = "Beta dispersion: CCA vs water"
)


#make PCoA plots with species and genus 
# Copy CCA metadata
meta_tax <- meta_CCA_df

# Add sample IDs from row names
meta_tax$SampleID <- rownames(meta_tax)

meta_tax$Species_plot <- dplyr::case_when(
  
  # No genus identification
  is.na(meta_tax$Host_Genus) |
    meta_tax$Host_Genus == "NF" ~
    "Unidentified",
  
  # Genus known, species not identified
  is.na(meta_tax$Host_Species) |
    meta_tax$Host_Species == "NF" ~
    paste0(meta_tax$Host_Genus, " sp."),
  
  # Genus + species both known
  TRUE ~
    paste(
      meta_tax$Host_Genus,
      meta_tax$Host_Species
    )
)

#adding genus
meta_tax$Genus_plot <- ifelse(
  is.na(meta_tax$Host_Genus) |
    meta_tax$Host_Genus == "NF",
  "Unidentified",
  meta_tax$Host_Genus
)

pcoa_tax_scores <- pcoa_scores %>%
  dplyr::select(
    SampleID,
    PCoA1,
    PCoA2
  ) %>%
  dplyr::left_join(
    meta_tax,
    by = "SampleID"
  )
pcoa_tax_scores[
  ,
  c(
    "SampleID",
    "CCA_Morphology",
    "Host_Genus",
    "Host_Species",
    "Species_plot",
    "Genus_plot"
  )
]
#coloring the existing CCA ordination by host taxonomic identity to see whether samples with similar host taxonomy occupy similar regions of community space

#PCoA colored by species and shape = morphology
ggplot2::ggplot(
  pcoa_tax_scores,
  ggplot2::aes(
    x = PCoA1,
    y = PCoA2,
    color = Species_plot,
    shape = CCA_Morphology
  )
) +
  ggplot2::geom_point(
    size = 4,
    alpha = 0.9
  ) +
  ggrepel::geom_text_repel(
    ggplot2::aes(label = SampleID),
    size = 3,
    show.legend = FALSE,
    max.overlaps = Inf
  ) +
  ggplot2::labs(
    title = "Bray-Curtis PCoA by CCA host species",
    x = paste0(
      "PCoA1 (",
      round(var_exp[1], 1),
      "%)"
    ),
    y = paste0(
      "PCoA2 (",
      round(var_exp[2], 1),
      "%)"
    ),
    color = "Host taxon",
    shape = "CCA morphology"
  ) +
  ggplot2::theme_classic()

#same thing at genus level
meta_CCA_df$Host_Genus[
  meta_CCA_df$Host_Genus == "Lithopyllum"
] <- "Lithophyllum"

ggplot2::ggplot(
  pcoa_tax_scores,
  ggplot2::aes(
    x = PCoA1,
    y = PCoA2,
    color = Genus_plot,
    shape = CCA_Morphology
  )
) +
  ggplot2::geom_point(
    size = 4,
    alpha = 0.9
  ) +
  ggrepel::geom_text_repel(
    ggplot2::aes(label = SampleID),
    size = 3,
    show.legend = FALSE,
    max.overlaps = Inf
  ) +
  ggplot2::labs(
    title = "Bray-Curtis PCoA by CCA host genus",
    x = paste0(
      "PCoA1 (",
      round(var_exp[1], 1),
      "%)"
    ),
    y = paste0(
      "PCoA2 (",
      round(var_exp[2], 1),
      "%)"
    ),
    color = "Host genus",
    shape = "CCA morphology"
  ) +
  ggplot2::theme_classic()
#morphology alone does not explain the major spread within the CCA samples 

####do encrusting samples differ from other encrusting 
meta_CCA_df$Lichenare_group <- ifelse(
  meta_CCA_df$Host_Species == "lichenare",
  "L_lithophyllum_lichenare",
  "Other"
)

encrusting_ids <- rownames(meta_CCA_df)[
  meta_CCA_df$CCA_Morphology == "encrusting"
]

beta_encrusting <- beta_CCA[
  encrusting_ids,
  ,
  drop = FALSE
]

meta_encrusting <- meta_CCA_df[
  encrusting_ids,
  ,
  drop = FALSE
]

table(meta_encrusting$Lichenare_group)

bray_encrusting <- vegan::vegdist(
  beta_encrusting,
  method = "bray"
)

permanova_lichenare <- vegan::adonis2(
  bray_encrusting ~ Lichenare_group,
  data = meta_encrusting,
  permutations = 999
)

permanova_lichenare


disp_lichenare <- vegan::betadisper(
  bray_encrusting,
  group = meta_encrusting$Lichenare_group,
  bias.adjust = TRUE
)

vegan::permutest(
  disp_lichenare,
  permutations = 999
)
