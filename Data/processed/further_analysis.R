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

############################################################
###################INDICATOR ANALYSIS#######################
############################################################

install.packages("indicspecies")
library(indicspecies)
library(permute)

lichenare_group <- factor(
  meta_encrusting$Lichenare_group,
  levels = c(
    "Other",
    "L_lithophyllum_lichenare"
  )
)

table(lichenare_group)

#removing taxa that are too rare to be informative
taxon_prevalence <- colSums(beta_encrusting > 0)
#beta_encrusting has relative abundance values, not raw read counts 

beta_encrusting_ind <- beta_encrusting[
  ,
  taxon_prevalence >= 3,
  drop = FALSE
]

dim(beta_encrusting_ind)
#this only keeps taxa found in at least 3 encrusting samples to help reduce noise --> got 12 encrusting samples x 558 bacterial taxa

#running actual indicator value test
set.seed(123)

indval_lichenare <- indicspecies::multipatt(
  beta_encrusting_ind,
  lichenare_group,
  func = "IndVal.g",
  duleg = TRUE,
  control = permute::how(nperm = 999)
)

summary(
  indval_lichenare,
  indvalcomp = TRUE,
  alpha = 0.05
)
#gave 16 bacterial taxa associated with lichenare and 28 with other encrusting samples 

#attaching taxonomy to the indicator taxa
tax_df <- as.data.frame(
  as(
    phyloseq::tax_table(ps_analysis),
    "matrix"
  )
)

tax_df$Taxon <- rownames(tax_df)

indval_results <- indval_lichenare$sign

indval_results$Taxon <- rownames(indval_results)

indval_results_tax <- dplyr::left_join(
  indval_results,
  tax_df,
  by = "Taxon"
)

colnames(indval_results_tax)

#pull out the 16 taxa significantlly associated with lichenare 
lichenare_indicators <- indval_results_tax %>%
  dplyr::filter(
    s.L_lithophyllum_lichenare == 1,
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
  dplyr::arrange(p.value)

lichenare_indicators

##################ugh ok the correction i did and deleted before is messing things up so remake the inval values table 
indval_raw <- indval_lichenare$sign

indval_raw$Taxon <- rownames(indval_raw)

indval_raw <- dplyr::left_join(
  indval_raw,
  tax_df,
  by = "Taxon"
)

lichenare_indicators_raw <- indval_raw %>%
  dplyr::filter(
    s.L_lithophyllum_lichenare == 1,
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
  dplyr::arrange(p.value)

lichenare_indicators_raw


##############trying to do branching vs encrusting indicator analysis##################
morph_group <- factor(
  meta_CCA_df$CCA_Morphology,
  levels = c("encrusting", "branching")
)

table(morph_group)

taxon_prevalence_morph <- colSums(beta_CCA > 0)

beta_CCA_ind <- beta_CCA[
  ,
  taxon_prevalence_morph >= 3,
  drop = FALSE
]

dim(beta_CCA_ind)

#taxa were retained if they occurred in at least 3 of the 14 CCA samples. removes extremely rare samples

set.seed(123)

indval_morph <- indicspecies::multipatt(
  beta_CCA_ind,
  morph_group,
  func = "IndVal.g",
  duleg = TRUE,
  control = permute::how(nperm = 999)
)

summary(
  indval_morph,
  indvalcomp = TRUE,
  alpha = 0.05
)
#found a smaller set of individual taxa associated with each morphology
#most interested in ones that have A and B = 1.000
#branching results are less secure just because there are only two samples 

#attach taxonomy
indval_morph_raw <- indval_morph$sign

indval_morph_raw$Taxon <- rownames(indval_morph_raw)

indval_morph_tax <- dplyr::left_join(
  indval_morph_raw,
  tax_df,
  by = "Taxon"
)

encrusting_indicators <- indval_morph_tax %>%
  dplyr::filter(
    s.encrusting == 1,
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
  dplyr::arrange(dplyr::desc(stat))

encrusting_indicators
#these taxa were specifically associated with encrusting CCA relative to branching CCA 

#want to make a stacked bar plot showing the relative abundance of just those indicator taxa across all CCA samples 
indicator_ids <- encrusting_indicators$Taxon

indicator_mat <- beta_CCA[
  ,
  indicator_ids,
  drop = FALSE
]

indicator_df <- as.data.frame(indicator_mat)

indicator_df$SampleID <- rownames(indicator_df)

indicator_long <- indicator_df %>%
  tidyr::pivot_longer(
    cols = -SampleID,
    names_to = "Taxon",
    values_to = "RelativeAbundance"
  )

indicator_tax_labels <- encrusting_indicators %>%
  dplyr::mutate(
    TaxonLabel = dplyr::case_when(
      !is.na(species) & species != "" &
        !is.na(genus) & genus != "" ~
        paste(genus, species),
      
      !is.na(genus) & genus != "" ~
        genus,
      
      TRUE ~ Taxon
    )
  ) %>%
  dplyr::select(
    Taxon,
    TaxonLabel
  )

indicator_long <- indicator_long %>%
  dplyr::left_join(
    indicator_tax_labels,
    by = "Taxon"
  ) %>%
  dplyr::left_join(
    meta_CCA_df %>%
      dplyr::mutate(
        SampleID = rownames(meta_CCA_df)
      ) %>%
      dplyr::select(
        SampleID,
        CCA_Morphology
      ),
    by = "SampleID"
  )


sample_order <- indicator_long %>%
  dplyr::distinct(
    SampleID,
    CCA_Morphology
  ) %>%
  dplyr::arrange(
    factor(
      CCA_Morphology,
      levels = c("encrusting", "branching")
    ),
    SampleID
  ) %>%
  dplyr::pull(SampleID)

indicator_long$SampleID <- factor(
  indicator_long$SampleID,
  levels = sample_order
)

sample_labels <- meta_CCA_df %>%
  dplyr::mutate(
    SampleID = rownames(meta_CCA_df),
    
    HostLabel = dplyr::case_when(
      
      # Genus and species identified
      !is.na(Host_Genus) &
        Host_Genus != "NF" &
        !is.na(Host_Species) &
        Host_Species != "NF" ~
        paste(Host_Genus, Host_Species),
      
      # Genus identified, species not identified
      !is.na(Host_Genus) &
        Host_Genus != "NF" ~
        paste0(Host_Genus, " sp."),
      
      # No taxonomic ID
      TRUE ~ "Unidentified"
    ),
    
    PlotLabel = paste0(
      SampleID,
      "\n",
      HostLabel
    )
  )

sample_label_vector <- setNames(
  sample_labels$PlotLabel,
  sample_labels$SampleID
)

ggplot2::ggplot(
  indicator_long,
  ggplot2::aes(
    x = SampleID,
    y = RelativeAbundance,
    fill = TaxonLabel
  )
) +
  ggplot2::geom_col() +
  ggplot2::facet_grid(
    ~ CCA_Morphology,
    scales = "free_x",
    space = "free_x"
  ) +
  ggplot2::scale_x_discrete(
    labels = sample_label_vector
  ) +
  ggplot2::scale_y_continuous(
    labels = scales::percent
  ) +
  ggplot2::labs(
    title = "Relative abundance of encrusting CCA indicator taxa",
    x = "Sample and host taxon",
    y = "Relative abundance",
    fill = "Indicator taxon"
  ) +
  ggplot2::theme_classic() +
  ggplot2::theme(
    axis.text.x = ggplot2::element_text(
      angle = 45,
      hjust = 1,
      size = 8
    )
  )
#the bars don't go to 100 because they show what fraction of the entire bacterial community is represented by these 9 candidate taxa
#this lets us see what is genuinely abundant or just statistically consistent at low abundance 


#getting rid of the EMU_taxon_78 and ordering things by how similar they are 

# ============================================================
# 8 ENCRUSTING INDICATOR TAXA
# Remove EMU_Taxon_78 and visualize relative abundance
# ============================================================


# ------------------------------------------------------------
# 1. Select the 8 indicator taxa
# ------------------------------------------------------------

indicator_ids_8 <- encrusting_indicators %>%
  dplyr::filter(Taxon != "EMU_Taxon_78") %>%
  dplyr::pull(Taxon)


# Pull their relative abundances from the CCA community table
# Rows = samples
# Columns = bacterial taxa

indicator_mat_8 <- beta_CCA[
  ,
  indicator_ids_8,
  drop = FALSE
]


# ------------------------------------------------------------
# 2. Order taxa by similarity in abundance pattern
# ------------------------------------------------------------

# IMPORTANT:
# taxa are COLUMNS, so do NOT transpose this matrix here

taxon_cor <- cor(
  indicator_mat_8,
  method = "spearman"
)

# Convert correlation to a distance
taxon_dist <- as.dist(1 - taxon_cor)

# Cluster taxa with similar abundance patterns
taxon_cluster <- hclust(
  taxon_dist,
  method = "average"
)

# Extract taxon order from clustering
taxon_order <- colnames(indicator_mat_8)[
  taxon_cluster$order
]

taxon_order

# ------------------------------------------------------------
# 3. Create readable bacterial taxon labels
# ------------------------------------------------------------

indicator_tax_labels_8 <- encrusting_indicators %>%
  dplyr::filter(Taxon != "EMU_Taxon_78") %>%
  dplyr::mutate(
    TaxonLabel = dplyr::case_when(
      
      # Genus + species available
      !is.na(genus) &
        genus != "" &
        !is.na(species) &
        species != "" ~
        paste(genus, species),
      
      # Only genus available
      !is.na(genus) &
        genus != "" ~
        genus,
      
      # Fallback
      TRUE ~ Taxon
    )
  ) %>%
  dplyr::select(
    Taxon,
    TaxonLabel
  )


# Match readable names to the similarity-based order

taxon_label_order <- indicator_tax_labels_8$TaxonLabel[
  match(
    taxon_order,
    indicator_tax_labels_8$Taxon
  )
]

taxon_label_order

# ------------------------------------------------------------
# 4. Convert to long dataframe
# ------------------------------------------------------------

indicator_df_8 <- as.data.frame(
  indicator_mat_8
)

indicator_df_8$SampleID <- rownames(
  indicator_df_8
)


indicator_long_8 <- indicator_df_8 %>%
  tidyr::pivot_longer(
    cols = -SampleID,
    names_to = "Taxon",
    values_to = "RelativeAbundance"
  ) %>%
  dplyr::left_join(
    indicator_tax_labels_8,
    by = "Taxon"
  )


# Apply similarity-based taxon order

indicator_long_8$TaxonLabel <- factor(
  indicator_long_8$TaxonLabel,
  levels = taxon_label_order
)


# ------------------------------------------------------------
# 5. Create sample labels with host taxonomy
# ------------------------------------------------------------

sample_labels <- meta_CCA_df %>%
  dplyr::mutate(
    SampleID = rownames(meta_CCA_df),
    
    HostLabel = dplyr::case_when(
      
      # Genus + species known
      !is.na(Host_Genus) &
        Host_Genus != "NF" &
        !is.na(Host_Species) &
        Host_Species != "NF" ~
        paste(
          Host_Genus,
          Host_Species
        ),
      
      # Only genus known
      !is.na(Host_Genus) &
        Host_Genus != "NF" ~
        paste0(
          Host_Genus,
          " sp."
        ),
      
      # No host identification
      TRUE ~
        "Unidentified"
    ),
    
    PlotLabel = paste0(
      SampleID,
      "\n",
      HostLabel
    )
  )


# Add morphology and host labels to abundance dataframe

indicator_long_8 <- indicator_long_8 %>%
  dplyr::left_join(
    sample_labels %>%
      dplyr::select(
        SampleID,
        CCA_Morphology,
        PlotLabel
      ),
    by = "SampleID"
  )

# ------------------------------------------------------------
# 6. Set sample order
# ------------------------------------------------------------

sample_order <- sample_labels %>%
  dplyr::filter(
    SampleID %in% rownames(beta_CCA)
  ) %>%
  dplyr::arrange(
    factor(
      CCA_Morphology,
      levels = c(
        "encrusting",
        "branching"
      )
    ),
    SampleID
  ) %>%
  dplyr::pull(SampleID)


indicator_long_8$SampleID <- factor(
  indicator_long_8$SampleID,
  levels = sample_order
)


# Named vector for x-axis labels

sample_label_vector <- setNames(
  sample_labels$PlotLabel,
  sample_labels$SampleID
)

table(
  indicator_long_8$TaxonLabel,
  useNA = "ifany"
)


# ------------------------------------------------------------
# 7. FACETED RELATIVE ABUNDANCE PLOT
# ------------------------------------------------------------

indicator_facet_plot <- ggplot2::ggplot(
  indicator_long_8,
  ggplot2::aes(
    x = SampleID,
    y = RelativeAbundance,
    fill = CCA_Morphology
  )
) +
  
  ggplot2::geom_col() +
  
  ggplot2::facet_wrap(
    ~ TaxonLabel,
    ncol = 2
  ) +
  
  ggplot2::scale_x_discrete(
    labels = sample_label_vector
  ) +
  
  ggplot2::scale_y_continuous(
    labels = scales::percent
  ) +
  
  ggplot2::labs(
    title = "Relative abundance of encrusting CCA indicator taxa",
    x = "Sample and host taxon",
    y = "Relative abundance",
    fill = "CCA morphology"
  ) +
  
  ggplot2::theme_classic() +
  
  ggplot2::theme(
    axis.text.x = ggplot2::element_text(
      angle = 45,
      hjust = 1,
      size = 6
    ),
    strip.text = ggplot2::element_text(
      face = "italic"
    )
  )


indicator_facet_plot


########################################
###doing this with water samples too####
########################################
meta_three <- as.data.frame(
  phyloseq::sample_data(ps_analysis)
)

meta_three$Group3 <- dplyr::case_when(
  meta_three$Environment == "Water" ~ "Water",
  meta_three$CCA_Morphology == "branching" ~ "Branching",
  meta_three$CCA_Morphology == "encrusting" ~ "Encrusting",
  TRUE ~ NA_character_
)

meta_three$Group3 <- factor(
  meta_three$Group3,
  levels = c(
    "Encrusting",
    "Branching",
    "Water"
  )
)

table(meta_three$Group3)


taxon_prevalence_three <- colSums(beta_rel > 0)

beta_three_ind <- beta_rel[
  ,
  taxon_prevalence_three >= 3,
  drop = FALSE
]

dim(beta_three_ind)


group_three <- meta_three$Group3

names(group_three) <- rownames(meta_three)

group_three <- group_three[
  rownames(beta_three_ind)
]

table(group_three)



set.seed(123)

indval_three <- indicspecies::multipatt(
  beta_three_ind,
  group_three,
  func = "IndVal.g",
  duleg = TRUE,
  control = permute::how(nperm = 999)
)

summary(
  indval_three,
  indvalcomp = TRUE,
  alpha = 0.05
)


# Pull original, uncorrected indicator results
indval_three_raw <- indval_three$sign

indval_three_raw$Taxon <- rownames(indval_three_raw)

# Attach bacterial taxonomy
indval_three_tax <- dplyr::left_join(
  indval_three_raw,
  tax_df,
  by = "Taxon"
)

# Extract only Encrusting indicators
encrusting_indicators_three <- indval_three_tax %>%
  dplyr::filter(
    s.Encrusting == 1,
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

encrusting_indicators_three

#stacked relative-abundance chart across all 17 biological samples
indicator_ids_three <- encrusting_indicators_three %>%
  dplyr::filter(Taxon != "EMU_Taxon_78") %>%
  dplyr::pull(Taxon)

# Pull relative abundance values from ALL 17 biological samples
indicator_mat_three <- beta_rel[
  ,
  indicator_ids_three,
  drop = FALSE
]

indicator_tax_labels_three <- encrusting_indicators_three %>%
  dplyr::filter(Taxon != "EMU_Taxon_78") %>%
  dplyr::mutate(
    TaxonLabel = dplyr::case_when(
      !is.na(genus) & genus != "" &
        !is.na(species) & species != "" ~
        paste(genus, species),
      
      !is.na(genus) & genus != "" ~
        genus,
      
      TRUE ~ Taxon
    )
  ) %>%
  dplyr::select(
    Taxon,
    TaxonLabel
  )



taxon_cor_three <- cor(
  indicator_mat_three,
  method = "spearman"
)

taxon_dist_three <- as.dist(
  1 - taxon_cor_three
)

taxon_cluster_three <- hclust(
  taxon_dist_three,
  method = "average"
)

taxon_order_three <- colnames(indicator_mat_three)[
  taxon_cluster_three$order
]

taxon_label_order_three <-
  indicator_tax_labels_three$TaxonLabel[
    match(
      taxon_order_three,
      indicator_tax_labels_three$Taxon
    )
  ]


indicator_df_three <- as.data.frame(
  indicator_mat_three
)

indicator_df_three$SampleID <- rownames(
  indicator_df_three
)

indicator_long_three <- indicator_df_three %>%
  tidyr::pivot_longer(
    cols = -SampleID,
    names_to = "Taxon",
    values_to = "RelativeAbundance"
  ) %>%
  dplyr::left_join(
    indicator_tax_labels_three,
    by = "Taxon"
  )

indicator_long_three$TaxonLabel <- factor(
  indicator_long_three$TaxonLabel,
  levels = taxon_label_order_three
)


# -------------------------------
# Rebuild meta_three as a plain data frame
# -------------------------------

meta_three_df <- data.frame(
  phyloseq::sample_data(ps_analysis),
  stringsAsFactors = FALSE,
  check.names = FALSE
)

meta_three_df$SampleID <- rownames(meta_three_df)

meta_three_df$Group3 <- dplyr::case_when(
  meta_three_df$Environment == "Water" ~ "Water",
  meta_three_df$CCA_Morphology == "branching" ~ "Branching",
  meta_three_df$CCA_Morphology == "encrusting" ~ "Encrusting",
  TRUE ~ NA_character_
)

meta_three_df$Group3 <- factor(
  meta_three_df$Group3,
  levels = c("Encrusting", "Branching", "Water")
)

# -------------------------------
# Create sample labels
# -------------------------------

sample_labels_three <- meta_three_df %>%
  dplyr::mutate(
    HostLabel = dplyr::case_when(
      Group3 == "Water" ~ "Water",
      
      !is.na(Host_Genus) &
        Host_Genus != "NF" &
        !is.na(Host_Species) &
        Host_Species != "NF" ~
        paste(Host_Genus, Host_Species),
      
      !is.na(Host_Genus) &
        Host_Genus != "NF" ~
        paste0(Host_Genus, " sp."),
      
      TRUE ~ "Unidentified"
    ),
    
    PlotLabel = paste0(
      SampleID,
      "\n",
      HostLabel
    )
  )

# -------------------------------
# Join labels back to plotting dataframe
# -------------------------------

indicator_long_three <- indicator_long_three %>%
  dplyr::left_join(
    sample_labels_three %>%
      dplyr::select(
        SampleID,
        Group3,
        PlotLabel
      ),
    by = "SampleID"
  )

# -------------------------------
# Order samples: Encrusting -> Branching -> Water
# -------------------------------

sample_order_three <- sample_labels_three %>%
  dplyr::arrange(Group3, SampleID) %>%
  dplyr::pull(SampleID)

indicator_long_three$SampleID <- factor(
  indicator_long_three$SampleID,
  levels = sample_order_three
)

sample_label_vector_three <- setNames(
  sample_labels_three$PlotLabel,
  sample_labels_three$SampleID
)



indicator_three_stacked <- ggplot2::ggplot(
  indicator_long_three,
  ggplot2::aes(
    x = SampleID,
    y = RelativeAbundance,
    fill = TaxonLabel
  )
) +
  ggplot2::geom_col() +
  ggplot2::facet_grid(
    ~ Group3,
    scales = "free_x",
    space = "free_x"
  ) +
  ggplot2::scale_x_discrete(
    labels = sample_label_vector_three
  ) +
  ggplot2::scale_y_continuous(
    labels = scales::percent
  ) +
  ggplot2::labs(
    title = "Relative abundance of encrusting CCA indicator taxa",
    subtitle = "Encrusting, branching, and water samples",
    x = "Sample and host taxon",
    y = "Relative abundance",
    fill = "Indicator taxon"
  ) +
  ggplot2::theme_classic() +
  ggplot2::theme(
    axis.text.x = ggplot2::element_text(
      angle = 45,
      hjust = 1,
      size = 7
    )
  )

indicator_three_stacked
