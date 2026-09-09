#  ooooooooo.                      oooo                                      .                  oooooo     oooo  .oooooo..o      ooooo                                                                                          .                                
#  `888   `Y88.                    `888                                    .o8                   `888.     .8'  d8P'    `Y8      `888'                                                                                        .o8                                
#   888   .d88'           .ooooo.   888   .oooo.   oooo    ooo  .oooo.   .o888oo  .oooo.          `888.   .8'   Y88bo.            888                     .oooo.o  .oooo.   oooo d8b ooo. .oo.  .oo.   .ooooo.  ooo. .oo.   .o888oo .ooooo.  .oooo.o   .oooo.   
#   888ooo88P'           d88' `"Y8  888  `P  )88b   `88.  .8'  `P  )88b    888   `P  )88b          `888. .8'     `"Y8888o.        888                   d88(  "8 `P  )88b  `888""8P `888P"Y88bP"Y88b  d88' `88b `888P"Y88b    888   d88' `88b d88(  "8 `P  )88b  
#   888                  888        888   .oP"888    `88..8'    .oP"888    888    .oP"888           `888.8'          `"Y88b       888                   `"Y88b.   .oP"888   888      888   888   888  888ooo888  888   888    888   888   888 `"Y88b.   .oP"888  
#   888         .o.      888   .o8  888  d8(  888     `888'    d8(  888    888 . d8(  888            `888'      oo     .d8P       888       o  .o.      o.  )88b d8(  888   888      888   888   888  888    .o  888   888    888 . 888   888 o.  )88b d8(  888  
#  o888o        Y8P      `Y8bod8P' o888o `Y888""8o     `8'     `Y888""8o   "888" `Y888""8o            `8'       8""88888P'       o888ooooood8  Y8P      8""888P' `Y888""8o d888b    o888o o888o o888o `Y8bod8P' o888o o888o   "888" `Y8bod8P' 8""888P' `Y888""8o    ALESSIA__AVESANI



# Packages
library(readxl)
library(tidyverse)
library(vegan)
library(ape)
library(decontam)
library(ggrepel)



# Original datasets
ASV_counts_Pcla_Lsar <- read_delim("ASV_counts_Pcla_Lsar.csv", 
                                   delim = ",", escape_double = FALSE, trim_ws = TRUE)

ASVs_taxonomy <- read_delim("C:/Users/AlessiaAvesani/OneDrive - unige.it/Desktop/Università/Dottorato/PROGETTI CORALLI/bioinformatic results/results/ASVs_taxonomy.tsv", 
                            delim = "\t", escape_double = FALSE, 
                            trim_ws = TRUE)

metadata_Pclav_Lsar <- read_excel("Pclavata + Lsarmentosa/metadata_Pclav_Lsar.xlsx", 
                                  col_types = c("text", "text", "text", 
                                                "text", "text", "date", "numeric", 
                                                "numeric", "numeric", "text", "numeric", 
                                                "text", "text", "numeric"))



# ┌────────────────────────────────────────────────────────────────────────────┐
# │                            CLEANING THE DATASET                            │
# └────────────────────────────────────────────────────────────────────────────┘

# 1 - Remove mitochondria and chloroplasts =======================================


# 1.1 - Quantify how many reads are Mitochondria and Chloroplast _________________

ASVs_mitochondria <- ASVs_taxonomy %>%
  filter(family == "Mitochondria")

ASVs_chloroplast <- ASVs_taxonomy %>%
  filter(order == "Chloroplast")

ASVs_countsXmitoechloro <- ASV_counts_Pcla_Lsar.csv #to avoid working on the original table, just in case

# 1) Create a vector with ASV names from the Mito and Chloro taxonomy tables
chloro_ids <- ASVs_chloroplast$ASVs
mito_ids   <- ASVs_mitochondria$ASVs

# 2) Rename ASVs in the dedicated count table
counts_labeled <- ASVs_countsXmitoechloro %>%
  mutate(ASVs = case_when(
    ASVs %in% mito_ids   ~ "Mitochondria",
    ASVs %in% chloro_ids ~ "Chloroplast",
    TRUE                 ~ ASVs
  ))

# 3) Aggregate (sum) rows that now have the same name (Mito/Chloro)
counts_agg <- counts_labeled %>%
  group_by(ASVs) %>%
  summarise(across(where(is.numeric), sum, na.rm = TRUE))

# 4) Relative abundance % per sample
relabund_pct <- counts_agg %>%
  mutate(across(where(is.numeric), ~ .x / sum(.x, na.rm = TRUE) * 100))

# 5) Final table only with Mitochondria and Chloroplast (in %)
mito_chloro_pct <- relabund_pct %>%
  filter(ASVs %in% c("Mitochondria", "Chloroplast"))

mito_chloro_pct

write.csv(mito_chloro_pct, "mito_chloro_pct_Pclav_Lsar.csv", row.names = FALSE)

# 6) Plot only Mito and Chloro
datasetXplot_mito_chloro <- mito_chloro_pct %>%
  pivot_longer(-ASVs, names_to = "Sample", values_to = "Percent")

ggplot(datasetXplot_mito_chloro, aes(x = Sample, y = Percent, fill = ASVs)) +
  geom_col() +
  labs(
    title = "Relative abundance (%) of Mitochondria and Chloroplast P. clavata VS L.sarmentosa",
    x = "Sample",
    y = "Relative abundance (%)",
    fill = NULL
  ) +
  theme_minimal() +
  theme(
    axis.text.x = element_text(angle = 90, vjust = 0.5, hjust = 1)
  )


# 1.2 - Remove Mitochondria and Chloroplast reads ________________________________

# Vector with ASVs to remove
remove_ASVs <- union(ASVs_mitochondria$ASVs, ASVs_chloroplast$ASVs)

length(remove_ASVs) #check how many ASVs there are

# Create a new ASV_counts table without Mitochondria and Chloroplast
ASVs_counts_Pcla_Lsar_noMitoChloro <- ASV_counts_Pcla_Lsar.csv %>%
  filter(!(ASVs %in% remove_ASVs))

nrow(ASVs_counts_Pcla_Lsar_noMitoChloro) #remaining ASVs


# 2 - Remove ASVs not present in any sample ======================================

# ==> NB: I do this because ASV_counts is a subset of a larger dataset 
# that included other samples with a microbial composition different from these 
# so, in my modest opinion, it is appropriate to perform this step 
# (not mandatory) to slim down the dataset

ASVs_counts_Pcla_Lsar_clean <- ASVs_counts_Pcla_Lsar_noMitoChloro %>%
  filter(rowSums(across(-ASVs), na.rm = TRUE) > 0)

nrow(ASVs_counts_Pcla_Lsar_clean) #remaining ASVs ==> I basically halved the dataset


# 3 - Remove singleton ASVs from ASV_counts ====================================

ASVs_counts_Pcla_Lsar_noSingletons <- ASVs_counts_Pcla_Lsar_clean %>%
  filter(rowSums(across(-ASVs), na.rm = TRUE) > 1)

nrow(ASVs_counts_Pcla_Lsar_noSingletons) #remaining ASVs


# 4 - Remove ASVs with no assigned phylum from ASV_counts ======================

ASVs_taxonomy_noPhylum <- ASVs_taxonomy %>%
  filter(is.na(phylum) | phylum == "")

ASVs_list_noPhylum <- ASVs_taxonomy_noPhylum %>%
  pull(ASVs) #create the ASV list with no assigned phylum

ASVs_counts_Pcla_Lsar_noSingletons_noPhylum <- ASVs_counts_Pcla_Lsar_noSingletons %>%
  filter(!(ASVs %in% ASVs_list_noPhylum))

nrow(ASVs_counts_Pcla_Lsar_noSingletons_noPhylum) #remaining ASVs


# 5 - Remove ASVs assigned to family Propionibacteriaceae ======================

ASVs_taxonomy_Propionibacteriaceae <- ASVs_taxonomy %>%
  filter(family == "Propionibacteriaceae")

ASVs_list_Propionibacteriaceae <- ASVs_taxonomy_Propionibacteriaceae %>%
  pull(ASVs) #create the ASV list assigned to family Propionibacteriaceae

ASVs_counts_Pcla_Lsar_noSingletons_noPhylum_noPropionibacteriaceae <- ASVs_counts_Pcla_Lsar_noSingletons_noPhylum %>%
  filter(!(ASVs %in% ASVs_list_Propionibacteriaceae))

nrow(ASVs_counts_Pcla_Lsar_noSingletons_noPhylum_noPropionibacteriaceae) #remaining ASVs


# 6 - Remove non-bacterial ASVs from ASV_counts ================================

ASVs_taxonomy_Eucarioti_Archea <- ASVs_taxonomy %>%
  filter(domain == "Eukaryota" | domain == "Archaea" | is.na(domain) | domain == "") #taxonomy only for eukaryotes and archaea

ASVs_list_Euk_Archaea <- ASVs_taxonomy_Eucarioti_Archea %>%
  pull(ASVs) #create the ASV list for eukaryotes and archaea

#clean ASV_counts without Eukaryota and Archaea
ASVs_counts_Pcla_Lsar_Bacteria_clean <- ASVs_counts_Pcla_Lsar_noSingletons_noPhylum_noPropionibacteriaceae %>%
  filter(!(ASVs %in% ASVs_list_Euk_Archaea))






# ============================== CHECK POINT ==============================

# ======================================================= #
#   CHECK POINT: How many reads are left in the samples   #
# ======================================================= #

reads_per_sample <- ASVs_counts_Pcla_Lsar_Bacteria_clean %>%
  select(-ASVs) %>%
  summarise(across(everything(), sum))

# ==> nicer format
reads_per_sample_post_rimozione_chloro_mito_nonBatteri <- reads_per_sample %>%
  pivot_longer(cols = everything(),
               names_to = "Sample",
               values_to = "Total_reads")

matrix_X_info_ASVREADS_before_rem_cont <- ASVs_counts_Pcla_Lsar_Bacteria_clean[, -1]

N_ASVs <- colSums(matrix_X_info_ASVREADS_before_rem_cont > 0, na.rm = TRUE)

asvs_per_sample <- data.frame(
  Sample = names(N_ASVs),
  N_ASVs = as.integer(N_ASVs)
)

reads_ASVs_per_sample_before_rem_cont <- merge(
  reads_per_sample_post_rimozione_chloro_mito_nonBatteri,
  asvs_per_sample,
  by = "Sample",
  all.x = TRUE
)

write.csv(reads_ASVs_per_sample_before_rem_cont, "reads_ASVs_per_sample_before_rem_cont.csv", row.names = FALSE)



# 7 - Remove possible contaminant ASVs ==========================================

# >>>>> PUT THE EXACT NAME of the DNA column in your metadata here <<<<<
Quantification <- "Quantification"

# 1) Align samples between counts and metadata
samples_counts <- setdiff(colnames(ASVs_counts_Pcla_Lsar_Bacteria_clean), "ASVs")

metadata_Pcla_Lsar_new_forcontaminants <- metadata_Pclav_Lsar %>%
  filter(IIT_ID %in% samples_counts) %>%
  distinct(IIT_ID, .keep_all = TRUE)

# 2) Create seqtab in the format required by decontam: samples x ASV
mat_asv_x_samp <- ASVs_counts_Pcla_Lsar_Bacteria_clean %>%
  column_to_rownames("ASVs") %>%
  as.matrix()

# reorder columns (samples) as in metadata
mat_asv_x_samp <- mat_asv_x_samp[, metadata_Pcla_Lsar_new_forcontaminants$IIT_ID, drop = FALSE]

# transposed: rows = samples, columns = ASV  (this is the decontam format)
seqtab <- t(mat_asv_x_samp)

# 3) Prepare DNA concentration (conc) aligned to seqtab rows
conc <- metadata_Pcla_Lsar_new_forcontaminants[[Quantification]]
names(conc) <- metadata_Pcla_Lsar_new_forcontaminants$IIT_ID

# conc must be numeric and > 0
conc <- as.numeric(conc)

# 4a) decontam frequency t = 0.2
# threshold = 0.1 (default) -> mild
# 0.2 -> more permissive (removes more), 0.5 -> very aggressive       ----------------> VALUE TO STUDY!!!!
contam_res02 <- isContaminant(seqtab, method = "frequency", conc = conc,
                              threshold = 0.2, normalize = TRUE, detailed = TRUE)

table(contam_res02$contaminant) # ==> NB: ASVs left after decontamination (False)

contam_asvs02 <- rownames(contam_res02)[contam_res02$contaminant]
length(contam_asvs02)

# 4b) decontam frequency t = 0.5
# threshold = 0.1 (default) -> mild
# 0.2 -> more permissive (removes more), 0.5 -> very aggressive       ----------------> VALUE TO STUDY!!!!
contam_res05 <- isContaminant(seqtab, method = "frequency", conc = conc,
                              threshold = 0.5, normalize = TRUE, detailed = TRUE)

table(contam_res05$contaminant) # ==> NB: ASVs left after decontamination (False)

contam_asvs05 <- rownames(contam_res05)[contam_res05$contaminant]
length(contam_asvs05)

# 5) Remove contaminated ASVs from your original count table (ASV x samples)

#threshold = 0.2
ASVs_counts_decontam_02 <- ASVs_counts_Pcla_Lsar_Bacteria_clean %>%
  filter(!(ASVs %in% contam_asvs))

#threshold = 0.5
ASVs_counts_decontam_05 <- ASVs_counts_Pcla_Lsar_Bacteria_clean %>%
  filter(!(ASVs %in% contam_asvs))


# 6) Useful outputs
# Table with p-value/probability (detailed=TRUE) and contaminant flag
contam_table02 <- contam_res02 %>%   #threshold = 0.2
  rownames_to_column("ASVs") %>%
  arrange(desc(contaminant), p)

contam_table02  #threshold = 0.2

contam_table05 <- contam_res05 %>%   #threshold = 0.5
  rownames_to_column("ASVs") %>%
  arrange(desc(contaminant), p)

contam_table05  #threshold = 0.5

# save contaminant list
write.csv(contam_table02, "decontam_frequency_results_Pcla_Lsar02.csv", row.names = FALSE)
write.csv(ASVs_counts_decontam_02, "ASVs_counts_decontam_02_Pcla_Lsar.csv", row.names = FALSE)

write.csv(contam_table05, "decontam_frequency_results_Pcla_Lsar05.csv", row.names = FALSE)
write.csv(ASVs_counts_decontam_05, "ASVs_counts_decontam_05_Pcla_Lsar.csv", row.names = FALSE)


# =============================================================== #
#   CHECK POINT: How many reads and ASVs are left in the samples  #
# =============================================================== #

############################### threshold = 0.2 ################################

reads_per_sample_postdecontaminazione02 <- ASVs_counts_decontam_02 %>%
  select(-ASVs) %>%
  summarise(across(everything(), sum))

# ==> nicer format
reads_per_sample_postdecontaminazione02 <- reads_per_sample_postdecontaminazione02 %>%
  pivot_longer(cols = everything(),
               names_to = "Sample",
               values_to = "Total_reads")

write.csv(reads_per_sample_postdecontaminazione02, "reads_per_sample_postdecontaminazione02.csv", row.names = FALSE)

# ASV x samples matrix (post-decontam), with ASVs as the first column
matrix_X_info_ASVREADS02 <- ASVs_counts_Pcla_Lsar_readytouse[, -1]

N_ASVs02 <- colSums(matrix_X_info_ASVREADS02 > 0, na.rm = TRUE)

asvs_per_sample02 <- data.frame(
  Sample = names(N_ASVs),
  N_ASVs = as.integer(N_ASVs)
)

reads_ASVs_per_sample_postdecontaminazione02 <- merge(
  reads_per_sample_postdecontaminazione02,
  asvs_per_sample,
  by = "Sample",
  all.x = TRUE
)

write.csv(reads_ASVs_per_sample_postdecontaminazione02, "reads_ASVs_per_sample_postdecontaminazione02.csv", row.names = FALSE)

############################### threshold = 0.5 ################################

reads_per_sample_postdecontaminazione05 <- ASVs_counts_decontam_05 %>%
  select(-ASVs) %>%
  summarise(across(everything(), sum))

# ==> nicer format
reads_per_sample_postdecontaminazione05 <- reads_per_sample_postdecontaminazione05 %>%
  pivot_longer(cols = everything(),
               names_to = "Sample",
               values_to = "Total_reads")

write.csv(reads_per_sample_postdecontaminazione05, "reads_per_sample_postdecontaminazione05.csv", row.names = FALSE)

# ASV x samples matrix (post-decontam), with ASVs as the first column
matrix_X_info_ASVREADS05 <- ASVs_counts_Pcla_Lsar_readytouse[, -1]

N_ASVs05 <- colSums(matrix_X_info_ASVREADS05 > 0, na.rm = TRUE)

asvs_per_sample05 <- data.frame(
  Sample = names(N_ASVs),
  N_ASVs = as.integer(N_ASVs)
)

reads_ASVs_per_sample_postdecontaminazione05 <- merge(
  reads_per_sample_postdecontaminazione05,
  asvs_per_sample,
  by = "Sample",
  all.x = TRUE
)

write.csv(reads_ASVs_per_sample_postdecontaminazione05, "reads_ASVs_per_sample_postdecontaminazione05.csv", row.names = FALSE)





##### IMPORTANT: UPDATED ASVs_COUNTS TABLE ###################

ASVs_counts_Pcla_Lsar_readytouse_withLS98 <- ASVs_counts_decontam_02
write.csv(ASVs_counts_Pcla_Lsar_readytouse_withLS98, "ASVs_counts_Pcla_Lsar_readytouse_withLS98.csv", row.names = FALSE)

metadata_Pclav_Lsar_noLS98 <- metadata_Pclav_Lsar %>% filter(IIT_ID != "LS98")
write.csv(metadata_Pclav_Lsar_noLS98, "metadata_Pclav_Lsar_noLS98.csv", row.names = FALSE)

ASVs_counts_Pcla_Lsar_readytouse <- ASVs_counts_Pcla_Lsar_readytouse_withLS98 %>% select(-any_of("LS98"))
write.csv(ASVs_counts_Pcla_Lsar_readytouse, "ASVs_counts_Pcla_Lsar_readytouse.csv", row.names = FALSE)


#####




# ┌────────────────────────────────────────────────────────────────────────────┐
# │                              RAREFACTION CURVE                             │
# └────────────────────────────────────────────────────────────────────────────┘

# 1) AFTER cleaning ------------------------------------------------------------

library(dplyr)
library(tibble)
library(vegan)
library(ggplot2)
library(purrr)


meta <- metadata_Pclav_Lsar  #metadata

mat <- ASVs_counts_Pcla_Lsar_readytouse_withLS98 %>%  #ASVs_counts
  select(-ASVs) %>%
  mutate(across(everything(), ~ suppressWarnings(as.numeric(.)))) %>%
  select(where(~ all(!is.na(.)))) %>%     # keep only truly numeric columns
  as.matrix() %>%
  t()

# observed-point data
df_pts <- data.frame(
  IIT_ID  = rownames(mat),
  reads   = rowSums(mat),
  asv_obs = specnumber(mat)
) %>%
  left_join(meta %>% select(IIT_ID, Specie), by = "IIT_ID")

# curves for each sample up to its OWN depth
# n_points = how many points you want along the curve (more = smoother, slower)
n_points <- 35
min_start <- 0  # starting point of the curve

df_curve <- map_dfr(df_pts$IIT_ID, function(id){
  counts <- mat[id, ]
  N <- sum(counts)
  
  # handle the case where a sample has few reads
  start <- min(min_start, N)
  depths <- unique(round(seq(start, N, length.out = n_points)))
  depths <- depths[depths > 0 & depths <= N]
  
  # vegan::rarefy wants a community matrix (rows = samples)
  m1 <- matrix(counts, nrow = 1)
  
  data.frame(
    IIT_ID  = id,
    depth   = depths,
    asv_exp = as.numeric(vegan::rarefy(m1, sample = depths))
  )
}) %>%
  left_join(meta %>% select(IIT_ID, Specie), by = "IIT_ID")

# palette
spec_levels <- c("Paramuricea_clavata",
                 "Leptogorgia_sarmentosa",
                 "Paramuricea_clavata_(malata)")

pal <- c("Paramuricea_clavata" = "#1b9e77",
         "Leptogorgia_sarmentosa" = "#7570b3",
         "Paramuricea_clavata_(malata)" = "#d95f02")

df_curve$Specie <- factor(df_curve$Specie, levels = spec_levels)
df_pts$Specie   <- factor(df_pts$Specie,   levels = spec_levels)

# plot curves per sample + observed point, faceted by Species
ggplot(df_curve, aes(x = depth, y = asv_exp, group = IIT_ID, color = Specie)) +
  geom_line(alpha = 0.35, linewidth = 0.6) +
  geom_point(data = df_pts, aes(x = reads, y = asv_obs, color = Specie),
             inherit.aes = FALSE, size = 1.6, alpha = 0.85) +
  geom_text_repel(data = df_pts,
                  aes(x = reads, y = asv_obs, label = IIT_ID, color = Specie),
                  inherit.aes = FALSE,
                  size = 2.6,
                  max.overlaps = Inf,
                  box.padding = 0.25,
                  point.padding = 0.15,
                  show.legend = FALSE) +
  scale_color_manual(values = pal, drop = FALSE) +
  facet_wrap(~Specie, scales = "free_x") +
  labs(x = "Reads",
       y = "ASVs",
       title = "Rarefaction curves per campione fino alla depth reale + punto osservato") +
  theme_bw() +
  theme(legend.position = "none")


# 1) BEFORE cleaning -----------------------------------------------------------

library(dplyr)
library(tibble)
library(vegan)
library(ggplot2)
library(purrr)


meta <- metadata_Pclav_Lsar  #metadata

mat <- ASV_counts_Pcla_Lsar %>%  #ASVs_counts
  select(-ASVs) %>%
  mutate(across(everything(), ~ suppressWarnings(as.numeric(.)))) %>%
  select(where(~ all(!is.na(.)))) %>%     # keep only truly numeric columns
  as.matrix() %>%
  t()

# observed-point data
df_pts <- data.frame(
  IIT_ID  = rownames(mat),
  reads   = rowSums(mat),
  asv_obs = specnumber(mat)
) %>%
  left_join(meta %>% select(IIT_ID, Specie), by = "IIT_ID")

# curves for each sample up to its OWN depth
# n_points = how many points you want along the curve (more = smoother, slower)
n_points <- 35
min_start <- 0  # starting point of the curve

df_curve <- map_dfr(df_pts$IIT_ID, function(id){
  counts <- mat[id, ]
  N <- sum(counts)
  
  # handle the case where a sample has few reads
  start <- min(min_start, N)
  depths <- unique(round(seq(start, N, length.out = n_points)))
  depths <- depths[depths > 0 & depths <= N]
  
  # vegan::rarefy wants a community matrix (rows = samples)
  m1 <- matrix(counts, nrow = 1)
  
  data.frame(
    IIT_ID  = id,
    depth   = depths,
    asv_exp = as.numeric(vegan::rarefy(m1, sample = depths))
  )
}) %>%
  left_join(meta %>% select(IIT_ID, Specie), by = "IIT_ID")

# palette
spec_levels <- c("Paramuricea_clavata",
                 "Leptogorgia_sarmentosa",
                 "Paramuricea_clavata_(malata)")

pal <- c("Paramuricea_clavata" = "#1b9e77",
         "Leptogorgia_sarmentosa" = "#7570b3",
         "Paramuricea_clavata_(malata)" = "#d95f02")

df_curve$Specie <- factor(df_curve$Specie, levels = spec_levels)
df_pts$Specie   <- factor(df_pts$Specie,   levels = spec_levels)

# plot curves per sample + observed point, faceted by Species
ggplot(df_curve, aes(x = depth, y = asv_exp, group = IIT_ID, color = Specie)) +
  geom_line(alpha = 0.35, linewidth = 0.6) +
  geom_point(data = df_pts, aes(x = reads, y = asv_obs, color = Specie),
             inherit.aes = FALSE, size = 1.6, alpha = 0.85) +
  geom_text_repel(data = df_pts,
                  aes(x = reads, y = asv_obs, label = IIT_ID, color = Specie),
                  inherit.aes = FALSE,
                  size = 2.6,
                  max.overlaps = Inf,
                  box.padding = 0.25,
                  point.padding = 0.15,
                  show.legend = FALSE) +
  scale_color_manual(values = pal, drop = FALSE) +
  facet_wrap(~Specie, scales = "free_x") +
  labs(x = "Reads",
       y = "ASVs",
       title = "Rarefaction curves per campione fino alla depth reale + punto osservato") +
  theme_bw() +
  theme(legend.position = "none")

##############


# ┌────────────────────────────────────────────────────────────────────────────┐
# │                   ALPHA-DIVERSITY - SAMPLES NOT MERGED                     │
# └────────────────────────────────────────────────────────────────────────────┘

#####

# ASV x samples matrix
matrix_ASVs_counts_reverse <- ASVs_counts_Pcla_Lsar_readytouse %>%
  column_to_rownames("ASVs") %>%   #set ASV IDs as rownames
  as.matrix()

# reorder matrix columns as in metadata
matrixdiversity <- matrix_ASVs_counts_reverse[, metadata_Pclav_Lsar_noLS98 $IIT_ID, drop = FALSE]

# invert rows and columns of the metadata matrix
matrix_ASV_a_diversity <- t(matrixdiversity)
nrow(matrixdiversity)

# ASV matrix with alpha diversity
alpha_df <- data.frame(
  IIT_ID = rownames(matrix_ASV_a_diversity),
  Observed = vegan::specnumber(matrix_ASV_a_diversity),
  Shannon  = vegan::diversity(matrix_ASV_a_diversity, index = "shannon"),
  Simpson  = vegan::diversity(matrix_ASV_a_diversity, index = "simpson")
) %>%
  left_join(metadata_Pclav_Lsar_noLS98, by = "IIT_ID")

alpha_df

# all alpha diversity plots together
alpha_long <- alpha_df %>%
  select(IIT_ID, Specie, Observed, Shannon, Simpson) %>%
  pivot_longer(cols = c(Observed, Shannon, Simpson),
               names_to = "Metric",
               values_to = "Value")


ggplot(alpha_long, aes(x = Specie, y = Value)) +
  geom_boxplot(outlier.shape = NA) +
  geom_jitter(width = 0.15, alpha = 0.6, size = 1.5) +
  #  geom_text_repel(
  #    aes(label = IIT_ID),
  #    size = 2.7,
  #    max.overlaps = 70,   # increase if it hides too many
  #    box.padding = 0.25,
  #    point.padding = 0.15,
  #    show.legend = FALSE
  #  ) +
  facet_wrap(~ Metric, scales = "free_y") +
  theme_minimal() +
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust = 1)) +
  labs(title = "Alpha diversity by species", x = NULL, y = NULL)

######



# ┌────────────────────────────────────────────────────────────────────────────┐
# │                    BETA-DIVERSITY - SAMPLES NOT MERGED                     │
# └────────────────────────────────────────────────────────────────────────────┘

# 1a - POST DECONTAMINATION: Bray Curtis ======================================

# ASV x samples matrix
matrix_ASVs_counts_reverse <- ASVs_counts_Pcla_Lsar_readytouse %>%
  column_to_rownames("ASVs") %>%   #set ASV IDs as rownames
  as.matrix()

# reorder matrix columns as in metadata
matrixdiversity <- matrix_ASVs_counts_reverse[, metadata_Pclav_Lsar_noLS98 $IIT_ID, drop = FALSE]

matrix_relativeabu_B_diversity <- sweep(matrixdiversity, 2, colSums(matrixdiversity), "/")          # ASV x samples (rel abund)
distance_bray_Bdiv  <- vegdist(t(matrix_relativeabu_B_diversity), method = "bray")      # distance between samples

pcoa_betadiv <- pcoa(distance_bray_Bdiv)
scores <- as.data.frame(pcoa_betadiv$vectors[, 1:2])
scores$IIT_ID <- rownames(scores)

scores <- scores %>%
  left_join(metadata_Pclav_Lsar_noLS98, by = "IIT_ID")

scores_plot <- scores %>%
  filter(!is.na(Axis.1), !is.na(Axis.2), !is.na(Specie))

#color palette that I like and that is clearly visible
pal_21 <- c("#66c2a5","#fc8d62","#8da0cb")

sp_levels <- sort(unique(scores_plot$Specie))
pal_named <- setNames(pal_21[seq_along(sp_levels)], sp_levels)

# % variance explained by the first axes
var_exp <- pcoa_betadiv$values$Relative_eig * 100

ggplot(scores_plot, aes(x = Axis.1, y = Axis.2)) +
  geom_point(aes(color = Specie), size = 5.0, alpha = 0.9) +
  scale_color_manual(values = pal_named) +
  theme_minimal() +
  geom_text_repel(
    aes(label = IIT_ID),
    size = 2.7,
    max.overlaps = 200,   # increase if it hides too many
    box.padding = 0.25,
    point.padding = 0.15,
    show.legend = FALSE
  ) +
  labs(
    title = "PCoA (Bray-Curtis) - post-decontaminazione 0.2",
    x = paste0("PCoA1 (", round(var_exp[1], 1), "%)"),
    y = paste0("PCoA2 (", round(var_exp[2], 1), "%)"),
    color = "Species"
  )


# 1b - POST DECONTAMINATION: Jaccard ==========================================

# ASV x samples matrix
matrix_ASVs_counts_reverse <- ASVs_counts_Pcla_Lsar_readytouse %>%
  column_to_rownames("ASVs") %>%   #set ASV IDs as rownames
  as.matrix()

# reorder matrix columns as in metadata
matrixdiversity <- matrix_ASVs_counts_reverse[, metadata_Pclav_Lsar_noLS98 $IIT_ID, drop = FALSE]

matrix_relativeabu_B_diversity <- sweep(matrixdiversity, 2, colSums(matrixdiversity), "/")          # ASV x samples (rel abund)
distance_bray_Bdiv  <- vegdist(t(matrix_relativeabu_B_diversity), method = "jaccard")      # distance between samples
?vegdist

pcoa_betadiv <- pcoa(distance_bray_Bdiv)
scores <- as.data.frame(pcoa_betadiv$vectors[, 1:2])
scores$IIT_ID <- rownames(scores)

scores <- scores %>%
  left_join(metadata_Pclav_Lsar_noLS98, by = "IIT_ID")

scores_plot <- scores %>%
  filter(!is.na(Axis.1), !is.na(Axis.2), !is.na(Specie))

#color palette that I like and that is clearly visible
pal_21 <- c("#66c2a5","#fc8d62","#8da0cb")

sp_levels <- sort(unique(scores_plot$Specie))
pal_named <- setNames(pal_21[seq_along(sp_levels)], sp_levels)

# % variance explained by the first axes
var_exp <- pcoa_betadiv$values$Relative_eig * 100

ggplot(scores_plot, aes(x = Axis.1, y = Axis.2)) +
  geom_point(aes(color = Specie), size = 5.0, alpha = 0.9) +
  scale_color_manual(values = pal_named) +
  theme_minimal() +
  #  geom_text_repel(
  #    aes(label = IIT_ID),
  #    size = 2.7,
  #    max.overlaps = 200,   # increase if it hides too many
  #   box.padding = 0.25,
  #    point.padding = 0.15,
  #    show.legend = FALSE
  #  ) +
  labs(
    title = "PCoA (Jaccard) - post-decontaminazione 0.2",
    x = paste0("PCoA1 (", round(var_exp[1], 1), "%)"),
    y = paste0("PCoA2 (", round(var_exp[2], 1), "%)"),
    color = "Species"
  )


# 2 - PRE DECONTAMINATION =====================================================

ASVs_counts_Pcla_Lsar_Bacteria_clean #ASVs_counts dataset

# ASV x samples matrix already created
matrix_ASVs_counts_reverse <- ASVs_counts_Pcla_Lsar_Bacteria_clean %>%
  column_to_rownames("ASVs") %>%   #set ASV IDs as rownames
  as.matrix()

# reorder matrix columns as in metadata
matrixdiversity <- matrix_ASVs_counts_reverse[, metadata_Pclav_Lsar_noLS98 $IIT_ID, drop = FALSE]


matrix_relativeabu_B_diversity <- sweep(matrixdiversity, 2, colSums(matrixdiversity), "/")          # ASV x samples (rel abund)
distance_bray_Bdiv  <- vegdist(t(matrix_relativeabu_B_diversity), method = "bray")      # distance between samples

pcoa_betadiv <- pcoa(distance_bray_Bdiv)
scores <- as.data.frame(pcoa_betadiv$vectors[, 1:2])
scores$IIT_ID <- rownames(scores)

scores <- scores %>%
  left_join(metadata_Pclav_Lsar_noLS98, by = "IIT_ID")

scores_plot <- scores %>%
  filter(!is.na(Axis.1), !is.na(Axis.2), !is.na(Specie))

ggplot(scores_plot, aes(x = Axis.1, y = Axis.2)) +
  geom_point(aes(color = Specie), size = 2, alpha = 0.9) +
  theme_minimal() +
  labs(title = "PCoA (Bray-Curtis)", x = "PCoA1", y = "PCoA2")

pal_21 <- c("#66c2a5","#fc8d62","#8da0cb")


sp_levels <- sort(unique(scores_plot$Specie))
pal_named <- setNames(pal_21[seq_along(sp_levels)], sp_levels)

#% variance explained by the first axes
var_exp <- pcoa_betadiv$values$Relative_eig * 100

ggplot(scores_plot, aes(x = Axis.1, y = Axis.2)) +
  geom_point(aes(color = Specie), size = 3, alpha = 0.9) +
  scale_color_manual(values = pal_named) +
  geom_text_repel(aes(label = Numero, color = Specie), size = 3, show.legend = FALSE) +
  theme_minimal() +
  labs(
    title = "PCoA (Bray-Curtis) - pre-decontaminazione",
    x = paste0("PCoA1 (", round(var_exp[1], 1), "%)"),
    y = paste0("PCoA2 (", round(var_exp[2], 1), "%)"),
    color = "Species"
  )

ggplot(scores_plot, aes(x = Axis.1, y = Axis.2)) +
  geom_point(aes(color = Specie), size = 3, alpha = 0.9) +
  scale_color_manual(values = pal_named) +
  geom_text_repel(
    aes(
      label = ifelse(Specie == "Paramuricea_clavata_(malata)", as.character(Numero), NA),
      color = Specie
    ),
    size = 3,
    show.legend = FALSE,
    na.rm = TRUE
  ) +
  theme_minimal() +
  labs(
    title = "PCoA (Bray-Curtis) - pre-decontaminazione",
    x = paste0("PCoA1 (", round(var_exp[1], 1), "%)"),
    y = paste0("PCoA2 (", round(var_exp[2], 1), "%)"),
    color = "Species"
  )

#species and number must be factors
scores_plot <- scores_plot %>%
  mutate(
    Specie = factor(Specie),
    Numero = factor(Numero)
  )

#species shape, number color 
ggplot(scores_plot, aes(x = Axis.1, y = Axis.2)) +
  geom_point(aes(shape = Specie, color = Numero), size = 3, alpha = 0.9) +
  theme_minimal() +
  labs(
    title = "PCoA (Bray-Curtis) - pre-decontaminazione",
    x = paste0("PCoA1 (", round(var_exp[1], 1), "%)"),
    y = paste0("PCoA2 (", round(var_exp[2], 1), "%)"),
    color = "Colony (Numero)",
    shape = "Species"
  )
#connection
ggplot(scores_plot, aes(x = Axis.1, y = Axis.2)) +
  geom_path(aes(group = Numero), alpha = 0.4) +           # connect replicates from the same colony
  geom_point(aes(color = Specie, shape = Specie), size = 3, alpha = 0.9) +
  scale_color_manual(values = pal_named) +
  theme_minimal() +
  labs(
    title = "PCoA (Bray-Curtis) - pre-decontaminazione",
    x = paste0("PCoA1 (", round(var_exp[1], 1), "%)"),
    y = paste0("PCoA2 (", round(var_exp[2], 1), "%)"),
    color = "Species",
    shape = "Species"
  )


#centroid plot
centroids_numero <- scores_plot %>%
  filter(!is.na(Numero), !is.na(Specie)) %>%
  group_by(Specie, Numero) %>%
  summarise(
    Axis.1 = mean(Axis.1, na.rm = TRUE),
    Axis.2 = mean(Axis.2, na.rm = TRUE),
    n_samples = n(),
    .groups = "drop"
  )

ggplot(centroids_numero, aes(x = Axis.1, y = Axis.2)) +
  geom_point(aes(color = Specie), size = 3.5, alpha = 0.95) +
  geom_text_repel(aes(label = Numero, color = Specie), size = 3, show.legend = FALSE) +
  scale_color_manual(values = pal_named) +
  theme_minimal() +
  labs(
    title = "PCoA (Bray-Curtis) - centroidi per colonia (Numero)",
    x = paste0("PCoA1 (", round(var_exp[1], 1), "%)"),
    y = paste0("PCoA2 (", round(var_exp[2], 1), "%)"),
    color = "Species"
  )
#####



# ┌────────────────────────────────────────────────────────────────────────────┐
# │          BETA-DIVERSITY - MERGED SAMPLES - POST-DECONTAMINATION             │
# └────────────────────────────────────────────────────────────────────────────┘

# 1 - MERGING ALL ==============================================================

#vector of samples present in the count table
samples_in_counts <- setdiff(colnames(ASVs_counts_Pcla_Lsar_readytouse), "ASVs")

#keep only metadata for samples that actually exist in the count table
map_samp2num <- metadata_Pclav_Lsar_noLS98 %>%
  filter(IIT_ID %in% samples_in_counts) %>%
  select(IIT_ID, Specie, Numero)

#merge by Numero
ASV_merged_byNumero <- ASVs_counts_Pcla_Lsar_readytouse %>%
  pivot_longer(
    cols = -ASVs,
    names_to = "IIT_ID",
    values_to = "Count"
  ) %>%
  left_join(map_samp2num, by = "IIT_ID") %>%
  group_by(ASVs, Specie, Numero) %>%
  summarise(Count = sum(Count, na.rm = TRUE), .groups = "drop") %>%
  mutate(Numero = as.character(Numero)) %>%
  pivot_wider(
    id_cols = c(ASVs),
    names_from = Numero,
    values_from = Count,
    values_fill = 0
  )



# 2 - B-DIVERSITY WITH MERGED SAMPLES ===========================================

metadata_byNumero <- map_samp2num %>%
  group_by(Numero) %>%
  summarise(
    Specie = first(Specie),
    n_samples_merged = n(),
    .groups = "drop"
  )

matrix_counts_num <- ASV_merged_byNumero %>%
  column_to_rownames("ASVs") %>%
  as.matrix()

matrix_rel_num <- sweep(matrix_counts_num, 2, colSums(matrix_counts_num), "/")
dist_bray_num  <- vegdist(t(matrix_rel_num), method = "bray")

pcoa_num <- pcoa(dist_bray_num)
scores_num <- as.data.frame(pcoa_num$vectors[,1:2])
scores_num$Numero <- rownames(scores_num)

scores_num <- scores_num %>%
  left_join(metadata_byNumero, by = "Numero")

var_exp_num <- pcoa_num$values$Relative_eig * 100

ggplot(scores_num, aes(Axis.1, Axis.2, color = Specie)) +
  geom_point(size = 3.5, alpha = 0.95) +
  scale_color_manual(values = pal_named) +
  theme_minimal() +
  labs(
    title = "PCoA (Bray-Curtis) - profili per colonia (Numero)",
    x = paste0("PCoA1 (", round(var_exp_num[1], 1), "%)"),
    y = paste0("PCoA2 (", round(var_exp_num[2], 1), "%)"),
    color = "Species"
  )

# 3 - MERGING ALL EXCEPT DISEASED P. CLAVATA ===================================

#samples present in the count table
samples_in_counts <- setdiff(colnames(ASVs_counts_Pcla_Lsar_readytouse), "ASVs")

#sample -> info map (only samples that exist in counts)
map_samp <- metadata_Pclav_Lsar_noLS98 %>%
  filter(IIT_ID %in% samples_in_counts) %>%
  select(IIT_ID, Specie, Numero)

stopifnot(!anyDuplicated(map_samp$IIT_ID))

#create merge_id: for "diseased" DO NOT merge (merge_id = IIT_ID), otherwise merge by Numero
map_samp <- map_samp %>%
  mutate(
    merge_id = if_else(Specie == "Paramuricea_clavata_(malata)",
                       IIT_ID,
                       as.character(Numero))
  )

#MERGE counts: sum reads per ASV within merge_id
ASVs_counts_merged_custom <- ASVs_counts_Pcla_Lsar_readytouse %>%
  pivot_longer(cols = -ASVs, names_to = "IIT_ID", values_to = "Count") %>%
  left_join(map_samp %>% select(IIT_ID, merge_id), by = "IIT_ID") %>%
  # if there are columns in counts that are NOT in metadata, keep them as single samples
  mutate(merge_id = if_else(is.na(merge_id), IIT_ID, merge_id)) %>%
  group_by(ASVs, merge_id) %>%
  summarise(Count = sum(Count, na.rm = TRUE), .groups = "drop") %>%
  pivot_wider(names_from = merge_id, values_from = Count, values_fill = 0)

#metadata collapsed by merge_id
metadata_merged_custom <- map_samp %>%
  group_by(merge_id) %>%
  summarise(
    Specie = first(Specie),
    Numero = first(Numero),
    n_samples_merged = n(),
    merged_rule = if_else(first(Specie) == "Paramuricea_clavata_(malata)", "kept_single", "merged_by_Numero"),
    .groups = "drop"
  )

#reorder columns as in metadata_merged_custom
ASVs_counts_merged_custom <- ASVs_counts_merged_custom %>%
  select(ASVs, all_of(metadata_merged_custom$merge_id))


# 4 - B-DIVERSITY WITH MERGED SAMPLES EXCEPT DISEASED P. CLAVATA ================

#ASV x unit (merge_id)
mat_counts <- ASVs_counts_merged_custom %>%
  column_to_rownames("ASVs") %>%
  as.matrix()

mat_rel <- sweep(mat_counts, 2, colSums(mat_counts), "/")   # ASV x unit (rel abund)
dist_bray <- vegdist(t(mat_rel), method = "bray")           # distances between units

pcoa_res <- pcoa(dist_bray)

#% explained variance
var_exp <- pcoa_res$values$Relative_eig * 100

scores <- as.data.frame(pcoa_res$vectors[, 1:2])
scores$merge_id <- rownames(scores)

scores_plot <- scores %>%
  left_join(metadata_merged_custom, by = "merge_id") %>%
  filter(!is.na(Axis.1), !is.na(Axis.2), !is.na(Specie))

pal_21 <- c("#66c2a5","#fc8d62","#8da0cb")
sp_levels <- sort(unique(scores_plot$Specie))
pal_named <- setNames(pal_21[seq_along(sp_levels)], sp_levels)

ggplot(scores_plot, aes(x = Axis.1, y = Axis.2)) +
  geom_point(aes(color = Specie), size = 3, alpha = 0.9) +
  scale_color_manual(values = pal_named) +
  theme_minimal() +
  labs(
    title = "PCoA (Bray-Curtis) - merged by colony, excluding P.clavata malata",
    x = paste0("PCoA1 (", round(var_exp[1], 1), "%)"),
    y = paste0("PCoA2 (", round(var_exp[2], 1), "%)"),
    color = "Species"
  )

# 5 - B-DIVERSITY WITH SPLIT NON-MERGED SAMPLES =================================

#IMPORT DATASET -----> for the future me: REDO ASV-counts if needed, and also the metadata tables!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

only_LS_ASVs_counts <- read_excel("only_LS_ASVs_counts.xlsx")
only_PC_ASVs_counts <- read_excel("only_PC_ASVs_counts.xlsx")
only_LS_metadata <- read_excel("only_LS_metadata.xlsx", 
                               col_types = c("text", "text", "text", 
                                             "text", "text", "date", "numeric", 
                                             "numeric", "text", "numeric", "text", 
                                             "text", "numeric"))
only_PC_metadata <- read_excel("only_PC_metadata.xlsx", 
                               col_types = c("text", "text", "text", 
                                             "text", "text", "date", "numeric", 
                                             "numeric", "text", "numeric", "text", 
                                             "text", "numeric"))
only_PC_metadata <- only_PC_metadata %>%
  rename(IIT_ID = `IIT-ID`)

# ASV x samples matrix
matrix_ASVs_counts_reverse <- only_PC_ASVs_counts %>%
  column_to_rownames("ASVs") %>%   #set ASV IDs as rownames
  as.matrix()

# reorder matrix columns as in metadata
matrixdiversity <- matrix_ASVs_counts_reverse[, only_PC_metadata $IIT_ID, drop = FALSE]

matrix_relativeabu_B_diversity <- sweep(matrixdiversity, 2, colSums(matrixdiversity), "/")          # ASV x samples (rel abund)
distance_bray_Bdiv  <- vegdist(t(matrix_relativeabu_B_diversity), method = "bray")      # distance between samples
?vegdist

pcoa_betadiv <- pcoa(distance_bray_Bdiv)
scores <- as.data.frame(pcoa_betadiv$vectors[, 1:2])
scores$IIT_ID <- rownames(scores)

scores <- scores %>%
  left_join(only_PC_metadata, by = "IIT_ID")

scores_plot <- scores_plot %>%
  mutate(Numero = as.factor(Numero))

#color palette that I like and that is clearly visible
pal_14 <- c(
  "#E41A1C", "#377EB8", "#4DAF4A", "#984EA3",
  "#FF7F00", "#FFFF33", "#A65628", "#F781BF",
  "#999999", "#66C2A5", "#FC8D62", "#8DA0CB",
  "#E78AC3", "#A6D854"
)


num_levels <- sort(unique(scores_plot$Numero))
pal_named_num <- setNames(pal_14[seq_along(num_levels)], num_levels)

# % variance explained by the first axes
var_exp <- pcoa_betadiv$values$Relative_eig * 100

ggplot(scores_plot, aes(x = Axis.1, y = Axis.2)) +
  geom_path(aes(group = Numero), alpha = 0.4) +
  geom_point(aes(color = Numero), size = 3, alpha = 0.9) +
  scale_color_manual(values = pal_named_num) +
  theme_minimal() +
  labs(
    title = "PCoA (Jaccard) - Paramuricea clavata",
    x = paste0("PCoA1 (", round(var_exp[1], 1), "%)"),
    y = paste0("PCoA2 (", round(var_exp[2], 1), "%)"),
    color = "Numero"
  )
