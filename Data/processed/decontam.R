if (!requireNamespace("BiocManager", quietly = TRUE)) {
  install.packages("BiocManager")
}

BiocManager::install("decontam")

library(phyloseq)
library(decontam)

packageVersion("decontam")

ps_clean <- readRDS(
  "processed/annabelle_16S_taxclean_predecontam_phyloseq.rds"
)

sample_names(ps_clean)
"barcode18" %in% sample_names(ps_clean)
"barcode19" %in% sample_names(ps_clean)

metadata <- data.frame(
  SampleID = sample_names(ps_clean),
  is_negative = FALSE
)

rownames(metadata) <- metadata$SampleID

# Mark the negative control
metadata["barcode18", "is_negative"] <- TRUE
sample_data(ps_clean) <- sample_data(metadata)

###decontam###
contam_prev <- isContaminant(
  ps_clean,
  method = "prevalence",
  neg = "is_negative",
  threshold = 0.1
)

table(contam_prev$contaminant)
sum(contam_prev$contaminant)
#this is saying nothing was flagged as a contaminant - prevalence looks for features 
#that are more prevalent in negative controls than in true samples --> since I only have 1 negative there's not very much
#replication and maybe this is ok?

seqtab <- as(otu_table(ps_clean), "matrix")

if(taxa_are_rows(ps_clean)) {
  seqtab <- t(seqtab)
}

dim(seqtab)
head(rownames(seqtab))
head(colnames(seqtab))
##------------------
  
library(phyloseq)
library(decontam)
library(readr)

ps_clean <- readRDS(
  "annabelle_16S_taxclean_predecontam_phyloseq.rds"
)

metadata <- read.csv(
  "CCA Project Metadata - Sheet1.csv",
  check.names = FALSE
)

head(metadata)
str(metadata)

ps_clean <- phyloseq::prune_samples(
  phyloseq::sample_names(ps_clean) %in% metadata$SampleID,
  ps_clean
)

metadata <- metadata[
  match(phyloseq::sample_names(ps_clean), metadata$SampleID),
  ,
  drop = FALSE
]

identical(
  metadata$SampleID,
  phyloseq::sample_names(ps_clean)
)

metadata$is_negative <- metadata$ControlType == "negative"
rownames(metadata) <- metadata$SampleID
sample_data(ps_clean) <- sample_data(metadata)

sample_data(ps_clean)[
  c("barcode18", "barcode19"),
  c("SampleType", "Quantification", "ControlType", "is_negative")
]

##ALESSIA DECONTAM WORKFLOW##
#going to use different thresholds and then see how many contaminants those each flag to decide which threshold is aggressive/mild enough for my data

# Extract abundance/count table
seqtab <- as(
  phyloseq::otu_table(ps_clean),
  "matrix"
)

# decontam wants samples as rows, taxa as columns
if (phyloseq::taxa_are_rows(ps_clean)) {
  seqtab <- t(seqtab)
}

conc <- as.numeric(
  phyloseq::sample_data(ps_clean)$Quantification
)

names(conc) <- phyloseq::sample_names(ps_clean)

identical(
  names(conc),
  rownames(seqtab)
)

any(is.na(conc))
any(conc <= 0)



contam_freq_01 <- isContaminant(
  seqtab,
  method = "frequency",
  conc = conc,
  threshold = 0.1,
  normalize = TRUE,
  detailed = TRUE
)
#3.1% flagged

contam_freq_02 <- isContaminant(
  seqtab,
  method = "frequency",
  conc = conc,
  threshold = 0.2,
  normalize = TRUE,
  detailed = TRUE
)
#7.8% flagged

contam_freq_05 <- isContaminant(
  seqtab,
  method = "frequency",
  conc = conc,
  threshold = 0.5,
  normalize = TRUE,
  detailed = TRUE
)
#24.1% flagged

#threshold default in decontam is 0.1 --> will make score histogram that tells me if the data is split any which way and could also choose threshold based on that 

hist(
  contam_freq_01$p,
  breaks = 100,
  main = "decontam frequency scores: 0–0.5",
  xlab = "decontam score",
  xlim = c(0, 0.5)
)

abline(v = 0.1, lty = 2)
abline(v = 0.2, lty = 2)

#my project is interested in microbial differences so real low-abundnace taxa could matter
#plot does not show a clear trough separation that would justify moving the cutoff up
#will use 0.1 as threshold because it is conservative --> 0.2 would remove an additional 56 taxa for kind of no reason
#the frequency-based decontam method was applied using post-PCR
#DNA concentration. the default threshold of 0.1 was retained because inspection of the score distribution
#did not reveal a clear data-driven breakpoint supporting higher threshold

neg <- as.logical(
  phyloseq::sample_data(ps_clean)$is_negative
)

names(neg) <- phyloseq::sample_names(ps_clean)

table(neg)
identical(
  names(neg),
  rownames(seqtab)
)

contam_prev_01 <- decontam::isContaminant(
  seqtab,
  method = "prevalence",
  neg = neg,
  threshold = 0.1,
  normalize = TRUE,
  detailed = TRUE
)
table(contam_prev_01$contaminant)

contam_combined_01 <- decontam::isContaminant(
  seqtab,
  method = "combined",
  conc = conc,
  neg = neg,
  threshold = 0.1,
  normalize = TRUE,
  detailed = TRUE
)

table(contam_combined_01$contaminant)


freq_ids <- rownames(contam_freq_01)[
  contam_freq_01$contaminant
]

prev_ids <- rownames(contam_prev_01)[
  contam_prev_01$contaminant
]

combined_ids <- rownames(contam_combined_01)[
  contam_combined_01$contaminant
]

length(freq_ids)
length(prev_ids)
length(combined_ids)

comparison <- data.frame(
  Taxon = rownames(contam_freq_01),
  Frequency = contam_freq_01$contaminant,
  Prevalence = contam_prev_01$contaminant,
  Combined = contam_combined_01$contaminant
)

table(
  Frequency = comparison$Frequency,
  Combined = comparison$Combined
)

# Flagged by both frequency and combined
both_ids <- comparison$Taxon[
  comparison$Frequency &
    comparison$Combined
]

# Frequency-only calls
frequency_only_ids <- comparison$Taxon[
  comparison$Frequency &
    !comparison$Combined
]

# Combined-only calls
combined_only_ids <- comparison$Taxon[
  !comparison$Frequency &
    comparison$Combined
]

length(both_ids)
length(frequency_only_ids)
length(combined_only_ids)


neg_counts <- seqtab["barcode18", ]

neg_present <- sort(
  neg_counts[neg_counts > 0],
  decreasing = TRUE
)

length(neg_present)

head(neg_present, 30)


sum(
  freq_ids %in% names(neg_present)
)

table(contam_prev_01$contaminant)

table(contam_combined_01$contaminant)
#only 3 taxa were really convincing under the combined method 
table(
  Frequency = comparison$Frequency,
  Combined = comparison$Combined
)
#1136 taxa are not contaminants under either prevalence or frequency 
#33 taxa are contaminants to frequency, but combined test doesn't have enough support 
#3 taxa are contaminants for frequency AND prevalence 
#the combined test didn't flag anything that the frequency method didn't already flag 
#combined method combines the two probabilities using Fisher's method --> feature can have strong frequency evidence and some supportive prevalence evidence without independently crossing the prevalence threshold

sum(freq_ids %in% names(neg_present))
#only have one negative control so it's hard to use the prevalence method, which is dependent on the negative control to find contamination
#of the 36 taxa identified as contaminants by the frequency method, 3 were actually detected in my negative control


##making a separate object just for fitting decontam:
ps_decontam_fit <- phyloseq::prune_samples(
  phyloseq::sample_names(ps_clean) != "barcode19",
  ps_clean
)
phyloseq::nsamples(ps_decontam_fit)
phyloseq::sample_names(ps_decontam_fit)
table(
  phyloseq::sample_data(ps_decontam_fit)$is_negative
)

#rerunning frequency correctly - wihtout barcode 18
contam_freq_01_final <- decontam::isContaminant(
  ps_decontam_fit,
  method = "frequency",
  conc = "Quantification",
  neg = "is_negative",
  threshold = 0.1,
  normalize = TRUE,
  detailed = TRUE
)
table(contam_freq_01_final$contaminant)

#prevalence:
contam_prev_01_final <- decontam::isContaminant(
  ps_decontam_fit,
  method = "prevalence",
  neg = "is_negative",
  threshold = 0.1,
  detailed = TRUE
)

#combined:
contam_combined_01_final <- decontam::isContaminant(
  ps_decontam_fit,
  method = "combined",
  conc = "Quantification",
  neg = "is_negative",
  threshold = 0.1,
  normalize = TRUE,
  detailed = TRUE
)

#checking:
table(contam_prev_01_final$contaminant)

table(contam_combined_01_final$contaminant)

table(contam_freq_01_final$contaminant)

table(
  Frequency = comparison$Frequency,
  Combined = comparison$Combined
)

sum(freq_ids %in% names(neg_present))
#prevalence still didn't flag any contamination 

# Taxa identified by the FINAL frequency analysis
freq_ids_final <- rownames(contam_freq_01_final)[
  contam_freq_01_final$contaminant
]

# Taxa identified by the FINAL prevalence analysis
prev_ids_final <- rownames(contam_prev_01_final)[
  contam_prev_01_final$contaminant
]

# Taxa identified by the FINAL combined analysis
combined_ids_final <- rownames(contam_combined_01_final)[
  contam_combined_01_final$contaminant
]
length(freq_ids_final)
length(prev_ids_final)
length(combined_ids_final)


#new comparison
comparison_final <- data.frame(
  Taxon = rownames(contam_freq_01_final),
  Frequency = contam_freq_01_final$contaminant,
  Prevalence = contam_prev_01_final$contaminant,
  Combined = contam_combined_01_final$contaminant
)
table(
  Frequency = comparison_final$Frequency,
  Combined = comparison_final$Combined
)
all(
  combined_ids_final %in% freq_ids_final
)#the three taxa are still contained in the frequency calls 

table(contam_freq_01_final$contaminant)

table(
  Frequency = comparison_final$Frequency,
  Combined = comparison_final$Combined
)

sum(
  freq_ids_final %in% names(neg_present)
)
#so with this CORRECTED COMBINATION (because code above for the combo wasn't right), 
#29 flagged by frequency 
#3 flagged by frequency and strong enough to be in combined 
#0 flagged by combined but not frequency 
#The frequency method is my primary method because I have DNA concentration measurements across the biological samples. 
#I retained decontam's default threshold of 0.1 because my score histogram did not show a clear breakpoint supporting a higher cutoff. 
#My single negative control was used as supporting evidence, 
#but prevalence-based detection was limited by the absence of replicate negative controls.

#want to save this before deleting contaminants:
# Taxonomy
tax_df <- as.data.frame(
  as(phyloseq::tax_table(ps_clean), "matrix")
)

# Negative-control abundances
neg_abundance <- seqtab["barcode18", ]

# Final contaminant table
contaminant_table_final <- data.frame(
  Taxon = freq_ids_final,
  FrequencyScore = contam_freq_01_final[
    freq_ids_final,
    "p"
  ],
  CombinedFlag = freq_ids_final %in% combined_ids_final,
  NegativeControlAbundance = neg_abundance[freq_ids_final],
  tax_df[freq_ids_final, , drop = FALSE],
  row.names = NULL
)

contaminant_table_final

write.csv(
  contaminant_table_final,
  "../decontam_final_contaminants_threshold_0.1.csv",
  row.names = FALSE
)

####remove those 32 taxa:
ps_decontam <- phyloseq::prune_taxa(
  !(phyloseq::taxa_names(ps_clean) %in% freq_ids_final),
  ps_clean
)
phyloseq::ntaxa(ps_clean)
phyloseq::ntaxa(ps_decontam)

saveRDS(
  ps_decontam,
  "annabelle_16S_final_decontam_phyloseq.rds"
)


# ============================================================
# CHECKPOINT: BEFORE VS AFTER DECONTAM
# ============================================================

qc_before_after <- data.frame(
  SampleID = phyloseq::sample_names(ps_clean),
  
  Counts_Before = phyloseq::sample_sums(ps_clean),
  
  Counts_After = phyloseq::sample_sums(ps_decontam)
)

# Extract abundance matrices
mat_before <- as(
  phyloseq::otu_table(ps_clean),
  "matrix"
)

mat_after <- as(
  phyloseq::otu_table(ps_decontam),
  "matrix"
)

# Put samples in rows if necessary
if (phyloseq::taxa_are_rows(ps_clean)) {
  mat_before <- t(mat_before)
}

if (phyloseq::taxa_are_rows(ps_decontam)) {
  mat_after <- t(mat_after)
}

qc_before_after$Taxa_Before <- rowSums(mat_before > 0)

qc_before_after$Taxa_After <- rowSums(mat_after > 0)

qc_before_after$Counts_Removed <-
  qc_before_after$Counts_Before -
  qc_before_after$Counts_After

qc_before_after$Taxa_Removed <-
  qc_before_after$Taxa_Before -
  qc_before_after$Taxa_After

qc_before_after$Percent_Counts_Removed <-
  100 *
  qc_before_after$Counts_Removed /
  qc_before_after$Counts_Before

qc_before_after

qc_before_after[
  order(
    qc_before_after$Percent_Counts_Removed,
    decreasing = TRUE
  ),
]

write.csv(
  qc_before_after,
  "../decontam_before_after_QC.csv",
  row.names = FALSE
)


##making a separate analysis object that doesn't have the controls for ecological analysis
ps_analysis <- phyloseq::prune_samples(
  phyloseq::sample_data(ps_decontam)$ControlType == "none",
  ps_decontam
)

saveRDS(
  ps_analysis,
  "../annabelle_16S_analysis_no_controls_phyloseq.rds"
)

# ============================================================
# MAKING THE RAREFACTION CURVES
# ============================================================
qc_before_after


