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
