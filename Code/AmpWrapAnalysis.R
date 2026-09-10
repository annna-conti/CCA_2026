library(phyloseq)
physeq <- readRDS("Data/Annabelle_phyloseq.rds")
physeq

#this makes sure that there will always be an original untouched copy 
physeq_raw <- physeq 

#to extract the ASV count matrix --> OTU = operational taxonomic unit, to classify groups of closely related microbes based on sequence similarity, can use as a 
#proxy for species 
ASV_counts <- as(otu_table(physeq), "matrix")

#ASVs as rows and samples as columns 
if(!taxa_are_rows(physeq)) {
  ASV_counts <- t(ASV_counts)
}

dim(ASV_counts) #this will tell you the number of rows and columns 
ASV_counts[1:5, 1:5]

#extracting taxonomy 
ASV_taxonomy <- as(tax_table(physeq), "matrix")

dim(ASV_taxonomy)
head(ASV_taxonomy)

colnames(ASV_taxonomy)

#finding anything labeled mitochondria in taxonomy
mito <- apply(
  ASV_taxonomy,
  1,
  function(x) any(grepl("Mitochondria",x, ignore.case = TRUE))
)

#finding anything labeled chloroplast
chloro <- apply(
  ASV_taxonomy, 
  1, 
  function(x) any(grepl("Chloroplast", x, ignore.case = TRUE))
)

sum(mito)
sum(chloro)

ASV_taxonomy[mito, , drop = FALSE]

ASV_taxonomy[chloro, , drop = FALSE]

sum(mito)
sum(chloro)

#combining mito and chloro to take them out 
organelle <- mito | chloro
sum(organelle)

#removing mito and chloro from ASV count table
ASV_counts_clean <- ASV_counts[!organelle, , drop = FALSE]
#removing mito and chloro from taxonomy table 
ASV_taxonomy_clean <- ASV_taxonomy[!organelle, , drop = FALSE]

dim(ASV_counts_clean)
dim(ASV_taxonomy_clean)

identical(
  rownames(ASV_counts_clean),
  rownames(ASV_taxonomy_clean)
)

#checking how many reads were removed 
#before removal:
sum(ASV_counts)
#after removal: 
sum(ASV_counts_clean)
#number of mito and chloro reads removed 
sum(ASV_counts) - sum(ASV_counts_clean)

#removing any ASVs with zero reads and singletons
ASV_counts_clean
ASV_taxonomy_clean

ASV_totals <- rowSums(ASV_counts_clean)
sum(ASV_totals ==0)
#there are no unused taxa so this returns 0 

dim(ASV_counts_clean)
dim(ASV_taxonomy_clean)

identical(
  rownames(ASV_counts_clean),
  rownames(ASV_taxonomy_clean)
)

ASV_totals <- rowSums(ASV_counts_clean)
singleton <- ASV_totals == 1
sum(singleton)

ASV_counts_clean[singleton, , drop = FALSE]
ASV_taxonomy_clean[singleton, , drop = FALSE]
#logical vector singleton did not select any rows 
sum(singleton)
#this confirms that there are no singleton ASVs in your current table 

keep_nonsingleton <- rowSums(ASV_counts_clean) > 1

ASV_counts_clean <- ASV_counts_clean[
  keep_nonsingleton,
  ,
  drop = FALSE
]

ASV_taxonomy_clean <- ASV_taxonomy_clean[
  keep_nonsingleton,
  ,
  drop = FALSE
]

min(rowSums(ASV_counts_clean))

#the outputs of these tell me that my phyloseq abundance table contains relative abundances, not read counts 
colSums(ASV_counts)

summary(colSums(ASV_counts))

sum(abs(ASV_counts - round(ASV_counts)) > 1e-8)

head(sort(ASV_counts[ASV_counts > 0]), 20)

# Re-create the cleaned tables starting from the originals
# Remove ONLY mitochondria and chloroplasts for now

organelle <- mito | chloro

ASV_counts_clean <- ASV_counts[
  !organelle,
  ,
  drop = FALSE
]

ASV_taxonomy_clean <- ASV_taxonomy[
  !organelle,
  ,
  drop = FALSE
]

dim(ASV_counts_clean)
dim(ASV_taxonomy_clean)

identical(
  rownames(ASV_counts_clean),
  rownames(ASV_taxonomy_clean)
)
#Alessia's code for singletons doesn't work for me because my table doesn't contain counts 
#Skip singleton stuff 

#Remove anything without assigned phylum for now instead
no_phylum <- is.na(ASV_taxonomy_clean[, "phylum"]) |
  trimws(ASV_taxonomy_clean[, "phylum"]) == ""
sum(no_phylum)
#this returned a 0 so there are no ASVs that are missing a phylum assignment 

#getting rid of Propionibacteriaceae
propionibacteriaceae <-
  !is.na(ASV_taxonomy_clean[, "family"]) &
  ASV_taxonomy_clean[, "family"] == "Propionibacteriaceae"

sum(propionibacteriaceae)

ASV_taxonomy_clean[
  propionibacteriaceae,
  ,
  drop = FALSE
]

ASV_taxonomy_clean[
  propionibacteriaceae,
  ,
  drop = FALSE
]

ASV_taxonomy_clean <- ASV_taxonomy_clean[
  !propionibacteriaceae,
  ,
  drop = FALSE
]

dim(ASV_counts_clean)
dim(ASV_taxonomy_clean)

identical(
  rownames(ASV_counts_clean),
  rownames(ASV_taxonomy_clean)
)
#the tables are not identical anymore so need to see what is mismatched
setdiff(
  rownames(ASV_counts_clean),
  rownames(ASV_taxonomy_clean)
)

setdiff(
  rownames(ASV_taxonomy_clean),
  rownames(ASV_counts_clean)
)

# REBUILD CLEAN TABLES USING OTU IDs

mito_ids <- rownames(ASV_taxonomy)[
  apply(
    ASV_taxonomy,
    1,
    function(x) any(grepl("Mitochondria", x, ignore.case = TRUE))
  )
]

chloro_ids <- rownames(ASV_taxonomy)[
  apply(
    ASV_taxonomy,
    1,
    function(x) any(grepl("Chloroplast", x, ignore.case = TRUE))
  )
]

prop_ids <- rownames(ASV_taxonomy)[
  !is.na(ASV_taxonomy[, "family"]) &
    ASV_taxonomy[, "family"] == "Propionibacteriaceae"
]

# Combine all OTUs to remove
remove_ids <- unique(
  c(mito_ids, chloro_ids, prop_ids)
)

remove_ids
length(remove_ids)

ASV_counts_clean <- ASV_counts[
  !(rownames(ASV_counts) %in% remove_ids),
  ,
  drop = FALSE
]

ASV_taxonomy_clean <- ASV_taxonomy[
  !(rownames(ASV_taxonomy) %in% remove_ids),
  ,
  drop = FALSE
]

ASV_taxonomy_clean <- ASV_taxonomy_clean[
  rownames(ASV_counts_clean),
  ,
  drop = FALSE
]

#CHECK
dim(ASV_counts_clean)
dim(ASV_taxonomy_clean)

identical(
  rownames(ASV_counts_clean),
  rownames(ASV_taxonomy_clean)
)

