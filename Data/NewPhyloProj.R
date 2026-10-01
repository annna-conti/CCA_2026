library(phyloseq)
ps_raw <- readRDS(
  "annabelle_16S_emu_taxon_count_phyloseq.rds"
)

ps_raw

counts_raw <- as( 
  otu_table(ps_raw),
  "matrix"
  )

if (!taxa_are_rows(otu_table(ps_raw))) {
  counts_raw <- t(counts_raw)
}

tax_raw <- as( 
  tax_table(ps_raw),
  "matrix"
  )

counts_raw <- t(counts_raw)

all(
  rownames(counts_raw %in% rownames(tax_raw)
)

tax_raw <- tax_raw[
  rownames(counts_raw),
  ,
  drop = FALSE
]

identical(
  rownames(counts_raw),
  rownames(tax_raw)
)

####################################
#cleaning taxonomy#
####################################

identical(
  rownames(counts_raw),
  rownames(tax_raw)
)

is_bacteria <- 
  !is.na(tax_raw[, "superkingdom"]) &
  trimws(tax_raw[, "superkingdom"]) == "Bacteria"

sum(!is_bacteria)

has_phylum <- 
  !is.na(tax_raw[, "phylum"]) &
  trimws(tax_raw[, "phylum"]) != ""

sum(!has_phylum)

is_propionibacteriaceae <- 
  !is.na(tax_raw[, "family"]) &
  tax_raw[, "family"] == "Propionibacteriaceae"

sum(is_propionibacteriaceae)


is_organelle <- apply(
  tax_raw,
  1,
  function(x) {
    any(
      grepl(
        "Mitochondria|Chloroplast",
        x,
        ignore.case = TRUE
      )
    )
  }
)

data.frame(
  category = c(
    "Non-bacteria",
    "No phylum",
    "Propionibacteriaceae",
    "Mitochondria/chloroplast"
  ),
  n_taxa = c(
    sum(!is_bacteria),
    sum(!has_phylum),
    sum(is_propionibacteriaceae),
    sum(is_organelle)
  )
)
#only need to remove 3 taxa total 

keep_taxa <- 
  is_bacteria &
  has_phylum &
  !is_propionibacteriaceae &
  !is_organelle

counts_taxclean <- counts_raw[
  keep_taxa,
  ,
  drop = FALSE
]

tax_clean <- tax_raw[
  keep_taxa,
  ,
  drop = FALSE
]

identical(
  rownames(counts_taxclean),
  rownames(tax_clean)
)

taxon_totals <- rowSums(counts_taxclean)

head(
  sort(taxon_totals),
  20
)

sum(taxon_totals == 1)
sum(taxon_totals <= 1)

####trying to get rid of singletons 
singleton_like <- taxon_totals <= (1 + 1e-6)
sum(singleton_like)
taxon_totals[singleton_like]

counts_clean <- counts_taxclean[
  !singleton_like,
  ,
  drop = FALSE
]

tax_clean <- tax_clean[
  !singleton_like,
  ,
  drop = FALSE
]
tax_clean <- tax_clean[
  rownames(counts_clean),
  ,
  drop = FALSE
]

identical(
  rownames(counts_clean),
  rownames(tax_clean)
)
