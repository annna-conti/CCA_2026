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

