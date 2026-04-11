setwd(file.path(Sys.getenv("THESIS_DIR"), "code/05_GO_annotation"))


# Load BLAST results
blast_best <- read.delim(
    "ensembl_blast_best.tsv",
    header = FALSE,
    sep = "\t",
    stringsAsFactors = FALSE
)

colnames(blast_best) <- c(
    "protein_id", "xtrop_protein_id",
    "pident", "length", "evalue", "bitscore"
)

# Load Xenopus tropicalis protein to gene mapping
xtrop_map <- read.delim(
    "xtrop_protein_to_gene.tsv",
    header = FALSE,
    sep = "\t",
    stringsAsFactors = FALSE
)

colnames(xtrop_map) <- c(
    "xtrop_protein_id",
    "xtrop_gene_id"
)

# Map protein IDs to Xenopus tropicalis gene IDs
protein_to_xtrop_gene <- merge(
    blast_best,
    xtrop_map,
    by = "xtrop_protein_id"
)

loc_to_xtrop_gene <- merge(
    map_loc_protein,
    protein_to_xtrop_gene,
    by = "protein_id"
)


save.image("03_map_xtropproteins_xtrop_genes.RData")
