setwd(file.path(Sys.getenv("THESIS_DIR"), "code/05_GO_annotation"))
load("03_map_xtropproteins_xtrop_genes.RData")

# 1. Retrieve GO annotations for ENSXETG

library(biomaRt)

# Connect to Ensembl
mart <- useEnsembl(
    biomart = "genes",
    dataset = "xtropicalis_gene_ensembl"
)

# Retrieve GO annotations
go_annot <- getBM(
    attributes = c(
        "ensembl_gene_id", # ENSXETG
        "go_id",
        "namespace_1003", # BP / MF / CC
        "name_1006" # GO term name
    ),
    mart = mart
)

# Remove entries without GO
go_annot <- go_annot[go_annot$go_id != "", ]

# 2. Build ENSXETG → GO mapping (TERM2GENE)
TERM2GENE <- unique(go_annot[, c("go_id", "ensembl_gene_id")])
colnames(TERM2GENE) <- c("GO", "GENE")

TERM2NAME <- unique(go_annot[, c("go_id", "name_1006")])
colnames(TERM2NAME) <- c("GO", "NAME")

# 3. Remove version suffix from Xenopus gene IDs
loc_to_xtrop_gene$xtrop_gene_id_clean <- sub("\\..*$", "", loc_to_xtrop_gene$xtrop_gene_id)

# Check that the formats now match
head(loc_to_xtrop_gene$xtrop_gene_id)
head(loc_to_xtrop_gene$xtrop_gene_id_clean)
head(TERM2GENE$GENE)

# 4. Check overlap
sum(loc_to_xtrop_gene$xtrop_gene_id_clean %in% TERM2GENE$GENE)

# 5. Map LOC_* → GO
loc_to_go <- merge(
    loc_to_xtrop_gene,
    TERM2GENE,
    by.x = "xtrop_gene_id_clean",
    by.y = "GENE"
)

# Final mapping
# LOC_* → GO


# 6. Save
save.image("04_GO_terms_mapping.RData")

# Save loc_to_go
write.csv(loc_to_go, "loc_to_go.csv", row.names = FALSE)
