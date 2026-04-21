setwd(file.path(Sys.getenv("THESIS_DIR"), "code/05_GO_annotation"))
load("03_map_xtropproteins_xtrop_genes.RData")

library(biomaRt)

## --------------------------------------------
## 1. Retrieve GO annotations (ENSXETG level)
## --------------------------------------------

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

## --------------------------------------------
## 2. Build GO mapping tables ENSXETG → GO (TERM2GENE)
## --------------------------------------------

TERM2GENE <- unique(go_annot[, c("go_id", "ensembl_gene_id")])
colnames(TERM2GENE) <- c("GO", "GENE")

TERM2NAME <- unique(go_annot[, c("go_id", "name_1006")])
colnames(TERM2NAME) <- c("GO", "NAME")

## Extended mapping with ontology
TERM2GENE_full <- unique(go_annot[, c(
    "go_id",
    "ensembl_gene_id",
    "namespace_1003"
)])
colnames(TERM2GENE_full) <- c("GO", "GENE", "ONTOLOGY")

TERM2NAME_full <- TERM2NAME

## --------------------------------------------
## 3. Remove version suffix from Xenopus gene IDs
## --------------------------------------------

loc_to_xtrop_gene$xtrop_gene_id_clean <- sub("\\..*$", "", loc_to_xtrop_gene$xtrop_gene_id)

# Check that the formats now match
head(loc_to_xtrop_gene$xtrop_gene_id)
head(loc_to_xtrop_gene$xtrop_gene_id_clean)
head(TERM2GENE$GENE)

# Overlap check
sum(loc_to_xtrop_gene$xtrop_gene_id_clean %in% TERM2GENE$GENE)

## --------------------------------------------
## 4. Map LOC_* → GO (pooled)
## --------------------------------------------
loc_to_go <- merge(
    loc_to_xtrop_gene,
    TERM2GENE,
    by.x = "xtrop_gene_id_clean",
    by.y = "GENE"
)

# Result LOC_* → GO

## --------------------------------------------
## 5. Map LOC_* → GO with ontology
## --------------------------------------------

loc_to_go_full <- merge(
    loc_to_xtrop_gene,
    TERM2GENE_full,
    by.x = "xtrop_gene_id_clean",
    by.y = "GENE"
)

## --------------------------------------------
## 6. Split by ontology
## --------------------------------------------

loc_to_go_BP <- loc_to_go_full[loc_to_go_full$ONTOLOGY == "biological_process", ]
loc_to_go_MF <- loc_to_go_full[loc_to_go_full$ONTOLOGY == "molecular_function", ]
loc_to_go_CC <- loc_to_go_full[loc_to_go_full$ONTOLOGY == "cellular_component", ]

## --------------------------------------------
## 7. Save
## --------------------------------------------
write.csv(loc_to_go_full, "loc_to_go_full.csv", row.names = FALSE)
write.csv(loc_to_go_BP, "loc_to_go_BP.csv", row.names = FALSE)
write.csv(loc_to_go_MF, "loc_to_go_MF.csv", row.names = FALSE)
write.csv(loc_to_go_CC, "loc_to_go_CC.csv", row.names = FALSE)


# 6. Save
save.image("04_GO_terms_mapping.RData")
