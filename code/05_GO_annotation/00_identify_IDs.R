# setwd(file.path(Sys.getenv("THESIS_DIR")))
# load("code/04_de/resultsInteraction.RData")
# load("code/04_de/resultsTemperatureEffects.RData")
# gtf <- read.delim("data/metadata/moor_frog_annotation.gtf", header = FALSE, sep = "\t", comment.char = "#")

library(stringr)
setwd(file.path(Sys.getenv("THESIS_DIR"), "code/05_GO_annotation"))
load("0_identify_IDs.RData")

# Define the gene universe
genes_to_annotate <- union(
    rownames(dds_group),
    rownames(dds_LRT_interaction)
)

length(unique(genes_to_annotate))

# Extract protein IDs
gtf_tx <- gtf[gtf$V3 == "transcript", ]

gene_id <- str_match(gtf_tx$V9, "gene_id ([^;]+)")[, 2]
protein_id <- str_match(gtf_tx$V9, "evidenceProteinID ([^:]+)")[, 2]

map_loc_protein <- unique(data.frame(
    gene_id = gene_id,
    protein_id = protein_id,
    stringsAsFactors = FALSE
))

map_loc_protein <- map_loc_protein[
    map_loc_protein$gene_id %in% genes_to_annotate,
]

map_loc_protein <- map_loc_protein[!is.na(map_loc_protein$protein_id), ]

# Get GO annotations for the protein IDs
proteins <- unique(map_loc_protein$protein_id)
length(proteins)

writeLines(proteins, "protein_ids.txt")
