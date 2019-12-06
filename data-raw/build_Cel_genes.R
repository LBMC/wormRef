requireNamespace('biomaRt', quietly = T)

# get gene ids and length info from the biomart
mart <- biomaRt::useMart("ensembl", dataset = "celegans_gene_ensembl", version = "Ensembl Genes 98")

Cel_genes <- biomaRt::getBM(attributes = c("ensembl_gene_id",        # WB gene ids
                                           "ensembl_transcript_id",  # Sequence names
                                           "external_gene_name",     # Public names
                                           "transcript_length"),     # Transcript lengths
                            mart = mart)

colnames(Cel_genes) <- c("wb_id", "transcript_name", "public_name", "transcript_length")

# save object to data
save('Cel_genes', file = "data/Cel_genes.RData")
rm(mart, Cel_genes)
