requireNamespace('biomaRt', quietly = T)

# get gene ids and length info from the biomart
mart <- biomaRt::useMart("parasite_mart", dataset = "wbps_gene", host = "https://parasite.wormbase.org", port = 443)

Cel_genes <- biomaRt::getBM(attributes = c("wormbase_gene",          # WB gene ids
                                           "wormbase_transcript",    # Transcript names
                                           "external_gene_id",       # Public names
                                           "wormbase_gseq",          # Sequence names
                                           "transcript_start",       # for transcript lengths
                                           "transcript_end"          # for transcript lengths
                                           ), 
                            filters = c("species_id_1010"), values = c("caelegprjna13758"), # filter C. elegans only
                            mart = mart)

colnames(Cel_genes) <- c("wb_id", "transcript_name", "public_name", "sequence_name", "tr_start", "tr_end")

# get transcript length
Cel_genes$transcript_length <- Cel_genes$tr_end - Cel_genes$tr_start
Cel_genes <- Cel_genes[,c("wb_id", "transcript_name", "sequence_name", "public_name", "transcript_length")]

# save object to data
save('Cel_genes', file = "data/Cel_genes.RData")
rm(mart, Cel_genes)
