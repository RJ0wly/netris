library(tidyverse)
counts <- read.csv("/home/rj/Desktop/NP137_rna_seq/count/F3899ReadsPerGene.out.tab",sep="",header=TRUE)
counts <- counts[-1:-3,]

#convert ensembl ID
ensg_point_look_up <- gsub("\\.[0-9]*$", "", counts$N_unmapped)
convert_ENSEMBL_id <- clusterProfiler::bitr(ensg_point_look_up,
                                           fromType = "ENSEMBL",
                                           toType = c("ENTREZID",
                                                      "SYMBOL",
                                                      "ENSEMBL"),
                                           OrgDb = "org.Hs.eg.db")
#changing rownames
rownames(counts) <- ensg_point_look_up

#ensembl list to search in the database

ensembl_list <- rownames(counts)
human <- useMart("ensembl", dataset="hsapiens_gene_ensembl")
gene_coords=getBM(attributes=c("hgnc_symbol","ensembl_gene_id", "start_position","end_position"), 
                  filters="ensembl_gene_id", 
                  values=ensembl_list, mart=human)
gene_coords$size=(gene_coords$end_position - gene_coords$start_position)/1000

#subsetting counts to only have ensemblID that have been find in the database
counts <- counts[convert_ENSEMBL_id$ENSEMBL,]
#adding a column with the ensembl ID
counts$ENSEMBLID <- convert_ENSEMBL_id$ENSEMBL
#adding a column with the gene symbol ID
counts$SYMBOL <- convert_ENSEMBL_id$SYMBOL
#changing rownames to the converted ensembl ID
rownames(counts) <- make.names(convert_ENSEMBL_id$ENSEMBL, unique=TRUE)
#combining symbol ID with size of the gene
df_ensembl_length <- as.data.frame(cbind("SYMBOL"=gene_coords$hgnc_symbol,"length_g"=gene_coords$size))
#merging df_ensembl_length with counts
df_counts_length <- merge(df_ensembl_length,counts,by="SYMBOL")
#converting counts to tpm
df_counts_length$length_g <- as.numeric(df_counts_length$length_g)
df_counts_length$rpk <- df_counts_length$X1774736.2 / df_counts_length$length_g
scale_factor <- sum(df_counts_length$rpk) / 1e6
df_counts_length$tpm <- df_counts_length$rpk / scale_factor
