base::lapply(c("tidymodels",
               "tidyverse",
               "ggplot2",
               "ggrepel",
               "stringr"),
             require, character.only=T)

files_to_merge <- list.files("/home/rj/Desktop/NP137_rna_seq/count")

vec_ID_patient <- vector()
for(i in seq_along(files_to_merge)){
  extract_name <- stringr::str_extract(files_to_merge[i],"\\D*\\d*")
  vec_ID_patient[i] <- extract_name
}


# for(i in seq_along(files_to_merge)){
#   path_file <- paste0("/home/rj/Desktop/NP137_rna_seq/count/")
#   assign(vec_ID_patient[i],read.csv(paste0(path_file,files_to_merge[i]),sep=""))
# }


df_list <- list()
for(i in seq(1:length(files_to_merge))){
  path_file <- paste0("/home/rj/Desktop/NP137_rna_seq/count/")
  tmp <- read.csv(paste0(path_file,files_to_merge[i]),sep="",header=TRUE)
  tmp <- tmp[-1:-3,]
  tmp_name <- vec_ID_patient[i]
  df_list[[tmp_name]] <- tmp
}

for(i in vec_ID_patient){
  tmp_df <- df_list[[i]]
  ensg_point_look_up <- gsub("\\.[0-9]*$", "", tmp_df$N_unmapped)
  convert_ENSEMBL_id <- clusterProfiler::bitr(ensg_point_look_up,
                                              fromType = "ENSEMBL",
                                              toType = c("ENTREZID",
                                                         "SYMBOL",
                                                         "ENSEMBL"),
                                              OrgDb = "org.Hs.eg.db")
  rownames(tmp_df) <- ensg_point_look_up
  df_list[[i]] <- tmp_df
}

#ensembl list to search in the database and compute length of gene
ensembl_list <- rownames(df_list[[vec_ID_patient[1]]])
human <- useMart("ensembl", dataset="hsapiens_gene_ensembl")
gene_coords=getBM(attributes=c("hgnc_symbol","ensembl_gene_id", "start_position","end_position"), 
                  filters="ensembl_gene_id", 
                  values=ensembl_list, mart=human)
gene_coords$size=(gene_coords$end_position - gene_coords$start_position)/1000

df_list_to_tpm = list()
for(i in vec_ID_patient){
  tmp_df <- df_list[[i]]
  tmp_df <- tmp_df[convert_ENSEMBL_id$ENSEMBL,]
  tmp_df$ENSEMBLID <- convert_ENSEMBL_id$ENSEMBL
  tmp_df$SYMBOL <- convert_ENSEMBL_id$SYMBOL
  rownames(tmp_df) <- make.names(convert_ENSEMBL_id$ENSEMBL, unique=TRUE)
  df_symbol_length <- as.data.frame(cbind("SYMBOL"=gene_coords$hgnc_symbol,"length_g"=gene_coords$size))
  df_counts_with_length <- merge(df_symbol_length,tmp_df,by="SYMBOL")
  df_list_to_tpm[[i]] <- df_counts_with_length
}

df_list_tpm = list()
for(i in vec_ID_patient){
  tmp_df <- df_list_to_tpm[[i]]
  tmp_df$length_g <- as.numeric(tmp_df$length_g)
  tmp_df$rpk <- tmp_df[,6] / tmp_df$length_g
  scale_factor <- sum(tmp_df$rpk) / 1e6
  tmp_df$tpm <- tmp_df$rpk / scale_factor
  df_list_tpm[[i]] <- tmp_df
}

tpm_gene_name <- tmp_df[,1]
# concat row tmp
df_tpm <- data.frame(row.names = make.names(tpm_gene_name, unique=TRUE))
for(i in vec_ID_patient){
  tmp_df_tpm <- df_list_tpm[[i]]
  df_tpm <- cbind(df_tpm,tmp_df_tpm[,9])
}
colnames(df_tpm) <- vec_ID_patient

saveRDS(df_tpm,"TPM_NP137.RDS")

rownames(tmp_df) <- make.names(convert_ENSEMBL_id$ENSEMBL, unique=TRUE)
# concat row counts
df <- data.frame(row.names = make.names(tmp_df[,1], unique=TRUE))
for(i in df_list){
  tmp_df <- i
  df <- cbind(df,tmp_df[,4])
}
colnames(df) <- vec_ID_patient

clinical_data <- read.csv("/home/rj/Desktop/NP137_rna_seq/clinic20210921.csv",
                          sep=";")

clinical_data$Cancer.type <- tolower(clinical_data$Cancer.type)

logi <- str_detect(clinical_data$Cancer.type,'endometrium')
clinical_data <- clinical_data[logi,]
# 
C3D1  <- clinical_data[clinical_data$treatment == "C3D1",]
# 
C1D1  <- clinical_data[clinical_data$treatment == "C1D1",]

intersect(C3D1$patient,C1D1$patient)
# 
# 
# 
# intersect(rownames(df),EMT_gene_list_cor)
# 
# df <- as.data.frame(t(df))
# ensg_point_look_up <- gsub("\\.[0-9]*$", "", colnames(df))
# 
# convert_ENSEMBL_id<- clusterProfiler::bitr(ensg_point_look_up, 
#                                            fromType = "ENSEMBL", 
#                                            toType = c("ENTREZID",
#                                                       "SYMBOL",
#                                                       "ENSEMBL"), 
#                                            OrgDb = "org.Hs.eg.db")
# 
# colnames(df) <- ensg_point_look_up
# 
# df <- df[,convert_ENSEMBL_id$ENSEMBL]
# # renaming column with SYMBOL
# colnames(df) <- convert_ENSEMBL_id$SYMBOL
# 
# df <- as.data.frame(t(df))
# 
# intersect(C3D1$patient,C1D1$patient)
# 
C3D1_pat <- C3D1[C3D1$patient %in% intersect(C3D1$patient,C1D1$patient),]
C3D1_pat$Sample.Id
C1D1_pat <- C1D1[C1D1$patient %in% intersect(C3D1$patient,C1D1$patient),]
C1D1_pat$Sample.Id

df_C1D1 <- df_tpm[,C1D1_pat$Sample.Id]
df_C1D3 <- df_tpm[,C3D1_pat$Sample.Id]

# 
# df_C1D1_EMT_sig <- df_C1D1[intersect(rownames(df),EMT_gene_list_cor),]
# df_C1D3_EMT_sig <- df_C1D3[intersect(rownames(df),EMT_gene_list_cor),]
# 
# mat1 <- as.matrix(df_C1D1_EMT_sig)
# mat2 <- as.matrix(df_C1D3_EMT_sig)
# Heatmap(mat,rect_gp = gpar(col = "white", lwd = 2), column_order = c(C3D1_pat$Sample.Id,C1D1_pat$Sample.Id))
# 
# 
# 
# 
# Heatmap(mat2,rect_gp = gpar(col ="white",lwd = 2),column_order=C3D1_pat$Sample.Id)
# Heatmap(mat1,rect_gp = gpar(col ="white",lwd = 2),column_order=C1D1_pat$Sample.Id)
# 
# 
