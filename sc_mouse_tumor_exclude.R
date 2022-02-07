library(Seurat)
library(SeuratObject)
library(tidyverse)
library(ggpubr)
library(rstatix)
library(biomaRt)
library(future)
library(DESeq2)


set.seed(500)

load_sc_and_create_SeuratObject <- function(path_to_data, 
                                            column_to_chose = 2,
                                            project_name,
                                            min.cells = 100,
                                            min.features = 200){
  
  data <- Seurat::Read10X(data.dir = path_to_data,gene.column = 2)
  
  SO_data <- SeuratObject::CreateSeuratObject(counts = data, 
                                              min.cells = min.cells, 
                                              min.features = min.features, 
                                              project = project_name, 
                                              assay = "RNA")
  return(SO_data)
}

convertHumanGeneList <- function(x){
  
  new_config <- httr::config(ssl_verifypeer = FALSE)
  human <- useMart("ensembl", dataset = "hsapiens_gene_ensembl")
  mouse <- useMart("ensembl", dataset = "mmusculus_gene_ensembl")
  
  genesV2 <- getLDS(attributes = c("hgnc_symbol"), filters = "hgnc_symbol", 
                    values = x, mart = human, attributesL = c("mgi_symbol"), martL = mouse, uniqueRows=T)
  
  
  return(genesV2)
}

seurat_obj_list <- readRDS("/home/owly/Desktop/single_cell_analysis/single_cell_mouse/tumor_cluster_without_G2M_S_cell.RDS")

seurat_obj_list <- base::lapply(X = seurat_obj_list, FUN = function(x) {
  x[["percent_mt"]] <- Seurat::PercentageFeatureSet(x, 
                                                    pattern = "^mt-")
  x[["percent_mt"]] <- x[["percent_mt"]]/100
  
  x@meta.data <- x@meta.data %>% dplyr::rename(nUMI = nCount_RNA,
                                               nGene = nFeature_RNA)
  return(x)
})


feature_scatter_plot_list <- base::lapply(X = seurat_obj_list,
                                          FUN = function(x){
                                            plot1 <- Seurat::FeatureScatter(x, 
                                                                            feature1 = "nUMI", 
                                                                            feature2 = "percent_mt")
                                            plot2 <- Seurat::FeatureScatter(x, 
                                                                            feature1 = "nUMI", 
                                                                            feature2 = "nGene")
                                            plot3 <- plot1 + plot2
                                          })
pdf(file = paste0(getwd(),"/results/feature_scatter_plot_list_tumor_exclude.pdf"),
    width = 10,
    height = 10)
ggpubr::ggarrange(plotlist = feature_scatter_plot_list,
                  ncol = 1,
                  nrow = length(feature_scatter_plot_list)
)
dev.off()
rm(feature_scatter_plot_list)

# change nFeature_RNA value in function of the results in feature_scatter_plot_list
seurat_obj_list <- base::lapply(X = seurat_obj_list, 
                                FUN = function(x) {
                                  x <- Seurat::NormalizeData(x,
                                                             normalization.method = "LogNormalize",
                                                             scale.factor = 1e6)
                                  x <- subset(x, 
                                              subset = nGene > 300 & nGene < 7500 & percent_mt < 0.20)
                                })


variable_feature_plot_list <- base::lapply(X = seurat_obj_list, 
                                           FUN = function(x) {
                                             x <- Seurat::FindVariableFeatures(x, selection.method = "vst")
                                             
                                             # Identify the 10 most highly variable genes
                                             top10 <- head(VariableFeatures(x), 10)
                                             
                                             # plot variable features with and without labels
                                             plot1 <- VariableFeaturePlot(x)
                                             plot2 <- LabelPoints(plot = plot1, 
                                                                  points = top10, 
                                                                  repel = TRUE)
                                             plot3 <- plot1 + plot2
                                           })

seurat_obj_list <- base::lapply(X = seurat_obj_list, 
                                FUN = function(x) {
                                  x <- Seurat::FindVariableFeatures(x, 
                                                                    selection.method = "vst")
                                })

pdf(file = paste0(getwd(),"/results/variable_feature_plot_list_tumor_exclude.pdf"),
    width = 15,
    height = 20)
ggpubr::ggarrange(plotlist = variable_feature_plot_list,
                  ncol = 1,
                  nrow = length(variable_feature_plot_list),
                  labels=c("not treated", 
                           "treated"),
                  label.x = 0.4)
dev.off()
rm(variable_feature_plot_list)

# Integrated Analysis part 1 -------------------------------------------------------------------

features <- SelectIntegrationFeatures(object.list = seurat_obj_list)

seurat_obj_anchors <- FindIntegrationAnchors(object.list = seurat_obj_list, 
                                             anchor.features = features,
                                             reduction = "cca")

seurat_obj_anchors_combined <- IntegrateData(anchorset = seurat_obj_anchors)

rm(seurat_obj_anchors)

DefaultAssay(seurat_obj_anchors_combined) <- "integrated"

seurat_obj_anchors_combined <- ScaleData(seurat_obj_anchors_combined, 
                                         verbose = FALSE)
seurat_obj_anchors_combined <- RunPCA(seurat_obj_anchors_combined, 
                                      npcs = 30, 
                                      verbose = FALSE)
rm(seurat_obj_list)

ElbowPlot(seurat_obj_anchors_combined)

seurat_obj_anchors_combined <- RunUMAP(seurat_obj_anchors_combined, 
                                       reduction = "pca", 
                                       dims = 1:5,
                                       n.neighbors = 960L,
                                       min.dist = 0.3)
seurat_obj_anchors_combined <- FindNeighbors(seurat_obj_anchors_combined, 
                                             reduction = "pca", 
                                             dims = 1:5)
seurat_obj_anchors_combined <- FindClusters(seurat_obj_anchors_combined, 
                                            resolution = 0.1)
#saveRDS(seurat_obj_anchors_combined,"tumor_Rshiny.RDS")
p1 <- DimPlot(seurat_obj_anchors_combined, 
              reduction = "umap", 
              group.by = "orig.ident") + ggtitle("")
p2 <- DimPlot(seurat_obj_anchors_combined, 
              reduction = "umap", 
              label = TRUE, 
              repel = TRUE)
p3 <- p1 + p2
plot(p3)

metrics <- seurat_obj_anchors_combined@meta.data
endo_mak <- read.table("/home/owly/Desktop/single_cell_analysis/single_cell_mouse/EMT_signature_endo.txt")
endo_mak <- str_to_lower(endo_mak$V1)
genes_vector <- c("NTN1","UNC5B","NEO1","EPCAM","CLDN7","CD14","CDH1","CLDN4",
                  "HOOK1","MUC1","GRHL2","TFF3","CDS1","ESRP1","ESRP2","MARVELD2","F11R","CTNND1",
                  "VIM","VCAM1","FBN1","SPARC","FN1","COL6A3","CDH2","ZEB1","ZEB2","THBS2","SULF1",
                  "PCOLCE","INHBA","CALD1")
genes_vector <- c(genes_vector,"COL1A1","COL1A2","COL3A1","COL5A1","COL5A2","COL6A1","COL6A2","COL10A1",endo_mak)
genes_vector <- str_to_lower(genes_vector)
data <- seurat_obj_anchors_combined@assays[["RNA"]]@data
data <- data %>% as.data.frame()
rownames(data) <- str_to_lower(rownames(data))
genes_vector <- intersect(rownames(data),genes_vector)
data <- as.data.frame(data[genes_vector,])
data_t <- as.data.frame(t(data[,rownames(metrics[metrics$orig.ident == "treated",])]))
data_nt <- as.data.frame(t(data[,rownames(metrics[metrics$orig.ident == "not treated",])]))

bs   <- 100
res1 <- lapply(1:bs,
               function(i){
                 cat("boostrap",i,"/",as.character(bs),"\n")
                 tpm_data_t <- data_t[sample(
                   rownames(data_t),
                   size    = nrow(data_nt),
                   replace = TRUE
                 ),]
                 stat_test <- sapply(1:length(genes_vector),
                                     function(j){
                                       wilcox.test(
                                         tpm_data_t[,j],
                                         data_nt[,j]
                                       )[["p.value"]]
                                     })
                 
                 val_var <- sapply(1:length(genes_vector),
                                   function(j){
                                     var(
                                       tpm_data_t[,j]
                                     )
                                   })
                 
                 val_sd <- sapply(1:length(genes_vector),
                                  function(j){
                                    sd(
                                      tpm_data_t[,j]
                                    )
                                  })
                 
                 tmp_mean <- apply(
                   tpm_data_t,
                   2,
                   mean
                 )
                 return(list(
                   "stat_test"=stat_test,
                   "mean"=tmp_mean,
                   "variance"=val_var,
                   "std deviation"=val_sd
                 ))
                 
                 
               })


sub_res1 <- lapply(1:bs,function(i){
  res1[[i]][["mean"]]
})

sub_res2 <- lapply(1:bs,function(i){
  res1[[i]][["stat_test"]]
})

sub_res3 <- lapply(1:bs,function(i){
  res1[[i]][["variance"]]
})

sub_res4 <- lapply(1:bs,function(i){
  res1[[i]][["std deviation"]]
})


val_sd <- sapply(1:length(genes_vector),
                 function(i){
                   sd(
                     data_nt[,i]
                   )})
val_var <- sapply(1:length(genes_vector),
                  function(i){
                    var(
                      data_nt[,i]
                    )})



#boostrap wilcox result
df_stats_test <- do.call(rbind, sub_res2)
df_var <- do.call(rbind,sub_res3)
df_std <- do.call(rbind,sub_res4)

mean_bs_wilcox_test_p_value <- apply(df_stats_test,2,mean)
mean_var <- apply(df_var,2,mean)
mean_std <- apply(df_std,2,mean)


df_bs_mean <- do.call(rbind, sub_res1)
mean_t <- apply(df_bs_mean,2,mean)
mean_nt <- apply(data_nt,2,mean)
df_compare <- cbind(
  "SYMBOL"                         = genes_vector,
  "not treated"                    = mean_nt,
  "treated"                        = mean_t,
  "wilcox test (p value)"          = mean_bs_wilcox_test_p_value,
  "Variance not treated"           = mean_var,
  "Variance treated"               = val_var,
  "Standard deviation not treated" = mean_std,
  "Standard deviation treated"     = val_sd
)

df_compare <- df_compare %>% as.data.frame()
#df_compare <- df_compare[-c(1,28),]


df_compare_mean <- apply(df_compare[,2:3],1,as.numeric)
mean_genes <- apply(df_compare_mean,2,mean)


test <- lapply(1:bs,function(i){
  tpm_data_t <- data_t[sample(
    rownames(data_t),
    size    = nrow(data_nt),
    replace = TRUE
  ),]
  
  val_percent <- sapply(1:length(mean_genes),
                        function(j){
                          length(which(tpm_data_t[,j] > mean_genes[j]))*100/nrow(tpm_data_t)
                        })
  return(val_percent)
})

df_percent_t <- do.call(rbind, test)
mean_percent <- apply(df_percent_t,2,mean)

percent_nt <- lapply(1:length(mean_genes),function(i){
  percent_t <- length(which(data_nt[,i] > mean_genes[i]))*100/nrow(data_nt)
})

df_percent_nt <- as.data.frame(do.call(rbind, percent_nt))

df_compare$`treated cell %` <- mean_percent
df_compare$`not treated cell %` <- df_percent_nt


write_excel_csv(as.data.frame(as.matrix(df_compare)),"mouse_boostrap_signature_enrichie.csv")

df_compare[,4] <- as.numeric(df_compare[,4])
keep_gene <- df_compare[df_compare$`wilcox test (p value)` < 0.05,]
keep_gene <- na.omit(keep_gene)
write_excel_csv(as.data.frame(as.matrix(keep_gene)),"mouse_boostrap_signature_enrichie_significative.csv")

pdf(file = paste0(getwd(),
                  "/results/umap_plot_analysis_tumor_exclude.pdf"),
    width = 15,
    height = 6)
plot(p3)
dev.off()

p4 <- DimPlot(seurat_obj_anchors_combined, 
              reduction = "umap", split.by = "orig.ident")

pdf(file = paste0(getwd(),"/results/umap_plot_analysis_split_by_orig_tumor_exclude.pdf"),
    width = 7,
    height = 6)
plot(p4)
dev.off()

library(ggrepel)

piechart <- function(data = data, title = "", logical_palette_color = F,palette_color = "Set1"){
  data <- as.data.frame(data)
  df <- data %>%
    arrange(desc(Var1)) %>%
    mutate(prop = Freq / sum(Freq) *100) %>%
    mutate(ypos = cumsum(prop) - 0.4*prop)
  
  df$prop <- round(df$prop,0)
  
  p <- ggplot(df, aes(x="", y=prop, fill=Var1,label=Var1)) +
    geom_bar(stat="identity", width=1, color="white") +
    coord_polar("y", start=0) +
    theme_void() +
    guides(fill = guide_legend(title = "Cell type"))
  if(logical_palette_color == TRUE){
    p + scale_fill_manual(values=palette_color)
  }else{
    p + scale_fill_brewer(palette = palette_color)
  }
  p <- p +  geom_text_repel(aes(y = ypos,
                                label = paste0(round(prop,2),"%")), 
                            color = "black", size=6) + 
    labs(title = as.character(title))
  
  return(
    list(
      "table" = df,
      "plot"  = p)
  )
}

palette_to_color <- c("#F87766D","#7CAE00","#00BFC4","#C77CF77")

metrics <- seurat_obj_anchors_combined@meta.data
nt <- metrics[metrics$orig.ident == "not treated",]
t  <- metrics[metrics$orig.ident == "treated",]
t <- metrics[sample(rownames(t),size= nrow(nt)),]
table_nt <- table(nt$seurat_clusters)
table_t <- table(t$seurat_clusters)

table_t_plot <- piechart(table_t,
                         "Treated",
                         logical_palette_color = TRUE,
                         palette_color = palette_to_color)
table_nt_plot <- piechart(table_nt,
                          "Not treated",
                          logical_palette_color = TRUE,
                          palette_color = palette_to_color)



cluster.markers <- FindMarkers(seurat_obj_anchors_combined, ident.1 = 1,ident.2 = 0, min.pct = 0.25)
cluster.markers <- cluster.markers[abs(cluster.markers$avg_log2FC) > 2 & cluster.markers$p_val_adj < 0.05,]

varnames <- rownames(seurat_obj_anchors_combined@assays[["RNA"]]@data)
data <- seurat_obj_anchors_combined@assays[["RNA"]]@data

for(i in rownames(cluster.markers)){
  png(file = paste0(paste0(paste0(getwd(),"/results/cluster_deg_0_1/features_plot_"),i),".png"),
      width = 1040,
      height = 555)
  p <- suppressWarnings(FeaturePlot(seurat_obj_anchors_combined, 
                                    features = c(i),
                                    pt.size = 1.2,
                                    min.cutoff = 0,
                                    max.cutoff = max(data[i,]) - quantile(data[i,],0.25),
                                    keep.scale = "NULL",
                                    split.by = "orig.ident",
                                    cols = c("grey", "red")) + theme(legend.position ="none"))
  plot(p)
  dev.off()
}

cluster.markers <- FindMarkers(seurat_obj_anchors_combined, ident.1 = 2,ident.2 = 1, min.pct = 0.25)
cluster.markers <- cluster.markers[abs(cluster.markers$avg_log2FC) > 2 & cluster.markers$p_val_adj < 0.05,]

data <- seurat_obj_anchors_combined@assays[["RNA"]]@data

for(i in rownames(cluster.markers)){
  png(file = paste0(paste0(paste0(getwd(),"/results/cluster_deg_1_2/features_plot_"),i),".png"),
      width = 1040,
      height = 555)
  p <- suppressWarnings(FeaturePlot(seurat_obj_anchors_combined, 
                                    features = c(i),
                                    pt.size = 1.2,
                                    min.cutoff = 0,
                                    max.cutoff = max(data[i,]) - quantile(data[i,],0.25),
                                    keep.scale = NULL,
                                    split.by = "orig.ident",
                                    cols = c("grey", "red")) + theme(legend.position ="none"))
  plot(p)
  dev.off()
}

# saveRDS(seurat_obj_anchors_combined_markers,"DEG_treated_vs_not_treated.RDS")
# 
# seurat_obj_anchors_combined_markers <- readRDS("/home/owly/Desktop/single_cell_analysis/single_cell_mouse/DEG_treated_vs_not_treated.RDS")
# 
seurat_obj_anchors_combined_markers %>%
  group_by(cluster) %>%
  top_n(n = 10, wt = avg_log2FC) -> top10
DoHeatmap(seurat_obj_anchors_combined, features = top10$gene)
# 
# pdf(file = paste0(getwd(),"/results/heatmap_markers_top10_analysis_tumor_exclude.pdf"),
#     width = 15,
#     height = 10)
# DoHeatmap(seurat_obj_anchors_combined, features = top10$gene) + NoLegend()
# dev.off()

# saveRDS(as.data.frame(seurat_obj_anchors_combined@assays[["RNA"]]@counts),"tumor_data.RDS")
# saveRDS(seurat_obj_anchors_combined@meta.data,"metrics_tumor_data.RDS")

# n_gene <- "Syne2"
# FeaturePlot(seurat_obj_anchors_combined, features = c(n_gene)) + ggtitle(n_gene)
# VlnPlot(seurat_obj_anchors_combined, 
#         features = c(n_gene), 
#         slot = "counts", 
#         log = TRUE,
#         pt.size = 0) + ggtitle(n_gene)
# 
# cluster_2_gene <- seurat_obj_anchors_combined_markers[seurat_obj_anchors_combined_markers$cluster == 2,]$gene
# 
# markers <- c("Epcam","Acta2","Pecam1","Ptprc","C1qa","Cd3d")
# nrow_feature_plot <- round(length(markers)/4)
# pdf(file = paste0(getwd(),"/results/umap_markers.pdf"),
#     width = 15,
#     height = 40)
# FeaturePlot(seurat_obj_anchors_combined, 
#             features = markers,
#             split.by ="orig.ident",
#             min.cutoff = 0.01)
# dev.off()
# 
# markers <- c("Vcam1","Itgb3","Itgav","Vim","Pten")
# nrow_feature_plot <- round(length(markers)/4)
# pdf(file = paste0(getwd(),"/results/umap_markers.pdf"),
#     width = 15,
#     height = 40)
# FeaturePlot(seurat_obj_anchors_combined, 
#             features = markers,
#             split.by ="orig.ident",
#             min.cutoff = 0.01)
# dev.off()
# 
# FeaturePlot(seurat_obj_anchors_combined, 
#             features = c("Ntn1","Itgb1","App","Adora2b","Unc50"),
#             split.by ="orig.ident",
#             min.cutoff = 0.01)
# 
# nrow_feature_plot <- round(length(markers)/4)
# pdf(file = paste0(getwd(),"/results/umap_markers_ntn1.pdf"),
#     width = 15,
#     height = 40)
# FeaturePlot(seurat_obj_anchors_combined, 
#             features = c("Ntn1","Itgb1","App","Adora2b","Unc50"),
#             split.by ="orig.ident",
#             min.cutoff = 0.01)
# dev.off()
# 
# pdf(file = paste0(getwd(),"/results/umap_markers_ntn1_V2.pdf"),
#     width = 15,
#     height = 40)
# FeaturePlot(seurat_obj_anchors_combined, 
#             features = c("Ntn1","Col1a1","Col6a3","Acta2","Vcam1"),
#             split.by ="orig.ident",
#             min.cutoff = 0.01)
# dev.off()
# 
# pdf(file = paste0(getwd(),"/results/umap_markers_col.pdf"),
#     width = 15,
#     height = 80)
# FeaturePlot(seurat_obj_anchors_combined, 
#             features = c("Col5a2","Col4a4","Col15a1","Col9a2","Col12a1",
#                          "Col18a1","Col14a1","Col11a2","Col4a6","Col4a5",
#                          "Col3a1","Col27a1","Col1a2","Col6a2","Col6a1","Col1a1"),
#             split.by ="orig.ident",
#             min.cutoff = 0.01)
# dev.off()
# 
# pdf(file = paste0(getwd(),"/results/umap_markers_cldn.pdf"),
#     width = 15,
#     height = 80)
# FeaturePlot(seurat_obj_anchors_combined, 
#             features = varname[grep(pattern = "Cldn",varname)][1:9],
#             split.by ="orig.ident",
#             min.cutoff = 0.01)
# dev.off()
# 
# pdf(file = paste0(getwd(),"/results/umap_markers_cam.pdf"),
#     width = 15,
#     height = 80)
# FeaturePlot(seurat_obj_anchors_combined, 
#             features = c("Vcam1","Bcam","Epcam","Ncam1","Icam1","Caml",
#                          "Icam2","Ticam1","Alcam","Pecam1","Ceacam1"),
#             split.by ="orig.ident",
#             min.cutoff = 0.01)
# dev.off()
# 
# varname[grep(pattern = "Unc",varname)]
# 
# pdf(file = paste0(getwd(),"/results/umap_markers_Itg.pdf"),
#     width = 15,
#     height = 80)
# FeaturePlot(seurat_obj_anchors_combined, 
#             features = varname[grep(pattern = "Itg",varname)],
#             split.by ="orig.ident",
#             min.cutoff = 0.01)
# dev.off()
# 
# pdf(file = paste0(getwd(),"/results/umap_markers_unc.pdf"),
#     width = 15,
#     height = 80)
# FeaturePlot(seurat_obj_anchors_combined, 
#             features = varname[grep(pattern = "Unc",varname)]
# ,
#             split.by ="orig.ident",
#             min.cutoff = 0.01)
# dev.off()
# 
# pdf(file = paste0(getwd(),"/results/umap_markers_krt.pdf"),
#     width = 15,
#     height = 80)
# FeaturePlot(seurat_obj_anchors_combined, 
#             features = varname[grep(pattern = "Krt",varname)]
#             ,
#             split.by ="orig.ident",
#             min.cutoff = 0.01)
# dev.off()
# 
# 
# pdf(file = paste0(getwd(),"/results/umap_markers_cdh.pdf"),
#     width = 15,
#     height = 80)
# FeaturePlot(seurat_obj_anchors_combined, 
#             features = varname[grep(pattern = "Cdh",varname)]
#             ,
#             split.by ="orig.ident",
#             min.cutoff = 0.01)
# dev.off()
# 
# pdf(file = paste0(getwd(),"/results/umap_markers_cdc.pdf"),
#     width = 15,
#     height = 120)
# FeaturePlot(seurat_obj_anchors_combined, 
#             features = varname[grep(pattern = "Cdc",varname)]
#             ,
#             split.by ="orig.ident",
#             min.cutoff = 0.01)
# dev.off()
# 
# pdf(file = paste0(getwd(),"/results/umap_markers_cdk.pdf"),
#     width = 15,
#     height = 120)
# FeaturePlot(seurat_obj_anchors_combined, 
#             features = varname[grep(pattern = "Cdk",varname)]
#             ,
#             split.by ="orig.ident",
#             min.cutoff = 0.01)

# DEG analysis ------------------------------------------------------------

metrics <- seurat_obj_anchors_combined@meta.data

endo_genes_emt <- read.table("/home/owly/Desktop/Project_TCGA_signature/scripts/EMT_signature_endo.txt", header = T)

endo_gene_epi <- endo_genes_emt[endo_genes_emt$EMT_status == "Epithelial",]$SYMBOL
endo_gene_mes <- endo_genes_emt[endo_genes_emt$EMT_status == "Mesenchymal",]$SYMBOL

convert_epi_genes <- convertHumanGeneList(endo_gene_epi)
convert_mes_genes <- convertHumanGeneList(endo_gene_mes)

convert_mus_epi_genes <- convert_epi_genes$MGI.symbol
convert_mus_mes_genes <- convert_mes_genes$MGI.symbol
convert_mus_genes <- c(convert_mus_epi_genes,convert_mus_mes_genes)



# # START boostrap ----------------------------------------------------------------

function_resample <- function(df,logical_margin,metadata,metadata_col,vec_categorial_var){
  if(logical_margin == 1){
    rowname_data_1 <- rownames(metadata[metadata_col == vec_categorial_var[1],])
    data_1 <- df[,rowname_data_1]
    rowname_data_2 <- rownames(metadata[metadata_col == vec_categorial_var[2],])
    data_2 <- sample(df[,rowname_data_2],length(data_1))

    id_to_select <- c(colnames(data_1),colnames(data_2))
    metadata <- metadata[id_to_select,]
    df <- as.data.frame(t(cbind(data_1,data_2)))

    return(
      list(
           "df"= df,
           "metadata"=metadata)
      )}else{
       stop()
  }
}


function_compute_score_mak <- function(df,vec1_var,vec2_var,metadata){
  mean_vec1_var <- apply(df[,intersect(vec1_var,colnames(df))],1,mean)
  mean_vec2_var <- apply(df[,intersect(vec2_var,colnames(df))],1,mean)
  metadata$score_mak <- mean_vec2_var-mean_vec1_var
  metadata$log10_score_mak <- log10(metadata$score_mak - (min(metadata$score_mak)-1))
  return(metadata)
}

boostrap_wilcox_test <- function(df){
  data <- function_resample(df,
                            1,
                            metrics,
                            metrics$orig.ident,
                            c("not treated","treated"))

  metrics_data <- function_compute_score_mak(data[["df"]],
                                             convert_mus_epi_genes,
                                             convert_mus_mes_genes,
                                             data[["metadata"]])

  stat_test <- metrics_data %>% wilcox_test(as.formula(paste("score_mak", "~ orig.ident"))) %>% adjust_pvalue(method = "bonferroni") %>% add_significance("p.adj")
  stat_test <- as.data.frame(stat_test)
  stat_test$mean_score_mak_not_treated <- mean(metrics_data[metrics_data$orig.ident == "not treated",]$score_mak)
  stat_test$sd_score_mak_not_treated <- sd(metrics_data[metrics_data$orig.ident == "not treated",]$score_mak)
  stat_test$var_score_mak_not_treated <- var(metrics_data[metrics_data$orig.ident == "not treated",]$score_mak)
  stat_test$mean_score_mak_treated <- mean(metrics_data[metrics_data$orig.ident == "treated",]$score_mak)
  stat_test$sd_score_mak_treated <- sd(metrics_data[metrics_data$orig.ident == "treated,"]$score_mak)
  stat_test$sd_score_mak_treated <- var(metrics_data[metrics_data$orig.ident == "treated",]$score_mak)

  return(stat_test)
}


NB=100
boot_list=list()
for(i in seq(1:NB)){
  cat("boostrap ",i,"\n")
  boot_list[[as.character(i)]] <- boostrap_wilcox_test(data_RNA)
}
# # END bootstrap -----------------------------------------------------------


# data_RNA_nt <- data_RNA[,rownames(metrics[metrics$orig.ident == "not treated",])]
# data_RNA_t <- sample(data_RNA[,rownames(metrics[metrics$orig.ident == "treated",])],length(data_RNA_nt))
# id_cells <- c(colnames(data_RNA_nt),colnames(data_RNA_t))
# metrics <- metrics[id_cells,]
# 
# data_RNA <- as.data.frame(t(cbind(data_RNA_nt,data_RNA_t)))
# mean_epi <- apply(data_RNA[,intersect(convert_mus_epi_genes,colnames(data_RNA))],1,mean)
# mean_mes <- apply(data_RNA[,intersect(convert_mus_mes_genes,colnames(data_RNA))],1,mean)
# metrics$score_mak <- mean_mes-mean_epi
# 
# median(metrics[metrics$orig.ident == "not treated",]$score_mak)
# sd(metrics[metrics$orig.ident == "not treated",]$score_mak)
# median(metrics[metrics$orig.ident == "treated",]$score_mak)
# sd(metrics[metrics$orig.ident == "treated",]$score_mak)
# 
# 
# wilcox_test(metrics,as.formula(paste("score_mak", "~ orig.ident")))
# min(metrics$score_mak)
# metrics$log10_score_mak <- log10(metrics$score_mak - (min(metrics$score_mak)-1))
# 
# rm(stat.test)
# 
#   stat.test <- metrics %>%
#     wilcox_test(as.formula(paste("log10_score_mak", "~ orig.ident"))) %>%
#     # wilcox_test(as.formula(paste(i, "~ treatment")),paired=TRUE)
#     adjust_pvalue(method = "bonferroni") %>%
#     add_significance("p.adj")
#   
#   bxp <- ggplot(data = metrics,aes_string(x="orig.ident", y = "log10_score_mak")) +
#     geom_boxplot(aes(fill = orig.ident)) +
#     scale_fill_manual(values = c("coral1","palegreen2"))
#   # Ajoutez des p-values sur les graphiques en box plot
#   stat.test <- stat.test %>%
#     add_xy_position(x = "orig.ident", dodge = 0.8)
#   stat.test$y.position = 1.5
#   
#   
#   # bxp <- bxp + stat_pvalue_manual(
#   #  stat.test,  label = "{p.adj}{p.adj.signif}", tip.length = 0
#   # )
#   
#   # Ajoutez 10 % d'espaces entre les étiquettes des p-values et la bordure du graphique
#   bxp <- bxp + stat_pvalue_manual(
#     stat.test,  label = "{p.adj.signif}", tip.length = .02, hide.ns = FALSE,
#     step.increase = 0.05) +
#     scale_y_continuous(expand = expansion(mult = c(0, 0.1))) +
#     theme(plot.title = element_text(hjust = 0.5,size=20),
#           panel.background = element_rect(fill="white", colour="black", size=0.5,
#                                           linetype="solid"),
#           panel.grid.major = element_line(size = 0.5, linetype = 0,
#                                           colour = "white"),
#           panel.grid.minor = element_line(size = 0.25, linetype = 'solid',
#                                           colour = "azure3"),
#           legend.text = element_text(size=15),
#           legend.title = element_text(size=20),
#           axis.text = element_text(size=12),
#           axis.title = element_text(size=15)) +
#     ggtitle("tumor cells treated vs not treated") +
#     ylab("score MAK") +
#     xlab("treatment") +
#     labs(fill="treatment")
  
  
suppressPackageStartupMessages(library(escape))
suppressPackageStartupMessages(library(dittoSeq))
suppressPackageStartupMessages(library(singscore))

# 
gene.sets <- list(
  "epithelial" = convert_epi_genes$MGI.symbol,
  "mesenchymal" = convert_mes_genes$MGI.symbol
)

ES <- enrichIt(obj = seurat_obj_anchors_combined,
               gene.sets = gene.sets,
               groups = 1000, cores = 4)

seurat_obj_anchors_combined <- AddMetaData(seurat_obj_anchors_combined, ES)

seurat_obj_anchors_combined@meta.data$active.idents <- seurat_obj_anchors_combined@active.ident
colorblind_vector <- colorRampPalette(rev(c("#0D0887FF", "#47039FFF",
                                            "#7301A8FF", "#9C179EFF", "#BD3786FF", "#D8576BFF",
                                            "#ED7953FF","#FA9E3BFF", "#FDC926FF", "#F0F921FF")))
DefaultAssay(seurat_obj_anchors_combined) <- "integrated"

dittoPlot(seurat_obj_anchors_combined, "epithelial", group.by = "orig.ident") +
  scale_fill_manual(values = c("coral1","cornflowerblue"))

dittoPlot(seurat_obj_anchors_combined, "mesenchymal", group.by = "orig.ident") +
  scale_fill_manual(values = c("coral1","cornflowerblue"))


dittoScatterHex(seurat_obj_anchors_combined,
                x.var = "mesenchymal",
                y.var = "epithelial",
                do.contour = TRUE,
                split.by = "orig.ident") +
  theme_classic() +
  scale_fill_gradientn(colors = rev(colorblind_vector(11))) +
  geom_vline(xintercept = 0, lty=2) +
  geom_hline(yintercept = 0, lty=2)
