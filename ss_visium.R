base::lapply(c("Seurat",
               "tidyverse",
               "SeuratObject",
               "SeuratData",
               "patchwork",
               "tidyverse",
               "singscore",
               "escape",
               "dittoSeq",
               "ggpubr",
               "ggrepel",
               "cowplot",
               "RColorBrewer",
               "ggrepel",
               "harmony",
               "rstatix"),
             require, character.only=T)

load_visium <- Seurat::Load10X_Spatial(data.dir   = "/home/owly/Desktop/analysis/visium/visium/folder_visium/synovial_sarcoma/NYESO_144",
                                       filename   = "filtered_feature_bc_matrix.h5",
                                       assay      = "Spatial",
                                       slice      = "NYESO144")

# violin plot of ncounts per spot
plot1 <- VlnPlot(load_visium, 
                 features = colnames(load_visium@meta.data)[2], 
                 pt.size  = 0.1) + NoLegend()
plot2 <- SpatialFeaturePlot(load_visium, 
                            features = colnames(load_visium@meta.data)[2]) + 
  theme(legend.position = "right")
display_p1_p2 <- wrap_plots(plot1, 
                            plot2)

# normalization with STC algorithm which builds regularized negative 
# binomial models of gene expression in order to account for technical artifacts 
# while preserving biological variance
normalized_data <- SCTransform(load_visium, 
                               assay   = "Spatial",
                               verbose = FALSE,
                               variable.features.n = 3000)

# feature plot 
spatial_features_plot <- SpatialFeaturePlot(normalized_data, 
                                            features = c("NTN1","UNC5B" ,"EPCAM","PTPRC","PECAM1","ACTA2"),
                                            alpha = c(0.05, 1))

SpatialFeaturePlot(normalized_data, 
                   features = c(
                     "KRT7"),
                   alpha = c(0.05, 1))

normalized_data <- RunPCA(normalized_data, 
                          assay   = "SCT", 
                          verbose = FALSE)
normalized_data <- FindNeighbors(normalized_data, 
                                 reduction = "pca", 
                                 dims      = 1:30)
# algorithm 4 is Leiden clustering algorithm (supposed 
# to be better than louvain algorithm)
normalized_data <- FindClusters(normalized_data, 
                                verbose    = FALSE,
                                algorithm  = 4,
                                resolution = 0.2)

normalized_data <- RunUMAP(normalized_data, 
                           reduction = "pca", 
                           dims      = 1:30)


normalized_data <- RunTSNE(normalized_data,
                           reduction = "pca",
                           dims =  1:30)

p1 <- DimPlot(normalized_data, 
              reduction = "umap", 
              label     = TRUE,
              pt.size = 2)
p2 <- SpatialDimPlot(normalized_data, 
                     label      = TRUE,
                     label.size = 3)
p3 <- p1 + p2

normalized_data_data_marker <- FindAllMarkers(normalized_data,
                                              test.use = "bimod",
                                              min.pct = 0.15,
                                              logfc.threshold = 0.15)

normalized_data_data_marker %>%
  group_by(cluster) %>%
  top_n(n = 10, wt = avg_log2FC) -> top10

p_hm <- DoHeatmap(normalized_data, 
                  features = top10$gene) + NoLegend()  +
  scale_fill_gradientn(colors = c("seagreen3", "black", "brown1"))


data <- as.data.frame(t(as.matrix(normalized_data@assays[["SCT"]]@data)))

normalized_data@meta.data$NTN1 <- data[,"NTN1"]

normalized_data@meta.data <- normalized_data@meta.data %>% mutate(NTN1_status = case_when(NTN1 > quantile(normalized_data@meta.data$NTN1,0.5) ~ "High",
                                                                                          TRUE ~ "Low"))
median(normalized_data@meta.data[normalized_data@meta.data$NTN1_status == "High",]$NTN1)

Idents(normalized_data) <- "NTN1_status"
normalized_data_data_marker <- FindAllMarkers(normalized_data,
                                              test.use = "bimod",
                                              min.pct = 0.15,
                                              logfc.threshold = 0.15)

normalized_data_data_marker %>%
  group_by(cluster) %>%
  top_n(n = 10, wt = avg_log2FC) -> top10

p_hm <- DoHeatmap(normalized_data, 
                  features = top10$gene) + NoLegend()  +
  scale_fill_gradientn(colors = c("seagreen3", "black", "brown1"))

gene.setsfull <- getGeneSets(library = "H")

ES_full <- enrichIt(obj = as.matrix(normalized_data@assays[["SCT"]]@counts),
                    gene.sets = gene.setsfull,
                    groups = 200, cores = 4)

normalized_data <- AddMetaData(normalized_data, ES_full)

gene_mak <- read.table(paste0(getwd(),"/signature_mak_1.csv"),sep=",",header=F)

gene_set_mak <- list("MAK_epithelial"= gene_mak[gene_mak$V2 == "E",]$V1,
                     "MAK_mesenchymal"= gene_mak[gene_mak$V2 == "M",]$V1)

ESmakfull <- enrichIt(obj = as.matrix(normalized_data@assays[["SCT"]]@counts),
                      gene.sets = gene_set_mak,
                      groups = 200, cores = 4 )

normalized_data <- AddMetaData(normalized_data, ESmakfull)


wilcox.test(normalized_data$HALLMARK_EPITHELIAL_MESENCHYMAL_TRANSITION[normalized_data$NTN1_status =="High"],
normalized_data$HALLMARK_EPITHELIAL_MESENCHYMAL_TRANSITION[normalized_data$NTN1_status =="Low"])

ggplot(normalized_data@meta.data,aes(x = NTN1_status, y = HALLMARK_EPITHELIAL_MESENCHYMAL_TRANSITION, fill=NTN1_status)) +
  geom_boxplot() +
  theme_light()

normalized_data@meta.data$score_MAK <-normalized_data@meta.data$MAK_mesenchymal - normalized_data@meta.data$MAK_epithelial

ggplot(normalized_data@meta.data,aes(x = NTN1_status, y = score_MAK, fill=NTN1_status)) +
  geom_boxplot() +
  theme_light()

i <- "HALLMARK_EPITHELIAL_MESENCHYMAL_TRANSITION"
stat_test <- normalized_data@meta.data %>%
  wilcox_test(as.formula(paste(as.character(i), "~", "NTN1_status"))) %>%
  add_significance()
stat_test <- stat_test %>% add_xy_position(x = "NTN1_status")

bxp <- ggboxplot(data = normalized_data@meta.data, 
                 x    = "NTN1_status", 
                 y    = i, 
                 fill = "NTN1_status")
p <- bxp + 
  stat_pvalue_manual(stat_test, label = "p.signif") +
  scale_y_continuous(expand = expansion(mult = c(0.05, 0.1))) +
  ylab(paste(i," expression")) + 
  xlab("") +
  scale_fill_discrete("NTN1 status") +
  ggtitle(paste(i))


saveRDS(normalized_data,"obj_NY144.RDS")

normalized_data <- readRDS("obj_NY144.RDS")

Idents(normalized_data) <- "seurat_clusters"
normalized_data_1_2 <- subset(normalized_data, idents = c(1, 2))
p1 <- SpatialDimPlot(normalized_data_1_2, crop = TRUE, label = TRUE)
p2 <- SpatialDimPlot(normalized_data_1_2, crop = FALSE, label = TRUE, 
                     pt.size.factor = 1, label.size = 3)
p1 + p2

p1 <- SpatialDimPlot(normalized_data_1_2, crop = TRUE,group.by = "NTN1_status",pt.size.factor = 1.2)
p2 <- SpatialFeaturePlot(normalized_data_1_2, 
                         features = c(
                           "KRT7"),
                         alpha = c(0.05, 1))
p3 <- SpatialFeaturePlot(normalized_data_1_2, 
                         features = c(
                           "CD24"),
                         alpha = c(0.05, 1))
p4 <- SpatialFeaturePlot(normalized_data_1_2, 
                         features = c(
                           "NPPC"),
                         alpha = c(0.05, 1))

p5 <- p1 + p2 + p3 + p4

SpatialFeaturePlot(normalized_data_1_2, 
                   features = c(
                     "NTN1"),
                   alpha = c(0.05, 1))

data <- as.data.frame(t(as.matrix(normalized_data_1_2@assays[["SCT"]]@data)))

normalized_data_1_2@meta.data$NTN1 <- data[,"NTN1"]

normalized_data_1_2@meta.data <- normalized_data_1_2@meta.data %>% mutate(NTN1_status = case_when(NTN1 > quantile(normalized_data_1_2@meta.data$NTN1,0.5) ~ "High",
                                                                                          TRUE ~ "Low"))
table(normalized_data_1_2@meta.data$NTN1_status)

Idents(normalized_data_1_2) <- "NTN1_status"
normalized_data_data_marker <- FindAllMarkers(normalized_data_1_2,
                                              test.use = "wilcox",
                                              min.pct = 0.15,
                                              logfc.threshold = 0.15)


gene.setsfull <- getGeneSets(library = "H")

ES_full <- enrichIt(obj = as.matrix(normalized_data_1_2@assays[["SCT"]]@counts),
                    gene.sets = gene.setsfull,
                    groups = 200, cores = 4)

normalized_data_1_2 <- AddMetaData(normalized_data_1_2, ES_full)

gene_mak <- read.table(paste0(getwd(),"/signature_mak_1.csv"),sep=",",header=F)

gene_set_mak <- list("MAK_epithelial"= gene_mak[gene_mak$V2 == "E",]$V1,
                     "MAK_mesenchymal"= gene_mak[gene_mak$V2 == "M",]$V1)

ESmakfull <- enrichIt(obj = as.matrix(normalized_data_1_2@assays[["SCT"]]@counts),
                      gene.sets = gene_set_mak,
                      groups = 200, cores = 4 )

normalized_data_1_2 <- AddMetaData(normalized_data_1_2, ESmakfull)

normalized_data_1_2@meta.data$score_MAK <- normalized_data_1_2@meta.data$MAK_mesenchymal - normalized_data_1_2@meta.data$MAK_epithelial

normalized_data_data_marker %>%
  group_by(cluster) %>%
  top_n(n = 5, wt = avg_log2FC) -> top10

p_hm <- DoHeatmap(normalized_data_1_2, 
                  features = top10$gene,
                  slot = "scale.data") + NoLegend()  +
  scale_fill_gradientn(colors = c("seagreen3", "black", "brown1"))


i <- "score_MAK"
stat_test <- normalized_data_1_2@meta.data %>%
  wilcox_test(as.formula(paste(as.character(i), "~", "NTN1_status"))) %>%
  add_significance()
stat_test <- stat_test %>% add_xy_position(x = "NTN1_status")

bxp <- ggboxplot(data = normalized_data_1_2@meta.data, 
                 x    = "NTN1_status", 
                 y    = i, 
                 fill = "NTN1_status")
p <- bxp + 
  stat_pvalue_manual(stat_test, label = "p.signif") +
  scale_y_continuous(expand = expansion(mult = c(0.1, 0.1))) +
  ylab(paste(i," expression")) + 
  xlab("") +
  scale_fill_discrete("NTN1 status") +
  ggtitle(paste(i))


load_visium <- Seurat::Load10X_Spatial(data.dir   = "/home/owly/Desktop/analysis/visium/visium/folder_visium/synovial_sarcoma/NYESO_181",
                                       filename   = "filtered_feature_bc_matrix.h5",
                                       assay      = "Spatial",
                                       slice      = "NYESO181")

# violin plot of ncounts per spot
plot1 <- VlnPlot(load_visium, 
                 features = colnames(load_visium@meta.data)[2], 
                 pt.size  = 0.1) + NoLegend()
plot2 <- SpatialFeaturePlot(load_visium, 
                            features = colnames(load_visium@meta.data)[2]) + 
  theme(legend.position = "right")
display_p1_p2 <- wrap_plots(plot1, 
                            plot2)

# normalization with STC algorithm which builds regularized negative 
# binomial models of gene expression in order to account for technical artifacts 
# while preserving biological variance
normalized_data <- SCTransform(load_visium, 
                               assay   = "Spatial",
                               verbose = FALSE,
                               variable.features.n = 3000)

# feature plot 
spatial_features_plot <- SpatialFeaturePlot(normalized_data, 
                                            features = c("NTN1","UNC5B" ,"EPCAM","PTPRC","PECAM1","ACTA2"),
                                            alpha = c(0.05, 1))

SpatialFeaturePlot(normalized_data, 
                   features = c("NTN1" ,"EPCAM","VCAN",),
                   alpha = c(0.05, 1))

normalized_data <- RunPCA(normalized_data, 
                          assay   = "SCT", 
                          verbose = FALSE)
normalized_data <- FindNeighbors(normalized_data, 
                                 reduction = "pca", 
                                 dims      = 1:30)
# algorithm 4 is Leiden clustering algorithm (supposed 
# to be better than louvain algorithm)
normalized_data <- FindClusters(normalized_data, 
                                verbose    = FALSE,
                                algorithm  = 4,
                                resolution = 0.2)

normalized_data <- RunUMAP(normalized_data, 
                           reduction = "pca", 
                           dims      = 1:30)


normalized_data <- RunTSNE(normalized_data,
                           reduction = "pca",
                           dims =  1:30)
Idents(normalized_data) <- "seurat_clusters"
p1 <- DimPlot(normalized_data, 
              reduction = "umap", 
              label     = TRUE,
              pt.size = 1.5)
p3 <- SpatialDimPlot(normalized_data, 
                     label      = TRUE,
                     label.size = 3)

p2 <- SpatialDimPlot(normalized_data,
                     group.by = "NTN1_status",
              pt.size = 1.5)
p4 <- p1 + p3 + p2

data <- as.data.frame(t(as.matrix(normalized_data@assays[["SCT"]]@data)))

normalized_data@meta.data$NTN1 <- data[,"NTN1"]

normalized_data@meta.data <- normalized_data@meta.data %>% 
  mutate(NTN1_status = case_when(NTN1 > quantile(normalized_data@meta.data$NTN1,0.05) ~ "High",
                                                                                          TRUE ~ "Low"))
table(normalized_data@meta.data $NTN1_status)
Idents(normalized_data) <- "seurat_clusters"
normalized_data_data_marker <- FindAllMarkers(normalized_data,
                                              test.use = "wilcox",
                                              min.pct = 0.15,
                                              logfc.threshold = 0.15)

normalized_data_data_marker %>%
  group_by(cluster) %>%
  top_n(n = 10, wt = avg_log2FC) -> top10

p_hm <- DoHeatmap(normalized_data, 
                  features = top10$gene) + NoLegend()  +
  scale_fill_gradientn(colors = c("seagreen3", "black", "brown1"))

SpatialDimPlot(normalized_data, crop = TRUE,group.by = "NTN1_status",pt.size.factor = 1.2)

gene.setsfull <- getGeneSets(library = "H")

ES_full <- enrichIt(obj = as.matrix(normalized_data@assays[["SCT"]]@counts),
                    gene.sets = gene.setsfull,
                    groups = 200, cores = 4)

normalized_data <- AddMetaData(normalized_data, ES_full)

gene_mak <- read.table(paste0(getwd(),"/signature_mak_1.csv"),sep=",",header=F)

gene_set_mak <- list("MAK_epithelial"= gene_mak[gene_mak$V2 == "E",]$V1,
                     "MAK_mesenchymal"= gene_mak[gene_mak$V2 == "M",]$V1)

ESmakfull <- enrichIt(obj = as.matrix(normalized_data@assays[["SCT"]]@counts),
                      gene.sets = gene_set_mak,
                      groups = 200, cores = 4 )

normalized_data <- AddMetaData(normalized_data, ESmakfull)


normalized_data@meta.data$score_MAK <- normalized_data@meta.data$MAK_mesenchymal - normalized_data@meta.data$MAK_epithelial

saveRDS(normalized_data,"obj_NY181.RDS")

normalized_data$HALLMARK_EPITHELIAL_MESENCHYMAL_TRANSITION
i <- "HALLMARK_EPITHELIAL_MESENCHYMAL_TRANSITION"
stat_test <- normalized_data@meta.data %>%
  wilcox_test(as.formula(paste(as.character(i), "~", "NTN1_status"))) %>%
  add_significance()
stat_test <- stat_test %>% add_xy_position(x = "NTN1_status")

bxp <- ggboxplot(data = normalized_data@meta.data, 
                 x    = "NTN1_status", 
                 y    = i, 
                 fill = "NTN1_status")
p <- bxp + 
  stat_pvalue_manual(stat_test, label = "p.signif") +
  scale_y_continuous(expand = expansion(mult = c(0.1, 0.1))) +
  ylab(paste(i," expression")) + 
  xlab("") +
  scale_fill_discrete("NTN1 status") +
  ggtitle(paste(i))
