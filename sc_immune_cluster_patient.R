library(Seurat)
library(SeuratObject)
library(tidyverse)
library(ggpubr)
library(rstatix)

merge_seurat_obj <- readRDS("/home/owly/Desktop/single_cell_analysis/datasinglecell/tumors_obj_cluster_immune_cluster.RDS")

seurat_obj_list <- Seurat::SplitObject(merge_seurat_obj, split.by = "orig.ident")

features <- SelectIntegrationFeatures(object.list = seurat_obj_list)

seurat_obj_anchors <- FindIntegrationAnchors(object.list = seurat_obj_list, 
                                             anchor.features = features,
                                             reduction = "cca")

seurat_obj_anchors_combined <- IntegrateData(anchorset = seurat_obj_anchors)

rm(seurat_obj_anchors)

DefaultAssay(seurat_obj_anchors_combined) <- "integrated"

seurat_obj_anchors_combined <- ScaleData(seurat_obj_anchors_combined, verbose = FALSE)
seurat_obj_anchors_combined <- RunPCA(seurat_obj_anchors_combined, npcs = 30, verbose = FALSE)

ElbowPlot(seurat_obj_anchors_combined)

seurat_obj_anchors_combined <- RunUMAP(seurat_obj_anchors_combined, reduction = "pca", dims = 1:15,
                                       n.neighbors = 50L,
                                       metric = "euclidean",
                                       min.dist = 0.2,
                                       local.connectivity = 2L)
seurat_obj_anchors_combined <- FindNeighbors(seurat_obj_anchors_combined, reduction = "pca", dims = 1:15)
seurat_obj_anchors_combined <- FindClusters(seurat_obj_anchors_combined, resolution = 0.2)

p1 <- DimPlot(seurat_obj_anchors_combined, reduction = "umap", group.by = "orig.ident")
p2 <- DimPlot(seurat_obj_anchors_combined, reduction = "umap", label = TRUE, repel = TRUE)
p3 <- p1 + p2
plot(p3)


varnames <- rownames(seurat_obj_anchors_combined@assays[["RNA"]]@counts)
varnames[grep("CD",varnames)]

marker_list <- list(
  "Macrophage_M1"             = c("CD86", "CD80", "CCR7", "ITGAX", "CXCL10"),
  "Macrophage_M2"             = c("CD163", "MRC1", "CXCR1", "CXCR2","FCER1A","FCER1G"),
  "LT_cytotoxic"              = c("CD8A","GZMB", "RUNX3","TBX21"),
  "LT_helper"                 = c("CD4","TBX21", "GATA3","RORC"),
  "LT_Reg"                    = c("FOXP3"),
  "PMN_MDSC"                  = c("ITGAM", "FUT4", "NLRP3",varnames[grep("CEACA",varnames)]),
  "NK"                        = c("FCGR3A","CD2", "B3GAT1", "CCR7", "NCAM1", "FGFBP2", "CCL3",varnames[grep("S1P",varnames)],"PRF1","IFNG"),
  "Dendritic_Cells_monocytes" = c("CADM1","CD5")
)

for( i in seq(marker_list)){
  cat(paste0(paste0("plotting ",names(marker_list)[i]),"\n"))
  tmp <- as.vector(marker_list[[i]])
  if(length(tmp) <= 2){
    height_plot   <- 4
    width_plot    <- 7
  }else if(length(tmp) > 5){
    width_plot    <- 10
    height_plot   <- 15
  }else if(length(tmp) > 9){
    width_plot    <- 100
    height_plot   <- 50
  }else{
    width_plot    <- 10
    height_plot   <- 10}
  pdf(file = paste0(paste0(paste0(getwd(),"/results/immune_compartment/umap_immune_plot_analysis_marker_",names(marker_list)[i])),".pdf"),
      width = width_plot,
      height = height_plot)
  plot_features <-suppressWarnings(FeaturePlot(seurat_obj_anchors_combined,
                                               features = tmp,
                                               min.cutoff = 0.01,
                                               max.cutoff = 2,
                                               pt.size = 0.30))
  plot(plot_features)
  dev.off()
}


pdf(file = paste0(paste0(paste0(getwd(),"/results/immune_compartment/umap_immune_plot_analysis_marker_",names(marker_list)[7])),".pdf"),
    width = 15,
    height = 20)
plot_features <-suppressWarnings(FeaturePlot(seurat_obj_anchors_combined,
                                             features = marker_list[[7]],
                                             min.cutoff = 0.01,
                                             pt.size = 0.30))
plot(plot_features)
dev.off()

pdf(file = paste0(paste0(paste0(getwd(),"/results/immune_compartment/umap_immune_plot_analysis_marker_",names(marker_list)[6])),".pdf"),
    width = 15,
    height = 20)
plot_features <-suppressWarnings(FeaturePlot(seurat_obj_anchors_combined,
                                             features = marker_list[[6]],
                                             min.cutoff = 0.01,
                                             pt.size = 0.30))
plot(plot_features)
dev.off()