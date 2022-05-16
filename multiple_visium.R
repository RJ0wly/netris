library(Seurat)
library(SeuratData)
library(ggplot2)
library(patchwork)
library(tidyverse)
library(ggpubr)
library(rstatix)
library(RColorBrewer)
main_path <- getwd()

source(paste0(main_path,"/utils.R"))

# create a result folder
if(dir.exists(paste0(main_path,"/results")) == FALSE){
  dir.create("results")
  
}

if(dir.exists(paste0(main_path,"/results/multiple_slices_results")) == FALSE){
  setwd(paste0(main_path,"/results"))
  dir.create("multiple_slices_results")
  setwd(main_path)
}

# list name of files in folder_visium and store path of this folder into a vector that contains those paths
list_files <- list.files(paste0(main_path,"/folder_visium"))

# extract names of patients and ID of treatment and put in a list
list_patient_filename <- list("ID_patient"   = unique(str_sub(list_files,1,6)),
                              "ID_treatment" = unique(str_sub(list_files,-4,-1)))

# create a list with the path of each visium files per patients
path_list_to_file_patient <- list()
for(i in 1:length(list_patient_filename[["ID_patient"]])){
  for(j in 1:length(list_patient_filename[["ID_treatment"]])){
    tmp_file_name <- paste0(list_patient_filename[["ID_patient"]][i],
                            "_",
                            list_patient_filename[["ID_treatment"]][j])


    path_list_to_file_patient[[as.character(list_patient_filename[["ID_patient"]][i])]][j] = paste0(main_path,
                                                                                                    "/folder_visium/",
                                                                                                    tmp_file_name,
                                                                                                    "/outs")

  }
}

result_folder_path <- paste0(main_path,"/results","/multiple_slices_results")

# Loop for that allow to create all plot for each visium
for(i in 1:length(path_list_to_file_patient)){

  sub_result_folder_path <- paste0(result_folder_path,"/",list_patient_filename[["ID_patient"]][i])

  # check if folder exist if not create it
  if(dir.exists(sub_result_folder_path) == FALSE){
    setwd(result_folder_path)
    dir.create(list_patient_filename[["ID_patient"]][i])
    setwd(sub_result_folder_path)
  }else{
    setwd(sub_result_folder_path)
  }

  for(j in 1:length(path_list_to_file_patient[i])){
    # load visium experience
    load_visium_1 <- Seurat::Load10X_Spatial(data.dir   = path_list_to_file_patient[[i]][1],
                                             filename   = "filtered_feature_bc_matrix.h5",
                                             assay      = "Spatial",
                                             slice      = "not treated")
    
    load_visium_1@meta.data <- load_visium_1@meta.data %>% 
      mutate(orig.ident = case_when(orig.ident == "SeuratProject" ~ "not treated"))

    load_visium_2 <- Seurat::Load10X_Spatial(data.dir   = path_list_to_file_patient[[i]][2],
                                             filename   = "filtered_feature_bc_matrix.h5",
                                             assay      = "Spatial",
                                             slice      = "treated")
    
    load_visium_2@meta.data <- load_visium_2@meta.data %>% 
      mutate(orig.ident = case_when(orig.ident == "SeuratProject" ~ "treated"))


    normalized_data_1 <- SCTransform(load_visium_1, assay = "Spatial", verbose = FALSE)
    normalized_data_2 <- SCTransform(load_visium_2, assay = "Spatial", verbose = FALSE)

    patient_merge_data <- merge(normalized_data_1, normalized_data_2)
    # feature plot
    spatial_features_plot <- SpatialFeaturePlot(patient_merge_data,
                                                features = c("NTN1","UNC5B" ,"EPCAM","PTPRC","PECAM1","ACTA2"),
                                                alpha = c(0.05, 1),
                                                min.cutoff = 0)
    pdf("spatial_features_plot_genes.pdf",
        height = 20,
        width = 16)
    plot(spatial_features_plot)
    dev.off()
    
    
    VariableFeatures(patient_merge_data) <- c(VariableFeatures(normalized_data_1), VariableFeatures(normalized_data_2))
    
    patient_merge_data <- RunPCA(patient_merge_data,
                              assay   = "SCT",
                              verbose = FALSE)
    patient_merge_data <- FindNeighbors(patient_merge_data,
                                     reduction = "pca",
                                     dims      = 1:30)
    # algorithm 4 is Leiden clustering algorithm (supposed
    # to be better than louvain algorithm)
    patient_merge_data <- FindClusters(patient_merge_data,
                                    verbose    = FALSE,
                                    algorithm  = 4,
                                    resolution = 0.2)

    patient_merge_data <- RunUMAP(patient_merge_data,
                               reduction = "pca",
                               dims      = 1:30)

    p1 <- SpatialDimPlot(patient_merge_data,
                         label      = TRUE,
                         label.size = 3)
    p2 <- DimPlot(patient_merge_data, reduction = "umap", group.by = c("ident", "orig.ident"))
    
    p3 <- p1 + p2
    
    pdf("spatial_umap_clustering", height = 16, width = 16)
    plot(p3)
    dev.off()

    # Characterizing each cluster with his top 10 most expressed markers
    patient_merge_data_marker <- FindAllMarkers(patient_merge_data,
                                                test.use = "negbinom",
                                                min.pct = 0.5,
                                                logfc.threshold = 0.7)

    write.csv(patient_merge_data_marker,file = "markers_per_clusters.csv")

    patient_merge_data_marker %>%
      group_by(cluster) %>%
      top_n(n = 10, wt = avg_log2FC) -> top10
    p_hm <- DoHeatmap(patient_merge_data, features = top10$gene) + NoLegend()  +
      scale_fill_gradientn(colors = c("brown1", "black", "green"))

    pdf("heatmap_top10_marker.pdf",width = 16,height = 12)
    plot(p_hm)
    dev.off()
    
    if(i == 1){
      p_epcam <- features_boxplot(vec_clusters = c(2,4),feature = "EPCAM")
      pdf("boxplot_epcam.pdf",width = 8,height = 7)
      plot(p_epcam)
      dev.off()
      
      p_tff3 <- features_boxplot(vec_clusters = c(2,4),feature = "TFF3")
      pdf("boxplot_tff3.pdf",width = 8,height = 7)
      plot(p_tff3)
      dev.off()
      
      p_pgr <- features_boxplot(vec_clusters = c(2,4),feature = "PGR")
      pdf("boxplot_pgr.pdf",width = 8,height = 7)
      plot(p_pgr)
      dev.off()
      
      
    }else if(i == 2){
      p_epcam <- features_boxplot(vec_clusters = c(6,1),feature = "EPCAM")
      pdf("boxplot_epcam.pdf",width = 8,height = 7)
      plot(p_epcam)
      dev.off()
      
      p_tff3 <- features_boxplot(vec_clusters = c(6,1),feature = "TFF3")
      pdf("boxplot_tff3.pdf",width = 8,height = 7)
      plot(p_tff3)
      dev.off()
      
      p_pgr <- features_boxplot(vec_clusters = c(6,1),feature = "PGR")
      pdf("boxplot_pgr.pdf",width = 8,height = 7)
      plot(p_pgr)
      dev.off()
    }
  }
  setwd(result_folder_path)
}


# test zone ---------------------------------------------------------------


# load_visium_1 <- Seurat::Load10X_Spatial(data.dir   = path_list_to_file_patient[[1]][1],
#                                          filename   = "filtered_feature_bc_matrix.h5",
#                                          assay      = "Spatial",
#                                          slice      = "not treated")
# 
# load_visium_1@meta.data <- load_visium_1@meta.data %>% 
#   mutate(orig.ident = case_when(orig.ident == "SeuratProject" ~ "not treated"))
# 
# load_visium_2 <- Seurat::Load10X_Spatial(data.dir   = path_list_to_file_patient[[1]][2],
#                                          filename   = "filtered_feature_bc_matrix.h5",
#                                          assay      = "Spatial",
#                                          slice      = "treated")
# 
# load_visium_2@meta.data <- load_visium_2@meta.data %>% 
#   mutate(orig.ident = case_when(orig.ident == "SeuratProject" ~ "treated"))
# 
# 
# normalized_data_1 <- SCTransform(load_visium_1, assay = "Spatial", verbose = FALSE)
# normalized_data_2 <- SCTransform(load_visium_2, assay = "Spatial", verbose = FALSE)
# 
# patient_merge_data <- merge(normalized_data_1, normalized_data_2)
# # feature plot
# spatial_features_plot <- SpatialFeaturePlot(patient_merge_data,
#                                             features = c("NTN1","UNC5B" ,"EPCAM","PTPRC","PECAM1","ACTA2"),
#                                             alpha = c(0.05, 1),
#                                             min.cutoff = 0)
# pdf("spatial_features_plot_genes.pdf",
#     height = 20,
#     width = 16)
# plot(spatial_features_plot)
# dev.off()
# 
# 
# VariableFeatures(patient_merge_data) <- c(VariableFeatures(normalized_data_1), VariableFeatures(normalized_data_2))
# 
# patient_merge_data <- RunPCA(patient_merge_data,
#                              assay   = "SCT",
#                              verbose = FALSE)
# patient_merge_data <- FindNeighbors(patient_merge_data,
#                                     reduction = "pca",
#                                     dims      = 1:30)
# # algorithm 4 is Leiden clustering algorithm (supposed
# # to be better than louvain algorithm)
# patient_merge_data <- FindClusters(patient_merge_data,
#                                    verbose    = FALSE,
#                                    algorithm  = 4,
#                                    resolution = 0.2)
# 
# patient_merge_data <- RunUMAP(patient_merge_data,
#                               reduction = "pca",
#                               dims      = 1:30)
# 
# p1 <- SpatialDimPlot(patient_merge_data,
#                      label      = TRUE,
#                      label.size = 3)
# p2 <- DimPlot(patient_merge_data, reduction = "umap", group.by = c("ident", "orig.ident"))
# 
# p3 <- FeaturePlot(patient_merge_data,features=("EPCAM"),reduction="umap",cols = rev(brewer.pal(n = 11, name = "RdYlBu")))
# 
# spatial_features_plot <- SpatialFeaturePlot(patient_merge_data,
#                                             features = c("EPCAM","TFF3"),
#                                             alpha = c(0.05, 1),
#                                             min.cutoff = 0)
# 
# p4 <- p2 + p3 + p1 + spatial_features_plot + p_epcam
