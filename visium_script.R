library(Seurat)
library(SeuratData)
library(ggplot2)
library(patchwork)
library(tidyverse)


main_path <- getwd()

# create a result folder
if(dir.exists(paste0(main_path,"/results")) == FALSE){
  dir.create("results")
}

# list name of files in folder_visium and store path of this folder into a vector that contains those paths
list_files <- list.files(paste0(main_path,"/folder_visium"))

path_to_image <- sapply(1:length(list_files),function(i){paste0(main_path,"/folder_visium/",list_files[i],"/outs")})

result_folder_path <- paste0(getwd(),"/results")

# Loop for that allow to create all plot for each visium
for(i in 1:length(path_to_image)){
  
  sub_result_folder_path <- paste0(result_folder_path,"/",list_files[i])
  
  # check if folder exist if not create it
  if(dir.exists(sub_result_folder_path) == FALSE){
    setwd(result_folder_path)
    dir.create(list_files[i])
    setwd(sub_result_folder_path)
  }else{
    setwd(sub_result_folder_path)
  }
  
  # load visium experience
  load_visium <- Seurat::Load10X_Spatial(data.dir   = path_to_image[i],
                                         filename   = "filtered_feature_bc_matrix.h5",
                                         assay      = paste0("Spatial_",str_sub(list_files[i],1,6)),
                                         slice      = str_sub(list_files[i],-4,-1)
  )
  
  # violin plot of ncounts per spot
  plot1 <- VlnPlot(load_visium, 
                   features = colnames(load_visium@meta.data)[2], 
                   pt.size  = 0.1) + NoLegend()
  plot2 <- SpatialFeaturePlot(load_visium, 
                              features = colnames(load_visium@meta.data)[2]) + 
    theme(legend.position = "right")
  display_p1_p2 <- wrap_plots(plot1, 
                              plot2)
  
  pdf("ncounts_plot.pdf",
      height = 6, 
      width = 16)
  plot(display_p1_p2)
  dev.off()
  
  # normalization with STC algorithm which builds regularized negative 
  # binomial models of gene expression in order to account for technical artifacts 
  # while preserving biological variance
  normalized_data <- SCTransform(load_visium, 
                                 assay   = paste0("Spatial_",str_sub(list_files[i],1,6)), 
                                 verbose = FALSE)
  
  # feature plot 
  spatial_features_plot <- SpatialFeaturePlot(normalized_data, 
                                              features = c("NTN1","UNC5B" ,"EPCAM","PTPRC","PECAM1","ACTA2"),
                                              alpha = c(0.05, 1))
  pdf("spatial_features_plot_genes.pdf", 
      height = 8, 
      width = 16)
  plot(spatial_features_plot)
  dev.off()
  
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
                                  resolution = 0.5)
  
  normalized_data <- RunUMAP(normalized_data, 
                             reduction = "pca", 
                             dims      = 1:30)
  

  p1 <- DimPlot(normalized_data, 
                reduction = "umap", 
                label     = TRUE,
                pt.size = 2)
  p2 <- SpatialDimPlot(normalized_data, 
                       label      = TRUE, 
                       label.size = 3)
  p3 <- p1 + p2
  
  pdf("spatial_features_clustering.pdf", height = 8, width = 16)
  plot(p3)
  dev.off()
  
  # Characterizing each cluster with his top 10 most expressed markers
  normalized_data_marker <- FindAllMarkers(normalized_data,
                 test.use = "wilcox",
                 min.pct = 0.5,
                 logfc.threshold = 0.7)
  
  write.csv(normalized_data_marker,file = "markers_per_clusters.csv")
  
  normalized_data_marker %>%
    group_by(cluster) %>%
    top_n(n = 10, wt = avg_log2FC) -> top10
  p_hm <- DoHeatmap(normalized_data, features = top10$gene) + NoLegend()  + 
    scale_fill_gradientn(colors = c("brown1", "black", "green"))
  
  pdf("heatmap_top10_marker.pdf",width = 16,height = 12)
  plot(p_hm)
  dev.off()
  
  setwd(result_folder_path)
}



# multiple slices ---------------------------------------------------------

setwd(main_path)

result_folder_path <- paste0(getwd(),"/results")

setwd(result_folder_path)

result_folder_path_2 <- paste0(result_folder_path,"/multiple_slices_results")
if(dir.exists(result_folder_path_2) == FALSE){
  dir.create("multiple_slices_results")
}

setwd(main_path)
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

result_folder_path_2 <- paste0(result_folder_path,"/multiple_slices_results")

# Loop for that allow to create all plot for each visium
for(i in 1:length(path_list_to_file_patient)){
  
  sub_result_folder_path <- paste0(result_folder_path_2,"/",list_patient_filename[["ID_patient"]][i])
  
  # check if folder exist if not create it
  if(dir.exists(sub_result_folder_path) == FALSE){
    setwd(result_folder_path_2)
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
                                             slice      = "not_treated")
    
    load_visium_2 <- Seurat::Load10X_Spatial(data.dir   = path_list_to_file_patient[[i]][2],
                                             filename   = "filtered_feature_bc_matrix.h5",
                                             assay      = "Spatial",
                                             slice      = "treated")
  
  
  normalized_data_1 <- SCTransform(load_visium_1, assay = "Spatial", verbose = FALSE)
  normalized_data_2 <- SCTransform(load_visium_2, assay = "Spatial", verbose = FALSE)
  
  patient_merge_data <- merge(normalized_data_1, normalized_data_2)
  
  # feature plot 
  spatial_features_plot <- SpatialFeaturePlot(patient_merge_data, 
                                              features = c("NTN1","UNC5B" ,"EPCAM","PTPRC","PECAM1","ACTA2"),
                                              alpha = c(0.05, 1))
  pdf("spatial_features_plot_genes.pdf", 
      height = 8, 
      width = 16)
  plot(spatial_features_plot)
  dev.off()
  
  normalized_data <- RunPCA(patient_merge_data, 
                            assay   = "SCT", 
                            verbose = FALSE)
  normalized_data <- FindNeighbors(patient_merge_data, 
                                   reduction = "pca", 
                                   dims      = 1:30)
  # algorithm 4 is Leiden clustering algorithm (supposed 
  # to be better than louvain algorithm)
  normalized_data <- FindClusters(patient_merge_data, 
                                  verbose    = FALSE,
                                  algorithm  = 4,
                                  resolution = 0.5)
  
  normalized_data <- RunUMAP(patient_merge_data, 
                             reduction = "pca", 
                             dims      = 1:30)
  
  
  p1 <- DimPlot(patient_merge_data, 
                reduction = "umap", 
                label     = TRUE,
                pt.size = 2)
  p2 <- SpatialDimPlot(patient_merge_data, 
                       label      = TRUE, 
                       label.size = 3)
  p3 <- p1 + p2
  
  pdf("spatial_features_clustering.pdf", height = 8, width = 16)
  plot(p3)
  dev.off()
  
  # Characterizing each cluster with his top 10 most expressed markers
  patient_merge_data_marker <- FindAllMarkers(patient_merge_data,
                                           test.use = "wilcox",
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
  
  }
  setwd(result_folder_path_2)
}
