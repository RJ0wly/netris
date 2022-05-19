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
               "rstatix",
               "ggplot2"),
             require, character.only=T)

set.seed(500)

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
                                             slice      = "bftrt_")
    
    load_visium_1@meta.data <- load_visium_1@meta.data %>% 
      mutate(orig.ident = case_when(orig.ident == "SeuratProject" ~ "before treatment"))
    
    load_visium_2 <- Seurat::Load10X_Spatial(data.dir   = path_list_to_file_patient[[i]][2],
                                             filename   = "filtered_feature_bc_matrix.h5",
                                             assay      = "Spatial",
                                             slice      = "aftrt_")
    
    load_visium_2@meta.data <- load_visium_2@meta.data %>% 
      mutate(orig.ident = case_when(orig.ident == "SeuratProject" ~ "after treatment"))
    
    patient_merge_data <- merge(load_visium_1, load_visium_2)
    normalized_data_1 <- SCTransform(load_visium_1, assay = "Spatial", verbose = FALSE)
    normalized_data_2 <- SCTransform(load_visium_2, assay = "Spatial", verbose = FALSE)
    
    patient_merge_data <- subset(patient_merge_data, subset = nFeature_Spatial > 500)
    
    patient_merge_data@meta.data$orig.ident <- factor(patient_merge_data@meta.data$orig.ident, 
                                                      levels =c("before treatment","after treatment"))
    
    patient_merge_data <- SCTransform(patient_merge_data, assay = "Spatial", verbose = FALSE)
    
    patient_merge_data@meta.data$orig.ident <- factor(patient_merge_data@meta.data$orig.ident, 
                                                      levels =c("before treatment","after treatment"))
    
    patient_merge_data@meta.data <- patient_merge_data@meta.data %>% 
      mutate(orig.ident = case_when(orig.ident == "after treatment" ~ "after treatment",
                                    TRUE ~ "before treatment"))
    
    patient_merge_data@meta.data$orig.ident <- factor(patient_merge_data@meta.data$orig.ident, 
                                                      levels =c("before treatment","after treatment"))
    
    # feature plot
    spatial_features_plot <- SpatialPlot(patient_merge_data,
                                                features = c("NTN1","UNC5B" ,"EPCAM","PTPRC","PECAM1","ACTA2"),
                                                alpha = c(0.05, 1),
                                                min.cutoff = 0)
    pdf("spatial_features_plot_genes.pdf",
        height = 20,
        width = 16)
    plot(spatial_features_plot)
    dev.off()
    # 
    patient_merge_data <- RunPCA(patient_merge_data,
                                 assay   = "SCT",
                                 verbose = FALSE)
    
    # Harmony algorithm after Run PCA in order to minimize the batch effect
    patient_merge_data <- patient_merge_data %>% 
      RunHarmony("orig.ident", plot_convergence = F,assay.use ="SCT")
    
    # merge harmony reduction in the seurat object
    harmony_embeddings <- Embeddings(patient_merge_data, 'harmony')
    
    DimPlot(object    = patient_merge_data, 
            reduction = "harmony", 
            pt.size   = 1, 
            group.by  = "orig.ident")
    
    patient_merge_data <- RunUMAP(object    = patient_merge_data, 
                                  reduction = "harmony",
                                  dims      = 1:30)
    patient_merge_data <- FindNeighbors(object    = patient_merge_data,
                                        reduction = "harmony",
                                        dims      = 1:30,
                                        k.param   = 50,
                                        n.trees   = 100)
    
    # algorithm 4 is Leiden clustering algorithm (supposed
    # to be better than louvain algorithm)
    if(i == 1){
      patient_merge_data <- FindClusters(patient_merge_data,
                                         verbose    = FALSE,
                                         algorithm  = 4,
                                         resolution = 0.5)
      
      patient_merge_data@meta.data$ident_clustering <- patient_merge_data@meta.data$SCT_snn_res.0.5
      
      patient_merge_data@meta.data <- patient_merge_data@meta.data %>% 
        mutate("ident_annotation" = case_when(ident_clustering == 1 ~ "CAFs",
                                            ident_clustering == 2 ~ "Tumoral",
                                            ident_clustering == 3 ~ "Tumoral",
                                            ident_clustering == 4 ~ "Vascularized stroma"))
    }else{
    patient_merge_data <- FindClusters(patient_merge_data,
                                       verbose    = FALSE,
                                       algorithm  = 4,
                                       resolution = 0.2)
    
    patient_merge_data@meta.data$ident_clustering <- patient_merge_data@meta.data$SCT_snn_res.0.2
    
    patient_merge_data@meta.data<- patient_merge_data@meta.data %>% 
      mutate("ident_annotation" = case_when(ident_clustering == 1 ~ "Tumoral",
                                          ident_clustering == 2 ~ "Vascularized stroma",
                                          ident_clustering == 3 ~ "CAFs",
                                          ident_clustering == 4 ~ "Vascularized stroma"))
    }
    
    p1 <- SpatialDimPlot(patient_merge_data,
                         label      = TRUE,
                         label.size = 3,
                         group.by = c("ident_annotation"))
    
    p2 <- DimPlot(patient_merge_data, reduction = "harmony", 
            group.by = c("ident_annotation"))
    
    p3 <- DimPlot(patient_merge_data, reduction = "harmony", 
                  group.by = c("orig.ident"))
    
    p4 <- p1 + p2 + p3
    
    pdf("spatial_umap_clustering.pdf", height = 16, width = 24)
    plot(p4)
    dev.off()
    
    Idents(patient_merge_data) <- "ident_annotation"
    KRT_gene_to_check <- rownames(patient_merge_data)[grep("^KRT[0-9]{2}",rownames(patient_merge_data))]
    gene_to_check <- c("EPCAM","CLDN5","CLDN4","DRAXIN","VIM")
    
    ftr_boxplot <- features_boxplot(patient_merge_data,c("Tumoral"),gene_to_check)
    
    pdf("Features_boxplot_stats.pdf",width = 16,height = 8*ceiling(length(gene_to_check)/2))
    plot(ftr_boxplot)
    dev.off()
  
    
    # ftr_boxplot <- features_boxplot(patient_merge_data,c("Tumoral"),KRT_gene_to_check)
    # pdf("Features_boxplot_KRT.pdf",width = 16,height = 8*ceiling(length(KRT_gene_to_check)/2))
    # plot(ftr_boxplot)
    # dev.off()
    
    # prepare multiple SCT object, by regressing the minimum median UMI
    # patient_merge_data <- PrepSCTFindMarkers(patient_merge_data,assay = "SCT")
    # Characterizing each cluster with his top 10 most expressed markers
    
    patient_merge_data_marker <- FindAllMarkers(patient_merge_data,
                                                test.use = "negbinom",
                                                min.pct = 0.5,
                                                logfc.threshold = 0.5)
    
    
    
    write.csv(patient_merge_data_marker,file = "markers_per_clusters.csv")
    
    
    patient_merge_data_marker %>%
      group_by(cluster) %>%
      top_n(n = 10, wt = avg_log2FC) -> top10
    
    p_hm <- DoHeatmap(patient_merge_data, 
                      features = top10$gene) + NoLegend()  +
      scale_fill_gradientn(colors = c("seagreen3", "black", "brown1"))
    
    pdf("heatmap_top10_marker.pdf",width = 16,height = 12)
    plot(p_hm)
    dev.off()
    
    gene.setsfull <- getGeneSets(library = "H")
    
    ES_full <- enrichIt(obj = patient_merge_data@assays[["Spatial"]]@counts,
                        gene.sets = gene.setsfull,
                        groups = 1000, cores = 4)
    
    patient_merge_data <- AddMetaData(patient_merge_data, ES_full)
    
    gene_mak <- read.table(paste0(main_path,"/signature_mak_1.csv"),sep="\t",header=F)
    
    gene_set_mak <- list("MAK_epithelial"= gene_mak[gene_mak$V2 == "E",]$V1,
                         "MAK_mesenchymal"= gene_mak[gene_mak$V2 == "M",]$V1)
    
    ESmakfull <- enrichIt(obj = patient_merge_data@assays[["Spatial"]]@counts,
                          gene.sets = gene_set_mak,
                          groups = 1000, cores = 4 )
    
    patient_merge_data <- AddMetaData(patient_merge_data, ESmakfull)
    
    # if(i == 1){
    #   cluster_to_get = 1
    # }
    # else{
    #   cluster_to_get = 1
    # }
    sub_metrics <- patient_merge_data@meta.data[patient_merge_data@meta.data$ident_annotation == "Tumoral",]
    
    sub_metrics$MAK_score <- sub_metrics$MAK_mesenchymal - sub_metrics$MAK_epithelial
    
    if(nrow(sub_metrics[sub_metrics$orig.ident == "not treated",]) <
       nrow(sub_metrics[sub_metrics$orig.ident == "treated",])){
      bigger_set <- sub_metrics[sub_metrics$orig.ident == "treated",]
      lesser_set <- sub_metrics[sub_metrics$orig.ident == "not treated",]
      
    }else{
      bigger_set <- sub_metrics[sub_metrics$orig.ident == "not treated",]
      lesser_set <- sub_metrics[sub_metrics$orig.ident == "treated",]
    }
    
    NES_list_name <- colnames(sub_metrics)[grep("^HALLMARK_.*",colnames(sub_metrics))]
    
    plot_HMRK <- lapply(NES_list_name,function(iterator){
      
      p <- ggplot(sub_metrics,aes_string(x = "orig.ident",
                                    y = iterator,
                                    fill = "orig.ident")) +
        geom_boxplot() +
        theme_light() +
        ggtitle(paste0("NES :",iterator))
      
      return(p)
    })
    
    splitPlots <- split(plot_HMRK, ceiling(seq_along(plot_HMRK)/2))
    
    pdf(file = "HMRK_plot.pdf", width = 16, height = 8)
    lapply(splitPlots, plotPlots)
    dev.off()
    
    p_HEMT <- ggplot(sub_metrics,aes(x = orig.ident, 
                                     y = HALLMARK_EPITHELIAL_MESENCHYMAL_TRANSITION,
                                     fill = orig.ident)) +
      geom_boxplot() +
      ylab("NES Hallmark EMT") +
      xlab("") +
      scale_fill_discrete("") +
      theme_light() +
      ggtitle("HALLMARK EMT")
    
    p_MAK <- ggplot(sub_metrics,aes(x = orig.ident, 
                                    y = MAK_score,
                                    fill = orig.ident)) +
      geom_boxplot() +
      ylab("NES score MAK") +
      xlab("") +
      scale_fill_discrete("") +
      theme_light() +
      ggtitle("Score MAK")
    
    p_NES <- p_HEMT + p_MAK
    
    pdf("NES_score_boxplot.pdf",width = 16,height = 8)
    plot(p_NES)
    dev.off()
    
    saveRDS(patient_merge_data,paste0("R_obj_visium",list_patient_filename[["ID_patient"]][i],".RDS"))
  }  
  setwd(result_folder_path)
}