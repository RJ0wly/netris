library(Seurat)
library(SeuratObject)
library(tidyverse)
library(ggpubr)
library(rstatix)
library(biomaRt)

set.seed(500)

merge_seurat_obj <- readRDS("/home/owly/Desktop/single_cell_analysis/datasinglecell/tumors_cluster_obj__cluster2_3_5.RDS")

seurat_obj_list <- Seurat::SplitObject(merge_seurat_obj, split.by = "orig.ident")

seurat_obj_list <- base::lapply(X = seurat_obj_list, 
                                FUN = function(x) {
                                  x <- Seurat::NormalizeData(x,
                                                             normalization.method = "LogNormalize",
                                                             scale.factor = 10000)}
)

seurat_obj_list <- base::lapply(X = seurat_obj_list, 
                                FUN = function(x) {
                                  x <- Seurat::FindVariableFeatures(x, selection.method = "vst")
                                })

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

seurat_obj_anchors_combined <- RunUMAP(seurat_obj_anchors_combined, reduction = "pca", dims = 1:10,
                                       n.neighbors = 40L,
                                       metric = "euclidean",
                                       min.dist = 0.1,
                                       local.connectivity = 2L)
seurat_obj_anchors_combined <- FindNeighbors(seurat_obj_anchors_combined, reduction = "pca", dims = 1:5)
seurat_obj_anchors_combined <- FindClusters(seurat_obj_anchors_combined, resolution = 0.1)

p1 <- DimPlot(seurat_obj_anchors_combined, reduction = "umap", group.by = "orig.ident")
p2 <- DimPlot(seurat_obj_anchors_combined, reduction = "umap", label = TRUE, repel = TRUE)
p3 <- p1 + p2
plot(p3)

genes_vector <- c("NTN1","UNC5B","NEO1","EPCAM","CLDN7","CD14","CDH1","CLDN4",
  "HOOK1","MUC1","GRHL2","TFF3","CDS1","ESRP1","ESRP2","MARVELD2","F11R","CTNND1",
  "VIM","VCAM1","FBN1","SPARC","FN1","COL6A3","CDH2","ZEB1","ZEB2","THBS2","SULF1",
  "PCOLCE","INHBA","CALD1")

# violin plot
for (i in genes_vector){
  png(file   = paste0(paste0(paste0(getwd(),"/results/result_endo/violon_plot/vln_plot_"),i),".png"),
      width  = 1040,
      height = 555)
  p <- suppressWarnings(VlnPlot(seurat_obj_anchors_combined, features = c(i),
          pt.size = 0.1))
  p <- p + ggtitle(i)
  plot(p)
  dev.off()
}


for (i in genes_vector){
  png(file   = paste0(paste0(paste0(getwd(),"/results/result_endo/feature_plot/feature_plot_"),i),".png"),
      width  = 1040,
      height = 555)
  p <- suppressWarnings(FeaturePlot(seurat_obj_anchors_combined, 
                                features = c(i),
                                pt.size = 0.7,
                                min.cutoff = 0.0,
                                keep.scale = "feature",
                                split.by = "orig.ident"))
  p + theme(legend.positon = "right")
  plot(p)
  dev.off()
}

# pie chart
library(ggrepel)
piechart <- function(data = data, title = "", logical_palette_color = F,palette_color = "Set1"){
  data <- as.data.frame(data)
  df <- data %>%
    arrange(desc(Var1)) %>%
    mutate(prop = Freq / sum(Freq) *100) %>%
    mutate(ypos = cumsum(prop)- 0.4*prop) %>%
    mutate(xpos = cumsum(prop) + 0.5*prop)
  
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

palette_to_color <- c("#F87766D","#ABA300","#00BE67","#00A9FF","#ED68ED")

#boostrap proportion
metrics <- seurat_obj_anchors_combined@meta.data
nt <- metrics[metrics$orig.ident == "not treated",]
t  <- metrics[metrics$orig.ident == "treated",]

res <- lapply(1:100,function(i){
  cat("boostrap",i,"/100\n")
  nt <- metrics[sample(rownames(nt),size= nrow(t),replace=TRUE),]
  table_nt <- table(nt$seurat_clusters)
})

# make dataframe from boostrap result
df <- do.call(rbind, res)
table_nt <- as.table(round(apply(df,2,mean)),0)

#nt <- metrics[metrics$orig.ident == "not treated",]
table_t <- table(t$seurat_clusters)

table_t_plot <- piechart(table_t,
                         "Treated",
                         logical_palette_color = TRUE,
                         palette_color = palette_to_color)
table_nt_plot <- piechart(table_nt,
                          "Not treated",
                          logical_palette_color = TRUE,
                          palette_color = palette_to_color)

png(file   = paste0(paste0(paste0(getwd(),"/results/result_endo/piechart/piechart_bs"),"not_treated"),".png"),
    width  = 1040,
    height = 555)
p <- table_nt_plot[["plot"]]

plot(p)
dev.off()


png(file   = paste0(paste0(paste0(getwd(),"/results/result_endo/piechart/piechart_bs"),"treated"),".png"),
    width  = 1040,
    height = 555)
p <- table_t_plot[["plot"]]
plot(p)
dev.off()

# FIND marker
seurat_obj_anchors_combined_markers <- FindAllMarkers(seurat_obj_anchors_combined, only.pos = TRUE, min.pct = 0.25, logfc.threshold = 0.25)


seurat_obj_anchors_combined_markers %>%
  group_by(cluster) %>%
  top_n(n = 10, wt = avg_log2FC) -> top10

png(file   = paste0(paste0(paste0(getwd(),"/results/result_endo/heatmap/heatmap"),""),".png"),
    width  = 1226,
    height = 704)
p <- DoHeatmap(seurat_obj_anchors_combined, features = top10$gene) + scale_fill_gradientn(colors = c("brown1", "black", "green"))
plot(p)
dev.off()


#print pdf
pdf(file = paste0(paste0(paste0(getwd(),"/results/result_endo/umap/umap_"),"plot3"),".pdf"),
    width = 10,
    height = 6)
plot(p3)
dev.off()

p4 <- DimPlot(seurat_obj_anchors_combined, reduction = "umap", split.by = "orig.ident")

pdf(file = paste0(paste0(paste0(getwd(),"/results/result_endo/umap/umap_"),"plot4"),".pdf"),
    width = 7,
    height = 6)
plot(p4)
dev.off()

genes_vector <- c(genes_vector,"COL1A1","COL1A2","COL3A1","COL5A1","COL5A2","COL6A1","COL6A2","COL10A1")
data <- seurat_obj_anchors_combined@assays[["RNA"]]@data
genes_vector <- intersect(rownames(data),genes_vector)
data <- as.data.frame(data[genes_vector,])
data_t <- as.data.frame(t(data[,rownames(metrics[metrics$orig.ident == "treated",])]))
data_nt <- as.data.frame(t(data[,rownames(metrics[metrics$orig.ident == "not treated",])]))

bs   <- 1000
res1 <- lapply(1:bs,
                        function(i){
  cat("boostrap",i,"/",as.character(bs),"\n")
  tpm_data_nt <- data_nt[sample(
                                rownames(data_nt),
                                size    = nrow(data_t),
                                replace = TRUE
                                ),]
    stat_test <- sapply(1:length(genes_vector),
                        function(j){
      wilcox.test(
                  tpm_data_nt[,j],
                  data_t[,j]
                 )[["p.value"]]
          })
      
      val_var <- sapply(1:length(genes_vector),
                        function(j){
        var(
          tpm_data_nt[,j]
          )
        })
      
      val_sd <- sapply(1:length(genes_vector),
                       function(j){
        sd(
          tpm_data_nt[,j]
           )
        })
    
    tmp_mean <- apply(
                      tpm_data_nt,
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
                    data_t[,i]
                   )})
val_var <- sapply(1:length(genes_vector),
                 function(i){
                   var(
                     data_t[,i]
                   )})
                   


#boostrap wilcox result
df_stats_test <- do.call(rbind, sub_res2)
df_var <- do.call(rbind,sub_res3)
df_std <- do.call(rbind,sub_res4)

mean_bs_wilcox_test_p_value <- apply(df_stats_test,2,mean)
mean_var <- apply(df_var,2,mean)
mean_std <- apply(df_std,2,mean)


df_bs_mean <- do.call(rbind, sub_res1)
mean_nt <- apply(df_bs_mean,2,mean)
mean_t <- apply(data_t,2,mean)
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
write_excel_csv(df_compare,"boostrap_signature_enrichie.csv")

print(df_compare[df_compare$`wilcox test (p value)` < 0.05,]$SYMBOL)

df_compare_mean <- apply(df_compare[,2:3],1,as.numeric)
mean_genes <- apply(df_compare_mean,2,mean)

test <- lapply(1:1000,function(i){
  tpm_data_nt <- data_nt[sample(
    rownames(data_nt),
    size    = nrow(data_t),
    replace = TRUE
  ),]
  
  val_percent <- sapply(1:length(mean_genes),
                  function(j){
                    length(which(tpm_data_nt[,j] > mean_genes[j]))*100/nrow(tpm_data_nt)
                    })
  return(val_percent)
})

df_percent_nt <- do.call(rbind, test)
mean_percent <- apply(df_percent_nt,2,mean)

percent_t <- lapply(1:length(mean_genes),function(i){
  percent_t <- length(which(data_t[,i] > mean_genes[i]))*100/nrow(data_t)
})

df_percent_t <- as.data.frame(do.call(rbind, percent_t))

df_compare$`not treated cell %` <- mean_percent
df_compare$`treated cell %` <- df_percent_t$V1

write_excel_csv(df_compare,"boostrap_signature_enrichie.csv")


max(data[,"COL1A1"]) - quantile(data[,"COL1A1"],0.95)


for(i in genes_vector[grep("COL",genes_vector)]){
png(file = paste0(paste0(paste0(getwd(),"/results/result_endo/feature_plot/features_v2_"),i),".png"),
      width = 1040,
      height = 555)
  p <- suppressWarnings(FeaturePlot(seurat_obj_anchors_combined, 
                               features = c(i),
                               pt.size = 1.2,
                               min.cutoff = 0,
                               max.cutoff = max(data[,i]) - quantile(data[,i],0.95),
                               keep.scale = "feature",
                               split.by = "orig.ident",
                               cols = c("grey", "red")) + theme(legend.position ="none"))
  plot(p)
  dev.off()
}




data <- seurat_obj_anchors_combined@assays[["RNA"]]@data
data <- data %>% as.data.frame() %>% t() %>% as.data.frame()

cluster.markers <- FindMarkers(seurat_obj_anchors_combined, ident.1 = 1,ident.2 = 0, min.pct = 0.25)
cluster.markers <- cluster.markers[order(cluster.markers[,"avg_log2FC"],decreasing = TRUE),]
cluster.markers <- cluster.markers[abs(cluster.markers$avg_log2FC) > 1,]

rownames(cluster.markers[grepl("COL",rownames(cluster.markers)),])

for(i in rownames(cluster.markers)){
png(file = paste0(paste0(paste0(getwd(),"/results/result_endo/DEG_feature_plot/features_plot_"),i),".png"),
      width = 1040,
      height = 555)
p <- suppressWarnings(FeaturePlot(seurat_obj_anchors_combined, 
                                  features = c(i),
                                  pt.size = 1.2,
                                  min.cutoff = 0,
                                  max.cutoff = max(data[,i]) - quantile(data[,i],0.25),
                                  keep.scale = "feature",
                                  split.by = "orig.ident",
                                  cols = c("grey", "red")) + theme(legend.position ="none"))
plot(p)
dev.off()
}

i <- "COL23A1"
suppressWarnings(FeaturePlot(seurat_obj_anchors_combined, 
                                  features = c(i),
                                  pt.size = 1.2,
                                  min.cutoff = 0,
                                  max.cutoff = max(data[,i]) - quantile(data[,i],0.25),
                                  keep.scale = "feature",
                                  split.by = "orig.ident",
                                  cols = c("grey", "red")) + theme(legend.position ="none"))


endo_mak <- read.table("/home/owly/Desktop/single_cell_analysis/single_cell_mouse/EMT_signature_endo.txt")

for(i in intersect(endo_mak$V1,colnames(data))){
  png(file = paste0(paste0(paste0(getwd(),"/results/result_endo/endo_mak_feature_plot/features_plot_"),i),".png"),
      width = 1040,
      height = 555)
  p <- suppressWarnings(FeaturePlot(seurat_obj_anchors_combined, 
                                    features = c(i),
                                    pt.size = 1.2,
                                    min.cutoff = 0,
                                    max.cutoff = max(data[,i]) - quantile(data[,i],0.25),
                                    keep.scale = NULL,
                                    split.by = "orig.ident",
                                    cols = c("grey", "red")) + theme(legend.position ="none"))
  plot(p)
  dev.off()
}

genes_vector <- c("NTN1","UNC5B","NEO1","EPCAM","CLDN7","CD14","CDH1","CLDN4",
                  "HOOK1","MUC1","GRHL2","TFF3","CDS1","ESRP1","ESRP2","MARVELD2","F11R","CTNND1",
                  "VIM","VCAM1","FBN1","SPARC","FN1","COL6A3","CDH2","ZEB1","ZEB2","THBS2","SULF1",
                  "PCOLCE","INHBA","CALD1")
genes_vector <- c(genes_vector,"COL1A1","COL1A2","COL3A1","COL5A1","COL5A2","COL6A1","COL6A2","COL10A1")
data <- seurat_obj_anchors_combined@assays[["RNA"]]@data
genes_vector <- intersect(rownames(data),genes_vector)
genes_vector<- c(intersect(endo_mak$V1,rownames(data)),genes_vector)
data <- as.data.frame(data[genes_vector,])
data_t <- as.data.frame(t(data[,rownames(metrics[metrics$orig.ident == "treated",])]))
data_nt <- as.data.frame(t(data[,rownames(metrics[metrics$orig.ident == "not treated",])]))

bs   <- 1000
res1 <- lapply(1:bs,
               function(i){
                 cat("boostrap",i,"/",as.character(bs),"\n")
                 tpm_data_nt <- data_nt[sample(
                   rownames(data_nt),
                   size    = nrow(data_t),
                   replace = TRUE
                 ),]
                 stat_test <- sapply(1:length(genes_vector),
                                     function(j){
                                       wilcox.test(
                                         tpm_data_nt[,j],
                                         data_t[,j]
                                       )[["p.value"]]
                                     })
                 
                 val_var <- sapply(1:length(genes_vector),
                                   function(j){
                                     var(
                                       tpm_data_nt[,j]
                                     )
                                   })
                 
                 val_sd <- sapply(1:length(genes_vector),
                                  function(j){
                                    sd(
                                      tpm_data_nt[,j]
                                    )
                                  })
                 
                 tmp_mean <- apply(
                   tpm_data_nt,
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
                     data_t[,i]
                   )})
val_var <- sapply(1:length(genes_vector),
                  function(i){
                    var(
                      data_t[,i]
                    )})



#boostrap wilcox result
df_stats_test <- do.call(rbind, sub_res2)
df_var <- do.call(rbind,sub_res3)
df_std <- do.call(rbind,sub_res4)

mean_bs_wilcox_test_p_value <- apply(df_stats_test,2,mean)
mean_var <- apply(df_var,2,mean)
mean_std <- apply(df_std,2,mean)


df_bs_mean <- do.call(rbind, sub_res1)
mean_nt <- apply(df_bs_mean,2,mean)
mean_t <- apply(data_t,2,mean)
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

test <- lapply(1:1000,function(i){
  tpm_data_nt <- data_nt[sample(
    rownames(data_nt),
    size    = nrow(data_t),
    replace = TRUE
  ),]
  
  val_percent <- sapply(1:length(mean_genes),
                        function(j){
                          length(which(tpm_data_nt[,j] > mean_genes[j]))*100/nrow(tpm_data_nt)
                        })
  return(val_percent)
})

df_percent_nt <- do.call(rbind, test)
mean_percent <- apply(df_percent_nt,2,mean)

percent_t <- lapply(1:length(mean_genes),function(i){
  percent_t <- length(which(data_t[,i] > mean_genes[i]))*100/nrow(data_t)
})

df_percent_t <- as.data.frame(do.call(rbind, percent_t))

df_compare$`not treated cell %` <- mean_percent
df_compare$`treated cell %` <- df_percent_t$V1

keep_gene <- df_compare[df_compare$`wilcox test (p value)` < 0.05,]$SYMBOL
write_excel_csv(df_compare,"boostrap_signature_enrichie_all.csv")

df_compare_2 <- df_compare[df_compare$SYMBOL %in% keep_gene,]
df_compare_2 <- df_compare_2[df_compare_2$treated > 0,]
write_excel_csv(df_compare_2,"boostrap_signature_enrichie_significative.csv")

s_genes <- cc.genes$s.genes
g2m_genes <- cc.genes$g2m.genes

seurat_obj_anchors_combined <- CellCycleScoring(seurat_obj_anchors_combined, s.features = s_genes, g2m.features = g2m_genes, set.ident = TRUE)
seurat_obj_anchors_combined <- RunPCA(seurat_obj_anchors_combined, features = c(s_genes, g2m_genes))
pdf(file = paste0(paste0(paste0(getwd(),"/results/result_endo/umap/umap_"),"score_cycle"),".pdf"),
    width = 7,
    height = 6)
DimPlot(seurat_obj_anchors_combined)
dev.off() 


# filter <- rownames(metrics[metrics$seurat_clusters == 2,])
# # immune cluster try
# # filter <- rownames(metrics[metrics$seurat_clusters == 0 |
# #                              metrics$seurat_clusters == 1 |
# #                              metrics$seurat_clusters == 4 | 
# #                              metrics$seurat_clusters == 6| 
# #                              metrics$seurat_clusters == 8,])
# merge_seurat_obj <- merge_seurat_obj[,filter]
# 
# seurat_obj_list <- Seurat::SplitObject(merge_seurat_obj, split.by = "orig.ident")
# 
# seurat_obj_list <- base::lapply(X = seurat_obj_list, 
#                                 FUN = function(x) {
#                                   x <- Seurat::NormalizeData(x,
#                                                              normalization.method = "LogNormalize",
#                                                              scale.factor = 10000)}
# )
# 
# seurat_obj_list <- base::lapply(X = seurat_obj_list, 
#                                 FUN = function(x) {
#                                   x <- Seurat::FindVariableFeatures(x, selection.method = "vst")
#                                 })
# 
# features <- SelectIntegrationFeatures(object.list = seurat_obj_list)
# 
# seurat_obj_anchors <- FindIntegrationAnchors(object.list = seurat_obj_list, 
#                                              anchor.features = features,
#                                              reduction = "cca")
# 
# seurat_obj_anchors_combined <- IntegrateData(anchorset = seurat_obj_anchors)
# 
# rm(seurat_obj_anchors)
# 
# DefaultAssay(seurat_obj_anchors_combined) <- "integrated"
# 
# seurat_obj_anchors_combined <- ScaleData(seurat_obj_anchors_combined, verbose = FALSE)
# seurat_obj_anchors_combined <- RunPCA(seurat_obj_anchors_combined, npcs = 30, verbose = FALSE)
# 
# ElbowPlot(seurat_obj_anchors_combined)
# 
# seurat_obj_anchors_combined <- RunUMAP(seurat_obj_anchors_combined, reduction = "pca", dims = 1:5,
#                                        n.neighbors = 40L,
#                                        metric = "euclidean",
#                                        min.dist = 0.1,
#                                        local.connectivity = 1L)
# seurat_obj_anchors_combined <- FindNeighbors(seurat_obj_anchors_combined, reduction = "pca", dims = 1:5)
# seurat_obj_anchors_combined <- FindClusters(seurat_obj_anchors_combined, resolution = 0.1)
# 
# p1 <- DimPlot(seurat_obj_anchors_combined, reduction = "umap", group.by = "orig.ident")
# p2 <- DimPlot(seurat_obj_anchors_combined, reduction = "umap", label = TRUE, repel = TRUE)
# p3 <- p1 + p2
# plot(p3)
# 
# DimPlot(seurat_obj_anchors_combined, reduction = "umap", split.by = "orig.ident")
# 
# seurat_obj_anchors_combined_markers <- FindAllMarkers(seurat_obj_anchors_combined, only.pos = TRUE, min.pct = 0.25, logfc.threshold = 0.25)
# 
# 
# seurat_obj_anchors_combined_markers %>%
#   group_by(cluster) %>%
#   top_n(n = 20, wt = avg_log2FC) -> top10
# DoHeatmap(seurat_obj_anchors_combined, features = top10$gene) + scale_fill_gradientn(colors = c("brown1", "black", "green"))
# 
# suppressWarnings(FeaturePlot(seurat_obj_anchors_combined, 
#                              features = c("NTN1"),
#                              pt.size = 2,
#                              min.cutoff = 0,
#                              #max.cutoff = 3,
#                              keep.scale = NULL,
#                              split.by = "orig.ident",
#                              cols = c("grey", "red")) + theme(legend.position ="none"))