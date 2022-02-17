base::lapply(c("Seurat",
               "tidyverse",
               "SeuratObject",
               "future",
               "rstatix",
               "singscore",
               "escape",
               "dittoSeq",
               "ggpubr",
               "ggrepel"),
             require, character.only=T)

set.seed(500)

merge_seurat_obj <- readRDS("/home/rj/Desktop/single_cell_analysis/datasinglecell/tumors_cluster_obj__cluster2_3_5.RDS")

seurat_obj_list <- Seurat::SplitObject(merge_seurat_obj, 
                                       split.by = "orig.ident")

seurat_obj_list <- base::lapply(X = seurat_obj_list, 
                                FUN = function(x) {
                                  x <- Seurat::NormalizeData(x,
                                                             normalization.method = "LogNormalize",
                                                             scale.factor = 10000)}
)

seurat_obj_list <- base::lapply(X = seurat_obj_list, 
                                FUN = function(x) {
                                  x <- Seurat::FindVariableFeatures(x, 
                                                                    selection.method = "vst",
                                                                    nfeature=5000)
                                })

features <- SelectIntegrationFeatures(object.list = seurat_obj_list,
                                      nfeature=5000)

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

ElbowPlot(seurat_obj_anchors_combined)

seurat_obj_anchors_combined <- RunUMAP(seurat_obj_anchors_combined, reduction = "pca", dims = 1:10,
                                       n.neighbors = 40L,
                                       metric = "euclidean",
                                       min.dist = 0.1,
                                       local.connectivity = 1L)

seurat_obj_anchors_combined <- FindNeighbors(seurat_obj_anchors_combined, 
                                             reduction = "pca", 
                                             dims = 1:5)
seurat_obj_anchors_combined <- FindClusters(seurat_obj_anchors_combined, 
                                            resolution = 0.1)

p1 <- DimPlot(seurat_obj_anchors_combined, reduction = "umap", 
              group.by = "orig.ident")
p2 <- DimPlot(seurat_obj_anchors_combined, reduction = "umap", 
              label = TRUE, 
              repel = TRUE)
p3 <- p1 + p2
plot(p3)

s_genes <- cc.genes$s.genes 
g2m_genes <- cc.genes$g2m.genes

seurat_obj_anchors_combined <- CellCycleScoring(seurat_obj_anchors_combined, 
                                                s.features = s_genes, 
                                                g2m.features = g2m_genes, 
                                                set.ident = TRUE)
seurat_obj_anchors_combined <- RunPCA(seurat_obj_anchors_combined, 
                                      features = c(s_genes, 
                                                   g2m_genes))
DimPlot(seurat_obj_anchors_combined,reduction="pca")

seurat_obj_anchors_combined <- ScaleData(seurat_obj_anchors_combined, 
                                         vars.to.regress = c("S.Score", "G2M.Score"), 
                                         features = rownames(seurat_obj_anchors_combined))
seurat_obj_anchors_combined <- RunPCA(seurat_obj_anchors_combined, 
                                      features = c(s_genes, 
                                                   g2m_genes))

DimPlot(seurat_obj_anchors_combined,reduction="pca")

seurat_obj_anchors_combined <- RunUMAP(seurat_obj_anchors_combined, reduction = "pca", dims = 1:10,
                                       n.neighbors = 40L,
                                       metric = "euclidean",
                                       min.dist = 0.1,
                                       local.connectivity = 1L)

seurat_obj_anchors_combined <- FindNeighbors(seurat_obj_anchors_combined, 
                                             reduction = "pca", 
                                             dims = 1:5)
seurat_obj_anchors_combined <- FindClusters(seurat_obj_anchors_combined, 
                                            resolution = 0.1)
p1 <- DimPlot(seurat_obj_anchors_combined, reduction = "umap", 
              group.by = "orig.ident")
p2 <- DimPlot(seurat_obj_anchors_combined, reduction = "umap", 
              label = TRUE, 
              repel = TRUE)
p3 <- DimPlot(seurat_obj_anchors_combined, reduction = "umap", 
              group.by = "Phase")
p4 <- p1 + p2 + p3
plot(p4)

metrics <- seurat_obj_anchors_combined@meta.data

endo_genes_emt <- read.table("/home/rj/Desktop/single_cell_analysis/datasinglecell/signature_mak.txt",
                             header = T)
endo_genes_emt$gene_status
endo_gene_epi <- endo_genes_emt[endo_genes_emt$gene_status == "E",]$SYMBOL
endo_gene_mes <- endo_genes_emt[endo_genes_emt$gene_status == "M",]$SYMBOL


gene.sets <- list(
  "epithelial" = endo_gene_epi,
  "mesenchymal" = endo_gene_mes
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
                x.var = "epithelial",
                y.var = "mesenchymal",
                do.contour = TRUE,
                split.by = "orig.ident") +
  theme_classic() +
  scale_fill_gradientn(colors = rev(colorblind_vector(11))) +
  geom_vline(xintercept = 0, lty=2) +
  geom_hline(yintercept = 0, lty=2)

cor_stats <- cor.test(ES$epithelial,ES$mesenchymal)

metrics <- seurat_obj_anchors_combined@meta.data
ggplot(metrics,aes(x=epithelial,y=mesenchymal,color=orig.ident)) + 
  geom_point(alpha = 0.5) + 
  geom_smooth(method='lm', color="brown3") +
  theme(plot.title = element_text(hjust = 0.5,size=15),
        panel.background = element_rect(fill="white", colour="black", size=0.5, 
                                        linetype="solid"),
        panel.grid.major = element_line(size = 0.5, linetype = 0,
                                        colour = "white"), 
        panel.grid.minor = element_line(size = 0.25, linetype = 'solid',
                                        colour = "azure3"),
        legend.text = element_text(size=15),
        legend.title = element_text(size=20),
        axis.text = element_text(size=12),
        axis.title = element_text(size=15)) +
  stat_cor(method = "pearson") +
  ggtitle("Correlation of Normalized Enrichment Score : Epithelial vs Mesenchymal Genes")


metrics$mak_score <- metrics$mesenchymal - metrics$epithelial

metrics_nt <- metrics[metrics$orig.ident == "not treated",]
metrics_t <- metrics[metrics$orig.ident == "treated",]

bs   <- 100
res1 <- lapply(1:bs,
               function(i){
                 cat("boostrap",i,"/",as.character(bs),"\n")
                 tpm_data_nt <- metrics_nt[sample(
                   rownames(metrics_nt),
                   size    = nrow(metrics_t),
                   replace = TRUE
                 ),]
                 stat_test_NES_epi <- 
                   wilcox.test(
                   tpm_data_nt[,"epithelial"],
                   metrics_t[,"epithelial"]
                 )[["p.value"]]
                                     
                 
                 stat_test_NES_mes<- 
                   wilcox.test(
                     tpm_data_nt[,"mesenchymal"],
                     metrics_t[,"mesenchymal"]
                   )[["p.value"]]
                 
                 stat_test_mak <- 
                   wilcox.test(
                     tpm_data_nt[,"mak_score"],
                     metrics_t[,"mak_score"]
                   )[["p.value"]]
                 
                 return(list(
                   "NES wilcox epi" = stat_test_NES_epi,
                   "NES wilcox mes" = stat_test_NES_mes,
                   "NES wilcox mak" = stat_test_mak
                 ))
               })
                 
                 
sub_res1 <- lapply(1:bs,function(i){
  res1[[i]][["NES wilcox epi"]]
})

sub_res2 <- lapply(1:bs,function(i){
  res1[[i]][["NES wilcox mes"]]
})

sub_res3 <- lapply(1:bs,function(i){
  res1[[i]][["NES wilcox mak"]]
})


epi <- do.call(rbind, sub_res1)
mes <- do.call(rbind,sub_res2)
mak  <- do.call(rbind,sub_res3)

GeomSplitViolin <- ggproto("GeomSplitViolin", GeomViolin, 
                           draw_group = function(self, data, ..., draw_quantiles = NULL) {
                             data <- transform(data, xminv = x - violinwidth * (x - xmin), xmaxv = x + violinwidth * (xmax - x))
                             grp <- data[1, "group"]
                             newdata <- plyr::arrange(transform(data, x = if (grp %% 2 == 1) xminv else xmaxv), if (grp %% 2 == 1) y else -y)
                             newdata <- rbind(newdata[1, ], newdata, newdata[nrow(newdata), ], newdata[1, ])
                             newdata[c(1, nrow(newdata) - 1, nrow(newdata)), "x"] <- round(newdata[1, "x"])
                             
                             if (length(draw_quantiles) > 0 & !scales::zero_range(range(data$y))) {
                               stopifnot(all(draw_quantiles >= 0), all(draw_quantiles <=
                                                                         1))
                               quantiles <- ggplot2:::create_quantile_segment_frame(data, draw_quantiles)
                               aesthetics <- data[rep(1, nrow(quantiles)), setdiff(names(data), c("x", "y")), drop = FALSE]
                               aesthetics$alpha <- rep(1, nrow(quantiles))
                               both <- cbind(quantiles, aesthetics)
                               quantile_grob <- GeomPath$draw_panel(both, ...)
                               ggplot2:::ggname("geom_split_violin", grid::grobTree(GeomPolygon$draw_panel(newdata, ...), quantile_grob))
                             }
                             else {
                               ggplot2:::ggname("geom_split_violin", GeomPolygon$draw_panel(newdata, ...))
                             }
                           })

geom_split_violin <- function(mapping = NULL, data = NULL, stat = "ydensity", position = "identity", ..., 
                              draw_quantiles = NULL, trim = TRUE, scale = "area", na.rm = FALSE, 
                              show.legend = NA, inherit.aes = TRUE) {
  layer(data = data, mapping = mapping, stat = stat, geom = GeomSplitViolin, 
        position = position, show.legend = show.legend, inherit.aes = inherit.aes, 
        params = list(trim = trim, scale = scale, draw_quantiles = draw_quantiles, na.rm = na.rm, ...))
}

ggplot(metrics, aes(x=seurat_clusters, y=epithelial, fill = orig.ident)) + 
  geom_split_violin() + 
  geom_jitter(width=0.25,alpha=0.05) +  
  theme(plot.title = element_text(hjust = 0.5,size=15),
        panel.background = element_rect(fill="white", colour="black", size=0.5, 
                                        linetype="solid"),
        panel.grid.major = element_line(size = 0.5, linetype = 0,
                                        colour = "white"), 
        panel.grid.minor = element_line(size = 0.25, linetype = 'solid',
                                        colour = "azure3"),
        legend.text = element_text(size=15),
        legend.title = element_text(size=20),
        axis.text = element_text(size=12),
        axis.title = element_text(size=15)) +
  labs(x="cluters",y="NES epithelial")

ggplot(metrics, aes(x=seurat_clusters, y=mesenchymal, fill = orig.ident)) + 
  geom_split_violin() + 
  geom_jitter(width=0.25,alpha=0.05) +  
  theme(plot.title = element_text(hjust = 0.5,size=15),
        panel.background = element_rect(fill="white", colour="black", size=0.5, 
                                        linetype="solid"),
        panel.grid.major = element_line(size = 0.5, linetype = 0,
                                        colour = "white"), 
        panel.grid.minor = element_line(size = 0.25, linetype = 'solid',
                                        colour = "azure3"),
        legend.text = element_text(size=15),
        legend.title = element_text(size=20),
        axis.text = element_text(size=12),
        axis.title = element_text(size=15)) +
  labs(x="cluters",y="NES mesenchymal")

metrics$score_mak <- metrics$mesenchymal -metrics$epithelial
ggplot(metrics, aes(x=seurat_clusters, y=mak_score, fill = orig.ident)) + 
  geom_split_violin() + 
  geom_jitter(width=0.25,alpha=0.05) +  
  theme(plot.title = element_text(hjust = 0.5,size=15),
        panel.background = element_rect(fill="white", colour="black", size=0.5, 
                                        linetype="solid"),
        panel.grid.major = element_line(size = 0.5, linetype = 0,
                                        colour = "white"), 
        panel.grid.minor = element_line(size = 0.25, linetype = 'solid',
                                        colour = "azure3"),
        legend.text = element_text(size=15),
        legend.title = element_text(size=20),
        axis.text = element_text(size=12),
        axis.title = element_text(size=15)) +
  labs(x="cluters",y="NES mak score")

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
#boostrap proportion


metrics_nt <- metrics[metrics$orig.ident == "not treated",]
metrics_t <- metrics[metrics$orig.ident == "treated",]

set.seed(400)
res <- lapply(1:100,function(i){
  tpm_data_nt <- metrics_nt[sample(
    rownames(metrics_nt),
    size    = nrow(metrics_t),
    replace = FALSE
  ),]
  table_nt <- table(tpm_data_nt$Phase)
})

# make dataframe from boostrap result
df <- do.call(rbind, res)
table_nt <- as.table(round(apply(df,2,mean)),0)

#nt <- metrics[metrics$orig.ident == "not treated",]
table_t <- table(metrics_t$Phase)

palette_to_color <- c("#F87766D","#7CAE00","#00BFC4","#C77CF77")
table_t_plot <- piechart(table_t,
                         "Treated",
                         logical_palette_color = TRUE,
                         palette_color = palette_to_color)
table_nt_plot <- piechart(table_nt,
                          "Not treated",
                          logical_palette_color = TRUE,
                          palette_color = palette_to_color)
table_nt_plot[["plot"]]
table_t_plot[["plot"]]

Result=cbind(table_nt,table_t)
res_prop <- chisq.test(Result)
print(res_prop)

nt <- cbind(Result[,1],rownames(Result))
t <- cbind(Result[,2],rownames(Result))
test <- rbind(nt,t)
test <- cbind(test,c("not treated","not treated","not treated","treated","treated","treated"))
colnames(test) <- c("freq","phase","treatment")
test <- as.data.frame(test)
test$freq <- as.numeric(test$freq)

res <- pairwise.prop.test(Result) %>%
  tidy()
res$xmin <- c(1,1,2)
res$xmax <- c(2,3,3)
res$ypos <- c(2*750,2*800,2*216)
res$signif <- sapply(res$p.value,function(x){
  if(as.numeric(x) < 0.05 & as.numeric(x) > 0.01){
    tmp_sig <- "*"
  }else if(as.numeric(x) < 0.01 & as.numeric(x) > 0.001){
    tmp_sig <- "**"
  }else if(as.numeric(x) < 0.001){
    tmp_sig <- "***"
  }else{
    tmp_sig <- "ns"
  }
})

# ggplot(test, aes(x=phase, y=freq)) +
#   geom_bar(stat="identity",aes(fill=treatment)) +
#   geom_signif(
#     xmin = as.integer(res$xmin),
#     xmax = as.integer(res$xmax),
#     y_position = 10 + res$ypos,
#     annotation = format(res$signif),
#     tip_length = 0.05)


library(uwot)

data <- as.data.frame(seurat_obj_anchors_combined@assays[["RNA"]]@data)
sub_data <- as.data.frame(t(data[endo_genes_emt$SYMBOL,]))
sub_data <- sub_data[,-grep("NA",colnames(sub_data))]

umap_data <- uwot::umap(X = sub_data,
           n_neighbors = 40,
           metric = "euclidean",
           scale = "none",
           n_trees = 50,
           min_dist = 0.3,
           n_epochs = 500,
           init = "pca")

umap_data  <- as.data.frame(umap_data)

ggplot(umap_data, aes(x=V1,y=V2,color = metrics$orig.ident)) + geom_point(size=2) + theme_light() +xlab("PC1") +ylab("PC2")

#saveRDS(seurat_obj_anchors_combined,"compute_tumoral_mak_cluster.RDS")
#seurat_obj_anchors_combined <- readRDS("compute_tumoral_mak_cluster.RDS")

metrics <- seurat_obj_anchors_combined@meta.data

metrics$mak_score <- metrics$mesenchymal - metrics$epithelial

metrics <- metrics %>% mutate(status_EMT = case_when(mak_score<=quantile(metrics$mak_score,0.10) ~ "Epithelial-like", 
                                          mak_score>=quantile(metrics$mak_score,0.90) ~ "Mesenchymal-like", 
                                          TRUE ~ "intermediate"))

ggplot(umap_data, aes(x=V1,y=V2,color = metrics$status_EMT)) + geom_point(size=1,alpha=0.7) + 
  theme_light() +
  xlab("PC1") +
  ylab("PC2") + 
  scale_color_discrete("EMT status")

meta_data_deseq <- rbind(metrics[metrics$status_EMT == "Epithelial-like",],metrics[metrics$status_EMT == "Mesenchymal-like",])
meta_data_deseq$status_EMT <- as.factor(meta_data_deseq$status_EMT)
data <- as.data.frame(seurat_obj_anchors_combined@assays[["RNA"]]@counts)
data <- data[seurat_obj_anchors_combined@assays[["integrated"]]@var.features,rownames(meta_data_deseq)]
data <- data+1

library(DESeq2)

dds <- DESeqDataSetFromMatrix(countData=data, 
                              colData=meta_data_deseq, 
                              design=~status_EMT)
dds <- DESeq(dds)
res <- results(dds)
result_des <- results(dds, tidy=TRUE)
top_result_des <- result_des[abs(result_des$log2FoldChange) > 1.5 & result_des$padj < 0.01,]


library(glmnet)
y <- meta_data_deseq$status_EMT
x <- data[top_result_des$row,]
fit <- glmnet(x=t(x), y=y,family="binomial")

cvfit <- cv.glmnet(x=t(x), y=y,family="binomial")
plot(cvfit)
tmp_coeffs <- coef(fit, s = cvfit$lambda.min)

variable_selected = tmp_coeffs@Dimnames[[1]][tmp_coeffs@i + 1]
variable_selected <- variable_selected[-1]


data <- as.data.frame(seurat_obj_anchors_combined@assays[["RNA"]]@data)
sub_data <- as.data.frame(t(data[variable_selected,]))
sub_data <- sub_data[,-grep("NA",colnames(sub_data))]

umap_data <- uwot::umap(X = sub_data,
                        n_neighbors = 40,
                        metric = "euclidean",
                        scale = "none",
                        n_trees = 50,
                        min_dist = 0.01,
                        n_epochs = 500,
                        init = "pca")

umap_data  <- as.data.frame(umap_data)

ggplot(umap_data, aes(x=V1,y=V2,color = metrics$orig.ident)) + geom_point(size=2) + theme_light() +xlab("PC1") +ylab("PC2")

metrics <- seurat_obj_anchors_combined@meta.data

metrics$mak_score <- metrics$mesenchymal - metrics$epithelial

metrics <- metrics %>% mutate(status_EMT = case_when(mak_score<=quantile(metrics$mak_score,0.25) ~ "epithelial-like", 
                                                     mak_score>=quantile(metrics$mak_score,0.75) ~ "mesenchymal-like", 
                                                     TRUE ~ "intermediate"))

ggplot(umap_data, aes(x=V1,y=V2,color = metrics$status_EMT)) + geom_point(size=1,alpha=0.7) + 
  theme_light() +
  xlab("PC1") +
  ylab("PC2") + 
  scale_color_discrete("EMT status")

ggplot(umap_data, aes(x=V1,y=V2,color = metrics$orig.ident)) + geom_point(size=1,alpha=0.7) + theme_light() +xlab("PC1") +ylab("PC2") + scale_color_discrete(("Treament"))


metric_nt <- metrics[metrics$orig.ident == "not treated",]
metric_t <- metrics[metrics$orig.ident == "treated",]

res <- lapply(1:100,function(i){
tpm_data_nt <- metric_nt[sample(
  rownames(metrics_nt),
  size    = nrow(metric_t),
  replace = TRUE
),]
table_nt <- table(tpm_data_nt$status_EMT)
})

df <- do.call(rbind, res)
table_nt <- as.table(round(apply(df,2,mean)),0)
                  
table_t <- table(metric_t[metric_t$orig.ident == "treated",]$status_EMT)

Result=cbind(table_nt,table_t)

res_prop <- prop.test(Result)
res_pair_prop <- pairwise.prop.test(Result) %>% tidy()

table_t_plot <- piechart(table_t,
                         "Treated",
                         logical_palette_color = TRUE,
                         palette_color = palette_to_color)
table_nt_plot <- piechart(table_nt,
                          "Not treated",
                          logical_palette_color = TRUE,
                          palette_color = palette_to_color)

table_nt_plot[["plot"]]
table_t_plot[["plot"]]
