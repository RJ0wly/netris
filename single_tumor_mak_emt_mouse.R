base::lapply(c("Seurat",
               "tidyverse",
               "SeuratObject",
               "future",
               "rstatix",
               "singscore",
               "escape",
               "dittoSeq",
               "ggpubr",
               "ggrepel",
               "future",
               "uwot",
               "caret",
               "DESeq2",
               "tidymodels"),
             require, character.only=T)
source("/home/rj/Desktop/single_cell_analysis/datasinglecell/utils.R")


#interactive()

plan("multicore", workers = 10)

options(future.globals.maxSize = 2000 * 1024^2)

set.seed(500)

merge_seurat_obj <- readRDS("/home/rj/Desktop/single_cell_analysis/single_cell_mouse/tumors_cluster_obj.RDS")

seurat_obj_list <- base::lapply(X = merge_seurat_obj, 
                                FUN = function(x) {
                                  x <- Seurat::NormalizeData(object               = x,
                                                             normalization.method = "LogNormalize",
                                                             scale.factor         = 10000)
                                })

seurat_obj_list <- base::lapply(X = seurat_obj_list, 
                                FUN = function(x) {
                                  x <- Seurat::FindVariableFeatures(object           = x, 
                                                                    selection.method = "vst",
                                                                    nfeature         = 5000)
                                })

features <- SelectIntegrationFeatures(object.list = seurat_obj_list,
                                      nfeature    = 5000)

seurat_obj_anchors <- FindIntegrationAnchors(object.list     = seurat_obj_list, 
                                             anchor.features = features,
                                             reduction       = "cca",
                                             n.trees         = 100,
                                             nn.method       = "annoy")

seurat_obj_anchors_combined <- IntegrateData(anchorset = seurat_obj_anchors)

rm(seurat_obj_anchors)

DefaultAssay(seurat_obj_anchors_combined) <- "integrated"

seurat_obj_anchors_combined <- ScaleData(object  = seurat_obj_anchors_combined, 
                                         verbose = FALSE)
seurat_obj_anchors_combined <- RunPCA(object  = seurat_obj_anchors_combined, 
                                      npcs    = 30, 
                                      verbose = FALSE)

ElbowPlot(seurat_obj_anchors_combined)

seurat_obj_anchors_combined <- RunUMAP(object      = seurat_obj_anchors_combined, 
                                       reduction   = "pca", 
                                       dims        = 1:10,
                                       n.neighbors = 500L,
                                       metric      = "euclidean",
                                       min.dist    = 0.1)

seurat_obj_anchors_combined <- FindNeighbors(object    = seurat_obj_anchors_combined, 
                                             reduction = "pca", 
                                             dims      = 1:6)
seurat_obj_anchors_combined <- FindClusters(object     = seurat_obj_anchors_combined, 
                                            resolution = 0.1)

p1 <- DimPlot(object    = seurat_obj_anchors_combined, 
              reduction = "umap", 
              group.by = "orig.ident")
p2 <- DimPlot(object    = seurat_obj_anchors_combined, 
              reduction = "umap", 
              label     = TRUE, 
              repel      = TRUE)
p3 <- p1 + p2
plot(p3)

# Cycle score

s_genes <- cc.genes$s.genes %>% str_to_lower() %>% str_to_title()
g2m_genes <- cc.genes$g2m.genes %>% str_to_lower() %>% str_to_title()


DefaultAssay(seurat_obj_anchors_combined) <- "RNA"

seurat_obj_anchors_combined <- ScaleData(object  = seurat_obj_anchors_combined, 
                                         verbose = FALSE)

seurat_obj_anchors_combined <- suppressWarnings(CellCycleScoring(seurat_obj_anchors_combined, 
                                                                 s.features = s_genes, 
                                                                 g2m.features = g2m_genes, 
                                                                 set.ident = TRUE))
seurat_obj_anchors_combined <- suppressWarnings(
  RunPCA(object = seurat_obj_anchors_combined, 
         features      = c(s_genes,g2m_genes)))

DimPlot(object    = seurat_obj_anchors_combined,
        reduction ="pca")

# if(interactive()){
#   
#   ask.regress.cycle.score <- readline(prompt="Would you like to regress the source of this heterogeneity [Y/N]: ")
#   
#   while(!(ask.regress.cycle.score %in% c("Y","Yes","yes","YES","N","No","no","NO"))){
#     message("Answer should be Y, Yes, yes, YES, N, No, no, NO")
#     ask.regress.cycle.score <- readline(prompt="Would you like to regress the source of this heterogeneity [Y/N]: ")
#   }
# }
# 
# if(ask.regress.cycle.score %in% c("Y","Yes","yes","YES")){
#  message("Regressing out cycle cell heterogeneity \n")

seurat_obj_anchors_combined <- ScaleData(object           = seurat_obj_anchors_combined, 
                                         vars.to.regress = c("S.Score", "G2M.Score"), 
                                         features        = rownames(seurat_obj_anchors_combined),
                                         verbose         = TRUE)
seurat_obj_anchors_combined <- suppressWarnings(
  RunPCA(object   = seurat_obj_anchors_combined, 
         features = c(s_genes, g2m_genes)))

DimPlot(object    = seurat_obj_anchors_combined,
        reduction = "pca")
# }else{
#     message("Next step of the analysis")
# }


DefaultAssay(seurat_obj_anchors_combined) <- "integrated"

seurat_obj_anchors_combined <- RunUMAP(object      = seurat_obj_anchors_combined, 
                                       reduction   = "pca", 
                                       dims        = 1:10,
                                       n.neighbors = 40L,
                                       metric      = "euclidean",
                                       min.dist    = 0.1)

seurat_obj_anchors_combined <- FindNeighbors(object    = seurat_obj_anchors_combined, 
                                             reduction = "pca", 
                                             dims      = 1:6)
seurat_obj_anchors_combined <- FindClusters(object     = seurat_obj_anchors_combined, 
                                            resolution = 0.1)
p1 <- DimPlot(object    = seurat_obj_anchors_combined, 
              reduction = "umap", 
              group.by  = "orig.ident")
p2 <- DimPlot(object    = seurat_obj_anchors_combined, 
              reduction = "umap", 
              label     = TRUE, 
              repel     = TRUE)
p3 <- DimPlot(object    = seurat_obj_anchors_combined, 
              reduction = "umap", 
              group.by  = "Phase")
p4 <- p1 + p2 + p3
plot(p4)


metrics <- seurat_obj_anchors_combined@meta.data

endo_genes_emt <- read.table("/home/rj/Desktop/single_cell_analysis/datasinglecell/signature_mak.txt",
                             header = T)

epi <- endo_genes_emt[endo_genes_emt$gene_status == "E",]$SYMBOL %>% str_to_lower() %>% str_to_title()
mes <- endo_genes_emt[endo_genes_emt$gene_status == "M",]$SYMBOL %>% str_to_lower() %>% str_to_title()

gene.sets <- list(
  "epithelial" = epi,
  "mesenchymal" = mes
)

ES <- enrichIt(obj = seurat_obj_anchors_combined,
               gene.sets = gene.sets,
               groups = 1000, cores = 10)
seurat_obj_anchors_combined <- AddMetaData(seurat_obj_anchors_combined, ES)

seurat_obj_anchors_combined@meta.data$active.idents <- seurat_obj_anchors_combined@active.ident
colorblind_vector <- colorRampPalette(rev(c("#0D0887FF", "#47039FFF",
                                            "#7301A8FF", "#9C179EFF", "#BD3786FF", "#D8576BFF",
                                            "#ED7953FF","#FA9E3BFF", "#FDC926FF", "#F0F921FF")))


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
  geom_point(alpha = 0.7) + 
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
  stat_cor(method = "spearman",show.legend = FALSE) +
  ggtitle("Correlation of Normalized Enrichment Score : Epithelial vs Mesenchymal Genes") +
  scale_color_discrete("Treatment")

metrics$mak_score <- metrics$mesenchymal - metrics$epithelial

metrics_nt <- metrics[metrics$orig.ident == "not treated",]
metrics_t <- metrics[metrics$orig.ident == "treated",]

bs_res_list <- boostrap_sampling(nb.sampling       = 100, 
                                 MARGIN            = 1,
                                 dataset.to.sample = metrics_t,
                                 size              = nrow(metrics_nt),
                                 logical.replace   = TRUE)


epi <- do.call(rbind, boostrap_comparaison_function(bs_res_list,metrics_t,"epithelial","wilcox.test","p.value")) %>% mean()
mes <- do.call(rbind, boostrap_comparaison_function(bs_res_list,metrics_t,"mesenchymal","wilcox.test","p.value")) %>% mean()
mak <- do.call(rbind, boostrap_comparaison_function(bs_res_list,metrics_t,"mak_score","wilcox.test","p.value")) %>% mean()

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

res_proportion <- lapply(1:length(bs_res_list),function(i){
  tmp_df <- bs_res_list[[i]]
  table_tmp_df <- table(tmp_df$Phase)
})

res_proportion_df <- do.call(rbind, res_proportion)

table_t <- as.table(round(apply(res_proportion_df,2,mean)),0)
table_nt <- table(metrics_nt$Phase)

palette_to_color <- c("#F87766D","#7CAE00","#00BFC4")

table_t_plot <- piechart(data                  = table_t,
                         title                 = "Treated",
                         logical.palette.color = TRUE,
                         palette.color         = palette_to_color)

table_nt_plot <- piechart(data                  = table_nt,
                          title                 = "Not treated",
                          logical.palette.color = TRUE,
                          palette.color         = palette_to_color)
table_nt_plot[["plot"]]
table_t_plot[["plot"]]

Result=cbind(table_nt,table_t)
res.prop <- chisq.test(Result)
print(res.prop)

nt <- cbind(Result[,1],rownames(Result))
t <- cbind(Result[,2],rownames(Result))
test <- rbind(nt,t)
test <- cbind(test,c("not treated","not treated","not treated","treated","treated","treated"))
colnames(test) <- c("freq","phase","treatment")
test <- as.data.frame(test)
test$freq <- as.numeric(test$freq)

res.pairwise.prop.test <- pairwise.prop.test(Result) %>%
  tidy()
# res$xmin <- c(1,1,2)
# res$xmax <- c(2,3,3)
# res$ypos <- c(2*750,2*800,2*216)
# res$signif <- sapply(res$p.value,function(x){
#   if(as.numeric(x) < 0.05 & as.numeric(x) > 0.01){
#     tmp_sig <- "*"
#   }else if(as.numeric(x) < 0.01 & as.numeric(x) > 0.001){
#     tmp_sig <- "**"
#   }else if(as.numeric(x) < 0.001){
#     tmp_sig <- "***"
#   }else{
#     tmp_sig <- "ns"
#   }
# })
# 
# ggplot(test, aes(x=phase, y=freq)) +
#   geom_bar(stat="identity",aes(fill=treatment)) +
#   geom_signif(
#     xmin = as.integer(res$xmin),
#     xmax = as.integer(res$xmax),
#     y_position = 10 + res$ypos,
#     annotation = format(res$signif),
#     tip_length = 0.05)



data <- as.data.frame(seurat_obj_anchors_combined@assays[["RNA"]]@data)
emt_genes <- endo_genes_emt$SYMBOL %>% str_to_lower() %>% str_to_title()
sub_data <- as.data.frame(t(data[emt_genes,]))
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

ggplot(umap_data, aes(x=V1,y=V2,color = metrics$orig.ident)) + geom_point(size=2) + 
  theme_light() +
  xlab("PC1") +
  ylab("PC2") +
  scale_color_discrete("Treatment")

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


dds <- DESeqDataSetFromMatrix(countData=data, 
                              colData=meta_data_deseq, 
                              design=~status_EMT)
dds <- DESeq(dds)
res <- results(dds)
result_des <- results(dds, tidy=TRUE)
top_result_des <- result_des[abs(result_des$log2FoldChange) > 1.5 & result_des$padj < 0.01,]


library(glmnet)
y <- meta_data_deseq$status_EMT
x <- as.data.frame(t(data[top_result_des$row,]))

dataset <- cbind(x,"status"=y)
## 75% of the sample size
dataset_split <- initial_split(dataset, prop = 0.70, strata = status)
data_train <- training(dataset_split)
x_train <- data_train[,1:848]
y_train <- data_train[,849]
data_test  <-  testing(dataset_split)
x_test <- data_test[,1:848]
y_test <- data_test[,849]


fit <- glmnet(x=x_train, y=y_train,family="binomial")

cv_fit <- cv.glmnet(as.matrix(x_train),
                    y_train,
                    alpha = 1, 
                    nfolds = 5,
                    family="binomial")
tmp_coeffs <- coef(fit, s = cv_fit$lambda.min)

best_lam <- cv_fit$lambda.min

variable_selected = tmp_coeffs@Dimnames[[1]][tmp_coeffs@i + 1]
variable_selected <- variable_selected[-1]


lasso_best <- glmnet(x_train, y_train, alpha = 1, lambda = best_lam,family ="binomial")
pred <- predict(lasso_best, s = best_lam, newx = as.matrix(x_test),type="class")

xtab <- table(y_test,as.factor(pred))

obj_caret_result <- caret::confusionMatrix(xtab, cutoff = 0.5,mode="everything")

data <- as.data.frame(seurat_obj_anchors_combined@assays[["RNA"]]@data)
sub_data <- as.data.frame(t(data[variable_selected,]))

umap_data <- uwot::umap(X = sub_data,
                        n_neighbors = 500,
                        metric = "euclidean",
                        scale = "none",
                        n_trees = 50,
                        min_dist = 0.01,
                        n_epochs = 500,
                        init = "pca",
                        n_components = 10)

umap_data  <- as.data.frame(umap_data)

ggplot(umap_data, aes(x=V2,y=V3,color = metrics$orig.ident)) + geom_point(size=2) + theme_light() +xlab("PC1") +ylab("PC2")

metrics <- metrics %>% mutate(status_EMT = case_when(mak_score<=quantile(metrics$mak_score,0.20) ~ "epithelial-like", 
                                                     mak_score>=quantile(metrics$mak_score,0.80) ~ "mesenchymal-like", 
                                                     TRUE ~ "intermediate"))
ggplot(umap_data, aes(x=V1,y=V2,color = metrics$status_EMT)) + geom_point(size=1,alpha=0.7) + 
  theme_light() +
  xlab("PC1") +
  ylab("PC2") + 
  scale_color_discrete("EMT status")

pca_rec <- recipe(~., data = sub_data) %>%
  step_pca(all_predictors(),threshold = .99)

pca_prep <- prep(pca_rec)
pca_data <- bake(pca_prep, sub_data)

ggplot(pca_data, aes(x=PC01,y=PC02,color = metrics$status_EMT)) + geom_point(size=1,alpha=0.7) + 
  theme_light() +
  xlab("PC1") +
  ylab("PC2") + 
  scale_color_discrete("EMT status")

write.table(variable_selected,"selected_genes_emt_deg_mak_signature_mouse.txt",quote=FALSE,row.names = FALSE,col.names = FALSE)


metric_nt <- metrics[metrics$orig.ident == "not treated",]
metric_t <- metrics[metrics$orig.ident == "treated",]

res <- lapply(1:100,function(i){
  tpm_data_t <- metric_t[sample(
    rownames(metric_t),
    size    = nrow(metric_nt),
    replace = TRUE
  ),]
  table_t <- table(tpm_data_t$status_EMT)
})

df <- do.call(rbind, res)
table_t <- as.table(round(apply(df,2,mean)),0)

table_nt <- table(metric_nt[metric_nt$orig.ident == "not treated",]$status_EMT)

Result=cbind(table_nt,table_t)

res_prop <- prop.test(Result)
res_pair_prop <- pairwise.prop.test(Result) %>% tidy()

table_t_plot <- piechart(table_t,
                         "Treated",
                         logical.palette.color = TRUE,
                         palette.color = palette_to_color)
table_nt_plot <- piechart(table_nt,
                          "Not treated",
                          logical.palette.color = TRUE,
                          palette.color = palette_to_color)

table_nt_plot[["plot"]]
table_t_plot[["plot"]]
