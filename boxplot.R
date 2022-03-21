metrics <- cbind(seurat_obj_anchors_combined@meta.data,data)

metrics_nt <- metrics[metrics$orig.ident == "not treated",]

metrics_nt <- metrics_nt[sample(rownames(metrics_nt),size = nrow(metrics_t)),]

metrics_t <- metrics[metrics$orig.ident == "treated",]

metrics <- rbind(metrics_nt,metrics_t)
metrics
my_comparisons <- list( c("not treated", "treated"))


plot_list_boxplot <- lapply(1:length(features_to_boxplot), function(i){
  ggboxplot(metrics, 
                      x ="orig.ident", 
                      y = as.character(features_to_boxplot[i]),
                      fill = "orig.ident")+ 
    stat_compare_means(comparisons = my_comparisons) + 
    stat_compare_means(label.y = max(metrics[,as.character(features_to_boxplot[i])])+1) + 
    xlab(as.character(features_to_boxplot[i])) + 
    ylab("RNA expression") + 
    theme(legend.position = "none")
})

cowplot::plot_grid(plotlist = plot_list_boxplot)
