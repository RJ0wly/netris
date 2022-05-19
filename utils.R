features_boxplot <- function(object,vec_clusters,feature){
  
  intersect_features <- intersect(rownames(object@assays[["SCT"]]@data[,]),feature)
  
  message(paste0(c("Features :",intersect_features),collapse=" "))

  if(length(feature) < 2){
  features_cluster <- as.data.frame(object@assays[["SCT"]]@data[intersect_features,
                                                                WhichCells(object, 
                                                                           idents = vec_clusters)])
  }else{
    features_cluster <- as.data.frame(object@assays[["SCT"]]@data[intersect_features,
                                                                  WhichCells(object, 
                                                                             idents = vec_clusters)])
    features_cluster <- features_cluster %>% t() %>% as.data.frame()
  }
  colnames(features_cluster) <- c(intersect_features)
  
  features_cluster$ID <- rownames(features_cluster)
  
  object@meta.data$ID <- rownames(object@meta.data)
  
  obj_df_stats_test <- merge(object@meta.data,features_cluster,by="ID")
  
  list_boxplot <- lapply(intersect_features,function(i){
      stat_test <- obj_df_stats_test %>%
        wilcox_test(as.formula(paste(as.character(i), "~", "orig.ident"))) %>%
        add_significance()
      stat_test <- stat_test %>% add_xy_position(x = "orig.ident")
      
      bxp <- ggboxplot(data = obj_df_stats_test, 
                       x    = "orig.ident", 
                       y    = i, 
                       fill = "orig.ident")
      p <- bxp + 
        stat_pvalue_manual(stat_test, label = "p.signif") +
        scale_y_continuous(expand = expansion(mult = c(0.05, 0.1))) +
        ylab(paste(i," expression")) + 
        xlab("") +
        scale_fill_discrete("") +
        ggtitle(paste(i))
      
      return(p)
    }) 
  ggpubr::ggarrange(
                    plotlist = list_boxplot,
                    ncol     = 2,
                    nrow     = ceiling(length(intersect_features)/2)
                    )
}


plotPlots <- function(x){
  out <- cowplot::plot_grid(plotlist = x, ncol = 2, nrow = 1)
  plot(out)
}
