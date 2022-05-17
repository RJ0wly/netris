features_boxplot <- function(vec_clusters,feature){
  features_cluster <- as.data.frame(patient_merge_data@assays[["SCT"]]@data[feature,WhichCells(patient_merge_data, idents = vec_clusters)])
  
  colnames(features_cluster) <- c(feature)
  
  features_cluster$ID <- rownames(features_cluster)
  
  patient_merge_data@meta.data$ID <- rownames(patient_merge_data@meta.data)
  
  df_compare <- merge(patient_merge_data@meta.data,features_cluster,by="ID")
  
  
  stat_test <- df_compare %>%
    wilcox_test(as.formula(paste(as.character(feature), "~", "orig.ident"))) %>%
    add_significance()
  stat_test <- stat_test %>% add_xy_position(x = "orig.ident")
  
  bxp <- ggboxplot(df_compare, x = "orig.ident", y = feature, fill = "orig.ident")
  p <- bxp + 
    stat_pvalue_manual(stat_test, label = "p.signif") +
    scale_y_continuous(expand = expansion(mult = c(0.05, 0.1))) +
    ylab(paste(feature," expression")) + 
    xlab("Treatment") +
    scale_fill_discrete("") +
    ggtitle(paste(feature," expression in tumoral clusters"))
  
  plot(p)
  return(p)
} 
