piechart <- function(data                  = data, 
                     title                 = "", 
                     logical.palette.color = F,
                     palette.color         = "Set1"){
  
  data <- as.data.frame(data)
  df <- data %>%
    arrange(desc(Var1)) %>%
    mutate(prop = Freq / sum(Freq) *100) %>%
    mutate(ypos = cumsum(prop) - 0.4*prop)
  
  df$prop <- round(df$prop,0)
  
  p <- ggplot(df, aes(x     = "", 
                      y     = prop, 
                      fill  = Var1,
                      label = Var1)) +
    geom_bar(
             stat           = "identity", 
             width          = 1, 
             color          = "white"
      ) +
    coord_polar(theta       = "y", #check if theta is the good argument
                start       = 0) +
    theme_void() +
    guides(fill             = guide_legend(
      title                 = "Cell type")
      )
  if(logical.palette.color == TRUE){
    p + scale_fill_manual(values = palette.color)
  }else{
    p + scale_fill_brewer(palette = palette.color)
  }
  p <- p +  geom_text_repel(
                  aes(y     = ypos,
                      label = paste0(round(prop,2),"%")), 
                      color = "black", 
                      size  = 6) + 
    labs(title  = as.character(title))
  
  return(
    list(
      "table" = df,
      "plot"  = p)
  )
}


GeomSplitViolin <- ggproto("GeomSplitViolin", 
                           GeomViolin, 
                           draw_group = function(self, 
                                                 data, 
                                                 ..., 
                                                 draw_quantiles = NULL) {
                             
                             data <- transform(data, 
                                               xminv = x - violinwidth * (x - xmin), 
                                               xmaxv = x + violinwidth * (xmax - x))
                             grp <- data[1, "group"]
                             newdata <- plyr::arrange(transform(data, 
                                                                x = if (grp %% 2 == 1) xminv else xmaxv), if (grp %% 2 == 1) y else -y)
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

geom_split_violin <- function(mapping = NULL, 
                              data = NULL, 
                              stat = "ydensity", 
                              position = "identity", ..., 
                              draw_quantiles = NULL, 
                              trim = TRUE, 
                              scale = "area", 
                              na.rm = FALSE, 
                              show.legend = NA, 
                              inherit.aes = TRUE) {
  layer(data = data, 
        mapping = mapping, 
        stat = stat, 
        geom = GeomSplitViolin, 
        position = position, 
        show.legend = show.legend, 
        inherit.aes = inherit.aes, 
        params = list(trim = trim, 
                      scale = scale, 
                      draw_quantiles = draw_quantiles, 
                      na.rm = na.rm, ...))
}


boostrap_sampling <- function(nb.sampling     = 100,
                              MARGIN,
                              dataset.to.sample,
                              size,
                              logical.replace = TRUE){
  #TO DO : add a FUN argument in order to manipulate data
  
  if(logical.replace == FALSE){
     if(MARGIN == 1 & length(rownames(metrics_nt)) < 100){
    stop("length of row to resample is inferior to the number of sampling")
  }else if(MARGIN == 2 & length(colnames(metrics_nt)) < 100){
    print("check colname length")
    stop("length of column to resample is inferior to the number of sampling")
  }
  }else if(MARGIN == 1){
    tmp_list <- lapply(1:nb.sampling,function(i){
      cat('boostrap sampling : ',i,"/",nb.sampling,"\n")
      tpm_data <- dataset.to.sample[sample(
        rownames(dataset.to.sample),
        size    = size,
        replace = logical.replace
      ),]
    })
  }else if(MARGIN == 2){
    tmp_list <- lapply(1:nb.sampling,function(i){
      cat('boostrap sampling : ',i,"/",nb.sampling,"\n")
      tpm_data <- dataset.to.sample[,sample(
        colnames(dataset.to.sample),
        size    = size,
        replace = logical.replace
      )]
    })
  }else{
        if(MARGIN == 1){
          print("lapply check")
          tmp_list <- lapply(1:nb.sampling,function(i){
            cat('boostrap sampling : ',i,"/",nb.sampling,"\n")
            tpm_data <- dataset.to.sample[sample(
              rownames(dataset.to.sample),
              size    = size,
              replace = logical.replace
            ),]
          })
        }else if(MARGIN == 2){
          print("lapply check")
          tmp_list <- lapply(1:nb.sampling,function(i){
            cat('boostrap sampling : ',i,"/",nb.sampling,"\n")
            tpm_data <- dataset.to.sample[,sample(
              colnames(dataset.to.sample),
              size    = size,
              replace = logical.replace
            )]
          })
        }else{
          stop("argument MARGIN indicates (1) rows and (2) columns.")
        }
    }
    return(tmp_list)
}

boostrap_comparaison_function <- function(list.to.use,
                                          dataframe.to.use,
                                          col.to.check,
                                          FUN = NULL,
                                          value.to.get = NULL){
  if(!is.null(FUN)){
    
  lapply(1:length(list.to.use),function(i){
    tmp_df       <- list.to.use[[i]]
    
    tmp_function <- do.call(what = FUN,
                            args=list(tmp_df[,col.to.check],
                                      dataframe.to.use[,col.to.check]))
    tmp_result   <- tmp_function[[value.to.get]]

  })
    
    }else{
    
    lapply(1:length(list.to.use),function(i){
      tmp_df     <- list.to.use[[i]]
      })
  }
}
