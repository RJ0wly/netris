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

change_column_name <- function(data = data,
                                old.name,
                                new.name){
    if(!is.character(old.name) | !is.character(new.name) ){
      stop("old.name or new.name is not a character type")
    }else{colnames(data)[which(names(data) == old.name)] <- new.name
  }
}

logistic_regression_bionomial <- function(training.data       = training.data,
                                          testing.data        = testing.data,
                                          var.name.to.predict = var.name.to.predict,
                                          classes.to.predict  = classes.to.predict,
                                          testing.response    = testing.response){
  
  dep <- deparse(substitute(var.name.to.predict))
  dep <- paste0("factor(", dep, ")")
  fmla <- paste0(dep, "~.")
  fmla <- as.formula(fmla)
  
  model              <- glm(formula = fmla, 
                             data   = training.data, 
                             family = "binomial")
  list_summary_model <- summary(model)
  probabilities      <- model %>% predict(testing.data, 
                                          type   = "response")
  
  predicted_classes  <- ifelse(probabilities> 0.5, classes.to.predict[1], classes.to.predict[2])
  accuracy           <- mean(predicted_classes == testing.response)
  
  plot_roc = roc(testing.response~ probabilities,
                 plot      = TRUE, 
                 print.auc = TRUE)
  
  cat("Model accuracy is",accuracy,"% : \n")
  print(list_summary_model[["coefficients"]])
  
  
  
  res <- list("model summary"  = list_summary_model,
              "predict classes" = predicted_classes,
              "predict obj"     = probabilities,
              "accuracy"        = accuracy,
              "roc curve"       = plot_roc)
}


counts_and_create_dataframe_from_a_list<- function(list.to.use){
  unlist_list <- unlist(list.to.use)
  counts_and_make_df <- as.data.frame(table(unlist_list))
  #add percent
  counts_and_make_df$percent <- counts_and_make_df$Freq*100/max(counts_and_make_df$Freq)
  return(counts_and_make_df)
}

bootstrap_lasso=function(x, 
                         y,
                         bs          = 10, 
                         kfold       = 10,
                         family      = "binomial",
                         standardize = FALSE,
                         maxit       = 10^5){
  rowx <- nrow(x)
  n    <- length(y)
  if (rowx != n){
    stop("The number of rows in x is not equal to the length of y!")
  }
  res <- lapply(1:bs,function(i){
    start_time= as.numeric(Sys.time())
    repeat{
      s <- sample(n, replace=TRUE)
      # repeat while it does not have at least two discrete value from each group
      if(length(table(y[s])) >= 2 & length(table(y[-s])) >= 2)
        break
    }
    # selecting row
    BoostrapX           <- as.matrix(x[s, ])
    colnames(BoostrapX) <- colnames(x)
    BoostrapY           <- y[s]
    
    # logistic regression lasso
    if(family == "binomial"){
      
      fit       <- glmnet(x      = BoostrapX, 
                          y      = BoostrapY,
                          family =  "binomial",
                          standardize = FALSE)
      
      cvfit     <- cv.glmnet(x            = BoostrapX,
                             y                = BoostrapY,
                             type.measure     = "auc",
                             nfolds           = kfold,
                             family           = "binomial",
                             alpha            = 1,
                             maxit            = maxit,
                             standardize      = standardize)
      
      tmp_coeffs <- coef(object  = fit, 
                         s       = cvfit$lambda.min)
      
      best_lam   <- cvfit$lambda.min
      
      variable_selected <- tmp_coeffs@Dimnames[[1]][tmp_coeffs@i + 1]
      variable_selected <- variable_selected[-1]
    }
    
    end_time     <- as.numeric(Sys.time())
    cat("Boostrap ", i,"/",bs, "time taken",end_time - start_time,"s\n")
    
    model_final <- list("fit" = fit,
                        "cvfit" = cvfit,
                        "selected variable" = variable_selected)
  })
  return(res)
}

reduction_plot <- function(object, 
                           method.of.reduction = c("pca","tsne","umap"),
                           list.to.color,
                           verbose = TRUE){
  res <- list()
  to_plot <- function(data){
    lapply(1:length(list.to.color),function(i){
      to_color <- list.to.color[[i]]
      ggplot(data, aes(x=V1,y=V2,color = to_color)) + 
        geom_point(size=1,alpha=0.7) +
        theme_light() +
        xlab("PC1") +
        ylab("PC2") +
        scale_color_discrete(paste0(names(list.to.color[i])))})
  }
  
  if("umap" %in% method.of.reduction){
    if(verbose){cat("umap running ... \n")}
    umap_data      <- uwot::umap(X = object[,sapply(object, class) == "numeric"],
                                 verbose = TRUE)
    umap_data      <- as.data.frame(umap_data)
    umap_list_plot <- to_plot(umap_data)
    res <- c(res, "umap" = list(umap_list_plot))

  }
  if("tsne" %in% method.of.reduction){
    if(verbose){cat("tsne running ... \n")}
    rtsne_data     <- Rtsne::Rtsne(X = as.matrix(object[,sapply(object, class) == "numeric"]),
          pca_center = FALSE,
          normalize = FALSE,
          check_duplicates = FALSE,
          verbose = TRUE)
    rtsne_data     <- as.data.frame(rtsne_data[["Y"]])
    tsne_list_plot <- to_plot(rtsne_data)
    res <- c(res, "tsne" = list(tsne_list_plot))
  }
  if("pca" %in% method.of.reduction){
    if(verbose){cat("pca running ... \n")}
    pca_rec       <- recipe(~.,
                        data = object[,sapply(object, class) == "numeric"]) %>%
      step_pca(all_predictors(),num_comp = 2)

    pca_prep      <- prep(pca_rec)
    pca_data      <- bake(pca_prep, object)
    colnames(pca_data) <- c("V1","V2")
    pca_list_plot <- to_plot(pca_data)
    res <- c(res, "pca" = list(pca_list_plot))
  }
  return(res)
}
