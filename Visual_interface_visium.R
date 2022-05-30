base::lapply(c("Seurat",
               "tidyverse",
               "SeuratObject",
               "SeuratData",
               "patchwork",
               "tidyverse",
               "ggpubr",
               "ggrepel",
               "cowplot",
               "RColorBrewer",
               "ggrepel",
               "harmony",
               "rstatix",
               "ggplot2",
               "shiny",
               "plotly",
               "shinyFeedback",
               "shinycssloaders",
               "STutility"),
             require, character.only=T)

main_path <- getwd()

# list name of files in folder_visium and store path of this folder into a vector that contains those paths
list_files <- list.files(paste0(main_path,"/folder_visium"))

# extract names of patients and ID of treatment and put in a list
patient_ID <- unique(str_sub(list_files,1,6))

ui <- fluidPage(
  sidebarLayout(
    sidebarPanel(
      selectizeInput("name_patient",
                     choices  = patient_ID,
                     label    = "Select patient",
                     multiple = F),
      br(),
      actionButton("load_obj_data","Load data"),
      br(),
      br(),
      selectizeInput("select_method_clustering",
                     choices  = c("Kmeans","Leiden"),
                     label    = "Select clustering method :",
                     multiple = F),
      br(),
      conditionalPanel(
        condition = "input.select_method_clustering == 'Leiden'",
        sliderInput("cluster_resolution_leiden", "Clustering resolution", 
                    value = 0.5, 
                    min   = 0, 
                    max   = 1,
                    step  = 0.1),
        br()
        ),
      conditionalPanel(
        condition = "input.select_method_clustering == 'Kmeans'",
        sliderInput("cluster_resolution_k_means", "Number of clusters :", 
                    value = 2, 
                    min   = 0, 
                    max   = 10,
                    step  = 1),
        br(),
        ),
      actionButton("cluster_action_button","Compute cluster"),
      br(),
      selectizeInput("name_gene",
                     choices = NULL,
                     label = "ID genes",
                     selected = NULL,
                     multiple = TRUE),
      br(),
      selectizeInput("violin_groupby",
                     choices = NULL,
                     label = "Select group of identification : ",
                     selected = NULL,
                     multiple = FALSE)
    ),
    mainPanel(
      textOutput("obj_text"),
      shinycssloaders::withSpinner(
        ui_element = plotOutput(outputId = "cluster_plot"),
        type       = 8,
        size       = 2,
        color      = "grey"),
      shinycssloaders::withSpinner(
        ui_element = plotOutput(outputId = "harmony_plot"),
        type       = 8,
        size       = 2,
        color      = "grey"),
      plotOutput("feature_plot"),
      shinycssloaders::withSpinner(
        ui_element = plotOutput(outputId = "violin_plot"),
        type       = 8,
        size       = 2,
        color      = "grey")
      
    )
  )
)

server <- function(input, output, session) {
  
  
  path_data_to_data <- eventReactive(input$load_obj_data, {
    path_to_load_data <- 
      paste0(main_path,
             "/results/multiple_slices_results/",
             as.character(input$name_patient),
             "/R_obj_visium",
             as.character(input$name_patient),
             ".RDS")
    
    
    return(path_to_load_data)
    
  })
  
  path_to_data_ST_utility <- eventReactive(input$load_obj_data,{
    path_to_load_data <- 
      paste0(main_path,
             "/results/multiple_slices_results/",
             as.character(input$name_patient),
             "/STutility_obj",
             as.character(input$name_patient),
             ".RDS")
  })
  
  
  
  obj_data <- reactive({
    readRDS(path_data_to_data())
  })
  
  obj_STutility <- reactive({
    readRDS(path_to_data_ST_utility())
  })
  
  observeEvent(input$load_obj_data,{
    updateSelectizeInput(session, "name_gene", choices = rownames(obj_data()@assays[["SCT"]]@counts), server = TRUE)
  })
  
  
  output$obj_text <- renderPrint({
    obj_data()
  })
  
  leiden_function <- function(object,res){
      data <- object
      print("Leiden occuring")
      data <- FindClusters(object     = data,
                           algorithm  = 4,
                           resolution = res,
                           verbose    = F)
      return(data)
  }
  
  k_means_function <- function(object,k){
    data <- object
    k_mean_res <- kmeans(t(data@assays[["SCT"]]@scale.data), 
                         k, 
                         iter.max = 100, 
                         nstart = 1)
    
    data@meta.data$k_means <- as.factor(k_mean_res[["cluster"]])
    return(data)
  }
  
  clust_obj <- eventReactive(input$cluster_action_button,{
    switch(input$select_method_clustering,
           Leiden = leiden_function(object = obj_data(), res = input$cluster_resolution_leiden),
           Kmeans = k_means_function(object = obj_data(), k = input$cluster_resolution_k_means)
    )
  })
  

  cluster_plot_kmeans <- function(object){
    Idents(object) <- "k_means"
    SpatialPlot(object        = object,
                   label      = TRUE,
                   label.size = 3)
  }
  
  cluster_plot_leiden <- function(object,res){
    Idents(object) <- paste0("SCT_snn_res.",as.character(res))
    SpatialDimPlot(object     = object,
                   label      = TRUE,
                   label.size = 3)
  }
  
  plot_cluster_type <- eventReactive(input$cluster_action_button,{
    req(clust_obj())
    switch(input$select_method_clustering,
           Leiden = cluster_plot_leiden(clust_obj(),input$cluster_resolution_leiden),
           Kmeans = cluster_plot_kmeans(clust_obj())
           )
  })
  
  output$cluster_plot <- renderPlot({plot_cluster_type()})
  
  plot_harmony_umap <- eventReactive(input$cluster_action_button,{
    req(clust_obj())
    switch(input$select_method_clustering,
           Leiden = DimPlot(object    = clust_obj(), 
                            reduction = "harmony",
                            group.by = paste0("SCT_snn_res.",
                                              as.character(input$cluster_resolution_leiden)))
           ,
           Kmeans = DimPlot(object    = clust_obj(), 
                            reduction = "harmony",
                            group.by = "k_means")
    )
  })
  
  output$harmony_plot <- renderPlot({plot_harmony_umap()})
  
  
  output$feature_plot <- renderPlot({
    req(input$name_gene)
    req(obj_STutility)
    
    ST.FeaturePlot(obj_STutility(), 
                   features = c(input$name_gene),
                   cols = c("dodgerblue2", "seagreen3", "khaki2", "darkgoldenrod2", "brown4"),
                   ncol = 2,
                   pt.size = 2)
  })

  
  violin_leiden <- function(object,res,gene_input,group_id){
    Idents(object) <- paste0("SCT_snn_res.",
                             as.character(res))
    VlnPlot(object    = object,
            group.by = group_id,
            features  = gene_input,
            pt.size   = 0)
    
  }
  
  violin_k_means <- function(object,gene_input,group_id){
    Idents(object) <- "k_means"
    VlnPlot(object    = object,
            group.by = group_id,
            features  = gene_input,
            pt.size   = 0)
  }
  
  observeEvent(input$cluster_action_button,{
    updateSelectizeInput(session, "violin_groupby", 
                         choices = names(clust_obj()@meta.data)[sapply(clust_obj()@meta.data, is.factor)], 
                         server = TRUE)
  })
  
  plot_violin_cluster_type <- reactive({
    req(input$name_gene)
    req(clust_obj())
    req(input$violin_groupby)
    
    switch(input$select_method_clustering,
           Leiden = violin_leiden(object     = clust_obj(),
                                  res        = input$cluster_resolution_leiden,
                                  gene_input = input$name_gene,
                                  group_id = input$violin_groupby),
           Kmeans = violin_k_means(object = clust_obj(),
                                   gene_input = input$name_gene,
                                   group_id = input$violin_groupby)
    )
  })
  
  output$violin_plot <- renderPlot({
    req(input$name_gene)
    req(clust_obj)
    req(input$violin_groupby)
    
    plot_violin_cluster_type()
  })
}

shinyApp(ui, server)
