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
               "shinycssloaders"),
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
      selectizeInput("select_method_clustering",
                     choices  = c("Hierarchical Clustering","K-means","Leiden: community clustering"),
                     label    = "Select clustering",
                     multiple = F),
      br(),
      conditionalPanel(
        condition = "input.select_method_clustering == 'Leiden: community clustering'",
      sliderInput("cluster_resolution_leiden", "Clustering resolution", 
                  value = 0.5, 
                  min   = 0, 
                  max   = 1,
                  step  = 0.1)
      ,
      br(),
      actionButton("cluster_resolution_action_button_leiden","Compute cluster")
      ),
      conditionalPanel(
        condition = "input.select_method_clustering == 'K-means'",
        sliderInput("cluster_resolution_k_means", "Number of clusters", 
                    value = 2, 
                    min   = 0, 
                    max   = 10,
                    step  = 1)
        ,
        br(),
        actionButton("cluster_resolution_action_button_k_means","Compute cluster")
      ),
      conditionalPanel(
        condition = "input.select_method_clustering == 'Hierarchical Clustering'",
        sliderInput("cluster_resolution_hierarchical_clustering", "Number of clusters :", 
                    value = 2, 
                    min   = 0, 
                    max   = 10,
                    step  = 1)
        ,
        br(),
        actionButton("cluster_resolution_action_button_hierachical_clustering","Compute cluster")
      )
      ,
      br(),
      selectizeInput("name_gene",
                     choices = NULL,
                     label = "ID genes",
                     selected = NULL,
                     multiple = TRUE)
    ),
    mainPanel(
      textOutput("obj_text"),
      plotOutput("cluster_plot"),
      plotOutput("harmony_plot"),
      plotOutput("feature_plot"),
      plotOutput("violin_plot")
      
      
    )
  )
)

server <- function(input, output, session) {
  
  
  path_data_to_data <- eventReactive( input$load_obj_data, {
    path_to_load_data <- 
      paste0(main_path,
             "/results/multiple_slices_results/",
             as.character(input$name_patient),
             "/R_obj_visium",
             as.character(input$name_patient),
             ".RDS")
    
    
    return(path_to_load_data)
    
  })
  
  
  
  obj_data <- reactive({
    
    readRDS(path_data_to_data())
    
  })
  
  observeEvent(input$load_obj_data,{
    updateSelectizeInput(session, "name_gene", choices = rownames(obj_data()@assays[["SCT"]]@counts), server = TRUE)
  })
  
  
  output$obj_text <- renderPrint({
    obj_data()
  })
  
  clust_obj <- eventReactive(
    input$cluster_resolution_action_button_leiden, {
      data <- obj_data()
      
      data <- FindClusters(object     = data,
                           algorithm  = 4,
                           resolution = input$cluster_resolution_leiden,
                           verbose    = F)
      
    })
  
  clust_obj <- eventReactive(
    input$cluster_resolution_action_button_hierarchical_clustering, {
      data <- obj_data()
      
      dist_mat <- dist(t(data@assays[["SCT"]]@scale.data), 
                       method = 'euclidean')
      hclust_avg <- hclust(dist_mat, 
                           method = 'centroid')
      cut_avg <- cutree(hclust_avg, 
                        k = input$cluster_resolution_hierarchical_clustering)
      data@meta.data$hierarchical_clustering <- cut_avg
      return(data)
    }
  )
  
  clust_obj <- eventReactive(
    input$cluster_resolution_action_button_k_means, {
      data <- obj_data()
      
      k_mean_res <- kmeans(t(data@assays[["SCT"]]@scale.data), 
                           input$cluster_resolution_k_means, 
                           iter.max = 10, 
                           nstart = 1)
      
      data@meta.data$k_means <- k_mean_res[["cluster"]]
      return(data)
    }
  )
  
  output$cluster_plot <- renderPlot({
    
    if(input$select_method_clustering == "Hierarchical clustering"){
      
      Idents(clust_obj()) <- "hierarchical_clustering"
      SpatialDimPlot(object    = clust_obj(),
                    label      = TRUE,
                    label.size = 3)
      
    }else if(input$select_method_clustering == "K-means"){
      Idents(clust_obj()) <- "k_means"
      SpatialDimPlot(object    = clust_obj(),
                     label      = TRUE,
                     label.size = 3)
    }else{
      SpatialDimPlot(object    = clust_obj(),
                     label      = TRUE,
                     label.size = 3)
    }
  })
  
  output$harmony_plot <- renderPlot({
    DimPlot(object    = clust_obj(), 
            reduction = "harmony")
  })
  
  
  output$feature_plot <- renderPlot({
    req(input$name_gene)
    req(clust_obj)
    
    SpatialPlot(object     = clust_obj(),
                features   = input$name_gene,
                alpha      = c(0.05, 1),
                min.cutoff = 0)
  })
  
  output$violin_plot <- renderPlot({
    req(input$name_gene)
    req(clust_obj)
    
    VlnPlot(object    = clust_obj(),
            features  = input$name_gene,
            pt.size   = 0)
  })
}

shinyApp(ui, server)
