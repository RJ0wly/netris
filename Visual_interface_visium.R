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
      sliderInput("cluster_resolution", "Clustering resolution", 
                  value = 0.5, 
                  min   = 0, 
                  max   = 1,
                  step  = 0.1),
      br(),
      actionButton("cluster_resolution_action_button","Compute cluster"),
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
    input$cluster_resolution_action_button, {
      data <- obj_data()
      
      data <- FindClusters(object     = data,
                           algorithm  = 4,
                           resolution = input$cluster_resolution,
                           verbose    = F)
      
    })
  
  output$cluster_plot <- renderPlot({
    SpatialDimPlot(object    = clust_obj(),
                  label      = TRUE,
                  label.size = 3)
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
