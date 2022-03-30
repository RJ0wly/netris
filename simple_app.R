library(shiny)
library(ggplot2)
library(Seurat)
library(SeuratObject)
library(RColorBrewer)
library(shinyFeedback)
library(shinysky)
library(org.Hs.eg.db)


human_gene_vector <- as.vector(unlist(as.list(org.Hs.egSYMBOL)[mappedkeys(org.Hs.egGO)]))
obj <- readRDS("/home/rj/Desktop/Shiny_single_cell_project/datasets/big_obj_seurat.RDS")
ui <- fluidPage(
  sidebarLayout(
  sidebarPanel(
  shinyFeedback::useShinyFeedback(),
  selectInput("dataset",label = "Dataset", choices = names(obj)),
  br(),
  selectizeInput("name_gene",
                 choices = NULL,
                 label = "ID genes",
                 selected = NULL,
                 multiple = TRUE)
  ),
  mainPanel(
    tabsetPanel(
      tabPanel("Clusters", plotOutput("plot_cluster")),
      tabPanel("Plot", plotOutput("plot")),
      tabPanel("Violon plot", plotOutput("violin_plot")),
      tabPanel("Heatmap",plotOutput("Heatmap_plot"))
    )
  )
  )
)

server <- function(input, output, session) {
  data <- reactive({
    req(input$dataset)

    tmp <- obj[[as.character(input$dataset)]]
  })
  
  updateSelectizeInput(session, "name_gene", choices = human_gene_vector, server = TRUE)
  
  vector_gene <- reactive({
    
    req(input$name_gene)
    exists <- input$name_gene %in% rownames(data()@assays[["RNA"]]@data)
    shinyFeedback::feedbackDanger("name_gene", !exists, "Unknown gene(s)")
    req(exists, cancelOutput = TRUE)
    
    as.vector(input$name_gene)
  })
  
  output$text <- renderPrint(summary(data()))
  
  output$plot <- renderPlot(
    FeaturePlot(object  = data(),
                features = c(as.character(vector_gene())),
                min.cutoff = 0.0,
                pt.size = 1,
                cols = rev(brewer.pal(n = 11, name = "RdYlBu")),
                keep.scale = "feature") + theme(legend.position = "None")
  )
  
  output$violin_plot <- renderPlot({
    VlnPlot(object = data(),
            assay = "integrated",
            slot = "data",
            features = c(as.character(vector_gene())),
            pt.size = 0,
            group.by = "seurat_clusters")
  })
  
  output$plot_cluster <- renderPlot({
    DimPlot(object  = data(), 
            reduction = "umap",
            pt.size = 1,
            group.by = "seurat_clusters") + ggtitle("Clusters")
  })

}


shinyApp(ui, server)
