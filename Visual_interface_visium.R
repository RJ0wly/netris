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
               "plotly"),
             require, character.only=T)

main_path <- getwd()

# list name of files in folder_visium and store path of this folder into a vector that contains those paths
list_files <- list.files(paste0(main_path,"/folder_visium"))

# extract names of patients and ID of treatment and put in a list
list_patient_filename <- list("ID_patient"   = unique(str_sub(list_files,1,6)),
                              "ID_treatment" = unique(str_sub(list_files,-4,-1)))

# create a list with the path of each visium files per patients
path_list_to_file_patient <- list()
for(i in 1:length(list_patient_filename[["ID_patient"]])){
  for(j in 1:length(list_patient_filename[["ID_treatment"]])){
    tmp_file_name <- paste0(list_patient_filename[["ID_patient"]][i],
                            "_",
                            list_patient_filename[["ID_treatment"]][j])
    
    
    path_list_to_file_patient[[as.character(list_patient_filename[["ID_patient"]][i])]][j] = paste0(main_path,
                                                                                                    "/folder_visium/",
                                                                                                    tmp_file_name,
                                                                                                    "/outs")
  }
}

names(path_list_to_file_patient)

ui <- fluidPage(
  sidebarLayout(
    sidebarPanel(
      selectizeInput("name_patient",
                     choices = names(path_list_to_file_patient),
                     label = "Select patient",
                     multiple = F),
      br(),
      actionButton("load_obj_data","Load data"),
      br(),
      sliderInput("cluster_resolution", "Clustering resolution", value = 0.5, min = 0, max = 1),
      br(),
      actionButton("cluster_resolution_action_button","Compute cluster")
    ),
    mainPanel(
      plotOutput("cluster_plot")
    )
  )
)

server <- function(input, output, session) {
  
  obj_data <- eventReactive( input$load_obj_data, {
    path_to_load_data <- 
      paste0(main_path,
             "/results/multiple_slices_results/",
             as.character(input$name_patient),
             "/R_obj_visium",
             as.character(input$name_patient),
             ".RDS")
    
    obj <- readRDS(path_to_load_data)
    
    return(obj)
  })
  
  obj_data <- eventReactive( input$cluster_resolution_action_button, {
    req(obj_data)
    req(input$cluster_resolution)
    
    obj_data()  <- FindClusters(obj_data(),
                                verbose    = FALSE,
                                algorithm  = 4,
                                resolution = input$cluster_resolution)
  })
  
  output$cluster_plot <- renderPlot({
              req(obj_data)
              req(input$cluster_resolution)
              
              SpatialDimPlot(obj_data(),
                         label      = TRUE,
                         label.size = 3,
                         group.by = c(colnames(obj_data()@meta.data)[length(obj_data()@meta.data)]))
  })
}

shinyApp(ui, server)

