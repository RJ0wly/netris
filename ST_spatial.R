base::lapply(c("Seurat",
               "tidyverse",
               "SeuratObject",
               "SeuratData",
               "patchwork",
               "tidyverse",
               "singscore",
               "escape",
               "dittoSeq",
               "ggpubr",
               "ggrepel",
               "cowplot",
               "RColorBrewer",
               "ggrepel",
               "harmony",
               "rstatix",
               "STutility"),
             require, character.only=T)
main_path <- getwd()

samples=c("/home/rj/Desktop/analysis/visium/visium/folder_visium/01_034_C1d1/outs/filtered_feature_bc_matrix.h5",
          "/home/rj/Desktop/analysis/visium/visium/folder_visium/01_034_C3d1/outs/filtered_feature_bc_matrix.h5")
spotfiles=c("/home/rj/Desktop/analysis/visium/visium/folder_visium/01_034_C1d1/outs/spatial/tissue_positions_list.csv",
            "/home/rj/Desktop/analysis/visium/visium/folder_visium/01_034_C3d1/outs/spatial/tissue_positions_list.csv")
imgs=c("/home/rj/Desktop/analysis/visium/visium/folder_visium/01_034_C1d1/outs/spatial/tissue_hires_image.png",
       "/home/rj/Desktop/analysis/visium/visium/folder_visium/01_034_C3d1/outs/spatial/tissue_hires_image.png")
json=c("/home/rj/Desktop/analysis/visium/visium/folder_visium/01_034_C1d1/outs/spatial/scalefactors_json.json",
       "/home/rj/Desktop/analysis/visium/visium/folder_visium/01_034_C3d1/outs/spatial/scalefactors_json.json")
sample_ID=c("before treatment","after treatment")

infoTable <- data.frame(samples,spotfiles,imgs,json,sample_ID)

se <- InputFromTable(infotable = infoTable, 
                     min.gene.count = 100, 
                     min.gene.spots = 5,
                     min.spot.count = 500,
                     platform =  "Visium")
se <- LoadImages(se, time.resolve = FALSE)

se <- SCTransform(se)

ST.FeaturePlot(se, features = c("EPCAM"), 
               cols = c("dodgerblue2", "seagreen3", "khaki2", "darkgoldenrod2", "brown4"), ncol = 2, pt.size = 2)

ImagePlot(se, method = "raster", annotate = FALSE)

FeatureOverlay(se, features = "EPCAM",
               sampleids = 1:2,
               cols = c("dodgerblue2", "seagreen3", "khaki2", "darkgoldenrod2", "brown4"),
               ncol = 2)

