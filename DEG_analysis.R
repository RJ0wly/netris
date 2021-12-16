base::lapply(c("tidymodels",
               "tidyverse",
               "ggplot2",
               "ggrepel",
               "recipes",
               "DESeq2",
               "here",
               "ggpubr",
               "clusterProfiler",
               "survival",
               "survminer",
               "lubridate",
               "DOSE",
               "enrichplot",
               "enrichR",
               "glmnet",
               "cowplot"),
             require, character.only=T)
source("Utils.R")


TCGA_analysis = lapply(c("CESC", "HNSC","HNSC","LUAD","LUNG","LUSC","OV","UCEC"),TCGA_DEG_analysis)

main_p <- ggarrange(TCGA_analysis[[1]][["plot"]],
                    TCGA_analysis[[2]][["plot"]],
                    TCGA_analysis[[3]][["plot"]],
                    TCGA_analysis[[4]][["plot"]],
                    TCGA_analysis[[5]][["plot"]],
                    TCGA_analysis[[6]][["plot"]],
                    TCGA_analysis[[7]][["plot"]],
                    labels = c("CESC", "HNSC","HNSC","LUAD","LUNG","LUSC","OV","UCEC"),
                    ncol = 2, nrow = 4)
png("TCGA_analysis.png",
    width=4000,
    height=2260)
plot(main_p)
dev.off()

main_p <- arrange_ggsurvplots(TCGA_analysis[[1]][["survplot"]],
                    TCGA_analysis[[2]][["survplot"]],
                    TCGA_analysis[[3]][["survplot"]],
                    TCGA_analysis[[4]][["survplot"]],
                    TCGA_analysis[[5]][["survplot"]],
                    TCGA_analysis[[6]][["survplot"]],
                    TCGA_analysis[[7]][["survplot"]],
                    labels = c("CESC", "HNSC","HNSC","LUAD","LUNG","LUSC","OV","UCEC"),
                    ncol = 2, nrow = 4)
png("TCGA_analysis_surv.png",
    width=4000,
    height=2260)
plot(main_p)
dev.off()





cor_mat <-cor(as.matrix(data2),method="pearson")
cor_mat_CDH1 <- as.data.frame(cor_mat['CDH1',])
cor_mat_FN1 <- as.data.frame(cor_mat['FN1',])
cor_mat_CDH2 <- as.data.frame(cor_mat['CDH2',])
cor_mat_NTN1 <- as.data.frame(cor_mat['NTN1',])

colnames(cor_mat_CDH1) <- c("correlation")
colnames(cor_mat_FN1) <- c("correlation")
colnames(cor_mat_CDH2) <- c("correlation")
colnames(cor_mat_NTN1) <- c("correlation")

cor_mat_CDH1 <- cbind(cor_mat_CDH1,rownames(cor_mat_CDH1))
test <- cor_mat_CDH1[cor_mat_CDH1$correlation > abs(0.5),]
test <- test %>% drop_na()

cor_mat_FN1 <- cbind(cor_mat_FN1,rownames(cor_mat_FN1))
test2 <- cor_mat_FN1[cor_mat_FN1$correlation > abs(0.5),]
test2 <- test2 %>% drop_na()


cor_mat_CDH2 <- cbind(cor_mat_CDH2,rownames(cor_mat_CDH2))
test3 <- cor_mat_CDH2[cor_mat_CDH2$correlation > abs(0.4),]
test3 <- test3 %>% drop_na()


cor_mat_NTN1 <- cbind(cor_mat_NTN1,rownames(cor_mat_NTN1))
test4 <- cor_mat_NTN1[cor_mat_NTN1$correlation > abs(0.4),]
test4 <- test3 %>% drop_na()

intersect(test$`rownames(cor_mat_CDH1)`,test3$`rownames(cor_mat_CDH2)`)

# test --------------------------------------------------------------------
data2 <-readRDS("/home/rj/Desktop/TCGA_project/Project_TCGA_signature/datasets/TCGA_UCEC.RDS")
data2 <- as.data.frame(t(data2))

# data4 <- read.csv("/home/rj/Downloads/TCGA.UCEC.sampleMap_RPPA_RBN.gz",sep="")
# rownames(data4) <- data4[,1]
# data4 <- data4[,-1]
# data4 <- as.data.frame(t(data4))
# 
# inter <-  intersect(rownames(data4),rownames(data))
# 
# data4 <- data4[inter,]
# data <- data[inter,]
# 
# data3 <- cbind(data4,"NTN1"=data[,"NTN1"])
# 
# cor_mat <-cor(as.matrix(data),method="pearson")
# cor_mat_NTN1 <- as.data.frame(cor_mat['NTN1',])
# colnames(cor_mat_NTN1) <- c("correlation")
# cor_mat_NTN1 <- cbind(cor_mat_NTN1,rownames(cor_mat_NTN1))
# test <- cor_mat_NTN1[cor_mat_NTN1$correlation > abs(0.45),]
# test <- test %>% drop_na()
# rownames(test)


file_name_path <- "/home/rj/Desktop/TCGA_project/Project_TCGA_signature/scripts/20530_genes_names.txt"
genes_name_list <- read.table(file_name_path)
file_path <- paste0(paste0(paste0("/home/rj/Desktop/TCGA_project/Project_TCGA_signature/datasets/",
                                  "TCGA_"),
                           "UCEC"),
                    ".RDS")
cat(paste0(as.character("UCEC"),"\n"))
data <- readRDS(file_path)
data <- as.data.frame(t(data))
data <- convert_gene_id_to_ensembl_from_symbol(data)
sample_ID <- rownames(data)
data <- as.data.frame(apply(data,2,as.numeric))
rownames(data) <- sample_ID
col_genes_names <- colnames(data)

score_NTN1 <- read.table("/home/rj/Desktop/TCGA_project/Project_TCGA_signature/scripts/enhancer_promoter_NTN1.txt")
intersect(colnames(data),score_NTN1$V)


test <- sapply(as.data.frame(t(data[,intersect(colnames(data),score_NTN1$V)])), function(df) (sum(df)))

test2 <- sapply(as.data.frame(test), function(df) (df-mean(df))/sd(df))

df_meta_data <- cbind(df_meta_data,test2)

# density plot
p_density <- ggplot2::ggplot(df_meta_data,aes(x=NTN1)) +
  ggplot2::geom_histogram(aes(y=..density..), colour="black", fill="white") +
  ggplot2::geom_density(alpha=.2, fill="#FF6666")

df_meta_data <- df_meta_data %>% mutate(netrin_status = case_when_character_type(
  var_of_interest = test,
  min = -1,
  max = 1))

# binding status with expression matrix
df <- cbind(data,"status"=df_meta_data$netrin_status)

df_tmp <- subset(df, status == "pos" | status == "neg")
status <- df_tmp$status

intersect_genes_list <- intersect(colnames(data),genes_name_list$V1)
# high counts variables
index <- colSums(df_tmp[,intersect_genes_list])
filter <- index > 20
RNA_data_h <- df_tmp[,intersect_genes_list][,filter]

#low counts variables
index <- colSums(df_tmp[,intersect_genes_list])
filter <- index < 20
RNA_data_l <- df_tmp[,intersect_genes_list][,filter]

# drop na

RNA_data_h <- RNA_data_h %>% drop_na()

# most variable column
var.all <- RNA_data_h %>% apply(2,var)%>% sort(decreasing=T)
var.all <- var.all %>% as.data.frame()
rowOfvarAll <- row.names(var.all)
mostvariable <- rowOfvarAll[1:10000]

RNA_data_h_10K <- RNA_data_h[,mostvariable]

rec <- recipes::recipe( ~ ., data = RNA_data_h_10K)

pca_trans <- rec %>%
  recipes::step_pca(all_numeric(),num_comp = 10)

pca_estimates <- recipes::prep(pca_trans, training = RNA_data_h_10K)
pca_data <- bake(pca_estimates, RNA_data_h_10K)

p_pca1 <-ggplot(pca_data,aes(x=PC01,y=PC02,colour=status)) + geom_point(size=3)
p_pca1 <- p_pca1 + theme(legend.position="top") + ggtitle("10K most variables genes")

# DEG
status_df <- as.data.frame(status)
meta_data <- cbind(rownames(RNA_data_h),status_df)

mat <- t(as.matrix(RNA_data_h))
meta_data$status <- factor(meta_data$status)

dds <- DESeqDataSetFromMatrix(countData=round(mat),colData=meta_data,design= ~ status)
dds <- estimateSizeFactors(dds)
dds <- DESeq(dds)
res <- results(dds)
resDF = as.data.frame(res)
resDF = cbind(resDF,"SYMBOL_gene"=rownames(resDF))
keep <- resDF[abs(resDF$log2FoldChange) > 1 & resDF$padj < 0.05 & resDF$baseMean > 1,]
resDF$genelabels <- factor(resDF$SYMBOL_gene, levels = rownames(keep))

# vulcano plot
p_vulcano <- ggplot2::ggplot(data=resDF,
                             aes(x=log2FoldChange,
                                 y=-log10(padj),
                                 label = genelabels )) +
  ggplot2::geom_point(alpha=0.2,
                      colour="dodgerblue2",
                      size=3) +
  ggrepel::geom_text_repel(col = "black",
                           na.rm = TRUE,
                           box.padding = unit(0.45, "lines"),
                           hjust = 1,
                           force = 10,
                           force_pull = 10,
                           max.overlaps = 100) +
  ggplot2::geom_hline(yintercept = 1.3,linetype = "dashed")+
  ggplot2::geom_vline(xintercept = 1,linetype = "dashed") +
  ggplot2::geom_vline(xintercept = -1,linetype = "dashed") +
  theme_bw()


rec <- recipes::recipe( ~ ., data = RNA_data_h[,rownames(keep)])

pca_trans <- rec %>%
  recipes::step_pca(all_numeric(),num_comp = 10)

pca_estimates <- recipes::prep(pca_trans, training = RNA_data_h[,rownames(keep)])
pca_data <- bake(pca_estimates, RNA_data_h[,rownames(keep)])

if(length(rownames(keep)) < 10){
  p_pca2 <- ggplot(pca_data,aes(x=PC1,y=PC2,colour=status)) + geom_point(size=3)
  p_pca2 <- p_pca2 + theme(legend.position="top") + ggtitle("most differentially expressed genes")
}else{
  p_pca2 <- ggplot(pca_data,aes(x=PC01,y=PC02,colour=status)) + geom_point(size=3)
  p_pca2 <- p_pca2 + theme(legend.position="top") + ggtitle("most differentially expressed genes")
}

# Survival analysis part

file_path_os <- paste0(paste0(paste0("/home/rj/Desktop/TCGA_project/Project_TCGA_signature/OS_datasets/",
                                     "OS_"),
                              "UCEC"),
                       ".txt")
OS_data <- read.csv(file_path_os,sep="")
OS_data$sample <- stringr::str_replace_all(OS_data$sample,"-",".")
OS_data$sample <- OS_data$sample %>% substring(1,15)
intersect_sample <- intersect(sample_ID,OS_data$sample)

df_meta_data <- cbind(df_meta_data,"sample"=sample_ID)

df_surv <- merge(df_meta_data,OS_data,by="sample")

cor(df_surv$NTN1,df_surv$OS.time)

fit <- survfit(Surv(OS.time,OS) ~ + netrin_status, data = df_surv)

p_surv <- ggsurvplot(fit,
                     pval = TRUE, conf.int = TRUE,
                     risk.table = TRUE, # Add risk table
                     risk.table.col = "strata", # Change risk table color by groups
                     linetype = "strata", # Change line type by groups
                     surv.median.line = "none", # Specify median survival
                     ggtheme = theme_bw(), # Change ggplot2 theme
                     fontsize = 9,
                     pval.size = 10,
                     censor.size = 5)

df_surv2 <- df_surv[df_surv$netrin_status == "pos" | df_surv$netrin_status == "neg",]

fit2 <- survfit(Surv(OS.time,OS) ~ + netrin_status, data = df_surv2)

p_surv2 <- ggsurvplot(fit2,
                      pval = TRUE, conf.int = TRUE,
                      risk.table = TRUE, # Add risk table
                      risk.table.col = "strata", # Change risk table color by groups
                      linetype = "strata", # Change line type by groups
                      surv.median.line = "none", # Specify median survival
                      ggtheme = theme_bw(), # Change ggplot2 theme
                      fontsize = 9,
                      pval.size = 10,
                      censor.size = 5)

main_p <- ggarrange(plotlist=list(p_pca1,
                    p_pca2,
                    p_density,
                    p_vulcano),
                    ncol = 2, nrow = 2)

splots <-list()
splots[[1]] <- p_surv
splots[[2]] <- p_surv2

survival_plot <- arrange_ggsurvplots(splots,print = FALSE,
                                     ncol = 2, nrow = 1, risk.table.height = 0.4)

obj <- list("plot"=main_p,"result"=keep,"surv_plot"=survival_plot)
ggsave("survplot_.pdf", survival_plot)
