library(tidyverse)
library(clusterProfiler)
library(biomaRt)
library(embed)
library(readxl)
library(ggplot2)
library(ggpubr)
library(rstatix)
library(survival)
library(survminer)
library(lubridate)



df <- readRDS("mat_exp_data.RDS")

meta_data <- read_excel("/home/rj/Desktop/analysis/synovial_sarcoma_valacker/SS -Rémy.xlsx")

meta_data <- meta_data[meta_data$n..ARN %in% colnames(df),]

df <- as.data.frame(t(df))

meta_data <- meta_data[meta_data$n..ARN %in% rownames(df),]
meta_data <- meta_data %>% mutate(Sub.type.sample =  case_when(`Sous.type.d.√©chantillon` == "M√©tastase" ~ 'Metastasis',
                                                               `Sous.type.d.√©chantillon` == "Meta" ~ 'Metastasis',
                                                               `Sous.type.d.√©chantillon` == "primitive tumor" ~ 'Primary',
                                                               `Sous.type.d.√©chantillon` == "Tumeur primitive" ~ 'Primary',
                                                               `Sous.type.d.√©chantillon` == "tumeur primitive" ~ 'Primary',
                                                               `Sous.type.d.√©chantillon` == "tumeur primitvie" ~ 'Primary',
                                                               `Sous.type.d.√©chantillon` == "Local reccurrence" ~ 'Local reccurrence'))

meta_data <- meta_data %>% mutate(Histotype =  case_when(Histotype == "NA" ~ 'Synovial sarcoma',
                                                         Histotype == "Synovial sarcoma" ~ 'Synovial sarcoma',
                                                         Histotype == "Synovial sarcoma - biphasic" ~ 'Synovial sarcoma - biphasic',
                                                         Histotype == "Synovial sarcoma - monophasic" ~ 'Synovial sarcoma - monophasic',
                                                         Histotype == "Synovial sarcoma - NOS" ~ 'Synovial sarcoma - NOS',
                                                         Histotype == "Synovial sarcoma - poorly differentiated" ~ 'Synovial sarcoma - poorly differentiated'))

meta_data <- meta_data %>% mutate(Classification =  case_when(Morphological.classification...62 == "." |
                                                                Morphological.classification...62 ==  "NA" | 
                                                                Morphological.classification...62 == "Other" ~"Other",
                                                              Morphological.classification...62 == "Pleomorphic" ~ "Pleomorphic",
                                                              Morphological.classification...62 == "Pseudo-epithelial" ~ "Pseudo-epithelial",
                                                              Morphological.classification...62 ==  "Round cells" ~ "Round cells",
                                                              Morphological.classification...62 == "Spindle" ~ "Spindle"))


meta_data <- meta_data %>% mutate(Grade.of.tumour =  case_when(Grade.of.tumour == "." |
                                                                 Grade.of.tumour ==  "NA" ~ "NA",
                                                               Grade.of.tumour == "2" ~ "2",
                                                               Grade.of.tumour == "3" ~ "3"))


names(meta_data)[names(meta_data) == "n..ARN"] <- "ID_patient"

# NTN1 --------------------------------------------------------------------

meta_data_NTN1 <- cbind('NTN1'=df[,'NTN1'],"ID_patient"=rownames(df))
meta_data_NTN1 <- meta_data_NTN1 %>% as.data.frame()
meta_data_NTN1$NTN1 <- as.numeric(meta_data_NTN1$NTN1)

meta_data_NTN1 <- meta_data_NTN1 %>% mutate("NTN1_status" = case_when(NTN1 >= quantile(meta_data_NTN1$NTN1,0.50) ~ "High",
                                                                     TRUE ~ "Low"
                                                                     
))

top_gene_deg <- read.csv("/home/rj/Desktop/analysis/visium/visium/visium_deg_low_vs_high.csv")
top_gene_deg <- intersect(top_gene_deg$gene,colnames(df))
df_features <- df[,top_gene_deg]
df_features$ID_patient <- rownames(df_features)

merge_data <- merge(meta_data_NTN1,df_features,by=c("ID_patient","NTN1"))

merge_data <- merge(merge_data,meta_data, by=c("ID_patient"))

sub_merge_data <- merge_data %>% filter(Sub.type.sample == "Primary")


# sub_merge_data <- sub_merge_data %>% mutate(KRT7_status=case_when(KRT7 > 2.5 ~ "Pos",
#                                                                   TRUE ~ "Neg"))
# 
# # KRT7 ~ 2.5 = biphasic
# 
# names(sub_merge_data)[names(sub_merge_data) == "OS Year"] <- "OS"
# names(sub_merge_data)[names(sub_merge_data) == "Vital status"] <- "status"
# 
sub_merge_data <- merge_data %>% filter(Sub.type.sample == "Primary") %>% filter(Histotype == "Synovial sarcoma - biphasic" | Histotype == "Synovial sarcoma - monophasic")
# 
# 
# survdiff(Surv(OS, status) ~ KRT7_status, data = sub_merge_data)
# 
# ggsurvplot(
#   fit = survfit(Surv(OS, status) ~ Histotype, data = sub_merge_data), 
#   xlab = "Year", 
#   ylab = "Overall survival probability")
# 
# sd <- survdiff(Surv(OS, status) ~ Histotype, data = sub_merge_data)
# 1 - pchisq(sd$chisq, length(sd$n) - 1)
# 
# hist(sub_merge_data$OS,na.rm = T)

# ggplot(sub_merge_data,aes(x = Site.of.tumour.category, y = NTN1_status)) + geom_bar(stat="identity",aes(fill=NTN1_status)) + theme_light()
# 
test <- lapply(top_gene_deg, function(i){
  stat_test <- sub_merge_data %>%
    wilcox_test(as.formula(paste(as.character(i), "~", "NTN1_status"))) %>%
    add_significance()
})

names(test) <- top_gene_deg

test <- do.call(rbind.data.frame, test)

test <- test[test$p < 0.05,]
names_genes <- test$.y.

test_2 <- lapply(names_genes, function(i){
  stat_test <- sub_merge_data %>%
    wilcox_test(as.formula(paste(as.character(i), "~", "NTN1_status"))) %>% add_significance()
  stat_test <- stat_test %>% add_xy_position(x = "NTN1_status")
  stat_test$y.position <- stat_test$y.position

  bxp <- ggboxplot(data = sub_merge_data,
                   x    = "NTN1_status",
                   y    = i,
                   fill = "NTN1_status")

  p <- bxp +
    stat_pvalue_manual(stat_test, label = "p") +
    scale_y_continuous(expand = expansion(mult = c(0.05, 0.1))) +
    ylab(paste(i," expression")) +
    xlab("") +
    scale_fill_discrete("NTN1 status") +
    ggtitle(paste(i))
})

pdf("plot_bxp_deg.pdf")
test_2
dev.off()

library(DESeq2)
library(tidyverse)
library(survminer)
library(survival)
library(lubridate)

names(sub_merge_data)[names(sub_merge_data) == "OS Year"] <- "OS"
names(sub_merge_data)[names(sub_merge_data) == "Vital status"] <- "status"

test_survdiff <- lapply(names_genes,function(i){
  res_cox_model <- coxph(as.formula(paste('Surv(OS, status) ~ strata(NTN1_status) +', i)), data = sub_merge_data)
  fit <- survfit(res_cox_model)
  plot_surv <- ggsurvplot(fit, conf.int = F, 
                          censor = FALSE ,
                          surv.median.line = "hv",
                          data = sub_merge_data) + ggtitle(i)
  list("plot"=plot_surv,"summary"=summary(res_cox_model))
  
})


test_survdiff <- lapply(c("MUC1","PCDH18"),function(i){
  res_cox_model <- coxph(as.formula(paste('Surv(OS, status) ~ strata(NTN1_status) +', i)), data = sub_merge_data)
  p_val <- round(summary(res_cox_model)[["coefficients"]][[5]],3)
  fit <- survfit(res_cox_model)
  plot_surv <- ggsurvplot(fit, conf.int = F, 
                          censor = FALSE ,
                          surv.median.line = "hv",
                          data = sub_merge_data,
                          risk.table = TRUE)
  
  plot_surv$table <- plot_surv$table +
    theme(plot.caption = element_text(hjust = 0)) +
    labs(caption = paste0("pvalue = ",p_val, " gene = ", i))
  return(plot_surv)
})



pdf("plot_surv.pdf")
test_survdiff
dev.off()
