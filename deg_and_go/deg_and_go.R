##### 3.4. Enrichment

library(enrichR)


### 1.1 Calculate DEG for differnet clusters , in o2 code: 3.4_DEG #### ---------------------
#!/usr/bin/Rscript

### 1. ## DEGs for 11 clusters

library(SingleCellExperiment)
library(Seurat)
library(tidyverse)
library(Matrix)
library(scales)
library(cowplot)
library(RCurl)
library(cowplot)
library(BiocManager)
library(multtest)
library(metap)
library(RColorBrewer)
library(Seurat)
library(SeuratDisk)
library(ggplot2)
library(patchwork)

options(future.globals.maxSize = 4000 * 1024^10)
skin.query <- qs_read("01_data/skin.query_use_v2.qs2")
#load("01_data/skin.ref.RData")


tmp <- Idents(object = skin.query)
tmp2 <- gsub("-.*","",as.vector(tmp))
Idents(object = skin.query) <- tmp2

DEG_11cluster <- FindAllMarkers(object = skin.query,
                                logfc.threshold = 0.25)  
write.table(DEG_11cluster , file= "03_results/markers_11cluster_16sam.txt", quote = F, sep = "\t")



#### ---------------------

### 1.2 DEG files

deg_allcelltype <- read.table(file="/n/groups/church/Lily/Skin_aging/scRNA_china/03_analysis_add12/03_results/markers_11cluster_16sam.txt", header = T)
deg_allcelltype_filter <- subset(deg_allcelltype, deg_allcelltype$p_val_adj < 0.05 & deg_allcelltype$avg_log2FC > 0.25 )

write.table(deg_allcelltype_filter, file= "03_results/deg_allcelltype_filter_16sam.txt", quote = F, sep = "\t")


deg_group1 <- subset(deg_allcelltype_filter, cluster=="Granular_cell")
deg_group2 <- subset(deg_allcelltype_filter, cluster=="Spinous_cell")
deg_group3 <- subset(deg_allcelltype_filter, cluster=="Basal_cell")
deg_group4 <- subset(deg_allcelltype_filter, cluster=="Endothelial_cell")
deg_group5 <- subset(deg_allcelltype_filter, cluster=="Fibroblast")
deg_group6 <- subset(deg_allcelltype_filter, cluster=="Hair_follicle")
deg_group7 <- subset(deg_allcelltype_filter, cluster=="Hair_follicle_(Sebocyte)")
deg_group8 <- subset(deg_allcelltype_filter, cluster=="Melanocyte")
deg_group9 <- subset(deg_allcelltype_filter, cluster=="Pericyte")
deg_group10 <- subset(deg_allcelltype_filter, cluster=="Immune_cell")
deg_group11 <- subset(deg_allcelltype_filter, cluster=="Langerhans")




### 1.2 Enrichr

library(enrichR)
setEnrichrSite("Enrichr")
dbs <- listEnrichrDbs()
dbs <- c("GO_Molecular_Function_2025", "GO_Cellular_Component_2025", "GO_Biological_Process_2025")
websiteLive <- TRUE


if (websiteLive) {
  enriched <- enrichr(deg_group1$gene, dbs)
}
if (websiteLive) enriched[["GO_Biological_Process_2015"]]
if (websiteLive) plotEnrich(enriched[[3]], showTerms = 20, numChar = 40, y = "Count", orderBy = "Adjusted.P.value")


enriched_group1 <- enrichr(deg_group1$gene, dbs)
save(enriched_group1, file = "enriched_group1.RData")
plotEnrich(enriched_group1[[3]], showTerms = 20, numChar = 40, y = "Count", orderBy = "Adjusted.P.value")
ggsave("enriched_group1.pdf",device = "pdf",width = 5, height = 8)

enriched_group2 <- enrichr(deg_group2$gene, dbs)
save(enriched_group2, file = "enriched_group2.RData")
plotEnrich(enriched_group2[[3]], showTerms = 20, numChar = 40, y = "Count", orderBy = "Adjusted.P.value")
ggsave("enriched_group2.pdf",device = "pdf",width = 5, height = 8)

enriched_group3 <- enrichr(deg_group3$gene, dbs)
save(enriched_group3, file = "enriched_group3.RData")
plotEnrich(enriched_group3[[3]], showTerms = 20, numChar = 40, y = "Count", orderBy = "Adjusted.P.value")
ggsave("enriched_group3.pdf",device = "pdf",width = 5, height = 8)

enriched_group4 <- enrichr(deg_group4$gene, dbs)
save(enriched_group4, file = "enriched_group4.RData")
plotEnrich(enriched_group4[[3]], showTerms = 20, numChar = 40, y = "Count", orderBy = "Adjusted.P.value")
ggsave("enriched_group4.pdf",device = "pdf",width = 5, height = 8)

enriched_group5 <- enrichr(deg_group5$gene, dbs)
save(enriched_group5, file = "enriched_group5.RData")
plotEnrich(enriched_group5[[3]], showTerms = 20, numChar = 40, y = "Count", orderBy = "Adjusted.P.value")
ggsave("enriched_group5.pdf",device = "pdf",width = 5, height = 8)

enriched_group6 <- enrichr(deg_group6$gene, dbs)
save(enriched_group6, file = "enriched_group6.RData")
plotEnrich(enriched_group6[[3]], showTerms = 20, numChar = 40, y = "Count", orderBy = "Adjusted.P.value")
ggsave("enriched_group6.pdf",device = "pdf",width = 5, height = 8)

enriched_group7 <- enrichr(deg_group7$gene, dbs)
save(enriched_group7, file = "enriched_group7.RData")
plotEnrich(enriched_group7[[3]], showTerms = 20, numChar = 40, y = "Count", orderBy = "Adjusted.P.value")
ggsave("enriched_group7.pdf",device = "pdf",width = 5, height = 8)

enriched_group8 <- enrichr(deg_group8$gene, dbs)
save(enriched_group8, file = "enriched_group8.RData")
plotEnrich(enriched_group8[[3]], showTerms = 20, numChar = 40, y = "Count", orderBy = "Adjusted.P.value")
ggsave("enriched_group8.pdf",device = "pdf",width = 5, height = 8)

enriched_group9 <- enrichr(deg_group9$gene, dbs)
save(enriched_group9, file = "enriched_group9.RData")
plotEnrich(enriched_group9[[3]], showTerms = 20, numChar = 40, y = "Count", orderBy = "Adjusted.P.value")
ggsave("enriched_group9.pdf",device = "pdf",width = 5, height = 8)

enriched_group10 <- enrichr(deg_group10$gene, dbs)
save(enriched_group10, file = "enriched_group10.RData")
plotEnrich(enriched_group10[[3]], showTerms = 20, numChar = 40, y = "Count", orderBy = "Adjusted.P.value")
ggsave("enriched_group10.pdf",device = "pdf",width = 5, height = 8)

enriched_group11 <- enrichr(deg_group11$gene, dbs)
save(enriched_group11, file = "enriched_group11.RData")
plotEnrich(enriched_group11[[3]], showTerms = 20, numChar = 40, y = "Count", orderBy = "Adjusted.P.value")
ggsave("enriched_group11.pdf",device = "pdf",width = 5, height = 8)



load("enriched_group1.RData")
load("enriched_group2.RData")
load("enriched_group3.RData")
load("enriched_group4.RData")
load("enriched_group5.RData")
load("enriched_group6.RData")
load("enriched_group7.RData")
load("enriched_group8.RData")
load("enriched_group9.RData")
load("enriched_group10.RData")
load("enriched_group11.RData")

fig1_go_data <- rbind(enriched_group1[[3]][1:5,], enriched_group2[[3]][1:5,], enriched_group3[[3]][1:5,], enriched_group4[[3]][1:5,], enriched_group5[[3]][1:5,],
                      enriched_group6[[3]][1:5,], enriched_group7[[3]][1:5,],enriched_group8[[3]][1:5,], enriched_group9[[3]][1:5,],enriched_group10[[3]][1:5,],
                      enriched_group11[[3]][1:5,])
fig1_go_data$group <- c(rep("Granular",5), rep("Spinous",5), rep("Basal",5), rep("Endothelial",5), rep("Fibroblast",5), 
                        rep("Hairfollicle",5), rep("Sebocyte",5), rep("Melanocyte",5), rep("Pericyte",5), rep("Immune",5), 
                        rep("Langerhans",5))


fig1_go_data_input <- data.frame(pathway= fig1_go_data$Term, value= -log10(fig1_go_data$Adjusted.P.value), group= fig1_go_data$group)
fig1_go_data_input$pathway <- gsub(" \\(.*)", "", fig1_go_data_input$pathway)

write.csv(fig1_go_data_input, file= "03_results/fig1_go_data_input_16sam.csv")



