### Reference -- mapping
library(Seurat)
library(sctransform)
library(SeuratWrappers)
library(patchwork)
library(Matrix)
library(stringr)
library(future)
library(rlang)
library(ggplot2)
library(pheatmap)
library(qs2)



# reference projection ----------------------------------------------------
load("/n/groups/church/Lily/Skin_aging/scRNA_china/03_analysis_add12/01_data/seurat_filtered_all16samples.RData")
DefaultAssay(filtered_seurat_final) <-"RNA"

# subset 16 levels
keep_levels <- levels(filtered_seurat_final)[1:16]

filtered_seurat_16 <- subset(filtered_seurat_final, idents = keep_levels)

levels(filtered_seurat_16)
table(filtered_seurat_16$orig.ident)

save(filtered_seurat_16, file = "01_data/seurat_filtered_all16samples_v4.RData")
DefaultAssay(filtered_seurat_16) <-"RNA"


### 6 sample object
skin.ref_sub <- readRDS("/n/groups/church/Lily/Skin_aging/scRNA_china/02_analysis/results/integrated11_seurat_6samples_sub.rds")
skin.ref_all <- readRDS("/n/groups/church/Lily/Skin_aging/scRNA_china/02_analysis/results/integrated11_seurat_6samples.rds")


meta_sub <- skin.ref_sub@meta.data
meta_all <- skin.ref_all@meta.data

skin.ref_sub$celltype <- skin.ref_sub@active.ident


skin.ref_sub <- RunUMAP(skin.ref_sub, dims = 1:40, reduction = "pca", return.model = TRUE)

# work on the skin.ref_all 

Idents(skin.ref_all) <- skin.ref_all$integrated_snn_res.0.6
skin.ref_all <- RenameIdents(object = skin.ref_all, 
                             cluster_annotation <- c(
                               "0"  = "Spinous cell-3",
                               "1"  = "Spinous cell-2",
                               "2"  = "Hair follicle-4",
                               "3"  = "Basal cell-6",
                               "4"  = "Endothelial cell-2",
                               "5"  = "Spinous cell-1",
                               "6"  = "Fibroblast-3",
                               "7"  = "Endothelial cell-3",
                               "8"  = "Fibroblast-2",
                               "9"  = "Basal cell-2",
                               "10" = "Basal cell-5",
                               "11" = "Pericyte-1",
                               "12" = "Hair follicle (Sebocyte)",
                               "13" = "Hair follicle-2",
                               "14" = "Granular cell",
                               "15" = "Basal cell-4",
                               "16" = "Immune cell-1",
                               "17" = "Basal cell-1",
                               "18" = "Fibroblast-1",
                               "19" = "Langerhans",
                               "20" = "Melanocyte",
                               "21" = "Fibroblast-4",
                               "22" = "Basal cell-7",
                               "23" = "Pericyte-2",
                               "24" = "Immune cell-2",
                               "25" = "Basal cell-3",
                               "26" = "Hair follicle-3",
                               "27" = "NA",
                               "28" = "Endothelial cell-4",
                               "29" = "Endothelial cell-5",
                               "30" = "Hair follicle-1",
                               "31" = "Endothelial cell-1"))

skin.ref_all$celltype <- skin.ref_all@active.ident
DimPlot(skin.ref_all, reduction = "umap", group.by = "celltype", label = T, label.size = 2, repel = TRUE) + ggtitle("Reference annotations") + NoLegend()
DimPlot(skin.ref_all, reduction = "umap", group.by = "integrated_snn_res.0.6", label = T, label.size = 2, repel = TRUE) + ggtitle("Reference annotations") + NoLegend()

save(skin.ref_all, file = "01_data/skin.ref_all_2025_v1.RData")

load("01_data/skin.ref_all_2025_v1.RData")


skin.ref <- skin.ref_all

skin.ref$celltype <- skin.ref@active.ident
DefaultAssay(skin.ref) <- "RNA"

skin.ref <- NormalizeData(skin.ref, verbose = F)
skin.ref <- FindVariableFeatures(skin.ref, verbose = F)
skin.ref <- ScaleData(skin.ref, verbose = F)
#skin.ref <- RunPCA(skin.ref, verbose = F)
skin.ref <- RunUMAP(skin.ref, dims = 1:40, reduction = "pca", return.model = TRUE)
qs_save(skin.ref, "01_data/skin.ref_all_2025_v2.qs2")



p1 <- DimPlot(skin.ref, reduction = "umap", group.by = "sample", label = T, label.size = 2, repel = TRUE) + ggtitle("Reference annotations") + NoLegend()
p2 <- DimPlot(skin.ref, reduction = "umap", group.by = "celltype", label = T, label.size = 2, repel = TRUE) + ggtitle("Reference annotations") + NoLegend()
p1 + p2
ggsave("Umap.reference_all.png", width=18, height=8, units = "in")
#qs_save(skin.ref, "skin.ref.qs2")




## predict cell type
load("01_data/seurat_filtered_all16samples_v4.RData")

seuset <- filtered_seurat_16 #[, 1:10000]
seuset <- filtered_seurat_16
seuset <- JoinLayers(seuset)
DefaultAssay(seuset) <-"RNA"
seuset <- NormalizeData(seuset, normalization.method = "LogNormalize", scale.factor = 10000)
seuset <- FindVariableFeatures(seuset, selection.method = "vst", nfeatures = 2000)
feature <- VariableFeatures(object = seuset)
seuset <- ScaleData(seuset, features = feature)
skin.query <- seuset

skin.anchors <- FindTransferAnchors(reference = skin.ref, query = skin.query, dims = 1:40, reference.reduction = "pca", reduction = "rpca", k.anchor = 20)

##  Unimodal UMAP Projection
LayerData(skin.ref[["RNA"]], layer = "counts") <- GetAssayData(skin.ref, layer = "data")

skin.query <- MapQuery(anchorset = skin.anchors, reference = skin.ref, query = skin.query,
                       refdata = list(cluster = "seurat_clusters", celltype = "celltype"), 
                       reference.reduction = "pca", reduction.model = "umap")

p3 <- DimPlot(skin.query, reduction = "ref.umap", group.by = "predicted.celltype", label = TRUE,
              label.size = 2, repel = TRUE, raster=FALSE) + NoLegend() + ggtitle("Query transferred labels")

p1 + p2 + p3
# qs_save(skin.query, "skin.query.qs2")
ggsave("Umap.reference_16sample.png", width=24, height=8, units = "in")

save(skin.ref, file = "01_data/skin.ref.16_09252025.RData")
save(skin.query, file = "01_data/skin.query.16_09252025.RData")


qs_save(skin.ref, "01_data/skin.ref.16_09252025.qs2")
qs_save(skin.query, "01_data/skin.query.16_09252025.qs2")


index <- skin.query$orig.ident == "age23b"
table(index)
head(names(skin.ref$celltype))
df <- table(skin.ref$celltype, skin.query$predicted.celltype[names(skin.ref$celltype)])
pheatmap(df, scale="row")


## subset, remove the cluster "NA"

skin.query_use <- subset(skin.query, subset = predicted.celltype != "NA")

qs_save(skin.query_use, "01_data/skin.query_use.qs2")










