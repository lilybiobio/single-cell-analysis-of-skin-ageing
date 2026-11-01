##### Public dataset : wang, basal cells. #####
# https://htmlpreview.github.io/?https://github.com/WangShuxiong/Human_Epi/blob/master/Main_Figure1_2.html

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

library(Seurat)
library(dplyr)
library(stringr)


## 1. Read the reference data

skin.ref <- qs_read("01_data/skin.ref.16_09252025.qs2")
FeaturePlot(skin.ref, features = "ATF3", order = F, slot = "data", reduction = "umap.rna", ncol = 2, max.cutoff = 10, cols = c("grey90", "darkred"), split.by = "age")

p1 <- DimPlot(skin.ref, reduction = "umap", group.by = "sample", label = T, label.size = 2, repel = TRUE) + ggtitle("Reference annotations") + NoLegend()
p2 <- DimPlot(skin.ref, reduction = "umap", group.by = "celltype", label = T, label.size = 2, repel = TRUE) + ggtitle("Reference annotations") + NoLegend()
p1 + p2


## 2. Read the wang data

datadir <- "/n/groups/church/Lily/Skin_aging/scRNA_china/03_analysis_add12/01_data/wang_basal/"

# 1) list the matrix files and extract donor IDs (e.g., Donor1..Donor5)
mat_files <- list.files(datadir, pattern = "matrix\\.mtx\\.gz$", full.names = FALSE)
m <- str_match(mat_files, "(GSM\\d+)_?(Donor\\d+)")   # cols: full, GSM..., DonorN
donors <- unique(m[,3])
donors <- donors[!is.na(donors)]

objs <- list()

# 2) read each donor (matrix, barcodes, genes/features)
for (d in donors) {
  mtx  <- list.files(datadir, pattern = paste0(d, ".*matrix\\.mtx\\.gz$"),   full.names = TRUE)
  bar  <- list.files(datadir, pattern = paste0(d, ".*barcodes\\.tsv\\.gz$"), full.names = TRUE)
  feat <- list.files(datadir, pattern = paste0(d, ".*(genes|features)\\.tsv\\.gz$"), full.names = TRUE)
  
  if (!(length(mtx)==1 && length(bar)==1 && length(feat)==1)) {
    stop("Missing or multiple files for ", d, "\nmtx: ", toString(mtx),
         "\nbar: ", toString(bar), "\nfeat: ", toString(feat))
  }
  
  mat <- ReadMtx(mtx = mtx, features = feat, cells = bar)
  so  <- CreateSeuratObject(mat, project = d, min.cells = 3, min.features = 200)
  so$donor <- d
  objs[[d]] <- so
}

# 3) merge donors (optional)
skin <- Reduce(function(x, y) merge(x, y), objs)

skin
table(skin$donor)


skin.query.wang <- skin


## 3. Reference -- mapping

seuset <- skin.query.wang
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

p3

ggsave("umap.reference_wang.pdf", width = 8, height = 8)

qs_save(skin.query, "01_data/skin.query.wang.qs2")



##### get the github data #######

library(R.matlab)

# Read the .mat file
wang_data <- readMat("/n/groups/church/Lily/Skin_aging/scRNA_china/03_analysis_add12/01_data/wang_basal/Data_Figure1_and_Figure2.mat")

# See what’s inside
names(data)

# 1. Extract the list from the MATLAB structure
library(Seurat)
library(Matrix)

raw <- wang_data$Data.Figure1.and.Figure2

# Turn the 10 × 1 × 1 array into a named list:
wang_list <- setNames(
  lapply(rownames(raw), function(nm) raw[nm, 1, 1][[1]]),
  rownames(raw)
)

str(wang_list, max.level = 1)  # just to sanity check


# 2. Build the Seurat object from counts
counts <- wang_list$data

# Optionally set gene names from allgenes if missing
genes <- unlist(wang_list$allgenes)
if (!is.null(genes) && length(genes) == nrow(counts)) {
  rownames(counts) <- make.unique(as.character(genes))
}

seu <- CreateSeuratObject(counts = counts, project = "wang_basal")


# 3. Add cluster labels + batch as metadata

n_cells <- ncol(seu)

if (length(wang_list$cluster.label) == n_cells) {
  seu$cluster_label <- factor(wang_list$cluster.label)
}

if (length(wang_list$batch.id) == n_cells) {
  seu$batch_id <- factor(wang_list$batch.id)
}



# 4. Add the latent embedding as a DimReduc

latent <- wang_list$latent

if (length(latent) == 2 * n_cells) {
  # cells × dims
  latent_mat <- matrix(latent, nrow = n_cells, ncol = 2, byrow = FALSE)
  colnames(latent_mat) <- c("latent_1", "latent_2")
  rownames(latent_mat) <- colnames(seu)
  
  seu[["latent"]] <- CreateDimReducObject(
    embeddings = latent_mat,
    key = "LAT_",
    assay = DefaultAssay(seu)
  )
}

DimPlot(seu, reduction = "latent", group.by = "cluster_label")

qs_save(seu, "01_data/skin.query.wang.figure1and2.qs2")





### ------------ Map our 16sample data to Wang's data ---------------##################

skin.ref.wang <- qs_read("01_data/skin.query.wang.figure1and2.qs2")

skin.query.16sam <- qs_read("01_data/skin.query.16_09252025.qs2")


# 1. reference data
skin.ref <- skin.ref.wang
skin.ref$celltype <- skin.ref@meta.data$cluster_label
DefaultAssay(skin.ref) <- "RNA"

skin.ref <- NormalizeData(skin.ref, verbose = F)
skin.ref <- FindVariableFeatures(skin.ref, verbose = F)
skin.ref <- ScaleData(skin.ref, verbose = F)
skin.ref <- RunPCA(skin.ref, verbose = F)
skin.ref <- RunUMAP(skin.ref, dims = 1:40, reduction = "pca", return.model = TRUE)

cluster_colors <- c(
  "1"  = "#6EC6FF",
  "2" = "#D95F02",
  "3"= "#E6AB02",
  "4" = "#66A61E",
  "5"    = "#FFD54F",
  "6"    = "#7570B3",
  "7"    = "#1E88E5"
)
DimPlot(skin.ref, group.by = "celltype", cols = cluster_colors)
ggsave("02_figures/review_fig/umap.only_wang.pdf", width = 6.5, height = 6.5)

p2 <- DimPlot(skin.ref, reduction = "umap", group.by = "celltype", label = T) + ggtitle("Reference annotations") + NoLegend()
qs_save(skin.ref, "01_data/skin.query.wang.ready.qs2")


# 2. query data

skin.query <- skin.query.16sam


skin.anchors <- FindTransferAnchors(reference = skin.ref, query = skin.query, dims = 1:40, reference.reduction = "pca", reduction = "rpca", k.anchor = 20)

##  Unimodal UMAP Projection
LayerData(skin.ref[["RNA"]], layer = "counts") <- GetAssayData(skin.ref, layer = "data")

skin.query <- MapQuery(anchorset = skin.anchors, reference = skin.ref, query = skin.query,
                       refdata = list(cluster = "seurat_clusters", celltype = "celltype"), 
                       reference.reduction = "pca", reduction.model = "umap")


p3 <- DimPlot(skin.query, reduction = "ref.umap", group.by = "predicted.celltype", label = TRUE,
              label.size = 2, repel = TRUE, raster=FALSE) + NoLegend() + ggtitle("Query transferred labels")

p2 + p3








### ------------ Combined mapping  ---------------##################

skin.ref$dataset <- "ref"
skin.query$dataset <- "query"

combined <- merge(skin.ref, y = skin.query)

DefaultAssay(combined) <- "RNA"  # or "SCT", just be consistent
combined <- NormalizeData(combined)
combined <- FindVariableFeatures(combined)
combined <- ScaleData(combined)
combined <- RunPCA(combined, npcs = 50)
combined <- RunUMAP(
  combined,
  dims = 1:30,
  reduction.name = "latent",
  reduction.key = "latent_"
)

DimPlot(combined, reduction = "latent", group.by = "dataset")
DimPlot(combined, reduction = "latent", group.by = "celltype")  # if annotated


ggsave("02_figures/review_fig/umap.combine_wang.pdf", width = 8, height = 8)
qs_save(combined, "01_data/skin.query.wang.combine.qs2")





### ------------  Map our data (kera + mela) to Wang's data     ---------------##################


## 1. Read the ref data

skin.ref <- qs_read("01_data/skin.query.wang.ready.qs2")


cluster_colors <- c(
  "1"  = "#6EC6FF",
  "2" = "#D95F02",
  "3"= "#E6AB02",
  "4" = "#66A61E",
  "5"    = "#FFD54F",
  "6"    = "#7570B3",
  "7"    = "#1E88E5"
)
DimPlot(skin.ref, group.by = "celltype", cols = cluster_colors, label.size = 2)

ggsave("02_figures/ref.our_to_wang_ref_pc10.pdf",width = 7,height = 7)



## 2. Read the query data (our data (kera + mela))

our_data <- qs_read("01_data/skin.query_use_v2.qs2")

Idents(our_data) <- "predicted.celltype"

kera_mela_labels <- c("Basal cell-1","Basal cell-2","Basal cell-3",
                         "Basal cell-4","Basal cell-5","Basal cell-6","Basal cell-7",
                         "Spinous cell-1","Spinous cell-2","Spinous cell-3",
                         "Granular cell", "Melanocyte")

sub_kera_mela <- subset(our_data, idents = kera_mela_labels)

skin.query <- sub_kera_mela



## 3. Reference - mapping 

skin.anchors <- FindTransferAnchors(reference = skin.ref, query = skin.query, dims = 1:40, reference.reduction = "pca", reduction = "rpca", k.anchor = 20)

##  Unimodal UMAP Projection
LayerData(skin.ref[["RNA"]], layer = "counts") <- GetAssayData(skin.ref, layer = "data")

skin.query <- MapQuery(anchorset = skin.anchors, reference = skin.ref, query = skin.query,
                       refdata = list(cluster = "seurat_clusters", celltype = "celltype"), 
                       reference.reduction = "pca", reduction.model = "umap")

p3 <- DimPlot(skin.query, reduction = "ref.umap", group.by = "predicted.celltype", label = TRUE,
              label.size = 2, repel = TRUE, raster=FALSE, cols = cluster_colors) + NoLegend() + ggtitle("Query transferred labels")


p4 <- DimPlot(skin.query, reduction = "ref.umap", group.by = "our_label", label = TRUE,
              label.size = 2, repel = TRUE, raster=FALSE, cols = custom_colors) + NoLegend() + ggtitle("Query transferred labels")


ggsave("02_figures/ref.our_to_wang_query.pdf",width = 7,height = 7)
ggsave("02_figures/ref.our_to_wang_our_color.pdf",width = 7,height = 7)


qs_save(skin.ref, "01_data/skin.ref.our_to_wang.qs2")
qs_save(skin.query, "01_data/skin.query.our_to_wang.qs2")

skin.ref <- qs_read("01_data/skin.ref.our_to_wang.qs2")
skin.query <- qs_read("01_data/skin.query.our_to_wang.qs2")

meta.data <- skin.query@meta.data
meta.data$our_label <- sub_kera_mela@meta.data$predicted.celltype
skin.query@meta.data <- meta.data


## 4. Calculate the label accuracy

library(dplyr)

# Suppose:
obj1 <- skin.query
obj2 <- sub_kera_mela

# 1. Extract meta.data from both
meta1 <- obj1@meta.data %>%
  select(predicted.celltype) %>%
  mutate(barcode = rownames(.))

meta2 <- obj2@meta.data %>%
  select(predicted.celltype) %>%
  mutate(barcode = rownames(.))

# 2. Merge by barcode (assuming barcodes overlap)
merged <- inner_join(meta1, meta2, by = "barcode",
                     suffix = c(".num", ".annot"))

# 3. Calculate composition of each numeric cluster
celltype_mapping <- merged %>%
  group_by(predicted.celltype.num, predicted.celltype.annot) %>%
  summarise(n = n(), .groups = "drop_last") %>%
  mutate(prop = n / sum(n)) %>%
  arrange(predicted.celltype.num, desc(prop))

# 4. Set custom order of annotated cell types ---
ordered_levels <- c(
  "Basal cell-1", "Basal cell-2", "Basal cell-3",
  "Basal cell-4", "Basal cell-5", "Basal cell-6", "Basal cell-7",
  "Spinous cell-1", "Spinous cell-2", "Spinous cell-3",
  "Granular cell", "Melanocyte"
)

celltype_mapping$predicted.celltype.annot <- factor(
  celltype_mapping$predicted.celltype.annot,
  levels = ordered_levels
)


print(celltype_mapping)
write.csv(celltype_mapping, "celltype_mapping_proportions.csv", row.names = FALSE)


# 5. plot
custom_colors <- c(
  "#EBC7E8", "#F0D9FF", "#A084CA", "#6867AC",
  "#BFACE0", "#B980F0", "#A267AC",
  "#FFC1F3", "#E6739F", "#BF5CAA",
  "#D59DC5", "#F07DEA"
)

# --- 8. Plot with fixed order and color ---
p <- ggplot(celltype_mapping, aes(x = factor(predicted.celltype.num),
                                  y = prop,
                                  fill = predicted.celltype.annot)) +
  geom_bar(stat = "identity", color = "white", linewidth = 0.3) +
  scale_fill_manual(values = custom_colors, drop = FALSE) +
  labs(
    x = "Numeric cluster (from obj1)",
    y = "Proportion of annotated cell types (from obj2)",
    fill = "Annotated cell type"
  ) +
  theme_bw(base_size = 12) +
  theme(
    axis.text.x = element_text(angle = 45, hjust = 1),
    panel.grid.minor = element_blank(),
    panel.grid.major.x = element_blank(),
    legend.position = "right",
    legend.key.size = unit(0.5, "cm"),
    legend.text = element_text(size = 10),
    panel.border = element_rect(color = "black", linewidth = 0.6)
  )

ggsave("02_figures/celltype_mapping_proportions.pdf", p, width = 6, height = 6)
ggsave("02_figures/celltype_mapping_proportions.png", p, width = 6, height = 6, dpi = 300)







#### Map our dataset (kera + mela) to Wang's dataset


### ------------ Map our 16sample data to Wang's data ---------------##################

skin.ref.wang <- qs_read("01_data/skin.query.wang.figure1and2.qs2")


# 1. reference data
skin.ref <- skin.ref.wang
skin.ref$celltype <- skin.ref@meta.data$cluster_label
DefaultAssay(skin.ref) <- "RNA"

skin.ref <- NormalizeData(skin.ref, verbose = F)
skin.ref <- FindVariableFeatures(skin.ref, verbose = F)
skin.ref <- ScaleData(skin.ref, verbose = F)
skin.ref <- RunPCA(skin.ref, verbose = F)
skin.ref <- RunUMAP(skin.ref, dims = 1:10, reduction = "pca", return.model = TRUE)

cluster_colors <- c(
  "1"  = "#6EC6FF",
  "2" = "#D95F02",
  "3"= "#E6AB02",
  "4" = "#66A61E",
  "5"    = "#FFD54F",
  "6"    = "#7570B3",
  "7"    = "#1E88E5"
)
DimPlot(skin.ref, group.by = "celltype", cols = cluster_colors)
ggsave("02_figures/review_fig/umap.only_wang_pc10.pdf", width = 6.5, height = 6.5)

p2 <- DimPlot(skin.ref, reduction = "umap", group.by = "celltype", label = T) + ggtitle("Reference annotations") + NoLegend()
qs_save(skin.ref, "01_data/skin.query.wang.ready.pc10.qs2")




## 2. Read the query data (our data (kera + mela))

our_data <- qs_read("01_data/skin.query_use_v2.qs2")

Idents(our_data) <- "predicted.celltype"

kera_mela_labels <- c("Basal cell-1","Basal cell-2","Basal cell-3",
                      "Basal cell-4","Basal cell-5","Basal cell-6","Basal cell-7",
                      "Spinous cell-1","Spinous cell-2","Spinous cell-3",
                      "Granular cell", "Melanocyte")

sub_kera_mela <- subset(our_data, idents = kera_mela_labels)

skin.query <- sub_kera_mela



## 3. Reference - mapping 

skin.anchors <- FindTransferAnchors(reference = skin.ref, query = skin.query, dims = 1:10, reference.reduction = "pca", reduction = "rpca", k.anchor = 20)

##  Unimodal UMAP Projection
LayerData(skin.ref[["RNA"]], layer = "counts") <- GetAssayData(skin.ref, layer = "data")

skin.query <- MapQuery(anchorset = skin.anchors, reference = skin.ref, query = skin.query,
                       refdata = list(cluster = "seurat_clusters", celltype = "celltype"), 
                       reference.reduction = "pca", reduction.model = "umap")

p3 <- DimPlot(skin.query, reduction = "ref.umap", group.by = "predicted.celltype", label = TRUE,
              label.size = 2, repel = TRUE, raster=FALSE, cols = cluster_colors) + NoLegend() + ggtitle("Query transferred labels")


meta.data <- skin.query@meta.data
meta.data$our_label <- sub_kera_mela@meta.data$predicted.celltype
skin.query@meta.data <- meta.data


custom_colors <- c(
  "#EBC7E8", "#F0D9FF", "#A084CA", "#6867AC",
  "#BFACE0", "#B980F0", "#A267AC",
  "#FFC1F3", "#E6739F", "#BF5CAA",
  "#D59DC5", "#F07DEA"
)


p4 <- DimPlot(skin.query, reduction = "ref.umap", group.by = "our_label", label = TRUE,
              label.size = 2, repel = TRUE, raster=FALSE, cols = custom_colors) + NoLegend() + ggtitle("Query transferred labels")


ggsave("02_figures/ref.our_to_wang_query_pc10.pdf",width = 7,height = 7)
ggsave("02_figures/ref.our_to_wang_our_color_pc10.pdf",width = 7,height = 7)


qs_save(skin.ref, "01_data/skin.ref.our_to_wang.qs2")
qs_save(skin.query, "01_data/skin.query.our_to_wang_pc10.qs2")

skin.ref <- qs_read("01_data/skin.ref.our_to_wang.qs2")
skin.query <- qs_read("01_data/skin.query.our_to_wang.qs2")



## 4. Calculate the label accuracy

library(dplyr)

# Suppose:
obj1 <- skin.query
obj2 <- sub_kera_mela

# 1. Extract meta.data from both
meta1 <- obj1@meta.data %>%
  select(predicted.celltype) %>%
  mutate(barcode = rownames(.))

meta2 <- obj2@meta.data %>%
  select(predicted.celltype) %>%
  mutate(barcode = rownames(.))

# 2. Merge by barcode (assuming barcodes overlap)
merged <- inner_join(meta1, meta2, by = "barcode",
                     suffix = c(".num", ".annot"))

# 3. Calculate composition of each numeric cluster
celltype_mapping <- merged %>%
  group_by(predicted.celltype.num, predicted.celltype.annot) %>%
  summarise(n = n(), .groups = "drop_last") %>%
  mutate(prop = n / sum(n)) %>%
  arrange(predicted.celltype.num, desc(prop))

# 4. Set custom order of annotated cell types ---
ordered_levels <- c(
  "Basal cell-1", "Basal cell-2", "Basal cell-3",
  "Basal cell-4", "Basal cell-5", "Basal cell-6", "Basal cell-7",
  "Spinous cell-1", "Spinous cell-2", "Spinous cell-3",
  "Granular cell", "Melanocyte"
)

celltype_mapping$predicted.celltype.annot <- factor(
  celltype_mapping$predicted.celltype.annot,
  levels = ordered_levels
)


print(celltype_mapping)
write.csv(celltype_mapping, "01_data/celltype_mapping_proportions_pc10.csv", row.names = FALSE)


# 5. plot
custom_colors <- c(
  "#EBC7E8", "#F0D9FF", "#A084CA", "#6867AC",
  "#BFACE0", "#B980F0", "#A267AC",
  "#FFC1F3", "#E6739F", "#BF5CAA",
  "#D59DC5", "#F07DEA"
)

# --- 8. Plot with fixed order and color ---
p <- ggplot(celltype_mapping, aes(x = factor(predicted.celltype.num),
                                  y = prop,
                                  fill = predicted.celltype.annot)) +
  geom_bar(stat = "identity", color = "white", linewidth = 0.3) +
  scale_fill_manual(values = custom_colors, drop = FALSE) +
  labs(
    x = "Numeric cluster (from obj1)",
    y = "Proportion of annotated cell types (from obj2)",
    fill = "Annotated cell type"
  ) +
  theme_bw(base_size = 12) +
  theme(
    axis.text.x = element_text(angle = 45, hjust = 1),
    panel.grid.minor = element_blank(),
    panel.grid.major.x = element_blank(),
    legend.position = "right",
    legend.key.size = unit(0.5, "cm"),
    legend.text = element_text(size = 10),
    panel.border = element_rect(color = "black", linewidth = 0.6)
  )

ggsave("02_figures/celltype_mapping_proportions_pc10.pdf", p, width = 6, height = 6)
ggsave("02_figures/celltype_mapping_proportions.png", p, width = 6, height = 6, dpi = 300)










