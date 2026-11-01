
# ----------- Load libraries -------------
library(Seurat)
library(Matrix)
library(stringr)
library(future)
options(future.globals.maxSize = 4000 * 1024^2)

# ----------- Step 1: Define 18 sample directories and sample IDs -------------
ages <- c(23, 24, 28, 46, 52, 55, 76, 85)
sample_dirs <- as.vector(sapply(ages, function(age) {
  c(paste0("age", age, "a_filtered_feature_bc_matrix"),
    paste0("age", age, "b_filtered_feature_bc_matrix"))
}))
data_root <- "/n/groups/church/Lily/Skin_aging/scRNA_china/03_analysis_add12/01_data/input_filtered_feature_addseq/"
dirs <- file.path(data_root, sample_dirs)
sample_ids <- str_replace(sample_dirs, "_filtered_feature_bc_matrix", "")



# ----------- Step 2: Read all samples into one Seurat object (v5 layers) -------------
seurat_list <- mapply(function(path, id) {
  mat <- Read10X(data.dir = path)
  CreateSeuratObject(counts = mat, min.features = 100, project = id)
}, dirs, sample_ids, SIMPLIFY = FALSE)


merged_seurat <- merge(
  x = seurat_list[[1]],
  y = seurat_list[-1],
  add.cell.ids = sample_ids
)

merged_seurat$sample <- sapply(strsplit(colnames(merged_seurat), "_"), `[`, 1)


# ----------- Step 3: Add QC metrics and filter cells -------------
merged_seurat$log10GenesPerUMI <- log10(merged_seurat$nFeature_RNA) / log10(merged_seurat$nCount_RNA)
merged_seurat$mitoRatio <- PercentageFeatureSet(merged_seurat, pattern = "^MT-") / 100

filtered_seurat <- subset(
  merged_seurat,
  subset = nCount_RNA >= 200 & nFeature_RNA >= 200 & nFeature_RNA <= 6000 & mitoRatio < 0.10
)

save(filtered_seurat, file = "01_data/seurat_filtered_all16samples_v2.RData")


# ----------- Step 4: Filter genes (≥3 cells expressing) -------------
# counts" layer 

layers <- grep("^counts", Layers(filtered_seurat[["RNA"]]), value = TRUE)

gene_filter_list <- list()

for (ly in layers) {
  counts_mat <- LayerData(filtered_seurat[["RNA"]], layer = ly)
  keep_genes <- rowSums(counts_mat > 0) >= 3
  gene_filter_list[[ly]] <- rownames(counts_mat)[keep_genes]
}

genes_to_keep <- Reduce(union, gene_filter_list)

filtered_seurat2 <- filtered_seurat[genes_to_keep, ]


# ----------- Step 5: Save filtered object -------------
save(filtered_seurat2, file = "01_data/seurat_filtered_all16samples_v3.RData")




