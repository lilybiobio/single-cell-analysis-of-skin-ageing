library(Seurat)
library(sctransform)


load("data/seurat_filtered_all16samples_v3.RData")
seuset <- filtered_seurat_final
seuset <- JoinLayers(seuset)
DefaultAssay(seuset) <-"RNA"
seuset <- NormalizeData(seuset, normalization.method = "LogNormalize", scale.factor = 10000)
seuset <- FindVariableFeatures(seuset, selection.method = "vst", nfeatures = 2000)
feature <- VariableFeatures(object = seuset)
seuset <- ScaleData(seuset, features = feature)


# PCA
PCs = 30
seuset <- RunPCA(object = seuset, features = feature, npcs = PCs)

## Cluster the cells
seuset <- FindNeighbors(object = seuset, dims = 1:PCs, return.neighbor = F) 
seuset <- FindClusters(object = seuset, resolution = 0.8, algorithm = 1)

## Look at cluster IDs of the first 5 cells
head(x = Idents(object = seuset), 5)

# Umap
seuset <- RunUMAP(seed.use = 42, object = seuset, dims = 1:PCs, n.neighbors = 30, # 1:15
                  reduction.name = 'umap.rna', reduction.key = 'rnaUMAP_', min.dist = 0.3 # 0.05 umap.rna.p3 umap.rna
)
DimPlot(object = seuset, reduction = "umap.rna",  pt.size = 1, label = T) # + scale_color_manual(values=celltypecolors)
ggsave("../RNA/rna.umap.png", width = 8, height = 6, units = "in") 
DimPlot(object = seuset,  reduction = "umap.rna", pt.size = 1, label = F, group.by="orig.ident", shuffle = T)
ggsave("../RNA/rna.umap.batch.png", width = 8, height = 6, units = "in")

## harmony
library(harmony)
seuset <- RunHarmony(seuset, "orig.ident")
seuset <- FindNeighbors(object = seuset, dims = 1:PCs, return.neighbor = F, reduction = "harmony") 
seuset <- FindClusters(object = seuset, resolution = 0.8, algorithm = 1)
seuset <- RunUMAP(seuset, reduction = "harmony", dims = 1:30, reduction.name = "umap.harmony")
DimPlot(object = seuset, reduction = "umap.harmony", pt.size = 1, label = F, group.by="orig.ident", shuffle = T)
ggsave("../RNA/rna.umap.batch.harmony.png", width = 8, height = 6, units = "in")

DimPlot(object = seuset, reduction = "umap.harmony", pt.size = 1, label = T, group.by = "seurat_clusters", shuffle = T)
ggsave("../RNA/rna.umap.harmony.png", width = 8, height = 6, units = "in")



