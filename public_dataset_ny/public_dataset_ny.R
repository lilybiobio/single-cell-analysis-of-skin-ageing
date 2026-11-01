#### 3.8 Public datasets  #########

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



## 1. Read the reference data

skin.ref <- qs_read("01_data/skin.ref.16_09252025.qs2")

p1 <- DimPlot(skin.ref, reduction = "umap", group.by = "sample", label = T, label.size = 2, repel = TRUE) + ggtitle("Reference annotations") + NoLegend()
p2 <- DimPlot(skin.ref, reduction = "umap", group.by = "celltype", label = T, label.size = 2, repel = TRUE) + ggtitle("Reference annotations") + NoLegend()
p1 + p2


## 2. Read the public data

skin.query.public <- qs_read("01_data/seuset_public.qs2")


## 3. Reference -- mapping

seuset <- skin.query.public
#seuset <- JoinLayers(seuset)
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

ggsave("Umap.reference_public.png", width=24, height=8, units = "in")

qs_save(skin.query, "01_data/skin.query.public.qs2")

skin.query <- qs_read("01_data/skin.query.public.qs2")




#### ------  Plot ATF3 in public dataset ------ ########

## 1. Plot in differnt cell types


library(Seurat)
library(dplyr)
library(readxl)
library(stringr)
library(ggplot2)
library(tidyr)
library(purrr)

## ---- 1) read Excel and prepare donor_age ----
# Expect the Excel to have a column 'donor_age' (numeric) and an ID to join on
# such as 'sample_barcode' or 'donor_id'.
sample_df <- read_excel("01_data/scrna_samplemetadata.xlsx")

# Clean names
names(sample_df) <- names(sample_df) |>
  str_replace_all("\\s+", "_") |>
  tolower()

# Ensure donor_age is numeric
if (!"donor_age" %in% names(sample_df)) {
  stop("Excel must contain a 'donor_age' column.")
}
sample_df <- sample_df |> mutate(donor_age = suppressWarnings(as.numeric(donor_age)))

## 1) collapse Excel to one row/key and join donor_age (robust to duplicates)
# choose join key automatically
md <- skin.query@meta.data
cn_md  <- tolower(colnames(md))
cn_xls <- tolower(colnames(sample_df))
pick_key <- function(cands) cands[which(cands %in% cn_md & cands %in% cn_xls)][1]
join_key <- pick_key(c("sample_barcode","donor_id"))
if (is.na(join_key)) stop("No common join key found (looked for 'sample_barcode' or 'donor_id').")

# clean + numeric age
names(sample_df) <- names(sample_df) |> str_replace_all("\\s+","_") |> tolower()
sample_df <- sample_df |> mutate(donor_age = suppressWarnings(as.numeric(donor_age)))

# collapse to one row per key (use median age if multiple)
sample_key <- sample_df |>
  filter(!is.na(.data[[join_key]])) |>
  group_by(.data[[join_key]]) |>
  summarise(donor_age = median(donor_age, na.rm = TRUE), .groups = "drop")

# join (many-to-one, won’t duplicate cells)
md2 <- md |>
  tibble::rownames_to_column("cell_id") |>
  left_join(sample_key, by = setNames(join_key, join_key)) |>
  tibble::column_to_rownames("cell_id")

skin.query@meta.data <- md2
# keep only cells with donor_age present
skin.query <- subset(skin.query, subset = !is.na(donor_age))

## 2) sanity checks on required meta columns
need_meta <- c("donor_age","cell_type","donor_id")
missing_meta <- setdiff(need_meta, colnames(skin.query@meta.data))
if (length(missing_meta)) {
  stop("Missing metadata columns in Seurat object after join: ",
       paste(missing_meta, collapse = ", "))
}

## 3) find ATF3 and the assay that actually contains it
assays <- Assays(skin.query)
find_atf3 <- function(obj) {
  hits <- lapply(assays, function(a) {
    rn <- rownames(obj[[a]])
    idx <- which(toupper(rn) == "ATF3")
    if (length(idx)) data.frame(assay=a, gene=rn[idx], stringsAsFactors=FALSE) else NULL
  })
  do.call(rbind, hits)
}
hit <- find_atf3(skin.query)
if (is.null(hit) || nrow(hit) == 0) stop("ATF3 not found in any assay rownames.")

assay_use <- hit$assay[1]
gene_use  <- hit$gene[1]
DefaultAssay(skin.query) <- assay_use

## 4) fetch per-cell data safely (FetchData will now find all variables)
vars_use <- c("donor_age","cell_type","donor_id", gene_use)
df <- FetchData(skin.query, vars = vars_use)
colnames(df)[colnames(df) == gene_use] <- "ATF3"

## 5) pseudobulk per donor × cell_type (reduces single-cell noise)
df_avg <- df |>
  group_by(cell_type, donor_id) |>
  summarise(
    donor_age = suppressWarnings(as.numeric(first(na.omit(donor_age)))),
    mean_expr = mean(ATF3, na.rm = TRUE),
    prop_expressing = mean(ATF3 > 0, na.rm = TRUE),
    n_cells = dplyr::n(),
    .groups = "drop"
  ) |>
  filter(!is.na(donor_age))






### 2. Put KCs together, plot KC


library(dplyr)
library(ggplot2)

## 1) Merge the three KC subtypes into a single group and aggregate per donor
kc_levels <- c("Bas KC", "Cyc KC", "Diff KC")

df_kc <- df_avg %>%
  filter(cell_type %in% kc_levels) %>%
  group_by(donor_id) %>%
  summarise(
    donor_age = suppressWarnings(as.numeric(first(na.omit(donor_age)))),
    # Weighted mean across KC subtypes (weights = number of cells per subtype)
    mean_expr = if (all(!is.na(n_cells))) {
      weighted.mean(mean_expr, w = n_cells, na.rm = TRUE)
    } else {
      mean(mean_expr, na.rm = TRUE)  # fallback if n_cells missing
    },
    n_cells_total = sum(n_cells, na.rm = TRUE),
    .groups = "drop"
  ) %>%
  filter(!is.na(donor_age))


## 2) Fit linear regression
kc_lm <- lm(mean_expr ~ donor_age, data = df_kc)
kc_sum <- summary(kc_lm)
r2_txt <- paste0("R² = ", signif(kc_sum$r.squared, 3))
pval_txt <- paste0("p = ", signif(coef(kc_sum)[2,4], 3))
slope_txt <- paste0("slope = ", signif(coef(kc_lm)[2], 3))

## 3) Regression plot
p_kc <- ggplot(df_kc, aes(x = donor_age, y = mean_expr)) +
  geom_point(size = 2, alpha = 0.85) +
  geom_smooth(method = "lm", se = TRUE, linewidth = 0.9) +
  labs(
    title = "ATF3 expression vs donor age — Keratinocytes (Bas KC + Cyc KC + Diff KC)",
    x = "Donor age",
    y = "Mean log-normalized ATF3 expression (per donor, merged KC)"
  ) +
  theme_bw(base_size = 12) +
  annotate("text",
           x = Inf, y = -Inf, hjust = 1.02, vjust = -0.8, size = 3.6,
           label = paste(r2_txt, pval_txt, slope_txt, sep = "   "))

print(p_kc)

## Optional: save
ggsave("ATF3_vs_donor_age_Keratinocytes_regression_v2.pdf", p_kc, width = 5, height = 6.5)



### plot for Nature style #### 


library(ggplot2)

# Nature-inspired palette
nature_blue <- "#2C7BB6"   # cool blue for regression line
nature_fill <- "#A6CEE3"   # light blue for confidence interval

p_kc <- ggplot(df_kc, aes(x = donor_age, y = mean_expr)) +
  geom_point(color = "black", size = 2.4, alpha = 0.9, stroke = 0.2) +
  geom_smooth(
    method = "lm", se = TRUE,
    color = nature_blue, fill = alpha(nature_fill, 0.3),
    linewidth = 1.3, linetype = "solid"
  ) +
  labs(
    title = "ATF3 expression versus donor age in keratinocytes",
    subtitle = "Combined basal, cycling, and differentiated KC populations",
    x = "Donor age (years)",
    y = "Mean log-normalized ATF3 expression per donor"
  ) +
  theme_classic(base_size = 12, base_family = "Arial") +
  theme(
    plot.title = element_text(face = "bold", size = 12, hjust = 0),
    plot.subtitle = element_text(size = 10, color = "black", hjust = 0),
    axis.title = element_text(size = 11, color = "black"),
    axis.text = element_text(size = 10, color = "black"),
    axis.line = element_line(linewidth = 0.6, color = "black"),  # Black axes
    panel.grid = element_blank(),
    panel.border = element_blank(),
    plot.margin = margin(6, 6, 6, 6, "pt"),
    text = element_text(family = "Arial")
  )

# Export publication-quality figure
ggsave(
  "ATF3_KC_Age_Nature_BlackAxis.pdf",
  p_kc,
  width = 3.5, height = 3,
  units = "in", useDingbats = FALSE
)




#### Public data + our dataset ###

skin.query <- qs_read("01_data/skin.query.public.qs2")

DimPlot(skin.query, group.by = "cell_type")



## run the code below
# add a label to each object so Harmony knows the batch
skin.ref$dataset   <- "ref"
skin.query$dataset <- "query"

# merge into one Seurat object
skin_merged <- merge(
  skin.ref,
  y = skin.query,
  add.cell.ids = c("ref", "query"),
  project = "skin_harmony"
)


# normalize and find HVGs
skin_merged <- NormalizeData(skin_merged)
skin_merged <- FindVariableFeatures(skin_merged, selection.method = "vst", nfeatures = 3000)

# scale & PCA
skin_merged <- ScaleData(skin_merged, features = VariableFeatures(skin_merged))
skin_merged <- RunPCA(skin_merged, features = VariableFeatures(skin_merged), npcs = 50)


## harmony
library(harmony)
seuset <- RunHarmony(seuset, "orig.ident")
seuset <- FindNeighbors(object = seuset, dims = 1:PCs, return.neighbor = F, reduction = "harmony") 
seuset <- FindClusters(object = seuset, resolution = 0.8, algorithm = 1)
seuset <- RunUMAP(seuset, reduction = "harmony", dims = 1:30, reduction.name = "umap.harmony")
# DimPlot(object = seuset, reduction = "umap.harmony", pt.size = 1, label = F, group.by="orig.ident", shuffle = T)
# ggsave("../RNA/rna.umap.batch.harmony.png", width = 8, height = 6, units = "in")


my_colors2 <- c(
  "Bas KC"     = "#A3C9E2",   # light sky blue
  "Cyc KC"     = "#8EE1DD",   # aqua pastel turquoise
  "Diff KC"    = "#7FAADE",   # fresh denim blue
  "Melano"     = "#CFB7E3",   # light lavender
  "HFE"        = "#9FCF9A",   # mint green
  "Ecc Gland"  = "#C5E6C3",   # pale green mint
  "Fib"        = "#7EC9A4",   # minty teal
  "EC"         = "#C1A6D8",   # pastel violet
  "LEC"        = "#9AB9E8",   # baby blue
  "Lym"        = "#8FB7F0",   # clean bright blue
  "Myel"       = "#F29AC2",   # soft candy pink
  "Peri"       = "#FFB3C6",   # light rose
  "Mast"       = "#E7C5DA",   # powdery mauve
  "APC"        = "#BA99C9",   # lilac purple
  "NA"         = "grey82"
)


DimPlot(object = seuset, reduction = "umap.harmony", pt.size = 1, label = T, group.by = "cell_type", cols = my_colors2 , shuffle = T)

ggsave("02_figures/public_plus_16sam_umap_harmony_celltype_label_v1.pdf", width = 7, height = 6)


qs_save(seuset, "01_data/public_and_our16sam.qs2")




### Plot 1: our data on the top 



# 1.1 set color

my_colors <- c(
  "Basal cell-1"             = "#EBC7E8",
  "Basal cell-2"             = "#F0D9FF",
  "Basal cell-3"             = "#A084CA",
  "Basal cell-4"             = "#6867AC",
  "Basal cell-5"             = "#BFACE0",
  "Basal cell-6"             = "#B980F0",
  "Basal cell-7"             = "#A267AC",
  "Spinous cell-1"           = "#FFC1F3",
  "Spinous cell-2"           = "#E6739F",
  "Spinous cell-3"           = "#BF5CAA",
  "Granular cell"            = "#D59DC5",
  "Melanocyte"               = "#F07DEA",
  "Langerhans"               = "#E89F71",
  "Hair follicle-1"          = "#A2B38B",
  "Hair follicle-2"          = "#94B4A4",
  "Hair follicle-3"          = "#EFFAD3",
  "Hair follicle-4"          = "#A8DDA8",
  "Hair follicle (Sebocyte)" = "#00AFC1",
  "Fibroblast-1"             = "#557B83",
  "Fibroblast-2"             = "#39AEA9",
  "Fibroblast-3"             = "#A2D5AB",
  "Fibroblast-4"             = "#E5EFC1",
  "Pericyte-1"               = "#F1D4D4",
  "Pericyte-2"               = "#FF9292",
  "Endothelial cell-1"       = "#D2DAFF",
  "Endothelial cell-2"       = "#B1B2FF",
  "Endothelial cell-3"       = "#AAC4FF",
  "Endothelial cell-4"       = "#9AD0EC",
  "Endothelial cell-5"       = "#9BA3EB",
  "Immune cell-1"            = "#FFD8A9",
  "Immune cell-2"            = "#EDCFA9"
)


# 1.2 force cell type follow the order
df$celltype <- trimws(as.character(df$celltype))

df$celltype <- factor(
  df$celltype,
  levels = names(my_colors)   # <- this order is now the legend/order
)

# 1.3 plot code

ggplot() +
  # grey background
  geom_point(
    data = subset(df, celltype == "NA" | is.na(celltype)),
    aes(umapharmony_1, umapharmony_2),
    color = "grey80",
    size  = 0.3
  ) +
  # colored cell types
  geom_point(
    data = subset(df, celltype != "NA" & !is.na(celltype)),
    aes(umapharmony_1, umapharmony_2, color = celltype),
    size = 0.5
  ) +
  scale_color_manual(values = my_colors, drop = FALSE) +
  theme_classic() +
  guides(color = guide_legend(override.aes = list(size = 4)))

ggsave("02_figures/public_plus_16sam_umap_ourup.pdf", width = 8, height = 5)



