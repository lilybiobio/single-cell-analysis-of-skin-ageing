#### Subcluster -- Fibroblast


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
#library(metap)
library(RColorBrewer)
library(Seurat)
#library(SeuratDisk)
library(ggplot2)
library(patchwork)
library(pheatmap)
library(RColorBrewer)
#library(viridis)
library(qs2)


###### 2025, fibroblast in 16 samples  ########

markers_fibroblast_clusterwise <- read.table("../02_analysis/results/markers/markers_fibroblast_clusterwise.txt")
skin.ref <- qs_read("01_data/skin.ref.16_09252025.qs2")
skin.query <- qs_read("01_data/skin.query.16_09252025.qs2")

## fibroblast changes in different cell types 

pap_genes <- c("APCDD1","WIF1","COL23A1","PDPN","NTN1","COL6A6","OLFML3","LUM","DCN")
ret_genes <- c("DPP4","SFRP2","MFAP5","FBLN2","ELN","COL1A1","COL1A2","MMP2","AEBP1")
fib <- subset(skin.ref, subset = celltype %in% c("Fibroblast-1", "Fibroblast-2", "Fibroblast-3", "Fibroblast-4"))
fib_16sam <- subset(skin.query, subset = predicted.celltype %in% c("Fibroblast-1", "Fibroblast-2", "Fibroblast-3", "Fibroblast-4"))

fib <- AddModuleScore(fib, features = list(pap_genes), name = "Papillary_Score")
fib <- AddModuleScore(fib, features = list(ret_genes), name = "Reticular_Score")




#  Give short, readable labels to match the figure idea
fib$fb_subtype <- recode(
  fib$celltype,
  "Fibroblast-1" = "Reticular",      # S.R in the screenshot
  "Fibroblast-2" = "Papillary",      # S.P
  "Fibroblast-3" = "Inflammatory",   # INF
  "Fibroblast-4" = "Myofibroblast"   # MES in that figure
)


# create the age groups 
fib$age_group <- case_when(
  fib$sample %in% c("age23a", "age23b") ~ "YOUNG",
  fib$sample %in% c("age52a", "age52b") ~ "MIDDLE",
  fib$sample %in% c("age85a", "age85b") ~ "OLD",
)

# Make it a factor for consistent plot order
fib$age_group <- factor(fib$age_group, levels = c("YOUNG", "MIDDLE", "OLD"))


# violin plot
library(dplyr)
library(ggplot2)
library(patchwork)
# ----------------------------------------------------------
# 0) Shift each score so the minimum = 0 (avoid negative values)
# ----------------------------------------------------------
shift_features <- c("Papillary_Score", "Reticular_Score")

for (feat in shift_features) {
  shifted_name <- paste0(feat, "_shifted")
  fib_yo[[shifted_name]] <- fib_yo[[feat]] - min(fib_yo[[feat]], na.rm = TRUE)
}

# ----------------------------------------------------------
# 1) Helper function: violin plot with per-(subtype × age_group) medians
# ----------------------------------------------------------
vln_with_median <- function(obj, feature, ylab = feature, dodge = 0.9, show_labels = TRUE) {
  # Compute medians for the shifted feature
  df_med <- FetchData(obj, vars = c("fb_subtype", "age_group", feature)) %>%
    dplyr::filter(age_group %in% c("YOUNG", "OLD")) %>%
    dplyr::group_by(fb_subtype, age_group) %>%
    dplyr::summarize(median_val = median(.data[[feature]], na.rm = TRUE), .groups = "drop")
  
  # Base violin (split by age_group)
  p <- VlnPlot(
    obj, features = feature, group.by = "fb_subtype", split.by = "age_group",
    pt.size = 0, assay = DefaultAssay(obj)
  ) +
    theme_classic(base_size = 12) +
    labs(x = NULL, y = ylab)
  
  # Overlay median point (+ optional numeric label)
  p +
    geom_point(
      data = df_med, inherit.aes = FALSE,
      aes(x = fb_subtype, y = median_val, shape = age_group),
      position = position_dodge(width = dodge),
      size = 2.3, stroke = 0.6, color = "black", fill = "white"
    ) +
    { if (show_labels)
      geom_text(
        data = df_med, inherit.aes = FALSE,
        aes(x = fb_subtype, y = median_val, label = sprintf("%.2f", median_val)),
        position = position_dodge(width = dodge),
        vjust = -0.8, size = 3.0, color = "black"
      ) else NULL } +
    scale_shape_manual(values = c(YOUNG = 21, OLD = 24)) +  # circle vs triangle
    guides(shape = "none")
}

# ----------------------------------------------------------
# 2) Build the two shifted panels
# ----------------------------------------------------------
v_pap_median  <- vln_with_median(fib_yo, "Papillary_Score_shifted", ylab = "Papillary signature (shifted)")
v_ret_median  <- vln_with_median(fib_yo, "Reticular_Score_shifted", ylab = "Reticular signature (shifted)")

# ----------------------------------------------------------
# 3) Combine them
# ----------------------------------------------------------
v_combined_shifted <- ( v_pap_median / v_ret_median) + plot_layout(heights = c(1,1))

v_combined_shifted

ggsave("fibroblast_violin_with_median.pdf",
       width = 7, height = 7, units = "in", dpi = 600)



#### median value table  ###-------------
library(dplyr)
library(purrr)

# Define the shifted score names
scores_shifted <- c("Papillary_Score_shifted", "Reticular_Score_shifted")

# Function to compute medians for one shifted feature
get_median_table_shifted <- function(feature) {
  FetchData(fib_yo, vars = c("fb_subtype", "age_group", feature)) %>%
    dplyr::filter(age_group %in% c("YOUNG", "OLD")) %>%
    dplyr::group_by(fb_subtype, age_group) %>%
    dplyr::summarize(median = median(.data[[feature]], na.rm = TRUE), .groups = "drop") %>%
    dplyr::mutate(Feature = feature)
}

# Combine both scores
median_table_shifted <- purrr::map_dfr(scores_shifted, get_median_table_shifted)

# Reorder columns
median_table_shifted <- median_table_shifted %>%
  dplyr::select(Feature, fb_subtype, age_group, median) %>%
  dplyr::arrange(Feature, fb_subtype, age_group)

# View the shifted-median table
median_table_shifted

write.csv(median_table_shifted, "fibroblast_violin_medians_shifted.csv", row.names = FALSE)
write.csv(summary_ttest_table, "fibroblast_violin_ttest_summary.csv", row.names = FALSE)


qs_save(fib, "01_data/fib.qs2")

qs_save(fib_yo, "01_data/fib_yo.qs2")




