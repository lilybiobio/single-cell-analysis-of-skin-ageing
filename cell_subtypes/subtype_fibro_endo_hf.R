#### 3.11 subtypes of fibroblast, endo, and HF

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
library(SeuratDisk)
library(ggplot2)
library(patchwork)
library(pheatmap)
library(RColorBrewer)
library(viridis)
library(qs2)

skin.query <- qs_read("01_data/skin.query_use_v2.qs2")

Idents(skin.query) <- "predicted.celltype"

sub_fib_16sam <- subset(skin.query, subset = predicted.celltype %in% c("Fibroblast-1", "Fibroblast-2", "Fibroblast-3", "Fibroblast-4"))
hf <- subset(obj, subset = predicted.celltype %in% c("Hair follicle-1", "Hair follicle-2", "Hair follicle-3", "Hair follicle-4"))


#### 1. Umaps for fibro, endo, and hf 

color_cluster=c('#557B83', '#39AEA9','#A2D5AB','#E5EFC1')

names(color_cluster)=c('Fibroblast-1',
                       'Fibroblast-2',
                       'Fibroblast-3',
                       'Fibroblast-4')

#group.by = 'orig.ident',

DimPlot(object = sub_fib_16sam, pt.size = 0.8,
        cols = color_cluster,
        label = TRUE)+xlim(-6.5,7) +ylim(8,17) 

ggsave("02_figures/sub_fib_16sam.pdf",device = "pdf",width = 12, height = 10)

qs_save(sub_fib_16sam, "01_data/sub_fib_16sam.qs2")




#### 2. Fibroblast proprotion change 

library(dplyr)
library(ggplot2)
library(stringr)   # for str_extract
library(scales)    # for percent_format

meta <- sub_fib_16sam@meta.data

# step1  Define a 4-subtype column from predicted.celltype
#    (Here I assume predicted.celltype already is "Fibroblast-1"..."Fibroblast-4")
meta <- meta %>%
  mutate(
    fb_subtype4 = predicted.celltype
  )

# 2) Optional: extract numeric age from "age23a" -> 23
meta <- meta %>%
  mutate(
    donor_age = as.numeric(str_extract(orig.ident, "\\d+"))
  )

# Quick check
table(meta$fb_subtype4)
head(meta[, c("orig.ident", "donor_age", "fb_subtype4")])


# step 2
prop_by_sample <- meta %>%
  group_by(orig.ident, fb_subtype4) %>%
  summarise(n = n(), .groups = "drop_last") %>%
  mutate(prop = n / sum(n)) %>%      # proportion within each sample
  ungroup()

prop_by_sample


# step 3

prop_subset <- prop_by_sample %>%
  filter(orig.ident %in% c("age23a", "age52a", "age85a"))

ggplot(prop_subset,
       aes(x = orig.ident, y = prop, fill = fb_subtype4)) +
  geom_col() +
  labs(
    x = "Sample",
    y = "Proportion of fibroblast subtypes",
    fill = "Subtype"
  ) +
  scale_y_continuous(labels = percent_format(accuracy = 1)) +
  theme_classic() +
  theme(
    axis.text.x = element_text(angle = 45, hjust = 1)
  )


ggsave("02_figures/fibro_proportion_batch1_arm.pdf", width = 6, height = 6.5)




