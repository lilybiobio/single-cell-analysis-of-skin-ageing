#!/usr/bin/env Rscript
script_file <- sub("^--file=", "", grep("^--file=", commandArgs(FALSE), value = TRUE)[1])
PACKAGE_ROOT <- dirname(dirname(normalizePath(script_file)))
source(file.path(PACKAGE_ROOT, "common", "analysis_helpers.R"))

# Recovered six-sample reconstruction, not the original FindAllMarkers tables.
out <- output_directory("pairwise_singlecell_degs")
cell_order <- c("Fibroblast", "Endothelial", "Hair follicle", "Keratinocyte", "Pericyte",
                "Melanocyte", "Immune", "Langerhans", "Sebocyte")
contrast_specs <- list(
  M_vs_Y = c(group1 = "Middle", group2 = "Young"),
  O_vs_Y = c(group1 = "Old", group2 = "Young"),
  O_vs_M = c(group1 = "Old", group2 = "Middle")
)

message("Reading six-sample integrated object...")
obj <- load_six_sample_object()
# Differential expression is performed on normalized RNA expression, not the
# batch-corrected integrated assay (whose negative corrected values can yield
# undefined fold changes in Seurat's data-slot fold-change calculation).
DefaultAssay(obj) <- "RNA"

md <- obj@meta.data
md$cell_name <- rownames(md)
md$broad_cell_type <- broad_cell_type(as.character(Idents(obj)))
md$sample_id <- as.character(md$orig.ident)
md$site_code <- sub("^age[0-9]+([ab])$", "\\1", md$sample_id)
md$age <- as.integer(sub("^age([0-9]+)[ab]$", "\\1", md$sample_id))
md$age_group <- ifelse(md$age == 23, "Young", ifelse(md$age == 52, "Middle", "Old"))

assay <- obj[["RNA"]]
single_cell_sets <- list()
single_cell_tables <- list()
single_cell_summary_parts <- list()
enrichment_parts <- list()
summary_index <- 1L
table_index <- 1L
enrichment_index <- 1L

for (cell_type in cell_order) {
  for (site_code in c("a", "b")) {
    aging_type <- ifelse(site_code == "a", "Combined aging", "Chronological aging")
    cell_site <- md$broad_cell_type == cell_type & md$site_code == site_code

    for (contrast_short in names(contrast_specs)) {
      groups <- contrast_specs[[contrast_short]]
      cells1 <- md$cell_name[cell_site & md$age_group == groups[["group1"]]]
      cells2 <- md$cell_name[cell_site & md$age_group == groups[["group2"]]]
      if (!length(cells1) || !length(cells2)) next

      message(
        cell_type, " / ", aging_type, " / ", contrast_short,
        " (", length(cells1), " vs ", length(cells2), " cells)"
      )
      markers <- FindMarkers(
        assay,
        cells.1 = cells1,
        cells.2 = cells2,
        test.use = "wilcox",
        slot = "data",
        logfc.threshold = 0,
        min.pct = 0,
        only.pos = FALSE,
        verbose = FALSE
      )
      markers$gene <- rownames(markers)
      markers$significant <- markers$p_val_adj < 0.05 & abs(markers$avg_log2FC) > 0.25
      markers$direction <- ifelse(markers$avg_log2FC > 0, "Up", "Down")
      markers$cell_type <- cell_type
      markers$aging_type <- aging_type
      markers$contrast <- contrast_short
      single_cell_tables[[table_index]] <- markers
      table_index <- table_index + 1L

      sig <- markers[markers$significant, , drop = FALSE]
      set_base <- paste(cell_type, aging_type, contrast_short, sep = "||")
      single_cell_sets[[paste0(set_base, "||Up")]] <- sig$gene[sig$avg_log2FC > 0]
      single_cell_sets[[paste0(set_base, "||Down")]] <- sig$gene[sig$avg_log2FC < 0]

      single_cell_summary_parts[[summary_index]] <- data.frame(
        cell_type = cell_type,
        aging_type = aging_type,
        contrast = contrast_short,
        single_cell_DEGs = nrow(sig),
        up_DEGs = sum(sig$avg_log2FC > 0),
        down_DEGs = sum(sig$avg_log2FC < 0),
        stringsAsFactors = FALSE
      )
      summary_index <- summary_index + 1L


    }
  }
}
write.csv(do.call(rbind, single_cell_tables), file.path(out, "pairwise_singlecell_degs.csv"), row.names = FALSE)
write.csv(do.call(rbind, single_cell_summary_parts), file.path(out, "pairwise_singlecell_summary.csv"), row.names = FALSE)
saveRDS(single_cell_sets, file.path(out, "pairwise_gene_sets.rds"))
write_session(out)
