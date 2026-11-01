#!/usr/bin/env Rscript
script_file <- sub("^--file=", "", grep("^--file=", commandArgs(FALSE), value = TRUE)[1])
PACKAGE_ROOT <- dirname(dirname(normalizePath(script_file)))
source(file.path(PACKAGE_ROOT, "common", "analysis_helpers.R"))

suppressPackageStartupMessages({library(clusterProfiler); library(org.Hs.eg.db); library(enrichplot)})
obj <- load_six_sample_object()
out <- output_directory("supplementary_basal_enrichment")
basal <- subset(obj, cells = colnames(obj)[grepl("^Basal cell-", obj$celltype)])
basal$basal_class <- ifelse(basal$celltype %in% BASAL_STEM_SUPPLEMENT, "Stem", "Other basal")
markers <- FindMarkers(basal, ident.1 = "Stem", ident.2 = "Other basal", group.by = "basal_class")
markers$gene <- rownames(markers)
write.csv(markers, file.path(out, "basal_stem_markers.csv"), row.names = FALSE)
# Historical supplementary definition: positive fold change and nominal p < .05.
genes <- markers$gene[markers$avg_log2FC > 0 & markers$p_val < 0.05]
enrich_bp <- function(genes) enrichGO(gene = unique(genes), OrgDb = org.Hs.eg.db,
  keyType = "SYMBOL", ont = "BP", pAdjustMethod = "BH", qvalueCutoff = 0.05)
result <- enrich_bp(genes)
write.csv(as.data.frame(result), file.path(out, "basal_stem_GO_BP.csv"), row.names = FALSE)
save_plot(enrichplot::dotplot(result, showCategory = 15) + theme_minimal(),
          out, "basal_stem_marker_enrichment", 7, 6)

# Within each site, compare sample identities in the BC3/4/6 population.
site_results <- lapply(c("a", "b"), function(site) {
  sub <- subset(obj, cells = colnames(obj)[obj$site == site & obj$celltype %in% BASAL_STEM_SUPPLEMENT])
  if (!"sample" %in% colnames(sub[[]])) stop("Historical sample metadata is required.")
  sub <- ScaleData(NormalizeData(sub, verbose = FALSE), verbose = FALSE)
  Idents(sub) <- "sample"
  m <- FindAllMarkers(sub, min.pct = 0.05)
  write.csv(m, file.path(out, paste0("basal_markers_site_", site, ".csv")), row.names = FALSE)
  # Use explicit gene symbols rather than possibly uniquified dataframe row names.
  e <- as.data.frame(enrich_bp(m$gene[m$p_val_adj < 0.05]))
  e$site <- rep(site, nrow(e))
  e
})
pathways <- do.call(rbind, site_results)
write.csv(pathways, file.path(out, "basal_pathways_by_site.csv"), row.names = FALSE)
# This clean panel shows the top 20 adjusted-p pathways per site. It is not
# claimed to recover the screenshot's final curated pathway list.
selected <- unique(unlist(lapply(site_results, function(x) head(x$ID[order(x$p.adjust)], 20))))
plot_data <- pathways[pathways$ID %in% selected, ]
p <- ggplot(plot_data, aes(x = -log10(p.adjust), y = reorder(Description, -log10(p.adjust)), color = site)) +
  geom_point(size = 3) + scale_color_manual(values = AGING_COLORS, labels = AGING_LABELS) +
  labs(x = expression(-log[10]("BH-adjusted p")), y = NULL, color = NULL) + theme_minimal()
save_plot(p, out, "basal_pathway_comparison", 9, 10)
write_session(out)
