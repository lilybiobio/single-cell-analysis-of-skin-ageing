#!/usr/bin/env Rscript
script_file <- sub("^--file=", "", grep("^--file=", commandArgs(FALSE), value = TRUE)[1])
PACKAGE_ROOT <- dirname(dirname(normalizePath(script_file)))
source(file.path(PACKAGE_ROOT, "common", "analysis_helpers.R"))

# Extracted from Figure1_S1_SSaging_ra.R. Uses the existing saved UMAP.
obj <- load_six_sample_object()
out <- output_directory("supplementary_umap")
if (!"umap" %in% names(obj@reductions)) stop("Input object has no saved UMAP.")
set.seed(1)
colors <- c(age23a = "lightblue", age23b = "#FFCCCC", age52a = "#4DBBD5FF",
            age52b = "#E64B35FF", age85a = "darkblue", age85b = "maroon")
p <- DimPlot(obj, reduction = "umap", group.by = "orig.ident", cols = colors, shuffle = TRUE)
save_plot(p, out, "umap_by_sample", 8, 6)
p <- DimPlot(obj, reduction = "umap", group.by = "celltype", label = TRUE, repel = TRUE)
save_plot(p, out, "umap_by_celltype", 11, 7)
write.csv(cbind(cell = colnames(obj), Embeddings(obj, "umap"), obj[[]][, c("orig.ident", "celltype")]),
          file.path(out, "umap_coordinates_and_metadata.csv"), row.names = FALSE)
write_session(out)
