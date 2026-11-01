#!/usr/bin/env Rscript
script_file <- sub("^--file=", "", grep("^--file=", commandArgs(FALSE), value = TRUE)[1])
PACKAGE_ROOT <- dirname(dirname(normalizePath(script_file)))
source(file.path(PACKAGE_ROOT, "common", "analysis_helpers.R"))

# Historical score = sum of significant positive/negative marker log2FCs.
obj <- load_six_sample_object()
out <- output_directory("perturbation_scores")
results <- list()
for (celltype in sort(unique(obj$celltype))) for (site in c("a", "b")) {
  sub <- subset(obj, cells = colnames(obj)[obj$celltype == celltype & obj$site == site])
  if (!"sample" %in% colnames(sub[[]])) stop("Historical sample metadata is required.")
  sub <- ScaleData(NormalizeData(sub, verbose = FALSE), verbose = FALSE)
  Idents(sub) <- "sample"
  # Fail on errors rather than silently converting failed analyses into zeros.
  m <- FindAllMarkers(sub, min.pct = 0.05)
  name <- gsub("[^A-Za-z0-9]+", "_", paste(celltype, site))
  write.csv(m, file.path(out, paste0(name, "_markers.csv")), row.names = FALSE)
  m <- m[m$p_val_adj < 0.05, ]
  results[[length(results) + 1]] <- data.frame(celltype = celltype, site = site,
    enrichment = sum(m$avg_log2FC[m$avg_log2FC > 0]), depletion = sum(m$avg_log2FC[m$avg_log2FC < 0]))
}
tab <- do.call(rbind, results)
write.csv(tab, file.path(out, "perturbation_scores_all_celltypes.csv"), row.names = FALSE)
long <- rbind(transform(tab, direction = "Enrichment", score = enrichment),
              transform(tab, direction = "Depletion", score = depletion))
p <- ggplot(long, aes(celltype, score, fill = site)) + geom_col(position = "dodge") + coord_flip() +
  scale_fill_manual(values = AGING_COLORS, labels = AGING_LABELS) + theme_classic() +
  labs(x = NULL, y = "Sum of significant marker log2 fold changes", fill = NULL)
save_plot(p, out, "perturbation_scores_all_celltypes", 8, 10)
write_session(out)
