#!/usr/bin/env Rscript
script_file <- sub("^--file=", "", grep("^--file=", commandArgs(FALSE), value = TRUE)[1])
PACKAGE_ROOT <- dirname(dirname(normalizePath(script_file)))
source(file.path(PACKAGE_ROOT, "common", "analysis_helpers.R"))

suppressPackageStartupMessages({library(dplyr); library(rstatix); library(ggpubr)})
obj <- load_six_sample_object()
out <- output_directory("supplementary_trends_and_pathways")
obj$trend_celltype <- broad_cell_type(obj$celltype, include_sebocyte_in_hair = TRUE)
keep_types <- c("Fibroblast", "Endothelial", "Hair follicle", "Keratinocyte", "Melanocyte")
obj <- subset(obj, cells = colnames(obj)[obj$trend_celltype %in% keep_types])
genes <- c("NOTCH4", "CXCL8", "CEACAM1", "SOCS3", "MEF2A", "PRMT9", "XAF1")
if (!all(genes %in% rownames(obj))) stop("Some configured trend genes are missing.")
expr <- GetAssayData(obj, assay = "RNA", slot = "data")[genes, , drop = FALSE]
# As in Seurat AverageExpression(data), average in linear normalized units.
linear_expression <- expm1(expr)
md <- obj[[]]
groups <- split(seq_len(ncol(obj)), interaction(md$age_group, md$site, md$trend_celltype, drop = TRUE))
trend <- do.call(rbind, lapply(groups, function(i) data.frame(gene = genes,
  expression = Matrix::rowMeans(linear_expression[, i, drop = FALSE]),
  age_group = md$age_group[i[1]], site = md$site[i[1]], celltype = md$trend_celltype[i[1]])))
write.csv(trend, file.path(out, "gene_trend_means.csv"), row.names = FALSE)
for (site in c("a", "b")) {
  p <- ggplot(trend[trend$site == site, ], aes(age_group, expression, group = gene, color = site)) +
    geom_line() + geom_point(size = 1.5) + facet_grid(gene ~ celltype, scales = "free_y") +
    scale_color_manual(values = AGING_COLORS, labels = AGING_LABELS) + theme_minimal() +
    labs(x = NULL, y = "Mean normalized expression", color = NULL)
  save_plot(p, out, paste0("gene_trends_site_", site), 9, 7)
}
set_dir <- file.path(PACKAGE_ROOT, "supplementary_trends_and_pathways", "gene_sets")
sets <- list()
for (path in list.files(set_dir, "\\.tsv$", full.names = TRUE)) {
  tab <- read.delim(path, check.names = FALSE)
  symbols <- tab[tab[[1]] == "GENE_SYMBOLS", 2]
  if (length(symbols) != 1) stop("Missing or ambiguous GENE_SYMBOLS in ", path)
  sets[[names(tab)[2]]] <- strsplit(symbols, ",", fixed = TRUE)[[1]]
}
sets$genage_human <- read.csv(file.path(set_dir, "genage_human.csv"))$symbol
for (name in names(sets)) {
  # Keep historical AddModuleScore defaults and seed; offline matching is explicit.
  retained <- intersect(sets[[name]], rownames(obj))
  write.csv(data.frame(gene = sets[[name]], present = sets[[name]] %in% retained),
            file.path(out, paste0(name, "_gene_membership.csv")), row.names = FALSE)
  if (!length(retained)) stop("No genes found for ", name)
  obj <- AddModuleScore(obj, features = list(retained), name = name, seed = 1, search = FALSE)
  tab <- data.frame(celltype = obj$celltype, sample_id = obj$orig.ident, age = obj$age_group,
                    site = obj$site, score = obj[[paste0(name, "1")]][, 1])
  tests <- tab %>% group_by(age, celltype) %>% wilcox_test(score ~ site) %>% ungroup()
  # One explicit correction family per gene set: all age-by-celltype tests.
  tests$p_adj <- p.adjust(tests$p, "bonferroni")
  write.csv(tests, file.path(out, paste0(name, "_wilcoxon.csv")), row.names = FALSE)
  p <- ggplot(tab, aes(age, score, fill = site)) + geom_violin(scale = "width", position = position_dodge(0.8)) +
    facet_wrap(~celltype, scales = "free_y", ncol = 5) +
    scale_fill_manual(values = AGING_COLORS, labels = AGING_LABELS) + theme_classic(base_size = 9) +
    labs(x = NULL, y = name, fill = NULL)
  save_plot(p, out, paste0(name, "_violin"), 12, 9)
}
write_session(out)
