#!/usr/bin/env Rscript

# Donor-level pseudobulk sensitivity analysis for reviewer comment 18.
#
# Study design in the supplied object:
#   age23a, age52a, age85a: arm / sun-exposed (combined aging)
#   age23b, age52b, age85b: back / sun-protected (chronological aging)
#
# There is one donor at each age. Age is therefore confounded with donor.
# The continuous-age models below are exploratory sensitivity analyses and
# must not be described as independently replicated age-group comparisons.

suppressPackageStartupMessages({
  library(Matrix)
  library(edgeR)
  library(ggplot2)
})

object_file <- Sys.getenv("SKIN_SIX_SAMPLE_RDS", "integrated11_seurat_6samples_sub.rds")
output_dir <- Sys.getenv("PSEUDOBULK_OUT", "generated")
cache_dir <- Sys.getenv("PSEUDOBULK_CACHE", "saved_results")
dir.create(output_dir, showWarnings = FALSE, recursive = TRUE)
for (name in c("pseudobulk_counts.rds", "pseudobulk_edgeR_all_results.csv", "cell_counts_by_celltype_sample.csv")) {
  source <- file.path(cache_dir, name)
  target <- file.path(output_dir, name)
  if (file.exists(source) && !file.exists(target)) file.copy(source, target)
}

sample_levels <- c("age23a", "age23b", "age52a", "age52b", "age85a", "age85b")
celltype_levels <- c(
  "Fibroblast", "Endothelial", "Hair follicle", "Keratinocyte", "Pericyte",
  "Melanocyte", "Immune", "Langerhans", "Sebocyte"
)
aging_levels <- c("Combined aging", "Chronological aging")

collapse_cell_type <- function(label) {
  ifelse(
    grepl("^Basal cell-|^Spinous cell-|^Granular cell$", label),
    "Keratinocyte",
    ifelse(
      label == "Hair follicle (Sebocyte)",
      "Sebocyte",
      ifelse(
        grepl("^Hair follicle-", label),
        "Hair follicle",
        ifelse(
          grepl("^Fibroblast-", label),
          "Fibroblast",
          ifelse(
            grepl("^Pericyte-", label),
            "Pericyte",
            ifelse(
              grepl("^Endothelial cell-", label),
              "Endothelial",
              ifelse(
                grepl("^Immune cell-", label),
                "Immune",
                ifelse(
                  label == "Melanocyte",
                  "Melanocyte",
                  ifelse(label == "Langerhans", "Langerhans", NA_character_)
                )
              )
            )
          )
        )
      )
    )
  )
}

counts_file <- file.path(output_dir, "pseudobulk_counts.rds")
cell_counts_file <- file.path(output_dir, "cell_counts_by_celltype_sample.csv")

if (!file.exists(counts_file)) {
  message("Reading Seurat object (large file; this can take about one minute)...")
  obj <- readRDS(object_file)
  labels <- as.character(obj@active.ident)
  broad_cell_type <- collapse_cell_type(labels)
  sample_id <- as.character(obj@meta.data$sample)

  if (anyNA(broad_cell_type)) {
    stop("Unmapped active identities: ", paste(unique(labels[is.na(broad_cell_type)]), collapse = ", "))
  }
  if (!setequal(unique(sample_id), sample_levels)) {
    stop("Unexpected sample IDs: ", paste(sort(unique(sample_id)), collapse = ", "))
  }

  cell_count_matrix <- as.data.frame.matrix(
    table(
      factor(broad_cell_type, levels = celltype_levels),
      factor(sample_id, levels = sample_levels)
    )
  )
  cell_count_matrix$cell_type <- rownames(cell_count_matrix)
  cell_count_matrix <- cell_count_matrix[, c("cell_type", sample_levels)]
  write.csv(cell_count_matrix, cell_counts_file, row.names = FALSE, quote = FALSE)

  group_levels <- as.vector(outer(celltype_levels, sample_levels, paste, sep = "__"))
  group <- factor(paste(broad_cell_type, sample_id, sep = "__"), levels = group_levels)
  membership <- sparseMatrix(
    i = seq_along(group),
    j = as.integer(group),
    x = 1,
    dims = c(length(group), length(group_levels)),
    dimnames = list(colnames(obj@assays$RNA@counts), group_levels)
  )
  message("Aggregating raw RNA counts by cell type and sample...")
  pb_counts <- obj@assays$RNA@counts %*% membership
  colnames(pb_counts) <- group_levels
  saveRDS(pb_counts, counts_file, compress = TRUE)
  rm(obj, membership)
  gc()
} else {
  message("Using existing pseudobulk count matrix: ", counts_file)
  pb_counts <- readRDS(counts_file)
}

run_one_model <- function(cell_type, site_code) {
  columns <- paste0(cell_type, "__age", c(23, 52, 85), site_code)
  y0 <- DGEList(counts = pb_counts[, columns, drop = FALSE])

  # Scaling makes the coefficient an estimated log2 fold change from age 23
  # to age 85 while using the middle-aged donor in the trend fit.
  design_data <- data.frame(age_scaled = c(0, 29 / 62, 1))
  design <- model.matrix(~age_scaled, data = design_data)
  keep <- filterByExpr(y0, design = design, min.count = 10)
  y <- calcNormFactors(y0[keep, , keep.lib.sizes = FALSE], method = "TMM")
  y <- estimateDisp(y, design, robust = TRUE)
  fit <- glmQLFit(y, design, robust = TRUE)

  qlf <- glmQLFTest(fit, coef = 2)
  qlf_table <- topTags(qlf, n = Inf, sort.by = "none")$table

  data.frame(
    gene = rownames(qlf_table),
    cell_type = cell_type,
    site_code = site_code,
    aging_type = ifelse(site_code == "a", "Combined aging", "Chronological aging"),
    logFC_old_vs_young = qlf_table$logFC,
    logCPM = qlf_table$logCPM,
    F = qlf_table$F,
    PValue = qlf_table$PValue,
    FDR = qlf_table$FDR,
    significant = qlf_table$FDR < 0.05 & abs(qlf_table$logFC) > 0.25,
    direction = ifelse(qlf_table$logFC > 0, "Upregulated", "Downregulated"),
    stringsAsFactors = FALSE
  )
}

all_results_file <- file.path(output_dir, "pseudobulk_edgeR_all_results.csv")
if (!file.exists(all_results_file)) {
  result_list <- list()
  index <- 1L
  for (cell_type in celltype_levels) {
    for (site_code in c("a", "b")) {
      message("Fitting ", cell_type, " / site ", site_code, "...")
      result_list[[index]] <- run_one_model(cell_type, site_code)
      index <- index + 1L
    }
  }
  all_results <- do.call(rbind, result_list)
  write.csv(all_results, all_results_file, row.names = FALSE, quote = FALSE)
} else {
  message("Using existing edgeR result table: ", all_results_file)
  all_results <- read.csv(all_results_file, stringsAsFactors = FALSE, check.names = FALSE)
}

base_grid <- expand.grid(
  cell_type = celltype_levels,
  aging_type = aging_levels,
  stringsAsFactors = FALSE
)

count_rows <- function(flag, output_name) {
  selected <- all_results[flag, , drop = FALSE]
  if (nrow(selected) == 0) {
    output <- base_grid
    output[[output_name]] <- 0L
    return(output)
  }
  output <- aggregate(gene ~ cell_type + aging_type, data = selected, FUN = length)
  names(output)[3] <- output_name
  output
}

merge_count <- function(target, source, output_name) {
  output <- merge(target, source, by = c("cell_type", "aging_type"), all.x = TRUE)
  output[[output_name]][is.na(output[[output_name]])] <- 0L
  output
}

summary_table <- base_grid
summary_table <- merge_count(
  summary_table,
  count_rows(all_results$significant, "total_DEGs"),
  "total_DEGs"
)
summary_table <- merge_count(
  summary_table,
  count_rows(all_results$significant & all_results$direction == "Upregulated", "upregulated"),
  "upregulated"
)
summary_table <- merge_count(
  summary_table,
  count_rows(all_results$significant & all_results$direction == "Downregulated", "downregulated"),
  "downregulated"
)
tested <- aggregate(gene ~ cell_type + aging_type, data = all_results, FUN = length)
names(tested)[3] <- "genes_tested"
summary_table <- merge(summary_table, tested, by = c("cell_type", "aging_type"), all.x = TRUE)
summary_table <- summary_table[
  order(
    match(summary_table$cell_type, celltype_levels),
    match(summary_table$aging_type, aging_levels)
  ),
]
write.csv(
  summary_table,
  file.path(output_dir, "pseudobulk_DEG_summary.csv"),
  row.names = FALSE,
  quote = FALSE
)

significant_results <- all_results[all_results$significant, , drop = FALSE]
significant_results <- significant_results[
  order(significant_results$aging_type, significant_results$cell_type, significant_results$FDR),
]
write.csv(
  significant_results,
  file.path(output_dir, "pseudobulk_significant_DEGs.csv"),
  row.names = FALSE,
  quote = FALSE
)

plot_summary <- summary_table
plot_summary$cell_type <- factor(plot_summary$cell_type, levels = celltype_levels)
plot_summary$aging_type <- factor(plot_summary$aging_type, levels = aging_levels)

plot_counts <- ggplot(
  plot_summary,
  aes(x = cell_type, y = total_DEGs, fill = aging_type)
) +
  geom_col(position = position_dodge(width = 0.74), width = 0.66) +
  scale_fill_manual(
    values = c("Combined aging" = "#2DB88F", "Chronological aging" = "#B48AC5")
  ) +
  scale_y_continuous(
    breaks = seq(0, 80, by = 20),
    limits = c(0, 82),
    expand = expansion(mult = c(0, 0))
  ) +
  labs(
    x = NULL,
    y = "Number of DEGs",
    fill = NULL
  ) +
  theme_classic(base_size = 10) +
  theme(
    legend.position = "top",
    legend.justification = "right",
    legend.text = element_text(size = 9),
    legend.key.size = grid::unit(0.62, "lines"),
    legend.spacing.x = grid::unit(0.12, "cm"),
    axis.text.x = element_text(angle = 28, hjust = 1, vjust = 1, size = 8.5, colour = "#333333"),
    axis.text.y = element_text(size = 9, colour = "#333333"),
    axis.title.y = element_text(size = 10, margin = margin(r = 7)),
    axis.line = element_line(linewidth = 0.45),
    axis.ticks = element_line(linewidth = 0.4),
    plot.margin = margin(4, 6, 4, 4)
  )
ggsave(
  file.path(output_dir, "Fig_pseudobulk_DEG_counts.png"),
  plot_counts,
  width = 6.2,
  height = 4.0,
  dpi = 320,
  bg = "white"
)
ggsave(
  file.path(output_dir, "Fig_pseudobulk_DEG_counts.pdf"),
  plot_counts,
  width = 6.2,
  height = 4.0,
  bg = "white"
)

message("Analysis complete. Outputs written to ", output_dir)
print(summary_table)
