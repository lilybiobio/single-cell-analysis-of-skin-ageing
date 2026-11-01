suppressPackageStartupMessages({library(Seurat); library(Matrix); library(ggplot2)})
AGING_COLORS <- c(a = "#2DB88F", b = "#B48AC5")
AGING_LABELS <- c(a = "Combined aging", b = "Chronological aging")
BASAL_STEM_SUPPLEMENT <- paste("Basal cell", c(3, 4, 6), sep = "-")

output_directory <- function(section) {
  root <- Sys.getenv("SKIN_OUTPUT_ROOT", file.path(PACKAGE_ROOT, "generated"))
  path <- file.path(root, section)
  dir.create(path, recursive = TRUE, showWarnings = FALSE)
  path
}
load_six_sample_object <- function() {
  path <- Sys.getenv("SKIN_SIX_SAMPLE_RDS")
  if (!nzchar(path) || !file.exists(path) || file.info(path)$size == 0)
    stop("Set SKIN_SIX_SAMPLE_RDS to the nonempty original six-sample Seurat RDS.")
  obj <- readRDS(path)
  obj$celltype <- as.character(Idents(obj))
  ids <- as.character(obj$orig.ident)
  expected <- as.vector(outer(c(23, 52, 85), c("a", "b"), function(a, b) paste0("age", a, b)))
  if (!setequal(unique(ids), expected)) stop("The six-sample entrypoint requires age23/52/85 a/b IDs.")
  obj$site <- sub("^age[0-9]+([ab])$", "\\1", ids)
  obj$age <- as.integer(sub("^age([0-9]+)[ab]$", "\\1", ids))
  obj$age_group <- factor(c("Young", "Middle", "Old")[match(obj$age, c(23, 52, 85))],
                          levels = c("Young", "Middle", "Old"))
  DefaultAssay(obj) <- "RNA"
  obj
}
save_plot <- function(plot, out, name, width = 7, height = 5) {
  for (ext in c("png", "pdf")) ggsave(file.path(out, paste0(name, ".", ext)), plot,
    width = width, height = height, dpi = 300, bg = "white")
}
write_session <- function(out) capture.output(sessionInfo(), file = file.path(out, "sessionInfo.txt"))
broad_cell_type <- function(labels, include_sebocyte_in_hair = FALSE) {
  result <- rep(NA_character_, length(labels))
  patterns <- c(Keratinocyte = "^Basal cell-|^Spinous cell-|^Granular cell$",
    Fibroblast = "^Fibroblast-", Endothelial = "^Endothelial cell-",
    `Hair follicle` = "^Hair follicle-", Pericyte = "^Pericyte-", Immune = "^Immune cell-")
  for (name in names(patterns)) result[grepl(patterns[[name]], labels)] <- name
  for (name in c("Melanocyte", "Langerhans")) result[labels == name] <- name
  result[labels == "Hair follicle (Sebocyte)"] <- if (include_sebocyte_in_hair) "Hair follicle" else "Sebocyte"
  result
}
