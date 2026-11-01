#!/usr/bin/env Rscript
script_file <- sub("^--file=", "", grep("^--file=", commandArgs(FALSE), value = TRUE)[1])
PACKAGE_ROOT <- dirname(dirname(normalizePath(script_file)))
source(file.path(PACKAGE_ROOT, "common", "analysis_helpers.R"))

# Replot the recovered marker tables without rerunning single-cell inference.
out <- output_directory("supplementary_deg")
input <- file.path(PACKAGE_ROOT, "supplementary_deg", "saved_results", "DEG_each_celltype")
files <- sort(list.files(input, "\\.txt$", full.names = TRUE))
if (!length(files)) stop("No saved marker tables.")
counts <- do.call(rbind, lapply(files, function(path) {
  d <- read.table(path, header = TRUE, check.names = FALSE)
  tokens <- strsplit(sub("\\.txt$", "", basename(path)), "_", fixed = TRUE)[[1]]
  type <- tokens[length(tokens) - 1]; site <- tail(tokens, 1)
  data.frame(source_file = basename(path), celltype = type, site = site,
    up = sum(d$p_val_adj < 0.05 & d$avg_log2FC > 0.25),
    down = sum(d$p_val_adj < 0.05 & d$avg_log2FC < -0.25),
    historical_down_lt_positive_025 = sum(d$p_val_adj < 0.05 & d$avg_log2FC < 0.25))
}))
write.csv(counts, file.path(out, "marker_row_counts.csv"), row.names = FALSE)
long <- rbind(transform(counts, direction = "Up", count = up),
              transform(counts, direction = "Down", count = -down))
p <- ggplot(long, aes(celltype, count, fill = site)) + geom_col(position = "dodge") +
  scale_fill_manual(values = c(arm = AGING_COLORS[["a"]], butt = AGING_COLORS[["b"]])) +
  labs(x = NULL, y = "Significant marker rows", fill = "Site") + theme_classic() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))
save_plot(p, out, "marker_counts", 9, 5)
write_session(out)
