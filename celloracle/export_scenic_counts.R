suppressPackageStartupMessages({
  library(qs2); library(Seurat); library(Matrix)
})

OUT_DIR <- Sys.getenv("SCENIC_OUT", "generated/scenic")
dir.create(OUT_DIR, showWarnings = FALSE, recursive = TRUE)

cat("[", format(Sys.time()), "] loading\n", sep="")
input <- Sys.getenv("SCENIC_INPUT_QS2")
if (!nzchar(input) || !file.exists(input)) stop("Set SCENIC_INPUT_QS2 to the original full-atlas object.")
obj <- qs_read(input)

# drop NA or "NA"-string celltype
ct <- as.character(obj$predicted.celltype)
keep <- !is.na(ct) & ct != "NA" & ct != ""
cat("dropping", sum(!keep), "NA celltype cells\n")
obj <- obj[, keep]
obj$celltype <- obj$predicted.celltype
DefaultAssay(obj) <- "RNA"

oid <- as.character(obj$orig.ident)
stopifnot(all(grepl("^age[0-9]+[ab]$", oid)))
obj$donor_age <- as.integer(sub("^age([0-9]+)[ab]$", "\\1", oid))
obj$site      <- sub("^age[0-9]+([ab])$", "\\1", oid)
obj$age_group <- cut(obj$donor_age, breaks = c(0,35,65,200),
                      labels = c("Y","M","O"), right = FALSE)

cat("full cells:", ncol(obj), "\n")
cts <- GetAssayData(obj, assay="RNA", layer="counts")
gene_ncell <- Matrix::rowSums(cts > 0)
# Retain genes detected in at least 0.5% of retained cells.
keep_gene <- gene_ncell >= ncol(obj) * 0.005
cat("genes kept:", sum(keep_gene), "of", length(keep_gene), "\n")
stopifnot("ATF3" %in% rownames(cts)[keep_gene])
full_cts <- cts[keep_gene, ]

stopifnot(all(full_cts@x == round(full_cts@x)))

# MTX + metadata
writeMM(full_cts, file.path(OUT_DIR, "full_counts.mtx"))
writeLines(rownames(full_cts), file.path(OUT_DIR, "full_genes.tsv"))
writeLines(colnames(full_cts), file.path(OUT_DIR, "full_barcodes.tsv"))
meta_out <- data.frame(
  barcode    = colnames(full_cts),
  celltype   = as.character(obj$celltype),
  age_group  = as.character(obj$age_group),
  site       = as.character(obj$site),
  orig.ident = as.character(obj$orig.ident),
  donor_age  = as.integer(obj$donor_age),
  stringsAsFactors = FALSE
)
write.table(meta_out, file.path(OUT_DIR, "full_meta.tsv"),
            sep = "\t", quote = FALSE, row.names = FALSE)
cat("[", format(Sys.time()), "] done\n", sep="")
