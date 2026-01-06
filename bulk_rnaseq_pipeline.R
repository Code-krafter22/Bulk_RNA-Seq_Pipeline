#!/usr/bin/env Rscript

suppressPackageStartupMessages({
  library(optparse)
  library(edgeR)
  library(limma)
  library(ggplot2)
  library(reshape2)
  library(data.table)
})

# ----------------------------
# Helpers: validation + IO
# ----------------------------
fail <- function(...) stop(sprintf(...), call. = FALSE)

is_zip_magic <- function(path) {
  con <- file(path, "rb")
  on.exit(close(con), add = TRUE)
  sig <- readBin(con, what = "raw", n = 4)
  identical(sig, as.raw(c(0x50, 0x4B, 0x03, 0x04)))  # "PK\003\004"
}

zip_looks_like_xlsx <- function(path) {
  # true when the ZIP contains the XLSX structure
  out <- tryCatch(utils::unzip(path, list = TRUE)$Name, error = function(e) character())
  if (length(out) == 0) return(FALSE)
  any(out == "xl/workbook.xml") || (any(grepl("^xl/", out)) && any(out == "[Content_Types].xml"))
}

# Extracts .zip to temp dir and returns a single extracted file path.
# If the input is actually an XLSX but mislabeled (e.g., .tsv that starts with PK),
# we copy it to a temp .xlsx path so readxl can read it.
materialize_input <- function(path, prefer_pattern = "\\.(csv|tsv|txt)$") {
  if (!file.exists(path)) fail("File not found: %s", path)

  ext <- tolower(tools::file_ext(path))

  # Normal Excel: read directly
  if (ext %in% c("xlsx", "xls")) return(path)

  # If it starts with PK, it's a ZIP container of some kind (zip or xlsx)
  if (ext == "zip" || (is_zip_magic(path) && ext != "gz")) {
    # If it looks like an Excel workbook but has the wrong extension, treat it as xlsx
    if (zip_looks_like_xlsx(path)) {
      tmp_xlsx <- tempfile("excel_", fileext = ".xlsx")
      file.copy(path, tmp_xlsx, overwrite = TRUE)
      return(tmp_xlsx)
    }

    # Otherwise treat it as a normal zip of text files
    td <- tempfile("unz_")
    dir.create(td)
    files <- utils::unzip(path, list = TRUE)$Name
    pick <- files[grepl(prefer_pattern, files, ignore.case = TRUE)]
    if (length(pick) == 0) pick <- files
    if (length(pick) == 0) fail("Zip is empty: %s", path)

    utils::unzip(path, files = pick[1], exdir = td)
    return(file.path(td, pick[1]))
  }

  # non-zip
  path
}

# Reads tabular data from csv/tsv/txt (auto delimiter) and supports .gz,
# PLUS Excel .xlsx/.xls (first sheet by default, configurable).
read_any_table <- function(path, sheet = 1) {
  real_path <- materialize_input(path)
  ext <- tolower(tools::file_ext(real_path))

  if (ext %in% c("xlsx", "xls")) {
    if (!requireNamespace("readxl", quietly = TRUE)) {
      fail(
        "Excel input detected (%s) but package 'readxl' is not installed.\nInstall one of:\n  - conda install -c conda-forge r-readxl\n  - R: install.packages('readxl')",
        real_path
      )
    }
    dt <- readxl::read_excel(real_path, sheet = sheet)
    dt <- as.data.frame(dt, check.names = FALSE)
  } else {
    dt <- data.table::fread(real_path, data.table = FALSE, check.names = FALSE)
  }

  if (nrow(dt) < 2 || ncol(dt) < 2) fail("File has too few rows/cols: %s", path)
  dt
}

# Detect feature column:
# - if user provides feature_col (index or name), use it
# - else pick the first column that is mostly non-numeric while others are numeric-ish
detect_feature_col <- function(df, feature_col = NULL) {
  if (!is.null(feature_col)) {
    if (is.numeric(feature_col)) {
      if (feature_col < 1 || feature_col > ncol(df)) fail("feature-col index out of range.")
      return(feature_col)
    } else {
      if (!feature_col %in% names(df)) fail("feature-col '%s' not in file.", feature_col)
      return(which(names(df) == feature_col)[1])
    }
  }

  # heuristic: choose first column with many non-numeric values
  is_num_like <- function(x) {
    # try numeric conversion; if lots convert, treat as numeric-like
    z <- suppressWarnings(as.numeric(as.character(x)))
    mean(!is.na(z)) > 0.9
  }

  numlike <- vapply(df, is_num_like, logical(1))
  # feature col = first NOT numeric-like
  if (any(!numlike)) return(which(!numlike)[1])

  # fallback: assume first col is feature
  1
}

read_counts_matrix <- function(path, feature_col = NULL, sheet = 1) {
  df <- read_any_table(path, sheet = sheet)
  fc <- detect_feature_col(df, feature_col)

  feature <- as.character(df[[fc]])
  if (any(is.na(feature) | feature == "")) fail("Feature ID column has missing/empty values.")
  feature <- make.unique(feature)

  # remove feature col from numeric matrix
  sample_df <- df[, -fc, drop = FALSE]
  if (ncol(sample_df) < 2) fail("Counts file must contain >=2 sample columns besides the feature column.")

  # coerce to numeric
  suppressWarnings({
    sample_num <- as.data.frame(lapply(sample_df, function(x) as.numeric(as.character(x))))
  })
  bad_cols <- names(sample_num)[colSums(!is.na(sample_num)) == 0]
  if (length(bad_cols) > 0) {
    fail("These sample columns are non-numeric (became all NA): %s", paste(bad_cols, collapse=", "))
  }

  sample_num[is.na(sample_num)] <- 0
  mat <- as.matrix(sample_num)
  rownames(mat) <- feature

  if (any(mat < 0, na.rm = TRUE)) fail("Counts contain negative values. Raw counts must be >= 0.")
  if (anyDuplicated(colnames(mat)) > 0) fail("Duplicate sample column names in counts file.")

  return(mat)
}

read_metadata_any <- function(path, sheet = 1) {
  meta <- read_any_table(path, sheet = sheet)
  meta
}

align_counts_meta <- function(counts_mat, meta, sample_col, group_col) {
  if (!sample_col %in% names(meta)) fail("Metadata missing sample column '%s'.", sample_col)
  if (!group_col %in% names(meta)) fail("Metadata missing group column '%s'.", group_col)

  # sample IDs
  meta[[sample_col]] <- as.character(meta[[sample_col]])
  if (anyDuplicated(meta[[sample_col]]) > 0) fail("Metadata sample IDs have duplicates.")

  # group
  if (any(is.na(meta[[group_col]])) || any(meta[[group_col]] == "")) {
    fail("Group column '%s' contains missing/empty values.", group_col)
  }

  # intersect samples
  samp_counts <- colnames(counts_mat)
  samp_meta   <- meta[[sample_col]]
  common <- intersect(samp_counts, samp_meta)

  if (length(common) < 2) {
    fail("Too few overlapping samples between counts (%d) and metadata (%d). Check sample naming.",
         length(samp_counts), length(samp_meta))
  }

  # warn on mismatches
  missing_in_counts <- setdiff(samp_meta, samp_counts)
  missing_in_meta   <- setdiff(samp_counts, samp_meta)
  if (length(missing_in_counts) > 0) message("WARNING: samples in metadata but not in counts: ", paste(missing_in_counts, collapse=", "))
  if (length(missing_in_meta) > 0) message("WARNING: samples in counts but not in metadata: ", paste(missing_in_meta, collapse=", "))

  meta2 <- meta[match(common, meta[[sample_col]]), , drop = FALSE]
  counts2 <- counts_mat[, common, drop = FALSE]

  # final sanity
  if (!identical(colnames(counts2), meta2[[sample_col]])) fail("Internal alignment error: counts/metadata order mismatch.")
  return(list(counts = counts2, meta = meta2))
}

make_design <- function(meta, group_col, ref_level, covariates = character()) {
  meta[[group_col]] <- factor(meta[[group_col]])
  if (!ref_level %in% levels(meta[[group_col]])) {
    fail("ref-level '%s' not found in %s levels: %s",
         ref_level, group_col, paste(levels(meta[[group_col]]), collapse=", "))
  }
  meta[[group_col]] <- relevel(meta[[group_col]], ref = ref_level)

  # covariate checks
  covariates <- covariates[covariates != ""]
  missing_cov <- setdiff(covariates, names(meta))
  if (length(missing_cov) > 0) fail("Covariates not found in metadata: %s", paste(missing_cov, collapse=", "))

  # if covariate exists but has NAs -> stop (safer)
  for (cv in covariates) {
    if (any(is.na(meta[[cv]]))) fail("Covariate '%s' has missing values. Impute/remove before running.", cv)
  }

  # build formula: ~ Condition + age + sex + batch ...
  rhs <- c(group_col, covariates)
  fml <- as.formula(paste("~", paste(rhs, collapse = " + ")))
  design <- model.matrix(fml, data = meta)
  return(list(design = design, meta = meta, formula = fml))
}

save_session_info <- function(outdir) {
  sink(file.path(outdir, "sessionInfo.txt"))
  print(sessionInfo())
  sink()
}

# ----------------------------
# CLI arguments
# ----------------------------
opt_list <- list(
  make_option(c("--counts"), type="character",
              help="Counts matrix TSV/CSV: rows=genes, cols=samples (first col gene)"),
  make_option(c("--meta"), type="character",
              help="Metadata TSV/CSV/XLSX with sample + group columns"),
  make_option(c("--outdir"), type="character", default="results",
              help="Output directory [default %default]"),
  make_option(c("--sample-col"), type="character", default="Sample",
              help="Metadata sample column [default %default]"),
  make_option(c("--feature-col"), type="character", default=NULL,
              help="Feature ID column (name like 'GeneID' OR index like '1'). If not set, auto-detected."),
  make_option(c("--group-col"), type="character", default="Condition",
              help="Metadata group column [default %default]"),
  make_option(c("--ref-level"), type="character", default="Control",
              help="Reference level for group [default %default]"),
  make_option(c("--covariates"), type="character", default="",
              help="Comma-separated covariates, e.g. age,sex,batch"),
  make_option(c("--filter-logcpm"), type="double", default=1,
              help="Filter threshold on log2CPM [default %default]"),
  make_option(c("--filter-min-samples"), type="integer", default=0,
              help="Min samples above threshold; 0 uses min group size [default %default]"),
  make_option(c("--alpha"), type="double", default=0.05,
              help="P-value cutoff [default %default]"),
  make_option(c("--lfc"), type="double", default=0.58,
              help="abs(logFC) cutoff [default %default]"),
  make_option(c("--counts-sheet"), type="integer", default=1,
              help="Excel sheet index for counts if counts is .xlsx/.xls [default %default]"),
  make_option(c("--meta-sheet"), type="integer", default=1,
              help="Excel sheet index for metadata if meta is .xlsx/.xls [default %default]")
)


opt <- optparse::parse_args(OptionParser(option_list = opt_list))

if (is.null(opt$counts) || is.null(opt$meta)) {
  fail("You must provide --counts and --meta.")
}

dir.create(opt$outdir, recursive = TRUE, showWarnings = FALSE)
dir.create(file.path(opt$outdir, "qc"), showWarnings = FALSE)
dir.create(file.path(opt$outdir, "tables"), showWarnings = FALSE)
dir.create(file.path(opt$outdir, "plots"), showWarnings = FALSE)

# ----------------------------
# Load + validate
# ----------------------------
feature_col <- opt$`feature-col`
if (!is.null(feature_col) && grepl("^[0-9]+$", feature_col)) feature_col <- as.integer(feature_col)

counts_mat <- read_counts_matrix(opt$counts, feature_col = feature_col, sheet = opt$`counts-sheet`)
meta       <- read_metadata_any(opt$meta, sheet = opt$`meta-sheet`)

aligned <- align_counts_meta(counts_mat, meta, opt$`sample-col`, opt$`group-col`)
counts_mat <- aligned$counts
meta <- aligned$meta

covars <- if (opt$covariates == "") character() else trimws(strsplit(opt$covariates, ",")[[1]])

des_out <- make_design(meta, opt$`group-col`, opt$`ref-level`, covariates = covars)
design <- des_out$design
meta <- des_out$meta

# group sanity
tbl <- table(meta[[opt$`group-col`]])
message("Samples per group:")
print(tbl)
if (length(tbl) < 2) fail("Need at least 2 groups in '%s'.", opt$`group-col`)
if (any(tbl < 2)) message("WARNING: One or more groups have <2 samples; results may be unstable.")

# ----------------------------
# QC plots: raw log2(count+1)
# ----------------------------
log_counts <- log2(counts_mat + 1)
log_counts_df <- as.data.frame(log_counts)
log_counts_df$Gene <- rownames(log_counts_df)
melted_counts <- reshape2::melt(log_counts_df, id.vars = "Gene",
                      variable.name = "Sample", value.name = "LogCount")

p_raw <- ggplot(melted_counts, aes(x = Sample, y = LogCount)) +
  geom_boxplot(outlier.size = 0.3) +
  theme_bw() +
  theme(axis.text.x = element_text(angle = 90, hjust = 1, size = 6)) +
  ylab("Log2(Count + 1)") +
  ggtitle("Boxplot: Raw Log-Counts")

ggsave(file.path(opt$outdir, "qc", "boxplot_raw_logcounts.png"), p_raw, width = 10, height = 5, dpi = 300)

# ----------------------------
# edgeR + filter + voom/limma
# ----------------------------
dge <- DGEList(counts = counts_mat, samples = meta)
dge <- edgeR::calcNormFactors(dge, method = "TMM")

log2cpm_raw <- cpm(dge, log = TRUE, prior.count = 1)

min_n <- if (opt$`filter-min-samples` > 0) opt$`filter-min-samples` else min(tbl)
keep <- rowSums(log2cpm_raw > opt$`filter-logcpm`) >= min_n

message(sprintf("Genes before filtering: %d", nrow(dge)))
message(sprintf("Genes kept after filtering: %d (%.2f%%)",
                sum(keep), 100 * sum(keep) / nrow(dge)))

dge_filt <- dge[keep, , keep.lib.sizes = FALSE]
dge_filt <- edgeR::calcNormFactors(dge_filt, method = "TMM")

log2cpm_filtered <- cpm(dge_filt, log = TRUE, prior.count = 1)
write.csv(log2cpm_filtered, file.path(opt$outdir, "tables", "normalized_log2cpm.csv"))

# Boxplot filtered logCPM
log_counts_df_after <- as.data.frame(log2cpm_filtered)
log_counts_df_after$Gene <- rownames(log_counts_df_after)
melted_counts_after <- reshape2::melt(log_counts_df_after, id.vars = "Gene",
                            variable.name = "Sample", value.name = "Log2CPM")

p_filt <- ggplot(melted_counts_after, aes(x = Sample, y = Log2CPM)) +
  geom_boxplot(outlier.size = 0.3) +
  theme_bw() +
  theme(axis.text.x = element_text(angle = 90, hjust = 1, size = 6)) +
  ylab("log2(CPM)") +
  ggtitle("Boxplot: Filtered log2-CPM")

ggsave(file.path(opt$outdir, "qc", "boxplot_filtered_log2cpm.png"), p_filt, width = 10, height = 5, dpi = 300)

# Voom plot
png(file.path(opt$outdir, "plots", "voom_mean_variance_trend.png"), width = 900, height = 600)
v <- voom(dge_filt, design, plot = TRUE)
dev.off()

fit <- lmFit(v, design)
fit <- eBayes(fit)

# ----------------------------
# DE results: by default compare each non-ref level vs ref
# ----------------------------
group_levels <- levels(meta[[opt$`group-col`]])
ref <- opt$`ref-level`
other_levels <- setdiff(group_levels, ref)

all_deg <- list()
all_sig <- list()

for (lvl in other_levels) {
  # coef name depends on model.matrix; for "~ Condition" it is like "ConditionAD"
  coef_name <- paste0(opt$`group-col`, lvl)
  if (!coef_name %in% colnames(design)) {
    # if covariates cause different naming, try to find a matching coefficient
    hit <- grep(paste0("^", opt$`group-col`, lvl, "$"), colnames(design), value = TRUE)
    if (length(hit) == 1) coef_name <- hit
  }
  if (!coef_name %in% colnames(design)) {
    message("Design columns are: ", paste(colnames(design), collapse=", "))
    fail("Could not find coefficient for level '%s'. Expected something like '%s'.", lvl, paste0(opt$`group-col`, lvl))
  }

  deg <- topTable(fit, coef = coef_name, number = Inf, sort.by = "P")
  deg$Gene <- rownames(deg)

  out_all <- file.path(opt$outdir, "tables", sprintf("DEG_%s_vs_%s_ALL.csv", lvl, ref))
  write.csv(deg, out_all, row.names = FALSE)

  deg_sig <- subset(deg, abs(logFC) > opt$lfc & P.Value < opt$alpha)
  out_sig <- file.path(opt$outdir, "tables", sprintf("DEG_%s_vs_%s_SIG.csv", lvl, ref))
  write.csv(deg_sig, out_sig, row.names = FALSE)

  # gene lists
  up <- deg_sig$Gene[deg_sig$logFC > 0]
  dn <- deg_sig$Gene[deg_sig$logFC < 0]
  write.table(up, file.path(opt$outdir, "tables", sprintf("Up_%s_vs_%s.txt", lvl, ref)),
              quote = FALSE, row.names = FALSE, col.names = FALSE)
  write.table(dn, file.path(opt$outdir, "tables", sprintf("Down_%s_vs_%s.txt", lvl, ref)),
              quote = FALSE, row.names = FALSE, col.names = FALSE)

  message(sprintf("[%s vs %s] Significant DEGs: %d", lvl, ref, nrow(deg_sig)))

  all_deg[[paste0(lvl, "_vs_", ref)]] <- deg
  all_sig[[paste0(lvl, "_vs_", ref)]] <- deg_sig
}

# Save metadata used + design
write.csv(meta, file.path(opt$outdir, "tables", "metadata_used.csv"), row.names = FALSE)
write.csv(design, file.path(opt$outdir, "tables", "design_matrix.csv"))

save_session_info(opt$outdir)
message("DONE. Results written to: ", normalizePath(opt$outdir))
