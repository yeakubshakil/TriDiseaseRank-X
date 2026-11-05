############################################################
# TriDiseaseRank-X — FINAL CONCISE DISCOVERY PIPELINE
# Stops at prioritized biomarkers.
#
# Workflow:
#   raw counts + metadata
#   -> limma-voom DEG
#   -> WGCNA
#   -> shared DEG + WGCNA integration
#   -> TD-CCHI custom GNN (pure R/torch)
#   -> 113-model ML
#   -> pathway convergence
#   -> five-layer TriDiseaseRank-X ranking
#
# Deferred: ANN, external validation, regulatory interpretation.
############################################################

# ==========================================================
# 0. PROJECT PATHS + SETTINGS
# ==========================================================

get_script_dir <- function() {
  args <- commandArgs(trailingOnly = FALSE)
  f <- grep("^--file=", args, value = TRUE)
  if (length(f)) return(dirname(normalizePath(sub("^--file=", "", f[1]), mustWork = FALSE)))
  getwd()
}

SCRIPT_DIR <- get_script_dir()
DATA_DIR <- Sys.getenv("TDRX_DATA_DIR", unset = SCRIPT_DIR)
OUT_DIR <- file.path(SCRIPT_DIR, "TriDiseaseRank_results")
dir.create(OUT_DIR, recursive = TRUE, showWarnings = FALSE)

sample_col <- "Samples"
condition_col <- "Condition"
control_label <- "Control"
case_label <- "Case"
candidate_covariates <- c("batch", "sex", "age")

# DEG
DEG_FDR <- 0.05
DEG_LOG2FC <- 1
FILTER_MIN_COUNT <- 10
USE_SVA <- TRUE

# WGCNA
WGCNA_TOP_VAR <- 5000L
WGCNA_NETWORK_TYPE <- "signed"
WGCNA_COR <- "pearson"       # change to "bicor" only if prespecified
WGCNA_MIN_MODULE <- 30L
WGCNA_DEEP_SPLIT <- 2L
WGCNA_MERGE_HEIGHT <- 0.25
WGCNA_MODULE_P <- 0.05
WGCNA_SUPPORT_MIN <- 3L       # strict support in all 3 datasets

# TD-CCHI
RUN_GNN <- TRUE
GNN_BOOTSTRAPS <- 100L
GNN_TOP_FRACTION <- 0.10
GNN_TOP_N_FOR_ML <- 20L
GNN_WARMUP_EPOCHS <- 250L
GNN_CLUSTER_EPOCHS <- 200L
GNN_PATIENCE <- 60L

# 113-ML
RUN_ML <- TRUE
ML_OUTER_FOLDS <- 5L
ML_INNER_FOLDS <- 4L
ML_TUNING_LENGTH <- 3L
ML_TOP_MODELS <- 20L
ML_WORKERS <- min(4L, max(1L, parallel::detectCores(logical = TRUE) - 1L))

# Pathway
RUN_PATHWAY <- TRUE

# Final fusion
FINAL_TOP_N <- 10L
PRIMARY_WEIGHTS <- c(DEG = 0.20, WGCNA = 0.20, GNN = 0.20, ML = 0.20, Pathway = 0.20)

analysis_list <- list(
  D1 = list(gse_id = "GSE172381", counts = "GSE172381_counts_data.csv", metadata = "GSE172381_metadata.csv"),
  D2 = list(gse_id = "GSE120852", counts = "GSE120852_counts_data.csv", metadata = "GSE120852_metadata.csv"),
  D3 = list(gse_id = "GSE161349", counts = "GSE161349_counts_data.csv", metadata = "GSE161349_metadata.csv")
)

# ==========================================================
# 1. PACKAGES + INPUT PREFLIGHT
# ==========================================================

required <- c("edgeR", "limma", "sva", "WGCNA", "data.table", "dplyr", "tidyr",
              "ggplot2", "matrixStats", "openxlsx", "torch")
missing <- required[!vapply(required, requireNamespace, logical(1), quietly = TRUE)]
if (length(missing)) stop("Install missing R packages: ", paste(missing, collapse = ", "))

suppressPackageStartupMessages({
  library(edgeR); library(limma); library(sva); library(WGCNA)
  library(data.table); library(dplyr); library(tidyr); library(ggplot2)
  library(matrixStats); library(openxlsx); library(torch)
})

options(stringsAsFactors = FALSE, scipen = 999)
set.seed(20260908)
try(WGCNA::allowWGCNAThreads(), silent = TRUE)

resolve_input <- function(name) {
  candidates <- unique(c(file.path(DATA_DIR, name), file.path(SCRIPT_DIR, name), file.path(getwd(), name), name))
  hit <- candidates[file.exists(candidates)]
  if (!length(hit)) stop("Input file not found: ", name, "\nChecked: ", paste(candidates, collapse = " | "))
  normalizePath(hit[1], mustWork = TRUE)
}

for (d in names(analysis_list)) {
  analysis_list[[d]]$counts_file <- resolve_input(analysis_list[[d]]$counts)
  analysis_list[[d]]$metadata_file <- resolve_input(analysis_list[[d]]$metadata)

  md <- data.table::fread(analysis_list[[d]]$metadata_file, data.table = FALSE)
  ct <- data.table::fread(analysis_list[[d]]$counts_file, nrows = 2L, data.table = FALSE, check.names = FALSE)

  if (!all(c(sample_col, condition_col) %in% names(md))) {
    stop(d, ": metadata must contain '", sample_col, "' and '", condition_col,
         "'. Found: ", paste(names(md), collapse = ", "))
  }
  if (ncol(ct) < 4L) stop(d, ": count file needs gene ID, gene_symbol and sample columns.")

  count_samples <- names(ct)[-(1:2)]
  meta_samples <- trimws(as.character(md[[sample_col]]))
  missing_samples <- setdiff(meta_samples, count_samples)
  if (length(missing_samples)) stop(d, ": metadata/count sample mismatch: ", paste(head(missing_samples, 10), collapse = ", "))

  bad_labels <- setdiff(unique(trimws(as.character(md[[condition_col]]))), c(control_label, case_label))
  if (length(bad_labels)) stop(d, ": unexpected Condition label(s): ", paste(bad_labels, collapse = ", "))

  cat(sprintf("[OK] %s (%s): %d metadata samples\n", d, analysis_list[[d]]$gse_id, nrow(md)))
}

# ==========================================================
# 2. HELPERS
# ==========================================================

clean_gene <- function(x) {
  x <- toupper(trimws(as.character(x)))
  x[x %in% c("", "NA", "N/A", "---", "--", ".", "NULL")] <- NA_character_
  x
}

minmax01 <- function(x) {
  x <- as.numeric(x); x[!is.finite(x)] <- NA_real_
  r <- range(x, na.rm = TRUE)
  if (!all(is.finite(r)) || r[1] == r[2]) return(rep(0, length(x)))
  y <- (x - r[1]) / (r[2] - r[1]); y[!is.finite(y)] <- 0
  pmin(pmax(y, 0), 1)
}

clamp01 <- function(x) pmin(pmax(as.numeric(x), 0), 1)

usable_covariates <- function(md) {
  x <- intersect(candidate_covariates, names(md))
  x[vapply(x, function(v) length(unique(na.omit(md[[v]]))) > 1L, logical(1))]
}

prepare_meta <- function(md) {
  names(md) <- trimws(names(md))
  md[[sample_col]] <- trimws(as.character(md[[sample_col]]))
  md[[condition_col]] <- trimws(as.character(md[[condition_col]]))
  md <- md[md[[condition_col]] %in% c(control_label, case_label), , drop = FALSE]
  md[[condition_col]] <- factor(md[[condition_col]], levels = c(control_label, case_label))

  covs <- usable_covariates(md)
  for (v in covs) {
    if (is.numeric(md[[v]]) || is.integer(md[[v]])) {
      if (anyNA(md[[v]])) md[[v]][is.na(md[[v]])] <- median(md[[v]], na.rm = TRUE)
    } else {
      z <- as.character(md[[v]]); z[is.na(z) | !nzchar(trimws(z))] <- "Unknown"; md[[v]] <- factor(z)
    }
  }
  list(metadata = md, covariates = covs)
}

save_pca <- function(expr, md, path, title) {
  keep <- matrixStats::rowVars(as.matrix(expr), na.rm = TRUE) > 0
  if (sum(keep) < 2L || ncol(expr) < 3L) return(invisible(NULL))
  pca <- prcomp(t(expr[keep, , drop = FALSE]), scale. = TRUE)
  ve <- 100 * pca$sdev^2 / sum(pca$sdev^2)
  df <- data.frame(PC1 = pca$x[,1], PC2 = pca$x[,2], Condition = md[[condition_col]])
  p <- ggplot(df, aes(PC1, PC2, color = Condition)) + geom_point(size = 3) + theme_bw() +
    labs(title = title, x = sprintf("PC1 (%.1f%%)", ve[1]), y = sprintf("PC2 (%.1f%%)", ve[2]))
  ggsave(path, p, width = 6.5, height = 5, dpi = 300)
}

# ==========================================================
# 3. DEG: RAW COUNTS -> TMM -> VOOM -> LIMMA (+ OPTIONAL SVA)
# ==========================================================

run_deg <- function(d, cfg) {
  cat("\n=== DEG: ", d, " (", cfg$gse_id, ") ===\n", sep = "")
  out <- file.path(OUT_DIR, paste0("Phase1_DEG_", d)); dir.create(out, recursive = TRUE, showWarnings = FALSE)

  raw <- data.table::fread(cfg$counts_file, data.table = FALSE, check.names = FALSE)
  md0 <- data.table::fread(cfg$metadata_file, data.table = FALSE, check.names = FALSE)
  mi <- prepare_meta(md0); md <- mi$metadata; covs <- mi$covariates

  gene_id <- as.character(raw[[1]]); gene_symbol <- clean_gene(raw[[2]])
  counts <- as.matrix(raw[, -(1:2), drop = FALSE]); suppressWarnings(storage.mode(counts) <- "double")
  rownames(counts) <- gene_id

  common <- md[[sample_col]][md[[sample_col]] %in% colnames(counts)]
  if (length(common) != nrow(md)) stop(d, ": not all metadata samples are present in counts.")
  counts <- counts[, common, drop = FALSE]; md <- md[match(common, md[[sample_col]]), , drop = FALSE]

  if (anyNA(counts) || any(!is.finite(counts))) stop(d, ": non-finite count values found.")
  if (any(counts < 0)) stop(d, ": negative values found; this is not raw count data.")
  if (mean(counts %% 1 == 0) < 0.95) warning(d, ": >5% non-integer values; verify raw-count input.")

  # unique gene IDs: retain highest-count row
  keep <- !is.na(rownames(counts)) & nzchar(rownames(counts)) & rowSums(counts) > 0
  counts <- counts[keep, , drop = FALSE]; gene_symbol <- gene_symbol[keep]; gene_id <- gene_id[keep]
  if (anyDuplicated(rownames(counts))) {
    ord <- order(rownames(counts), -rowSums(counts)); counts <- counts[ord, , drop = FALSE]
    gene_symbol <- gene_symbol[ord]; gene_id <- gene_id[ord]
    u <- !duplicated(rownames(counts)); counts <- counts[u, , drop = FALSE]; gene_symbol <- gene_symbol[u]; gene_id <- gene_id[u]
  }
  annot <- data.frame(gene_id = rownames(counts), gene_symbol = gene_symbol, stringsAsFactors = FALSE)

  dge <- edgeR::DGEList(counts)
  design_formula <- reformulate(c(condition_col, covs))
  design <- model.matrix(design_formula, md); colnames(design) <- make.names(colnames(design))

  keep <- edgeR::filterByExpr(dge, group = md[[condition_col]], min.count = FILTER_MIN_COUNT)
  dge <- edgeR::calcNormFactors(dge[keep, , keep.lib.sizes = FALSE], method = "TMM")

  v0 <- limma::voom(dge, design, plot = FALSE)
  final_design <- design; n_sv <- 0L

  if (USE_SVA && ncol(v0$E) - qr(design)$rank >= 3L) {
    null_formula <- if (length(covs)) reformulate(covs) else ~1
    mod0 <- model.matrix(null_formula, md)
    n_sv <- tryCatch(as.integer(sva::num.sv(v0$E, design, method = "be")), error = function(e) 0L)
    if (is.finite(n_sv) && n_sv > 0L) {
      sv <- tryCatch(sva::sva(v0$E, mod = design, mod0 = mod0, n.sv = n_sv)$sv, error = function(e) NULL)
      if (!is.null(sv)) {
        candidate_design <- cbind(design, sv); colnames(candidate_design) <- c(colnames(design), paste0("SV", seq_len(ncol(sv))))
        if (qr(candidate_design)$rank == ncol(candidate_design)) final_design <- candidate_design else n_sv <- 0L
      } else n_sv <- 0L
    }
  }

  v <- limma::voom(dge, final_design, plot = FALSE)
  fit <- limma::eBayes(limma::lmFit(v, final_design), robust = TRUE)

  coef_target <- paste0(make.names(condition_col), make.names(case_label))
  hit <- which(colnames(final_design) == coef_target)
  if (length(hit) != 1L) hit <- grep(paste0(make.names(condition_col), ".*", make.names(case_label)), colnames(final_design))
  if (length(hit) != 1L) stop(d, ": cannot identify Case-vs-Control coefficient. Design: ", paste(colnames(final_design), collapse = ", "))

  tab <- limma::topTable(fit, coef = hit, number = Inf, adjust.method = "BH", sort.by = "P")
  tab$gene_id <- rownames(tab)
  tab$gene_symbol <- annot$gene_symbol[match(tab$gene_id, annot$gene_id)]
  tab$gene_symbol[is.na(tab$gene_symbol)] <- tab$gene_id[is.na(tab$gene_symbol)]
  tab$gene_symbol <- clean_gene(tab$gene_symbol)
  tab$padj <- tab$adj.P.Val
  tab$DEG_raw <- -log10(pmax(tab$padj, 1e-300)) * abs(tab$logFC)
  tab$S_DEG <- minmax01(tab$DEG_raw)
  tab <- tab[!is.na(tab$gene_symbol), , drop = FALSE]
  tab <- tab[order(tab$padj, -abs(tab$logFC)), ]; tab <- tab[!duplicated(tab$gene_symbol), ]
  sig <- tab[tab$padj < DEG_FDR & abs(tab$logFC) > DEG_LOG2FC, , drop = FALSE]

  write.csv(tab, file.path(out, "DEG_all.csv"), row.names = FALSE)
  write.csv(sig, file.path(out, "DEG_significant.csv"), row.names = FALSE)
  save_pca(v$E, md, file.path(out, "PCA_voom.png"), paste(d, cfg$gse_id, "voom PCA"))

  vol <- ggplot(tab, aes(logFC, -log10(pmax(padj, 1e-300)))) +
    geom_point(aes(color = padj < DEG_FDR & abs(logFC) > DEG_LOG2FC), alpha = .65, size = 1.2) +
    geom_vline(xintercept = c(-DEG_LOG2FC, DEG_LOG2FC), linetype = 2) +
    geom_hline(yintercept = -log10(DEG_FDR), linetype = 2) + theme_bw() +
    labs(title = paste(d, cfg$gse_id), x = "log2 fold change", y = "-log10 adjusted P", color = "Significant")
  ggsave(file.path(out, "Volcano.png"), vol, width = 7, height = 5.5, dpi = 300)

  cat("Samples:", ncol(v$E), "| genes:", nrow(v$E), "| significant DEGs:", nrow(sig), "| SVs:", n_sv, "\n")
  list(expr = v$E, metadata = md, deg_all = tab, sig_degs = sig, n_sv = n_sv)
}

deg_results <- lapply(names(analysis_list), function(d) run_deg(d, analysis_list[[d]])); names(deg_results) <- names(analysis_list)

# ==========================================================
# 4. WGCNA ON THE SAME VOOM EXPRESSION
# ==========================================================

run_wgcna <- function(d, deg_obj) {
  cat("\n=== WGCNA: ", d, " ===\n", sep = "")
  out <- file.path(OUT_DIR, paste0("Phase2_WGCNA_", d)); dir.create(out, recursive = TRUE, showWarnings = FALSE)

  expr <- deg_obj$expr; md <- deg_obj$metadata
  vars <- matrixStats::rowVars(expr, na.rm = TRUE); names(vars) <- rownames(expr)
  genes <- names(sort(vars[is.finite(vars) & vars > 0], decreasing = TRUE))
  genes <- head(genes, min(WGCNA_TOP_VAR, length(genes)))
  datExpr <- t(expr[genes, , drop = FALSE])

  gsg <- WGCNA::goodSamplesGenes(datExpr, verbose = 0)
  datExpr <- datExpr[gsg$goodSamples, gsg$goodGenes, drop = FALSE]
  md <- md[match(rownames(datExpr), md[[sample_col]]), , drop = FALSE]
  if (nrow(datExpr) < 15L) warning(d, ": WGCNA sample size is small; interpret modules cautiously.")

  trait <- as.numeric(md[[condition_col]] == case_label)
  powers <- c(1:10, seq(12, 30, 2))
  cor_fnc <- if (WGCNA_COR == "bicor") "bicor" else "cor"
  cor_opts <- if (WGCNA_COR == "bicor") list(use = "p", maxPOutliers = 0.10) else list(use = "p")

  sft <- WGCNA::pickSoftThreshold(datExpr, powerVector = powers, networkType = WGCNA_NETWORK_TYPE,
                                  corFnc = cor_fnc, corOptions = cor_opts, verbose = 0)
  r2 <- -sign(sft$fitIndices[,3]) * sft$fitIndices[,2]
  ok <- which(r2 >= 0.80)
  power <- if (length(ok)) sft$fitIndices[min(ok),1] else sft$fitIndices[which.max(r2),1]

  net <- WGCNA::blockwiseModules(
    datExpr, power = power, maxBlockSize = max(5000L, ncol(datExpr)),
    networkType = WGCNA_NETWORK_TYPE, TOMType = WGCNA_NETWORK_TYPE,
    corType = WGCNA_COR, minModuleSize = WGCNA_MIN_MODULE,
    deepSplit = WGCNA_DEEP_SPLIT, mergeCutHeight = WGCNA_MERGE_HEIGHT,
    reassignThreshold = 0, numericLabels = FALSE, pamRespectsDendro = FALSE,
    verbose = 0
  )

  colors <- net$colors; names(colors) <- colnames(datExpr)
  MEs <- WGCNA::orderMEs(net$MEs)
  mt_cor <- as.numeric(WGCNA::cor(MEs, trait, use = "p"))
  mt_p <- as.numeric(WGCNA::corPvalueStudent(matrix(mt_cor), nrow(datExpr)))
  module_names <- sub("^ME", "", colnames(MEs))
  module_table <- data.frame(Module = module_names, Correlation = mt_cor, P_value = mt_p,
                             FDR = p.adjust(mt_p, "BH"), stringsAsFactors = FALSE)
  module_table <- module_table[module_table$Module != "grey", ]
  sig_modules <- module_table$Module[module_table$P_value < WGCNA_MODULE_P]

  gs <- as.numeric(WGCNA::cor(datExpr, trait, use = "p")); names(gs) <- colnames(datExpr)
  mm_mat <- as.matrix(WGCNA::cor(datExpr, MEs, use = "p")); rownames(mm_mat) <- colnames(datExpr)
  own_mm <- vapply(colnames(datExpr), function(g) {
    col <- paste0("ME", colors[g]); if (col %in% colnames(mm_mat)) mm_mat[g, col] else NA_real_
  }, numeric(1))

  gene_id <- colnames(datExpr)
  symbol_map <- deg_obj$deg_all$gene_symbol[match(gene_id, deg_obj$deg_all$gene_id)]
  symbol_map[is.na(symbol_map)] <- gene_id[is.na(symbol_map)]

  gtab <- data.frame(
    gene_id = gene_id, gene_symbol = clean_gene(symbol_map), Module = unname(colors[gene_id]),
    MM = own_mm, GS = gs[gene_id], stringsAsFactors = FALSE
  )
  gtab$Significant_Module <- gtab$Module %in% sig_modules
  gtab$WGCNA_raw <- abs(gtab$MM) * abs(gtab$GS)
  gtab$S_WGCNA <- minmax01(gtab$WGCNA_raw)
  gtab <- gtab[!is.na(gtab$gene_symbol), ]; gtab <- gtab[order(-gtab$S_WGCNA), ]; gtab <- gtab[!duplicated(gtab$gene_symbol), ]
  sig_genes <- gtab[gtab$Significant_Module & gtab$Module != "grey", , drop = FALSE]

  write.csv(module_table, file.path(out, "module_trait.csv"), row.names = FALSE)
  write.csv(gtab, file.path(out, "WGCNA_all_genes.csv"), row.names = FALSE)
  write.csv(sig_genes, file.path(out, "WGCNA_significant_module_genes.csv"), row.names = FALSE)

  cat("Soft power:", power, "| modules:", nrow(module_table), "| significant modules:", length(sig_modules),
      "| supported genes:", nrow(sig_genes), "\n")
  list(all_genes = gtab, sig_module_genes = sig_genes, modules = module_table, soft_power = power)
}

wgcna_results <- lapply(names(deg_results), function(d) run_wgcna(d, deg_results[[d]])); names(wgcna_results) <- names(deg_results)

# ==========================================================
# 5. SHARED DEG + WGCNA CANDIDATE INTEGRATION
# ==========================================================

cat("\n=== DEG + WGCNA integration ===\n")
datasets <- names(analysis_list)
shared_deg <- Reduce(intersect, lapply(datasets, function(d) unique(clean_gene(deg_results[[d]]$sig_degs$gene_symbol))))
shared_deg <- shared_deg[!is.na(shared_deg)]
if (!length(shared_deg)) stop("No three-way shared significant DEGs under current thresholds.")

cand <- data.frame(Gene = shared_deg, stringsAsFactors = FALSE)
for (d in datasets) {
  de <- deg_results[[d]]$deg_all %>% select(gene_symbol, logFC, padj, S_DEG) %>% distinct(gene_symbol, .keep_all = TRUE)
  names(de) <- c("Gene", paste0("logFC_", d), paste0("padj_", d), paste0("S_DEG_", d))
  wg <- wgcna_results[[d]]$sig_module_genes %>% select(gene_symbol, Module, MM, GS, S_WGCNA) %>% distinct(gene_symbol, .keep_all = TRUE)
  names(wg) <- c("Gene", paste0("Module_", d), paste0("MM_", d), paste0("GS_", d), paste0("S_WGCNA_", d))
  cand <- cand %>% left_join(de, by = "Gene") %>% left_join(wg, by = "Gene")
}

sd_cols <- paste0("S_DEG_", datasets); sw_cols <- paste0("S_WGCNA_", datasets); lfc_cols <- paste0("logFC_", datasets)
cand$S_DEG <- rowMeans(as.matrix(cand[, sd_cols]), na.rm = TRUE)
cand$WGCNA_support_count <- rowSums(!is.na(cand[, sw_cols, drop = FALSE]))
for (x in sw_cols) cand[[x]][is.na(cand[[x]])] <- 0
cand$S_WGCNA <- rowMeans(as.matrix(cand[, sw_cols, drop = FALSE]))
cand$Same_direction <- apply(as.matrix(cand[, lfc_cols]), 1, function(x) all(sign(x) == sign(x[1])))
cand <- cand[cand$WGCNA_support_count >= WGCNA_SUPPORT_MIN, , drop = FALSE]
if (!nrow(cand)) stop("No candidates meet WGCNA_SUPPORT_MIN = ", WGCNA_SUPPORT_MIN, ". Change only as a prespecified sensitivity analysis.")
cand <- cand[order(-rowMeans(cbind(cand$S_DEG, cand$S_WGCNA))), ]

phase3_dir <- file.path(OUT_DIR, "Phase3_DEG_WGCNA_candidates"); dir.create(phase3_dir, recursive = TRUE, showWarnings = FALSE)
write.csv(cand, file.path(phase3_dir, "DEG_WGCNA_candidates.csv"), row.names = FALSE)
cat("Shared DEGs:", length(shared_deg), "| final DEG-WGCNA candidates:", nrow(cand), "\n")
if (nrow(cand) < 3L) stop("TD-CCHI requires at least 3 candidate genes.")

# ==========================================================
# 6. EXPORT NORMALIZED MODEL INPUTS
# ==========================================================

model_dir <- file.path(OUT_DIR, "ModelInputs"); dir.create(model_dir, recursive = TRUE, showWarnings = FALSE)
model_inputs <- setNames(lapply(datasets, function(d) {
  expr <- deg_results[[d]]$expr
  id <- rownames(expr)
  sym <- deg_results[[d]]$deg_all$gene_symbol[match(id, deg_results[[d]]$deg_all$gene_id)]
  sym[is.na(sym)] <- id[is.na(sym)]
  efile <- file.path(model_dir, paste0(d, "_expression.csv"))
  mfile <- file.path(model_dir, paste0(d, "_metadata.csv"))
  write.csv(data.frame(probe_id = id, gene_symbol = clean_gene(sym), expr, check.names = FALSE), efile, row.names = FALSE)
  md <- deg_results[[d]]$metadata
  write.csv(data.frame(Samples = md[[sample_col]], Condition = as.character(md[[condition_col]])), mfile, row.names = FALSE)
  list(expression = efile, metadata = mfile)
}), datasets)

candidate_file <- file.path(model_dir, "TD_CCHI_candidates.csv")
write.csv(data.frame(Gene = cand$Gene), candidate_file, row.names = FALSE)

# ==========================================================
# 7. TD-CCHI CUSTOM GNN — PURE R/torch, PAIRWISE ACROSS DATASETS
# ==========================================================

GNN_SCRIPT <- file.path(SCRIPT_DIR, "TD_CCHI_GNN_R.R")
phase4_dir <- file.path(OUT_DIR, "Phase4_TD_CCHI_GNN")
dir.create(phase4_dir, recursive = TRUE, showWarnings = FALSE)

if (!RUN_GNN) stop("RUN_GNN must be TRUE for the final discovery pipeline.")
if (!file.exists(GNN_SCRIPT)) stop("Pure-R GNN module not found beside master script: ", GNN_SCRIPT)
source(GNN_SCRIPT)

run_gnn_pair <- function(a, b) {
  pair <- paste0(a, "_vs_", b)
  od <- file.path(phase4_dir, pair)
  dir.create(od, recursive = TRUE, showWarnings = FALSE)

  fit <- run_td_cchi_r(
    cohort1_expression = model_inputs[[a]]$expression,
    cohort1_metadata = model_inputs[[a]]$metadata,
    cohort2_expression = model_inputs[[b]]$expression,
    cohort2_metadata = model_inputs[[b]]$metadata,
    candidate_gene_file = candidate_file,
    cohort1_name = a,
    cohort2_name = b,
    output_dir = od,
    bootstraps = GNN_BOOTSTRAPS,
    top_fraction = GNN_TOP_FRACTION,
    warmup_epochs = GNN_WARMUP_EPOCHS,
    cluster_epochs = GNN_CLUSTER_EPOCHS,
    patience = GNN_PATIENCE,
    seed = 20260908L + match(a, datasets) * 1000L + match(b, datasets) * 100L,
    device = "auto",
    verbose = TRUE
  )

  x <- fit$results
  x$Pair <- pair
  x$Pair_Top <- x$Rank <= max(1L, ceiling(nrow(x) * GNN_TOP_FRACTION))
  x
}

pairs <- combn(datasets, 2, simplify = FALSE)
gnn_long <- bind_rows(lapply(pairs, function(p) run_gnn_pair(p[1], p[2])))
write.csv(gnn_long, file.path(phase4_dir, "TD_CCHI_pairwise_long.csv"), row.names = FALSE)

gnn <- gnn_long %>% group_by(Gene) %>% summarise(
  S_GNN = mean(CCHS, na.rm = TRUE), GNN_SD = sd(CCHS, na.rm = TRUE),
  GNN_PairTopSupport = sum(Pair_Top, na.rm = TRUE),
  GNN_C = mean(Centrality_C, na.rm = TRUE), GNN_I = mean(InterCluster_I, na.rm = TRUE),
  GNN_S = mean(Stability_S, na.rm = TRUE), GNN_R = mean(Reproducibility_R, na.rm = TRUE), .groups = "drop"
) %>% arrange(desc(S_GNN), desc(GNN_PairTopSupport))
write.csv(gnn, file.path(phase4_dir, "TD_CCHI_aggregate.csv"), row.names = FALSE)

gnn_hubs <- head(gnn, min(GNN_TOP_N_FOR_ML, nrow(gnn)))
gnn_hub_file <- file.path(phase4_dir, "TD_CCHI_hubs_for_113ML.csv")
write.csv(data.frame(Gene = gnn_hubs$Gene), gnn_hub_file, row.names = FALSE)
cat("TD-CCHI pure-R hubs entering ML:", nrow(gnn_hubs), "\n")

# ==========================================================
# 8. 113-MODEL ML — INTERNAL NESTED CV PER DATASET
# ==========================================================

ML_SCRIPT <- file.path(SCRIPT_DIR, "TriDiseaseRankX_113ML_Internal.R")
phase5_dir <- file.path(OUT_DIR, "Phase5_113ML"); dir.create(phase5_dir, recursive = TRUE, showWarnings = FALSE)
if (!RUN_ML) stop("RUN_ML must be TRUE for final fusion.")
if (!file.exists(ML_SCRIPT)) stop("113-ML module not found beside master script: ", ML_SCRIPT)
source(ML_SCRIPT)

ml_results <- setNames(lapply(seq_along(datasets), function(i) {
  d <- datasets[i]; cat("\n=== 113-ML: ", d, " ===\n", sep = "")
  run_113ml_internal(
    expression_file = model_inputs[[d]]$expression,
    metadata_file = model_inputs[[d]]$metadata,
    candidate_gene_file = gnn_hub_file,
    dataset_name = d,
    output_dir = file.path(phase5_dir, d),
    outer_folds = ML_OUTER_FOLDS, inner_folds = ML_INNER_FOLDS,
    tuning_len = ML_TUNING_LENGTH, workers = ML_WORKERS,
    consensus_top_models = ML_TOP_MODELS, seed = 20260908L + i * 10000L
  )
}), datasets)

ml_long <- bind_rows(lapply(datasets, function(d) transform(ml_results[[d]]$gene_scores, Dataset = d)))
ml <- ml_long %>% select(Gene, Dataset, S_ML_Disease) %>%
  pivot_wider(names_from = Dataset, values_from = S_ML_Disease, names_prefix = "S_ML_")
for (nm in paste0("S_ML_", datasets)) { if (!nm %in% names(ml)) ml[[nm]] <- 0; ml[[nm]][is.na(ml[[nm]])] <- 0 }
ml$S_ML <- rowMeans(as.matrix(ml[, paste0("S_ML_", datasets), drop = FALSE]))
ml$ML_DatasetSupport <- rowSums(as.matrix(ml[, paste0("S_ML_", datasets), drop = FALSE]) > 0)
ml <- ml[order(-ml$S_ML), ]
write.csv(ml_long, file.path(phase5_dir, "113ML_gene_scores_long.csv"), row.names = FALSE)
write.csv(ml, file.path(phase5_dir, "113ML_aggregate.csv"), row.names = FALSE)

# ==========================================================
# 9. PATHWAY ENRICHMENT + CONVERGENCE
# ==========================================================

PATHWAY_SCRIPT <- file.path(SCRIPT_DIR, "TriDiseaseRankX_Pathway.R")
phase6_dir <- file.path(OUT_DIR, "Phase6_Pathway"); dir.create(phase6_dir, recursive = TRUE, showWarnings = FALSE)
if (!RUN_PATHWAY) stop("RUN_PATHWAY must be TRUE for final fusion.")
if (!file.exists(PATHWAY_SCRIPT)) stop("Pathway module not found beside master script: ", PATHWAY_SCRIPT)
source(PATHWAY_SCRIPT)

pathway_gene_lists <- setNames(lapply(datasets, function(d) {
  intersect(unique(clean_gene(deg_results[[d]]$sig_degs$gene_symbol)),
            unique(clean_gene(wgcna_results[[d]]$sig_module_genes$gene_symbol)))
}), datasets)
pathway_universe <- setNames(lapply(datasets, function(d) unique(na.omit(clean_gene(deg_results[[d]]$deg_all$gene_symbol)))), datasets)

pathway_result <- run_tdrx_pathway(
  disease_gene_lists = pathway_gene_lists,
  universe_gene_lists = pathway_universe,
  candidate_genes = cand$Gene,
  output_dir = phase6_dir,
  padj_cutoff = 0.05, pvalue_cutoff = 0.05, qvalue_cutoff = 0.20,
  min_gs_size = 5L, max_gs_size = 500L,
  run_reactome = TRUE, run_wikipathways = TRUE, run_disease_ontology = TRUE,
  top_terms_plot = 20L
)
pathway <- pathway_result$gene_scores

# ==========================================================
# 10. FIVE-LAYER TRIDISEASERANK-X FUSION
# ==========================================================

cat("\n=== Final TriDiseaseRank-X fusion ===\n")
phase7_dir <- file.path(OUT_DIR, "Phase7_Prioritized_Biomarkers"); dir.create(phase7_dir, recursive = TRUE, showWarnings = FALSE)

base <- cand %>% transmute(Gene, S_DEG = clamp01(S_DEG), S_WGCNA = clamp01(S_WGCNA), WGCNA_support_count) %>%
  inner_join(gnn, by = "Gene") %>% inner_join(ml, by = "Gene") %>%
  left_join(pathway %>% select(any_of(c("Gene", "S_Pathway", "PathwayCount", "DatasetPathwaySupport",
                                        "DatabaseSupport", "ConvergentPathwayCount", "SharedAllDatasetPathwayCount"))), by = "Gene")

base$S_GNN <- clamp01(base$S_GNN); base$S_ML <- clamp01(base$S_ML)
base$S_Pathway <- clamp01(ifelse(is.na(base$S_Pathway), 0, base$S_Pathway))

for (layer in names(PRIMARY_WEIGHTS)) {
  score_col <- paste0("S_", layer)
  base[[paste0("Contribution_", layer)]] <- PRIMARY_WEIGHTS[layer] * base[[score_col]]
}
base$TriDiseaseRankX_Score <- rowSums(base[, paste0("Contribution_", names(PRIMARY_WEIGHTS)), drop = FALSE])

rank_pct <- function(x) if (length(x) <= 1L) rep(1, length(x)) else 1 - (rank(-x, ties.method = "average") - 1) / (length(x) - 1)
rmat <- cbind(rank_pct(base$S_DEG), rank_pct(base$S_WGCNA), rank_pct(base$S_GNN), rank_pct(base$S_ML), rank_pct(base$S_Pathway))
base$RankFusion_Sensitivity <- rowMeans(rmat)
cmat <- as.matrix(base[, c("S_DEG", "S_WGCNA", "S_GNN", "S_ML", "S_Pathway")])
labels <- c("DEG", "WGCNA", "TD-CCHI GNN", "113-ML", "Pathway")
base$Strongest_Evidence <- labels[max.col(cmat, ties.method = "first")]
base$Weakest_Evidence <- labels[max.col(-cmat, ties.method = "first")]
base <- base %>% arrange(desc(TriDiseaseRankX_Score), desc(RankFusion_Sensitivity), Gene) %>% mutate(Final_Rank = row_number())

final_panel <- head(base, min(FINAL_TOP_N, nrow(base)))
write.csv(base, file.path(phase7_dir, "TriDiseaseRankX_full_ranking.csv"), row.names = FALSE)
write.csv(final_panel, file.path(phase7_dir, "TriDiseaseRankX_top_biomarkers.csv"), row.names = FALSE)
write.csv(data.frame(Gene = final_panel$Gene), file.path(phase7_dir, "TriDiseaseRankX_locked_panel_for_future_ANN.csv"), row.names = FALSE)
write.csv(data.frame(Layer = names(PRIMARY_WEIGHTS), Weight = as.numeric(PRIMARY_WEIGHTS)),
          file.path(phase7_dir, "TriDiseaseRankX_weights.csv"), row.names = FALSE)

openxlsx::write.xlsx(
  list(Final_Ranking = base, Top_Biomarkers = final_panel, DEG_WGCNA = cand,
       GNN = gnn, ML = ml, Pathway = pathway, Pathway_Convergence = pathway_result$convergence),
  file.path(phase7_dir, "TriDiseaseRankX_results.xlsx"), overwrite = TRUE
)

capture.output(sessionInfo(), file = file.path(phase7_dir, "sessionInfo.txt"))
cat("\n============================================================\n")
cat("TriDiseaseRank-X completed successfully.\n")
cat("Final ranking: ", file.path(phase7_dir, "TriDiseaseRankX_full_ranking.csv"), "\n", sep = "")
cat("Top biomarkers: ", file.path(phase7_dir, "TriDiseaseRankX_top_biomarkers.csv"), "\n", sep = "")
cat("ANN input: ", file.path(phase7_dir, "TriDiseaseRankX_locked_panel_for_future_ANN.csv"), "\n", sep = "")
cat("============================================================\n")
