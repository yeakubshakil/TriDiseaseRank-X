################################################################################
# PATCH: R-torch indexing hardened — all class/index labels are 1-based.
# TD-CCHI v4.1 — PURE R / torch PORT
# Transcriptome-Driven Cross-Cluster Hub Identification
#
# Faithful Stage-1 port of the supplied Python TD-CCHI core.
# External CRC-specific validation/adapters are intentionally omitted here.
#
# Primary score:
#   CCHS = (C + I + S + R) / 4
# where
#   C = intra-cluster centrality
#   I = inter-cluster participation/connectivity
#   S = bootstrap selection stability
#   R = cross-cohort representation reproducibility
#
# Required packages:
#   torch, data.table, ggplot2
################################################################################

if (!requireNamespace("torch", quietly = TRUE)) {
  stop("Package 'torch' is required. Install with install.packages('torch').")
}
if (!requireNamespace("data.table", quietly = TRUE)) {
  stop("Package 'data.table' is required.")
}
if (!requireNamespace("ggplot2", quietly = TRUE)) {
  stop("Package 'ggplot2' is required.")
}

suppressPackageStartupMessages({
  library(torch)
  library(data.table)
  library(ggplot2)
})

# ==============================================================================
# 1. CONFIGURATION
# ==============================================================================

td_cchi_config <- function(
    seed = 20260812L,
    edge_proj_dim = 24L,
    gat_hidden_per_head = 16L,
    gat_heads = 4L,
    embedding_dim = 32L,
    dropout = 0.25,
    warmup_epochs = 250L,
    cluster_epochs = 200L,
    learning_rate = 1e-3,
    weight_decay = 1e-4,
    align_temperature = 0.20,
    lambda_align = 0.75,
    lambda_cluster = 0.30,
    lambda_gate_sparse = 0.005,
    lambda_recon_sparse = 1e-4,
    patience = 60L,
    min_delta = 1e-5,
    k_min = 2L,
    k_max = 8L,
    kmeans_n_init = 30L,
    n_bootstrap = 100L,
    top_fraction = 0.10,
    top_n_report = 25L
) {
  list(
    seed = as.integer(seed),
    edge_proj_dim = as.integer(edge_proj_dim),
    gat_hidden_per_head = as.integer(gat_hidden_per_head),
    gat_heads = as.integer(gat_heads),
    embedding_dim = as.integer(embedding_dim),
    dropout = as.numeric(dropout),
    warmup_epochs = as.integer(warmup_epochs),
    cluster_epochs = as.integer(cluster_epochs),
    learning_rate = as.numeric(learning_rate),
    weight_decay = as.numeric(weight_decay),
    align_temperature = as.numeric(align_temperature),
    lambda_align = as.numeric(lambda_align),
    lambda_cluster = as.numeric(lambda_cluster),
    lambda_gate_sparse = as.numeric(lambda_gate_sparse),
    lambda_recon_sparse = as.numeric(lambda_recon_sparse),
    patience = as.integer(patience),
    min_delta = as.numeric(min_delta),
    k_min = as.integer(k_min),
    k_max = as.integer(k_max),
    kmeans_n_init = as.integer(kmeans_n_init),
    n_bootstrap = as.integer(n_bootstrap),
    top_fraction = as.numeric(top_fraction),
    top_n_report = as.integer(top_n_report)
  )
}

set_td_cchi_seed <- function(seed) {
  set.seed(seed)
  torch::torch_manual_seed(as.integer(seed))
  invisible(NULL)
}

get_td_cchi_device <- function(requested = c("auto", "cpu", "cuda")) {
  requested <- match.arg(requested)
  if (requested == "cpu") return(torch_device("cpu"))
  if (requested == "cuda") {
    if (!cuda_is_available()) stop("CUDA requested but no CUDA device is available.")
    return(torch_device("cuda"))
  }
  if (cuda_is_available()) torch_device("cuda") else torch_device("cpu")
}

# ==============================================================================
# 2. INPUT / PREPROCESSING
# ==============================================================================

clean_td_gene <- function(x) {
  x <- trimws(as.character(x))
  x[x %in% c("", "NA", "N/A", "---", "--", "NULL", ".")] <- NA_character_
  x <- sub("[[:space:]]*(///|//|;|\\|).*$", "", x)
  trimws(x)
}

read_td_candidates <- function(path) {
  if (!file.exists(path)) stop("Candidate file not found: ", path)
  ext <- tolower(tools::file_ext(path))
  if (ext == "csv") {
    x <- data.table::fread(path, data.table = FALSE)
  } else if (ext %in% c("txt", "tsv")) {
    x <- data.table::fread(path, data.table = FALSE)
  } else if (ext %in% c("xlsx", "xls")) {
    if (!requireNamespace("readxl", quietly = TRUE)) stop("Package 'readxl' is required for Excel candidate files.")
    x <- readxl::read_excel(path)
  } else {
    stop("Candidate file must be CSV/TSV/TXT/XLSX/XLS.")
  }
  g <- if ("Gene" %in% names(x)) x$Gene else x[[1]]
  g <- unique(na.omit(clean_td_gene(g)))
  if (length(g) < 3L) stop("TD-CCHI requires at least 3 candidate genes.")
  g
}

read_td_expression <- function(expression_file, metadata, genes, sample_col = "Samples") {
  if (!file.exists(expression_file)) stop("Expression file not found: ", expression_file)
  x <- data.table::fread(expression_file, data.table = FALSE, check.names = FALSE)
  if (!all(c("probe_id", "gene_symbol") %in% names(x))) {
    stop("Expression file must contain 'probe_id' and 'gene_symbol'.")
  }
  
  samples <- trimws(as.character(metadata[[sample_col]]))
  miss <- setdiff(samples, names(x))
  if (length(miss)) stop("Expression file is missing metadata samples: ", paste(head(miss, 8), collapse = ", "))
  
  sym <- clean_td_gene(x$gene_symbol)
  keep <- !is.na(sym) & sym %in% genes
  if (!any(keep)) stop("None of the candidate genes are present in the expression file.")
  
  dat <- as.data.frame(x[keep, samples, drop = FALSE], check.names = FALSE)
  dat[] <- lapply(dat, function(z) suppressWarnings(as.numeric(z)))
  mat <- as.matrix(dat)
  finite <- apply(mat, 1, function(z) all(is.finite(z)))
  mat <- mat[finite, , drop = FALSE]
  sym <- sym[keep][finite]
  
  # Collapse duplicate symbols by arithmetic mean, matching the Python implementation.
  split_idx <- split(seq_along(sym), sym)
  out <- t(vapply(split_idx, function(ii) {
    if (length(ii) == 1L) mat[ii, ] else colMeans(mat[ii, , drop = FALSE])
  }, numeric(ncol(mat))))
  
  out <- out[genes[genes %in% rownames(out)], , drop = FALSE]
  missing_genes <- setdiff(genes, rownames(out))
  if (length(missing_genes)) {
    stop(length(missing_genes), " candidate genes are absent after annotation collapse. Examples: ",
         paste(head(missing_genes, 10), collapse = ", "))
  }
  out <- out[genes, , drop = FALSE]
  colnames(out) <- samples
  storage.mode(out) <- "double"
  if (any(!is.finite(out))) stop("Non-finite values remain in candidate expression matrix.")
  out
}

maybe_log2_td <- function(expr) {
  q <- as.numeric(stats::quantile(expr, probs = c(0, .25, .50, .75, .99, 1), na.rm = TRUE, names = FALSE))
  needs <- q[5] > 100 || ((q[6] - q[1] > 50) && q[2] > 0)
  if (!needs) return(list(expr = expr, applied = FALSE))
  if (min(expr, na.rm = TRUE) < 0) stop("Expression appears to require log2 but contains negative values.")
  list(expr = log2(expr + 1), applied = TRUE)
}

read_td_metadata <- function(path, sample_col = "Samples", condition_col = "Condition") {
  if (!file.exists(path)) stop("Metadata file not found: ", path)
  md <- data.table::fread(path, data.table = FALSE, check.names = FALSE)
  names(md) <- trimws(names(md))
  if (!all(c(sample_col, condition_col) %in% names(md))) {
    stop("Metadata must contain '", sample_col, "' and '", condition_col, "'.")
  }
  md[[sample_col]] <- trimws(as.character(md[[sample_col]]))
  md[[condition_col]] <- trimws(as.character(md[[condition_col]]))
  if (anyDuplicated(md[[sample_col]])) stop("Duplicated metadata sample IDs detected.")
  if (!all(md[[condition_col]] %in% c("Control", "Case"))) {
    stop("TD-CCHI R core expects Condition values exactly 'Control' and 'Case'.")
  }
  md
}

has_complete_pairing_td <- function(meta) {
  if (!all(c("Samples", "Condition", "Pair_ID") %in% names(meta))) return(FALSE)
  if (anyNA(meta$Pair_ID)) return(FALSE)
  tab <- table(as.character(meta$Pair_ID), as.character(meta$Condition))
  if (!all(c("Control", "Case") %in% colnames(tab))) return(FALSE)
  nrow(tab) >= 2L && all(tab[, "Control"] == 1L) && all(tab[, "Case"] == 1L) && all(rowSums(tab) == 2L)
}

# ==============================================================================
# 3. TRANSCRIPTOMIC FEATURE REPRESENTATION
# ==============================================================================

row_spearman_td <- function(expr) {
  cmat <- suppressWarnings(stats::cor(t(expr), method = "spearman", use = "pairwise.complete.obs"))
  cmat[!is.finite(cmat)] <- 0
  cmat <- pmax(pmin(cmat, 1), -1)
  diag(cmat) <- 1
  cmat
}

build_td_features <- function(expr) {
  if (!is.matrix(expr)) expr <- as.matrix(expr)
  if (ncol(expr) < 3L) stop("At least 3 biological profiles per cohort are required.")
  
  mn <- rowMeans(expr)
  sdv <- apply(expr, 1, stats::sd)
  med <- apply(expr, 1, stats::median)
  mad_raw <- apply(expr, 1, function(z) stats::median(abs(z - stats::median(z))))
  qs <- t(apply(expr, 1, stats::quantile, probs = c(.10, .25, .75, .90), names = FALSE, type = 7))
  iqr <- qs[, 3] - qs[, 2]
  summaries <- cbind(mn, sdv, mad_raw, iqr, qs[, 1], qs[, 2], med, qs[, 3], qs[, 4])
  
  corr <- row_spearman_td(expr)
  fp <- corr
  diag(fp) <- 0
  feat <- cbind(summaries, fp)
  feat[!is.finite(feat)] <- 0
  storage.mode(feat) <- "double"
  list(features = feat, corr = corr)
}

fit_td_scaler <- function(f1, f2) {
  x <- rbind(f1, f2)
  center <- colMeans(x)
  scale <- apply(x, 2, stats::sd)
  scale[!is.finite(scale) | scale == 0] <- 1
  list(center = center, scale = scale)
}

apply_td_scaler <- function(x, scaler) {
  x <- sweep(x, 2, scaler$center, "-")
  x <- sweep(x, 2, scaler$scale, "/")
  x[!is.finite(x)] <- 0
  x
}

# ==============================================================================
# 4. EDGE GATE + WEIGHTED GAT
# ==============================================================================

EdgeGateR <- nn_module(
  "EdgeGateR",
  initialize = function(feature_dim, projection_dim = 24L) {
    self$project <- nn_sequential(
      nn_linear(feature_dim, projection_dim),
      nn_gelu(),
      nn_layer_norm(projection_dim)
    )
    self$mlp <- nn_sequential(
      nn_linear(2L * projection_dim + 1L, projection_dim),
      nn_gelu(),
      nn_linear(projection_dim, 1L)
    )
  },
  forward = function(x, base_rel) {
    u <- self$project(x)
    ui <- u$unsqueeze(2)
    uj <- u$unsqueeze(1)
    pair <- torch_cat(list(torch_abs(ui - uj), ui * uj, base_rel$unsqueeze(3)), dim = 3)
    gate <- torch_sigmoid(self$mlp(pair)$squeeze(3))
    gate <- 0.5 * (gate + gate$t())
    n <- gate$size(1)
    eye <- torch_eye(n, dtype = torch_float(), device = gate$device)
    gate <- gate * (1 - eye) + eye
    adaptive <- base_rel * gate
    adaptive <- adaptive * (1 - eye) + eye
    list(adaptive = adaptive, gate = gate)
  }
)

WeightedGATLayerR <- nn_module(
  "WeightedGATLayerR",
  initialize = function(in_dim, head_dim, heads, dropout) {
    self$head_dim <- as.integer(head_dim)
    self$heads <- as.integer(heads)
    self$output_dim <- as.integer(head_dim * heads)
    
    self$W <- nn_linear(in_dim, self$output_dim, bias = FALSE)
    self$a_src <- nn_parameter(torch_empty(c(self$heads, self$head_dim)))
    self$a_dst <- nn_parameter(torch_empty(c(self$heads, self$head_dim)))
    self$edge_gamma_raw <- nn_parameter(torch_zeros(self$heads))
    self$skip <- nn_linear(in_dim, self$output_dim, bias = FALSE)
    self$norm <- nn_layer_norm(self$output_dim)
    self$dropout <- nn_dropout(dropout)
    
    nn_init_xavier_uniform_(self$W$weight)
    nn_init_xavier_uniform_(self$skip$weight)
    nn_init_xavier_uniform_(self$a_src)
    nn_init_xavier_uniform_(self$a_dst)
  },
  forward = function(x, edge_weight, return_attention = FALSE) {
    n <- x$size(1)
    h <- self$W(x)$reshape(c(n, self$heads, self$head_dim))
    
    e_src <- torch_einsum("nhd,hd->nh", list(h, self$a_src))
    e_dst <- torch_einsum("nhd,hd->nh", list(h, self$a_dst))
    scores <- e_src$unsqueeze(2) + e_dst$unsqueeze(1)
    scores <- nnf_leaky_relu(scores, negative_slope = 0.2)
    
    gamma <- nnf_softplus(self$edge_gamma_raw)$reshape(c(1, 1, self$heads))
    scores <- scores + gamma * torch_log(edge_weight$clamp_min(1e-8))$unsqueeze(3)
    mask <- edge_weight <= 0
    scores <- scores$masked_fill(mask$unsqueeze(3), -Inf)
    
    attn <- nnf_softmax(scores, dim = 2)
    attn_drop <- self$dropout(attn)
    
    out <- torch_einsum("ijh,jhd->ihd", list(attn_drop, h))$reshape(c(n, self$output_dim))
    out <- self$norm(self$skip(x) + out)
    out <- nnf_gelu(out)
    
    if (return_attention) return(list(out = out, attention = attn))
    out
  }
)

CentroidHolderR <- nn_module(
  "CentroidHolderR",
  initialize = function(k, embedding_dim) {
    self$mu <- nn_parameter(torch_zeros(c(as.integer(k), as.integer(embedding_dim))))
  }
)

TDCCHIModelR <- nn_module(
  "TDCCHIModelR",
  initialize = function(feature_dim, cfg, k_values) {
    self$cfg <- cfg
    self$k_values <- as.integer(k_values)
    self$edge_gate <- EdgeGateR(feature_dim, cfg$edge_proj_dim)
    self$gat1 <- WeightedGATLayerR(feature_dim, cfg$gat_hidden_per_head, cfg$gat_heads, cfg$dropout)
    hidden_dim <- cfg$gat_hidden_per_head * cfg$gat_heads
    self$gat2 <- WeightedGATLayerR(hidden_dim, cfg$embedding_dim, 1L, cfg$dropout)
    self$centroids <- nn_module_list(lapply(self$k_values, function(k) CentroidHolderR(k, cfg$embedding_dim)))
  },
  encode = function(x, base_rel, return_attention = FALSE) {
    eg <- self$edge_gate(x, base_rel)
    adaptive <- eg$adaptive
    gate <- eg$gate
    h <- self$gat1(x, adaptive)
    
    if (return_attention) {
      z <- self$gat2(h, adaptive, return_attention = TRUE)
      return(list(h = z$out, adaptive = adaptive, gate = gate, attention = z$attention$squeeze(3)))
    }
    h <- self$gat2(h, adaptive)
    list(h = h, adaptive = adaptive, gate = gate)
  },
  soft_assign = function(h, k) {
    idx <- match(as.integer(k), self$k_values)
    if (is.na(idx)) stop("Centroids for K=", k, " are unavailable.")
    mu <- self$centroids[[idx]]$mu
    dist2 <- torch_cdist(h, mu, p = 2)$pow(2)
    q <- 1 / (1 + dist2)
    q / q$sum(dim = 2, keepdim = TRUE)$clamp_min(1e-12)
  }
)

# ==============================================================================
# 5. LOSSES / TRAINING
# ==============================================================================

offdiag_mask_td <- function(n, device) {
  !torch_eye(n, dtype = torch_bool(), device = device)
}

graph_reconstruction_loss_td <- function(h, base_rel, lambda_sparse) {
  pred <- torch_sigmoid(torch_mm(h, h$t()) / sqrt(max(1, h$size(2))))
  mask <- offdiag_mask_td(h$size(1), h$device)
  pv <- pred$masked_select(mask)
  rv <- base_rel$masked_select(mask)
  nnf_mse_loss(pv, rv) + lambda_sparse * pv$mean()
}

gene_infonce_td <- function(h1, h2, temperature) {
  a <- nnf_normalize(h1, p = 2, dim = 2)
  b <- nnf_normalize(h2, p = 2, dim = 2)
  logits <- torch_mm(a, b$t()) / temperature
  n <- as.integer(h1$size(1))
  
  # IMPORTANT: R torch CrossEntropy targets MUST be 1-based.
  # Do not use the Python/PyTorch convention 0:(N-1), and avoid torch_arange()
  # here so that endpoint semantics cannot accidentally create an invalid target.
  labels <- torch_tensor(
    as.integer(seq_len(n)),
    dtype = torch_long(),
    device = h1$device
  )
  
  0.5 * (
    nnf_cross_entropy(logits, target = labels) +
      nnf_cross_entropy(logits$t(), target = labels)
  )
}

dec_target_td <- function(q) {
  weight <- q$pow(2) / q$sum(dim = 1, keepdim = TRUE)$clamp_min(1e-12)
  p <- weight / weight$sum(dim = 2, keepdim = TRUE)$clamp_min(1e-12)
  p$detach()
}

kl_batchmean_td <- function(log_q, p) {
  nnf_kl_div(log_q, p, reduction = "sum") / log_q$size(1)
}

multi_resolution_cluster_loss_td <- function(model, h1, h2, k_values) {
  losses <- lapply(k_values, function(k) {
    q1 <- model$soft_assign(h1, k)
    q2 <- model$soft_assign(h2, k)
    p1 <- dec_target_td(q1)
    p2 <- dec_target_td(q2)
    dec1 <- kl_batchmean_td(torch_log(q1$clamp_min(1e-12)), p1)
    dec2 <- kl_batchmean_td(torch_log(q2$clamp_min(1e-12)), p2)
    cons <- 0.5 * (
      kl_batchmean_td(torch_log(q1$clamp_min(1e-12)), q2$detach()) +
        kl_batchmean_td(torch_log(q2$clamp_min(1e-12)), q1$detach())
    )
    dec1 + dec2 + 0.5 * cons
  })
  torch_stack(losses)$mean()
}

initialize_centroids_td <- function(model, h1, h2, k_values, cfg, device) {
  mean_h <- (as_array(h1$detach()$cpu()) + as_array(h2$detach()$cpu())) / 2
  for (j in seq_along(k_values)) {
    k <- k_values[j]
    set.seed(cfg$seed + k)
    km <- stats::kmeans(mean_h, centers = k, nstart = cfg$kmeans_n_init, iter.max = 100)
    centers <- torch_tensor(km$centers, dtype = torch_float(), device = device)
    with_no_grad({ model$centroids[[j]]$mu$copy_(centers) })
  }
  invisible(NULL)
}

clone_state_td <- function(model) {
  lapply(model$state_dict(), function(x) x$detach()$clone())
}

restore_state_td <- function(model, state) {
  model$load_state_dict(state)
  invisible(NULL)
}

train_td_cchi_r_core <- function(x1, r1, x2, r2, cfg, k_values, verbose = TRUE) {
  model <- TDCCHIModelR(x1$size(2), cfg, k_values)
  model$to(device = x1$device)
  history <- list()
  
  base_losses <- function() {
    z1 <- model$encode(x1, r1)
    z2 <- model$encode(x2, r2)
    h1 <- z1$h; h2 <- z2$h
    l_graph <- 0.5 * (
      graph_reconstruction_loss_td(h1, r1, cfg$lambda_recon_sparse) +
        graph_reconstruction_loss_td(h2, r2, cfg$lambda_recon_sparse)
    )
    l_align <- gene_infonce_td(h1, h2, cfg$align_temperature)
    mask <- offdiag_mask_td(z1$gate$size(1), z1$gate$device)
    l_sparse <- 0.5 * (z1$gate$masked_select(mask)$mean() + z2$gate$masked_select(mask)$mean())
    list(h1 = h1, h2 = h2, graph = l_graph, align = l_align, sparse = l_sparse)
  }
  
  # Warmup
  optimizer <- optim_adamw(model$parameters, lr = cfg$learning_rate, weight_decay = cfg$weight_decay)
  best_loss <- Inf; best_state <- clone_state_td(model); bad <- 0L
  model$train()
  
  for (epoch in seq_len(cfg$warmup_epochs)) {
    optimizer$zero_grad()
    L <- base_losses()
    loss <- L$graph + cfg$lambda_align * L$align + cfg$lambda_gate_sparse * L$sparse
    loss$backward()
    nn_utils_clip_grad_norm_(model$parameters, max_norm = 5)
    optimizer$step()
    lv <- loss$item()
    history[[length(history) + 1L]] <- data.frame(
      phase = "warmup", epoch = epoch, total = lv,
      graph = L$graph$item(), alignment = L$align$item(), cluster = 0,
      gate_sparse = L$sparse$item()
    )
    if (lv < best_loss - cfg$min_delta) {
      best_loss <- lv; best_state <- clone_state_td(model); bad <- 0L
    } else bad <- bad + 1L
    if (verbose && (epoch == 1L || epoch %% 50L == 0L || epoch == cfg$warmup_epochs || bad >= cfg$patience)) {
      cat(sprintf("  warmup %4d/%d loss=%.4f graph=%.4f align=%.4f patience=%d\n",
                  epoch, cfg$warmup_epochs, lv, L$graph$item(), L$align$item(), bad))
    }
    if (bad >= cfg$patience) break
  }
  restore_state_td(model, best_state)
  
  model$eval()
  with_no_grad({
    z1 <- model$encode(x1, r1)
    z2 <- model$encode(x2, r2)
    initialize_centroids_td(model, z1$h, z2$h, k_values, cfg, x1$device)
  })
  
  # Cluster phase
  optimizer <- optim_adamw(model$parameters, lr = cfg$learning_rate * 0.5, weight_decay = cfg$weight_decay)
  best_loss <- Inf; best_state <- clone_state_td(model); bad <- 0L
  model$train()
  
  for (epoch in seq_len(cfg$cluster_epochs)) {
    optimizer$zero_grad()
    L <- base_losses()
    l_cluster <- multi_resolution_cluster_loss_td(model, L$h1, L$h2, k_values)
    loss <- L$graph + cfg$lambda_align * L$align + cfg$lambda_cluster * l_cluster + cfg$lambda_gate_sparse * L$sparse
    loss$backward()
    nn_utils_clip_grad_norm_(model$parameters, max_norm = 5)
    optimizer$step()
    lv <- loss$item()
    history[[length(history) + 1L]] <- data.frame(
      phase = "cluster", epoch = epoch, total = lv,
      graph = L$graph$item(), alignment = L$align$item(), cluster = l_cluster$item(),
      gate_sparse = L$sparse$item()
    )
    if (lv < best_loss - cfg$min_delta) {
      best_loss <- lv; best_state <- clone_state_td(model); bad <- 0L
    } else bad <- bad + 1L
    if (verbose && (epoch == 1L || epoch %% 50L == 0L || epoch == cfg$cluster_epochs || bad >= cfg$patience)) {
      cat(sprintf("  cluster %4d/%d loss=%.4f graph=%.4f align=%.4f cluster=%.4f patience=%d\n",
                  epoch, cfg$cluster_epochs, lv, L$graph$item(), L$align$item(), l_cluster$item(), bad))
    }
    if (bad >= cfg$patience) break
  }
  restore_state_td(model, best_state)
  model$eval()
  list(model = model, history = do.call(rbind, history))
}

# ==============================================================================
# 6. C / I / R COMPONENTS
# ==============================================================================

bounded01_td <- function(x) pmax(0, pmin(1, ifelse(is.finite(x), x, 0)))

cluster_assignments_td <- function(model, h, k_values) {
  out <- vector("list", length(k_values)); names(out) <- as.character(k_values)
  with_no_grad({
    for (j in seq_along(k_values)) {
      q <- as_array(model$soft_assign(h, k_values[j])$cpu())
      out[[j]] <- max.col(q, ties.method = "first")
    }
  })
  out
}

compute_CI_single_td <- function(adjacency, assignments) {
  n <- nrow(adjacency)
  A <- adjacency
  diag(A) <- 0
  C <- numeric(n); I <- numeric(n)
  
  for (labels in assignments) {
    clusters <- sort(unique(labels))
    total_strength <- rowSums(A) + 1e-12
    
    for (cl in clusters) {
      members <- which(labels == cl)
      if (!length(members)) next
      within <- rowSums(A[members, members, drop = FALSE])
      mx <- max(within) + 1e-12
      C[members] <- C[members] + within / mx
    }
    
    for (i in seq_len(n)) {
      fr <- vapply(clusters, function(cl) sum(A[i, labels == cl]) / total_strength[i], numeric(1))
      I[i] <- I[i] + (1 - sum(fr^2))
    }
  }
  
  den <- max(1L, length(assignments))
  list(C = bounded01_td(C / den), I = bounded01_td(I / den))
}

compute_R_td <- function(h1, h2, model, k_values) {
  cos <- as_array(nnf_cosine_similarity(h1, h2, dim = 2)$detach()$cpu())
  emb <- pmax(0, pmin(1, (cos + 1) / 2))
  cr <- numeric(length(emb))
  
  with_no_grad({
    for (k in k_values) {
      q1 <- as_array(model$soft_assign(h1, k)$cpu())
      q2 <- as_array(model$soft_assign(h2, k)$cpu())
      cr <- cr + (1 - 0.5 * rowSums(abs(q1 - q2)))
    }
  })
  cr <- pmax(0, pmin(1, cr / max(1L, length(k_values))))
  bounded01_td(0.5 * emb + 0.5 * cr)
}

compute_components_td <- function(model, x1, r1, x2, r2, k_values) {
  model$eval()
  z <- with_no_grad({
    list(
      z1 = model$encode(x1, r1, return_attention = TRUE),
      z2 = model$encode(x2, r2, return_attention = TRUE)
    )
  })
  ass1 <- cluster_assignments_td(model, z$z1$h, k_values)
  ass2 <- cluster_assignments_td(model, z$z2$h, k_values)
  a1 <- as_array(z$z1$adaptive$detach()$cpu())
  a2 <- as_array(z$z2$adaptive$detach()$cpu())
  ci1 <- compute_CI_single_td(a1, ass1)
  ci2 <- compute_CI_single_td(a2, ass2)
  
  C <- bounded01_td(0.5 * (ci1$C + ci2$C))
  I <- bounded01_td(0.5 * (ci1$I + ci2$I))
  R <- compute_R_td(z$z1$h, z$z2$h, model, k_values)
  
  list(
    C = C, I = I, R = R,
    h1 = as_array(z$z1$h$detach()$cpu()),
    h2 = as_array(z$z2$h$detach()$cpu()),
    a1 = a1, a2 = a2,
    gate1 = as_array(z$z1$gate$detach()$cpu()),
    gate2 = as_array(z$z2$gate$detach()$cpu()),
    assignments1 = ass1, assignments2 = ass2
  )
}

# ==============================================================================
# 7. BOOTSTRAP STABILITY
# ==============================================================================

bootstrap_expr_td <- function(expr, meta) {
  if (has_complete_pairing_td(meta)) {
    pairs <- unique(as.character(meta$Pair_ID))
    sampled <- sample(pairs, length(pairs), replace = TRUE)
    idx <- integer(0)
    for (p in sampled) {
      for (cond in c("Control", "Case")) {
        hit <- which(as.character(meta$Pair_ID) == p & as.character(meta$Condition) == cond)
        if (length(hit) != 1L) stop("Invalid 1:1 pair structure during bootstrap.")
        idx <- c(idx, match(as.character(meta$Samples[hit]), colnames(expr)))
      }
    }
    return(expr[, idx, drop = FALSE])
  }
  
  idx <- integer(0)
  for (cond in c("Control", "Case")) {
    s <- as.character(meta$Samples[as.character(meta$Condition) == cond])
    if (length(s) < 2L) stop("Stratified bootstrap requires at least 2 samples per condition.")
    sampled <- sample(s, length(s), replace = TRUE)
    idx <- c(idx, match(sampled, colnames(expr)))
  }
  expr[, idx, drop = FALSE]
}

bootstrap_stability_td <- function(model, expr1, meta1, expr2, meta2, scaler, cfg, k_values, device, genes, verbose = TRUE) {
  n <- nrow(expr1)
  top_k <- max(1L, ceiling(cfg$top_fraction * n))
  counts <- numeric(n)
  log_rows <- vector("list", cfg$n_bootstrap)
  set.seed(cfg$seed + 99173L)
  model$eval()
  
  for (b in seq_len(cfg$n_bootstrap)) {
    e1 <- bootstrap_expr_td(expr1, meta1)
    e2 <- bootstrap_expr_td(expr2, meta2)
    f1 <- build_td_features(e1)
    f2 <- build_td_features(e2)
    x1 <- torch_tensor(apply_td_scaler(f1$features, scaler), dtype = torch_float(), device = device)
    x2 <- torch_tensor(apply_td_scaler(f2$features, scaler), dtype = torch_float(), device = device)
    r1 <- torch_tensor(abs(f1$corr), dtype = torch_float(), device = device)
    r2 <- torch_tensor(abs(f2$corr), dtype = torch_float(), device = device)
    
    comp <- compute_components_td(model, x1, r1, x2, r2, k_values)
    interim <- (comp$C + comp$I + comp$R) / 3
    selected <- order(interim, decreasing = TRUE)[seq_len(top_k)]
    counts[selected] <- counts[selected] + 1
    log_rows[[b]] <- data.frame(
      bootstrap = b,
      top_k = top_k,
      selected_gene_indices = paste(selected, collapse = ";"),
      selected_genes = paste(genes[selected], collapse = ";"),
      stringsAsFactors = FALSE
    )
    if (verbose && (b == 1L || b %% 10L == 0L || b == cfg$n_bootstrap)) {
      cat(sprintf("  bootstrap %4d/%d\n", b, cfg$n_bootstrap))
    }
  }
  list(S = bounded01_td(counts / cfg$n_bootstrap), log = do.call(rbind, log_rows))
}

# ==============================================================================
# 8. ENTROPY SENSITIVITY
# ==============================================================================

entropy_weights_td <- function(Y) {
  Y <- as.matrix(Y); n <- nrow(Y); p <- ncol(Y)
  if (n < 2L) return(rep(1 / p, p))
  div <- numeric(p)
  for (j in seq_len(p)) {
    z <- pmax(Y[, j], 0)
    if (diff(range(z)) < 1e-12 || sum(z) <= 1e-12) next
    pr <- z / sum(z)
    ent <- -sum(pr * log(pr + 1e-12)) / log(n)
    div[j] <- max(0, 1 - ent)
  }
  if (sum(div) <= 1e-12) rep(1 / p, p) else div / sum(div)
}

# ==============================================================================
# 9. PUBLIC RUNNER
# ==============================================================================

run_td_cchi_r <- function(
    cohort1_expression,
    cohort1_metadata,
    cohort2_expression,
    cohort2_metadata,
    candidate_gene_file,
    cohort1_name = "Cohort1",
    cohort2_name = "Cohort2",
    output_dir = "TD_CCHI_R_results",
    bootstraps = 100L,
    top_fraction = 0.10,
    warmup_epochs = 250L,
    cluster_epochs = 200L,
    patience = 60L,
    k_min = 2L,
    k_max = 8L,
    seed = 20260812L,
    device = c("auto", "cpu", "cuda"),
    verbose = TRUE
) {
  device <- match.arg(device)
  cfg <- td_cchi_config(
    seed = seed,
    warmup_epochs = warmup_epochs,
    cluster_epochs = cluster_epochs,
    n_bootstrap = bootstraps,
    top_fraction = top_fraction,
    patience = patience,
    k_min = k_min,
    k_max = k_max
  )
  if (!(cfg$top_fraction > 0 && cfg$top_fraction <= 1)) stop("top_fraction must be in (0,1].")
  set_td_cchi_seed(cfg$seed)
  dev <- get_td_cchi_device(device)
  dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)
  
  genes <- read_td_candidates(candidate_gene_file)
  n_genes <- length(genes)
  k_values <- seq.int(max(2L, cfg$k_min), min(cfg$k_max, n_genes - 1L))
  if (!length(k_values)) stop("Candidate count is too small for the requested K range.")
  
  if (verbose) {
    cat("\n", strrep("=", 72), "\n", sep = "")
    cat("TD-CCHI v4.1 — PURE R/torch\n")
    cat(strrep("=", 72), "\n", sep = "")
    cat("Device:", as.character(dev), "| candidates:", n_genes, "| K:", paste(range(k_values), collapse = ".."), "\n")
  }
  
  meta1 <- read_td_metadata(cohort1_metadata)
  meta2 <- read_td_metadata(cohort2_metadata)
  expr1 <- read_td_expression(cohort1_expression, meta1, genes)
  expr2 <- read_td_expression(cohort2_expression, meta2, genes)
  lg1 <- maybe_log2_td(expr1); expr1 <- lg1$expr
  lg2 <- maybe_log2_td(expr2); expr2 <- lg2$expr
  
  # Preserve metadata order exactly.
  expr1 <- expr1[, meta1$Samples, drop = FALSE]
  expr2 <- expr2[, meta2$Samples, drop = FALSE]
  
  f1 <- build_td_features(expr1)
  f2 <- build_td_features(expr2)
  scaler <- fit_td_scaler(f1$features, f2$features)
  
  x1 <- torch_tensor(apply_td_scaler(f1$features, scaler), dtype = torch_float(), device = dev)
  x2 <- torch_tensor(apply_td_scaler(f2$features, scaler), dtype = torch_float(), device = dev)
  r1 <- torch_tensor(abs(f1$corr), dtype = torch_float(), device = dev)
  r2 <- torch_tensor(abs(f2$corr), dtype = torch_float(), device = dev)
  
  if (verbose) cat("Feature dimension:", x1$size(2), "\nTraining adaptive GAT...\n")
  trained <- train_td_cchi_r_core(x1, r1, x2, r2, cfg, k_values, verbose)
  model <- trained$model
  core <- compute_components_td(model, x1, r1, x2, r2, k_values)
  
  if (verbose) cat("Estimating bootstrap stability...\n")
  boot <- bootstrap_stability_td(model, expr1, meta1, expr2, meta2, scaler, cfg, k_values, dev, genes, verbose)
  
  Y <- cbind(core$C, core$I, boot$S, core$R)
  colnames(Y) <- c("C", "I", "S", "R")
  cchs <- rowMeans(Y)
  ew <- entropy_weights_td(Y)
  cchs_entropy <- as.numeric(Y %*% ew)
  
  results <- data.frame(
    Gene = genes,
    CCHS = cchs,
    CCHS_EntropySensitivity = cchs_entropy,
    Centrality_C = core$C,
    InterCluster_I = core$I,
    Stability_S = boot$S,
    Reproducibility_R = core$R,
    stringsAsFactors = FALSE
  )
  results$Rank <- rank(-results$CCHS, ties.method = "min")
  results$EntropySensitivity_Rank <- rank(-results$CCHS_EntropySensitivity, ties.method = "min")
  results <- results[order(-results$CCHS, -results$Stability_S), , drop = FALSE]
  rownames(results) <- NULL
  
  write.csv(results, file.path(output_dir, "TD_CCHI_ranked_hub_genes.csv"), row.names = FALSE)
  write.csv(trained$history, file.path(output_dir, "TD_CCHI_training_history.csv"), row.names = FALSE)
  write.csv(boot$log, file.path(output_dir, "TD_CCHI_bootstrap_log.csv"), row.names = FALSE)
  write.csv(data.frame(Component = c("C", "I", "S", "R"), Primary_EqualWeight = rep(.25, 4), Entropy_Sensitivity_Weight = ew),
            file.path(output_dir, "TD_CCHI_component_weights.csv"), row.names = FALSE)
  write.csv(expr1, file.path(output_dir, paste0(cohort1_name, "_gene_expression_used.csv")))
  write.csv(expr2, file.path(output_dir, paste0(cohort2_name, "_gene_expression_used.csv")))
  write.csv(f1$corr, file.path(output_dir, paste0(cohort1_name, "_signed_spearman.csv")))
  write.csv(f2$corr, file.path(output_dir, paste0(cohort2_name, "_signed_spearman.csv")))
  write.csv(core$a1, file.path(output_dir, paste0(cohort1_name, "_adaptive_adjacency.csv")))
  write.csv(core$a2, file.path(output_dir, paste0(cohort2_name, "_adaptive_adjacency.csv")))
  write.csv(core$h1, file.path(output_dir, paste0(cohort1_name, "_gene_embeddings.csv")))
  write.csv(core$h2, file.path(output_dir, paste0(cohort2_name, "_gene_embeddings.csv")))
  
  cluster_rows <- list(); z <- 0L
  for (cohort in c(cohort1_name, cohort2_name)) {
    ass <- if (identical(cohort, cohort1_name)) core$assignments1 else core$assignments2
    for (j in seq_along(ass)) {
      z <- z + 1L
      cluster_rows[[z]] <- data.frame(Cohort = cohort, K = as.integer(names(ass)[j]), Gene = genes, Cluster = ass[[j]])
    }
  }
  write.csv(do.call(rbind, cluster_rows), file.path(output_dir, "TD_CCHI_multi_resolution_clusters.csv"), row.names = FALSE)
  
  # Save the complete R torch module rather than a Python state dict.
  torch_save(model, file.path(output_dir, "TD_CCHI_locked_model_R.pt"))
  saveRDS(list(center = scaler$center, scale = scaler$scale, genes = genes, k_values = k_values, config = cfg),
          file.path(output_dir, "TD_CCHI_preprocessing_R.rds"))
  
  # Concise publication figure.
  top <- head(results, min(cfg$top_n_report, nrow(results)))
  p <- ggplot(top, aes(x = reorder(Gene, CCHS), y = CCHS)) +
    geom_col() + coord_flip() + theme_bw(base_size = 11) +
    labs(x = NULL, y = "Cross-Cluster Hub Score (CCHS)", title = paste(cohort1_name, "vs", cohort2_name))
  ggsave(file.path(output_dir, "TD_CCHI_top_ranked_hubs.png"), p, width = 7.2,
         height = max(5, 0.26 * nrow(top)), dpi = 300)
  
  audit <- list(
    implementation = "Pure R torch port",
    cohort1 = cohort1_name,
    cohort2 = cohort2_name,
    candidate_gene_count = n_genes,
    profiles_cohort1 = ncol(expr1),
    profiles_cohort2 = ncol(expr2),
    log2_applied_cohort1 = lg1$applied,
    log2_applied_cohort2 = lg2$applied,
    k_values = k_values,
    device = as.character(dev),
    primary_weights = c(C = .25, I = .25, S = .25, R = .25),
    entropy_sensitivity_weights = setNames(ew, c("C", "I", "S", "R")),
    config = cfg
  )
  saveRDS(audit, file.path(output_dir, "TD_CCHI_analysis_audit_R.rds"))
  
  if (verbose) {
    cat("TD-CCHI R complete. Top hubs:\n")
    print(head(results, cfg$top_n_report), row.names = FALSE)
    cat("Results:", normalizePath(output_dir, mustWork = FALSE), "\n")
  }
  
  invisible(list(
    results = results,
    model = model,
    history = trained$history,
    bootstrap = boot$log,
    scaler = scaler,
    config = cfg,
    audit = audit
  ))
}

################################################################################
# End of TD_CCHI_GNN_R.R
################################################################################
