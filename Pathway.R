############################################################
# TriDiseaseRank-X pathway/convergence evidence module
# Purpose:
#   - disease/dataset-specific enrichment from DEG + WGCNA-supported genes
#   - GO BP/MF/CC, KEGG; optional Reactome, WikiPathways, Disease Ontology
#   - pathway convergence across datasets
#   - gene-level S_Pathway for final TriDiseaseRank-X rank fusion
#
# Important:
#   Pathway evidence is calculated from dataset-specific DEG-WGCNA gene sets,
#   not from the final ranked biomarker list. This reduces post-selection
#   circularity and permits pathway convergence even when exact genes differ.
############################################################

run_tdrx_pathway <- function(
    disease_gene_lists,
    universe_gene_lists,
    candidate_genes,
    output_dir,
    padj_cutoff = 0.05,
    pvalue_cutoff = 0.05,
    qvalue_cutoff = 0.20,
    min_gs_size = 5L,
    max_gs_size = 500L,
    run_reactome = TRUE,
    run_wikipathways = TRUE,
    run_disease_ontology = TRUE,
    top_terms_plot = 15L
) {
  required <- c(
    "clusterProfiler", "org.Hs.eg.db", "dplyr", "tidyr",
    "data.table", "ggplot2", "openxlsx"
  )
  missing <- required[!vapply(required, requireNamespace, logical(1), quietly = TRUE)]
  if (length(missing)) {
    stop("Missing pathway package(s): ", paste(missing, collapse = ", "))
  }

  dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)
  tables_dir <- file.path(output_dir, "tables")
  plots_dir <- file.path(output_dir, "plots")
  disease_dir <- file.path(output_dir, "by_dataset")
  dir.create(tables_dir, recursive = TRUE, showWarnings = FALSE)
  dir.create(plots_dir, recursive = TRUE, showWarnings = FALSE)
  dir.create(disease_dir, recursive = TRUE, showWarnings = FALSE)

  disease_names <- names(disease_gene_lists)
  if (is.null(disease_names) || length(disease_names) < 2L) {
    stop("disease_gene_lists must be a named list containing at least two datasets.")
  }
  if (!all(disease_names %in% names(universe_gene_lists))) {
    stop("universe_gene_lists must contain the same dataset names as disease_gene_lists.")
  }

  clean_symbol <- function(x) {
    x <- toupper(trimws(as.character(x)))
    x[x %in% c("", "NA", "N/A", "---", "--", ".", "NULL")] <- NA_character_
    x
  }

  minmax01 <- function(x) {
    x <- as.numeric(x)
    x[!is.finite(x)] <- NA_real_
    rng <- range(x, na.rm = TRUE)
    if (!all(is.finite(rng)) || rng[1] == rng[2]) return(rep(0, length(x)))
    y <- (x - rng[1]) / (rng[2] - rng[1])
    y[!is.finite(y)] <- 0
    pmin(pmax(y, 0), 1)
  }

  safe_enrich <- function(expr, label) {
    tryCatch(expr, error = function(e) {
      warning(label, " failed: ", conditionMessage(e))
      NULL
    })
  }

  map_symbols <- function(symbols) {
    symbols <- unique(clean_symbol(symbols))
    symbols <- symbols[!is.na(symbols)]
    if (!length(symbols)) {
      return(data.frame(SYMBOL = character(), ENTREZID = character()))
    }
    suppressMessages({
      x <- clusterProfiler::bitr(
        symbols,
        fromType = "SYMBOL",
        toType = "ENTREZID",
        OrgDb = org.Hs.eg.db::org.Hs.eg.db
      )
    })
    if (is.null(x) || !nrow(x)) {
      return(data.frame(SYMBOL = character(), ENTREZID = character()))
    }
    x$SYMBOL <- clean_symbol(x$SYMBOL)
    x$ENTREZID <- as.character(x$ENTREZID)
    x[!is.na(x$SYMBOL) & !is.na(x$ENTREZID), , drop = FALSE]
  }

  maps <- lapply(disease_names, function(d) map_symbols(disease_gene_lists[[d]]))
  names(maps) <- disease_names
  universe_maps <- lapply(disease_names, function(d) map_symbols(universe_gene_lists[[d]]))
  names(universe_maps) <- disease_names

  global_map <- dplyr::bind_rows(maps, universe_maps) %>%
    dplyr::distinct(ENTREZID, SYMBOL, .keep_all = TRUE)
  entrez_to_symbol <- global_map$SYMBOL
  names(entrez_to_symbol) <- global_map$ENTREZID

  convert_gene_ids <- function(gene_id_string) {
    if (is.na(gene_id_string) || !nzchar(as.character(gene_id_string))) return(character())
    ids <- trimws(unlist(strsplit(as.character(gene_id_string), "/", fixed = TRUE)))
    ids <- ids[nzchar(ids)]
    syms <- unname(entrez_to_symbol[ids])
    syms[is.na(syms)] <- ids[is.na(syms)]
    unique(clean_symbol(syms)[!is.na(clean_symbol(syms))])
  }

  optional_reactome <- run_reactome && requireNamespace("ReactomePA", quietly = TRUE)
  optional_dose <- run_disease_ontology && requireNamespace("DOSE", quietly = TRUE)
  optional_msigdbr <- run_wikipathways && requireNamespace("msigdbr", quietly = TRUE)

  wp_term2gene <- NULL
  if (optional_msigdbr) {
    wp <- tryCatch(
      msigdbr::msigdbr(
        species = "Homo sapiens",
        collection = "C2",
        subcollection = "CP:WIKIPATHWAYS"
      ),
      error = function(e1) {
        tryCatch(
          msigdbr::msigdbr(
            species = "Homo sapiens",
            category = "C2",
            subcategory = "CP:WIKIPATHWAYS"
          ),
          error = function(e2) NULL
        )
      }
    )
    if (!is.null(wp) && nrow(wp)) {
      gene_col <- if ("ncbi_gene" %in% colnames(wp)) "ncbi_gene" else if ("entrez_gene" %in% colnames(wp)) "entrez_gene" else NA_character_
      if (!is.na(gene_col) && "gs_name" %in% colnames(wp)) {
        wp_term2gene <- wp %>%
          dplyr::transmute(term = as.character(gs_name), gene = as.character(.data[[gene_col]])) %>%
          dplyr::filter(!is.na(term), !is.na(gene), nzchar(gene)) %>%
          dplyr::distinct()
      }
    }
  }

  prepare_table <- function(obj, dataset, database) {
    if (is.null(obj)) return(data.frame())
    df <- tryCatch(as.data.frame(obj), error = function(e) data.frame())
    if (!nrow(df)) return(data.frame())
    if (!"p.adjust" %in% colnames(df)) return(data.frame())
    df <- df[is.finite(df$p.adjust) & df$p.adjust < padj_cutoff, , drop = FALSE]
    if (!nrow(df)) return(data.frame())

    if (!"ID" %in% colnames(df)) df$ID <- seq_len(nrow(df))
    if (!"Description" %in% colnames(df)) df$Description <- as.character(df$ID)
    if (!"geneID" %in% colnames(df)) df$geneID <- ""
    if (!"qvalue" %in% colnames(df)) df$qvalue <- NA_real_

    df$Dataset <- dataset
    df$Database <- database
    df$Pathway_ID <- as.character(df$ID)
    df$Pathway_Description <- as.character(df$Description)
    df$Pathway_Key <- paste(database, df$Pathway_ID, sep = "::")
    df$MinusLog10Padj <- -log10(pmax(as.numeric(df$p.adjust), 1e-300))
    df$gene_symbols <- vapply(df$geneID, function(x) paste(convert_gene_ids(x), collapse = "/"), character(1))
    df
  }

  all_tables <- list()
  database_run <- character()

  for (d in disease_names) {
    out_d <- file.path(disease_dir, d)
    dir.create(out_d, recursive = TRUE, showWarnings = FALSE)

    gene_map <- maps[[d]]
    universe_map <- universe_maps[[d]]
    genes <- unique(as.character(gene_map$ENTREZID))
    universe <- unique(as.character(universe_map$ENTREZID))

    data.table::fwrite(
      data.frame(gene_symbol = clean_symbol(disease_gene_lists[[d]])),
      file.path(out_d, paste0(d, "_DEG_WGCNA_genes_for_pathway.csv"))
    )
    data.table::fwrite(gene_map, file.path(out_d, paste0(d, "_symbol_to_entrez.csv")))

    if (length(genes) < min_gs_size) {
      warning(d, ": fewer than ", min_gs_size, " mapped genes; enrichment skipped.")
      next
    }

    common_args_go <- list(
      gene = genes,
      OrgDb = org.Hs.eg.db::org.Hs.eg.db,
      keyType = "ENTREZID",
      pAdjustMethod = "BH",
      pvalueCutoff = pvalue_cutoff,
      qvalueCutoff = qvalue_cutoff,
      minGSSize = min_gs_size,
      maxGSSize = max_gs_size,
      readable = FALSE
    )
    if (length(universe)) common_args_go$universe <- universe

    for (ont in c("BP", "MF", "CC")) {
      args <- common_args_go
      args$ont <- ont
      obj <- safe_enrich(do.call(clusterProfiler::enrichGO, args), paste(d, "GO", ont))
      tab <- prepare_table(obj, d, paste0("GO_", ont))
      all_tables[[paste(d, "GO", ont, sep = "_")]] <- tab
      database_run <- union(database_run, paste0("GO_", ont))
      data.table::fwrite(tab, file.path(out_d, paste0(d, "_GO_", ont, ".csv")))
    }

    kegg_args <- list(
      gene = genes,
      organism = "hsa",
      keyType = "kegg",
      pAdjustMethod = "BH",
      pvalueCutoff = pvalue_cutoff,
      qvalueCutoff = qvalue_cutoff,
      minGSSize = min_gs_size,
      maxGSSize = max_gs_size
    )
    if (length(universe)) kegg_args$universe <- universe
    kegg <- safe_enrich(do.call(clusterProfiler::enrichKEGG, kegg_args), paste(d, "KEGG"))
    kegg_tab <- prepare_table(kegg, d, "KEGG")
    all_tables[[paste0(d, "_KEGG")]] <- kegg_tab
    database_run <- union(database_run, "KEGG")
    data.table::fwrite(kegg_tab, file.path(out_d, paste0(d, "_KEGG.csv")))

    if (optional_reactome) {
      react_args <- list(
        gene = genes,
        organism = "human",
        pAdjustMethod = "BH",
        pvalueCutoff = pvalue_cutoff,
        qvalueCutoff = qvalue_cutoff,
        minGSSize = min_gs_size,
        maxGSSize = max_gs_size,
        readable = FALSE
      )
      if (length(universe)) react_args$universe <- universe
      react <- safe_enrich(do.call(ReactomePA::enrichPathway, react_args), paste(d, "Reactome"))
      react_tab <- prepare_table(react, d, "Reactome")
      all_tables[[paste0(d, "_Reactome")]] <- react_tab
      database_run <- union(database_run, "Reactome")
      data.table::fwrite(react_tab, file.path(out_d, paste0(d, "_Reactome.csv")))
    }

    if (!is.null(wp_term2gene)) {
      wp_args <- list(
        gene = genes,
        TERM2GENE = wp_term2gene,
        pAdjustMethod = "BH",
        pvalueCutoff = pvalue_cutoff,
        qvalueCutoff = qvalue_cutoff,
        minGSSize = min_gs_size,
        maxGSSize = max_gs_size
      )
      if (length(universe)) wp_args$universe <- universe
      wp_obj <- safe_enrich(do.call(clusterProfiler::enricher, wp_args), paste(d, "WikiPathways"))
      wp_tab <- prepare_table(wp_obj, d, "WikiPathways")
      all_tables[[paste0(d, "_WikiPathways")]] <- wp_tab
      database_run <- union(database_run, "WikiPathways")
      data.table::fwrite(wp_tab, file.path(out_d, paste0(d, "_WikiPathways.csv")))
    }

    if (optional_dose) {
      do_obj <- safe_enrich(
        DOSE::enrichDO(
          gene = genes,
          ont = "DO",
          pvalueCutoff = pvalue_cutoff,
          pAdjustMethod = "BH",
          qvalueCutoff = qvalue_cutoff,
          minGSSize = min_gs_size,
          maxGSSize = max_gs_size,
          readable = FALSE
        ),
        paste(d, "DiseaseOntology")
      )
      do_tab <- prepare_table(do_obj, d, "DiseaseOntology")
      all_tables[[paste0(d, "_DiseaseOntology")]] <- do_tab
      database_run <- union(database_run, "DiseaseOntology")
      data.table::fwrite(do_tab, file.path(out_d, paste0(d, "_DiseaseOntology.csv")))
    }
  }

  enrichment <- dplyr::bind_rows(all_tables)
  if (!nrow(enrichment)) {
    warning("No significant pathway enrichment results were produced. S_Pathway will be 0 for all candidates.")
    scores <- data.frame(
      Gene = clean_symbol(candidate_genes),
      PathwayCount = 0L,
      DatasetPathwaySupport = 0L,
      DatabaseSupport = 0L,
      ConvergentPathwayCount = 0L,
      MeanMinusLog10Padj = 0,
      S_PathwayCount = 0,
      S_DatasetPathwaySupport = 0,
      S_DatabaseSupport = 0,
      S_PathwaySignificance = 0,
      S_PathwayConvergence = 0,
      S_Pathway = 0,
      stringsAsFactors = FALSE
    )
    data.table::fwrite(scores, file.path(tables_dir, "TriDiseaseRankX_gene_pathway_scores.csv"))
    return(list(gene_scores = scores, enrichment = enrichment, convergence = data.frame(), edges = data.frame()))
  }

  data.table::fwrite(enrichment, file.path(tables_dir, "TriDiseaseRankX_all_significant_enrichment.csv"))

  presence <- enrichment %>%
    dplyr::distinct(Database, Pathway_ID, Pathway_Key, Pathway_Description, Dataset) %>%
    dplyr::group_by(Database, Pathway_ID, Pathway_Key, Pathway_Description) %>%
    dplyr::summarise(
      DatasetCount = dplyr::n_distinct(Dataset),
      SupportedDatasets = paste(sort(unique(Dataset)), collapse = ";"),
      .groups = "drop"
    ) %>%
    dplyr::mutate(
      ConvergenceClass = dplyr::case_when(
        DatasetCount == length(disease_names) ~ "Shared_all_datasets",
        DatasetCount >= 2L ~ "Shared_pairwise_or_multi",
        TRUE ~ "Dataset_specific"
      )
    ) %>%
    dplyr::arrange(dplyr::desc(DatasetCount), Database, Pathway_Description)

  data.table::fwrite(presence, file.path(tables_dir, "TriDiseaseRankX_pathway_convergence.csv"))

  edges <- enrichment %>%
    dplyr::select(Dataset, Database, Pathway_Key, Pathway_Description, p.adjust, MinusLog10Padj, gene_symbols) %>%
    tidyr::separate_rows(gene_symbols, sep = "/") %>%
    dplyr::mutate(Gene = toupper(trimws(gene_symbols))) %>%
    dplyr::filter(!is.na(Gene), nzchar(Gene)) %>%
    dplyr::select(-gene_symbols) %>%
    dplyr::distinct(Dataset, Database, Pathway_Key, Gene, .keep_all = TRUE) %>%
    dplyr::left_join(presence %>% dplyr::select(Pathway_Key, DatasetCount, ConvergenceClass), by = "Pathway_Key")

  data.table::fwrite(edges, file.path(tables_dir, "TriDiseaseRankX_gene_pathway_edges.csv"))

  max_db <- max(1L, length(unique(enrichment$Database)))
  candidate_genes <- unique(clean_symbol(candidate_genes))
  candidate_genes <- candidate_genes[!is.na(candidate_genes)]

  gene_raw <- edges %>%
    dplyr::filter(Gene %in% candidate_genes) %>%
    dplyr::group_by(Gene) %>%
    dplyr::summarise(
      PathwayCount = dplyr::n_distinct(Pathway_Key),
      DatasetPathwaySupport = dplyr::n_distinct(Dataset),
      DatabaseSupport = dplyr::n_distinct(Database),
      ConvergentPathwayCount = dplyr::n_distinct(Pathway_Key[DatasetCount >= 2L]),
      SharedAllDatasetPathwayCount = dplyr::n_distinct(Pathway_Key[DatasetCount == length(disease_names)]),
      MeanMinusLog10Padj = mean(MinusLog10Padj, na.rm = TRUE),
      BestAdjustedP = min(p.adjust, na.rm = TRUE),
      SupportedDatasets = paste(sort(unique(Dataset)), collapse = ";"),
      SupportedDatabases = paste(sort(unique(Database)), collapse = ";"),
      .groups = "drop"
    )

  scores <- data.frame(Gene = candidate_genes, stringsAsFactors = FALSE) %>%
    dplyr::left_join(gene_raw, by = "Gene") %>%
    dplyr::mutate(
      PathwayCount = dplyr::coalesce(as.integer(PathwayCount), 0L),
      DatasetPathwaySupport = dplyr::coalesce(as.integer(DatasetPathwaySupport), 0L),
      DatabaseSupport = dplyr::coalesce(as.integer(DatabaseSupport), 0L),
      ConvergentPathwayCount = dplyr::coalesce(as.integer(ConvergentPathwayCount), 0L),
      SharedAllDatasetPathwayCount = dplyr::coalesce(as.integer(SharedAllDatasetPathwayCount), 0L),
      MeanMinusLog10Padj = dplyr::coalesce(as.numeric(MeanMinusLog10Padj), 0),
      BestAdjustedP = dplyr::coalesce(as.numeric(BestAdjustedP), 1),
      SupportedDatasets = dplyr::coalesce(SupportedDatasets, ""),
      SupportedDatabases = dplyr::coalesce(SupportedDatabases, ""),
      S_PathwayCount = minmax01(log1p(PathwayCount)),
      S_DatasetPathwaySupport = pmin(pmax(DatasetPathwaySupport / length(disease_names), 0), 1),
      S_DatabaseSupport = pmin(pmax(DatabaseSupport / max_db, 0), 1),
      S_PathwaySignificance = minmax01(MeanMinusLog10Padj),
      S_PathwayConvergence = minmax01(log1p(ConvergentPathwayCount)),
      S_Pathway = rowMeans(cbind(
        S_PathwayCount,
        S_DatasetPathwaySupport,
        S_DatabaseSupport,
        S_PathwaySignificance,
        S_PathwayConvergence
      )),
      Pathway_Rank = rank(-S_Pathway, ties.method = "min")
    ) %>%
    dplyr::arrange(Pathway_Rank, Gene)

  data.table::fwrite(scores, file.path(tables_dir, "TriDiseaseRankX_gene_pathway_scores.csv"))

  top_plot <- presence %>%
    dplyr::left_join(
      enrichment %>%
        dplyr::group_by(Pathway_Key) %>%
        dplyr::summarise(MinPadj = min(p.adjust, na.rm = TRUE), .groups = "drop"),
      by = "Pathway_Key"
    ) %>%
    dplyr::arrange(dplyr::desc(DatasetCount), MinPadj)

  top_plot <- utils::head(top_plot, min(as.integer(top_terms_plot), nrow(top_plot))) %>%
    dplyr::mutate(
      MinusLog10Padj = -log10(pmax(MinPadj, 1e-300)),
      Label = factor(Pathway_Description, levels = rev(unique(Pathway_Description)))
    )

  if (nrow(top_plot)) {
    p <- ggplot2::ggplot(top_plot, ggplot2::aes(x = MinusLog10Padj, y = Label)) +
      ggplot2::geom_point(ggplot2::aes(size = DatasetCount, shape = Database), alpha = 0.85) +
      ggplot2::theme_bw(base_size = 11) +
      ggplot2::labs(
        title = "Top convergent enriched pathways",
        x = "-log10(min adjusted P)",
        y = "Pathway",
        size = "Dataset support",
        shape = "Database"
      )
    ggplot2::ggsave(file.path(plots_dir, "Top_convergent_pathways.pdf"), p, width = 10, height = 7)
    ggplot2::ggsave(file.path(plots_dir, "Top_convergent_pathways.png"), p, width = 10, height = 7, dpi = 300)
  }

  wb <- openxlsx::createWorkbook()
  openxlsx::addWorksheet(wb, "Gene_Pathway_Scores")
  openxlsx::writeData(wb, "Gene_Pathway_Scores", scores)
  openxlsx::addWorksheet(wb, "Pathway_Convergence")
  openxlsx::writeData(wb, "Pathway_Convergence", presence)
  openxlsx::addWorksheet(wb, "Gene_Pathway_Edges")
  openxlsx::writeData(wb, "Gene_Pathway_Edges", edges)
  openxlsx::addWorksheet(wb, "All_Enrichment")
  openxlsx::writeData(wb, "All_Enrichment", enrichment)
  openxlsx::saveWorkbook(
    wb,
    file.path(output_dir, "TriDiseaseRankX_pathway_convergence_results.xlsx"),
    overwrite = TRUE
  )

  report <- c(
    "TriDiseaseRank-X pathway evidence module",
    paste0("Datasets: ", paste(disease_names, collapse = ", ")),
    paste0("Candidate genes scored: ", nrow(scores)),
    paste0("Significant enrichment rows: ", nrow(enrichment)),
    paste0("Unique pathway keys: ", nrow(presence)),
    paste0("Databases contributing evidence: ", paste(sort(unique(enrichment$Database)), collapse = ", ")),
    "S_Pathway = equal mean of five bounded components:",
    "  pathway-count support, dataset support, database support, enrichment significance, convergence support.",
    "Enrichment backgrounds are dataset-specific tested-gene universes where available."
  )
  writeLines(report, file.path(output_dir, "TriDiseaseRankX_pathway_report.txt"))

  list(
    gene_scores = scores,
    enrichment = enrichment,
    convergence = presence,
    edges = edges,
    databases = sort(unique(enrichment$Database))
  )
}
