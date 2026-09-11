options(stringsAsFactors = FALSE)

# FIX (KEGG 'Timeout of 60 seconds was reached'): raise R's default 60 s
# download timeout for ALL remote reads (KEGG REST, GEO, STRING). rest.kegg.jp
# can be slow or intermittently unreachable from some networks/regions.
options(timeout = max(600, getOption("timeout")))

# install only missing packages
install_missing <- function(pkgs, installer, ...) {
  missing_pkgs <- pkgs[!vapply(pkgs, requireNamespace, logical(1), quietly = TRUE)]
  if (length(missing_pkgs)) {
    message("Installing ", length(missing_pkgs), " package(s): ",
            paste(missing_pkgs, collapse = ", "))
    installer(missing_pkgs, ...)
  } else {
    message("All ", length(pkgs), " packages already installed.")
  }
  invisible(missing_pkgs)
}

cores <- max(1L, parallel::detectCores() - 1L)

# BiocManager
if (!requireNamespace("BiocManager", quietly = TRUE))
  install.packages("BiocManager", repos = "https://cloud.r-project.org")

# CRAN
cran_pkgs <- c(
  "pheatmap", "ggplot2", "ggfortify", "VennDiagram", "igraph",
  "dplyr", "tidyr", "stringr", "RColorBrewer", "pROC",
  "data.table", "cowplot", "scales"
)
install_missing(cran_pkgs, install.packages,
                repos = "https://cloud.r-project.org",
                Ncpus = cores)

# Bioconductor
bioc_pkgs <- c(
  "GEOquery", "limma", "GO.db", "impute",
  "ragene11sttranscriptcluster.db", "AnnotationDbi",
  "clusterProfiler", "org.Rn.eg.db", "ReactomePA",
  "STRINGdb", "ComplexHeatmap"
)
install_missing(bioc_pkgs, BiocManager::install,
                update = FALSE, ask = FALSE,
                Ncpus = cores)
# WGCNA if missing
if (!requireNamespace("WGCNA", quietly = TRUE)) install.packages("WGCNA")

BiocManager::install("preprocessCore")
# Library load
suppressPackageStartupMessages({
  library(GEOquery)
  library(limma)
  library(pheatmap)
  library(ggplot2)
  library(VennDiagram)
  library(ggfortify)
  library(WGCNA)
  library(igraph)
  library(dplyr)
  library(tidyr)
  library(stringr)
  library(RColorBrewer)
  library(pROC)
  library(data.table)
  library(cowplot)
  library(clusterProfiler)
  library(org.Rn.eg.db)
  library(ReactomePA)
  library(STRINGdb)
  library(ComplexHeatmap)
  library(ragene11sttranscriptcluster.db)
  library(AnnotationDbi)
})
# Output directories 
out_dir   <- "figures"    
results_dir <- "results"

dir.create(out_dir,     showWarnings = FALSE, recursive = TRUE)
dir.create(results_dir, showWarnings = FALSE, recursive = TRUE)

# Helper function: save ggplot base plot to PNG
save_png <- function(filename, plot_expr, width = 1800, height = 1500, res = 300) {
  filepath <- file.path(out_dir, filename)
  png(filepath, width = width, height = height, res = res)
  if (is.function(plot_expr)) plot_expr() else print(plot_expr)
  dev.off()
  message(sprintf("Saved: %s", filepath))
  invisible(filepath)
}

# A single seed for full reproducibility (Reviewer 2 Q2 / Q4)
set.seed(20240901)



# 1. Download GSE58438 (Reviewer 2 Q1: explicit platform / sample reporting)
gse_id <- "GSE58438"

# caching: getGEO() downloads to tempdir(); set destdir = for persistence.
gse <- getGEO(gse_id, GSEMatrix = TRUE, getGPL = TRUE)
eset <- gse[[1]]

# metadata block printed + saved to results Reviewer 2 Q1 
platform_id <- annotation(eset)
platform_title <- tryCatch(
  pData(eset)$platform_title[1],
  error = function(e) NA_character_
)
sample_info <- data.frame(
  GSM         = sampleNames(eset),
  Title       = pData(eset)$title,
  Group_orig  = pData(eset)$characteristics_ch1.1,
  stringsAsFactors = FALSE
)
write.csv(sample_info, file.path(results_dir, "GSE58438_sample_metadata.csv"),
          row.names = FALSE)

meta_log <- c(
  paste("Dataset                  :", gse_id),
  paste("Platform (GEO)           :", platform_id),
  paste("Platform title           :", platform_title),
  paste("Total samples in series  :", ncol(eset)),
  paste("Samples used in analysis :", "9 (5 Control + 4 AKI)"),
  paste("GSM accessions used      :",
        paste(c("GSM1411057","GSM1411058","GSM1411059","GSM1411060","GSM1411061",
                "GSM1411067","GSM1411068","GSM1411069","GSM1411070"), collapse=", ")),
  paste("Normalization            : values from GEO Matrix are already log2-RMA-scaled"),
  paste("log2-transform check     : max(exprSet)=~16, min(exprSet)=~2 (consistent with log2)"),
  paste("limma model              : ~ 0 + group (no intercept, group-means parameterization)"),
  paste("limma contrast           : AKI - Control"),
  paste("DEG thresholds (primary) : |log2FC| > 2 and raw P.Value < 0.05 (matches paper, 8,596 DEGs)"),
  paste("DEG thresholds (sens.)   : |log2FC| > 2 and BH-FDR (adj.P.Val) < 0.05"),
  paste("Multiple-testing method  : Benjamini-Hochberg (BH), reported in 'adj.P.Val'"),
  paste("R version                :", R.version.string)
)
writeLines(meta_log, file.path(results_dir, "analysis_metadata.txt"))
cat(meta_log, sep = "\n")


## 2. Subsetting

samples <- c("GSM1411057","GSM1411058","GSM1411059","GSM1411060","GSM1411061",
             "GSM1411067","GSM1411068","GSM1411069","GSM1411070")

available <- intersect(samples, colnames(eset))
missing   <- setdiff(samples, available)

if (length(missing) > 0)
  stop(sprintf("Missing samples in the downloaded ExpressionSet: %s",
               paste(missing, collapse = ", ")))

eset <- eset[, samples]

# Expression matrix
exprSet <- Biobase::exprs(eset)

# Reviewer 2 Q1: verify log2 scale. RMA-processed Affymetrix ST arrays
# from GEO are already log2-transformed; we sanity-check here.
log2_check <- list(
  min_value      = round(min(exprSet, na.rm = TRUE), 3),
  max_value      = round(max(exprSet, na.rm = TRUE), 3),
  median_value   = round(median(exprSet, na.rm = TRUE), 3),
  looks_log2     = (max(exprSet, na.rm = TRUE) < 25)  
)
print(log2_check)

## 3. Differential expression with limma (Reviewer 2 Q1: explicit model + FDR)
group  <- factor(c(rep("Control", 5), rep("AKI", 4)),
                 levels = c("Control", "AKI"))
design <- model.matrix(~ 0 + group)
colnames(design) <- levels(group)

fit  <- lmFit(exprSet, design)
contrast.matrix <- makeContrasts(AKIvsControl = AKI - Control, levels = design)
fit2 <- contrasts.fit(fit, contrast.matrix)
fit2 <- eBayes(fit2)

# Full DEG table with BOTH raw p and BH-FDR (Reviewer 2 Q1)
deg <- topTable(fit2, coef = "AKIvsControl", number = Inf,
                adjust.method = "BH")
deg$ProbeID <- rownames(deg)

# Primary filter (matches paper exactly): |logFC|>2 & P<0.05 
deg_filtered <- deg[abs(deg$logFC) > 2 & deg$P.Value < 0.05, ]

# Sensitivity analysis filter: |logFC|>2 & BH-FDR<0.05 
deg_fdr <- deg[abs(deg$logFC) > 2 & deg$adj.P.Val < 0.05, ]

# Persist DEG outputs (Reviewer 2 Minor 7: data availability)
write.csv(deg,          file.path(results_dir, "DEG_full_all_genes.csv"),
          row.names = FALSE)
write.csv(deg_filtered, file.path(results_dir, "DEG_filtered_logFC2_p05.csv"),
          row.names = FALSE)
write.csv(deg_fdr,      file.path(results_dir, "DEG_sensitivity_FDR05.csv"),
          row.names = FALSE)

# Reviewer 2 Q1 summary table 
deg_summary <- data.frame(
  Filter = c("|logFC|>2 & P<0.05  (paper, primary)",
             "|logFC|>2 & BH-FDR<0.05  (sensitivity)"),
  Total_DEGs      = c(nrow(deg_filtered), nrow(deg_fdr)),
  Upregulated     = c(sum(deg_filtered$logFC > 0),
                      sum(deg_fdr$logFC > 0)),
  Downregulated   = c(sum(deg_filtered$logFC < 0),
                      sum(deg_fdr$logFC < 0))
)
write.csv(deg_summary, file.path(results_dir, "DEG_summary.csv"),
          row.names = FALSE)
print(deg_summary)

# Up / down lists
deg_up   <- deg_filtered[deg_filtered$logFC > 0, ]
deg_down <- deg_filtered[deg_filtered$logFC < 0, ]

# 4. Annotation columns for heatmaps
annotation_col <- data.frame(Group = group, row.names = colnames(exprSet))


# =============================================================================
# 5. Figure 1A — Heatmap of DEGs (PNG)
# =============================================================================
plot_genes <- intersect(rownames(deg_filtered), rownames(exprSet))
stopifnot(length(plot_genes) >= 2)

save_png("Fig1A_DEG_heatmap.png", {
  pheatmap(exprSet[plot_genes, ],
           scale = "row",
           annotation_col = annotation_col,
           show_rownames = FALSE,
           show_colnames = TRUE,
           main = "Figure 1A. DEG heatmap (|logFC|>2, P<0.05)")
}, width = 2000, height = 2400)


# =============================================================================
# 6. Figure 1B — Volcano plot (PNG)
# =============================================================================
deg$threshold <- as.factor(abs(deg$logFC) > 2 & deg$P.Value < 0.05)
volcano <- ggplot(deg, aes(x = logFC, y = -log10(P.Value), color = threshold)) +
  geom_point(alpha = 0.6, size = 1.2) +
  scale_color_manual(values = c("grey", "red"),
                     labels = c("NS", "|logFC|>2 & P<0.05")) +
  geom_vline(xintercept = c(-2, 2), linetype = "dashed", color = "grey40") +
  geom_hline(yintercept = -log10(0.05), linetype = "dashed", color = "grey40") +
  theme_minimal(base_size = 12) +
  labs(title = "Figure 1B. Volcano plot (AKI vs Control)",
       x = "log2 Fold Change", y = "-log10(P value)", color = "Threshold")
save_png("Fig1B_volcano.png", volcano, width = 2000, height = 1800)


# =============================================================================
# 7. Figure 1C — PCA plot (PNG)
# =============================================================================
pca_res <- prcomp(t(exprSet), scale. = TRUE)
pca_df  <- data.frame(PC1 = pca_res$x[, 1], PC2 = pca_res$x[, 2], Group = group)
var_explained <- round(100 * pca_res$sdev^2 / sum(pca_res$sdev^2), 1)
pca_plot <- ggplot(pca_df, aes(PC1, PC2, color = Group)) +
  geom_point(size = 4, alpha = 0.85) +
  stat_ellipse(level = 0.95, type = "t") +
  scale_color_manual(values = c("Control" = "#377EB8", "AKI" = "#E41A1C")) +
  theme_minimal(base_size = 12) +
  labs(title = "Figure 1C. PCA (AKI vs Control)",
       x = sprintf("PC1 (%.1f%%)", var_explained[1]),
       y = sprintf("PC2 (%.1f%%)", var_explained[2]))
save_png("Fig1C_PCA.png", pca_plot, width = 1800, height = 1500)


# =============================================================================
# 8. Figure 1D — Up vs Down Venn (PNG)
#    (Up and Down are mutually exclusive by construction; we visualize this
#    fact for figure completeness, mirroring the original Figure 1D.)
# =============================================================================
save_png("Fig1D_UpVsDown_Venn.png", {
  grid.newpage()
  venn.plot <- draw.pairwise.venn(
    area1     = nrow(deg_up),
    area2     = nrow(deg_down),
    cross.area = 0,
    category  = c("Upregulated", "Downregulated"),
    fill      = c("#E41A1C", "#377EB8"),
    alpha     = 0.5,
    cat.pos   = c(0, 0),
    cat.dist  = c(0.05, 0.05),
    cat.cex   = 1.2,
    cex       = 1.2
  )
}, width = 1500, height = 1500)


# =============================================================================
# 9. Probe-to-gene symbol mapping (rat Gene 1.1 ST array)
# =============================================================================
probe_ids_all <- rownames(deg)
sym_all <- mapIds(ragene11sttranscriptcluster.db,
                  keys  = probe_ids_all,
                  column = "SYMBOL",
                  keytype = "PROBEID",
                  multiVals = "first")
deg$GeneSymbol <- sym_all[deg$ProbeID]

# Upregulated-only table for STRING submission (reviewer-noted choice)
deg_up_anno <- deg_filtered[deg_filtered$logFC > 0, ]
deg_up_anno$GeneSymbol <- mapIds(ragene11sttranscriptcluster.db,
                                 keys  = rownames(deg_up_anno),
                                 column = "SYMBOL",
                                 keytype = "PROBEID",
                                 multiVals = "first")

# Remove unmapped probes before STRING submission (Reviewer 2 Q2: transparency)
deg_up_anno_clean <- deg_up_anno[!is.na(deg_up_anno$GeneSymbol) &
                                   deg_up_anno$GeneSymbol != "", ]
write.csv(deg_up_anno_clean,
          file.path(results_dir, "Upregulated_DEGs_for_STRING.csv"),
          row.names = FALSE)

# Also export the upregulated probe table (original behavior preserved)
up_genes <- data.frame(
  ProbeID    = rownames(deg_up_anno),
  GeneSymbol = deg_up_anno$GeneSymbol,
  logFC      = deg_up_anno$logFC,
  P.Value    = deg_up_anno$P.Value,
  adj.P.Val  = deg_up_anno$adj.P.Val     # Reviewer 2 Q1: now also FDR
)
write.csv(up_genes,
          file.path(results_dir, "Upregulated_genes_with_probes.csv"),
          row.names = FALSE)


# =============================================================================
# 10. Figure 2 — Top-50 genes heatmap (PNG)
# =============================================================================
top_genes <- head(rownames(deg_filtered[order(abs(deg_filtered$logFC),
                                              decreasing = TRUE), ]), 50)
top_genes <- intersect(top_genes, rownames(exprSet))

# Add gene-symbol row labels for readability
top_anno <- data.frame(
  ProbeID = top_genes,
  Symbol = mapIds(ragene11sttranscriptcluster.db,
                  keys = top_genes, column = "SYMBOL",
                  keytype = "PROBEID", multiVals = "first")
)
save_png("Fig2_top50_heatmap.png", {
  pheatmap(exprSet[top_genes, ],
           scale = "row",
           annotation_col = annotation_col,
           labels_row = ifelse(is.na(top_anno$Symbol), top_anno$ProbeID,
                               top_anno$Symbol),
           show_rownames = TRUE,
           show_colnames = TRUE,
           main = "Figure 2. Top-50 DEG heatmap")
}, width = 2000, height = 2400)


# =============================================================================
# 11. GSEA — missing in original (Reviewer 2 Q1)
#     Implemented with clusterProfiler: Hallmark + KEGG + Reactome.
#     Uses all genes ranked by limma t-statistic (not filtered).
# =============================================================================
rank_list <- deg$t
names(rank_list) <- deg$GeneSymbol
rank_list <- sort(rank_list[!is.na(names(rank_list))], decreasing = TRUE)

# Map SYMBOL -> ENTREZID for KEGG / Reactome
eg <- bitr(names(rank_list), fromType = "SYMBOL",
           toType = "ENTREZID", OrgDb = org.Rn.eg.db)
rank_entrez <- rank_list[eg$SYMBOL]
names(rank_entrez) <- eg$ENTREZID
rank_entrez <- sort(rank_entrez[!duplicated(names(rank_entrez))],
                    decreasing = TRUE)

# --- KEGG GSEA (cites KEGG; Editor Q3 reminder printed) ---
cat("NOTE: KEGG is © Kanehisa Laboratories. Cite Kanehisa et al., 2021",
    "(Nucleic Acids Res. 41:D586-590) and Kanehisa & Goto, 2000.",
    "Access date:", format(Sys.Date(), "%Y-%m-%d"), "\n")

# --- Deterministic, duplicate-free ranked vector ----------------------------
o <- order(-abs(rank_entrez), names(rank_entrez))
rank_entrez <- rank_entrez[o]
rank_entrez <- rank_entrez[!duplicated(names(rank_entrez))]  # 1 Entrez = 1 entry
rank_entrez <- rank_entrez[order(-rank_entrez, names(rank_entrez))]  # tie-break

cat(sprintf("[GSEA] versions -> DOSE %s | clusterProfiler %s | ReactomePA %s | fgsea %s\n",
            as.character(packageVersion("DOSE")),
            as.character(packageVersion("clusterProfiler")),
            as.character(packageVersion("ReactomePA")),
            tryCatch(as.character(packageVersion("fgsea")),
                     error = function(e) "not installed")))
cat("[GSEA] NOTE: on the gson/enrichit-based Bioconductor stack (2026+),\n",
    "      gseKEGG()/gsePathway() no longer expose an 'eps' argument. The\n",
    "      'P-values are less than 1e-10' notice is then EXPECTED and harmless:\n",
    "      it flags only pathways whose multilevel p-value estimate hit the\n",
    "      1e-10 floor (i.e., the strongest hits). On older stacks eps=0 is\n",
    "      applied automatically via the formals() checks below.\n")

# --- Build args; add eps/seed only if this clusterProfiler version has them --
gsea_args <- list(geneList      = rank_entrez,
                  organism      = "rno",
                  keyType       = "kegg",
                  nPerm         = 1000,
                  minGSSize     = 10,
                  maxGSSize     = 500,
                  pvalueCutoff  = 0.25,
                  pAdjustMethod = "BH",
                  verbose       = FALSE)
fmls <- names(formals(clusterProfiler::gseKEGG))
if ("eps"  %in% fmls) gsea_args$eps  <- 0   # exact tail estimation
if ("seed" %in% fmls) gsea_args$seed <- TRUE
if (!("eps" %in% fmls))
  message("[GSEA] clusterProfiler has no 'eps' formal (expected on the ",
          "gson/enrichit stack) -> the p<1e-10 notice cannot be suppressed ",
          "via arguments; it is harmless (see version banner above).")

set.seed(42)
gsea_kegg <- NULL
for (kegg_att in 1:3) {
  gsea_kegg <- tryCatch(do.call(clusterProfiler::gseKEGG, gsea_args),
                        error = function(e) {
                          message(sprintf("KEGG GSEA attempt %d failed: %s",
                                          kegg_att, e$message))
                          NULL })
  if (!is.null(gsea_kegg)) break
  if (kegg_att < 3) {
    message("  retrying in 30 s (rest.kegg.jp can be slow or unreachable)...")
    Sys.sleep(30)
  }
}
if (is.null(gsea_kegg))
  message("[KEGG] KEGG GSEA skipped after 3 failed attempts (see the exact ",
          "error messages above for the cause -- connectivity, blocked rest.kegg.jp, ",
          "or an argument mismatch on this package version). ",
          "All other pipeline sections are unaffected.")

if (!is.null(gsea_kegg) && nrow(gsea_kegg) > 0) {
  save_png("Fig_GSEA_KEGG_dotplot.png", {
    print(dotplot(gsea_kegg, showCategory = 20,
                  title = "GSEA — KEGG pathways (AKI vs Control)"))
  }, width = 2000, height = 1800)
  write.csv(as.data.frame(gsea_kegg),
            file.path(results_dir, "GSEA_KEGG_results.csv"),
            row.names = FALSE)
}

# --- Reactome GSEA (FIX: eps/seed now applied here too, same formals() check
#     as the KEGG call -- the original call never received eps = 0, which is
#     why the 'P-values are less than 1e-10' warning fired a 2nd time) ---
react_args <- list(geneList      = rank_entrez,
                   organism      = "rat",
                   nPerm         = 1000,
                   minGSSize     = 10,
                   maxGSSize     = 500,
                   pvalueCutoff  = 0.25,
                   pAdjustMethod = "BH",
                   verbose       = FALSE)
fmls_r <- names(formals(ReactomePA::gsePathway))
if ("eps"  %in% fmls_r) react_args$eps  <- 0   # exact tail estimation
if ("seed" %in% fmls_r) react_args$seed <- TRUE
if (!("eps" %in% fmls_r))
  message("[GSEA] ReactomePA has no 'eps' formal (expected on the ",
          "gson/enrichit stack) -> the p<1e-10 notice is harmless.")
set.seed(42)
gsea_react <- tryCatch(do.call(ReactomePA::gsePathway, react_args),
                       error = function(e) { message("Reactome GSEA failed: ", e$message); NULL })
if (!is.null(gsea_react) && nrow(gsea_react) > 0) {
  save_png("Fig_GSEA_Reactome_dotplot.png", {
    print(dotplot(gsea_react, showCategory = 20,
                  title = "GSEA — Reactome pathways (AKI vs Control)"))
  }, width = 2000, height = 1800)
  write.csv(as.data.frame(gsea_react),
            file.path(results_dir, "GSEA_Reactome_results.csv"),
            row.names = FALSE)
}

# =============================================================================
# 12. PPI network via STRINGdb (Reviewer 2 Q2: documented parameters)
#     - Organism: Rattus norvegicus (taxid 10116)
#     - Score threshold: medium confidence 400 (STRING default)
#     - Only upregulated genes submitted (documented design choice)
#     Requires from upstream: results_dir, out_dir, save_png(),
#     deg_up_anno_clean, exprSet, group
# =============================================================================
suppressPackageStartupMessages({
  library(httr); library(igraph); library(dplyr); library(tidyr); library(tibble)
  library(pheatmap); library(pROC)
})

## ============================================================================
## PURE HELPERS  (no side effects; save this block as pipeline_helpers.R to
##                run the offline unit tests at the bottom of this answer)
## ============================================================================

# Normalize edge tables from STRINGdb (from/to/combined_score) OR the REST
# tsv API (#node1_string_id / node2_string_id / score).
normalize_string_edges <- function(edges) {
  stopifnot(is.data.frame(edges))
  cn <- colnames(edges)
  find_col <- function(patterns) {
    for (p in patterns) {
      hit <- cn[grepl(p, cn)]
      if (length(hit) > 0) return(hit[1])
    }
    NA_character_
  }
  from_col <- find_col(c("^from$", "^#?node1_string_id$", "^stringId_A$"))
  to_col   <- find_col(c("^to$",   "^node2_string_id$",   "^stringId_B$"))
  sc_col   <- find_col(c("^combined_score$", "^score$"))   # never nscore/fscore/...
  if (is.na(from_col) || is.na(to_col) || is.na(sc_col))
    stop("Unrecognized STRING column names: ", paste(cn, collapse = ", "))
  out <- data.frame(
    from           = as.character(edges[[from_col]]),
    to             = as.character(edges[[to_col]]),
    combined_score = suppressWarnings(as.numeric(edges[[sc_col]])),
    stringsAsFactors = FALSE
  )
  out[stats::complete.cases(out), , drop = FALSE]
}

# Collapse A->B / B->A duplicates, keeping the HIGHEST combined_score.
dedup_undirected <- function(edges) {
  key <- paste(pmin(edges$from, edges$to), pmax(edges$from, edges$to), sep = "|")
  o     <- order(-edges$combined_score)
  edges <- edges[o, , drop = FALSE]
  key   <- key[o]
  edges[!duplicated(key), c("from", "to", "combined_score"), drop = FALSE]
}

filter_within <- function(edges, keep_ids) {
  keep_ids <- unique(as.character(keep_ids))
  edges[edges$from %in% keep_ids & edges$to %in% keep_ids, , drop = FALSE]
}

fetch_string_rest <- function(ids, species = 10116, score = 400, chunk = 250) {
  ids <- unique(as.character(ids))
  if (length(ids) == 0) return(NULL)
  out <- list()
  for (start in seq(1, length(ids), by = chunk)) {
    end   <- min(start + chunk - 1, length(ids))
    batch <- ids[start:end]
    cat(sprintf("  [REST] ids %d-%d / %d\n", start, end, length(ids)))
    resp <- tryCatch(
      httr::POST("https://string-db.org/api/tsv/interaction-partners",
                 body = list(identifiers     = paste(batch, collapse = "\r\n"),
                             species         = species,
                             required_score  = score,
                             limit           = 1000,   # partners per input protein
                             caller_identity = "rat_AKI_DEG_pipeline"),
                 encode = "form", httr::timeout(180)),
      error = function(e) { message("  [REST] request failed: ", e$message); NULL })
    if (is.null(resp)) next
    if (httr::status_code(resp) != 200) {
      message("  [REST] HTTP ", httr::status_code(resp)); next
    }
    txt <- httr::content(resp, as = "text", encoding = "UTF-8")
    df  <- tryCatch(read.delim(text = txt, header = TRUE, sep = "\t",
                               check.names = FALSE, quote = "",
                               stringsAsFactors = FALSE),
                    error = function(e) NULL)
    if (is.null(df) || ncol(df) < 3 || nrow(df) == 0) next  # error pages, empty
    out[[length(out) + 1L]] <- df
    Sys.sleep(0.5)  # be polite to the public API
  }
  if (length(out) == 0) return(NULL)
  unique(do.call(rbind, out))
}

predict_venn <- function(mirdb_file, mirwalk_file, targetscan_file,
                         out_csv, out_png) {
  stopifnot(file.exists(mirdb_file), file.exists(mirwalk_file),
            file.exists(targetscan_file))
  mirdb      <- read.csv(mirdb_file,      stringsAsFactors = FALSE)$miRNA
  mirwalk    <- read.csv(mirwalk_file,    stringsAsFactors = FALSE)$miRNA
  targetscan <- read.csv(targetscan_file, stringsAsFactors = FALSE)$miRNA
  v <- VennDiagram::venn.diagram(
    x = list(miRDB = mirdb, miRWalk = mirwalk, TargetScan = targetscan),
    filename = NULL, fill = c("#E41A1C", "#377EB8", "#4DAF4A"),
    alpha = 0.5, cat.cex = 1.1, cex = 1.1, main = "miRNA 3-database overlap")
  grDevices::png(out_png, width = 1800, height = 1500, res = 300)
  grid::grid.draw(v); grDevices::dev.off()
  common <- Reduce(intersect, list(mirdb, mirwalk, targetscan))
  utils::write.csv(data.frame(miRNA = common), out_csv, row.names = FALSE)
  message(sprintf("Saved: %s ; common miRNAs: %d", out_png, length(common)))
  invisible(common)
}

run_mirna_mrna_spearman <- function(ct_df, gene_cols, mirna_cols,
                                    out_csv, out_png, min_n = 5) {
  stopifnot(all(c(gene_cols, mirna_cols) %in% colnames(ct_df)))
  grid_ <- expand.grid(Gene = gene_cols, miRNA = mirna_cols,
                       stringsAsFactors = FALSE)
  res <- do.call(rbind, lapply(seq_len(nrow(grid_)), function(i) {
    g <- grid_$Gene[i]; m <- grid_$miRNA[i]
    x <- suppressWarnings(as.numeric(ct_df[[g]]))
    y <- suppressWarnings(as.numeric(ct_df[[m]]))
    ok <- is.finite(x) & is.finite(y)
    if (sum(ok) < min_n || length(unique(x[ok])) < 2 || length(unique(y[ok])) < 2) {
      return(data.frame(Gene = g, miRNA = m, rho = NA_real_,
                        p.value = NA_real_, n = sum(ok)))
    }
    ct <- suppressWarnings(stats::cor.test(x[ok], y[ok],
                                           method = "spearman", exact = FALSE))
    data.frame(Gene = g, miRNA = m, rho = unname(ct$estimate),
               p.value = ct$p.value, n = sum(ok))
  }))
  res$adj.p.value <- stats::p.adjust(res$p.value, method = "BH")
  utils::write.csv(res, out_csv, row.names = FALSE)
  
  mat <- as.matrix(tibble::column_to_rownames(
    tidyr::pivot_wider(res, names_from = miRNA, values_from = rho), "Gene"))
  if (nrow(mat) >= 2 && ncol(mat) >= 2) {
    grDevices::png(out_png, width = 2000, height = 1500, res = 300)
    pheatmap::pheatmap(mat, cluster_rows = TRUE, cluster_cols = TRUE,
                       color = grDevices::colorRampPalette(
                         c("#377EB8", "white", "#E41A1C"))(50),
                       main = "Spearman correlation: mRNA vs miRNA")
    grDevices::dev.off()
    message("Saved: ", out_png)
  } else message("Heatmap skipped (need >=2 genes and >=2 miRNAs).")
  invisible(res)
}

# NOTE (Reviewer 2 Q5): exploratory only; same cohort as qPCR, no independent
# clinical validation. direction="<" => marker expected HIGHER in levels[2].
exploratory_roc <- function(df, marker, group_levels = c("Sham", "IRI")) {
  stopifnot(marker %in% colnames(df), "group" %in% colnames(df))
  d <- df[, c("group", marker)]
  x <- suppressWarnings(as.numeric(d[[marker]]))
  d <- d[is.finite(x), , drop = FALSE]
  d$pred <- x[is.finite(x)]
  d$grp  <- factor(d$group, levels = group_levels)
  d <- d[!is.na(d$grp), , drop = FALSE]
  tab <- table(d$grp)
  if (any(tab < 2))
    stop("Need >=2 samples per group (got ",
         paste(names(tab), as.integer(tab), sep = "=", collapse = ", "), ")")
  roc_obj <- pROC::roc(response = d$grp, predictor = d$pred,
                       levels = group_levels, direction = "<", quiet = TRUE)
  ci_obj <- as.numeric(pROC::ci.auc(roc_obj, conf.level = 0.95))
  best <- as.data.frame(pROC::coords(roc_obj, "best",
                                     ret = c("threshold", "sensitivity", "specificity"),
                                     transpose = FALSE))
  spec   <- best$specificity[1]
  lr_pos <- if (is.finite(spec) && spec < 1) best$sensitivity[1] / (1 - spec) else Inf
  x0 <- d$pred[d$grp == group_levels[1]]  # controls
  x1 <- d$pred[d$grp == group_levels[2]]  # cases
  # P(AUC == 0.5) == Mann-Whitney p-value on the marker between groups
  p_mw <- tryCatch(stats::wilcox.test(x1, x0, exact = FALSE)$p.value,
                   error = function(e) NA_real_)
  data.frame(marker = marker, n = nrow(d),
             n_controls = length(x0), n_cases = length(x1),
             AUC = as.numeric(roc_obj$auc),
             AUC_lower95 = ci_obj[1], AUC_upper95 = ci_obj[3],
             p_value = p_mw, cutoff = best$threshold[1],
             sensitivity = best$sensitivity[1], specificity = spec,
             LR_plus = lr_pos, stringsAsFactors = FALSE)
}

# FIX (fatal crash): the two lines previously pasted here --
#     geneList <- geneList[order(-geneList, names(geneList))]
#     egob     <- GSEA(geneList, ..., eps = 0)
# aborted the whole pipeline with:
#     Error in eval(ei, envir) : object 'geneList' not found
# because (1) no object named `geneList` exists in this script -- the ranked
# vector is `rank_entrez` (Section 11) -- and (2) a literal `...` is not valid
# inside a top-level call. Both fixes they intended are ALREADY implemented in
# Section 11 and need no code here:
#   (a) Ties warning  -> `rank_entrez` is deterministically tie-broken by gene
#       name via order(-rank_entrez, names(rank_entrez)) BEFORE the GSEA calls,
#       so tied genes keep a fixed, reproducible order. The fgsea notice itself
#       cannot be removed by sorting (ties exist in the t-stat VALUES) and is
#       safe to document and ignore.
#   (b) p<1e-10 warning -> `eps = 0` is now passed to BOTH gseKEGG() and
#       gsePathway() via formals() checks (see gsea_args / react_args below in
#       Section 11).
## ========================== END PURE HELPERS ================================

# ---- 12a. Initialize STRINGdb (multi-version fallback; cache ~GB on 1st run)
init_string_db <- function(species = 10116, score_threshold = 400, cache_dir) {
  dir.create(cache_dir, showWarnings = FALSE, recursive = TRUE)
  for (v in c("12.1", "12.0", "12", "11.5", "11.0")) {
    db <- tryCatch(
      STRINGdb$new(version = v, species = species,
                   score_threshold = score_threshold,
                   input_directory = cache_dir),
      error = function(e) {
        message("  [STRING] version ", v, " unavailable (", conditionMessage(e), ")")
        NULL
      })
    if (!is.null(db)) { message("[STRING] Using STRINGdb v", v); return(db) }
  }
  stop("[STRING] Could not initialise STRINGdb with any known version.")
}
string_db <- init_string_db(10116, 400, file.path(results_dir, "STRING_cache"))
# NOTE: do NOT call set_background() — species is pinned in the constructor.

# ---- 12b. Map upregulated gene symbols -> STRING_id
sym_all <- unique(na.omit(as.character(deg_up_anno_clean$GeneSymbol)))
map_in  <- unique(deg_up_anno_clean[!is.na(deg_up_anno_clean$GeneSymbol),
                                    c("GeneSymbol", "logFC"), drop = FALSE])
mapped <- string_db$map(map_in, "GeneSymbol", removeUnmappedRows = TRUE)
write.csv(mapped, file.path(results_dir, "STRING_mapped_upregulated.csv"),
          row.names = FALSE)

unmapped_syms <- setdiff(sym_all, mapped$GeneSymbol)
write.csv(data.frame(GeneSymbol = unmapped_syms),
          file.path(results_dir, "STRING_unmapped_symbols.csv"), row.names = FALSE)

# ---- 12c. Diagnostics (FIX: length(), not nrow(), on a vector)
n_submitted <- length(sym_all)
n_mapped    <- length(unique(na.omit(mapped$STRING_id)))
cat(sprintf("[STRING] Submitted symbols : %d\n", n_submitted))
cat(sprintf("[STRING] Mapped to STRING_id: %d (%.1f%%)\n",
            n_mapped, 100 * n_mapped / max(n_submitted, 1)))
if (n_mapped / max(n_submitted, 1) < 0.5)
  warning("[STRING] <50% mapping rate — check species taxid / identifier type.")
print(head(mapped$STRING_id, 5))

string_ids <- unique(mapped$STRING_id[!is.na(mapped$STRING_id) &
                                        nzchar(mapped$STRING_id)])
stopifnot(length(string_ids) > 0)
cat(sprintf("[STRING] Unique STRING_ids to query: %d\n", length(string_ids)))

# ---- 12d. Chunked STRINGdb query
chunk_size <- 50
chunks <- split(string_ids, ceiling(seq_along(string_ids) / chunk_size))
cat(sprintf("[STRING] Querying interactions in %d chunks (~%d ids/chunk)...\n",
            length(chunks), chunk_size))
string_edges_list <- list()
for (i in seq_along(chunks)) {
  if (i %% 10 == 0 || i == length(chunks))
    cat(sprintf("  chunk %d/%d\n", i, length(chunks)))
  edges_i <- tryCatch(string_db$get_interactions(chunks[[i]]),
                      error = function(e) {
                        message("  chunk ", i, " failed: ", e$message); NULL })
  if (!is.null(edges_i) && nrow(edges_i) > 0)
    string_edges_list[[length(string_edges_list) + 1L]] <- edges_i
}
string_edges_stringdb <- if (length(string_edges_list))
  unique(do.call(rbind, string_edges_list)) else NULL

# ---- 12e/f. Choose source, normalize, 12g. filter within input set + dedup
raw_edges <- if (is.null(string_edges_stringdb) || nrow(string_edges_stringdb) == 0) {
  message("[STRING] STRINGdb returned no edges -> REST POST fallback.")
  fetch_string_rest(string_ids)
} else string_edges_stringdb
if (is.null(raw_edges) || nrow(raw_edges) == 0)
  stop("[STRING] No interactions from STRINGdb or REST. ",
       "Check internet access / STRING service status.")

string_edges_raw <- normalize_string_edges(raw_edges)
string_edges     <- dedup_undirected(filter_within(string_edges_raw, string_ids))
cat(sprintf("[STRING] Raw edges: %d | normalized: %d | between input proteins: %d\n",
            nrow(raw_edges), nrow(string_edges_raw), nrow(string_edges)))
covered <- unique(c(string_edges$from, string_edges$to))
cat(sprintf("[STRING] Coverage: %d / %d input proteins in >=1 edge\n",
            length(covered), length(string_ids)))
write.csv(string_edges,
          file.path(results_dir, "STRING_edge_list_upregulated.csv"),
          row.names = FALSE)

# ---- 12h. igraph + centralities + reproducible hub lists --------------------
# IMPORTANT: clusterProfiler/DOSE/ggtree export an S4 generic `simplify` that
# masks igraph's. Every igraph call below is fully qualified on purpose.
stopifnot(nrow(string_edges) >= 2)

g <- igraph::graph_from_data_frame(string_edges, directed = FALSE)
g <- igraph::simplify(g, remove.multiple = TRUE, remove.loops = TRUE,
                      edge.attr.comb = list(combined_score = "max"))
igraph::V(g)$STRING_id <- igraph::V(g)$name   # preserve IDs before renaming
id2sym <- setNames(mapped$GeneSymbol, mapped$STRING_id)
igraph::V(g)$name <- ifelse(igraph::V(g)$name %in% names(id2sym),
                            id2sym[igraph::V(g)$name], igraph::V(g)$name)
cat(sprintf("[STRING] igraph: %d vertices, %d edges\n",
            igraph::vcount(g), igraph::ecount(g)))

centralities <- data.frame(
  Gene        = igraph::V(g)$name,
  Degree      = as.numeric(igraph::degree(g)),
  Betweenness = as.numeric(igraph::betweenness(g, normalized = TRUE)),
  Closeness   = suppressWarnings(
    as.numeric(igraph::closeness(g, normalized = TRUE))),
  Eigenvector = as.numeric(igraph::eigen_centrality(g)$vector),
  stringsAsFactors = FALSE)
centralities <- centralities[order(-centralities$Degree,
                                   -centralities$Betweenness,
                                   centralities$Gene), ]
rownames(centralities) <- NULL
write.csv(centralities,
          file.path(results_dir, "hub_gene_centrality_metrics.csv"),
          row.names = FALSE)

hub90 <- head(centralities$Gene, min(90, nrow(centralities)))
hub10 <- head(centralities$Gene, min(10, nrow(centralities)))
writeLines(hub90, file.path(results_dir, "hub_genes_90.txt"))
writeLines(hub10, file.path(results_dir, "hub_genes_10.txt"))
# =============================================================================
# 13. Hub-gene supplementary figures (unchanged logic, package-safe)
# =============================================================================
plot_top_n_bar <- function(metric, n, fill_col) {
  df <- centralities %>%
    arrange(desc(.data[[metric]])) %>%
    slice_head(n = n) %>%
    mutate(Gene = factor(Gene, levels = rev(Gene)))
  ggplot(df, aes(x = .data[[metric]], y = Gene)) +
    geom_col(fill = fill_col, color = "black", width = 0.7) +
    theme_minimal(base_size = 11) +
    labs(title = sprintf("Top-%d hub genes by %s", n, metric),
         x = metric, y = NULL)
}
metric_cols <- c(Degree = "#E41A1C", Betweenness = "#377EB8",
                 Closeness = "#4DAF4A", Eigenvector = "#984EA3")
for (m in names(metric_cols))
  for (n in c(20, 30))
    save_png(sprintf("HubFig_top%d_%s.png", n, m),
             plot_top_n_bar(m, n, metric_cols[[m]]),
             width = 1800, height = 1500)

# =============================================================================
# 14. Hub-of-hub refinement (FIX: filter to the 10 hub genes only)
# =============================================================================
if (!requireNamespace("VennDiagram", quietly = TRUE))
  BiocManager::install("VennDiagram")
suppressPackageStartupMessages(library(VennDiagram))
suppressMessages(futile.logger::flog.threshold(futile.logger::ERROR,
                                               name = "VennDiagramLogger"))

hub10_mapped <- string_db$map(data.frame(GeneSymbol = hub10,
                                         stringsAsFactors = FALSE),
                              "GeneSymbol", removeUnmappedRows = TRUE)
miss <- setdiff(hub10, hub10_mapped$GeneSymbol)
if (length(miss))
  message("[STRING] Hub genes not mapped: ", paste(miss, collapse = ", "))
hub10_ids <- unique(hub10_mapped$STRING_id)
raw10 <- tryCatch(string_db$get_interactions(hub10_ids), error = function(e) NULL)
if (is.null(raw10) || nrow(raw10) == 0)
  raw10 <- fetch_string_rest(hub10_ids, chunk = max(length(hub10_ids), 1))

if (is.null(raw10)) {
  message("[STRING] No interactions among the 10 hubs; writing empty table.")
  hub10_edges <- data.frame(from = character(), to = character(),
                            combined_score = numeric(),
                            from_symbol = character(), to_symbol = character())
} else {
  hub10_edges <- dedup_undirected(
    filter_within(normalize_string_edges(raw10), hub10_ids))
  id2sym10 <- setNames(hub10_mapped$GeneSymbol, hub10_mapped$STRING_id)
  hub10_edges$from_symbol <- unname(id2sym10[hub10_edges$from])
  hub10_edges$to_symbol   <- unname(id2sym10[hub10_edges$to])
}
write.csv(hub10_edges, file.path(results_dir, "STRING_edge_list_hub10.csv"),
          row.names = FALSE)
cat(sprintf("[STRING] Hub-10 internal edges: %d\n", nrow(hub10_edges)))

# Final 6 hub genes retained for experimental validation (paper-curated; rat
# ortholog symbols). Kept hardcoded for manuscript traceability.
hub6 <- c("BRCA1", "RAD51", "MCM7", "EXO1", "RFC3", "TP53")
writeLines(hub6, file.path(results_dir, "hub_genes_6_validated.txt"))

# =============================================================================
# 15. miRNA target prediction — Venn across miRDB / miRWalk / TargetScan
#     (predict_venn() defined in the PURE HELPERS block above)
# =============================================================================
# Template — un-comment and supply the three exported TSVs to regenerate Fig 4D:
# common_mirnas <- predict_venn(
#   "mirdb_export.tsv", "mirwalk_export.tsv", "targetscan_export.tsv",
#   out_csv = file.path(results_dir, "miRNA_intersection.csv"),
#   out_png = file.path(out_dir, "Fig4D_miRNA_Venn.png"))

gene_mirna_table <- data.frame(
  Gene = c(rep("BRCA1", 3), rep("RAD51", 2), rep("EXO1", 1),
           rep("MCM7", 1), rep("RFC3", 1), rep("TP53", 1)),
  miRNA = c("rno-miR-103-3p", "rno-miR-107-3p", "rno-miR-101b-3p",
            "rno-miR-92a-2-5p", "rno-miR-935",
            "rno-miR-103-3p", "rno-miR-107-3p",
            "rno-miR-935", "rno-miR-103-3p"),
  Predicted_by  = "miRDB, miRWalk, TargetScan",
  Evidence_type = "predicted (3-database overlap)",
  stringsAsFactors = FALSE)
write.csv(gene_mirna_table,
          file.path(results_dir, "gene_miRNA_prediction_table.csv"),
          row.names = FALSE)

# =============================================================================
# 16. Spearman miRNA-mRNA correlation (uses helper; NA-safe, ties-safe)
# =============================================================================
# Template — un-comment when the qPCR Ct matrix is loaded:
# run_mirna_mrna_spearman(qpcr_ct,
#   gene_cols  = c("BRCA1","RAD51","MCM7","EXO1","RFC3","TP53"),
#   mirna_cols = c("rno-miR-103-3p","rno-miR-107-3p","rno-miR-101b-3p",
#                  "rno-miR-92a-2-5p","rno-miR-935"),
#   out_csv = file.path(results_dir, "miRNA_mRNA_spearman.csv"),
#   out_png = file.path(out_dir, "Fig_miRNA_mRNA_corr_heatmap.png"))

# =============================================================================
# 17. Independent validation: GSE9943 (Reviewer 1 Q1)  [REWRITTEN - was FATAL]
# -----------------------------------------------------------------------------
# Ground truth verified against the GEO SOFT file (checked 2026-09-06):
#   * 12 samples = 2 strains x 2 groups x 3 replicates: Brown Norway (BN) and
#     Sprague-Dawley (SD) rats, each with a Control group (renal pedicles NOT
#     clamped, i.e. the sham analog) and an I/R group (45-min bilateral
#     ischemia + 24-h reperfusion).
#   * characteristics_ch1 records ONLY the strain ("Tissue from <strain> rat
#     strain") and is even WRONG for 2 samples (GSM251588 and GSM280411 are
#     labeled "Brown Norway" although their titles say Sprague Dawley).
#   * The treatment group appears ONLY in Sample_title, abbreviated "I/R",
#     which the old grepl("iri|ischem|reperfus|clamp") could never match, and
#     "control" never appears in characteristics_ch1 at all -> the old
#     auto-detection hit its own stop() and aborted the whole pipeline here.
#   * Platform = GPL2996 (ABI Rat Genome Survey, spotted oligo array) -- NOT
#     Affymetrix, so the old ragene10sttranscriptcluster.db fallback could
#     never fire. Only ~24% of probes carry a GEO "Gene Symbol"; the rest are
#     recovered via GB_LIST (first GenBank accession) -> org.Rn.eg.db.
#   * Strain caveat: the discovery cohort GSE58438 used WISTAR rats (Speir et
#     al., PLoS One 2015; PMID 25970334). GSE9943 provides BN and SD only, so
#     we validate in BOTH strains: SD as the primary arm (closest common
#     outbred strain) and BN as an I/R-resistant sensitivity arm; both overlap
#     tables are written.
# The whole section is wrapped in tryCatch(): a validation failure now
# degrades to a message instead of aborting WGCNA (Section 18) and beyond.
# =============================================================================

## Ground-truth GSM -> "strain|group" map (transcribed from GEO Sample_title)
gse9943_map <- c(
  GSM251560 = "BN|Control", GSM251561 = "SD|Control",
  GSM251586 = "BN|Control", GSM251588 = "SD|Control",
  GSM251591 = "BN|IR",      GSM251593 = "BN|IR",
  GSM251594 = "SD|IR",      GSM251595 = "SD|IR",
  GSM280410 = "BN|Control", GSM280411 = "SD|Control",
  GSM280416 = "BN|IR",      GSM280417 = "SD|IR")

val_strain <- "SD"   # primary validation arm used for the Venn figure ("BN" to switch)

## Helper: infer strain/group from the sample title, then reconcile with the
## ground-truth map (map wins on conflict; discrepancies are logged).
infer_group_9943 <- function(gsm, title) {
  hard <- unname(gse9943_map[gsm])
  t    <- tolower(ifelse(is.na(title), "", title))
  strn <- if (grepl("brown norway|\\bbn\\b", t)) "BN" else
    if (grepl("sprague|\\bsd\\b",       t)) "SD" else NA_character_
  grp  <- if (grepl("control|sham",            t)) "Control" else
    if (grepl("i/r|\\bir\\b|iri|ischem", t)) "IR"      else NA_character_
  out  <- c(strn, grp)
  if (!is.na(hard)) {
    hv <- strsplit(hard, "|", fixed = TRUE)[[1]]
    if (is.na(out[1])) out[1] <- hv[1]
    if (is.na(out[2])) out[2] <- hv[2]
    if (!identical(out, hv))
      message(sprintf("[GSE9943] %s: title infers %s/%s but map says %s/%s -> keeping map",
                      gsm, out[1], out[2], hv[1], hv[2]))
  }
  out
}

tryCatch(local({
  
  GSE_VAL <- "GSE9943"
  gse_val <- tryCatch(getGEO(GSE_VAL, GSEMatrix = TRUE, getGPL = TRUE),
                      error = function(e) { message(GSE_VAL, " download failed: ",
                                                    e$message); NULL })
  if (is.null(gse_val)) stop("GEO download failed", call. = FALSE)
  
  eset_val  <- gse_val[[1]]
  pheno_val <- Biobase::pData(eset_val)
  expr_val  <- Biobase::exprs(eset_val)
  stopifnot(ncol(expr_val) == nrow(pheno_val))
  
  ## Audit dump BEFORE any inference (Reviewer 2 traceability)
  write.csv(pheno_val, file.path(results_dir, "GSE9943_pData_full.csv"),
            row.names = TRUE)
  
  ## Group labels: parse titles + ground-truth map (see block comment above)
  gsm_ids <- sampleNames(eset_val)
  titles  <- as.character(pheno_val$title)
  inf2    <- lapply(seq_along(gsm_ids),
                    function(i) infer_group_9943(gsm_ids[i], titles[i]))
  strains <- vapply(inf2, function(x) x[1], character(1))
  grps    <- vapply(inf2, function(x) x[2], character(1))
  names(strains) <- names(grps) <- gsm_ids
  
  val_table <- data.frame(GSM = gsm_ids, Title = titles,
                          Strain = unname(strains), Group = unname(grps),
                          stringsAsFactors = FALSE)
  write.csv(val_table, file.path(results_dir, "GSE9943_group_assignment.csv"),
            row.names = FALSE)
  cat("[GSE9943] Group assignment (from Sample_title + GEO map):\n")
  print(val_table)
  
  if (any(is.na(strains)) || any(is.na(grps)))
    stop("Unassigned strain/group for: ",
         paste(gsm_ids[is.na(strains) | is.na(grps)], collapse = ", "), call. = FALSE)
  
  cat(sprintf("[GSE9943] SD: %d Control / %d IR | BN: %d Control / %d IR\n",
              sum(strains == "SD" & grps == "Control"),
              sum(strains == "SD" & grps == "IR"),
              sum(strains == "BN" & grps == "Control"),
              sum(strains == "BN" & grps == "IR")))
  
  ## FIX (kept): limma needs complete rows; do NOT coerce negatives to NA.
  expr_val <- expr_val[complete.cases(expr_val), , drop = FALSE]
  
  ## Probe -> SYMBOL. GPL2996 = ABI Rat Genome Survey: use the GEO "Gene
  ## Symbol" column first, then fall back to GB_LIST (first GenBank
  ## accession) mapped through org.Rn.eg.db.
  fdat    <- tryCatch(Biobase::fData(eset_val), error = function(e) NULL)
  sym_lut <- setNames(rep(NA_character_, nrow(expr_val)), rownames(expr_val))
  if (!is.null(fdat) && "Gene Symbol" %in% colnames(fdat)) {
    s <- as.character(fdat[rownames(expr_val), "Gene Symbol"])
    s[s %in% c("null", "", "NA")] <- NA_character_
    sym_lut <- setNames(s, rownames(expr_val))
  }
  if (!is.null(fdat) && "GB_LIST" %in% colnames(fdat)) {
    na_ids <- names(sym_lut)[is.na(sym_lut)]
    if (length(na_ids)) {
      gb  <- as.character(fdat[na_ids, "GB_LIST"])
      acc <- rep(NA_character_, length(gb))
      vg  <- !is.na(gb) & nzchar(gb)
      acc[vg] <- trimws(sub("^([^,;[:space:]]+).*$", "\\1", gb[vg]))
      acc[acc %in% c("null", "")] <- NA_character_
      ok  <- !is.na(acc)
      eg2 <- if (any(ok)) suppressMessages(tryCatch(
        AnnotationDbi::select(org.Rn.eg.db, keys = acc[ok],
                              keytype = "ACCNUM", columns = "SYMBOL"),
        error = function(e) NULL)) else NULL
      if (!is.null(eg2) && nrow(eg2) > 0) {
        eg2  <- eg2[!duplicated(eg2$ACCNUM) & !is.na(eg2$SYMBOL), , drop = FALSE]
        fill <- eg2$SYMBOL[match(acc, eg2$ACCNUM)]
        names(fill) <- na_ids
        keep <- !is.na(fill)
        sym_lut[names(fill)[keep]] <- fill[keep]
      }
    }
  }
  cat(sprintf("[GSE9943] probes with a gene symbol: %d / %d (%.1f%%)\n",
              sum(!is.na(sym_lut)), length(sym_lut),
              100 * sum(!is.na(sym_lut)) / length(sym_lut)))
  
  ## Discovery DEG symbols (|logFC|>2 & P<0.05), from Sections 3 + 9.
  ## (FIX: the old code overlapped against symbols of ALL tested probes via
  ##  `deg`, which inflated the discovery list; now restricted to DEGs.)
  sym_disc  <- setNames(deg$GeneSymbol, deg$ProbeID)
  disc_syms <- na.omit(unique(unname(sym_disc[rownames(deg_filtered)])))
  if (!length(disc_syms))
    stop("Discovery DEG symbol list is empty; check Sections 3 and 9.", call. = FALSE)
  
  ## Per-strain validation DEGs + overlap (SD = primary, BN = sensitivity)
  val_summary <- list()
  for (st in c("SD", "BN")) {
    sel <- names(strains)[strains == st & grps %in% c("Control", "IR")]
    gst <- table(droplevels(factor(grps[sel])))
    if (length(sel) < 4 || min(gst) < 2) {
      message("[GSE9943] strain ", st, ": insufficient replicates (", length(sel),
              " samples); skipped."); next
    }
    g_st <- factor(grps[sel], levels = c("Control", "IR"))
    e_st <- expr_val[, sel, drop = FALSE]
    
    d_st <- model.matrix(~ 0 + g_st)
    colnames(d_st) <- levels(g_st)
    f_st <- limma::lmFit(e_st, d_st)
    cont <- setNames(rep(0, ncol(d_st)), colnames(d_st))
    cont["IR"] <- 1; cont["Control"] <- -1          # IR - Control
    f_st <- limma::eBayes(limma::contrasts.fit(f_st, cont))
    dv   <- limma::topTable(f_st, number = Inf, adjust.method = "BH")
    
    dvf <- dv[abs(dv$logFC) > 2 & dv$P.Value < 0.05, , drop = FALSE]
    dvf$GeneSymbol <- unname(sym_lut[rownames(dvf)])
    write.csv(dv,  file.path(results_dir,
                             sprintf("GSE9943_%s_DEG_full.csv", st)), row.names = TRUE)
    write.csv(dvf, file.path(results_dir,
                             sprintf("GSE9943_%s_DEG_filtered.csv", st)), row.names = TRUE)
    
    vs <- na.omit(unique(dvf$GeneSymbol))
    ov <- intersect(disc_syms, vs)
    val_summary[[st]] <- data.frame(
      Strain = st, n_samples = length(sel),
      n_DEG_probes = nrow(dvf), n_DEG_symbols = length(vs),
      overlap_with_GSE58438 = length(ov), stringsAsFactors = FALSE)
    write.csv(data.frame(Validated_gene = ov),
              file.path(results_dir,
                        sprintf("GSE58438_vs_GSE9943_%s_overlap.csv", st)),
              row.names = FALSE)
    cat(sprintf("[validation/%s] discovery: %d | GSE9943-%s: %d symbols | overlap: %d\n",
                st, length(disc_syms), st, length(vs), length(ov)))
    
    ## Primary Venn figure uses the strain-matched arm (val_strain)
    if (st == val_strain && length(disc_syms) > 0 && length(vs) > 0) {
      save_png("Fig_GSE9943_validation_venn.png", {
        grid::grid.newpage()
        grid::grid.draw(draw.pairwise.venn(
          area1 = length(disc_syms), area2 = length(vs),
          cross.area = length(ov),
          category = c("GSE58438 (discovery)",
                       sprintf("GSE9943 %s (validation)", st)),
          fill = c("#E41A1C", "#377EB8"), alpha = 0.5,
          cat.pos = c(-30, 30), cat.dist = c(0.05, 0.05),
          cat.cex = 1.2, cex = 1.2))
      }, width = 1800, height = 1500)
    }
  }
  
  if (length(val_summary)) {
    val_sum_df <- do.call(rbind, val_summary)
    write.csv(val_sum_df,
              file.path(results_dir, "GSE9943_validation_summary.csv"),
              row.names = FALSE)
    print(val_sum_df)
  }
}), error = function(e) {
  message("[GSE9943] Validation section failed (pipeline continues): ",
          conditionMessage(e))
})

# =============================================================================
# 18. WGCNA module discovery (FIXES: MEs coloring, trait alignment, blocks)
# =============================================================================
suppressPackageStartupMessages(library(WGCNA))
enableWGCNAThreads()

datExpr <- t(as.matrix(exprSet))
if (is.null(rownames(datExpr)))
  rownames(datExpr) <- sprintf("S%02d", seq_len(nrow(datExpr)))
gsg <- goodSamplesGenes(datExpr, verbose = 0)
if (!gsg$allOK) {
  cat(sprintf("[WGCNA] goodSamplesGenes removed %d samples / %d genes\n",
              sum(!gsg$goodSamples), sum(!gsg$goodGenes)))
  datExpr <- datExpr[gsg$goodSamples, gsg$goodGenes, drop = FALSE]
}
# FIX: explicit levels so Control=0/AKI=1 is guaranteed; align with samples.
trait <- as.numeric(factor(group, levels = c("Control", "AKI"))) - 1
if (!gsg$allOK) trait <- trait[gsg$goodSamples]
stopifnot(length(trait) == nrow(datExpr))

powers <- 1:20
sft <- pickSoftThreshold(datExpr, powerVector = powers,
                         networkType = "unsigned", verbose = 3)
chosen_power <- 6   # value used in the paper (documented, Reviewer 2 Minor 7)
r2_at_power <- sft$fitIndices$SFT.R.sq[sft$fitIndices$Power == chosen_power]
if (!is.na(r2_at_power) && r2_at_power < 0.80)
  message(sprintf("[WGCNA] NOTE: scale-free R^2 at power %d is %.2f (<0.80); ",
                  chosen_power, r2_at_power),
          "sft$powerEstimate = ", sft$powerEstimate,
          " — consider reporting both.")

save_png("WGCNA_soft_threshold.png", {
  par(mfrow = c(1, 2))
  plot(sft$fitIndices$Power, -sign(sft$fitIndices$slope) * sft$fitIndices$SFT.R.sq,
       xlab = "Soft Threshold (power)",
       ylab = "Scale Free Topology Model Fit (signed R^2)",
       type = "n", main = "Soft-threshold selection")
  text(sft$fitIndices$Power, -sign(sft$fitIndices$slope) * sft$fitIndices$SFT.R.sq,
       labels = powers, col = "red")
  plot(sft$fitIndices$Power, sft$fitIndices$mean.k.,
       xlab = "Soft Threshold (power)", ylab = "Mean Connectivity",
       type = "l", main = "Mean connectivity decay")
}, width = 2200, height = 1100)

net <- blockwiseModules(datExpr,
                        power = chosen_power, TOMType = "unsigned",
                        minModuleSize = 30, reassignThreshold = 0,
                        mergeCutHeight = 0.25, numericLabels = TRUE,
                        pamRespectsDendro = FALSE,
                        saveTOMs = TRUE,   # large files; set FALSE if disk-bound
                        saveTOMFileBase = file.path(results_dir, "GSE58438_TOM"),
                        verbose = 3)

moduleColors <- labels2colors(net$colors)
# FIX: blockwiseModules returns numeric labels (ME1, ...); relabel to colors
# so MEs[, "MEblue"] works (previous version errored here).
MEs <- net$MEs
colnames(MEs) <- paste0("ME", labels2colors(
  as.numeric(gsub("^ME", "", colnames(MEs)))))
MEs <- orderMEs(as.data.frame(MEs))

save_png("WGCNA_dendro_colors.png", {                    # FIX: plot ALL blocks
  n_blocks <- length(net$dendrograms)
  par(mfrow = c(n_blocks, 1))
  for (b in seq_len(n_blocks))
    plotDendroAndColors(net$dendrograms[[b]], moduleColors[net$blockGenes[[b]]],
                        "Module colors", dendroLabels = FALSE, hang = 0.03,
                        addGuide = TRUE, guideHang = 0.05,
                        main = sprintf("Block %d", b))
}, width = 2200, height = 800 * length(net$dendrograms) + 300)

trait_mat <- matrix(trait, ncol = 1, dimnames = list(NULL, "AKI_vs_Control"))
moduleTraitCor    <- cor(as.matrix(MEs), trait_mat, use = "p")
moduleTraitPvalue <- corPvalueStudent(moduleTraitCor, nSamples = nrow(datExpr))
textMat <- paste(signif(moduleTraitCor, 2), "\n(",
                 signif(moduleTraitPvalue, 1), ")", sep = "")
dim(textMat) <- dim(moduleTraitCor)

save_png("WGCNA_module_trait_heatmap.png", {
  par(mar = c(4, 8, 2, 1))
  labeledHeatmap(Matrix = moduleTraitCor, xLabels = colnames(trait_mat),
                 yLabels = rownames(moduleTraitCor),
                 ySymbols = rownames(moduleTraitCor), colorLabels = FALSE,
                 colors = blueWhiteRed(50), textMatrix = textMat,
                 setStdMargins = FALSE, cex.text = 0.8,
                 main = "Module-trait relationships")
}, width = 1400, height = 2000)

module_of_interest <- "blue"
me_col <- paste0("ME", module_of_interest)
moduleGenes <- moduleColors == module_of_interest
if (sum(moduleGenes) > 0 && me_col %in% colnames(MEs)) {
  kME_blue <- as.numeric(cor(datExpr[, moduleGenes, drop = FALSE],
                             MEs[[me_col]], use = "p"))
  n_hub <- min(10, sum(moduleGenes))
  hubGenes <- colnames(datExpr)[moduleGenes][order(-kME_blue)][seq_len(n_hub)]
  writeLines(hubGenes, file.path(results_dir, "WGCNA_blue_hub10.txt"))
  cat("Top WGCNA blue-module hub genes (probe IDs):\n"); print(hubGenes)
} else message("[WGCNA] Module 'blue' absent in this dataset; skipped.")

# =============================================================================
# 19. ROC analysis — exploratory (Reviewer 2 Q5). helper defined above;
#     FIX: previous roc.test(roc_obj) call was invalid (needs two ROC curves);
#     AUC p-value now via Mann-Whitney (equivalent to H0: AUC = 0.5).
# =============================================================================
# Template call (un-comment when qPCR data loaded):
# roc_results <- do.call(rbind, lapply(
#   c("BRCA1","RAD51","MCM7","EXO1","RFC3","TP53"),
#   function(m) exploratory_roc(qpcr_df, m)))
# write.csv(roc_results, file.path(results_dir, "ROC_exploratory_summary.csv"),
#           row.names = FALSE)

# =============================================================================
# 20. Session info + artifact manifest
# =============================================================================
writeLines(capture.output(sessionInfo()),
           file.path(results_dir, "sessionInfo.txt"))
manifest <- list.files(out_dir, pattern = "\\.png$", full.names = FALSE)
results_manifest <- list.files(results_dir, full.names = FALSE)
writeLines(c("# Figures", paste0("- ", manifest), "",
             "# Result tables", paste0("- ", results_manifest)),
           file.path(results_dir, "MANIFEST.txt"))
cat("\n=== Pipeline finished ===\n")
cat("PNG figures -> ", normalizePath(out_dir), "\n")
cat("CSV tables  -> ", normalizePath(results_dir), "\n")
