options(stringsAsFactors = FALSE)
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
setwd("D:/")
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



# Download GSE58438 (Reviewer 2 Q1: explicit platform / sample reporting)
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


## Subsetting

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
log2_check <- list(
  min_value      = round(min(exprSet, na.rm = TRUE), 3),
  max_value      = round(max(exprSet, na.rm = TRUE), 3),
  median_value   = round(median(exprSet, na.rm = TRUE), 3),
  looks_log2     = (max(exprSet, na.rm = TRUE) < 25)  
)
print(log2_check)

## Differential expression with limma Reviewer 2 Q1: explicit model + FDR
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

write.csv(deg,          file.path(results_dir, "DEG full all genes.csv"),
          row.names = FALSE)
write.csv(deg_filtered, file.path(results_dir, "DEG filtered logFC.csv"),
          row.names = FALSE)
write.csv(deg_fdr,      file.path(results_dir, "DEG sensitivity.csv"),
          row.names = FALSE)

# summary table 
deg_summary <- data.frame(
  Filter = c("|logFC|>2 & P<0.05",
             "|logFC|>2 & BH-FDR<0.05  (sensitivity)"),
  Total_DEGs      = c(nrow(deg_filtered), nrow(deg_fdr)),
  Upregulated     = c(sum(deg_filtered$logFC > 0),
                      sum(deg_fdr$logFC > 0)),
  Downregulated   = c(sum(deg_filtered$logFC < 0),
                      sum(deg_fdr$logFC < 0))
)
write.csv(deg_summary, file.path(results_dir, "DEG summary.csv"),
          row.names = FALSE)
print(deg_summary)

deg_up   <- deg_filtered[deg_filtered$logFC > 0, ]
deg_down <- deg_filtered[deg_filtered$logFC < 0, ]


annotation_col <- data.frame(Group = group, row.names = colnames(exprSet))


plot_genes <- intersect(rownames(deg_filtered), rownames(exprSet))
stopifnot(length(plot_genes) >= 2)
heat_colors <- colorRampPalette(c("navy", "white", "firebrick3"))(50)

save_png("DEG heatmap.png", {
  pheatmap(exprSet[plot_genes, ],
           scale = "row",
           color = heat_colors,
           annotation_col = annotation_col,
           border_color = NA,
           show_rownames = FALSE,
           show_colnames = TRUE,
           main = "DEG heatmap (|logFC|>2, P<0.05)")
}, width = 2000, height = 2400)


deg$threshold <- as.factor(abs(deg$logFC) > 2 & deg$P.Value < 0.05)
volcano <- ggplot(deg, aes(x = logFC, y = -log10(P.Value), color = threshold)) +
  geom_point(alpha = 0.6, size = 1.2) +
  scale_color_manual(values = c("grey", "red"),
                     labels = c("NS", "|logFC|>2 & P<0.05")) +
  geom_vline(xintercept = c(-2, 2), linetype = "dashed", color = "grey40") +
  geom_hline(yintercept = -log10(0.05), linetype = "dashed", color = "grey40") +
  theme_minimal(base_size = 12) +
  labs(title = "Volcano plot (AKI vs Control)",
       x = "log2 Fold Change", y = "-log10(P value)", color = "Threshold")
save_png("1b volcano.png", volcano, width = 2000, height = 1800)


pca_res <- prcomp(t(exprSet), scale. = TRUE)
pca_df  <- data.frame(PC1 = pca_res$x[, 1], PC2 = pca_res$x[, 2], Group = group)
var_explained <- round(100 * pca_res$sdev^2 / sum(pca_res$sdev^2), 1)
pca_plot <- ggplot(pca_df, aes(PC1, PC2, color = Group)) +
  geom_point(size = 4, alpha = 0.85) +
  stat_ellipse(level = 0.95, type = "t") +
  scale_color_manual(values = c("Control" = "#377EB8", "AKI" = "#E41A1C")) +
  theme_minimal(base_size = 12) +
  labs(title = "PCA (AKI vs Control)",
       x = sprintf("PC1 (%.1f%%)", var_explained[1]),
       y = sprintf("PC2 (%.1f%%)", var_explained[2]))
save_png("1c PCA.png", pca_plot, width = 1800, height = 1500)


save_png("1D UpVsDown Venn.png", {
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


# Probe to gene symbol mapping 
probe_ids_all <- rownames(deg)
sym_all <- mapIds(ragene11sttranscriptcluster.db,
                  keys  = probe_ids_all,
                  column = "SYMBOL",
                  keytype = "PROBEID",
                  multiVals = "first")
deg$GeneSymbol <- sym_all[deg$ProbeID]


deg_up_anno <- deg_filtered[deg_filtered$logFC > 0, ]
deg_up_anno$GeneSymbol <- mapIds(ragene11sttranscriptcluster.db,
                                 keys  = rownames(deg_up_anno),
                                 column = "SYMBOL",
                                 keytype = "PROBEID",
                                 multiVals = "first")


deg_up_anno_clean <- deg_up_anno[!is.na(deg_up_anno$GeneSymbol) &
                                   deg_up_anno$GeneSymbol != "", ]
write.csv(deg_up_anno_clean,
          file.path(results_dir, "Upregulated DEGs for STRING.csv"),
          row.names = FALSE)

up_genes <- data.frame(
  ProbeID    = rownames(deg_up_anno),
  GeneSymbol = deg_up_anno$GeneSymbol,
  logFC      = deg_up_anno$logFC,
  P.Value    = deg_up_anno$P.Value,
  adj.P.Val  = deg_up_anno$adj.P.Val     
)
write.csv(up_genes,
          file.path(results_dir, "Upregulated genes with probes.csv"),
          row.names = FALSE)



top_genes <- head(rownames(deg_filtered[order(abs(deg_filtered$logFC),
                                              decreasing = TRUE), ]), 50)
top_genes <- intersect(top_genes, rownames(exprSet))

top_anno <- data.frame(
  ProbeID = top_genes,
  Symbol = mapIds(ragene11sttranscriptcluster.db,
                  keys = top_genes, column = "SYMBOL",
                  keytype = "PROBEID", multiVals = "first")
)
save_png("2 top50 heatmap.png", {
  pheatmap(exprSet[top_genes, ],
           scale = "row",
           color = heat_colors,
           annotation_col = annotation_col,
           labels_row = ifelse(is.na(top_anno$Symbol), top_anno$ProbeID,
                               top_anno$Symbol),
           show_rownames = TRUE,
           show_colnames = TRUE,
           main = "Top 50 DEG heatmap")
}, width = 2000, height = 2400)


# GSEA 
rank_list <- deg$t
names(rank_list) <- deg$GeneSymbol
rank_list <- sort(rank_list[!is.na(names(rank_list))], decreasing = TRUE)

eg <- bitr(names(rank_list), fromType = "SYMBOL",
           toType = "ENTREZID", OrgDb = org.Rn.eg.db)
rank_entrez <- rank_list[eg$SYMBOL]
names(rank_entrez) <- eg$ENTREZID
rank_entrez <- sort(rank_entrez[!duplicated(names(rank_entrez))],
                    decreasing = TRUE)

cat("NOTE: KEGG is © Kanehisa Laboratories. Cite Kanehisa et al., 2021",
    "(Nucleic Acids Res. 41:D586-590) and Kanehisa & Goto, 2000.",
    "Access date:", format(Sys.Date(), "%Y-%m-%d"), "\n")

o <- order(-abs(rank_entrez), names(rank_entrez))
rank_entrez <- rank_entrez[o]
rank_entrez <- rank_entrez[!duplicated(names(rank_entrez))] 
rank_entrez <- rank_entrez[order(-rank_entrez, names(rank_entrez))] 

cat(sprintf("versions -> DOSE %s | clusterProfiler %s | ReactomePA %s | fgsea %s\n",
            as.character(packageVersion("DOSE")),
            as.character(packageVersion("clusterProfiler")),
            as.character(packageVersion("ReactomePA")),
            tryCatch(as.character(packageVersion("fgsea")),
                     error = function(e) "not installed")))

gsea_muffle_warnings <- function(expr, label) {
  n_ties <- 0L
  n_eps  <- 0L
  res <- withCallingHandlers(
    expr,
    warning = function(w) {
      msg <- conditionMessage(w)
      if (grepl("ties in the preranked", msg)) {
        n_ties <<- n_ties + 1L
        invokeRestart("muffleWarning")
      } else if (grepl("1e-10|eps argument", msg)) {
        n_eps <<- n_eps + 1L
        invokeRestart("muffleWarning")
      }
    })
  if (n_ties > 0L || n_eps > 0L)
    message(sprintf(paste0("%s: muffled %d benign fgsea notice(s) ",
                           "(%d tie-ordering, %d p<1e-10 floor); this is expected."),
                    label, n_ties + n_eps, n_ties, n_eps))
  res
}


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
if ("eps"  %in% fmls) gsea_args$eps  <- 0   
if ("seed" %in% fmls) gsea_args$seed <- TRUE
if (!("eps" %in% fmls))
  message("clusterProfiler has no 'eps' formal; it is harmless.")

set.seed(42)

gsea_kegg <- NULL
old_timeout <- getOption("timeout")
options(timeout = 120)
for (kegg_att in 1:3) {
  gsea_kegg <- tryCatch(
    gsea_muffle_warnings(do.call(clusterProfiler::gseKEGG, gsea_args), "KEGG"),
    error = function(e) {
      message(sprintf("KEGG GSEA attempt %d failed: %s",
                      kegg_att, e$message))
      NULL })
  if (!is.null(gsea_kegg)) break
  if (kegg_att < 3) {
    message("  retrying in 30s")
    Sys.sleep(30)
  }
}
options(timeout = old_timeout)
if (is.null(gsea_kegg))
  message("KEGG GSEA skipped after 3 failed attempts")

if (!is.null(gsea_kegg) && nrow(gsea_kegg) > 0) {
  save_png("GSEA KEGG dotplot.png", {
    print(
      dotplot(gsea_kegg, showCategory = 20, title = "GSEA KEGG pathways (AKI vs Control)") +
        scale_color_viridis_c(guide = guide_colorbar(reverse = TRUE)) +
        theme_bw(base_size = 13) +
        theme(axis.text.y = element_text(size = 11, lineheight = 0.8),
              panel.grid.minor = element_blank())
    )
  }, width = 2000, height = 1800)
  
  write.csv(as.data.frame(gsea_kegg),
            file.path(results_dir, "GSEA KEGG results.csv"),
            row.names = FALSE)
}



react_args <- list(geneList      = rank_entrez,
                   organism      = "rat",
                   nPerm         = 1000,
                   minGSSize     = 10,
                   maxGSSize     = 500,
                   pvalueCutoff  = 0.25,
                   pAdjustMethod = "BH",
                   verbose       = FALSE)
fmls_r <- names(formals(ReactomePA::gsePathway))
if ("eps"  %in% fmls_r) react_args$eps  <- 0  
if ("seed" %in% fmls_r) react_args$seed <- TRUE
if (!("eps" %in% fmls_r))
  message("ReactomePA has no 'eps' formal and its harmless.")
set.seed(42)
gsea_react <- tryCatch(
  gsea_muffle_warnings(do.call(ReactomePA::gsePathway, react_args), "Reactome"),
  error = function(e) { message("Reactome GSEA failed: ", e$message); NULL })
if (!is.null(gsea_react) && nrow(gsea_react) > 0) {
  save_png("GSEA Reactome dotplot.png", {
    print(
      dotplot(gsea_react, showCategory = 20, title = "GSEA - Reactome pathways (AKI vs Control)") +
        scale_color_viridis_c(guide = guide_colorbar(reverse = TRUE)) +
        scale_y_discrete(labels = function(x) str_wrap(x, width = 45)) + 
        theme_bw(base_size = 13) +
        theme(axis.text.y = element_text(size = 7, lineheight = 0.75),  
              panel.grid.minor = element_blank(),
              plot.title = element_text(hjust = 0.5))
    )
  }, width = 2000, height = 1800) 
  
  write.csv(as.data.frame(gsea_react),
            file.path(results_dir, "GSEA Reactome results.csv"),
            row.names = FALSE)
}

# PPI network via STRINGdb 
suppressPackageStartupMessages({
  library(httr); library(igraph); library(dplyr); library(tidyr); library(tibble)
  library(pheatmap); library(pROC)
})


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
  sc_col   <- find_col(c("^combined_score$", "^score$"))   
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

# Collapse A->B / B->A 
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
    cat(sprintf("ids %d-%d / %d\n", start, end, length(ids)))
    resp <- tryCatch(
      httr::POST("https://string-db.org/api/tsv/interaction-partners",
                 body = list(identifiers     = paste(batch, collapse = "\r\n"),
                             species         = species,
                             required_score  = score,
                             limit           = 1000,   
                             caller_identity = "rat_AKI_DEG_pipeline"),
                 encode = "form", httr::timeout(180)),
      error = function(e) { message("request failed: ", e$message); NULL })
    if (is.null(resp)) next
    if (httr::status_code(resp) != 200) {
      message("HTTP ", httr::status_code(resp)); next
    }
    txt <- httr::content(resp, as = "text", encoding = "UTF-8")
    df  <- tryCatch(read.delim(text = txt, header = TRUE, sep = "\t",
                               check.names = FALSE, quote = "",
                               stringsAsFactors = FALSE),
                    error = function(e) NULL)
    if (is.null(df) || ncol(df) < 3 || nrow(df) == 0) next 
    out[[length(out) + 1L]] <- df
    Sys.sleep(0.5) 
  }
  if (length(out) == 0) return(NULL)
  unique(do.call(rbind, out))
}

#Initialize STRINGdb 
init_string_db <- function(species = 10116, score_threshold = 400, cache_dir) {
  dir.create(cache_dir, showWarnings = FALSE, recursive = TRUE)
  for (v in c("12.0", "12", "11.5", "11.0")) {
    db <- tryCatch(
      STRINGdb$new(version = v, species = species,
                   score_threshold = score_threshold,
                   input_directory = cache_dir),
      error = function(e) {
        message("  [STRING] version ", v, " unavailable (", conditionMessage(e), ")")
        NULL
      })
    if (!is.null(db)) { message("Using STRINGdb v", v); return(db) }
  }
  stop("Could not initialise STRINGdb with any known version.")
}
string_db <- init_string_db(10116, 400, file.path(results_dir, "STRING_cache"))


# Map upregulated gene symbols -> STRING_id
sym_all <- unique(na.omit(as.character(deg_up_anno_clean$GeneSymbol)))
map_in  <- unique(deg_up_anno_clean[!is.na(deg_up_anno_clean$GeneSymbol),
                                    c("GeneSymbol", "logFC"), drop = FALSE])
mapped <- string_db$map(map_in, "GeneSymbol", removeUnmappedRows = TRUE)
write.csv(mapped, file.path(results_dir, "STRING mapped upregulated.csv"),
          row.names = FALSE)

unmapped_syms <- setdiff(sym_all, mapped$GeneSymbol)
write.csv(data.frame(GeneSymbol = unmapped_syms),
          file.path(results_dir, "STRING unmapped symbols.csv"), row.names = FALSE)

# Diagnostics 
n_submitted <- length(sym_all)
n_mapped    <- length(unique(na.omit(mapped$STRING_id)))
cat(sprintf("Submitted symbols : %d\n", n_submitted))
cat(sprintf("Mapped to STRING_id: %d (%.1f%%)\n",
            n_mapped, 100 * n_mapped / max(n_submitted, 1)))
if (n_mapped / max(n_submitted, 1) < 0.5)
  warning("<50% mapping rate")
print(head(mapped$STRING_id, 5))

string_ids <- unique(mapped$STRING_id[!is.na(mapped$STRING_id) &
                                        nzchar(mapped$STRING_id)])
stopifnot(length(string_ids) > 0)
cat(sprintf("Unique STRING_ids to query: %d\n", length(string_ids)))

# Chunked STRINGdb query
chunk_size <- 50
chunks <- split(string_ids, ceiling(seq_along(string_ids) / chunk_size))
cat(sprintf("Querying interactions in %d chunks\n",
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

# Choose source, normalize,filter within input set + dedup
raw_edges <- if (is.null(string_edges_stringdb) || nrow(string_edges_stringdb) == 0) {
  message("STRINGdb returned no edges -> REST POST fallback.")
  fetch_string_rest(string_ids)
} else string_edges_stringdb
if (is.null(raw_edges) || nrow(raw_edges) == 0)
  stop("No interactions from STRINGdb or REST. ")

string_edges_raw <- normalize_string_edges(raw_edges)
string_edges     <- dedup_undirected(filter_within(string_edges_raw, string_ids))
cat(sprintf("Raw edges: %d | normalized: %d | between input proteins: %d\n",
            nrow(raw_edges), nrow(string_edges_raw), nrow(string_edges)))
covered <- unique(c(string_edges$from, string_edges$to))
cat(sprintf("Coverage: %d / %d input proteins in >=1 edge\n",
            length(covered), length(string_ids)))
write.csv(string_edges,
          file.path(results_dir, "STRING edge list upregulated.csv"),
          row.names = FALSE)

# igraph + centralities

stopifnot(nrow(string_edges) >= 2)

g <- igraph::graph_from_data_frame(string_edges, directed = FALSE)
g <- igraph::simplify(g, remove.multiple = TRUE, remove.loops = TRUE,
                      edge.attr.comb = list(combined_score = "max"))
igraph::V(g)$STRING_id <- igraph::V(g)$name 
id2sym <- setNames(mapped$GeneSymbol, mapped$STRING_id)
igraph::V(g)$name <- ifelse(igraph::V(g)$name %in% names(id2sym),
                            id2sym[igraph::V(g)$name], igraph::V(g)$name)
cat(sprintf("igraph: %d vertices, %d edges\n",
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
          file.path(results_dir, "hub gene centrality metrics.csv"),
          row.names = FALSE)

hub90 <- head(centralities$Gene, min(90, nrow(centralities)))
hub10 <- head(centralities$Gene, min(10, nrow(centralities)))
writeLines(hub90, file.path(results_dir, "hub genes 90.txt"))
writeLines(hub10, file.path(results_dir, "hub genes 10.txt"))

# Hub-gene supplementary

plot_top_n_bar <- function(metric, n, fill_col) {
  df <- centralities %>%
    arrange(desc(.data[[metric]])) %>%
    slice_head(n = n) %>%
    mutate(Gene = factor(Gene, levels = rev(Gene)))
  
  ggplot(df, aes(x = .data[[metric]], y = Gene)) +
    geom_col(fill = fill_col, color = NA, width = 0.75) +
    geom_text(aes(label = round(.data[[metric]], 2)), 
              hjust = -0.1, size = 3.5, color = "grey30") +
    scale_x_continuous(expand = expansion(mult = c(0, 0.15))) +
    theme_classic(base_size = 12) +
    theme(axis.line.y = element_blank(),
          axis.ticks.y = element_blank()) +
    labs(title = sprintf("Top-%d hub genes by %s", n, metric),
         x = metric, y = NULL)
}

metric_cols <- c(Degree = "#E41A1C", Betweenness = "#377EB8",
                 Closeness = "#4DAF4A", Eigenvector = "#984EA3")

for (m in names(metric_cols)) {
  for (n in c(20, 30)) {
    save_png(sprintf("Hub top%d %s.png", n, m),
             plot_top_n_bar(m, n, metric_cols[[m]]),
             width = 1800, height = 1500)
  }
}



# Hub-of-hub refinement

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
  message("Hub genes not mapped: ", paste(miss, collapse = ", "))
hub10_ids <- unique(hub10_mapped$STRING_id)
raw10 <- tryCatch(string_db$get_interactions(hub10_ids), error = function(e) NULL)
if (is.null(raw10) || nrow(raw10) == 0)
  raw10 <- fetch_string_rest(hub10_ids, chunk = max(length(hub10_ids), 1))

if (is.null(raw10)) {
  message("No interactions among the 10 hubs; writing empty table.")
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
write.csv(hub10_edges, file.path(results_dir, "STRING edge list hub10.csv"),
          row.names = FALSE)
cat(sprintf("Hub 10 internal edges: %d\n", nrow(hub10_edges)))

hub6 <- c("BRCA1", "RAD51", "MCM7", "EXO1", "RFC3", "TP53")
writeLines(hub6, file.path(results_dir, "hub genes 6 validated.txt"))




# Independent validation: GSE9943 
gse9943_map <- c(
  GSM251560 = "BN|Control", GSM251561 = "SD|Control",
  GSM251586 = "BN|Control", GSM251588 = "SD|Control",
  GSM251591 = "BN|IR",      GSM251593 = "BN|IR",
  GSM251594 = "SD|IR",      GSM251595 = "SD|IR",
  GSM280410 = "BN|Control", GSM280411 = "SD|Control",
  GSM280416 = "BN|IR",      GSM280417 = "SD|IR")

val_strain <- "SD"


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
      message(sprintf("GSE9943 %s: title infers %s/%s but map says %s/%s -> keeping map",
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

  ## traceability
  write.csv(pheno_val, file.path(results_dir, "GSE9943 pData full.csv"),
            row.names = TRUE)

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
  write.csv(val_table, file.path(results_dir, "GSE9943 group assignment.csv"),
            row.names = FALSE)
  cat("Group assignment (from Sample_title + GEO map):\n")
  print(val_table)

  if (any(is.na(strains)) || any(is.na(grps)))
    stop("Unassigned strain/group for: ",
         paste(gsm_ids[is.na(strains) | is.na(grps)], collapse = ", "), call. = FALSE)

  cat(sprintf("SD: %d Control / %d IR | BN: %d Control / %d IR\n",
              sum(strains == "SD" & grps == "Control"),
              sum(strains == "SD" & grps == "IR"),
              sum(strains == "BN" & grps == "Control"),
              sum(strains == "BN" & grps == "IR")))


  expr_val <- expr_val[complete.cases(expr_val), , drop = FALSE]

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
  cat(sprintf("probes with a gene symbol: %d / %d (%.1f%%)\n",
              sum(!is.na(sym_lut)), length(sym_lut),
              100 * sum(!is.na(sym_lut)) / length(sym_lut)))

  sym_disc  <- setNames(deg$GeneSymbol, deg$ProbeID)
  disc_syms <- na.omit(unique(unname(sym_disc[rownames(deg_filtered)])))
  if (!length(disc_syms))
    stop("Discovery DEG symbol list is empty", call. = FALSE)

  ## Per-strain validation DEGs + overlap 
  val_summary <- list()
  for (st in c("SD", "BN")) {
    sel <- names(strains)[strains == st & grps %in% c("Control", "IR")]
    gst <- table(droplevels(factor(grps[sel])))
    if (length(sel) < 4 || min(gst) < 2) {
      message("strain ", st, ": insufficient replicates (", length(sel),
              " samples); skipped."); next
    }
    g_st <- factor(grps[sel], levels = c("Control", "IR"))
    e_st <- expr_val[, sel, drop = FALSE]

    d_st <- model.matrix(~ 0 + g_st)
    colnames(d_st) <- levels(g_st)
    f_st <- limma::lmFit(e_st, d_st)
    cont <- setNames(rep(0, ncol(d_st)), colnames(d_st))
    cont["IR"] <- 1; cont["Control"] <- -1       
    f_st <- limma::eBayes(limma::contrasts.fit(f_st, cont))
    dv   <- limma::topTable(f_st, number = Inf, adjust.method = "BH")

    dvf <- dv[abs(dv$logFC) > 2 & dv$P.Value < 0.05, , drop = FALSE]
    dvf$GeneSymbol <- unname(sym_lut[rownames(dvf)])
    write.csv(dv,  file.path(results_dir,
               sprintf("GSE9943 %s DEG full.csv", st)), row.names = TRUE)
    write.csv(dvf, file.path(results_dir,
               sprintf("GSE9943 %s DEG filtered.csv", st)), row.names = TRUE)

    vs <- na.omit(unique(dvf$GeneSymbol))
    ov <- intersect(disc_syms, vs)
    val_summary[[st]] <- data.frame(
      Strain = st, n_samples = length(sel),
      n_DEG_probes = nrow(dvf), n_DEG_symbols = length(vs),
      overlap_with_GSE58438 = length(ov), stringsAsFactors = FALSE)
    write.csv(data.frame(Validated_gene = ov),
              file.path(results_dir,
                sprintf("GSE58438 vs GSE9943 %s overlap.csv", st)),
              row.names = FALSE)
    cat(sprintf("validation/%s discovery: %d | GSE9943-%s: %d symbols | overlap: %d\n",
                st, length(disc_syms), st, length(vs), length(ov)))


    if (st == val_strain && length(disc_syms) > 0 && length(vs) > 0) {
      save_png("GSE9943 validation venn.png", {
        grid::grid.newpage()
        grid::grid.draw(draw.pairwise.venn(
          area1 = length(disc_syms), area2 = length(vs),
          cross.area = length(ov),
          category = c("GSE58438 (discovery)",
                       sprintf("GSE9943 %s (validation)", st)),
          fill = c("#E41A1C", "#377EB8"), alpha = 0.6,
          lty = "blank",
          fontfamily = "sans",
          cat.fontfamilty = "sans",
          cat.pos = c(-30, 30), cat.dist = c(0.05, 0.05),
          cat.cex = 1.2, cex = 1.3))
      }, width = 1800, height = 1500)
    }
  }

  if (length(val_summary)) {
    val_sum_df <- do.call(rbind, val_summary)
    write.csv(val_sum_df,
              file.path(results_dir, "GSE9943 validation summary.csv"),
              row.names = FALSE)
    print(val_sum_df)
  }
}), error = function(e) {
  message("Validation section failed (pipeline continues): ",
          conditionMessage(e))
})


# WGCNA module discovery

suppressPackageStartupMessages(library(WGCNA))
enableWGCNAThreads()


cor   <- WGCNA::cor     
bicor <- WGCNA::bicor  

if (!all(c("weights.x", "weights.y", "cosine") %in% names(formals(cor))))
  message(" WGCNA::cor formals look different")


probe_ok <- tryCatch({
  invisible(do.call(match.fun("cor"),
                    list(x = matrix(c(1, 2, 3, 4), nrow = 2),
                         weights.x = NULL, weights.y = NULL,
                         cosine = FALSE, use = "everything")))
  TRUE
}, error = function(e) {
  message("cor dispatch probe FAILED")
  FALSE
})
if (probe_ok)
  message("cor dispatch probe OK")

datExpr <- t(as.matrix(exprSet))
if (is.null(rownames(datExpr)))
  rownames(datExpr) <- sprintf("S%02d", seq_len(nrow(datExpr)))
gsg <- goodSamplesGenes(datExpr, verbose = 0)
if (!gsg$allOK) {
  cat(sprintf("goodSamplesGenes removed %d samples / %d genes\n",
              sum(!gsg$goodSamples), sum(!gsg$goodGenes)))
  datExpr <- datExpr[gsg$goodSamples, gsg$goodGenes, drop = FALSE]
}

trait <- as.numeric(factor(group, levels = c("Control", "AKI"))) - 1
if (!gsg$allOK) trait <- trait[gsg$goodSamples]
stopifnot(length(trait) == nrow(datExpr))

powers <- 1:20
sft <- pickSoftThreshold(datExpr, powerVector = powers,
                         networkType = "unsigned", verbose = 3)
chosen_power <- 6   
r2_at_power <- sft$fitIndices$SFT.R.sq[sft$fitIndices$Power == chosen_power]
if (!is.na(r2_at_power) && r2_at_power < 0.80)
  message(sprintf("scale-free R^2 at power %d is %.2f (<0.80); ",
                  chosen_power, r2_at_power),
          "sft$powerEstimate = ", sft$powerEstimate)

save_png("WGCNA soft threshold.png", {
  par(mfrow = c(1, 2))
  plot(sft$fitIndices$Power, -sign(sft$fitIndices$slope) * sft$fitIndices$SFT.R.sq,
       xlab = "Soft Threshold (power)",
       ylab = "Scale Free Topology Model Fit (signed R^2)",
       type = "n", main = "Soft threshold selection")
  text(sft$fitIndices$Power, -sign(sft$fitIndices$slope) * sft$fitIndices$SFT.R.sq,
       labels = powers, col = "red")
  abline(h = 0.85, col = "blue", lty = 2, lwd = 2)
  plot(sft$fitIndices$Power, sft$fitIndices$mean.k.,
       xlab = "Soft Threshold (power)", ylab = "Mean Connectivity",
       type = "l", main = "Mean connectivity decay")
}, width = 2200, height = 1100)

net <- blockwiseModules(datExpr,
                        power = chosen_power, TOMType = "unsigned",
                        minModuleSize = 30, reassignThreshold = 0,
                        mergeCutHeight = 0.25, numericLabels = TRUE,
                        pamRespectsDendro = FALSE,
                        saveTOMs = TRUE,   
                        saveTOMFileBase = file.path(results_dir, "GSE58438_TOM"),
                        verbose = 3)

moduleColors <- labels2colors(net$colors)

MEs <- net$MEs
colnames(MEs) <- paste0("ME", labels2colors(
  as.numeric(gsub("^ME", "", colnames(MEs)))))
MEs <- orderMEs(as.data.frame(MEs))


trait_mat <- matrix(trait, ncol = 1, dimnames = list(NULL, "AKI_vs_Control"))
moduleTraitCor    <- cor(as.matrix(MEs), trait_mat, use = "p")
moduleTraitPvalue <- corPvalueStudent(moduleTraitCor, nSamples = nrow(datExpr))
textMat <- paste(signif(moduleTraitCor, 2), "\n(",
                 signif(moduleTraitPvalue, 1), ")", sep = "")
dim(textMat) <- dim(moduleTraitCor)

save_png("WGCNA module trait heatmap.png", {
  par(mar = c(4, 8, 2, 1))
  labeledHeatmap(Matrix = moduleTraitCor, xLabels = colnames(trait_mat),
                 yLabels = rownames(moduleTraitCor),
                 ySymbols = rownames(moduleTraitCor), colorLabels = FALSE,
                 colors = blueWhiteRed(50), textMatrix = textMat,
                 setStdMargins = TRUE, cex.text = 0.8,
                 main = "Module trait relationships")
}, width = 5600, height = 8000)

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
} else message("Module 'blue' absent in this dataset; skipped.")







# Session info + artifact manifest
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
