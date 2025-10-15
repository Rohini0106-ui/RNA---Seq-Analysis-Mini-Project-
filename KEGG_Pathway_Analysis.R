# Clean, single-file script for KEGG + GO enrichment using Entrez IDs

# Set working directory 
setwd("C:/Users/LENOVO/Downloads") 

# Save this as run_enrichment.R and run in your working directory.

# ---- Setup ----
# (Install only if missing; comment out the install lines after first successful run)
if (!requireNamespace("BiocManager", quietly = TRUE)) install.packages("BiocManager")
if (!requireNamespace("clusterProfiler", quietly = TRUE)) BiocManager::install("clusterProfiler")
if (!requireNamespace("org.Hs.eg.db", quietly = TRUE)) BiocManager::install("org.Hs.eg.db")
if (!requireNamespace("pathview", quietly = TRUE)) BiocManager::install("pathview")
if (!requireNamespace("enrichplot", quietly = TRUE)) BiocManager::install("enrichplot")

# Load libraries
library(clusterProfiler)
library(org.Hs.eg.db)
library(pathview)
library(enrichplot)
library(ggplot2)

# ---- Parameters ----
input_file <- "Annotated_Significant_DEGs_WithEnrichrHits.csv"  # ensure this file exists in WD
entrez_col <- "entrezgene_id"    # adjust if your file uses a different column name
lfc_col <- "log2FoldChange"
outdir <- "enrichment_output"
plots_dir <- file.path(outdir, "plots")

dir.create(outdir, showWarnings = FALSE, recursive = TRUE)
dir.create(plots_dir, showWarnings = FALSE, recursive = TRUE)

# ---- Read & validate input ----
if (!file.exists(input_file)) stop("Input file not found: ", input_file)
deg_data <- read.csv(input_file, stringsAsFactors = FALSE, check.names = FALSE)

# Check columns
if (!(entrez_col %in% colnames(deg_data))) {
  stop("Entrez column '", entrez_col, "' not found. Available columns: ", paste(colnames(deg_data), collapse = ", "))
}
if (!(lfc_col %in% colnames(deg_data))) {
  stop("log2FoldChange column '", lfc_col, "' not found. Available columns: ", paste(colnames(deg_data), collapse = ", "))
}

# Filter out rows with missing IDs or LFC
deg_data_filtered <- deg_data[!is.na(deg_data[[entrez_col]]) & deg_data[[entrez_col]] != "" &
                                !is.na(deg_data[[lfc_col]]) & deg_data[[lfc_col]] != "", ]

if (nrow(deg_data_filtered) == 0) stop("No rows left after filtering out missing entrez IDs / log2FoldChange.")

# Make sure Entrez IDs are character (no factors or decimals)
# Many pipelines write Entrez IDs as numbers — convert to character without scientific notation
deg_data_filtered[[entrez_col]] <- as.character(as.integer(as.numeric(deg_data_filtered[[entrez_col]])))

# Remove duplicates (keep the one with largest absolute LFC)
deg_data_filtered <- deg_data_filtered[order(-abs(as.numeric(deg_data_filtered[[lfc_col]]))), ]
deg_data_filtered <- deg_data_filtered[!duplicated(deg_data_filtered[[entrez_col]]), ]

# Create named numeric vector geneList (values = log2FC, names = Entrez)
geneList <- as.numeric(deg_data_filtered[[lfc_col]])
names(geneList) <- deg_data_filtered[[entrez_col]]

# Basic checks
cat("Number of genes in geneList:", length(geneList), "\n")
cat("First 6 entries (name = Entrez ID -> value = log2FC):\n")
print(head(geneList))

# ---- KEGG enrichment ----
# clusterProfiler's enrichKEGG expects Entrez IDs. Use names(geneList) directly.
kegg_enrich <- tryCatch({
  enrichKEGG(gene = names(geneList), organism = "hsa", pvalueCutoff = 0.05)
}, error = function(e) {
  message("KEGG enrichment failed: ", e$message)
  NULL
})

if (!is.null(kegg_enrich) && length(kegg_enrich) > 0) {
  kegg_df <- as.data.frame(kegg_enrich)
  write.csv(kegg_df, file = file.path(outdir, "KEGG_Enrichment_Results.csv"), row.names = FALSE)
  png(file.path(plots_dir, "KEGG_dotplot.png"), width = 1000, height = 800)
  print(dotplot(kegg_enrich, showCategory = 15) + ggtitle("KEGG Pathway Enrichment"))
  dev.off()
  
  # generate pathway diagrams for top 3 pathways (if any)
  top_pathways <- head(kegg_df$ID, 3)
  for (pid in top_pathways) {
    # pathview will create files in working dir; set prefix to keep organized
    tryCatch({
      pathview(gene.data = geneList, pathway.id = pid, species = "hsa", out.suffix = paste0("kegg_", pid))
    }, error = function(e) message("pathview failed for ", pid, ": ", e$message))
  }
} else {
  message("No KEGG results to save/plot.")
}

# ---- GO enrichment (BP, MF, CC) ----
# clusterProfiler::enrichGO expects a vector of gene IDs (ENTREZID) OR use bitr conversion if you have other IDs.
go_results <- list()
ontologies <- c(BP = "BP", MF = "MF", CC = "CC")

for (ont in ontologies) {
  message("Running GO ", ont, " ...")
  res <- tryCatch({
    enrichGO(gene = names(geneList),
             OrgDb = org.Hs.eg.db,
             keyType = "ENTREZID",
             ont = ont,
             pAdjustMethod = "BH",
             pvalueCutoff = 0.05,
             readable = TRUE)
  }, error = function(e) {
    message("enrichGO error for ", ont, ": ", e$message)
    NULL
  })
  go_results[[ont]] <- res
  if (!is.null(res) && length(res) > 0) {
    write.csv(as.data.frame(res), file = file.path(outdir, paste0("GO_", ont, "_enrichment.csv")), row.names = FALSE)
    png(file.path(plots_dir, paste0("GO_", ont, "_dotplot.png")), width = 1000, height = 800)
    print(dotplot(res, showCategory = 20) + ggtitle(paste("GO", ont, "Enrichment")))
    dev.off()
    png(file.path(plots_dir, paste0("GO_", ont, "_barplot.png")), width = 1000, height = 800)
    print(barplot(res, showCategory = 10, title = paste("GO", ont)))
    dev.off()
  } else {
    message("No GO results for ontology: ", ont)
  }
}

# ---- Combined GO three-ontology barplot (top terms) ----
# Build a combined dataframe of top terms (if available)
combine_top_terms <- function(res_obj, ont_label) {
  if (is.null(res_obj) || length(res_obj) == 0) return(NULL)
  df <- as.data.frame(res_obj)
  topn <- min(10, nrow(df))
  df <- df[1:topn, c("Description", "Count", "p.adjust")]
  df$Ontology <- ont_label
  colnames(df)[1] <- "Term"
  df
}

bp_df <- combine_top_terms(go_results[["BP"]], "BP")
cc_df <- combine_top_terms(go_results[["CC"]], "CC")
mf_df <- combine_top_terms(go_results[["MF"]], "MF")
combined_df <- do.call(rbind, list(bp_df, cc_df, mf_df))
if (!is.null(combined_df) && nrow(combined_df) > 0) {
  png(file.path(plots_dir, "GO_Three_Ontologies.png"), width = 1200, height = 900)
  p <- ggplot(combined_df, aes(x = Count, y = reorder(Term, Count), fill = Ontology)) +
    geom_bar(stat = "identity", position = "dodge") +
    labs(title = "GO Results: Three Ontologies", x = "Gene Count", y = "Term") +
    theme_bw() +
    facet_wrap(~Ontology, scales = "free_y", ncol = 1)
  print(p)
  dev.off()
}

# ---- Save session info for reproducibility ----
writeLines(capture.output(sessionInfo()), file.path(outdir, "sessionInfo.txt"))

cat("Enrichment run complete. Results and plots saved in:", outdir, "\n")
