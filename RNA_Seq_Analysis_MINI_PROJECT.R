#set working directory 
setwd("C:/Users/LENOVO/Downloads") 

# Load libraries
library(dplyr)
library(DESeq2)
library(tidyverse)
library(GEOquery)
library(EnhancedVolcano)
library(pheatmap)
library(RColorBrewer)
library(clusterProfiler)
library(org.Hs.eg.db)
library(AnnotationDbi)
library(enrichplot)
library(ggplot2)
library(biomaRt)

list.files()

# Step 1: Load raw counts & metadata
counts <- read.table("Merged_Feature_CountsFinal.txt", header=TRUE, row.names=1, sep="\t", check.names=FALSE)
metadata <- read.csv("metadata.csv", header=TRUE, row.names=1)
all(colnames(counts) %in% rownames(metadata))  # Should return TRUE

# Step 2: Create DESeq2 data set & pre-filter
dds <- DESeqDataSetFromMatrix(countData=counts,
                              colData=metadata,
                              design=~condition)

# Step 3: Run DESeq2 analysis
dds <- DESeq(dds)
res <- results(dds, contrast=c("condition", "T2D", "Non_T2D"))
head(res)
write.csv(as.data.frame(res), "DESeq2_results.csv")

# Step 4: MA plot & volcano plot
plotMA(res)

EnhancedVolcano(res,
                lab = rownames(res),
                x = 'log2FoldChange',
                y = 'padj',
                pCutoff = 0.05,
                FCcutoff = 1,
                title = 'Volcano Plot of DEGs')
# Order results by adjusted p-value (padj)
res_sorted <- res[order(res$padj), ]

# Step 5: Filter significant genes
res_df <- as.data.frame(res_sorted)
sig_genes <- subset(res_df, padj < 0.05 & abs(log2FoldChange) > 1)
write.csv(sig_genes, "Significant_DEGs.csv", row.names = TRUE)

# Step 6: Separate upregulated and downregulated genes
upregulated_genes <- subset(res_df, padj < 0.05 & log2FoldChange > 0)
write.csv(upregulated_genes, "Upregulated_genes.csv", row.names = TRUE)

downregulated_genes <- subset(res_df, padj < 0.05 & log2FoldChange < 0)
write.csv(downregulated_genes, "Downregulated_genes.csv", row.names = TRUE)

# Step 7: Normalized counts heatmap for significant genes
vsd <- vst(dds, blind = FALSE)
norm_counts <- assay(vsd)
sig_gene_names <- rownames(sig_genes)
sig_norm_counts <- norm_counts[sig_gene_names,]

annotation_col <- data.frame(Type = factor(meta_data$condition))
rownames(annotation_col) <- colnames(sig_norm_counts)

# Example: if  normalized counts have 6 samples (3 control, 3 T2D)
sample_info <- data.frame(
  condition = c(rep("Control", 3), rep("T2D", 3))
)
rownames(sample_info) <- colnames(sig_norm_counts)

# Assign it as annotation_col
annotation_col <- sample_info

# Now plot
pheatmap(sig_norm_counts,
         cluster_rows = TRUE,
         cluster_cols = TRUE,
         annotation_col = annotation_col,
         show_rownames = TRUE,
         show_colnames = FALSE,
         color = colorRampPalette(rev(brewer.pal(9, "RdBu")))(255),
         main = "Heatmap of Significant DEGs")

# Step 8: PCA plot for sample clustering
pca <- prcomp(t(norm_counts))
percent_var <- (pca$sdev^2) / sum(pca$sdev^2) * 100
pca_data <- as.data.frame(pca$x)
pca_data$Type <- metadata$condition  # Use "condition" here, not "Type"

ggplot(pca_data, aes(PC1, PC2, color = Type)) +
  geom_point(size = 3) +
  labs(
    title = "PCA of RNA-seq samples",
    x = sprintf("PC1 (%.2f%%)", percent_var[1]),
    y = sprintf("PC2 (%.2f%%)", percent_var[2])
  ) +
  theme_minimal()

# Step 9: Sample distance heatmap
rld <- rlog(dds, blind = TRUE)
sampleDists <- dist(t(assay(rld)))
sampleDistMatrix <- as.matrix(sampleDists)
colors <- colorRampPalette(rev(brewer.pal(9, "Blues")))(255)

pheatmap(sampleDistMatrix,
         clustering_distance_rows = sampleDists,
         clustering_distance_cols = sampleDists,
         col = colors,
         main = "Sample-to-Sample Distance Heatmap")

# Step 10: Dispersion plot
plotDispEsts(dds, main = "Dispersion Plot")

# Step 11: Gene annotation using biomaRt
ensembl <- useMart("ensembl", dataset = "hsapiens_gene_ensembl")
annot <- getBM(attributes = c('ensembl_gene_id', 'entrezgene_id', 'hgnc_symbol', 'description'),
               filters = 'entrezgene_id',
               values = sig_gene_names,
               mart = ensembl)
sig_genes_df <- as.data.frame(sig_genes)
sig_genes_df <- rownames_to_column(sig_genes_df, var = "entrezgene_id")
annot$entrezgene_id <- as.character(annot$entrezgene_id)
sig_genes_annot <- left_join(sig_genes_df, annot, by = "entrezgene_id")
head(sig_genes_annot)
write.csv(sig_genes_annot, "Annotated_Significant_DEGs.csv", row.names = FALSE)

library(clusterProfiler)
library(org.Hs.eg.db)
library(enrichplot)

# Extract valid Entrez IDs (remove NAs and duplicates)
entrez_ids <- unique(na.omit(sig_genes_annot$entrezgene_id))

# Ensure they are character strings
entrez_ids <- as.character(entrez_ids)

# --- GO Enrichment ---
ego_bp <- enrichGO(
  gene = entrez_ids,
  OrgDb = org.Hs.eg.db,
  keyType = "ENTREZID",
  ont = "BP",
  pAdjustMethod = "BH",
  pvalueCutoff = 0.05,
  qvalueCutoff = 0.05,
  readable = TRUE
)

ego_mf <- enrichGO(
  gene = entrez_ids,
  OrgDb = org.Hs.eg.db,
  keyType = "ENTREZID",
  ont = "MF",
  pAdjustMethod = "BH",
  pvalueCutoff = 0.05,
  qvalueCutoff = 0.05,
  readable = TRUE
)

ego_cc <- enrichGO(
  gene = entrez_ids,
  OrgDb = org.Hs.eg.db,
  keyType = "ENTREZID",
  ont = "CC",
  pAdjustMethod = "BH",
  pvalueCutoff = 0.05,
  qvalueCutoff = 0.05,
  readable = TRUE
)

# --- KEGG Enrichment ---
kegg_enrich <- enrichKEGG(
  gene = entrez_ids,
  organism = "hsa",
  pvalueCutoff = 0.05
)

# Step 13: Visualization of enrichment results
dotplot(ego_bp, showCategory = 20) + ggtitle("GO: Biological Process")
dotplot(ego_mf, showCategory = 20) + ggtitle("GO: Molecular Function")
dotplot(ego_cc, showCategory = 20) + ggtitle("GO: Cellular Component")
dotplot(kegg_enrich, showCategory = 20) + ggtitle("KEGG Pathway Enrichment")

getwd()
list.files()

# Step 14: Load Enrichr clinical overlap results (manually exported CSV)
enrichr_table <- read.csv("Annotated_Significant_DEGs.csv")

colnames(enrichr_table)

# 1 Create a hit_count column: 1 if padj < 0.05, else 0
enrichr_table$hit_count <- ifelse(!is.na(enrichr_table$padj) & enrichr_table$padj < 0.05, 1, 0)

# 2 Replace any NAs in hit_count with 0 (just to be safe)
enrichr_table$hit_count[is.na(enrichr_table$hit_count)] <- 0

# 3️ Sort the table by hit_count decreasing
enrichr_table <- enrichr_table[order(-enrichr_table$hit_count), ]

# 4️ Optional: view the top rows
head(enrichr_table)

# Replace NA with 0 and sort decreasing by overlap count
# Replace NA in hit_count with 0 (just in case)
enrichr_table$hit_count[is.na(enrichr_table$hit_count)] <- 0

# Sort decreasing by hit_count
enrichr_table <- enrichr_table[order(-enrichr_table$hit_count), ]

# View top rows
head(enrichr_table)

# Save sorted annotation
write.csv(enrichr_table, "Final_Annotated_DEGs_With_DiabetesHits_Sorted.csv", row.names=FALSE)

# Step 15: Biomarker prioritization - rank by overlap and fold change magnitude
enrichr_table$biomarker_rank_score <- enrichr_table$hit_count * abs(enrichr_table$log2FoldChange)

# Sort decreasing by biomarker_rank_score
enrichr_table <- enrichr_table[order(-enrichr_table$biomarker_rank_score), ]

# View top-ranked biomarkers
head(enrichr_table)

library(ggplot2)

# 16 Select top 10 biomarkers by hit_count or biomarker_rank_score
top_biomarkers <- head(enrichr_table[order(-enrichr_table$biomarker_rank_score), ], 10)

# Create the barplot
ggplot(top_biomarkers, aes(x = reorder(hgnc_symbol, hit_count), y = hit_count)) +
  geom_bar(stat = "identity", fill = "steelblue") +
  coord_flip() +
  labs(
    title = "Top 10 Biomarker Candidates by Diabetes Study Hit Count",
    x = "Gene Symbol",
    y = "Number of Diabetes GEO Overlaps"
  ) +
  theme_minimal()

library(dplyr)

#17 Select top 10 genes by biomarker_rank_score
top10_genes <- head(enrichr_table[order(-enrichr_table$biomarker_rank_score), ], 10)

# Create summary table with selected columns
summary_table <- dplyr::select(top10_genes, hgnc_symbol, log2FoldChange, padj, hit_count, description)

# Print summary table to console
print(summary_table)

# Export summary table to CSV
write.csv(summary_table, "Top_10_Biomarker_Genes_Summary.csv", row.names = FALSE)
