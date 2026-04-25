# ============================================================
# RNA-seq Differential Expression Analysis Pipeline
# Tool: DESeq2 | Author: Your Name | Date: 2026
# ============================================================

# --- 1. Install & Load Libraries ---
if (!requireNamespace("BiocManager", quietly = TRUE)) install.packages("BiocManager")
BiocManager::install(c("DESeq2", "AnnotationDbi", "org.Hs.eg.db",
                        "EnhancedVolcano", "clusterProfiler", "biomaRt"), ask = FALSE)
install.packages(c("tidyverse", "pheatmap", "RColorBrewer", "ggplot2"), repos = "https://cloud.r-project.org")

library(DESeq2)
library(tidyverse)
library(pheatmap)
library(RColorBrewer)
library(EnhancedVolcano)
library(clusterProfiler)
library(org.Hs.eg.db)
library(AnnotationDbi)

# --- 2. Load Count Matrix ---
# Replace with your actual count matrix path
# Rows = genes, Columns = samples
counts_path <- "data/raw_counts.csv"   # e.g., from featureCounts or STAR

if (!file.exists(counts_path)) {
  message("Simulating example count data...")
  set.seed(42)
  n_genes   <- 2000
  n_samples <- 6
  count_matrix <- matrix(
    rnbinom(n_genes * n_samples, mu = 500, size = 1),
    nrow = n_genes,
    dimnames = list(
      paste0("GENE_", seq_len(n_genes)),
      paste0("Sample_", seq_len(n_samples))
    )
  )
  dir.create("data", showWarnings = FALSE)
  write.csv(count_matrix, counts_path, row.names = TRUE)
}

count_matrix <- read.csv(counts_path, row.names = 1) |> as.matrix()
cat("Count matrix dimensions:", dim(count_matrix), "\n")

# --- 3. Sample Metadata ---
col_data <- data.frame(
  sample    = colnames(count_matrix),
  condition = factor(rep(c("Control", "Treatment"), each = ncol(count_matrix) / 2)),
  row.names = colnames(count_matrix)
)
cat("Sample metadata:\n")
print(col_data)

# --- 4. Build DESeq2 Object ---
dds <- DESeqDataSetFromMatrix(
  countData = count_matrix,
  colData   = col_data,
  design    = ~ condition
)

# --- 5. Pre-filtering: remove low-count genes ---
keep <- rowSums(counts(dds) >= 10) >= 3
dds  <- dds[keep, ]
cat("Genes after filtering:", nrow(dds), "\n")

# --- 6. Run DESeq2 ---
dds    <- DESeq(dds)
res    <- results(dds, contrast = c("condition", "Treatment", "Control"),
                  alpha = 0.05)
res_df <- as.data.frame(res) |> rownames_to_column("gene_id")
cat("DESeq2 summary:\n")
summary(res)

# --- 7. Significant Genes ---
sig_genes <- res_df |>
  filter(!is.na(padj), padj < 0.05, abs(log2FoldChange) > 1) |>
  arrange(padj)
cat("Significant DEGs:", nrow(sig_genes), "\n")

dir.create("results", showWarnings = FALSE)
write.csv(sig_genes, "results/significant_DEGs.csv", row.names = FALSE)
write.csv(res_df,    "results/all_DEGs.csv",          row.names = FALSE)

# --- 8. Visualisations ---
dir.create("plots", showWarnings = FALSE)

## 8a. PCA Plot
vsd <- vst(dds, blind = FALSE)
pca_data <- plotPCA(vsd, intgroup = "condition", returnData = TRUE)
pct_var  <- round(100 * attr(pca_data, "percentVar"))
png("plots/PCA_plot.png", width = 800, height = 600)
ggplot(pca_data, aes(PC1, PC2, color = condition, label = name)) +
  geom_point(size = 4, alpha = 0.85) +
  geom_text(vjust = -0.8, size = 3) +
  xlab(paste0("PC1: ", pct_var[1], "% variance")) +
  ylab(paste0("PC2: ", pct_var[2], "% variance")) +
  theme_bw(base_size = 14) +
  ggtitle("PCA — Variance Stabilised Counts") +
  scale_color_brewer(palette = "Set1")
dev.off()

## 8b. Volcano Plot
png("plots/Volcano_plot.png", width = 900, height = 700)
EnhancedVolcano(res_df,
  lab       = res_df$gene_id,
  x         = "log2FoldChange",
  y         = "padj",
  title     = "Volcano Plot: Treatment vs Control",
  pCutoff   = 0.05,
  FCcutoff  = 1,
  pointSize = 2.5,
  labSize   = 3)
dev.off()

## 8c. Heatmap of top 50 DEGs
top50 <- head(sig_genes$gene_id, 50)
mat   <- assay(vsd)[top50, ]
mat   <- mat - rowMeans(mat)
png("plots/Heatmap_top50_DEGs.png", width = 900, height = 1000)
pheatmap(mat,
  annotation_col = col_data["condition"],
  color          = colorRampPalette(rev(brewer.pal(11, "RdBu")))(100),
  main           = "Top 50 DEGs — Z-score",
  fontsize_row   = 7,
  show_rownames  = TRUE)
dev.off()

## 8d. MA Plot
png("plots/MA_plot.png", width = 800, height = 600)
plotMA(res, main = "MA Plot — DESeq2 Results", ylim = c(-5, 5), alpha = 0.05)
dev.off()

# --- 9. Gene Ontology Enrichment ---
gene_list <- sig_genes |>
  filter(log2FoldChange > 1) |>
  pull(gene_id)

# Map gene IDs to Entrez (works best with ENSEMBL or SYMBOL IDs)
# Here we skip mapping since genes are simulated
# In real data: use bitr() from clusterProfiler
# ego <- enrichGO(gene = entrez_ids, OrgDb = org.Hs.eg.db, ont = "BP", readable = TRUE)
# dotplot(ego)

cat("\n✅ Pipeline complete. Check results/ and plots/ directories.\n")
