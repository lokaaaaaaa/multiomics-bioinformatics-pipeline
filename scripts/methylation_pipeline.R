# ============================================================
# DNA Methylation Analysis Pipeline
# Tools: minfi, ChAMP, limma | Array: Illumina 450K / EPIC
# ============================================================

if (!requireNamespace("BiocManager", quietly = TRUE)) install.packages("BiocManager")
BiocManager::install(c("minfi", "ChAMP", "limma", "IlluminaHumanMethylation450kanno.ilmn12.hg19",
                        "IlluminaHumanMethylation450kmanifest", "missMethyl",
                        "DMRcate", "bumphunter"), ask = FALSE)
install.packages(c("tidyverse", "pheatmap", "RColorBrewer", "ggplot2"), repos = "https://cloud.r-project.org")

library(minfi)
library(limma)
library(ChAMP)
library(missMethyl)
library(DMRcate)
library(tidyverse)
library(pheatmap)
library(RColorBrewer)
library(IlluminaHumanMethylation450kanno.ilmn12.hg19)

dir.create("results", showWarnings = FALSE)
dir.create("plots",   showWarnings = FALSE)

# ── 1. Load IDAT Files (or simulate beta matrix) ───────────────
idat_dir <- "data/idat"

if (dir.exists(idat_dir) && length(list.files(idat_dir, pattern = ".idat")) > 0) {
  message("Loading IDAT files from: ", idat_dir)
  targets <- read.metharray.sheet(idat_dir)
  rgSet   <- read.metharray.exp(targets = targets)
} else {
  message("No IDAT files found — simulating beta matrix for demonstration...")
  set.seed(123)
  n_cpgs   <- 5000
  n_ctrl   <- 3
  n_treat  <- 3
  n_total  <- n_ctrl + n_treat

  # Simulate beta values (0–1 bounded)
  beta_sim <- matrix(
    rbeta(n_cpgs * n_total, shape1 = 2, shape2 = 5),
    nrow = n_cpgs,
    dimnames = list(
      paste0("cg", formatC(seq_len(n_cpgs), width = 8, flag = "0")),
      c(paste0("Ctrl_", seq_len(n_ctrl)), paste0("Treat_", seq_len(n_treat)))
    )
  )

  # Add ~500 truly differentially methylated CpGs
  dm_idx <- sample(n_cpgs, 500)
  beta_sim[dm_idx, (n_ctrl + 1):n_total] <-
    pmin(pmax(beta_sim[dm_idx, (n_ctrl + 1):n_total] +
              runif(length(dm_idx) * n_treat, 0.15, 0.4), 0), 1)

  beta_matrix <- beta_sim
  pd <- data.frame(
    Sample_Name = colnames(beta_matrix),
    Sample_Group = factor(c(rep("Control", n_ctrl), rep("Treatment", n_treat))),
    row.names    = colnames(beta_matrix)
  )
}

# ── 2. Quality Control (with real data only) ──────────────────
# If using minfi pipeline:
# MSet        <- preprocessRaw(rgSet)
# qcReport(rgSet, sampNames = targets$Sample_Name, pdf = "results/QCReport.pdf")
# detP        <- detectionP(rgSet)

# ── 3. Normalisation ──────────────────────────────────────────
if (exists("rgSet")) {
  # Functional normalisation (recommended for differential methylation)
  GRset    <- preprocessFunnorm(rgSet)
  beta_matrix <- getBeta(GRset)
  pd       <- pData(GRset)
  pd$Sample_Group <- factor(pd$Sample_Group)
}

cat("Beta matrix: ", nrow(beta_matrix), "CpGs ×", ncol(beta_matrix), "samples\n")

# ── 4. Filter CpGs ────────────────────────────────────────────
# Remove CpGs with low variance (uninformative)
cpg_var  <- apply(beta_matrix, 1, var)
beta_filt <- beta_matrix[cpg_var > quantile(cpg_var, 0.2), ]
cat("After variance filtering:", nrow(beta_filt), "CpGs remain\n")

# M-value transformation for statistical testing
M_vals <- log2(beta_filt / (1 - beta_filt) + 1e-6)

# ── 5. Differential Methylation — limma ───────────────────────
group    <- pd$Sample_Group
design   <- model.matrix(~ group)
fit      <- lmFit(M_vals, design)
fit2     <- eBayes(fit)
dmp      <- topTable(fit2, coef = 2, number = Inf, sort.by = "p")
dmp$delta_beta <- rowMeans(beta_filt[rownames(dmp), group == levels(group)[2]]) -
                  rowMeans(beta_filt[rownames(dmp), group == levels(group)[1]])

sig_dmps <- dmp |> filter(adj.P.Val < 0.05, abs(delta_beta) > 0.1)
cat("Significant DMPs:", nrow(sig_dmps), "\n")

write.csv(dmp,      "results/all_DMPs.csv",         row.names = TRUE)
write.csv(sig_dmps, "results/significant_DMPs.csv", row.names = TRUE)

# ── 6. DMR Analysis with bumphunter ───────────────────────────
# (Requires genomic coordinates from annotation)
# For real 450K/EPIC data:
# cpg_annot <- getAnnotation(IlluminaHumanMethylation450kanno.ilmn12.hg19)
# chr <- cpg_annot[rownames(M_vals), "chr"]
# pos <- cpg_annot[rownames(M_vals), "pos"]
# designB <- model.matrix(~ group)
# dmrs <- bumphunter(M_vals, design = designB, chr = chr, pos = pos,
#                    B = 1000, type = "Beta", cutoff = 0.2)

# ── 7. Visualisations ─────────────────────────────────────────
## 7a. Density plot of beta values
png("plots/beta_density.png", width = 800, height = 600)
d_ctrl  <- density(rowMeans(beta_filt[, group == levels(group)[1]]))
d_treat <- density(rowMeans(beta_filt[, group == levels(group)[2]]))
plot(d_ctrl,  col = "#1d3557", lwd = 2, main = "Beta Value Distribution",
     xlab = "Beta value", xlim = c(0, 1))
lines(d_treat, col = "#e63946", lwd = 2)
legend("topright", legend = levels(group), col = c("#1d3557", "#e63946"), lwd = 2)
dev.off()

## 7b. Volcano plot for DMPs
dmp_plot <- dmp |>
  mutate(sig = ifelse(adj.P.Val < 0.05 & abs(delta_beta) > 0.1,
                      ifelse(delta_beta > 0, "Hyper", "Hypo"), "NS"))

png("plots/methylation_volcano.png", width = 800, height = 600)
ggplot(dmp_plot, aes(delta_beta, -log10(P.Value), colour = sig)) +
  geom_point(size = 1.2, alpha = 0.6) +
  scale_colour_manual(values = c(Hyper = "#e63946", Hypo = "#457b9d", NS = "#adb5bd")) +
  geom_vline(xintercept = c(-0.1, 0.1), lty = 2, colour = "black") +
  geom_hline(yintercept = -log10(0.05), lty = 2, colour = "black") +
  labs(title = "Differential Methylation Volcano Plot",
       x = "Δβ (Treatment − Control)",
       y = "-log₁₀ p-value",
       colour = "Status") +
  theme_bw(base_size = 13)
dev.off()

## 7c. Heatmap of top 100 DMPs
top100 <- head(rownames(sig_dmps), 100)
mat    <- beta_filt[top100, ]
annot  <- data.frame(Group = group, row.names = colnames(mat))
cols   <- list(Group = c(Control = "#457b9d", Treatment = "#e63946"))

png("plots/methylation_heatmap.png", width = 900, height = 1200)
pheatmap(mat,
  annotation_col  = annot,
  annotation_colors = cols,
  color           = colorRampPalette(c("#1d3557", "white", "#e63946"))(100),
  show_rownames   = FALSE,
  main            = "Top 100 Differentially Methylated CpGs",
  clustering_distance_rows = "euclidean",
  clustering_distance_cols = "euclidean")
dev.off()

## 7d. MDS plot
png("plots/methylation_MDS.png", width = 700, height = 600)
plotMDS(M_vals,
        col    = ifelse(group == levels(group)[1], "#1d3557", "#e63946"),
        pch    = 16,
        cex    = 1.5,
        main   = "MDS Plot — M-values")
legend("topright", legend = levels(group),
       col = c("#1d3557", "#e63946"), pch = 16)
dev.off()

cat("\n✅ Methylation pipeline complete. See results/ and plots/.\n")
