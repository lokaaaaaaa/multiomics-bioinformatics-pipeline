# Multi-Omics Analysis Pipeline

> **Reproducible, open-source pipelines for transcriptomics, proteomics, methylation analysis, and AI/ML-driven omics integration**

[![License: MIT](https://img.shields.io/badge/License-MIT-yellow.svg)](LICENSE)
[![Python 3.10+](https://img.shields.io/badge/Python-3.10%2B-blue.svg)](https://python.org)
[![R 4.3+](https://img.shields.io/badge/R-4.3%2B-276DC3.svg)](https://r-project.org)

---

## Project Overview

This repository implements four end-to-end omics analysis modules, each independently reproducible and production-ready:

| # | Module | Tools | Language |
|---|--------|-------|----------|
| 1 | **RNA-seq / Transcriptomics** | STAR, featureCounts, DESeq2 | Bash + R |
| 2 | **Proteomics (LFQ)** | MaxQuant → custom DA pipeline | Python |
| 3 | **DNA Methylation** | minfi, ChAMP, limma, bumphunter | R |
| 4 | **AI/ML Multi-Omics Integration** | Autoencoder, Random Forest, SHAP | Python |

---

## Repository Structure

```
bioinformatics_project/
├── 01_transcriptomics/
│   ├── upstream_rnaseq.sh          # QC → Trim → Align → Count
│   └── rnaseq_deseq2_pipeline.R    # DESeq2 DE analysis + plots
├── 02_proteomics/
│   └── proteomics_pipeline.py      # LFQ DA analysis + plots
├── 03_methylation/
│   └── methylation_pipeline.R      # 450K/EPIC DMP/DMR analysis
├── 04_ml_omics/
│   └── multiomics_ml_pipeline.py   # Autoencoder + RF + SHAP
├── data/                            # (gitignored raw data)
├── results/                         # Output tables (CSV)
├── plots/                           # Output figures (PNG)
├── environment.yml                  # Conda environment
├── run_all.sh                       # Master execution script
└── README.md
```

---

## Quick Start

### 1. Clone and Set Up Environment
```bash
git clone https://github.com/YOUR_USERNAME/bioinformatics_project.git
cd bioinformatics_project

conda env create -f environment.yml
conda activate bioinformatics_omics
```

### 2. Run All Pipelines
```bash
bash run_all.sh
```

Or run each module individually (see below).

---

## 🔬 Module Details

### Module 1 — RNA-seq Pipeline
```bash
# Upstream (requires raw FASTQ files in data/fastq/)
bash 01_transcriptomics/upstream_rnaseq.sh

# Differential expression (works with simulated data too)
Rscript 01_transcriptomics/rnaseq_deseq2_pipeline.R
```
**Outputs:** `results/significant_DEGs.csv`, `plots/PCA_plot.png`, `plots/Volcano_plot.png`, `plots/Heatmap_top50_DEGs.png`

---

### Module 2 — Proteomics Pipeline
```bash
# Place MaxQuant proteinGroups.txt in data/ (or runs on simulated data)
python 02_proteomics/proteomics_pipeline.py
```
**Outputs:** `results/proteomics_significant.csv`, `plots/proteomics_volcano.png`, `plots/proteomics_heatmap.png`

---

### Module 3 — Methylation Pipeline
```bash
# Place 450K IDAT files in data/idat/ (or runs on simulated data)
Rscript 03_methylation/methylation_pipeline.R
```
**Outputs:** `results/significant_DMPs.csv`, `plots/methylation_volcano.png`, `plots/methylation_heatmap.png`

---

### Module 4 — AI/ML Multi-Omics Integration
```bash
python 04_ml_omics/multiomics_ml_pipeline.py
```
**Outputs:** `results/cross_validation_AUC.csv`, `results/latent_embeddings.csv`, `results/feature_importance_SHAP.csv`, `plots/shap_summary.png`, `plots/roc_curve.png`, `plots/latent_space.png`

---

## 📊 Methods Summary

### Transcriptomics
- Read trimming: **Trimmomatic** (adapter removal, quality filtering)
- Alignment: **STAR** (2-pass mode, splice-aware)
- Quantification: **featureCounts** (gene-level counts)
- DE analysis: **DESeq2** (negative binomial GLM, Wald test, BH correction)

### Proteomics
- Quantification: **MaxQuant** LFQ (or input from any label-free workflow)
- Normalisation: Median normalisation on log₂ intensities
- Missing value imputation: Column-wise median
- DA analysis: Welch's t-test + BH FDR correction

### Methylation
- Array processing: **minfi** (functional normalisation)
- DMP detection: **limma** (empirical Bayes moderated t-test)
- DMR detection: **bumphunter** / **DMRcate**

### AI/ML Integration
- **Autoencoder**: Unsupervised dimensionality reduction across omics layers (PyTorch)
- **Random Forest**: Supervised classification with feature selection (scikit-learn)
- **SHAP**: Feature explainability and omics contribution analysis
- **Visualisation**: PCA, t-SNE, UMAP for latent space exploration

---

## 📈 Expected Outputs

Each module generates:
- Filtered & normalised data matrices
- Statistically tested results tables (CSV)
- Publication-quality figures (PNG/PDF)
- Summary statistics in console output

---

## 🔧 Data Requirements

| Module | Required Input | Default |
|--------|----------------|---------|
| RNA-seq upstream | `data/fastq/*.fastq.gz` | Simulated |
| RNA-seq DE | `data/raw_counts.csv` | Simulated |
| Proteomics | `data/proteinGroups.txt` | Simulated |
| Methylation | `data/idat/*.idat` | Simulated |
| ML integration | Outputs from modules 1–3 | Simulated |

---

## 📖 References

- Love MI et al. (2014) DESeq2. *Genome Biology*
- Aryee MJ et al. (2014) minfi. *Bioinformatics*
- Cox J, Mann M (2008) MaxQuant. *Nature Biotechnology*
- Lundberg SM, Lee SI (2017) SHAP. *NeurIPS*
- Argelaguet R et al. (2020) MOFA+. *Genome Biology*

---

## 📄 License

MIT License — see [LICENSE](LICENSE) for details.

---

##  Author

Lokaveenasri — Biotechnology Graduate 
📧 dlokaveenasri@gmail.com 
🔗 [LinkedIn](https://www.linkedin.com/in/lokaveenasri-d-5884b1281/) | [GitHub](https://github.com/lokaaaaaaa)
