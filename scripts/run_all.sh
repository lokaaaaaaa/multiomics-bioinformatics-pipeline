#!/usr/bin/env bash
# ============================================================
# run_all.sh — Execute all four omics pipelines in order
# ============================================================
set -euo pipefail
echo "================================================================"
echo "  Multi-Omics Pipeline — Master Execution Script"
echo "================================================================"

echo ""
echo "[1/4] RNA-seq Differential Expression (DESeq2)..."
Rscript 01_transcriptomics/rnaseq_deseq2_pipeline.R

echo ""
echo "[2/4] Proteomics Differential Abundance..."
python 02_proteomics/proteomics_pipeline.py

echo ""
echo "[3/4] DNA Methylation Analysis (minfi + limma)..."
Rscript 03_methylation/methylation_pipeline.R

echo ""
echo "[4/4] AI/ML Multi-Omics Integration..."
python 04_ml_omics/multiomics_ml_pipeline.py

echo ""
echo "================================================================"
echo "  ALL PIPELINES COMPLETE"
echo "  Results  → results/"
echo "  Figures  → plots/"
echo "================================================================"
