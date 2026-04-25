#!/usr/bin/env bash
# ============================================================
# Upstream RNA-seq Pipeline: QC → Trim → Align → Count
# Tools: FastQC, Trimmomatic, STAR, featureCounts (Subread)
# ============================================================
set -euo pipefail

# ── Configuration ─────────────────────────────────────────────
THREADS=8
GENOME_DIR="references/star_index"
GTF="references/genome.gtf"
FASTQ_DIR="data/fastq"
TRIM_DIR="data/trimmed"
ALIGN_DIR="data/aligned"
COUNT_DIR="data/counts"
QC_DIR="results/qc"
ADAPTERS="/opt/Trimmomatic/adapters/TruSeq3-PE.fa"   # adjust path

mkdir -p "$TRIM_DIR" "$ALIGN_DIR" "$COUNT_DIR" "$QC_DIR"

# ── Step 1: FastQC on raw reads ────────────────────────────────
echo "==> [1/5] Running FastQC on raw reads..."
fastqc -t "$THREADS" -o "$QC_DIR" "$FASTQ_DIR"/*.fastq.gz
multiqc "$QC_DIR" -o "$QC_DIR/multiqc_raw"

# ── Step 2: Trim with Trimmomatic (PE mode) ────────────────────
echo "==> [2/5] Trimming adapters with Trimmomatic..."
for R1 in "$FASTQ_DIR"/*_R1_001.fastq.gz; do
  SAMPLE=$(basename "$R1" _R1_001.fastq.gz)
  R2="$FASTQ_DIR/${SAMPLE}_R2_001.fastq.gz"

  trimmomatic PE -threads "$THREADS" \
    "$R1" "$R2" \
    "$TRIM_DIR/${SAMPLE}_R1_paired.fastq.gz"   "$TRIM_DIR/${SAMPLE}_R1_unpaired.fastq.gz" \
    "$TRIM_DIR/${SAMPLE}_R2_paired.fastq.gz"   "$TRIM_DIR/${SAMPLE}_R2_unpaired.fastq.gz" \
    ILLUMINACLIP:"$ADAPTERS":2:30:10:2:True \
    LEADING:3 TRAILING:3 SLIDINGWINDOW:4:15 MINLEN:36 \
    2>> "results/trimmomatic_${SAMPLE}.log"

  echo "  Trimmed: $SAMPLE"
done

# ── Step 3: Build STAR genome index (run once) ─────────────────
if [ ! -d "$GENOME_DIR" ]; then
  echo "==> [3a] Building STAR genome index..."
  mkdir -p "$GENOME_DIR"
  STAR --runMode genomeGenerate \
       --genomeDir "$GENOME_DIR" \
       --genomeFastaFiles references/genome.fa \
       --sjdbGTFfile "$GTF" \
       --sjdbOverhang 149 \
       --runThreadN "$THREADS"
fi

# ── Step 4: Align with STAR ────────────────────────────────────
echo "==> [4/5] Aligning reads with STAR..."
for R1 in "$TRIM_DIR"/*_R1_paired.fastq.gz; do
  SAMPLE=$(basename "$R1" _R1_paired.fastq.gz)
  R2="$TRIM_DIR/${SAMPLE}_R2_paired.fastq.gz"

  STAR --runThreadN "$THREADS" \
       --genomeDir "$GENOME_DIR" \
       --readFilesIn "$R1" "$R2" \
       --readFilesCommand zcat \
       --outSAMtype BAM SortedByCoordinate \
       --outSAMattributes NH HI AS NM MD \
       --outFileNamePrefix "$ALIGN_DIR/${SAMPLE}_" \
       --quantMode GeneCounts \
       --outFilterMultimapNmax 20 \
       --alignSJoverhangMin 8 \
       --alignSJDBoverhangMin 1 \
       --outFilterMismatchNmax 999 \
       --outFilterMismatchNoverReadLmax 0.04 \
       --alignIntronMin 20 \
       --alignIntronMax 1000000 \
       --alignMatesGapMax 1000000

  samtools index "$ALIGN_DIR/${SAMPLE}_Aligned.sortedByCoord.out.bam"
  echo "  Aligned: $SAMPLE"
done

# ── Step 5: Count reads with featureCounts ─────────────────────
echo "==> [5/5] Counting reads with featureCounts..."
BAM_FILES=$(ls "$ALIGN_DIR"/*_Aligned.sortedByCoord.out.bam | tr '\n' ' ')

featureCounts \
  -T "$THREADS" \
  -p \
  -a "$GTF" \
  -o "$COUNT_DIR/raw_counts.txt" \
  $BAM_FILES

# Tidy count table (remove featureCounts metadata columns)
cut -f1,7- "$COUNT_DIR/raw_counts.txt" | \
  sed '1d' > "$COUNT_DIR/raw_counts_clean.txt"

# Run MultiQC on all results
multiqc . -o "results/multiqc_final"

echo ""
echo "✅ Upstream pipeline complete!"
echo "   Count matrix → $COUNT_DIR/raw_counts_clean.txt"
echo "   Proceed with: Rscript 01_transcriptomics/rnaseq_deseq2_pipeline.R"
