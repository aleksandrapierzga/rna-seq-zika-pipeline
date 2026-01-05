#!/bin/bash

# =========================================================
# RNA-seq pipeline (SRA -> QC -> trimming -> STAR mapping -> BAM -> featureCounts)
# Notes:
# - Downloads 8 SRA samples (limited to 1,000,000 reads per sample)
# - Runs FastQC, trims reads with Trimmomatic, maps with STAR
# - Converts SAM to sorted/indexed BAM and generates gene-level counts (featureCounts)
# - Prepares outputs for DESeq2 analysis in R
# =========================================================

# Project directories
PROJECT_DIR="$HOME/rna_seq_project"
READS_DIR="$PROJECT_DIR/reads"
REFS_DIR="$PROJECT_DIR/refs"
OUTPUT_DIR="$PROJECT_DIR/mapping_results"
STAR_INDEX_DIR="$REFS_DIR/STAR_index"
THREADS=4

# Create required directories
mkdir -p "$READS_DIR" "$REFS_DIR" "$OUTPUT_DIR" "$PROJECT_DIR/fastqc_reports" "$PROJECT_DIR/trimmed_reads" "$PROJECT_DIR/deseq2_input"

echo "=== Downloading data from SRA ==="
# 8 SRR samples; download up to 1 million reads per sample
for SRR in SRR3191542 SRR3191543 SRR3191544 SRR3191545 SRR3194428 SRR3194429 SRR3194430 SRR3194431; do
    fasterq-dump --split-files --max-reads 1000000 "$SRR" -O "$READS_DIR"
done

echo "=== Quality control with FastQC ==="
# Generate QC reports for all FASTQ files
fastqc "$READS_DIR"/*.fastq -o "$PROJECT_DIR/fastqc_reports/"

echo "=== Trimming low-quality reads (Trimmomatic) ==="
# Paired-end trimming for files matching *_1.fastq and *_2.fastq
for FILE in "$READS_DIR"/*_1.fastq; do
    BASE=$(basename "$FILE" _1.fastq)
    trimmomatic PE -threads "$THREADS" \
        "$READS_DIR/${BASE}_1.fastq" "$READS_DIR/${BASE}_2.fastq" \
        "$PROJECT_DIR/trimmed_reads/${BASE}_1_paired.fastq" "$PROJECT_DIR/trimmed_reads/${BASE}_1_unpaired.fastq" \
        "$PROJECT_DIR/trimmed_reads/${BASE}_2_paired.fastq" "$PROJECT_DIR/trimmed_reads/${BASE}_2_unpaired.fastq" \
        ILLUMINACLIP:TruSeq3-PE.fa:2:30:10 LEADING:3 TRAILING:3 SLIDINGWINDOW:4:15 MINLEN:36
done

echo "=== Building STAR genome index (hg19) ==="
# Requires reference FASTA: $REFS_DIR/hg19.fa
mkdir -p "$STAR_INDEX_DIR"
STAR --runThreadN "$THREADS" \
     --runMode genomeGenerate \
     --genomeDir "$STAR_INDEX_DIR" \
     --genomeFastaFiles "$REFS_DIR/hg19.fa" \
     --genomeSAindexNbases 9

echo "=== STAR mapping ==="
# Map paired reads using STAR
for FILE in "$PROJECT_DIR/trimmed_reads"/*_1_paired.fastq; do
    BASE=$(basename "$FILE" _1_paired.fastq)
    STAR --runThreadN "$THREADS" \
         --genomeDir "$STAR_INDEX_DIR" \
         --readFilesIn "$PROJECT_DIR/trimmed_reads/${BASE}_1_paired.fastq" "$PROJECT_DIR/trimmed_reads/${BASE}_2_paired.fastq" \
         --outFileNamePrefix "$OUTPUT_DIR/${BASE}_"
done

echo "=== Converting SAM -> BAM, sorting and indexing ==="
# Convert STAR outputs to sorted/indexed BAM files
for SAM_FILE in "$OUTPUT_DIR"/*.sam; do
    BASE=$(basename "$SAM_FILE" .sam)
    samtools view -S -b "$SAM_FILE" > "$OUTPUT_DIR/${BASE}.bam"
    samtools sort "$OUTPUT_DIR/${BASE}.bam" -o "$OUTPUT_DIR/${BASE}_sorted.bam"
    samtools index "$OUTPUT_DIR/${BASE}_sorted.bam"
done

echo "=== Counting reads with featureCounts ==="
# Requires annotation GTF: $REFS_DIR/hg19.gtf
featureCounts -T "$THREADS" \
              -a "$REFS_DIR/hg19.gtf" \
              -o "$OUTPUT_DIR/gene_counts.txt" \
              "$OUTPUT_DIR"/*_sorted.bam

echo "=== Preparing files for DESeq2 in R ==="
# Copy counts and BAM files to a dedicated directory for downstream analysis
cp "$OUTPUT_DIR/gene_counts.txt" "$PROJECT_DIR/deseq2_input/"
cp "$OUTPUT_DIR"/*.bam "$PROJECT_DIR/deseq2_input/"

echo "=== Results are in: $OUTPUT_DIR ==="
