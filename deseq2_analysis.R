# =========================================================
# RNA-seq differential expression analysis (DESeq2)
# - Reads BAM files and a GTF annotation
# - Counts reads per gene using featureCounts (Rsubread)
# - Runs DESeq2 (Zika vs Mock)
# - Generates PCA plot and heatmap of top variable genes
# =========================================================

# Libraries used
library(Rsubread)    # Read counting (featureCounts)
library(DESeq2)      # Differential expression analysis
library(edgeR)       # Alternative DE approach (not used below, kept for reference)
library(tximport)    # Import transcript quantification (not used below, kept for reference)
library(pheatmap)    # Heatmap visualization
library(ggplot2)     # PCA visualization
library(rtracklayer) # Working with GTF files

# -----------------------------
# Paths (CHANGE to your setup)
# -----------------------------
# NOTE: Replace Windows paths with your local paths or use relative paths in a repo.
gtf_file <- "Homo_sapiens.GRCh37.87.chr.gtf/Homo_sapiens.GRCh37.87.chr.gtf"
bam_dir  <- "bam bai(1)/bam bai"

# -----------------------------
# Load GTF annotation
# -----------------------------
gtf_data <- import(gtf_file, format = "gtf")
head(gtf_data)

# Keep only gene features
gtf_genes <- gtf_data[gtf_data$type == "gene"]
head(gtf_genes)

# -----------------------------
# Load BAM files
# -----------------------------
bam_files <- list.files(bam_dir, pattern = "\\.bam$", full.names = TRUE)

# Count reads for genes (all BAMs together)
counts <- featureCounts(
  files = bam_files,
  annot.ext = gtf_file,
  isGTFAnnotationFile = TRUE,
  GTF.featureType = "gene",
  GTF.attrType = "gene_name",
  useMetaFeatures = TRUE,
  isPairedEnd = TRUE,
  nthreads = 4
)

# -----------------------------
# Separate paired-end (PE) and single-end (SE) BAMs
# -----------------------------
bam_pe <- c(
  "bam bai(1)/bam bai/SRR3191542.bam",
  "bam bai(1)/bam bai/SRR3191543.bam",
  "bam bai(1)/bam bai/SRR3191544.bam",
  "bam bai(1)/bam bai/SRR3191545.bam"
)

bam_se <- c(
  "bam bai(1)/bam bai/SRR3194428.bam",
  "bam bai(1)/bam bai/SRR3194429.bam",
  "bam bai(1)/bam bai/SRR3194430.bam",
  "bam bai(1)/bam bai/SRR3194431.bam"
)

# Count reads for paired-end BAMs
counts_pe <- featureCounts(
  files = bam_pe,
  annot.ext = gtf_file,
  isGTFAnnotationFile = TRUE,
  GTF.featureType = "gene",
  GTF.attrType = "gene_name",
  useMetaFeatures = TRUE,
  isPairedEnd = TRUE,
  nthreads = 4
)

# Count reads for single-end BAMs
counts_se <- featureCounts(
  files = bam_se,
  annot.ext = gtf_file,
  isGTFAnnotationFile = TRUE,
  GTF.featureType = "gene",
  GTF.attrType = "gene_name",
  useMetaFeatures = TRUE,
  isPairedEnd = FALSE,
  nthreads = 4
)

# Merge counts (PE + SE) into one matrix
final_counts <- cbind(counts_pe$counts, counts_se$counts)

# -----------------------------
# Build DESeq2 dataset
# -----------------------------
dds <- DESeqDataSetFromMatrix(
  countData = final_counts,
  colData = data.frame(condition = c(rep("Zika", 4), rep("Mock", 4))),
  design = ~ condition
)

# Run DESeq2
dds <- DESeq(dds)
res <- results(dds)
head(res)

# -----------------------------
# PCA plot
# -----------------------------
vsd <- vst(dds)
pcaData <- plotPCA(vsd, intgroup = "condition", returnData = TRUE)

ggplot(pcaData, aes(PC1, PC2, color = condition)) +
  geom_point(size = 4) +
  theme_minimal() +
  ggtitle("RNA-seq PCA")

# -----------------------------
# Heatmap of top variable genes
# -----------------------------
# Convert colData to data.frame for pheatmap annotation
annotation_col <- as.data.frame(colData(dds))
rownames(annotation_col) <- colnames(assay(vsd))
annotation_col$condition <- as.factor(annotation_col$condition)

# Calculate gene-wise variance and select top 25 most variable genes
gene_vars <- apply(assay(vsd), 1, var)
topVarGenes <- order(gene_vars, decreasing = TRUE)
topVarGenesSubset <- head(topVarGenes, 25)

# Heatmap
pheatmap(
  assay(vsd)[topVarGenesSubset, ],
  scale = "row",
  annotation_col = annotation_col[, "condition", drop = FALSE],
  fontsize_row = 8,
  fontsize_col = 10,
  angle_col = 45,
  clustering_distance_rows = "euclidean",
  clustering_distance_cols = "euclidean",
  clustering_method = "complete"
)
