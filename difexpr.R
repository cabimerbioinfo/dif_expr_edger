# Load required packages
setwd("/workspace")
library(Rsubread)
library(edgeR)
library(biomaRt)
args <- commandArgs(trailingOnly = TRUE)

# Check if the required arguments are provided
if (length(args) < 10) {
  stop("Usage: edger_de.R <cond1> <cond2> <gtf> <normalization method (TMM, RLE, upperquartile)> <fdr th> <log2fc th> <output_dir> <bams> <SE or PE (true or false)> <spike-in (true or false)> <spike-in p>")
}

condition_names <- c(args[1], args[2])

# Input GTF:
gtfFile <- args[3]

# Normalization method:
normalization_method <- args[4]

# Thresholds:
fdr_threshold <- as.numeric(args[5])
log2fc_threshold <- as.numeric(args[6])
output_dir <- args[7]

# Extract bam files
bamFiles <- args[8:(length(args)-3)]
sampleNames <- tools::file_path_sans_ext(basename(bamFiles))
print("Sample names:")
print(sampleNames)

# Assign reads to genes with featureCounts
countData <- featureCounts(
  files = bamFiles,
  annot.ext = gtfFile,
  isGTFAnnotationFile = TRUE,
  GTF.featureType = "exon",
  GTF.attrType = "gene_id",
  strandSpecific = 2,
  isPairedEnd = as.logical(args[length(args)-2])  # Convert third-to-last argument to logical
)
print("Reads counted")

# Create DataFrame with condition information
conditions <- rep(condition_names, each = length(sampleNames) / length(condition_names))
colData <- data.frame(
  condition = factor(conditions),
  row.names = sampleNames
)
print("Condition DataFrame:")
print(colData)
print("DataFrame built")

# Create DGEList object
dge <- DGEList(counts = countData$counts, group = colData$condition)
print("DGEList object created")

# Normalization
if (normalization_method == "TMM") {
  dge <- calcNormFactors(dge, method = "TMM")
} else if (normalization_method == "RLE") {
  dge <- calcNormFactors(dge, method = "RLE")
} else if (normalization_method == "upperquartile") {
  dge <- calcNormFactors(dge, method = "upperquartile")
} else {
  stop("Invalid normalization method. Use one of: TMM, RLE, upperquartile")
}
print(paste("Normalization method used:", normalization_method))

# If spike-in is used
if (as.logical(args[length(args)-1])) {
  # Spike-in specific processing
  gene_names <- rownames(countData$counts)
  rows_with_dm <- grep(args[length(args)], gene_names)

  # Calculate normalization factors using spike-ins
  dge <- calcNormFactors(dge, method = normalization_method, lib.size = rows_with_dm)
}
print("Normalization factors calculated")

# Save normalization factors to a file
norm_factors_file <- paste0(output_dir, condition_names[1], "_", condition_names[2], "_norm_factors.txt")
write.table(dge$samples$norm.factors, file = norm_factors_file, sep = "\t", quote = FALSE, col.names = NA)
print(paste("Normalization factors saved to:", norm_factors_file))

# Estimate dispersion
dge <- estimateDisp(dge)
print("Dispersion estimated")

# Perform differential expression analysis
et <- exactTest(dge, pair = condition_names)
results <- et$table
print("Differential expression analysis completed")

# Calculate FDR
results$FDR <- p.adjust(results$PValue, method = "BH")
results$log10FDR <- -log10(results$FDR)
print("FDR calculated")

# Write normalized counts to file
norm_counts <- cpm(dge, normalized.lib.sizes = TRUE)
norm_counts <- as.data.frame(norm_counts)
norm_counts$ensembl_id <- rownames(norm_counts)
print("Normalized counts calculated")

# Retrieve gene symbols
ensembl <- useMart("ensembl", dataset = "hsapiens_gene_ensembl")
gene_symbols <- getBM(attributes = c("ensembl_gene_id", "external_gene_name"),
                      filters = "ensembl_gene_id",
                      values = norm_counts$ensembl_id,
                      mart = ensembl)
colnames(gene_symbols) <- c("ensembl_id", "gene_symbol")

# Merge gene symbols with normalized counts
norm_counts <- merge(norm_counts, gene_symbols, by = "ensembl_id", all.x = TRUE)

# Reorder columns to have gene symbol first
norm_counts <- norm_counts[, c("ensembl_id", "gene_symbol", setdiff(colnames(norm_counts), c("ensembl_id", "gene_symbol")))]

# Write normalized counts with gene symbols to file
normc <- paste0(output_dir, condition_names[1], "_", condition_names[2], "_normalizedcounts.tsv")
write.table(norm_counts, file = normc, sep = "\t", quote = FALSE, row.names = FALSE)
print(paste("Normalized counts with gene symbols saved to:", normc))

# Filter genes
results_df <- as.data.frame(results)
clean_results_df <- results_df[!is.na(results_df$FDR) & !is.na(results_df$logFC), ]
print("Filtered results for NA values")

ensembl_id <- rownames(clean_results_df)
clean_results_df$ensembl_id <- ensembl_id
df_symbols <- merge(clean_results_df, gene_symbols, by = "ensembl_id", all.x = TRUE)
print("Gene symbols retrieved and merged with results")

# Write edgeR results to file
edger_results <- paste0(output_dir, condition_names[1], "_", condition_names[2], "_difexpression.tsv")
write.table(df_symbols, file = edger_results, sep = "\t", quote = FALSE, row.names = FALSE)
print(paste("edgeR results saved to:", edger_results))

# Filter by thresholds
upregulated_genes <- df_symbols[df_symbols$FDR < fdr_threshold & df_symbols$logFC > log2fc_threshold, ]
downregulated_genes <- df_symbols[df_symbols$FDR < fdr_threshold & df_symbols$logFC < -log2fc_threshold, ]

# Save upregulated genes to a TSV file
upr <- paste0(output_dir, condition_names[2], "_upregulated.tsv")
write.table(upregulated_genes, file = upr, sep = "\t", quote = FALSE, row.names = FALSE)
print(paste("Upregulated genes saved to:", upr))

# Save downregulated genes to a TSV file
downr <- paste0(output_dir, condition_names[2], "_downregulated.tsv")
write.table(downregulated_genes, file = downr, sep = "\t", quote = FALSE, row.names = FALSE)
print(paste("Downregulated genes saved to:", downr))
