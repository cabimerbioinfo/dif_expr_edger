setwd("/workspace")
library(csaw)
library(edgeR)
library(dplyr)
# Get command line arguments
args <- commandArgs(trailingOnly = TRUE)

# Check if the correct arguments are provided
if (length(args) != 9) {
  stop("Usage: script.R <mutant> <wt> <param> <bed_file> <logFC_threshold> <FDR_threshold> <use_spikein> <spikein_pattern> <output_dir>")
}

# Assign arguments to variables
mutant <- args[1]
wt <- args[2]
param <- args[3]
bed_file <- args[4]
logFC_threshold <- as.numeric(args[5])
FDR_threshold <- as.numeric(args[6])
use_spikein <- as.logical(args[7])
spikein_pattern <- args[8]
output_dir<-args[9]



# Construct the comparison name
comp <- paste0(mutant, "_vs_", wt)

# List BAM files
bamFiles_wt <- list.files(pattern = paste0("^", wt, ".*\\.bam$"))
bamFiles_mutant <- list.files(pattern = paste0("^", mutant, ".*\\.bam$"))

# Check if BAM files are found
if (length(bamFiles_wt) != 2 || length(bamFiles_mutant) != 2) {
  stop("No BAM files found for ", wt, " or ", mutant)
}

# Combine BAM files into a single vector
bamFiles <- c(bamFiles_wt, bamFiles_mutant)
bam_names <- tools::file_path_sans_ext(basename(bamFiles))
print("Counting reads in large windows...")

# Count reads in large windows
large_bins_counts <- windowCounts(bamFiles, width = 10000, bin = TRUE, param = readParam(pe = param))

if (use_spikein) {
  print("Using spike-in normalization...")

  # Identify spike-in and endogenous data based on chromosome names
  is_spikein <- grepl(spikein_pattern, seqnames(rowRanges(large_bins_counts)))
  spike_data <- large_bins_counts[is_spikein,]
  endog_data <- large_bins_counts[!is_spikein,]

# Normalize endogenous data using spike-in data
  norm_factors <- normFactors(spike_data, se.out = endog_data)
  library_sizes <- endog_data$totals
  final_factors <- norm_factors$norm.factors * library_sizes
  perMillion_factors <- (final_factors / 1000000)^-1
  names(perMillion_factors) <- bam_names

  # Write normalization factors to disk
  write.table(t(perMillion_factors), file = paste0(output_dir,comp, "_normfactors_bamcov.txt"), sep = "\t", quote = FALSE,
              row.names = FALSE)

} else {
  print("Calculating normalization factors without spike-in...")

  # Compute normalization factors
  library_sizes <- large_bins_counts$totals
  norm_factors <- normFactors(large_bins_counts)
  final_factors <- norm_factors$norm.factors * library_sizes
  perMillion_factors <- (final_factors / 1000000)^-1
  names(perMillion_factors) <- bam_names

  # Write normalization factors to disk
  write.table(t(perMillion_factors), file = paste0(output_dir,comp, "_normfactors_bamcov.txt"), sep = "\t", quote = FALSE,
              row.names = FALSE)
}

print("Performing differential binding analysis...")

# Read the BED file provided by the user
merged <- read.table(file = bed_file, header = FALSE)
colnames(merged) <- c("seqnames", "start", "end")
merged <- GRanges(merged)

peak_counts <- regionCounts(bamFiles, merged, param = readParam(pe = param))
TMM_counts <- assay(peak_counts) * perMillion_factors
colnames(TMM_counts) <- bam_names
merged <- as.data.frame(merged)
TMM_counts <- cbind(merged[, 1:3], TMM_counts)
write.table(TMM_counts, file = paste0(output_dir,comp, "_TMM_counts.tsv"), row.names = FALSE, quote = FALSE, sep = "\t")

print("Preparing design for comparisons...")

# Design for comparisons
grouping <- factor(c(rep(wt, 2), rep(mutant, 2)))
design <- model.matrix(~0 + grouping)
colnames(design) <- levels(grouping)
design

print("Performing differential analysis...")

# Differential binding analysis
db <- asDGEList(peak_counts, norm.factors = norm_factors$norm.factors)
db <- estimateDisp(db, design)
fit_db <- glmQLFit(db, design, robust = TRUE)

contrast <- paste0(mutant, "-", wt)
comp_db <- makeContrasts(contrast, levels = design)
comp_db <- glmLRT(fit_db, contrast = comp_db)
comp_table <- comp_db$table

# Calculate adjusted p-values (FDR)
comp_table$FDR <- p.adjust(comp_table$PValue, method = "BH")

# Combine the merged data with the results
comp_table <- data.frame(cbind(merged[, 1:3], comp_table))
write.table(comp_table, file = paste0(comp, "_db.tsv"), row.names = FALSE, col.names = TRUE, sep = "\t")
comp_table <- read.table(file = paste0(comp, "_db.tsv"), header = TRUE)

print("Filtering results based on thresholds...")

# Filter results based on given thresholds
mutant_increase <- comp_table %>% filter(logFC >= logFC_threshold, FDR <= FDR_threshold)
mutant_decrease <- comp_table %>% filter(logFC <= -logFC_threshold, FDR <= FDR_threshold)
mutant_nochange <- comp_table %>% filter(abs(logFC) < logFC_threshold | FDR > FDR_threshold)

print("Writing results to output files...")

# Write results to output files
write.table(mutant_increase, file = paste0(output_dir,mutant, "_enriched.tsv"), row.names = FALSE, quote = FALSE, sep = "\t")
write.table(mutant_decrease, file = paste0(output_dir,wt, "_enriched.tsv"), row.names = FALSE, quote = FALSE, sep = "\t")
write.table(mutant_nochange, file = paste0(output_dir,comp, "_nochange.tsv"), row.names = FALSE, quote = FALSE, sep = "\t")

print("Process completed.")
