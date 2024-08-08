This Docker container is based on Rocker 4.4.0 and includes the following R packages: Rsubread, edgeR, csaw, dplyr, biomaRt, and data.table. The container facilitates differential expression/binding analysis using edgeR for RNA-seq/ChIP-seq data, allowing users to choose significance thresholds, normalization methods, and the use of spike-in controls. It includes two key scripts: difexpr.R for RNA-seq and difenrich.R for ChIP-seq data. 

**DE analysis with difexpr.R**
*Usage*: docker run --rm -v $(pwd):/workspace dif_expr_edger:v1 /scripts/difexpr.R <cond1> <cond2> </workspace/path_to_gtf_file> <normalization method (TMM, RLE, upperquartile)> <FDR threshold> <log2fc threshold> </workspace/output_dir> </workspace/path_to_bams_cond1> </workspace/path_to_bams_cond2> <PE (true) or SE (false)> <spike-in (true or false)> <spike-in pattern>

**DB analysis with difenrich.R**
*Usage*: docker run --rm -v $(pwd):/workspace dif_expr_edger:v1 /scripts/difenrich.R <mutant> <wt> <param> <bed_file> <logFC_threshold> <FDR_threshold> <use_spikein> <spikein_pattern> </workspace/output_dir>
