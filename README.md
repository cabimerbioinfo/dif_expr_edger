This Docker container is based on Rocker 4.4.0 and includes the following R packages: Rsubread, edgeR, csaw, dplyr, biomaRt, and data.table. The container facilitates differential expression analysis using edgeR, allowing users to choose significance thresholds, normalization methods, and the use of spike-in controls. We are also working on adding scripts for differential enrichment analysis of ChIP-Seq data.

*Usage*: edger_de.R <normalization method (TMM, RLE, upperquartile)> <output_dir> <PE (true) or SE (false)> <spike-in (true or false)>
