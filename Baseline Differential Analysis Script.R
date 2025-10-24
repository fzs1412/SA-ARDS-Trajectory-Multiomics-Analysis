```r
# ==============================================================================
# Baseline Differential Analysis Script (Step 1)
# Core Comparison: "Non-resolving" (class 2) vs "Rapid-resolving" (class 1) 
# gene expression differences on Day 1
# log2FoldChange > 0 indicates higher expression in "Non-resolving" patients
# Version updated based on final feedback
# ==============================================================================

# --- Step 1: Load Required R Packages ---
# Ensure all required packages are installed
packages_to_load <- c(
  "dplyr", "data.table", "tibble", "stringr", "DESeq2", 
  "EnhancedVolcano", "pheatmap", "clusterProfiler", "org.Hs.eg.db"
)

for (pkg in packages_to_load) {
  if (!requireNamespace(pkg, quietly = TRUE)) {
    if (pkg %in% c("clusterProfiler", "org.Hs.eg.db", "DESeq2", "EnhancedVolcano")) {
      if (!requireNamespace("BiocManager", quietly = TRUE)) install.packages("BiocManager")
      BiocManager::install(pkg)
    } else {
      install.packages(pkg)
    }
  }
  library(pkg, character.only = TRUE)
}

# Set your main working directory
setwd('d:/OneDrive/R/25SXFX/25ARDS')


# --- Step 2: Data Loading and Preparation ---

# 2.1 Load clinical data and filter correct subtypes
coldata_raw <- fread("./data/data331day1mice.csv") %>%
  dplyr::select(-no, -ID) %>%
  filter(class == 1 | class == 2) %>%
  mutate(class = as.factor(class)) %>%
  column_to_rownames("SampleName")

# 2.2 Load total expression profile data
countdata_path <- "D:/OneDrive/R/01sepsis/OMIX006457/combined_raw_counts.csv"
if (!file.exists(countdata_path)) {
  stop("Error: Expression profile file not found. Please check path: ", countdata_path)
}
combined_counts <- fread(countdata_path, data.table = FALSE)
rownames(combined_counts) <- combined_counts[, 1]
combined_counts <- combined_counts[, -1]

# 2.3 Filter Day 1 samples and match with clinical data
day1_samples <- colnames(combined_counts)[!grepl("d3$|d5$", colnames(combined_counts))]
countdata_d1 <- combined_counts[, day1_samples]

common_samples <- intersect(colnames(countdata_d1), rownames(coldata_raw))
cat("Number of common samples between clinical data and expression profile: ", length(common_samples), "\n")

countdata_final <- countdata_d1[, common_samples]
coldata_final <- coldata_raw[common_samples, ]

stopifnot(all(colnames(countdata_final) == rownames(coldata_final)))


# --- Step 3: DESeq2 Differential Expression Analysis ---

# [Corrected]: Preprocess row names of expression matrix to ensure pure ENSEMBL IDs
countdata_for_deseq <- countdata_final
# Extract pure ENSEMBL ID from 'ENSEMBL|SYMBOL' format
rownames(countdata_for_deseq) <- str_split_i(rownames(countdata_for_deseq), "\\|", 1)
# Remove potential duplicate ENSEMBL IDs, keep the first occurrence
countdata_for_deseq <- countdata_for_deseq[!duplicated(rownames(countdata_for_deseq)), ]

dds <- DESeqDataSetFromMatrix(countData = countdata_for_deseq,  # Use matrix with corrected row names
                              colData   = coldata_final,
                              design    = ~ class)

# Filter low-expression genes
keep <- rowSums(counts(dds)) >= 10
dds <- dds[keep,]
cat("Remaining genes after filtering low-expression genes: ", nrow(dds), "\n")

dds_analyzed <- DESeq(dds)
cat("DESeq2 analysis completed.\n")

# Compare class 2 (Non-resolving) vs class 1 (Rapid-resolving)
res <- results(dds_analyzed, contrast = c("class", "2", "1"))
summary(res)

# [Corrected]: Create unified output directories
dir.create("./data", showWarnings = FALSE)
dir.create("./figure", showWarnings = FALSE)

# [Corrected]: Save results to ./data/ folder
fwrite(as.data.frame(res) %>% rownames_to_column("GeneID"), "./data/DESeq2_baseline_class2_vs_1.csv")
saveRDS(res, file = "./data/baseline_DESeq2_results_object.rds")
saveRDS(dds_analyzed, file = "./data/baseline_dds_analyzed_object.rds")


# --- Step 4: Gene ID Conversion and Preparation ---
count_lc_map_df     <- fread("D:/OneDrive/R/01sepsis/OMIX006457/OMIX006457-02.csv", data.table = FALSE)
original_gene_names <- count_lc_map_df[, 1]

library(stringr)
library(dplyr)
library(tibble) 


initial_gene_map <- data.frame(
  ENSEMBL = str_split_i(original_gene_names, "\\|", 1),
  SYMBOL  = str_split_i(original_gene_names, "\\|", 2),
  stringsAsFactors = FALSE
) %>% filter(!is.na(ENSEMBL) & ENSEMBL != "")

# Merge mappings and resolve many-to-one relationships (keep ENSEMBL with smallest padj for each SYMBOL)
res_df_for_map <- as.data.frame(res) %>%
  rownames_to_column("ENSEMBL") %>%
  left_join(initial_gene_map, by = "ENSEMBL")

gene_map_cleaned <- res_df_for_map %>%
  mutate(SYMBOL = ifelse(is.na(SYMBOL) | SYMBOL == "", ENSEMBL, SYMBOL)) %>%
  group_by(SYMBOL) %>%
  slice_min(order_by = padj, n = 1, with_ties = FALSE) %>%
  ungroup() %>%
  select(ENSEMBL, SYMBOL)

gene_map_cleaned$SYMBOL <- make.unique(gene_map_cleaned$SYMBOL)
cat("Created unique ID-SYMBOL mapping table containing", nrow(gene_map_cleaned), "genes.\n")


# --- Step 5: Visualization - Volcano Plot ---

res_df_for_volcano <- as.data.frame(res) %>%
  rownames_to_column("ENSEMBL") %>%
  left_join(gene_map_cleaned, by = "ENSEMBL")

genes_to_label <- res_df_for_volcano %>%
  filter(padj < 0.05 & abs(log2FoldChange) > 1) %>%
  arrange(padj) %>%
  head(20) %>%
  pull(SYMBOL)

p_volcano <- EnhancedVolcano(
  res_df_for_volcano,
  lab = res_df_for_volcano$SYMBOL,
  selectLab = genes_to_label,
  x = 'log2FoldChange',
  y = 'padj',
  xlab = expression(Log[2]~'Fold Change'),
  ylab = expression(-Log[10]~'Adjusted P-value'),
  axisLabSize = 14,
  pCutoff = 0.05,
  FCcutoff = 1.0,
  pointSize = 2.0,
  labSize = 4.0,
  max.overlaps = Inf,
  title = "Volcano Plot: Non-resolving (class 2) vs. Rapid-resolving (class 1)",
  subtitle = "Day 1 Baseline",
  col = c('grey30', '#2166AC', '#2166AC', '#B2182B'),
  colAlpha = 0.8,
  drawConnectors = TRUE,
  widthConnectors = 0.5
)

# [Corrected]: Save plot to ./figure/ folder
ggsave("./figure/Volcano_Plot_Class2_vs_1.pdf", p_volcano, width=8, height=8)
cat("Volcano plot saved to ./figure/Volcano_Plot_Class2_vs_1.pdf\n")


# --- Step 6: Visualization - Heatmap ---

vsd <- vst(dds_analyzed, blind = FALSE)
vst_matrix <- assay(vsd)

degs_for_heatmap <- res_df_for_volcano %>%
  filter(padj < 0.05 & abs(log2FoldChange) > 1)

cat("Number of significant genes for heatmap: ", nrow(degs_for_heatmap), "\n")

if (nrow(degs_for_heatmap) > 1) {
  heatmap_matrix <- vst_matrix[degs_for_heatmap$ENSEMBL, ]
  # Safely convert row names to unique SYMBOLS
  heatmap_symbols <- gene_map_cleaned$SYMBOL[match(rownames(heatmap_matrix), gene_map_cleaned$ENSEMBL)]
  rownames(heatmap_matrix) <- make.unique(ifelse(is.na(heatmap_symbols), rownames(heatmap_matrix), heatmap_symbols))
  
  annotation_col <- data.frame(
    Group = factor(coldata_final$class, 
                   levels = c("1", "2"),
                   labels = c("Rapid-resolving (C1)", "Non-resolving (C2)")),
    row.names = colnames(heatmap_matrix)
  )
  
  pheatmap_plot <- pheatmap(
    heatmap_matrix,
    scale = "row",
    cluster_rows = TRUE,
    cluster_cols = TRUE,
    show_colnames = FALSE,
    show_rownames = TRUE,
    fontsize_row = 8,
    annotation_col = annotation_col,
    annotation_colors = list(
      Group = c("Rapid-resolving (C1)" = "#2166AC", "Non-resolving (C2)" = "#B2182B")
    ),
    border_color = 'white',
    main = "Heatmap of Differentially Expressed Genes (Day 1)",
    color = colorRampPalette(c("navy", "white", "firebrick3"))(50)
  )
  
  # [Corrected]: Save plot to ./figure/ folder
  ggsave("./figure/Heatmap_Class2_vs_1.pdf", pheatmap_plot, width=8, height=10, device = cairo_pdf)
  cat("Heatmap saved to ./figure/Heatmap_Class2_vs_1.pdf\n")
} else {
  cat("Insufficient number of significant differentially expressed genes, skipping heatmap generation.\n")
}

cat("\n===== Baseline Analysis (Step 1) Completed Successfully =====\n")
```