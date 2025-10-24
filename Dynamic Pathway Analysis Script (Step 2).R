```r
# ==============================================================================
# Dynamic Pathway Analysis Script (Step 2)
# Comparing Intra-Phenotype Dynamic Gene Changes (Day5 vs Day1) between 
# "Rapid-resolving" (class 1) and "Non-resolving" (class 2)
# Version updated based on feedback
# ==============================================================================
# Create output directories if they don't exist
if (!dir.exists("./table_dynamic")) {
  dir.create("./table_dynamic", recursive = TRUE, showWarnings = FALSE)
}
if (!dir.exists("./figure")) {
  dir.create("./figure", recursive = TRUE, showWarnings = FALSE)
}

setwd('d:/OneDrive/R/25SXFX/25ARDS')

# --- Step 0: Load Required R Packages ---
# Ensure all required packages are installed
packages_to_load <- c(
  "dplyr", "data.table", "tibble", "DESeq2", "EnhancedVolcano",
  "pheatmap", "ggvenn", "clusterProfiler", "org.Hs.eg.db", "ggplot2",
  "Cairo", "purrr"
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


# --- Step 1: Define a Reusable Paired Analysis Function ---
# This function completes the entire workflow from data preparation to DESeq2 analysis for a specified phenotype (class)

perform_paired_deseq2 <- function(target_class, coldata_path, countdata_path) {
  
  cat(paste0("\n===== Starting Analysis for Class: ", target_class, " =====\n"))
  
  # 1.1 Load and filter clinical data
  coldata_all <- fread(coldata_path)
  coldata <- coldata_all %>%
    filter(class == target_class) %>%
    as.data.frame()
  
  if (nrow(coldata) == 0) {
    stop(paste("Error: No samples for class", target_class, "found in clinical data."))
  }
  
  # 1.2 Load expression profile data
  # Use fread for faster reading and retain original row names
  countdata_all <- fread(countdata_path, data.table = FALSE)
  rownames(countdata_all) <- countdata_all[, 1]
  countdata_all <- countdata_all[, -1]
  
  # 1.3 Identify patients with both Day1 and Day5 samples in the specified class
  # Extract base IDs (remove _d1/_d3/_d5 suffixes)
  coldata$base_id <- sub("_(d1|d3|d5)$", "", coldata$SampleName)
  count_colnames_df <- data.frame(full_name = colnames(countdata_all))
  count_colnames_df$base_id <- sub("_(d1|d3|d5)$", "", count_colnames_df$full_name)
  
  # Filter base IDs with both Day1 and Day5 data
  paired_patients <- count_colnames_df %>%
    filter(grepl("_d1$|_d5$", full_name)) %>%
    group_by(base_id) %>%
    filter(any(grepl("_d1$", full_name)) & any(grepl("_d5$", full_name))) %>%
    pull(base_id) %>%
    unique()
  
  # Intersect with patients in the current class
  final_paired_patients <- intersect(paired_patients, coldata$base_id)
  
  if (length(final_paired_patients) == 0) {
    cat(paste("Warning: No patients with both Day1 and Day5 samples found in class", target_class, ". Skipping this group.\n"))
    return(NULL)
  }
  cat(paste("Found", length(final_paired_patients), "paired samples in Class", target_class, ".\n"))
  
  # 1.4 Prepare final countdata and coldata
  final_sample_names <- count_colnames_df %>%
    filter(base_id %in% final_paired_patients & grepl("_d1$|_d5$", full_name)) %>%
    pull(full_name)
  
  countdata_final <- countdata_all[, final_sample_names]
  
  coldata_final <- data.frame(
    SampleName = colnames(countdata_final),
    row.names = colnames(countdata_final)
  )
  coldata_final$patient <- sub("_(d1|d5)$", "", coldata_final$SampleName)
  coldata_final$time <- factor(ifelse(grepl("_d1$", coldata_final$SampleName), "d1", "d5"), levels = c("d1", "d5"))
  
  stopifnot(all(table(coldata_final$patient) == 2))  # Ensure each patient has 2 time points
  
  # 1.5 Create gene ID mapping table
  # Correction: Create mapping from original row names of combined_counts (not filtered countdata_final)
  original_gene_map_df <- data.frame(
    ENSEMBL_VERSION = rownames(countdata_all),
    stringsAsFactors = FALSE
  )
  
  gene_map <- original_gene_map_df %>%
    mutate(ENSEMBL = gsub("\\..*$", "", ENSEMBL_VERSION)) %>%
    left_join(
      # [Corrected]: Explicitly use AnnotationDbi::select to avoid conflict with dplyr::select
      AnnotationDbi::select(org.Hs.eg.db, keys = .$ENSEMBL, columns = "SYMBOL", keytype = "ENSEMBL"),
      by = "ENSEMBL"
    ) %>%
    filter(!is.na(SYMBOL) & SYMBOL != "") %>%
    distinct(ENSEMBL, .keep_all = TRUE)
  
  # Prepare matrix for DESeq2 (row names = ENSEMBL IDs without version)
  countdata_deseq <- countdata_final
  rownames(countdata_deseq) <- gsub("\\..*$", "", rownames(countdata_deseq))
  countdata_deseq <- countdata_deseq[rownames(countdata_deseq) %in% gene_map$ENSEMBL, ]
  
  # --- 1.6 Run DESeq2 Paired Analysis ---
  dds <- DESeqDataSetFromMatrix(countData = countdata_deseq,
                                colData = coldata_final,
                                design = ~ patient + time)
  # Filter low-expression genes
  keep <- rowSums(counts(dds)) >= 10
  dds <- dds[keep,]
  
  dds <- DESeq(dds)
  res <- results(dds, contrast = c("time", "d5", "d1"))  # Explicit comparison: Day5 vs Day1
  
  # Organize results and merge with gene names
  res_df <- as.data.frame(res) %>%
    rownames_to_column("ENSEMBL") %>%
    left_join(gene_map, by = "ENSEMBL")
  
  cat(paste0("===== Analysis for Class ", target_class, " Completed =====\n"))
  
  return(res_df)
}

# --- Step 2: Analyze the Two Phenotypes Separately ---

# Define file paths
coldata_file <- "./data/data331day1mice.csv"
countdata_file <- "D:/OneDrive/R/01sepsis/OMIX006457/combined_raw_counts.csv"

# [Corrected] Analyze "Rapid-resolving" (class 1)
results_resolving <- perform_paired_deseq2(
  target_class = 1,
  coldata_path = coldata_file,
  countdata_path = countdata_file
)

# [Corrected] Analyze "Non-resolving" (class 2)
results_non_resolving <- perform_paired_deseq2(
  target_class = 2,
  coldata_path = coldata_file,
  countdata_path = countdata_file
)

# Save results
save(results_resolving, results_non_resolving, 
     file = "./data/Dynamic_Pathway_Analysis_Step2_DESeq2_Results.RData")
results_file <- "./data/Dynamic_Pathway_Analysis_Step2_DESeq2_Results.RData"

# Load saved results
load(results_file)

# ==============================================================================
# Dynamic Pathway Analysis Script (Step 2 - GSEA Follow-up)
# ==============================================================================

# --- Step 1: Load Additional Packages for GSEA ---
packages_to_load <- c(
  "dplyr", "data.table", "tibble", "stringr", "clusterProfiler", 
  "org.Hs.eg.db", "ggplot2", "Cairo", "purrr", "fgsea",
  "msigdbr", "gridExtra", "officer", "flextable"
)

for (pkg in packages_to_load) {
  if (!requireNamespace(pkg, quietly = TRUE)) {
    if (pkg %in% c("clusterProfiler", "org.Hs.eg.db", "fgsea", "msigdbr")) {
      if (!requireNamespace("BiocManager", quietly = TRUE)) install.packages("BiocManager")
      BiocManager::install(pkg)
    } else {
      install.packages(pkg)
    }
  }
  library(pkg, character.only = TRUE)
}


# --- Step 2: Define GSEA Analysis Function ---
# This function takes a DESeq2 results data frame and performs complete GSEA analysis

run_gsea_analysis <- function(deseq_results, class_name, pathways_db) {
  
  cat(paste0("\n\n===== Starting GSEA Analysis for '", class_name, "' =====\n"))
  
  # 2.1 Prepare pre-ranked gene list for GSEA
  if (!all(c("stat", "SYMBOL") %in% names(deseq_results))) {
    stop("Error: Input DESeq2 results missing 'stat' or 'SYMBOL' column.")
  }
  
  ranks_df <- deseq_results %>%
    filter(!is.na(stat) & !is.na(SYMBOL) & SYMBOL != "") %>%
    group_by(SYMBOL) %>%
    slice_max(order_by = abs(stat), n = 1, with_ties = FALSE) %>%
    ungroup()
  
  ranks <- ranks_df$stat
  names(ranks) <- ranks_df$SYMBOL
  ranks <- sort(ranks, decreasing = TRUE)
  
  cat(paste("Created ranked list with", length(ranks), "genes for", class_name, ".\n"))
  
  # 2.2 Run GSEA
  fgsea_res <- fgsea(pathways = pathways_db, stats = ranks, minSize = 15, maxSize = 500, nPermSimple = 10000)
  fgsea_res_sorted <- fgsea_res %>% as_tibble() %>% arrange(padj, desc(NES))
  cat("GSEA analysis completed.\n")
  
  # 2.3 [Updated]: Visualize and prepare for saving
  # Create output directory
  dir.create("./figure", showWarnings = FALSE, recursive = TRUE)
  
  # Filter most significant upregulated and downregulated pathways
  top_path_up <- fgsea_res_sorted %>% filter(NES > 0) %>% slice_head(n = 1)
  top_path_down <- fgsea_res_sorted %>% filter(NES < 0) %>% slice_tail(n = 1)
  
  # Generate enrichment plots
  p_up <- if (nrow(top_path_up) > 0) {
    plotEnrichment(pathways_db[[top_path_up$pathway]], ranks) + 
      labs(title = top_path_up$pathway, 
           subtitle = sprintf("NES = %.2f, padj = %.2e", top_path_up$NES, top_path_up$padj))
  } else { NULL }
  
  p_down <- if (nrow(top_path_down) > 0) {
    plotEnrichment(pathways_db[[top_path_down$pathway]], ranks) + 
      labs(title = top_path_down$pathway, 
           subtitle = sprintf("NES = %.2f, padj = %.2e", top_path_down$NES, top_path_down$padj))
  } else { NULL }
  
  plot_list <- compact(list(p_up, p_down))
  if (length(plot_list) > 0) {
    # [Updated]: Save directly to ./figure folder
    ggsave(
      filename = paste0("./figure/GSEA_Top_Pathways_", class_name, ".pdf"), 
      plot = grid.arrange(grobs = plot_list, ncol = 1), 
      width = 8, 
      height = 5 * length(plot_list)
    )
    cat(paste0("GSEA enrichment plots saved to ./figure/\n"))
  }
  
  return(fgsea_res_sorted)
}


# --- Step 3: Main Analysis Workflow ---

# 3.1 Prepare pathway database (Hallmarks recommended)
h_gene_sets <- msigdbr(
  species = "Homo sapiens", 
  collection = "H"
) %>%
  dplyr::select(gs_name, gene_symbol) %>%
  as.data.frame() %>%
  split(.$gs_name) %>%
  lapply(function(x) x$gene_symbol)

cat("Loaded", length(h_gene_sets), "Hallmarks pathways.\n")

# 3.2 Run GSEA for "Rapid-resolving" (Class 1)
gsea_resolving <- run_gsea_analysis(
  deseq_results = results_resolving,
  class_name = "Rapid_Resolving_C1",
  pathways_db = h_gene_sets
)

# 3.3 Run GSEA for "Non-resolving" (Class 2)
gsea_non_resolving <- run_gsea_analysis(
  deseq_results = results_non_resolving,
  class_name = "Non_Resolving_C2",
  pathways_db = h_gene_sets
)

# --- [Updated]: Step 4: Merge GSEA Results into a Single Word Document ---

cat("\n\n===== Generating Word Report for GSEA Results =====\n")
# Create output directory
dir.create("./table", showWarnings = FALSE, recursive = TRUE)

# 1. Prepare data for "Rapid-resolving"
up_resolving <- gsea_resolving %>% 
  filter(NES > 0) %>% 
  slice_head(n = 10) %>% 
  dplyr::select(pathway, NES, pval, padj, size) %>% 
  mutate(across(where(is.numeric), ~ signif(., 3)))

down_resolving <- gsea_resolving %>% 
  filter(NES < 0) %>% 
  slice_tail(n = 10) %>% 
  arrange(padj, NES) %>% 
  dplyr::select(pathway, NES, pval, padj, size) %>% 
  mutate(across(where(is.numeric), ~ signif(., 3)))

# 2. Prepare data for "Non-resolving"
up_non_resolving <- gsea_non_resolving %>% 
  filter(NES > 0) %>% 
  slice_head(n = 10) %>% 
  dplyr::select(pathway, NES, pval, padj, size) %>% 
  mutate(across(where(is.numeric), ~ signif(., 3)))

down_non_resolving <- gsea_non_resolving %>% 
  filter(NES < 0) %>% 
  slice_tail(n = 10) %>% 
  arrange(padj, NES) %>% 
  dplyr::select(pathway, NES, pval, padj, size) %>% 
  mutate(across(where(is.numeric), ~ signif(., 3)))

# 3. Create and write to Word document
doc <- read_docx()
doc <- body_add_par(doc, "Supplementary Table: Dynamic GSEA Results (Day5 vs Day1)", style = "heading 1")
doc <- body_add_break(doc)

# Add results for "Rapid-resolving"
doc <- body_add_par(doc, "Rapid-Resolving Phenotype (Class 1)", style = "heading 2")
if (nrow(up_resolving) > 0) {
  doc <- body_add_par(doc, "Top 10 Upregulated Pathways", style = "heading 3")
  doc <- body_add_flextable(doc, flextable(up_resolving) %>% autofit())
}
if (nrow(down_resolving) > 0) {
  doc <- body_add_par(doc, "Top 10 Downregulated Pathways", style = "heading 3")
  doc <- body_add_flextable(doc, flextable(down_resolving) %>% autofit())
}
doc <- body_add_break(doc)

# Add results for "Non-resolving"
doc <- body_add_par(doc, "Non-Resolving Phenotype (Class 2)", style = "heading 2")
if (nrow(up_non_resolving) > 0) {
  doc <- body_add_par(doc, "Top 10 Upregulated Pathways", style = "heading 3")
  doc <- body_add_flextable(doc, flextable(up_non_resolving) %>% autofit())
}
if (nrow(down_non_resolving) > 0) {
  doc <- body_add_par(doc, "Top 10 Downregulated Pathways", style = "heading 3")
  doc <- body_add_flextable(doc, flextable(down_non_resolving) %>% autofit())
}

# 4. Save Word document
output_doc_file <- "./table/Supp_Table_GSEA_Dynamic_Results.docx"
print(doc, target = output_doc_file)
cat("GSEA results