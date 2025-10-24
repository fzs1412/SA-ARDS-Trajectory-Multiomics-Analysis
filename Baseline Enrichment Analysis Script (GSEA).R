```r
# ==============================================================================
# Baseline Enrichment Analysis Script (GSEA)
# Using GSEA to analyze pathway differences between "Non-resolving" (class 2) vs "Rapid-resolving" (class 1)
# ==============================================================================

# --- Step 1: Load Required R Packages ---
# Ensure all required packages are installed
packages_to_load <- c(
  "dplyr", "tibble", "data.table", "clusterProfiler", "org.Hs.eg.db",
  "fgsea", "ggplot2", "msigdbr", "gridExtra", "officer", "flextable",
  "Cairo", "purrr"  # purrr package provides the compact() function
)

for (pkg in packages_to_load) {
  if (!requireNamespace(pkg, quietly = TRUE)) {
    if (pkg %in% c("clusterProfiler", "org.Hs.eg.db", "msigdbr")) {
      if (!requireNamespace("BiocManager", quietly = TRUE)) install.packages("BiocManager")
      BiocManager::install(pkg)
    } else {
      install.packages(pkg)
    }
  }
  library(pkg, character.only = TRUE)
}

# Set network timeout to prevent msigdbr data download failures
options(timeout = 300)
# Increase expression nesting limit to handle complex pathway names
options(expressions = 5e5)

# Set your main working directory
setwd('d:/OneDrive/R/25SXFX/25ARDS')


# --- Step 2: Prepare Pre-ranked Gene List for GSEA ---

# 2.1 Load previously saved Day 1 DESeq2 analysis results
results_file <- "./data/baseline_DESeq2_results_object.rds"
if (!file.exists(results_file)) {
  stop("Error: '", results_file, "' file not found. Please ensure the Step 1 baseline analysis script ran successfully.")
}
res <- readRDS(results_file)
cat("Successfully loaded DESeq2 results file.\n")

# 2.2 Convert results to data frame and remove NA values
res_df <- as.data.frame(res) %>%
  rownames_to_column("ENSEMBL") %>%
  filter(!is.na(stat))  # GSEA uses the 'stat' column for ranking

# 2.3 Convert ENSEMBL IDs to ENTREZ IDs
ensembl_to_entrez <- bitr(res_df$ENSEMBL,
                          fromType = "ENSEMBL",
                          toType = "ENTREZID",
                          OrgDb = org.Hs.eg.db)

# 2.4 Create final ranked list
# This is a named vector with ENTREZ IDs as names and 'stat' values as content
ranks <- res_df %>%
  inner_join(ensembl_to_entrez, by = "ENSEMBL") %>%
  # When one ENTREZ ID maps to multiple ENSEMBL IDs, keep the one with the largest absolute 'stat' value
  group_by(ENTREZID) %>%
  slice_max(order_by = abs(stat), n = 1, with_ties = FALSE) %>%
  ungroup() %>%
  { setNames(.$stat, .$ENTREZID) } %>%
  sort(decreasing = TRUE)  # Must be sorted in descending order

cat("Successfully created ranked list containing", length(ranks), "genes for GSEA analysis.\n")


# --- Step 3: Prepare Pathway Gene Sets (Hallmarks and GO:BP recommended) ---

# Option A: Hallmarks gene sets (fewer pathways, more focused, highly recommended)
h_gene_sets <- msigdbr(species = "Homo sapiens", category = "H") %>%
  dplyr::select(gs_name, entrez_gene) %>%
  as.data.frame() %>%
  split(.$gs_name) %>%
  lapply(function(x) as.character(x$entrez_gene))  # Ensure IDs are character type

cat("Loaded", length(h_gene_sets), "Hallmarks pathways.\n")

# Option B: GO Biological Processes (GO:BP) gene sets (more pathways, comprehensive)
gobp_gene_sets <- msigdbr(species = "Homo sapiens", category = "C5", subcategory = "GO:BP") %>%
  dplyr::select(gs_name, entrez_gene) %>%
  as.data.frame() %>%
  split(.$gs_name) %>%
  lapply(function(x) as.character(x$entrez_gene))

cat("Loaded", length(gobp_gene_sets), "GO:BP pathways.\n")

# ----【Select pathway set to use here, default is Hallmarks】----
pathways_to_use <- h_gene_sets
cat("---- Currently selected pathway set: Hallmarks ----\n")
# To use GO:BP, uncomment the line below and comment out the line above
# pathways_to_use <- gobp_gene_sets
# cat("---- Currently selected pathway set: GO:BP ----\n")


# --- Step 4: Run GSEA Analysis ---
fgsea_res <- fgsea(
  pathways = pathways_to_use,
  stats    = ranks,
  minSize  = 15,
  maxSize  = 500,
  nPermSimple = 10000  # Increase permutations for more stable p-values
)

# Convert to tibble and sort
fgsea_res_sorted <- fgsea_res %>%
  as_tibble() %>%
  arrange(padj, desc(NES))

cat("--- GSEA analysis completed. Pathways with smallest adjusted p-values are shown below ---\n")
print(head(fgsea_res_sorted %>% dplyr::select(-leadingEdge, -ES), 10))


# --- Step 5: Visualization - Generate Publication-Quality Plots (Enrichment Plots) ---

# 1. Filter most significant upregulated and downregulated pathways (sorted by padj)
fgsea_df <- fgsea_res_sorted %>% filter(!is.na(padj))
top_pathway_up <- fgsea_df %>% filter(NES > 0) %>% slice_head(n = 1)
top_pathway_down <- fgsea_df %>% filter(NES < 0) %>% slice_head(n = 1)

# 2. Plot upregulated pathway
p_up <- if (nrow(top_pathway_up) > 0) {
  plotEnrichment(pathways_to_use[[top_pathway_up$pathway]], ranks) +
    labs(
      title = top_pathway_up$pathway,
      subtitle = sprintf("NES = %.2f, padj = %.2e", top_pathway_up$NES, top_pathway_up$padj)
    ) +
    theme(plot.title = element_text(face = "bold", size = 12))
} else { NULL }

# 3. Plot downregulated pathway
p_down <- if (nrow(top_pathway_down) > 0) {
  plotEnrichment(pathways_to_use[[top_pathway_down$pathway]], ranks) +
    labs(
      title = top_pathway_down$pathway,
      subtitle = sprintf("NES = %.2f, padj = %.2e", top_pathway_down$NES, top_pathway_down$padj)
    ) +
    theme(plot.title = element_text(face = "bold", size = 12))
} else { NULL }

# 4. Combine plots and save as PDF
output_plot_file <- "./figure/GSEA_Baseline_Enrichment.pdf"
plot_list <- compact(list(p_up, p_down))  # Remove potential NULLs

if (length(plot_list) > 0) {
  CairoPDF(output_plot_file, width = 8, height = 5 * length(plot_list))
  grid.arrange(grobs = plot_list, ncol = 1)
  dev.off()
  cat("GSEA enrichment plots saved to:", output_plot_file, "\n")
} else {
  cat("No significant enriched pathways found for plotting.\n")
}


# --- Step 6: Generate Supplementary Material Table (Word Document) ---

# 1. Filter top 10 upregulated and downregulated pathways
fgsea_table_data <- fgsea_res_sorted %>%
  dplyr::select(pathway, NES, pval, padj, size) %>%
  mutate(across(where(is.numeric), ~ signif(., 3)))  # Keep 3 significant digits

up_top10 <- fgsea_table_data %>% filter(NES > 0) %>% slice_head(n = 10)
down_top10 <- fgsea_table_data %>% filter(NES < 0) %>% arrange(padj) %>% slice_head(n = 10)  # Downregulated also sorted by padj

# 2. Create and write to Word document
doc <- read_docx()
doc <- body_add_par(doc, "Supplementary Table: GSEA Results for Day 1 (class 2 vs class 1)", style = "heading 1")

if (nrow(up_top10) > 0) {
  doc <- body_add_par(doc, "Top 10 Enriched Pathways (Upregulated in Non-resolving, class 2)", style = "heading 2")
  doc <- body_add_flextable(doc, flextable(up_top10) %>% autofit())
  doc <- body_add_break(doc)
}

if (nrow(down_top10) > 0) {
  doc <- body_add_par(doc, "Top 10 Enriched Pathways (Downregulated in Non-resolving, class 2)", style = "heading 2")
  doc <- body_add_flextable(doc, flextable(down_top10) %>% autofit())
}

output_doc_file <- "./table/Supp_Table_GSEA_Baseline.docx"
print(doc, target = output_doc_file)

cat("GSEA results table saved to:", output_doc_file, "\n")

cat("\n===== GSEA Baseline Enrichment Analysis Completed Successfully =====\n")
```