

# -------------------------- Key Step: Ensure All Required Packages Are Loaded --------------------------
required_packages <- c("dplyr", "readxl", "tidyr", "tibble", "limma", "data.table")  # Add limma and data.table
for (pkg in required_packages) {
  if (!requireNamespace(pkg, quietly = TRUE)) {  # Check if package is installed
    if (pkg == "limma") {
      # limma is a Bioconductor package and requires installation via BiocManager
      if (!requireNamespace("BiocManager", quietly = TRUE)) {
        install.packages("BiocManager")
      }
      BiocManager::install("limma")
    } else {
      install.packages(pkg)  # Install other packages using standard method
    }
  }
  library(pkg, character.only = TRUE)  # Load the package
}

# 1. Load Olink long-format data (NPX values)
olink_long_data <- read.delim("D:/OneDrive/R/01sepsis/OMIX006238-01.csv", sep = ";", stringsAsFactors = FALSE)
cat("Olink data loaded successfully: ", nrow(olink_long_data), "rows\n")

# 2. Load sample information mapping table
sample_info <- read_excel("D:/OneDrive/R/01sepsis/Sample_Info.xlsxOLINK.xlsx")
cat("Sample mapping table loaded successfully: ", nrow(sample_info), "rows\n")

# 3. Standardize column types for merging and remove potential whitespace
olink_long_data$SampleID <- trimws(as.character(olink_long_data$SampleID))
sample_info$Sample <- trimws(as.character(sample_info$Sample)) 

# 5. Merge data (using correct column name matching)
olink_mapped <- inner_join(olink_long_data, sample_info, by = c("SampleID" = "Sample"))
cat("[Checkpoint 1]: Successfully matched ", nrow(olink_mapped), " rows of data.\n")

# Rename final sample column to 'FinalSampleName' (adjust based on actual Excel column name; assuming merged column is SampleID.y here)
# If unsure about merged column names, run: colnames(olink_mapped) to check all column names
olink_mapped <- olink_mapped %>% rename(FinalSampleName = SampleID.y)
cat("Number of unique matched samples: ", length(unique(olink_mapped$FinalSampleName)), "\n\n")


# -------------------------- Step 3: Data Reshaping (Corrected) --------------------------
# 1. Check if merged data is empty
if (nrow(olink_mapped) == 0) {
  stop("Error: No matching sample IDs found. Please check raw data based on diagnostic results!")
}

# 2. Reshape from long to wide format (dplyr is loaded, so select() works normally)
olink_data_wide <- olink_mapped %>%
  dplyr::select(Assay, FinalSampleName, NPX) %>%  # Explicitly use dplyr::select to avoid package conflicts
  dplyr::filter(!is.na(FinalSampleName) & !is.na(Assay) & !is.na(NPX)) %>%
  dplyr::group_by(Assay, FinalSampleName) %>%
  dplyr::summarise(NPX = mean(NPX, na.rm = TRUE), .groups = "drop") %>%  # .groups="drop" clears grouping residue
  tidyr::pivot_wider(names_from = FinalSampleName, values_from = NPX) %>%  # Explicitly use tidyr::pivot_wider
  tibble::column_to_rownames("Assay")  # Explicitly use tibble::column_to_rownames


# 1. Load clinical data (focus on Class 1 vs Class 2 from data331)
coldata_raw <- fread("./data/data331.csv") %>%
  dplyr::select(-no) %>%
  filter(class == 1 | class == 2) %>%
  mutate(class = as.factor(class)) %>%
  tibble::column_to_rownames("SampleName")

# 2. Filter Day 1 samples
day1_samples_olink <- colnames(olink_data_wide)[!grepl("d3$|d5$", colnames(olink_data_wide))]
olink_data_d1 <- olink_data_wide[, day1_samples_olink]

# 3. Identify common samples
common_samples_prot <- intersect(colnames(olink_data_d1), rownames(coldata_raw))
cat("[Checkpoint 3]: Found ", length(common_samples_prot), " common samples between reshaped Olink Day1 data and clinical data.\n\n")

# -------------------------- Step 4: Final Data Alignment and Preparation --------------------------

# 4.1 Filter final protein expression matrix and clinical data based on common samples
prot_matrix_final <- olink_data_d1[, common_samples_prot]
coldata_prot_final <- coldata_raw[common_samples_prot, ]

# 4.2 Ultimate check: Ensure sample order is completely consistent
stopifnot(all(colnames(prot_matrix_final) == rownames(coldata_prot_final)))
cat("[Checkpoint 4]: Final protein matrix and clinical data successfully aligned.\n")

# 5.1 Create design matrix
# This matrix tells limma which group each sample belongs to
design <- model.matrix(~ 0 + class, data = coldata_prot_final)
colnames(design) <- c("class1", "class2")  # Rename columns to match group names

# 5.2 Fit linear model
# Olink NPX data is already log2-scaled and can be used directly in limma
fit <- lmFit(prot_matrix_final, design)

# 5.3 Create contrast matrix
# Define our comparison of interest: class1 vs class2
cont.matrix <- makeContrasts(class1 - class2, levels = design)

# 5.4 Apply contrast
fit2 <- contrasts.fit(fit, cont.matrix)

# 5.5 Calculate statistics using empirical Bayes method
fit2 <- eBayes(fit2)
cat("[Checkpoint 5]: limma differential expression analysis completed.\n\n")

# 5.6 Extract complete differential expression results table
# coef=1 indicates we only care about the first contrast (class1 - class2)
# number=Inf indicates extracting results for all proteins
top_proteins <- topTable(fit2, coef = 1, number = Inf, sort.by = "P")
top_proteins <- tibble::rownames_to_column(top_proteins, "Protein")  # Convert row names (protein names) to a column


# -------------------------- Step 6: Result Visualization (Optimized Version) --------------------------

# 6.1 Volcano Plot - Optimized Version
cat("Generating optimized volcano plot...\n")

# Set p-value and logFC thresholds
logFC_threshold_strict <- 1  # |logFC| > 1
logFC_threshold_loose <- 0.585  # |logFC| > 1.5-fold
p_value_threshold <- 0.05  # Use raw P-value as adj.P.Val may be too strict

# Proteins to label on the plot: P < 0.05 and |logFC| > 0.585
proteins_to_label <- top_proteins %>%
  filter(P.Value < p_value_threshold & abs(logFC) > logFC_threshold_loose) %>%
  pull(Protein)

p_volcano <- EnhancedVolcano(top_proteins,
                             lab = top_proteins$Protein,
                             selectLab = proteins_to_label,  # Label more potentially relevant proteins
                             x = 'logFC',
                             y = 'P.Value',  # Use P.Value for more dynamic range
                             title = 'Differentially Expressed Proteins (Class 1 vs Class 2)',
                             subtitle = 'Olink Day 1 Data',
                             pCutoff = p_value_threshold,
                             FCcutoff = logFC_threshold_loose,
                             pointSize = 3.0,
                             labSize = 4.5,
                             colAlpha = 0.8,
                             legendPosition = 'bottom',
                             drawConnectors = TRUE,
                             widthConnectors = 0.75) +
  ggplot2::labs(y = expression(-Log[10]~'(P-value)'))  # Explicit Y-axis label

print(p_volcano)


# 6.2 Heatmap - Optimized Version
cat("Generating optimized heatmap...\n")

# Filter proteins for heatmap: Select top 40 most significant proteins (sorted by p-value)
top_n_proteins <- top_proteins %>%
  arrange(P.Value) %>%
  head(40) %>%
  pull(Protein)

if (length(top_n_proteins) > 1) {
  # Extract expression data for these proteins
  heatmap_matrix <- prot_matrix_final[top_n_proteins, ]
  
  # Create column annotations (sample group information)
  annotation_col <- data.frame(
    Subtype = factor(coldata_prot_final$class, labels = c("Class 1", "Class 2")),
    row.names = rownames(coldata_prot_final)
  )
  
  # Draw heatmap
  p_heatmap <- pheatmap(heatmap_matrix,
                        scale = "row",  # Z-score normalization by row
                        annotation_col = annotation_col,
                        show_colnames = FALSE,
                        fontsize_row = 8,
                        main = "Heatmap of Top 40 Differentially Expressed Proteins")
  print(p_heatmap)
} else {
  cat("Warning: Insufficient number of proteins to generate heatmap.\n")
}


# -------------------------- Step 7: Save Results --------------------------
cat("Saving analysis results...\n")
dir.create("./results_olink", showWarnings = FALSE)  # Create results folder

# 7.1 Save differential protein list
fwrite(top_proteins, "./results_olink/DE_proteins_class1_vs_class2.csv")

# 7.2 Save volcano plot
ggsave("./results_olink/volcano_plot_class1_vs_class2.pdf", plot = p_volcano, width = 10, height = 10)
ggsave("./results_olink/volcano_plot_class1_vs_class2.tiff", plot = p_volcano, width = 10, height = 10, dpi = 300, compression = "lzw")


# 7.3 Save heatmap
if (length(top_n_proteins) > 1) {
  # Save as PDF
  pdf("./results_olink/heatmap_class1_vs_class2.pdf", width = 8, height = 10)
  pheatmap(heatmap_matrix,
           scale = "row",
           annotation_col = annotation_col,
           show_colnames = FALSE,
           fontsize_row = 8,
           main = "Heatmap of Top 40 Differentially Expressed Proteins")
  dev.off()
  
  # Save as TIFF
  tiff("./results_olink/heatmap_class1_vs_class2.tiff", width = 8, height = 10, units = "in", res = 300, compression = "lzw")
  pheatmap(heatmap_matrix,
           scale = "row",
           annotation_col = annotation_col,
           show_colnames = FALSE,
           fontsize_row = 8,
           main = "Heatmap of Top 40 Differentially Expressed Proteins")
  dev.off()
}

cat("[Complete]: All Olink proteomics analysis results saved to 'results_olink' folder.\n")
```