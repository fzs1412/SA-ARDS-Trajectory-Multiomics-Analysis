```r
# ==============================================================================
# Olink Proteomics Advanced Analysis Plan - Analysis Section
# Strategy: Strengthen validation of transcriptomics and anchor biological mechanisms
# ==============================================================================
setwd('d:/OneDrive/R/25SXFX/25ARDS')

# -------------------------- Key Step: Load All Required Packages --------------------------
required_packages <- c("dplyr", "readxl", "tidyr", "tibble", "limma", "data.table")  # Add limma and data.table
for (pkg in required_packages) {
  if (!requireNamespace(pkg, quietly = TRUE)) {  # Check if package is installed
    if (pkg == "limma") {
      # limma is a Bioconductor package, install via BiocManager
      if (!requireNamespace("BiocManager", quietly = TRUE)) {
        install.packages("BiocManager")
      }
      BiocManager::install("limma")
    } else {
      install.packages(pkg)  # Install CRAN packages via standard method
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

# 5. Merge data (use correct column name matching)
olink_mapped <- inner_join(olink_long_data, sample_info, by = c("SampleID" = "Sample"))
cat("[Checkpoint 1]: ", nrow(olink_mapped), "rows of data matched successfully.\n")

# Rename final sample column to 'FinalSampleName' (adjust based on actual Excel column name; assume merged column is SampleID.y here)
# If unsure about merged column names, run: colnames(olink_mapped) to check all column names
olink_mapped <- olink_mapped %>% rename(FinalSampleName = SampleID.y)
cat("Number of unique matched samples: ", length(unique(olink_mapped$FinalSampleName)), "\n\n")


# -------------------------- Step 3: Data Reshaping (Corrected) --------------------------
# 1. Check if merged data is empty
if (nrow(olink_mapped) == 0) {
  stop("Error: No matched sample IDs found. Please verify raw data based on diagnostic results!")
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
cat("[Checkpoint 3]: ", length(common_samples_prot), "common samples found between reshaped Olink Day1 data and clinical data.\n\n")

# -------------------------- Step 4: Final Data Alignment and Preparation --------------------------

# 4.1 Filter final protein expression matrix and clinical data based on common samples
prot_matrix_final <- olink_data_d1[, common_samples_prot]
coldata_prot_final <- coldata_raw[common_samples_prot, ]

# 4.2 Final consistency check: Ensure sample order is identical
stopifnot(all(colnames(prot_matrix_final) == rownames(coldata_prot_final)))
cat("[Checkpoint 4]: Final protein matrix and clinical data aligned successfully.\n")

# 5.1 Create design matrix
# This matrix tells limma which group each sample belongs to
design <- model.matrix(~ 0 + class, data = coldata_prot_final)
colnames(design) <- c("class1", "class2")  # Rename columns to match group names

# 5.2 Fit linear model
# Olink NPX data is already log2-scaled and can be used directly in limma
fit <- lmFit(prot_matrix_final, design)

# 5.3 Create contrast matrix
# Define the comparison of interest: class1 vs class2
cont.matrix <- makeContrasts(class1 - class2, levels = design)

# 5.4 Apply contrast
fit2 <- contrasts.fit(fit, cont.matrix)

# 5.5 Calculate statistics using empirical Bayes method
fit2 <- eBayes(fit2)
cat("[Checkpoint 5]: limma differential expression analysis completed.\n\n")

# 5.6 Extract complete differential expression results table
# coef=1 means we only care about the first contrast (class1 - class2)
# number=Inf means extract results for all proteins
top_proteins <- topTable(fit2, coef = 1, number = Inf, sort.by = "P")
top_proteins <- tibble::rownames_to_column(top_proteins, "Protein")  # Convert row names to column


# -------------------------- Section 0: Setup --------------------------
# Ensure all required packages are loaded
required_packages <- c("dplyr", "readxl", "tidyr", "tibble", "data.table", 
                       "limma", "ggplot2", "ggrepel", "pheatmap", "corrplot", "psych", "EnhancedVolcano")
for (pkg in required_packages) {
  # For Bioconductor packages, use BiocManager
  if (pkg %in% c("limma", "EnhancedVolcano") && !requireNamespace(pkg, quietly = TRUE)) {
    if (!requireNamespace("BiocManager", quietly = TRUE)) install.packages("BiocManager")
    BiocManager::install(pkg)
  } else if (!requireNamespace(pkg, quietly = TRUE)) {
    install.packages(pkg)
  }
  library(pkg, character.only = TRUE)
}

# Set working directory
setwd('d:/OneDrive/R/25SXFX/25ARDS')

# Create a dedicated output directory for this advanced analysis
dir.create("./results_olink_advanced", showWarnings = FALSE)
cat("--- Analysis setup complete. Outputs will be saved in 'results_olink_advanced' ---\n\n")


# -------------------------- Section 1: Data Preparation --------------------------
# Load and process data
olink_long_data <- read.delim("D:/OneDrive/R/01sepsis/OMIX006238-01.csv", sep = ";", stringsAsFactors = FALSE)
sample_info <- read_excel("D:/OneDrive/R/01sepsis/Sample_Info.xlsxOLINK.xlsx")
olink_long_data$SampleID <- trimws(as.character(olink_long_data$SampleID))
sample_info$Sample <- trimws(as.character(sample_info$Sample))

olink_mapped <- inner_join(olink_long_data, sample_info, by = c("SampleID" = "Sample")) %>%
  rename(FinalSampleName = SampleID.y)

olink_data_wide <- olink_mapped %>%
  dplyr::select(Assay, FinalSampleName, NPX) %>%
  dplyr::filter(!is.na(FinalSampleName) & !is.na(Assay) & !is.na(NPX)) %>%
  dplyr::group_by(Assay, FinalSampleName) %>%
  dplyr::summarise(NPX = mean(NPX, na.rm = TRUE), .groups = "drop") %>%
  tidyr::pivot_wider(names_from = FinalSampleName, values_from = NPX) %>%
  tibble::column_to_rownames("Assay")

# In this analysis, we compare class 1 (Rapid-resolving) vs class 3 (Non-resolving)
coldata_raw <- fread("./data/data331day1mice.csv") %>%
  dplyr::select(-no, -ID) %>%
  filter(class == 1 | class == 3) %>%  # Filtering for Class 1 and 3
  mutate(class = factor(class, levels = c("1", "3"))) %>% 
  tibble::column_to_rownames("SampleName")

day1_samples_olink <- colnames(olink_data_wide)[!grepl("d3$|d5$", colnames(olink_data_wide))]
olink_data_d1 <- olink_data_wide[, day1_samples_olink]
common_samples_prot <- intersect(colnames(olink_data_d1), rownames(coldata_raw))
prot_matrix_final <- olink_data_d1[, common_samples_prot]
coldata_prot_final <- coldata_raw[common_samples_prot, ]

# Final consistency check
stopifnot(all(colnames(prot_matrix_final) == rownames(coldata_prot_final)))
cat("--- Section 1: Data preparation and alignment complete. N =", length(common_samples_prot), " ---\n\n")


# -------------------------- Section 2: Full Differential Expression Analysis (Foundation) --------------------------
# This limma analysis is the foundation for all subsequent steps.

# 2.1 Create design matrix
# We use class 1 as the baseline reference group
coldata_prot_final$class <- factor(coldata_prot_final$class, levels = c("1", "3"))
design <- model.matrix(~ class, data = coldata_prot_final)
# R will automatically create two columns: '(Intercept)' for class 1, and 'class3' for the effect of class 3.
# We must use 'class3' as the coefficient to get the desired comparison.
cat("Design matrix columns:", colnames(design), "\n")

# 2.2 Fit linear model
fit <- lmFit(prot_matrix_final, design)
fit2 <- eBayes(fit)

# 2.3 Extract complete differential expression results table
# CORRECTED: Use the correct coefficient name 'class3' from the design matrix.
# Using a wrong name like "class" would cause a 'subscript out of bounds' error.
# A positive logFC value means the protein is more highly expressed in "Non-resolving" (class 3) than in "Rapid-resolving" (class 1).
full_results <- topTable(fit2, coef = "class3", number = Inf, sort.by = "P") %>%
  tibble::rownames_to_column("Protein")

# 2.4 Create and save results directory
dir.create("./results_olink_advanced", showWarnings = FALSE)
fwrite(full_results, "./results_olink_advanced/Supp_Table_Full_DE_Protein_Results_C3vsC1.csv")

cat("--- Section 2: Basic differential expression analysis completed ---\n")
cat("Top protein by p-value:", full_results$Protein[1], 
    " (logFC =", round(full_results$logFC[1], 2), 
    ", P.Value =", format.pval(full_results$P.Value[1], digits = 2), ")\n\n")


# -------------------------- Section 3: Goal 1 - Validate Key Transcriptomic Pathways (Core Step) --------------------------
cat("--- Section 3: Execute core task - Hypothesis-driven pathway validation ---\n")

# 3.1 Define key pathways and proteins to validate
# !!! Note: This list is pre-populated based on your previous GSEA results (Interferon, IL-6 pathways, etc.)
# Please check if these proteins are present in your Olink panel and add/remove as needed
pathway_protein_map <- list(
  `Interferon_Response` = c("CXCL9", "CXCL10", "CXCL11", "IDO1", "STAT1", "GBP2"),
  `IL6_JAK_STAT3_Signaling` = c("IL6", "IL6R", "SOCS3", "STAT3", "OSM"),
  `Neutrophil_Related` = c("MMP8", "ELANE", "MPO", "CXCL1", "CXCL8", "CXCL6")
)

# 3.2 Filter differential results for these target proteins
validation_proteins <- unlist(pathway_protein_map)
validation_results <- full_results %>% filter(Protein %in% validation_proteins)

# 3.3 Visualize validation results
if (nrow(validation_results) > 0) {
  # Create matrix for heatmap
  validation_heatmap_matrix <- validation_results %>%
    dplyr::select(Protein, logFC, P.Value) %>%
    mutate(`-log10(P.Value)` = -log10(P.Value)) %>%
    dplyr::select(Protein, logFC, `-log10(P.Value)`) %>%
    tibble::column_to_rownames("Protein")
  
  # Generate heatmap object
  p_validation_heatmap <- pheatmap(
    validation_heatmap_matrix,
    main = "Validation of Key Pathways at Protein Level\n(Non-resolving vs. Rapid-resolving)",
    cluster_cols = FALSE,
    fontsize_row = 10,
    angle_col = 0,
    legend = TRUE,
    display_numbers = TRUE,
    number_format = "%.2f",
    color = colorRampPalette(c("blue", "white", "red"))(50)
  )
  
  # Save heatmap as TIFF format
  tiff("./results_olink_advanced/Figure_A_Pathway_Validation_Heatmap.tiff", width = 6, height = 8, units = "in", res = 300, compression = "lzw")
  print(p_validation_heatmap)
  dev.off()
  
  cat("Pathway validation heatmap saved. The plot shows logFC and significance levels of preselected proteins.\n")
} else {
  cat("Warning: None of the validation proteins specified in the list were found in Olink results. Please check protein names.\n")
}
cat("--- Core Task 1 completed ---\n\n")


# -------------------------- Section 4: Goal 2 - Prepare Data for Network Analysis --------------------------
cat("--- Section 4: Execute Task 2 - Screen candidate proteins for network analysis ---\n")

# Filter proteins with p-value < 0.1 as candidates (captures proteins with differential trends)
network_candidates <- full_results %>%
  filter(P.Value < 0.1) %>%
  dplyr::select(Protein, logFC)

# Save list for subsequent upload to STRING database
fwrite(network_candidates, "./results_olink_advanced/For_STRING_DB_Network_Analysis_p_le_0.1.csv")

cat(nrow(network_candidates), "candidate proteins (P.Value < 0.1) saved.\n")
cat("Next step: Upload the 'For_STRING_DB_Network_Analysis_p_le_0.1.csv' file to the STRING database website (string-db.org) for protein-protein interaction network analysis.\n")
cat("--- Task 2 completed ---\n\n")


# -------------------------- Section 4: Goal 2 - Identify Core Biological Mechanisms (Protein Network Analysis) --------------------------
cat("--- Section 4: Execute Task 2 - Protein-Protein Interaction Network Analysis ---\n")

# 4.1 Install and load required network analysis packages
network_packages <- c("STRINGdb", "igraph", "ggraph")
for (pkg in network_packages) {
  if (!requireNamespace(pkg, quietly = TRUE)) {
    if (!requireNamespace("BiocManager", quietly = TRUE)) install.packages("BiocManager")
    BiocManager::install(pkg)
  }
  library(pkg, character.only = TRUE)
}

# 4.2 Filter candidate proteins for network analysis (P.Value < 0.1)
network_candidates <- full_results %>%
  filter(P.Value < 0.1)

if (nrow(network_candidates) > 2) {
  # 4.3 Initialize STRINGdb object
  # score_threshold=400 (medium confidence), species=9606 (Homo sapiens)
  string_db <- STRINGdb$new(version = "11.5", species = 9606, score_threshold = 