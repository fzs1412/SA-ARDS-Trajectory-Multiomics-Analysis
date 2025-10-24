
# -------------------------- Key Step: Ensure All Required Packages Are Loaded --------------------------
required_packages <- c("dplyr", "readxl", "tidyr", "tibble", "limma", "data.table")  # Add limma and data.table
for (pkg in required_packages) {
  if (!requireNamespace(pkg, quietly = TRUE)) {  # Check if package is installed
    if (pkg == "limma") {
      # limma is a Bioconductor package; install via BiocManager
      if (!requireNamespace("BiocManager", quietly = TRUE)) {
        install.packages("BiocManager")
      }
      BiocManager::install("limma")
    } else {
      install.packages(pkg)  # Install other packages via standard method
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

# Rename final sample column to 'FinalSampleName' (adjust based on actual Excel column name; assume merged column is SampleID.y here)
# If unsure about merged column names, run: colnames(olink_mapped) to check all column names
olink_mapped <- olink_mapped %>% rename(FinalSampleName = SampleID.y)
cat("Number of unique matched samples: ", length(unique(olink_mapped$FinalSampleName)), "\n\n")


# -------------------------- Step 3: Data Reshaping (Corrected) --------------------------
# 1. Check if merged data is empty
if (nrow(olink_mapped) == 0) {
  stop("Error: No matching sample IDs found. Please verify raw data based on diagnostic results!")
}

# 2. Reshape from long to wide format (dplyr is loaded, so select() works normally)
olink_data_wide <- olink_mapped %>%
  dplyr::select(Assay, FinalSampleName, NPX) %>%  # Explicitly use dplyr::select to avoid package conflicts
  dplyr::filter(!is.na(FinalSampleName) & !is.na(Assay) & !is.na(NPX)) %>%
  dplyr::group_by(Assay, FinalSampleName) %>%
  dplyr::summarise(NPX = mean(NPX, na.rm = TRUE), .groups = "drop") %>%  # .groups="drop" clears grouping residue
  tidyr::pivot_wider(names_from = FinalSampleName, values_from = NPX) %>%  # Explicitly use tidyr::pivot_wider
  tibble::column_to_rownames("Assay")  # Explicitly use tibble::column_to_rownames


# 1. Load clinical data (Class 2 = 36 samples, Class 1 = 53 samples; filter for Class 1 or 2)
coldata_raw <- fread("./data/data331day1mice.csv") %>%
  dplyr::select(-no) %>%
  mutate(class = as.factor(class)) %>%
  tibble::column_to_rownames("SampleName")

# 2. Filter Day 1 samples
day1_samples_olink <- colnames(olink_data_wide)[!grepl("d3$|d5$", colnames(olink_data_wide))]
olink_data_d1 <- olink_data_wide[, day1_samples_olink]

# 3. Identify common samples
common_samples_prot <- intersect(colnames(olink_data_d1), rownames(coldata_raw))
cat("[Checkpoint 3]: Found ", length(common_samples_prot), " common samples between reshaped Olink Day1 data and clinical data.\n\n")

# -------------------------- Step 3: Prepare Final Data for LASSO --------------------------
cat("--- Section 3: Prepare Data for LASSO Analysis ---\n")

# 3.1 Integrate data
prot_data_t <- t(olink_data_d1) %>% as.data.frame()
common_samples_all <- intersect(rownames(coldata_raw), rownames(prot_data_t))
coldata_for_model <- coldata_raw[common_samples_all, ]
prot_data_for_model <- prot_data_t[common_samples_all, ]

stopifnot(all(rownames(prot_data_for_model) == rownames(coldata_for_model)))

# 3.2 Create outcome variable
stopifnot(all(rownames(prot_data_for_model) == rownames(coldata_for_model)))
full_modeling_data <- cbind(coldata_for_model, prot_data_for_model)
cat("Successfully integrated all clinical and protein data. Total samples: ", nrow(full_modeling_data), "\n")

# 3.2 Create outcome variable
full_modeling_data$outcome <- ifelse(full_modeling_data$class == "2", 1, 0)
cat("Outcome variable distribution (1=Class 2, 0=Others):", paste(names(table(full_modeling_data$outcome)), table(full_modeling_data$outcome), collapse = ", "), "\n")

coldata_for_model$outcome <- ifelse(coldata_for_model$class == "2", 1, 0)
coldata_for_model <- coldata_for_model %>%
  select(-class)
prot_data_for_model <- full_modeling_data[, 75:438]

# Check data structure
str(coldata_for_model)
names(coldata_for_model)
table(coldata_for_model$outcome)


# -------------------------- Step 4: LASSO Feature Selection --------------------------
library(glmnet)
library(ggplot2)

# Function: Perform LASSO and return selected variables
perform_lasso <- function(data, outcome_col = ncol(data)) {
  # Extract predictor variables and outcome variable
  X <- as.matrix(data[, -outcome_col])
  y <- as.vector(data[, outcome_col])
  
  # Perform cross-validated LASSO
  cv_lasso <- cv.glmnet(X, y, alpha = 1, family = "binomial")
  
  # Build final model using optimal lambda
  best_lambda <- cv_lasso$lambda.min
  final_model <- glmnet(X, y, alpha = 1, lambda = best_lambda, family = "binomial")
  
  # Extract variables with non-zero coefficients
  coefficients <- coef(final_model)
  non_zero_vars <- rownames(coefficients)[which(coefficients != 0)][-1]  # Exclude intercept term
  
  # Visualize cross-validation results
  plot(cv_lasso, main = "LASSO Cross-Validation")
  
  return(list(selected_vars = non_zero_vars, model = final_model, cv_result = cv_lasso))
}

# LASSO for clinical feature selection
cat("\n--- Performing LASSO for Clinical Features ---\n")
clinical_data <- coldata_for_model
clinical_lasso <- perform_lasso(clinical_data)
cat("Clinical variables selected by LASSO:", paste(clinical_lasso$selected_vars, collapse = ", "), "\n")
## Example output: "mapmin", "pf", "rrmax", "k", "paco", "urine"

# LASSO for protein feature selection
cat("\n--- Performing LASSO for Protein Features ---\n")
protein_data <- cbind(prot_data_for_model, outcome = full_modeling_data$outcome)
protein_lasso <- perform_lasso(protein_data)
cat("Protein variables selected by LASSO:", paste(protein_lasso$selected_vars, collapse = ", "), "\n")
## Example output: "AGER", "NTF3", "PRDX5", "PSIP1"

# LASSO for integrated feature selection (clinical + protein)
cat("\n--- Performing LASSO for Integrated Features ---\n")
integrated_data <- full_modeling_data %>% select(-class)  # Ensure original 'class' column is removed
integrated_lasso <- perform_lasso(integrated_data)
cat("Integrated variables selected by LASSO:", paste(integrated_lasso$selected_vars, collapse = ", "), "\n")
## Example output: "paco", "urine", "fluidout", "pf", "AGER", "NTF3", "PRDX5", "PSIP1", "HCLS1", "ITGB6"


# ==============================================================================
# Olink Proteomics Advanced Analysis Plan - Final Model Building & Comparison
# ==============================================================================

# -------------------------- Step 7: Model Building & Performance Comparison (No Data Splitting) --------------------------
# Core Objective: Build models using pre-selected variables on all available data and evaluate goodness-of-fit

cat("\n--- Section 7: Start Final Model Building & Evaluation ---\n")

# 7.1 Install and load required packages
modeling_packages <- c("pROC", "caret", "dplyr", "glmnet", "tidyr")
for (pkg in modeling_packages) {
  if (!requireNamespace(pkg, quietly = TRUE)) {
    install.packages(pkg)
  }
  library(pkg, character.only = TRUE)
}

# 7.2 Define pre-selected variable signatures
clinical_signature <- clinical_lasso$selected_vars  # Obtained from LASSO results
protein_signature <- protein_lasso$selected_vars    # Obtained from LASSO results
integrated_signature <- integrated_lasso$selected_vars  # Obtained from LASSO results

# 7.3 Set up repeated cross-validation
set.seed(123)  # Ensure reproducibility
cv_folds <- createMultiFolds(full_modeling_data$outcome, k = 10, times = 10)

# 7.4 Initialize lists to store validation results
results_a <- data.frame(pred = numeric(), obs = integer())
results_b <- data.frame(pred = numeric(), obs = integer())
results_c <- data.frame(pred = numeric(), obs = integer())
auc_values <- list(
  model_a = numeric(),
  model_b = numeric(),
  model_c = numeric()
)

cat("--- Starting 10-fold Cross-Validation with 10 Repeats ---\n")
for (i in 1:length(cv_folds)) {
  # Split data
  train_indices <- cv_folds[[i]]
  train_data <- full_modeling_data[train_indices, ]
  test_data <- full_modeling_data[-train_indices, ]
  
  # --- Model A (Clinical-only Marker Model) ---
  if (length(clinical_signature) > 0 && all(clinical_signature %in% colnames(train_data))) {
    formula_a <- as.formula(paste("outcome ~", paste(clinical_signature, collapse = "+")))
    model_a <- glm(formula_a, data = train_data, family = "binomial")
    pred_a <- predict(model_a, newdata = test_data, type = "response")
    results_a <- rbind(results_a, data.frame(pred = pred_a, obs = test_data$outcome))
    
    # Calculate AUC
    if (length(unique(test_data$outcome)) == 2) {
      current_roc <- roc(test_data$outcome, pred_a, quiet = TRUE)
      auc_values$model_a <- c(auc_values$model_a, auc(current_roc))
    }
  }
  
  # --- Model B (Protein-only Marker Model) ---
  if (length(protein_signature) > 0 && all(protein_signature %in% colnames(train_data))) {
    formula_b <- as.formula(paste("outcome ~", paste(protein_signature, collapse = "+")))
    model_b <- glm(formula_b, data = train_data, family = "binomial")
    pred_b <- predict(model_b, newdata = test_data, type = "response")
    results_b <- rbind(results_b, data.frame(pred = pred_b, obs = test_data$outcome))
    
    # Calculate AUC
    if (length(unique(test_data$outcome)) == 2) {
      current_roc <- roc(test_data$outcome, pred_b, quiet = TRUE)
      auc_values$model_b <- c(auc_values$model_b, auc(current_roc))
    }
  }
  
  # --- Model C (Integrated Marker Model) ---
  if (length(integrated_signature) > 0 && all(integrated_signature %in% colnames(train_data))) {
    formula_c <- as.formula(paste("outcome ~", paste(integrated_signature, collapse = "+")))
    model_c <- glm(formula_c, data = train_data, family = "binomial")
    pred_c <- predict(model_c, newdata = test_data, type = "response")
    results_c <- rbind(results_c, data.frame(pred = pred_c, obs = test_data$outcome))
    
    # Calculate AUC
    if (length(unique(test_data$outcome)) == 2) {
      current_roc <- roc(test_data$outcome, pred_c, quiet = TRUE)
      auc_values$model_c <- c(auc_values$model_c, auc(current_roc))
    }
  }
  
  # Print progress
  if (i %% 10 == 0) cat("Completed ", i, "/", length(cv_folds), " folds...\n")
}
cat("--- Cross-Validation Completed ---\n\n")


# 7.5 Model Performance Comparison & Visualization (Based on Cross-Validation Results)
cat("--- Performing Model Performance Comparison ---\n")
dir.create("./results_olink", showWarnings = FALSE) 
tiff("./results_olink/Figure_Final_ROC_Comparison_CrossValidated.tiff", 
     width = 8, height = 8, units = "in", res = 300, compression = "lzw")

# Function to plot ROC with confidence interval
plot_roc_with_ci <- function(roc_data, color, label, add = FALSE) {
  roc_obj <- roc(roc_data$obs, roc_data$pred, quiet = TRUE)
  plot(roc_obj, col = color, lwd = 2, add = add,
       main = ifelse(!add, "Cross-Validated ROC Curve Comparison", ""),
       xlab = "1 - Specificity", ylab = "Sensitivity")
  return(roc_obj)
}

# Plot reference line
plot(NA, xlim = c(0, 1), ylim = c(0, 1), 
     xlab = "1 - Specificity", ylab = "Sensitivity",
     main = "Cross-Validated ROC Curve Comparison with 95% CI")
abline(a = 0, b = 1, lty = 2, col = "gray")

# Plot ROC curves for each model
legend_text <- c()
legend_colors <- c()

if (nrow(results_a) > 0) {
  roc_a <- plot_roc_with_ci(results_a, "blue", "Clinical Model", add = TRUE)
  ci_a <- ci.auc(roc_a)
  legend_text <- c(legend_text, sprintf("Clinical: AUC=%.3f (%.3f-%.3f)",
                                        auc(roc_a), ci_a[1], ci_a[3]))
  legend_colors <- c(legend_colors, "blue")
}

if (nrow(results_b) > 0) {
  roc_b <- plot_roc_with_ci(results_b, "darkorange", "Protein Model", add = TRUE)
  ci_b <- ci.auc(roc_b)
  legend_text <- c(legend_text, sprintf("Protein: AUC=%.3f (%.3f-%.3f)",
                                        auc(roc_b), ci_b[1], ci_b[3]))
  legend_colors <- c(legend_colors, "darkorange")
}

if (nrow(results_c) > 0) {
  roc_c <- plot_roc_with_ci(results_c, "green", "Integrated Model", add = TRUE)
  ci_c <- ci.auc(roc_c)
  legend_text <- c(legend_text, sprintf("Integrated: AUC=%.3f (%.3f-%.3f)",
                                        auc(roc_c), ci_c[1], ci_c[3]))
  legend_colors <- c(legend_colors, "green")
}

# Add legend
legend("bottomright", legend = legend_text, col = legend_colors, lwd = 2, bty = "n")
dev.off()

cat("Cross-validated ROC comparison plot saved to 'results_olink/Figure_Final_ROC_Comparison_CrossValidated.tiff'\n")

# 7.6 Statistical Comparison of Models (Based on Cross-Validated ROC Objects)
cat("\n【Model AUC and 95% Confidence Interval】\n")
if (length(auc_values$model_a) > 0) {
  ci_a <- calculate_ci(auc_values$model_a)
  cat(sprintf("Model A (Clinical): Mean AUC = %.3f, 95%%CI = [%.3f, %.3f], Based on %d valid folds\n", 
              ci_a$mean, ci_a$lower, ci_a$upper, length(auc_values$model_a)))
} else {
  cat("Model A (Clinical): No valid AUC data\n")
}

if (length(auc_values$model_b) > 0) {
  ci_b <- calculate_ci(auc_values$model_b)
  cat(sprintf("Model B (Protein): Mean AUC = %.3f, 95%%CI = [%.3f, %.3f], Based on %d valid folds\n", 
              ci_b$mean, ci_b$lower, ci_b$upper, length(auc_values$model_b)))
} else {
  cat("Model B (Protein): No valid AUC data\n")
}

if (length(auc_values$model_c) > 0) {
  ci_c <- calculate_ci(auc_values$model_c)
  cat(sprintf("Model C (Integrated): Mean AUC = %.3f, 95%%CI = [%.3f, %.3f], Based on %d valid folds\n", 
              ci_c$mean, ci_c$lower, ci_c$upper, length(auc_values$model_c)))
} else {
  cat("Model C (Integrated): No valid AUC data\n")
}

if (exists("roc_c") && exists("roc_a")) {
  roc_test_cv <- roc.test(roc_c, roc_a)
  cat("\n【Core Comparison】Statistical comparison of cross-validated AUC between Model C (Integrated) vs. Model A (Clinical):\n")
  print(roc_test_cv)
}

cat("\n***** The final predictive model development and validation process has been fully executed *****\n")