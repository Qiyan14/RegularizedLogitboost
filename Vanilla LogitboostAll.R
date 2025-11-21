# Load libraries
library(caret)
library(doParallel)
library(writexl)
library(readr)
library(dplyr)
library(glue)
library(pROC)
library(caTools)
library(parallelly)  
library(foreach)

# Convert outcome variable to factor for classification
data.dfPH$condL <- factor(data.dfPH$condL)
levels(data.dfPH$condL) <- make.names(levels(data.dfPH$condL))

# Use all biomarkers in the model
all_biomarkers <- biomarkers
print(paste("Using all", length(all_biomarkers), "biomarkers in the model"))

# Custom summary function using ROC analysis (for cross-validation)
customSummary <- function(data, lev = NULL, model = NULL) {
  if (length(lev) > 2) {
    stop(paste("Your outcome has", length(lev), "levels. The twoClassSummary() function isn't appropriate."))
  }
  
  requireNamespaceQuietStop <- function(pkg) {
    if (!requireNamespace(pkg, quietly = TRUE)) {
      stop(paste("Package", pkg, "needed for this function to work. Please install it."))
    }
  }
  requireNamespaceQuietStop("pROC")
  
  if (!all(levels(data[, "pred"]) == lev)) {
    stop("levels of observed and predicted data do not match")
  }
  
  sens <- sensitivity(data[, "pred"], data[, "obs"], lev[1])
  spec <- specificity(data[, "pred"], data[, "obs"], lev[2])
  
  accuracy <- mean(c(sens, spec))
  
  out <- c(sens, spec, accuracy)
  names(out) <- c("Sens", "Spec", "Accuracy")
  out
}

# Define control for repeated 7-fold cross-validation (7 repeats)
control <- trainControl(
  method = "repeatedcv",
  number = 7,
  repeats = 5,
  summaryFunction = customSummary,
  classProbs = TRUE,
  allowParallel = TRUE,
  savePredictions = "final"
)

# Define iterations to test - this is the only parameter for regular LogitBoost
iterations_to_test <- c(5, 10, 15, 25, 50, 75, 100, 150, 200, 500)

# Define base columns for results - simplified for regular LogitBoost
base_cols <- c("Accuracy", "Sensitivity", "Specificity", "Precision", 
               "NPV", "Threshold", "AUC", "nIter")

# Function to create results dataframe for all biomarkers model
create_all_biomarkers_results_df <- function(iterations_to_test) {
  total_rows <- length(iterations_to_test)
  
  # Create dataframe with model configuration and performance metrics
  df <- data.frame(matrix(NA, nrow = total_rows, 
                          ncol = length(base_cols) + 1),  # +1 for model description
                   stringsAsFactors = FALSE)
  colnames(df) <- c("Model_Description", base_cols)
  
  # Set column types
  df$Model_Description <- as.character(df$Model_Description)
  df[base_cols] <- lapply(df[base_cols], as.numeric)
  
  return(df)
}

# Calculate total number of iterations for progress tracking
total_iterations <- length(iterations_to_test)
print(paste("Total iterations to process:", total_iterations))

# Check if confounders are in the dataset
valid_confounders <- confounders[confounders %in% colnames(data.dfPH)]
if (length(valid_confounders) < length(confounders)) {
  missing_confounders <- setdiff(confounders, valid_confounders)
  warning(glue("The following confounders were not found in the dataset and will be ignored: {paste(missing_confounders, collapse=', ')}"))
  confounders <- valid_confounders
}

# Set up parallel processing using parallelly (more robust)
numCores <- parallelly::availableCores() - 1
cl <- parallelly::makeClusterPSOCK(numCores)
registerDoParallel(cl)

# Export necessary functions and data to the cluster
clusterEvalQ(cl, {
  library(caret)
  library(pROC)
  library(dplyr)
  library(glue)
})

# Export the custom functions and necessary variables
clusterExport(cl, c("customSummary", "data.dfPH", "all_biomarkers", 
                    "confounders", "control"))

# Prepare data with all biomarkers and confounders
subset_data <- data.dfPH[, c(all_biomarkers, confounders, "condL")]

# Create formula with all biomarkers and confounders
formula_str <- paste("condL ~", paste(c(all_biomarkers, confounders), collapse = " + "))
formula_obj <- as.formula(formula_str)

print(paste("Model formula:", formula_str))
print("Starting LogitBoost analysis with all biomarkers and confounders...")

# Process all iterations in parallel
all_results <- foreach(iter_val = iterations_to_test, .combine = 'rbind',
                       .packages = c("caret", "pROC", "dplyr", "glue")) %dopar% {
                         
                         # Try to fit model with this iteration parameter
                         tryCatch({
                           # Create a row for results
                           result_row <- data.frame(matrix(NA, nrow = 1, 
                                                           ncol = 1 + length(base_cols)),
                                                    stringsAsFactors = FALSE)
                           colnames(result_row) <- c("Model_Description", base_cols)
                           
                           # Store model description
                           result_row[1, "Model_Description"] <- paste0("All_", length(all_biomarkers), "_Biomarkers_with_Confounders")
                           
                           # Store iteration parameter
                           result_row[1, "nIter"] <- iter_val
                           
                           # Create tuning parameter grid for caret (just nIter for regular LogitBoost)
                           tuning_param <- data.frame(nIter = iter_val)
                           
                           # Train LogitBoost model (using built-in "LogitBoost" method)
                           set.seed(123)
                           logitboost_model <- train(
                             formula_obj,
                             data = subset_data,
                             method = "LogitBoost",  # Standard LogitBoost method in caret
                             tuneGrid = tuning_param,
                             metric = "Accuracy",
                             trControl = control
                           )
                           
                           # Get the factor levels from the outcome
                           lev <- levels(subset_data$condL)
                           pred_subset <- logitboost_model$pred
                           
                           if (length(unique(pred_subset$obs)) < 2 || !all(lev %in% colnames(pred_subset))) {
                             results <- rep(NA, 7)  # 7 is the number of performance metrics
                           } else {
                             roc_obj <- suppressMessages(roc(pred_subset$obs, pred_subset[[lev[2]]], quiet = TRUE))
                             
                             Metrics <- coords(roc_obj, "best", best.method = "youden", 
                                               ret = c("threshold", "accuracy", "sensitivity", "specificity", "precision", "npv"))
                             
                             results <- c(
                               Metrics["accuracy"],
                               Metrics["sensitivity"],
                               Metrics["specificity"],
                               Metrics["precision"],
                               Metrics["npv"],
                               Metrics["threshold"],
                               as.numeric(roc_obj$auc)
                             )
                           }
                           
                           # Store performance metrics
                           result_row[1, c("Accuracy", "Sensitivity", 
                                           "Specificity", "Precision", "NPV", 
                                           "Threshold", "AUC")] <- results
                           
                           return(result_row)
                           
                         }, error = function(e) {
                           warning(glue("Error processing all biomarkers with nIter={iter_val}: {e$message}"))
                           return(NULL)
                         })
                       }

# Stop the cluster and register sequential processing
tryCatch({
  parallel::stopCluster(cl)
}, error = function(e) {
  message("Error in stopping the cluster: ", e$message)
})

registerDoSEQ()

# Round all numeric values in results
if (!is.null(all_results)) {
  all_results <- all_results %>% 
    mutate_if(is.numeric, ~ round(., 5))
}

# Prepare export list
export_list <- list()

if (!is.null(all_results) && nrow(all_results) > 0) {
  # Sort by nIter for consistency
  all_results <- all_results %>% arrange(nIter)
  
  # Split into chunks if necessary (though unlikely with just iterations)
  if (nrow(all_results) > 1000000) {
    chunks <- split(all_results, ceiling(seq_len(nrow(all_results)) / 1000000))
    for (j in seq_along(chunks)) {
      export_list[[paste0("Results_All_Biomarkers_Part", j)]] <- chunks[[j]]
    }
  } else {
    export_list[["Results_All_Biomarkers"]] <- all_results
  }
  
  # Add a summary sheet with biomarker list
  biomarker_summary <- data.frame(
    Biomarker_Number = 1:length(all_biomarkers),
    Biomarker_Name = all_biomarkers,
    stringsAsFactors = FALSE
  )
  export_list[["Biomarker_List"]] <- biomarker_summary
  
  # Add confounders list if any
  if (length(confounders) > 0) {
    confounder_summary <- data.frame(
      Confounder_Number = 1:length(confounders),
      Confounder_Name = confounders,
      stringsAsFactors = FALSE
    )
    export_list[["Confounder_List"]] <- confounder_summary
  }
  
  # Add parameter summary
  param_summary <- data.frame(
    Parameter = c("Total_Biomarkers", "Total_Confounders", "CV_Folds", "CV_Repeats", "Iterations_Tested"),
    Value = c(length(all_biomarkers), length(confounders), 7, 7, length(iterations_to_test)),
    stringsAsFactors = FALSE
  )
  export_list[["Analysis_Parameters"]] <- param_summary
  
  # Write all data frames to an Excel file
  write_xlsx(export_list, path = "Results_StandardLogitBoost_AllBiomarkersSim.10.3.25.xlsx")
  
  print(paste("Analysis complete! Results saved for", length(all_biomarkers), "biomarkers"))
  print(paste("Total iterations tested:", nrow(all_results)))
  print("Results summary:")
  print(all_results)
} else {
  print("No results to export - check for errors in model fitting")
}

# Print final summary

print("ANALYSIS SUMMARY")

print(paste("Total biomarkers used:", length(all_biomarkers)))
print(paste("Total confounders used:", length(confounders)))
print(paste("Cross-validation: 7-fold repeated 7 times"))
print(paste("Method: Standard LogitBoost"))
print(paste("Iterations tested:", paste(iterations_to_test, collapse = ", ")))
