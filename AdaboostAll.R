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
library(ada)  

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

# Define iterations to test 
iter_to_test <- c(10, 25, 50, 100, 200, 500)

# Define comprehensive tuning parameters for AdaBoost
tuning_params <- expand.grid(
  iter = iter_to_test,           # Number of boosting iterations
  maxdepth = c(2, 3, 4),  # Maximum depth of trees
  nu = c(0.1, 0.5, 1.0)         # Learning rate (shrinkage parameter)
)

# Define base columns including tuning parameters
base_cols <- c("Accuracy", "Sensitivity", "Specificity", "Precision", 
               "NPV", "Threshold", "AUC", "iter", "maxdepth", "nu")

# Function to create results dataframe for all biomarkers model
create_all_biomarkers_results_df <- function(tuning_params) {
  total_rows <- nrow(tuning_params)
  
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
total_iterations <- nrow(tuning_params)
print(paste("Total iterations to process:", total_iterations))

# Check if confounders are in the dataset
valid_confounders <- confounders[confounders %in% colnames(data.dfPH)]
if (length(valid_confounders) < length(confounders)) {
  missing_confounders <- setdiff(confounders, valid_confounders)
  warning(glue("The following confounders were not found in the dataset and will be ignored: {paste(missing_confounders, collapse=', ')}"))
  confounders <- valid_confounders
}

# Set up parallel processing using parallelly 
numCores <- parallelly::availableCores() - 1
cl <- parallelly::makeClusterPSOCK(numCores)
registerDoParallel(cl)

# Export necessary functions and data to the cluster
clusterEvalQ(cl, {
  library(caret)
  library(pROC)
  library(dplyr)
  library(glue)
  library(ada)
})

# Export the custom functions and necessary variables
clusterExport(cl, c("customSummary", "data.dfPH", "all_biomarkers", "confounders", "control"))

# Prepare data with all biomarkers and confounders
subset_data <- data.dfPH[, c(all_biomarkers, confounders, "condL")]

# Create formula with all biomarkers and confounders
formula_str <- paste("condL ~", paste(c(all_biomarkers, confounders), collapse = " + "))
formula_obj <- as.formula(formula_str)

print(paste("Model formula:", formula_str))
print(paste("Total features in model:", length(all_biomarkers) + length(confounders)))

# Process all parameter combinations in parallel
all_results <- foreach(param_idx = 1:nrow(tuning_params), .combine = 'rbind',
                       .packages = c("caret", "pROC", "dplyr", "glue", "ada")) %dopar% {
                         
                         # Get the parameters from tuning_params
                         params <- tuning_params[param_idx, ]
                         
                         # Try to fit model with these parameters
                         tryCatch({
                           # Create a row for results
                           result_row <- data.frame(matrix(NA, nrow = 1, 
                                                           ncol = 1 + length(base_cols)),
                                                    stringsAsFactors = FALSE)
                           colnames(result_row) <- c("Model_Description", base_cols)
                           
                           # Store model description
                           result_row[1, "Model_Description"] <- paste0("All_", length(all_biomarkers), "_Biomarkers_Plus_", length(confounders), "_Confounders")
                           
                           # Store tuning parameters
                           result_row[1, "iter"] <- params$iter
                           result_row[1, "maxdepth"] <- params$maxdepth
                           result_row[1, "nu"] <- params$nu
                           
                           # Train AdaBoost model
                           set.seed(123)
                           ada_model <- train(
                             formula_obj,
                             data = subset_data,
                             method = "ada",
                             tuneGrid = params,
                             metric = "Accuracy",
                             trControl = control
                           )
                           
                           # Get the factor levels from the outcome
                           lev <- levels(subset_data$condL)
                           pred_subset <- ada_model$pred
                           
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
                           warning(glue("Error processing all biomarkers with iter={params$iter}, maxdepth={params$maxdepth}, nu={params$nu}: {e$message}"))
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
  # Split into chunks if necessary
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
  
  # Add a summary sheet with confounder list
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
    Parameter = c("iter", "maxdepth", "nu"),
    Values_Tested = c(
      paste(unique(tuning_params$iter), collapse = ", "),
      paste(unique(tuning_params$maxdepth), collapse = ", "),
      paste(unique(tuning_params$nu), collapse = ", ")
    ),
    stringsAsFactors = FALSE
  )
  export_list[["Parameter_Summary"]] <- param_summary
  
  # Add a summary of the best performing models
  if (nrow(all_results) > 0) {
    # Best by AUC
    best_auc_idx <- which.max(all_results$AUC)
    best_auc <- all_results[best_auc_idx, ]
    
    # Best by Accuracy
    best_accuracy_idx <- which.max(all_results$Accuracy)
    best_accuracy <- all_results[best_accuracy_idx, ]
    
    # Top 10 by AUC
    top_10_auc <- all_results[order(-all_results$AUC), ][1:min(10, nrow(all_results)), ]
    
    summary_results <- rbind(
      cbind(Rank = "Best_AUC", best_auc),
      cbind(Rank = "Best_Accuracy", best_accuracy)
    )
    
    export_list[["Best_Models_Summary"]] <- summary_results
    export_list[["Top_10_AUC"]] <- cbind(Rank = 1:nrow(top_10_auc), top_10_auc)
  }
  
  # Write all data frames to an Excel file
  write_xlsx(export_list, path = "Results_AdaBoost_BiomarkersSim.10.3.25.xlsx")
  
  print(paste("Analysis complete! Results saved for", length(all_biomarkers), "biomarkers and", length(confounders), "confounders"))
  print(paste("Total parameter combinations tested:", nrow(all_results)))
  
  # Print summary statistics
  if (nrow(all_results) > 0) {
    cat("\n=== SUMMARY STATISTICS ===\n")
    cat("Best AUC:", max(all_results$AUC, na.rm = TRUE), "\n")
    cat("Best Accuracy:", max(all_results$Accuracy, na.rm = TRUE), "\n")
    cat("Mean AUC:", round(mean(all_results$AUC, na.rm = TRUE), 4), "\n")
    cat("Mean Accuracy:", round(mean(all_results$Accuracy, na.rm = TRUE), 4), "\n")
    
    # Show best parameter combinations
    best_auc_params <- all_results[which.max(all_results$AUC), ]
    cat("\nBest AUC Parameters:\n")
    cat("  iter:", best_auc_params$iter, "\n")
    cat("  maxdepth:", best_auc_params$maxdepth, "\n")
    cat("  nu:", best_auc_params$nu, "\n")
  }
} else {
  print("No results to export - check for errors in model fitting")
}