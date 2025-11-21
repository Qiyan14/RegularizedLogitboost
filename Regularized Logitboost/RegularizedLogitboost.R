# Qi Yan
# 11.21.2025
# Core code to run LogitBoost with regularization and adaptive learning rate

# Load necessary libraries
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

# [LogitBoostTesting function 
LogitBoostTesting <- function(
    xlearn, ylearn,
    nIter       = ncol(xlearn),
    regType     = "none",
    lambda      = 0.1,
    alpha       = NA,
    nu          = 0.1,      # base learning rate
    gamma       = 0.05,     # inverse‐decay rate
    adaptLR     = FALSE,    # if TRUE, use inverse linear decay
    subsample   = 0.7,
    mtry        = 0.2,
    early.stop  = 15,
    K           = 20
) {
  # ——— Prep outcome ———
  if (is.factor(ylearn)) {
    lablist        <- levels(ylearn)
    ylearn_numeric <- as.numeric(ylearn) - 1
  } else {
    ylearn_numeric <- ylearn
    lablist        <- sort(unique(ylearn_numeric))
  }
  if (length(unique(ylearn_numeric)) != 2 ||
      !all(unique(ylearn_numeric) %in% c(0,1))) {
    stop("Outcome must be binary (0/1)")
  }
  nLearn <- nrow(xlearn)
  nFeat  <- ncol(xlearn)
  ylearn <- ylearn_numeric
  
  # ——— Sparse support for very high‑dim ———
  if (!inherits(xlearn, "dgCMatrix") && nFeat > 500) {
    if (requireNamespace("Matrix", quietly=TRUE)) {
      xlearn <- as(as.matrix(xlearn), "dgCMatrix")
    }
  }
  
  # ——— Regularization setup ———
  if (regType == "none")      lambda <- 0
  if (regType == "l1")        alpha  <- 1
  if (regType == "l2")        alpha  <- 0
  if (regType == "elastic" && is.na(alpha)) alpha <- 0.5
  
  # ——— Initialize ———
  f       <- rep(0, nLearn)
  p       <- rep(0.5, nLearn)
  Stump   <- matrix(0, nIter, 4,
                    dimnames = list(NULL, c("feature","threshold","sign","coef")))
  best_oob <- -Inf
  no_imp   <- 0
  
  for (m in seq_len(nIter)) {
    # a) subsample rows & feats
    rows <- sample(nLearn, floor(subsample * nLearn), FALSE)
    feats <- sample(nFeat, max(1, floor(mtry * nFeat)), FALSE)
    
    # b) working weights & response
    w <- pmax(p * (1 - p), 1e-24)
    z <- (ylearn - p) / w
    w <- w / sum(w)
    S_total <- sum(w[rows] * z[rows])
    
    # c) top‑K feature preselection
    feature_scores <- sapply(feats, function(iFeat) {
      xi <- xlearn[rows, iFeat]
      if (var(xi) == 0) return(NA)
      abs(cor(as.numeric(xi), z[rows]))
    })
    K_eff    <- min(K, length(feats))
    top_feats <- feats[order(-feature_scores, na.last=NA)][seq_len(K_eff)]
    
    # d) stump search
    best_impr  <- -Inf
    bestFeat   <- NA
    bestThresh <- NA
    bestS      <- NA
    for (iFeat in top_feats) {
      xi <- xlearn[rows, iFeat]
      if (var(xi) == 0) next
      ord   <- order(xi)
      xi_o  <- xi[ord]
      wz_o  <- w[rows][ord] * z[rows][ord]
      csum  <- cumsum(wz_o)
      S_cand <- 2 * csum - S_total
      idx   <- which(diff(xi_o) != 0)
      if (!length(idx)) next
      S_u   <- S_cand[idx]
      t_u   <- xi_o[idx]
      lambda_eff <- lambda * (1 + m / nIter)
      impr <- (pmax(abs(S_u) - lambda_eff*alpha, 0)^2) /
        (1 + lambda_eff*(1 - alpha))
      if (!length(impr)) next
      i_max <- which.max(impr)
      if (!length(i_max)) next
      if (impr[i_max] > best_impr) {
        best_impr  <- impr[i_max]
        bestFeat   <- iFeat
        bestThresh <- t_u[i_max]
        bestS      <- S_u[i_max]
      }
    }
    
    # e) compute coef with inverse linear decay
    if (adaptLR) {
      nu_m <- nu / (1 + gamma * (m - 1))
    } else {
      nu_m <- nu
    }
    a <- if (best_impr <= 0 || is.na(bestFeat)) {
      0
    } else {
      raw <- sign(bestS) *
        (max(abs(bestS) - lambda_eff*alpha, 0) /
           (1 + lambda_eff*(1 - alpha)))
      nu_m * raw
    }
    
    # f) record stump
    if (is.na(bestFeat)) {
      Stump[m, ] <- c(1, 0, 0, 0)
    } else {
      Stump[m, ] <- c(bestFeat, bestThresh, sign(a), a)
    }
    
    # g) update predictions
    if (!is.na(bestFeat) && bestFeat <= ncol(xlearn)) {
      f <- f + a * (2 * (xlearn[,bestFeat] <= bestThresh) - 1)
    }
    p <- 1 / (1 + exp(-pmax(pmin(f, 500), -500)))
    
    # h) OOB early stopping
    if (!is.null(early.stop)) {
      oob <- setdiff(seq_len(nLearn), rows)
      if (length(oob)) {
        p_oob     <- 1 / (1 + exp(-pmax(pmin(f[oob], 500), -500)))
        pred_oob  <- ifelse(p_oob > 0.5, 1, 0)
        acc       <- mean(pred_oob == ylearn[oob])
        if (acc > best_oob + 1e-6) {
          best_oob <- acc; no_imp <- 0
        } else {
          no_imp <- no_imp + 1
        }
        if (no_imp >= early.stop) {
          Stump <- Stump[1:m, , drop=FALSE]
          break
        }
      }
    }
  }
  
  # ——— Return model object ———
  object <- list(
    Stump     = Stump,
    lablist   = lablist,
    regType   = regType,
    lambda    = lambda,
    alpha     = alpha,
    nu        = nu,
    gamma     = gamma,
    adaptLR   = adaptLR,
    subsample = subsample,
    mtry      = mtry
  )
  class(object) <- "LogitBoostTesting"
  object
}

# [predict.LogitBoostTesting function 
predict.LogitBoostTesting <- function(object, xtest,
                                      type = c("class","raw"),
                                      nIter = NA, ...) {
  type    <- match.arg(type)
  Stump   <- object$Stump
  lablist <- object$lablist
  
  if (is.na(nIter)) nIter <- nrow(Stump) else
    nIter <- min(nIter, nrow(Stump))
  
  nTest  <- nrow(xtest)
  # Binary case
  f <- numeric(nTest)
  for (i in seq_len(nIter)) {
    feat <- Stump[i, "feature"]
    thr  <- Stump[i, "threshold"]
    coef <- Stump[i, "coef"]
    if (!is.na(feat) && feat > 0 && feat <= ncol(xtest)) {
      f <- f + coef * (2 * (xtest[, feat] <= thr) - 1)
    }
  }
  f   <- pmax(pmin(f, 500), -500)
  pr  <- 1 / (1 + exp(-f))
  Prob <- cbind(1 - pr, pr)
  colnames(Prob) <- lablist
  
  if (type == "raw") return(Prob)
  
  # class labels
  ord    <- ifelse(Prob[,2] > Prob[,1], 2, 1)
  labels <- lablist[ord]
  # ties → NA
  ties <- Prob[,1] == Prob[,2]
  labels[ties] <- NA
  labels
}

# [LogitBoostRegularizedModel caret definition 
LogitBoostRegularizedModel <- list(
  label    = "Boosted Logistic Regression with Adaptive LR",
  library  = NULL,
  loop     = function(grid) {
    loop      <- grid[which.max(grid$nIter), , drop = FALSE]
    submodels <- grid[-which.max(grid$nIter), , drop = FALSE]
    list(loop = loop, submodels = list(submodels))
  },
  type     = "Classification",
  parameters = data.frame(
    parameter = c('nIter','regType','lambda','alpha',
                  'nu','gamma','adaptLR','subsample','mtry'),
    class     = c('numeric','character','numeric','numeric',
                  'numeric','numeric','logical','numeric','numeric'),
    label     = c('# Boosting Iterations',
                  'Regularization Type ("none","l1","l2","elastic")',
                  'Penalty Strength',
                  'Mixing Parameter (for elastic)',
                  'Base Learning Rate (nu)',
                  'Inverse‐Decay Rate (γ)',
                  'Use Inverse‐Linear Decay (TRUE/FALSE)',
                  'Row Subsample Fraction',
                  'Feature Subsample Fraction'),
    stringsAsFactors = FALSE
  ),
  grid = function(x, y, len = NULL, search = "grid") {
    if (search == "grid") {
      data.frame(
        nIter     = 1 + ((1:len) * 10),
        regType   = rep("none", len),
        lambda    = rep(0.1, len),
        alpha     = rep(NA,  len),
        nu        = rep(0.1, len),
        gamma     = rep(0.05, len),
        adaptLR   = rep(TRUE, len),
        subsample = rep(0.7, len),
        mtry      = rep(0.5, len),
        stringsAsFactors = FALSE
      )
    } else {
      data.frame(
        nIter     = unique(sample(1:100, size = len, replace = TRUE)),
        regType   = rep("none",        len),
        lambda    = runif(len, 0, 1),
        alpha     = runif(len, 0, 1),
        nu        = runif(len, 0, 1),
        gamma     = runif(len, 0, 1),
        adaptLR   = sample(c(TRUE,FALSE), len, replace = TRUE),
        subsample = runif(len, 0.5, 1),
        mtry      = runif(len, 0.1, 1),
        stringsAsFactors = FALSE
      )
    }
  },
  fit = function(x, y, wts, param, lev, last, classProbs, ...) {
    if (!is.factor(y)) y <- factor(y, levels = lev)
    model <- LogitBoostTesting(
      as.matrix(x), y,
      nIter     = param$nIter,
      regType   = param$regType,
      lambda    = param$lambda,
      alpha     = param$alpha,
      nu        = param$nu,
      gamma     = param$gamma,
      adaptLR   = param$adaptLR,
      subsample = param$subsample,
      mtry      = param$mtry
    )
    model$obsLevels <- lev
    model
  },
  predict = function(modelFit, newdata, submodels = NULL) {
    out <- predict.LogitBoostTesting(modelFit, newdata, "class")
    if (!is.null(submodels)) {
      tmp <- out
      out <- vector("list", nrow(submodels) + 1)
      out[[1]] <- tmp
      for (i in seq_len(nrow(submodels))) {
        out[[i+1]] <- predict.LogitBoostTesting(
          modelFit, newdata,
          nIter = submodels$nIter[i]
        )
      }
    }
    out
  },
  prob = function(modelFit, newdata, submodels = NULL) {
    probs <- predict.LogitBoostTesting(modelFit, newdata, "raw")
    probs <- t(apply(probs, 1, function(z) z / sum(z)))
    if (!is.null(submodels)) {
      out <- list(probs)
      for (i in seq_len(nrow(submodels))) {
        tmp <- predict.LogitBoostTesting(
          modelFit, newdata,
          type  = "raw",
          nIter = submodels$nIter[i]
        )
        out[[i+1]] <- as.data.frame(t(apply(tmp, 1, function(z) z / sum(z))),
                                    stringsAsFactors = TRUE)
      }
      return(out)
    }
    probs
  },
  predictors = function(x, ...) {
    if (!is.null(x$xNames)) unique(x$xNames[x$Stump[, "feature"]])
    else NA
  },
  levels = function(x) x$obsLevels,
  tags   = c("Ensemble Model","Boosting","Adaptive LR","Feature Selection"),
  sort   = function(x) x[order(x[,1]), ]
)

# customSummary function for caret
customSummary <- function(data, lev = NULL, model = NULL) {
  if (length(lev) > 2) {
    stop(paste("Your outcome has", length(lev),
               "levels. TwoClassSummary isn't appropriate."))
  }
  if (!requireNamespace("pROC", quietly = TRUE)) {
    stop("Package pROC needed for this function to work. Please install it.")
  }
  if (!all(levels(data[, "pred"]) == lev)) {
    stop("levels of observed and predicted data do not match")
  }
  
  # Confusion matrix counts
  cm <- table(data[, "pred"], data[, "obs"])
  tp <- cm[lev[2], lev[2]]
  tn <- cm[lev[1], lev[1]]
  fp <- cm[lev[2], lev[1]]
  fn <- cm[lev[1], lev[2]]
  
  sensitivity <- if ((tp + fn) > 0) tp / (tp + fn) else 0
  specificity <- if ((tn + fp) > 0) tn / (tn + fp) else 0
  accuracy    <- (tp + tn) / (tp + tn + fp + fn)
  
  c(Sens = sensitivity, Spec = specificity, Accuracy = accuracy)
}


# Main analysis code

# Data.dfPH, biomarkers, confounders from the StartingData Code
data.dfPH$condL <- factor(data.dfPH$condL)
levels(data.dfPH$condL) <- make.names(levels(data.dfPH$condL))
all_biomarkers <- biomarkers
print(paste("Using all", length(all_biomarkers), "biomarkers in the model"))

# Convert factor predictors to numeric/dummies
factor_vars <- sapply(data.dfPH[, c(all_biomarkers, confounders)], is.factor)
if (any(factor_vars)) {
  vars <- names(factor_vars)[factor_vars]
  warning(paste("Converting factor variables:", paste(vars, collapse = ", ")))
  for (var in vars) {
    if (nlevels(data.dfPH[[var]]) == 2) {
      data.dfPH[[var]] <- as.numeric(data.dfPH[[var]]) - 1
    } else {
      mm <- model.matrix(~ . - 1, data = data.dfPH[var])
      colnames(mm) <- paste0(var, "_", colnames(mm))
      data.dfPH <- cbind(data.dfPH, mm)
      data.dfPH[[var]] <- NULL
      if (var %in% all_biomarkers) {
        all_biomarkers <- c(setdiff(all_biomarkers, var), colnames(mm))
      }
      if (var %in% confounders) {
        confounders <- c(setdiff(confounders, var), colnames(mm))
      }
    }
  }
}

# Control for repeated 7-fold cross-validation
control <- trainControl(
  method          = "repeatedcv",
  number          = 7,
  repeats         = 5,
  summaryFunction = customSummary,
  classProbs      = TRUE,
  allowParallel   = TRUE,
  savePredictions = "final"
)

# UPDATED: Expanded tuning values to include adaptive learning rate parameters
iterations_to_test  <- c(10, 25, 50, 100,200,500)
lambda_values       <- c(0.001, 0.01, 0.1)
nu_values           <- c(0.01, 0.05, 0.1, 0.5)  # Base learning rates
gamma_values        <- c(0.01, 0.05, 0.1, 0.5)  # Inverse-decay rates
adaptLR_values      <- c(TRUE, FALSE)            # Test both adaptive and fixed
subsample_values    <- c(0.3, 0.5, 0.7, 1.0)
mtry_values         <- c(0.3, 0.5, 0.7)



# Result columns 
base_cols <- c(
  "Accuracy", "Sensitivity", "Specificity", "Precision", "NPV", "AUC",
  "Feature_Coefficients", "Num_Nonzero_Coefficients", "Threshold",
  "nIter", "regType", "lambda", "alpha", "nu", "gamma", "adaptLR",
  "subsample", "mtry"
)

# Define regType/alpha combos
parameters <- list(
  l1            = list(regType="l1",      alpha=1),
  l2            = list(regType="l2",      alpha=0),
  elastic_0.20  = list(regType="elastic", alpha=0.2),
  elastic_0.50  = list(regType="elastic", alpha=0.5),
  elastic_0.80  = list(regType="elastic", alpha=0.8)
)
param_df <- do.call(rbind, lapply(names(parameters), function(key) {
  data.frame(
    param_key = key,
    regType   = parameters[[key]]$regType,
    alpha     = parameters[[key]]$alpha,
    stringsAsFactors = FALSE
  )
}))

# Calculate total iterations including new parameters
total_iterations <- length(iterations_to_test) *
  length(lambda_values) *
  nrow(param_df) *
  length(nu_values) *
  length(gamma_values) *
  length(adaptLR_values) *
  length(subsample_values) *
  length(mtry_values)
print(paste("Total iterations to process:", total_iterations))

# Subset data and feature names
subset_data <- data.dfPH[, c(all_biomarkers, confounders, "condL")]
feature_names <- colnames(subset_data[, c(all_biomarkers, confounders)])

# Data quality checks 
print("=== DATA QUALITY CHECKS ===")
print(paste("Dimensions:", nrow(subset_data), "x", ncol(subset_data)))
print(paste("Outcome levels:", paste(levels(subset_data$condL), collapse=", ")))
print(paste("Outcome counts:", paste(table(subset_data$condL), collapse=" | ")))

# Remove constant predictors
constant_vars <- sapply(subset_data[, c(all_biomarkers, confounders)],
                        function(x) var(x, na.rm=TRUE) == 0)
if (any(constant_vars)) {
  const_names <- names(constant_vars)[constant_vars]
  warning(paste("Removing constant vars:", paste(const_names, collapse=", ")))
  all_biomarkers <- setdiff(all_biomarkers, const_names)
  confounders    <- setdiff(confounders,    const_names)
  subset_data    <- subset_data[, !names(subset_data) %in% const_names]
  feature_names  <- colnames(subset_data[, c(all_biomarkers, confounders)])
}

# Handle missing data
if (any(is.na(subset_data))) {
  subset_data <- na.omit(subset_data)
  if (nrow(subset_data) < 0.8 * nrow(data.dfPH)) {
    warning(">20% data lost to NA; consider imputation.")
  }
}

# Formula
formula_str <- paste("condL ~", paste(c(all_biomarkers, confounders), collapse=" + "))
formula_obj <- as.formula(formula_str)
print(paste("Model formula:", formula_str))

# Helper function for feature importance 
get_feature_importance <- function(stump_matrix, feature_names) {
  feats   <- stump_matrix[, "feature"]
  coefs   <- as.numeric(stump_matrix[, "coef"])
  out     <- setNames(rep(0, length(feature_names)), feature_names)
  for (j in seq_along(feats)) {
    if (!is.na(feats[j]) && feats[j] > 0 && feats[j] <= length(feature_names)) {
      nm <- feature_names[feats[j]]
      out[nm] <- out[nm] + abs(coefs[j])
    }
  }
  # Only keep non-zero coefficients, and sort descending
  out <- out[out != 0]
  if (length(out) > 0) {
    out <- sort(out, decreasing = TRUE)
  }
  return(out)
}

# Parallel backend
numCores <- parallelly::availableCores() - 1
cl <- parallelly::makeClusterPSOCK(numCores)
registerDoParallel(cl)
clusterEvalQ(cl, {
  library(caret); library(pROC); library(dplyr); library(glue)
})
clusterExport(cl, c("LogitBoostTesting","predict.LogitBoostTesting",
                    "LogitBoostRegularizedModel","customSummary",
                    "subset_data","all_biomarkers","confounders",
                    "control","param_df","formula_obj","feature_names",
                    "get_feature_importance"))

print(paste("Using", numCores, "cores for parallel processing"))

# Grid search with adaptive learning rate parameters
all_results <- foreach(param_idx    = 1:nrow(param_df), .combine = 'rbind') %:%
  foreach(nu_val         = nu_values,        .combine = 'rbind') %:%
  foreach(gamma_val      = gamma_values,     .combine = 'rbind') %:%
  foreach(adaptLR_val    = adaptLR_values,   .combine = 'rbind') %:%
  foreach(sub_val        = subsample_values, .combine = 'rbind') %:%
  foreach(mtry_val       = mtry_values,      .combine = 'rbind') %:%
  foreach(lambda_val     = lambda_values,    .combine = 'rbind') %:%
  foreach(iter_val       = iterations_to_test,.combine = 'rbind',
          .packages      = c("caret","pROC","dplyr","glue")) %dopar% {
            
            regType_val <- param_df$regType[param_idx]
            alpha_val   <- param_df$alpha[param_idx]
            
            tryCatch({
              # Prepare result row with updated columns
              result_row <- data.frame(matrix(NA, nrow=1, ncol=1+length(base_cols)),
                                       stringsAsFactors=FALSE)
              colnames(result_row) <- c("Model_Description", base_cols)
              
              result_row$Model_Description <- paste0(
                "All_", length(all_biomarkers), "_Biomarkers_Plus_",
                length(confounders), "_Confounders"
              )
              result_row$nIter     <- iter_val
              result_row$regType   <- regType_val
              result_row$lambda    <- lambda_val
              result_row$alpha     <- alpha_val
              result_row$nu        <- nu_val
              result_row$gamma     <- gamma_val    
              result_row$adaptLR   <- adaptLR_val  
              result_row$subsample <- sub_val
              result_row$mtry      <- mtry_val
              
              # Tuning parameters 
              tuning_param <- data.frame(
                nIter     = iter_val,
                regType   = regType_val,
                lambda    = lambda_val,
                alpha     = alpha_val,
                nu        = nu_val,
                gamma     = gamma_val,    
                adaptLR   = adaptLR_val,  
                subsample = sub_val,
                mtry      = mtry_val,
                stringsAsFactors = FALSE
              )
              
              # Running caret
              set.seed(123)
              logitboost_model <- train(
                formula_obj,
                data      = subset_data,
                method    = LogitBoostRegularizedModel,
                tuneGrid  = tuning_param,
                metric    = "Accuracy",
                trControl = control
              )
              
              lev         <- levels(subset_data$condL)
              pred_subset <- logitboost_model$pred
              
              if (length(unique(pred_subset$obs)) < 2 ||
                  !all(lev %in% colnames(pred_subset))) {
                metrics <- rep(NA, 7)
              } else {
                roc_obj <- suppressMessages(
                  roc(pred_subset$obs, pred_subset[[lev[2]]], quiet = TRUE)
                )
                M <- coords(
                  roc_obj, "best", best.method = "youden",
                  ret = c("threshold","accuracy","sensitivity",
                          "specificity","precision","npv")
                )
                metrics <- c(M["accuracy"], M["sensitivity"], M["specificity"],
                             M["precision"], M["npv"], M["threshold"],
                             as.numeric(roc_obj$auc))
              }
              
              result_row[, c("Accuracy","Sensitivity","Specificity",
                             "Precision","NPV","Threshold","AUC")] <- metrics
              
              # Extract nonzero feature coefficients and their count
              fitted_model <- logitboost_model$finalModel
              feature_importance <- get_feature_importance(
                fitted_model$Stump,
                feature_names
              )
              result_row$Feature_Coefficients <- if (length(feature_importance) > 0) {
                paste(
                  paste(names(feature_importance), round(feature_importance, 4), sep = ":"),
                  collapse = "; "
                )
              } else {
                ""
              }
              result_row$Num_Nonzero_Coefficients <- length(feature_importance)
              
              result_row
              
            }, error = function(e) {
              warning(glue(
                "Error regType={regType_val}, lambda={lambda_val}, ",
                "alpha={alpha_val}, nu={nu_val}, gamma={gamma_val}, ",
                "adaptLR={adaptLR_val}, subsample={sub_val}, ",
                "mtry={mtry_val}, iter={iter_val}: {e$message}"
              ))
              NULL
            })
          }

# Stop cluster
tryCatch({ parallel::stopCluster(cl) }, error = function(e) {
  message("Error stopping cluster: ", e$message)
})
registerDoSEQ()

# Post‐process results
if (!is.null(all_results)) {
  all_results <- all_results %>% mutate_if(is.numeric, ~ round(., 5))
  all_results$regType <- as.character(all_results$regType)
  all_results$adaptLR <- as.logical(all_results$adaptLR)  # Ensure proper logical type
}

# Prepare export
export_list <- list()
if (!is.null(all_results) && nrow(all_results) > 0) {
  if (nrow(all_results) > 1e6) {
    chunks <- split(all_results, ceiling(seq_len(nrow(all_results))/1e6))
    for (j in seq_along(chunks)) {
      export_list[[paste0("Results_All_Biomarkers_Part", j)]] <- chunks[[j]]
    }
  } else {
    export_list[["Results_All_Biomarkers"]] <- all_results
  }
  
  biomarker_summary <- data.frame(
    Biomarker_Number = seq_along(all_biomarkers),
    Biomarker_Name   = all_biomarkers,
    stringsAsFactors = FALSE
  )
  export_list[["Biomarker_List"]] <- biomarker_summary
  
  if (length(confounders) > 0) {
    confounder_summary <- data.frame(
      Confounder_Number = seq_along(confounders),
      Confounder_Name   = confounders,
      stringsAsFactors   = FALSE
    )
    export_list[["Confounder_List"]] <- confounder_summary
  }
  
  # Best models analysis 
  best_auc      <- all_results[which.max(all_results$AUC), ]
  best_accuracy <- all_results[which.max(all_results$Accuracy), ]
  
  # Additional analysis: Compare adaptive vs non-adaptive performance
  adaptive_results <- all_results[all_results$adaptLR == TRUE, ]
  fixed_results    <- all_results[all_results$adaptLR == FALSE, ]
  
  if (nrow(adaptive_results) > 0 && nrow(fixed_results) > 0) {
    adaptive_summary <- data.frame(
      Learning_Rate_Type = "Adaptive",
      Mean_AUC = round(mean(adaptive_results$AUC, na.rm = TRUE), 4),
      Best_AUC = round(max(adaptive_results$AUC, na.rm = TRUE), 4),
      Mean_Accuracy = round(mean(adaptive_results$Accuracy, na.rm = TRUE), 4),
      Best_Accuracy = round(max(adaptive_results$Accuracy, na.rm = TRUE), 4),
      Count = nrow(adaptive_results),
      stringsAsFactors = FALSE
    )
    
    fixed_summary <- data.frame(
      Learning_Rate_Type = "Fixed",
      Mean_AUC = round(mean(fixed_results$AUC, na.rm = TRUE), 4),
      Best_AUC = round(max(fixed_results$AUC, na.rm = TRUE), 4),
      Mean_Accuracy = round(mean(fixed_results$Accuracy, na.rm = TRUE), 4),
      Best_Accuracy = round(max(fixed_results$Accuracy, na.rm = TRUE), 4),
      Count = nrow(fixed_results),
      stringsAsFactors = FALSE
    )
    
    lr_comparison <- rbind(adaptive_summary, fixed_summary)
    export_list[["Learning_Rate_Comparison"]] <- lr_comparison
  }
  
  # Parameter sensitivity analysis 
  if (nrow(adaptive_results) > 0) {
    # Best gamma values analysis
    gamma_performance <- adaptive_results %>%
      group_by(gamma) %>%
      summarise(
        Mean_AUC = round(mean(AUC, na.rm = TRUE), 4),
        Best_AUC = round(max(AUC, na.rm = TRUE), 4),
        Mean_Accuracy = round(mean(Accuracy, na.rm = TRUE), 4),
        Count = n(),
        .groups = 'drop'
      ) %>%
      arrange(desc(Mean_AUC))
    
    export_list[["Gamma_Performance_Analysis"]] <- gamma_performance
    
    # Best nu values analysis
    nu_performance <- adaptive_results %>%
      group_by(nu) %>%
      summarise(
        Mean_AUC = round(mean(AUC, na.rm = TRUE), 4),
        Best_AUC = round(max(AUC, na.rm = TRUE), 4),
        Mean_Accuracy = round(mean(Accuracy, na.rm = TRUE), 4),
        Count = n(),
        .groups = 'drop'
      ) %>%
      arrange(desc(Mean_AUC))
    
    export_list[["Nu_Performance_Analysis"]] <- nu_performance
  }
  
  # Overall best models summary
  summary_results <- rbind(
    cbind(Metric="Best_AUC", best_auc),
    cbind(Metric="Best_Accuracy", best_accuracy)
  )
  export_list[["Best_Models_Summary"]] <- summary_results
  
  # Save filename
  write_xlsx(export_list,
             path="Results_LogitBoostRegularizedRepeatingestingCompSim.xlsx")
  
  message(glue(
    "Analysis complete! Results saved for {length(all_biomarkers)} biomarkers ",
    "and {length(confounders)} confounders. Tested {nrow(all_results)} parameter sets ",
    "including adaptive learning rate variations."
  ))
  
  cat("\n=== SUMMARY STATISTICS ===\n")
  cat("Best AUC:     ", max(all_results$AUC,      na.rm = TRUE), "\n")
  cat("Best Accuracy:", max(all_results$Accuracy, na.rm = TRUE), "\n")
  cat("Mean AUC:     ", round(mean(all_results$AUC,      na.rm = TRUE), 4), "\n")
  cat("Mean Accuracy:", round(mean(all_results$Accuracy, na.rm = TRUE), 4), "\n")
  
  # Print adaptive vs fixed learning rate comparison
  if (exists("lr_comparison")) {
    cat("\n=== ADAPTIVE vs FIXED LEARNING RATE COMPARISON ===\n")
    print(lr_comparison)
    
    if (nrow(adaptive_results) > 0 && nrow(fixed_results) > 0) {
      # Statistical test for difference
      auc_diff <- mean(adaptive_results$AUC, na.rm = TRUE) - 
        mean(fixed_results$AUC, na.rm = TRUE)
      cat("\nMean AUC Difference (Adaptive - Fixed):", round(auc_diff, 4), "\n")
      
      acc_diff <- mean(adaptive_results$Accuracy, na.rm = TRUE) - 
        mean(fixed_results$Accuracy, na.rm = TRUE)
      cat("Mean Accuracy Difference (Adaptive - Fixed):", round(acc_diff, 4), "\n")
    }
  }
  
  # Print top gamma values if available
  if (exists("gamma_performance")) {
    cat("\n=== TOP GAMMA VALUES FOR ADAPTIVE LEARNING RATE ===\n")
    print(head(gamma_performance, 3))
  }
  
  # Print top nu values if available
  if (exists("nu_performance")) {
    cat("\n=== TOP NU VALUES FOR ADAPTIVE LEARNING RATE ===\n")
    print(head(nu_performance, 3))
  }
  
} else {
  warning("No results to export – please check for errors in model fitting.")
}

# Function to analyze learning rate decay patterns
analyze_learning_rate_decay <- function(nu, gamma, max_iter = 100) {
  iterations <- 1:max_iter
  adaptive_rates <- nu / (1 + gamma * (iterations - 1))
  fixed_rates <- rep(nu, max_iter)
  
  decay_data <- data.frame(
    Iteration = iterations,
    Adaptive_LR = adaptive_rates,
    Fixed_LR = fixed_rates,
    Decay_Ratio = adaptive_rates / nu
  )
  
  return(decay_data)
}
# 
# # Visualize learning rate decay for different gamma values
# cat("\n=== LEARNING RATE DECAY EXAMPLES ===\n")
# for (gamma_val in c(0.01, 0.05, 0.1, 0.2)) {
#   decay_example <- analyze_learning_rate_decay(nu = 0.1, gamma = gamma_val, max_iter = 50)
#   final_lr <- tail(decay_example$Adaptive_LR, 1)
#   cat(paste0("Gamma = ", gamma_val, 
#              ", Final LR at iter 50: ", round(final_lr, 4),
#              " (", round(100 * final_lr / 0.1, 1), "% of initial)\n"))
# }