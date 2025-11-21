#3 simulation is the real one

#91% diagnostic accuracy
#7x7 cross fold, repeat


#91% diagnostic accuracy 
#7x5 cross fold
# 200 features, 100 ppl
# ==============================================================================
# REFINED 95% ACCURACY - MINIMAL CHANGES TO YOUR WORKING 90% VERSION
# ==============================================================================
generate_refined_separable_data <- function(n = 100, p_biomarkers = 20, seed = 42) {
  set.seed(seed)
  X <- matrix(rnorm(n * p_biomarkers), n, p_biomarkers)
  biomarkers <- paste0("bio", 1:p_biomarkers)
  colnames(X) <- biomarkers
  
  # START WITH YOUR WORKING 90% VERSION, MAKE MINIMAL TARGETED IMPROVEMENTS
  score <- rep(0, n)
  
  # Slightly increase the strongest effects from your working version
  score <- score + 18 * (X[,1] > 0)           # Was 15, now 18 (20% increase)
  score <- score - 15 * (X[,2] > 0.5)         # Was 12, now 15 (25% increase)  
  score <- score + 12 * (X[,3] < -0.3)        # Was 10, now 12 (20% increase)
  score <- score - 10 * (X[,4] > 0.7)         # Was 8, now 10 (25% increase)
  score <- score + 8 * (X[,5] < 0)            # Was 6, now 8 (33% increase)
  
  # Keep the same amplification as your working version
  score <- score + ifelse(score > 0, 5, -5)   # Same as before
  
  # Make probabilities slightly more extreme than your working version
  probs <- plogis(score)
  probs <- ifelse(probs > 0.5, 0.97, 0.03)    # Was 0.95/0.05, now 0.97/0.03
  
  # Same outcome generation as your working version
  y <- ifelse(probs > 0.5, 1, 0)
  
  # Same balance logic as your working version
  if(sum(y) < 20 || sum(y) > 80) {
    target_pos <- 40
    pos_candidates <- which(score > median(score))
    y <- rep(0, n)
    y[sample(pos_candidates, target_pos)] <- 1
  }
  
  # Same data frame creation as your working version
  data.dfPH <- data.frame(X)
  data.dfPH$y <- y
  data.dfPH$cond <- data.dfPH$y
  names(data.dfPH)[names(data.dfPH) == "y"] <- "condL"
  data.dfPH$condL[data.dfPH$condL == 0] <- "H"
  data.dfPH$condL[data.dfPH$condL == 1] <- "P"
  data.dfPH$condL <- factor(data.dfPH$condL, levels = c("H", "P"))
  
  return(list(
    data.dfPH = data.dfPH,
    biomarkers = biomarkers,
    confounders = NULL,
    important_biomarkers = paste0("bio", 1:5),
    noise_biomarkers = paste0("bio", 6:20)
  ))
}

# ==============================================================================
# ALTERNATIVE: JUST ADD ONE MORE STRONG RULE TO YOUR WORKING VERSION
# ==============================================================================
generate_one_more_rule_data <- function(n = 100, p_biomarkers = 20, seed = 42) {
  set.seed(seed)
  X <- matrix(rnorm(n * p_biomarkers), n, p_biomarkers)
  biomarkers <- paste0("bio", 1:p_biomarkers)
  colnames(X) <- biomarkers
  
  # EXACT COPY of your working 90% version
  score <- rep(0, n)
  score <- score + 15 * (X[,1] > 0)           # bio1 > 0 → P (huge effect)
  score <- score - 12 * (X[,2] > 0.5)         # bio2 > 0.5 → H
  score <- score + 10 * (X[,3] < -0.3)        # bio3 < -0.3 → P
  score <- score - 8 * (X[,4] > 0.7)          # bio4 > 0.7 → H
  score <- score + 6 * (X[,5] < 0)            # bio5 < 0 → P
  
  # ADD JUST ONE MORE STRONG RULE
  score <- score + 7 * (X[,6] > 0.4)          # bio6 > 0.4 → P (new rule)
  
  score <- score + ifelse(score > 0, 5, -5)   # Same amplification
  
  probs <- plogis(score)
  probs <- ifelse(probs > 0.5, 0.95, 0.05)    # Same as your working version
  y <- ifelse(probs > 0.5, 1, 0)
  
  # Same balance logic
  if(sum(y) < 20 || sum(y) > 80) {
    target_pos <- 40
    pos_candidates <- which(score > median(score))
    y <- rep(0, n)
    y[sample(pos_candidates, target_pos)] <- 1
  }
  
  # Same data frame creation
  data.dfPH <- data.frame(X)
  data.dfPH$y <- y
  data.dfPH$cond <- data.dfPH$y
  names(data.dfPH)[names(data.dfPH) == "y"] <- "condL"
  data.dfPH$condL[data.dfPH$condL == 0] <- "H"
  data.dfPH$condL[data.dfPH$condL == 1] <- "P"
  data.dfPH$condL <- factor(data.dfPH$condL, levels = c("H", "P"))
  
  return(list(
    data.dfPH = data.dfPH,
    biomarkers = biomarkers,
    confounders = NULL,
    important_biomarkers = paste0("bio", 1:6),  # Now includes bio6
    noise_biomarkers = paste0("bio", 7:20)
  ))
}

# ==============================================================================
# MOST CONSERVATIVE: JUST TIGHTEN THE THRESHOLDS
# ==============================================================================
generate_tighter_thresholds_data <- function(n = 100, p_biomarkers = 20, seed = 42) {
  set.seed(seed)
  X <- matrix(rnorm(n * p_biomarkers), n, p_biomarkers)
  biomarkers <- paste0("bio", 1:p_biomarkers)
  colnames(X) <- biomarkers
  
  # Keep exact same coefficients as your working version
  score <- rep(0, n)
  score <- score + 15 * (X[,1] > -0.1)        # Easier threshold: was >0, now >-0.1
  score <- score - 12 * (X[,2] > 0.4)         # Easier threshold: was >0.5, now >0.4
  score <- score + 10 * (X[,3] < -0.2)        # Easier threshold: was <-0.3, now <-0.2
  score <- score - 8 * (X[,4] > 0.6)          # Easier threshold: was >0.7, now >0.6
  score <- score + 6 * (X[,5] < 0.1)          # Easier threshold: was <0, now <0.1
  
  score <- score + ifelse(score > 0, 5, -5)   # Same amplification
  
  probs <- plogis(score)
  probs <- ifelse(probs > 0.5, 0.95, 0.05)
  y <- ifelse(probs > 0.5, 1, 0)
  
  if(sum(y) < 20 || sum(y) > 80) {
    target_pos <- 40
    pos_candidates <- which(score > median(score))
    y <- rep(0, n)
    y[sample(pos_candidates, target_pos)] <- 1
  }
  
  data.dfPH <- data.frame(X)
  data.dfPH$y <- y
  data.dfPH$cond <- data.dfPH$y
  names(data.dfPH)[names(data.dfPH) == "y"] <- "condL"
  data.dfPH$condL[data.dfPH$condL == 0] <- "H"
  data.dfPH$condL[data.dfPH$condL == 1] <- "P"
  data.dfPH$condL <- factor(data.dfPH$condL, levels = c("H", "P"))
  
  return(list(
    data.dfPH = data.dfPH,
    biomarkers = biomarkers,
    confounders = NULL,
    important_biomarkers = paste0("bio", 1:5),
    noise_biomarkers = paste0("bio", 6:20)
  ))
}

# ==============================================================================
# USAGE - TRY THESE IN ORDER OF CONSERVATISM
# ==============================================================================

# Option 1: Most conservative - just tighter thresholds
cat("=== Option 1: Tighter Thresholds (Most Conservative) ===\n")
sim1 <- generate_tighter_thresholds_data(n = 100, p_biomarkers = 200, seed = 42)
data.dfPH <- sim1$data.dfPH
biomarkers <- sim1$biomarkers
confounders <- sim1$confounders
cat("Outcome distribution:", table(data.dfPH$condL), "\n")
# 
# # Option 2: Add one more rule
# cat("\n=== Option 2: One Additional Rule ===\n")
# sim2 <- generate_one_more_rule_data(n = 100, p_biomarkers = 20, seed = 42)
# cat("Outcome distribution:", table(sim2$data.dfPH$condL), "\n")
# 
# # Option 3: Refined version with slight increases
# cat("\n=== Option 3: Refined Version (Slight Increases) ===\n")
# sim3 <- generate_refined_separable_data(n = 100, p_biomarkers = 20, seed = 42)
# cat("Outcome distribution:", table(sim3$data.dfPH$condL), "\n")
# 
# # Quick test function
# test_simple_rules <- function(data) {
#   rule1 <- mean((data$bio1 > -0.1) == (data$condL == "P"))
#   rule2 <- mean((data$bio2 <= 0.4) == (data$condL == "H")) 
#   cat("Rule 1 accuracy:", round(rule1, 3), "\n")
#   cat("Rule 2 accuracy:", round(rule2, 3), "\n")
# }

# cat("\n=== Quick Rule Tests ===\n")
# test_simple_rules(sim1$data.dfPH)