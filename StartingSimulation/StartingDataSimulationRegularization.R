#Qi Yan
# 11.21.2025
# This code is for generating the simulating data


generate_rule_based_simulation <- function(n = 100, p_biomarkers = 20, seed = 42) {
  set.seed(seed)
  
  # Create a matrix of biomarker measurements drawn from a standard normal distribution.
  # Each row is one sample; each column is one biomarker.
  X <- matrix(rnorm(n * p_biomarkers), n, p_biomarkers)
  biomarkers <- paste0("bio", 1:p_biomarkers)
  colnames(X) <- biomarkers
  
  # Initialize a rule-based score for each sample.
  # This score will later determine the sample's class label.
  score <- rep(0, n)
  
  # Apply threshold-based rules:
  # For each rule, if the condition is satisfied, the score is adjusted.
  score <- score + 15 * (X[,1] > -0.1)
  score <- score - 12 * (X[,2] > 0.4)
  score <- score + 10 * (X[,3] < -0.2)
  score <- score - 8  * (X[,4] > 0.6)
  score <- score + 6  * (X[,5] < 0.1)
  
  # Amplify the separation between positive and negative scores.
  # Positive scores increase; negative scores decrease further.
  score <- score + ifelse(score > 0, 5, -5)
  
  # Convert scores into probabilities using a logistic function,
  # then push probabilities close to either 0 or 1 for near-deterministic labels.
  probs <- plogis(score)
  probs <- ifelse(probs > 0.5, 0.95, 0.05)
  
  # Generate class labels based on these probabilities.
  y <- ifelse(probs > 0.5, 1, 0)
  
  # If the outcome distribution becomes too imbalanced,
  # enforce a fixed number of positive samples chosen from high-scoring individuals.
  if(sum(y) < 20 || sum(y) > 80) {
    target_pos <- 40
    pos_candidates <- which(score > median(score))
    y <- rep(0, n)
    y[sample(pos_candidates, target_pos)] <- 1
  }
  
  # Build the output dataset:
  #  - biomarker values
  #  - binary label (cond)
  #  - human-readable label ("H" or "P")
  data.dfPH <- data.frame(X)
  data.dfPH$y <- y
  data.dfPH$cond <- data.dfPH$y
  names(data.dfPH)[names(data.dfPH) == "y"] <- "condL"
  data.dfPH$condL[data.dfPH$condL == 0] <- "H"
  data.dfPH$condL[data.dfPH$condL == 1] <- "P"
  data.dfPH$condL <- factor(data.dfPH$condL, levels = c("H", "P"))
  
  # Return the dataset and metadata describing which biomarkers drive the rules.
  return(list(
    data.dfPH = data.dfPH,
    biomarkers = biomarkers,
    confounders = NULL,
    important_biomarkers = paste0("bio", 1:5),
    noise_biomarkers = paste0("bio", 6:20)
  ))
}
# Run the simulation
cat("=== Option 1:===\n")
sim1 <- generate_rule_based_simulation(n = 100, p_biomarkers = 200, seed = 42)
data.dfPH <- sim1$data.dfPH
biomarkers <- sim1$biomarkers
confounders <- sim1$confounders
cat("Outcome distribution:", table(data.dfPH$condL), "\n")
# 