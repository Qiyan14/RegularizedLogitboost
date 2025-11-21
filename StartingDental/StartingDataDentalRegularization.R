###########
# Qi Yan
# 11.20.2025
# This file is used to prepare the data set for the dental data

# Load Library
library(VIM)
library(dbplyr)
library(dplyr)
library(readxl)

# Read in data

data.df1 <- read_excel("C:/Users/qiyan/OneDrive - University of Kentucky/MergedForZhang12JUL2025_from-29NOV2023.xlsx")
#data.df1 <- read_excel("C:/Users/qya231/OneDrive - University of Kentucky/MergedForZhang12JUL2025_from-29NOV2023.xlsx")


# Remove unneeded columns
data.df1 <- data.df1[, !(names(data.df1) %in% c("ID...28", "Group...29", "ID...30"))]

names(data.df1) <- gsub("-V1", "", names(data.df1))

names(data.df1) <- gsub("\\.\\.\\..*$", "", names(data.df1))

data.df1 <- data.df1[!(data.df1$ID %in% c("DM021")), ]

# Create a copy for processing
data.df2 <- data.df1

# Recode Visit.1 categories
names(data.df2) <- gsub("Visit 1", "Visit.1", names(data.df2))

data.df2$Visit.1[data.df2$Visit.1 %in% c("GG", "LG")] <- "G"
data.df2$Visit.1[data.df2$Visit.1 %in% c("GP", "LP")] <- "P"
data.df2$Visit.1[data.df2$Visit.1 %in% c("GG/LP", "LP")] <- NA
data.df2$Visit.1[!data.df2$Visit.1 %in% c("G", "P", "H")] <- NA

# Remove specific rows and convert character columns to numeric
data.df2 <- data.df2[-(86:92), ]
data.df2 <- data.df2[-(86:92), ]
char_columns <- sapply(data.df2[, -2, drop = FALSE], is.character)
for (col in names(char_columns[char_columns][-1])) {
  suppressWarnings(data.df2[[col]] <- as.numeric(data.df2[[col]]))
}

# Factor conversions
data.df2$M1F2 <- as.factor(data.df2$M1F2)
data.df2$RaceW1H2B3 <- as.factor(data.df2$RaceW1H2B3)
data.df2$Visit.1 <- as.factor(data.df2$Visit.1)

# Filter to relevant groups
data.df2 <- subset(data.df2, Visit.1 %in% c("G", "H", "P"))

# Create a copy for imputation
data.df5 <- data.df2

# kNN Imputation
names(data.df5) <- gsub("-", ".", names(data.df5))
cat('kNN_vars <- c("', paste(names(data.df2), collapse = '", "'), '")\n', sep = "")
# Define variables for kNN imputation
kNN_vars <- c("Age", "M1F2", "RaceW1H2B3", "A1C", "BMI", "Total_sites",  "MMP.9", "S100A8", "BAFF", "IFN.alpha", "TIMP.1", "IL1b", "IL6", "MIP1a", "MMP8", "Adiponectin", "Resistin", "PGE2", "Otu001", "Otu002", "Otu003", "Otu004", "Otu005", "Otu006", "Otu007", "Otu008", "Otu009", "Otu010", "Otu011", "Otu012", "Otu013", "Otu014", "Otu015", "Otu016", "Otu017", "Otu018", "Otu019", "Otu020", "Otu021", "Otu022", "Otu023", "Otu024", "Otu025", "Otu026", "Otu027", "Otu028", "Otu029", "Otu031", "Otu033", "Otu034", "Otu035", "Otu038", "Otu039", "Otu042", "Otu044", "Otu045", "Otu046", "Otu047", "Otu048", "Otu049", "Otu050", "Otu051", "Otu052", "Otu053", "Otu054", "Otu055", "Otu056", "Otu057", "Otu058", "Otu059", "Otu062", "Otu063", "Otu064", "Otu065", "Otu066", "Otu067", "Otu068", "Otu069", "Otu070", "Otu071", "Otu073", "Otu074", "Otu078", "Otu080", "Otu083", "Otu085", "Otu087", "Otu091", "Otu092", "Otu094", "Otu095", "Otu097", "Otu098", "Otu099", "Otu105", "Otu106", "Otu110", "Otu111", "Otu112", "Otu114", "Otu115", "Otu119", "Otu123", "Otu125", "Otu129", "Otu130", "Otu135", "Otu136", "Otu137", "Otu146", "Otu147", "Otu150", "Otu152", "Otu156", "Otu158", "Otu163", "Otu164", "Otu165", "Otu171", "Otu173", "Otu185", "Otu252")
# Impute and remove _imp columns
data.df5 <- kNN(data.df5, variable = kNN_vars)

# Group labels
data.df5$condL <- data.df5$Visit.1
data.df5$cond <- ifelse(data.df5$Visit.1 == "P", 1,
                        ifelse(data.df5$Visit.1 == "G", 1, 0))
data.df5$cond <- as.factor(data.df5$cond)
data.df5$cond <- relevel(data.df5$cond, ref = "0")
data.df5$condL <- as.factor(data.df5$condL)

data.df5 <- data.df5[, !grepl("_imp", names(data.df5))]


# Add small constant to all numeric values to account for measurement limitations
data.df5 <- data.df5 %>%
  mutate(across(where(is.numeric), ~ . + (1e-06 / 2)))

# Variables to remove
remove_vars <- c("Age", "A1C", "M1F2", "BMI", "RaceW1H2B3", "Total_sites")

# Updated kNN_vars without the removed variables
kNN_vars <- setdiff(kNN_vars, remove_vars)

data.df5 <- data.df5[, !(names(data.df5) %in% c("Visit", "Group...29", "ID...30"))]


# Log2 transformation for selected columns
data.df5[, kNN_vars] <- lapply(data.df5[, kNN_vars], log2)

# Standardize using stats from group H only
stats_df <- data.frame(Column = character(), Mean = numeric(), SD = numeric(), stringsAsFactors = FALSE)

for (col in kNN_vars) {
  mean_value <- mean(data.df5[data.df5$condL == "G", col], na.rm = TRUE)
  sd_value <- sd(data.df5[data.df5$condL == "G", col], na.rm = TRUE)
  
  data.df5[[col]] <- (data.df5[[col]] - mean_value) / sd_value
  
  stats_df <- rbind(stats_df, data.frame(Column = col, Mean = mean_value, SD = sd_value))
}

# Final column name fixes

colnames(data.df5)[colnames(data.df5) == "Otu156.Ff.59.3."] <- "Otu156"


# FINAL STEP: split into GH, PH, GP
data.dfPH <- subset(data.df5, Visit.1 != "G")
data.dfGP <- subset(data.df5, Visit.1 != "H")
data.dfGH <- subset(data.df5, Visit.1 != "P")

# Check
summary(data.dfPH)
summary(data.dfGH$condL)
summary(data.dfGP$condL)
# Biomarkers to be used in the model
biomarkers <- c("MMP.9", "S100A8", "BAFF", "IFN.alpha", "TIMP.1", "IL1b", "IL6", "MIP1a", "MMP8", "Adiponectin", "Resistin", "PGE2", "Otu001", "Otu002", "Otu003", "Otu004", "Otu005", "Otu006", "Otu007", "Otu008", "Otu009", "Otu010", "Otu011", "Otu012", "Otu013", "Otu014", "Otu015", "Otu016", "Otu017", "Otu018", "Otu019", "Otu020", "Otu021", "Otu022", "Otu023", "Otu024", "Otu025", "Otu026", "Otu027", "Otu028", "Otu029", "Otu031", "Otu033", "Otu034", "Otu035", "Otu038", "Otu039", "Otu042", "Otu044", "Otu045", "Otu046", "Otu047", "Otu048", "Otu049", "Otu050", "Otu051", "Otu052", "Otu053", "Otu054", "Otu055", "Otu056", "Otu057", "Otu058", "Otu059", "Otu062", "Otu063", "Otu064", "Otu065", "Otu066", "Otu067", "Otu068", "Otu069", "Otu070", "Otu071", "Otu073", "Otu074", "Otu078", "Otu080", "Otu083", "Otu085", "Otu087", "Otu091", "Otu092", "Otu094", "Otu095", "Otu097", "Otu098", "Otu099", "Otu105", "Otu106", "Otu110", "Otu111", "Otu112", "Otu114", "Otu115", "Otu119", "Otu123", "Otu125", "Otu129", "Otu130", "Otu135", "Otu136", "Otu137", "Otu146", "Otu147", "Otu150", "Otu152", "Otu156", "Otu158", "Otu163", "Otu164", "Otu165", "Otu171", "Otu173", "Otu185", "Otu252")
# Test Biomarkers to check if code works
#biomarkers <- c("MMP.9", "S100A8", "BAFF","IFN.alpha", "TIMP.1", "IL1b", "IL6", "MIP1a", "MMP8", "Adiponectin")


# Define biomarkers and confounders
confounders <- c("Age", "BMI", "Total_sites")  # Confounders to be used in the model
#confounders <- NULL