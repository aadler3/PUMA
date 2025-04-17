# PUMA Analysis Utility Functions
# Author: PUMA Research Team
# Date: 2025

# Load necessary packages
load_packages <- function() {
  packages <- c(
    "dplyr", "flextable", "ggplot2", "collapse", "labelled", 
    "compareGroups", "gee", "geepack", "pROC", "ROCR", 
    "plotROC", "gridExtra"
  )
  
  invisible(lapply(packages, function(pkg) {
    if (!require(pkg, character.only = TRUE)) {
      install.packages(pkg)
      library(pkg, character.only = TRUE)
    }
  }))
}

# Set data pathway
get_data_path <- function() {
  "~/Library/CloudStorage/Box-Box/Puma Analysis/Data/"
}

# Load PUMA data
load_puma_data <- function(path = get_data_path()) {
  file_path <- file.path(path, "PUMAextract_13Mar2025.csv")
  read.csv(file_path)
}

# Format model results for output
format_gee_results <- function(model, exponentiate = TRUE) {
  model_summary <- summary(model)
  
  # Calculate z values and p-values
  z_values <- model_summary$coefficients[, "Robust z"]
  p_values <- 2 * pnorm(-abs(z_values))
  
  # Get coefficients and standard errors
  coefficients <- model_summary$coefficients[, 1]
  robust_ses <- model_summary$coefficients[, "Robust S.E."]
  
  if (exponentiate) {
    # For Poisson models (RR)
    estimates <- exp(coefficients)
    lower_ci <- exp(coefficients - 1.96 * robust_ses)
    upper_ci <- exp(coefficients + 1.96 * robust_ses)
  } else {
    # For linear models
    estimates <- coefficients
    lower_ci <- coefficients - 1.96 * robust_ses
    upper_ci <- coefficients + 1.96 * robust_ses
  }
  
  # Create results data frame
  results <- data.frame(
    Term = rownames(model_summary$coefficients),
    Estimate = round(estimates, 3),
    Robust_SE = round(robust_ses, 3),
    z_value = round(z_values, 3),
    p_value = round(p_values, 3),
    CI_Lower = round(lower_ci, 3),
    CI_Upper = round(upper_ci, 3)
  )
  
  return(results)
} 