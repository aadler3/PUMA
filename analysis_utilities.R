# Analysis utilities for PUMA study
# This file contains functions for analyzing PUMA data

#' Run descriptive statistics on adherence scores by binary biomarker
#' 
#' @param data The dataset
#' @param adherence_var Name of the adherence variable
#' @param biomarker_var Name of the binary biomarker variable
#' @return A list with t-test results and summary statistics
summarize_adherence_by_biomarker <- function(data, adherence_var = "wilson_adh", biomarker_var) {
  # Run t-test
  formula <- as.formula(paste(adherence_var, "~", biomarker_var))
  t_test_result <- t.test(formula, data = data)
  
  # Generate summaries
  pos_data <- data %>% filter(get(biomarker_var) == 1)
  neg_data <- data %>% filter(get(biomarker_var) == 0)
  
  pos_mean <- mean(pos_data[[adherence_var]], na.rm = TRUE)
  pos_sd <- sd(pos_data[[adherence_var]], na.rm = TRUE)
  pos_n <- sum(!is.na(pos_data[[adherence_var]]))
  
  neg_mean <- mean(neg_data[[adherence_var]], na.rm = TRUE)
  neg_sd <- sd(neg_data[[adherence_var]], na.rm = TRUE)
  neg_n <- sum(!is.na(neg_data[[adherence_var]]))
  
  overall_mean <- mean(data[[adherence_var]], na.rm = TRUE)
  overall_sd <- sd(data[[adherence_var]], na.rm = TRUE)
  overall_n <- sum(!is.na(data[[adherence_var]]))
  
  # Return results
  return(list(
    t_test = t_test_result,
    positive = list(mean = pos_mean, sd = pos_sd, n = pos_n),
    negative = list(mean = neg_mean, sd = neg_sd, n = neg_n),
    overall = list(mean = overall_mean, sd = overall_sd, n = overall_n)
  ))
}

#' Run GEE model and extract results including relative risks and confidence intervals
#' 
#' @param data The dataset
#' @param outcome The outcome variable (binary)
#' @param predictors Character vector of predictor variables
#' @param id_var The ID variable for clustering
#' @param corstr Correlation structure (default: "exchangeable")
#' @param family Model family (default: poisson with log link)
#' @param round_digits Number of digits to round results
#' @return A data frame with model results
run_gee_model <- function(data, outcome, predictors, 
                          id_var = "ptid", 
                          corstr = "exchangeable",
                          family = poisson(link = "log"),
                          round_digits = 3) {
  
  # Create formula
  formula_str <- paste(outcome, "~", paste(predictors, collapse = " + "))
  formula <- as.formula(formula_str)
  
  # Run GEE model
  gee_model <- gee(formula,
                   id = data[[id_var]],
                   family = family,
                   corstr = corstr,
                   data = data)
  
  # Extract model summary
  model_summary <- summary(gee_model)
  
  # Calculate p-values from z-values
  z_values <- model_summary$coefficients[, "Robust z"]
  p_values <- 2 * pnorm(-abs(z_values))
  
  # Calculate relative risks and confidence intervals
  coefficients <- model_summary$coefficients[, 1]
  relative_risks <- exp(coefficients)
  robust_ses <- model_summary$coefficients[, "Robust S.E."]
  lower_ci <- exp(coefficients - 1.96 * robust_ses)
  upper_ci <- exp(coefficients + 1.96 * robust_ses)
  
  # Create results data frame
  terms <- rownames(model_summary$coefficients)
  
  results <- data.frame(
    Term = terms,
    Estimate = round(relative_risks, round_digits),
    Robust_SE = round(robust_ses, round_digits),
    z_value = round(z_values, round_digits),
    p_value = round(p_values, round_digits),
    CI_Lower = round(lower_ci, round_digits),
    CI_Upper = round(upper_ci, round_digits)
  )
  
  return(results)
}

#' Filter data for GEE analysis
#' 
#' @param data The dataset
#' @param outcome The outcome variable
#' @param predictors The predictor variables
#' @param id_var The ID variable
#' @param additional_vars Additional variables to include (optional)
#' @return A filtered dataset
prepare_gee_data <- function(data, outcome, predictors, id_var = "ptid", additional_vars = NULL) {
  all_vars <- c(outcome, predictors, id_var, additional_vars)
  
  # Filter complete cases for the required variables
  filtered_data <- data %>%
    filter(complete.cases(across(all_of(all_vars)))) %>%
    select(all_of(all_vars))
  
  return(filtered_data)
}

#' Run GEE model for continuous outcomes and extract results
#' 
#' @param data The dataset
#' @param outcome The continuous outcome variable
#' @param predictors Character vector of predictor variables
#' @param id_var The ID variable for clustering
#' @param corstr Correlation structure (default: "exchangeable")
#' @param family Model family (default: gaussian with identity link)
#' @param round_digits Number of digits to round results
#' @return A data frame with model results
run_gee_continuous_model <- function(data, outcome, predictors, 
                                     id_var = "ptid", 
                                     corstr = "exchangeable",
                                     family = gaussian(link = "identity"),
                                     round_digits = 3) {
  
  # Create formula
  formula_str <- paste(outcome, "~", paste(predictors, collapse = " + "))
  formula <- as.formula(formula_str)
  
  # Run GEE model
  gee_model <- gee(formula,
                   id = data[[id_var]],
                   family = family,
                   corstr = corstr,
                   data = data)
  
  # Extract model summary
  model_summary <- summary(gee_model)
  
  # Calculate p-values from z-values
  z_values <- model_summary$coefficients[, "Robust z"]
  p_values <- 2 * pnorm(-abs(z_values))
  
  # Extract coefficients and calculate confidence intervals
  coefficients <- model_summary$coefficients[, 1]
  robust_ses <- model_summary$coefficients[, "Robust S.E."]
  lower_ci <- coefficients - 1.96 * robust_ses
  upper_ci <- coefficients + 1.96 * robust_ses
  
  # Create results data frame
  terms <- rownames(model_summary$coefficients)
  
  results <- data.frame(
    Term = terms,
    Coefficient = round(coefficients, round_digits),
    Robust_SE = round(robust_ses, round_digits),
    z_value = round(z_values, round_digits),
    p_value = round(p_values, round_digits),
    CI_Lower = round(lower_ci, round_digits),
    CI_Upper = round(upper_ci, round_digits)
  )
  
  return(results)
}