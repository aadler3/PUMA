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

#' Run univariable GEE models for multiple predictors and outcomes
#' 
#' @param data The dataset
#' @param predictors Character vector of predictor variables
#' @param outcomes Character vector of outcome variables
#' @param id_var The ID variable for clustering
#' @param family Model family (default: poisson with log link)
#' @param corstr Correlation structure (default: "exchangeable")
#' @param print_progress Whether to print progress during analysis
#' @param round_digits Number of digits to round results
#' @return A data frame with all univariable model results
run_univariable_analysis <- function(data, predictors, outcomes, 
                                     id_var = "ptid",
                                     family = poisson(link = "log"),
                                     corstr = "exchangeable",
                                     print_progress = TRUE,
                                     round_digits = 4) {
  
  # Create a data frame to store results
  results <- data.frame(
    outcome = character(),
    predictor = character(),
    estimate = numeric(),
    exp_estimate = numeric(),
    std_error = numeric(),
    p_value = numeric(),
    significant = logical(),
    stringsAsFactors = FALSE
  )
  
  # Loop through each outcome and predictor combination
  for (outcome in outcomes) {
    for (predictor in predictors) {
      
      # Create formula for the model
      formula_str <- paste(outcome, "~", predictor)
      formula_obj <- as.formula(formula_str)
      
      # Try to fit the model, handle potential errors
      tryCatch({
        # Fit GEE model
        model <- geeglm(formula_obj, 
                        data = data,
                        family = family,
                        id = data[[id_var]],
                        corstr = corstr)
        
        # Extract model summary
        model_summary <- summary(model)
        
        # Get coefficient, std error and p-value for the predictor
        coef_data <- coef(model_summary)
        
        # Skip intercept, get row for predictor
        if (nrow(coef_data) > 1) {
          predictor_row <- coef_data[2, ]
          
          # Determine if result is statistically significant
          is_significant <- predictor_row[4] < 0.05
          
          # Store results
          results <- rbind(results, data.frame(
            outcome = outcome,
            predictor = predictor,
            estimate = round(predictor_row[1], round_digits),
            exp_estimate = round(exp(predictor_row[1]), round_digits),
            std_error = round(predictor_row[2], round_digits),
            p_value = round(predictor_row[4], round_digits),
            significant = is_significant,
            stringsAsFactors = FALSE
          ))
          
          # Print progress if requested
          if (print_progress) {
            cat("Outcome:", outcome, "| Predictor:", predictor, 
                "| Estimate:", round(predictor_row[1], round_digits),
                "| Exp(Estimate):", round(exp(predictor_row[1]), round_digits),
                "| p-value:", round(predictor_row[4], round_digits),
                "| Significant:", is_significant, "\n")
          }
        }
      }, error = function(e) {
        if (print_progress) {
          cat("Error in model for outcome:", outcome, "and predictor:", predictor, "-", e$message, "\n")
        }
      })
    }
  }
  
  # Return the results data frame
  return(results)
}

#' Summarize univariable analysis results
#'
#' @param univariable_results Results from run_univariable_analysis
#' @param significance_threshold P-value threshold for significance
#' @return A list with significant predictors by outcome
summarize_univariable_results <- function(univariable_results, significance_threshold = 0.05) {
  # Create a list to store significant predictors by outcome
  significant_predictors <- list()
  
  # Get unique outcomes
  outcomes <- unique(univariable_results$outcome)
  
  # For each outcome, find significant predictors
  for (outcome in outcomes) {
    outcome_results <- univariable_results[univariable_results$outcome == outcome, ]
    significant <- outcome_results[outcome_results$p_value < significance_threshold, ]
    
    if (nrow(significant) > 0) {
      significant_predictors[[outcome]] <- significant$predictor
      
      # Print summary
      cat("\nSignificant predictors for", outcome, ":\n")
      for (i in 1:nrow(significant)) {
        cat("- ", significant$predictor[i], 
            "(p =", significant$p_value[i], 
            ", estimate =", significant$estimate[i], 
            ", exp(estimate) =", significant$exp_estimate[i], ")\n")
      }
    } else {
      significant_predictors[[outcome]] <- character(0)
      cat("\nNo significant predictors found for", outcome, "\n")
    }
  }
  
  return(significant_predictors)
}

#' Run a GEE model for a specific study group
#' 
#' @param data The dataset
#' @param outcome The outcome variable
#' @param predictor The main predictor variable
#' @param group_var The group variable name
#' @param group_value The specific group value to analyze
#' @param id_var ID variable for clustering
#' @return A data frame with model results formatted exactly as in original code
run_group_specific_gee <- function(data, outcome, predictor, group_var, group_value,
                                   id_var = "ptid") {
  
  # Filter data for this group
  group_data <- data[data[[group_var]] == group_value, ]
  
  # Create formula
  formula_str <- paste(outcome, "~", predictor)
  formula <- as.formula(formula_str)
  
  # Run GEE model exactly as in original code
  group_gee <- gee(formula,
                   id = group_data[[id_var]],
                   family = poisson(link = "log"),
                   corstr = "exchangeable",
                   data = group_data)
  
  # Extract summary
  summary_group <- summary(group_gee)
  
  # Calculate p-values from z-values exactly as in original code
  z_group <- summary_group$coefficients[, "Robust z"]
  p_group <- 2*pnorm(-abs(z_group))
  
  # Calculate RR and 95% CIs exactly as in original code
  coefficients_group <- summary(group_gee)$coefficients[, 1]
  relative_risks_group <- exp(summary_group$coefficients[, 1])
  robust_ses_group <- summary_group$coefficients[, "Robust S.E."]
  
  # The original code calculates CI differently - using the exact same approach
  lower_ci_group <- exp(coefficients_group - 1.96 * robust_ses_group)
  upper_ci_group <- exp(coefficients_group + 1.96 * robust_ses_group)
  
  names_group <- rownames(summary_group$coefficients)
  
  # Create results data frame exactly as in original code
  group_results <- data.frame(
    Term = names_group,
    Estimate = round(exp(coefficients_group), 3), 
    Robust_SE = round(robust_ses_group, 3),
    z_value = round(z_group, 3),
    p_value = round(p_group, 3),
    CI_Lower = round(lower_ci_group, 3),
    CI_Upper = round(upper_ci_group, 3)
  )
  
  return(group_results)
}

#' Run interaction test between study groups
#' 
#' @param data The dataset
#' @param outcome The outcome variable
#' @param predictor The main predictor variable
#' @param group_var The group variable name
#' @param id_var ID variable for clustering
#' @return A data frame with interaction model results
run_group_interaction_test <- function(data, outcome, predictor, group_var,
                                       id_var = "ptid") {
  
  # Create interaction formula exactly as in original code
  formula_str <- paste(outcome, "~", predictor, "+", group_var, "+", 
                       group_var, "*", predictor)
  formula <- as.formula(formula_str)
  
  # Run GEE model
  lrt_model <- gee(formula,
                   id = data[[id_var]],
                   family = poisson(link = "log"),
                   corstr = "exchangeable",
                   data = data)
  
  # Extract summary
  summary_lrt <- summary(lrt_model)
  
  # Calculate p-values exactly as in original code
  z_lrt <- summary_lrt$coefficients[, "Robust z"]
  p_lrt <- 2*pnorm(-abs(z_lrt))
  
  # Calculate estimates exactly as in original code
  coefficients_lrt <- summary(lrt_model)$coefficients[, 1]
  relative_risks_lrt <- exp(summary_lrt$coefficients[, 1])
  robust_ses_lrt <- summary_lrt$coefficients[, "Robust S.E."]
  lower_ci_lrt <- exp(coefficients_lrt - 1.96 * robust_ses_lrt)
  upper_ci_lrt <- exp(coefficients_lrt + 1.96 * robust_ses_lrt)
  names_lrt <- rownames(summary_lrt$coefficients)
  
  # Create results data frame exactly as in original code
  lrt_results <- data.frame(
    Term = names_lrt,
    Estimate = round(exp(coefficients_lrt), 3), 
    Robust_SE = round(robust_ses_lrt, 3),
    z_value = round(z_lrt, 3),
    p_value = round(p_lrt, 3),
    CI_Lower = round(lower_ci_lrt, 3),
    CI_Upper = round(upper_ci_lrt, 3)
  )
  
  return(lrt_results)
}