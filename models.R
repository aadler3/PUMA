# PUMA Analysis - Statistical Models
# Author: PUMA Research Team
# Date: 2025

source("utils.R")

# Run GEE model for hair analysis
run_hair_model <- function(hair_gee_data, scale_method = "sd", include_group = TRUE) {
  if (scale_method == "sd") {
    predictor <- "wilson_adh_scaled_sd"
  } else if (scale_method == "10") {
    predictor <- "wilson_adh10" 
  } else {
    predictor <- "wilson_adh"
  }
  
  formula <- if (include_group) {
    as.formula(paste("detect_tfv ~", predictor, "+ group_new"))
  } else {
    as.formula(paste("detect_tfv ~", predictor))
  }
  
  model <- gee(formula,
               id = ptid,
               family = poisson(link = "log"),
               corstr = "exchangeable",
               data = hair_gee_data)
  
  return(model)
}

# Run GEE model for urine analysis
run_urine_model <- function(ut_gee_data, scale_method = "sd", include_group = TRUE) {
  if (scale_method == "sd") {
    predictor <- "wilson_adh_scaled_sd"
  } else if (scale_method == "10") {
    predictor <- "wilson_adh10" 
  } else {
    predictor <- "wilson_adh"
  }
  
  formula <- if (include_group) {
    as.formula(paste("ut_result ~", predictor, "+ group_new"))
  } else {
    as.formula(paste("ut_result ~", predictor))
  }
  
  model <- gee(formula,
               id = ptid,
               family = poisson(link = "log"),
               corstr = "exchangeable",
               data = ut_gee_data)
  
  return(model)
}

# Run adjusted models with covariates
run_adjusted_model <- function(data, outcome, predictors, id_var = "ptid") {
  # Create formula from outcome and predictor variables
  formula_str <- paste(outcome, "~", paste(predictors, collapse = " + "))
  formula_obj <- as.formula(formula_str)
  
  # Run the GEE model
  model <- gee(formula_obj,
               id = data[[id_var]],
               family = poisson(link = "log"),
               corstr = "exchangeable",
               data = data)
  
  return(model)
}

# Run model selection for urine data
run_urine_model_selection <- function(ut_multivar) {
  models <- list()
  
  # Start with minimum model
  models$min_model <- geeglm(ut_result ~ group_new,
                         id = ptid,
                         family = poisson(link = "log"),
                         corstr = "exchangeable",
                         data = ut_multivar)
  
  # Add in Wilson ADH
  models$wilson_model <- geeglm(ut_result ~ group_new + wilson_adh_scaled_sd,
                            id = ptid,
                            family = poisson(link = "log"),
                            corstr = "exchangeable",
                            data = ut_multivar)
  
  # Add in employment status 
  models$emp_model <- geeglm(ut_result ~ group_new + wilson_adh_scaled_sd + emp_coded,
                         id = ptid,
                         family = poisson(link = "log"),
                         corstr = "exchangeable",
                         data = ut_multivar)
  
  # Remove emp and add in age
  models$age_model <- geeglm(ut_result ~ group_new + wilson_adh_scaled_sd + age,
                         id = ptid,
                         family = poisson(link = "log"),
                         corstr = "exchangeable",
                         data = ut_multivar)
  
  # Add all in
  models$full_model <- geeglm(ut_result ~ group_new + emp_coded + age + wilson_adh_scaled_sd,
                          id = ptid,
                          family = poisson(link = "log"),
                          corstr = "exchangeable",
                          data = ut_multivar)
  
  # Calculate QIC for all models
  qic_values <- sapply(models, QIC)
  
  # Find the best model
  best_model_name <- names(which.min(qic_values))
  
  results <- list(
    models = models,
    qic_values = qic_values,
    best_model_name = best_model_name,
    best_model = models[[best_model_name]]
  )
  
  return(results)
}

# Run model selection for hair data
run_hair_model_selection <- function(hr_multivar) {
  models <- list()
  
  # Start with minimum model
  models$min_model <- geeglm(detect_tfv ~ group_new,
                         id = ptid,
                         family = poisson(link = "log"),
                         corstr = "exchangeable",
                         data = hr_multivar)
  
  # Add in Wilson ADH
  models$wilson_model <- geeglm(detect_tfv ~ group_new + wilson_adh_scaled_sd,
                            id = ptid,
                            family = poisson(link = "log"),
                            corstr = "exchangeable",
                            data = hr_multivar)
  
  # Add in edu
  models$edu_model <- geeglm(detect_tfv ~ group_new + wilson_adh_scaled_sd + edu,
                         id = ptid,
                         family = poisson(link = "log"),
                         corstr = "exchangeable",
                         data = hr_multivar)
  
  # Just edu by itself
  models$edu_model2 <- geeglm(detect_tfv ~ group_new + edu,
                          id = ptid,
                          family = poisson(link = "log"),
                          corstr = "exchangeable",
                          data = hr_multivar)
  
  # Calculate QIC for all models
  qic_values <- sapply(models, QIC)
  
  # Find the best model
  best_model_name <- names(which.min(qic_values))
  
  results <- list(
    models = models,
    qic_values = qic_values,
    best_model_name = best_model_name,
    best_model = models[[best_model_name]]
  )
  
  return(results)
}

# Run univariable analysis
run_univariable_analysis <- function(data, predictors, outcomes) {
  # Create a data frame to store results
  results <- data.frame(
    outcome = character(),
    predictor = character(),
    estimate = numeric(),
    std_error = numeric(),
    p_value = numeric(),
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
        # Fit GEE model with Poisson regression
        model <- geeglm(formula_obj, 
                      data = data,
                      family = poisson(link = "log"),
                      id = ptid,
                      corstr = "exchangeable")
        
        # Extract model summary
        model_summary <- summary(model)
        
        # Get coefficient, std error and p-value for the predictor
        coef_data <- coef(model_summary)
        
        # Skip intercept, get row for predictor
        if (nrow(coef_data) > 1) {
          predictor_row <- coef_data[2, ]
          
          # Store results
          results <- rbind(results, data.frame(
            outcome = outcome,
            predictor = predictor,
            estimate = predictor_row[1],
            std_error = predictor_row[2],
            p_value = predictor_row[4],
            stringsAsFactors = FALSE
          ))
        }
      }, error = function(e) {
        cat("Error in model for outcome:", outcome, "and predictor:", predictor, 
            "-", e$message, "\n")
      })
    }
  }
  
  return(results)
}

# Run ROC analysis
run_roc_analysis <- function(data, outcome_var, predictor_var) {
  # Create ROC object
  roc_obj <- roc(data[[outcome_var]], data[[predictor_var]])
  
  # Get coordinates for different cutoffs
  coords_all <- coords(roc_obj, "all", ret = c("threshold", "sensitivity", "specificity"))
  
  # Calculate AUC
  auc_value <- auc(roc_obj)
  
  # Find optimal cutoff using Youden's index
  best_cutoff <- coords(roc_obj, "best", ret = "threshold")
  
  # Create a sequence of cutoffs
  cutoffs <- seq(min(data[[predictor_var]], na.rm = TRUE), 
                max(data[[predictor_var]], na.rm = TRUE), 
                length.out = 20)
  
  # Pre-allocate vectors
  sensitivity_values <- numeric(length(cutoffs))
  specificity_values <- numeric(length(cutoffs))
  
  # Calculate sensitivity and specificity for each cutoff
  for(i in 1:length(cutoffs)) {
    cutoff <- cutoffs[i]
    # Get coordinates for this specific threshold
    coords_at_cutoff <- coords(roc_obj, cutoffs[i], input = "threshold")
    
    # Store values
    sensitivity_values[i] <- coords_at_cutoff$sensitivity
    specificity_values[i] <- coords_at_cutoff$specificity
  }
  
  # Create final dataframe
  results_df <- data.frame(
    cutoff = cutoffs,
    sensitivity = sensitivity_values,
    specificity = specificity_values
  )
  
  return(list(
    roc_obj = roc_obj,
    auc = auc_value,
    best_cutoff = best_cutoff,
    coords_all = coords_all,
    results_df = results_df
  ))
} 