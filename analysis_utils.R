# Analysis utilities for PUMA study
# This file contains functions for analyzing PUMA data

# Load required packages for analysis
library(dplyr)
library(geepack)
library(pROC)
library(car)
library(compareGroups)
library(collapse)

#' Run t-test for comparison between groups
#' 
#' @param data Data frame containing variables to test
#' @param outcome_var Name of the outcome variable
#' @param group_var Name of the grouping variable
#' @param na.rm Whether to remove NA values, default is TRUE
#' @return List containing t-test results, means and SDs
#' @export
compare_groups_ttest <- function(data, outcome_var, group_var, na.rm = TRUE) {
  # Run t-test
  t_result <- t.test(formula = as.formula(paste(outcome_var, "~", group_var)), 
                     data = data, na.rm = na.rm)
  
  # Calculate mean and SD for each group
  results <- data %>%
    group_by(across(all_of(group_var))) %>%
    summarize(
      mean = mean(get(outcome_var), na.rm = na.rm),
      sd = sd(get(outcome_var), na.rm = na.rm),
      n = sum(!is.na(get(outcome_var)))
    )
  
  return(list(
    t_test = t_result,
    summary = results,
    p_value = t_result$p.value
  ))
}

#' Run chi-square test for categorical variables
#'
#' @param data Data frame containing variables
#' @param var1 First variable name
#' @param var2 Second variable name
#' @return Chi-square test results
#' @export
run_chisq_test <- function(data, var1, var2) {
  # Create table
  table_data <- table(data[[var1]], data[[var2]])
  
  # Run chi-square test
  chisq_result <- chisq.test(table_data)
  
  return(chisq_result)
}

#' Create detailed chi-square analysis for all visits
#'
#' @param wide_data Data in wide format with visit columns
#' @param outcome_var Base name of outcome variable (e.g., "ut_result")
#' @param group_var Grouping variable (default: "group_new")
#' @param visit_count Number of visits to analyze
#' @return List of chi-square test results
#' @export
run_visit_chisq_analysis <- function(wide_data, outcome_var, group_var = "group_new", visit_count = 5) {
  results <- list()
  
  # Group comparisons for each visit
  for (i in 1:visit_count) {
    var_name <- paste0(outcome_var, ".", i)
    
    # Check if the variable exists in the data
    if (var_name %in% names(wide_data)) {
      # Create table
      table_data <- table(wide_data[[group_var]], wide_data[[var_name]])
      
      # Run chi-square test
      chisq_result <- chisq.test(table_data)
      
      results[[paste0("visit", i, "_by_group")]] <- chisq_result
      
      # Run test on just the outcome distribution
      observed_data <- table(wide_data[[var_name]])
      if (length(observed_data) > 1) {  # Need at least 2 categories
        overall_result <- chisq.test(observed_data)
        results[[paste0("visit", i, "_overall")]] <- overall_result
      }
    }
  }
  
  return(results)
}

#' Create a GEE model for longitudinal data
#'
#' @param data Longitudinal dataset in long format
#' @param formula Model formula as a string
#' @param id_var Subject identifier variable
#' @param family Model family (default: binomial with logit link)
#' @param corstr Correlation structure (default: independence)
#' @return GEE model object
#' @export
run_gee_model <- function(data, formula, id_var = "ptid", 
                          family = binomial(link = "logit"), 
                          corstr = "independence") {
  # Convert formula from string to formula
  model_formula <- as.formula(formula)
  
  # Run GEE model
  gee_model <- geeglm(
    formula = model_formula,
    family = family,
    id = as.formula(paste("~", id_var)),
    corstr = corstr,
    data = data
  )
  
  return(gee_model)
}

#' Test interaction terms in a GEE model
#'
#' @param gee_model GEE model object
#' @param visit_levels Visit levels to test (default: 2:5)
#' @param group_var Grouping variable name (default: "group_new")
#' @param group_value Value of the intervention group (default: "Intervention")
#' @return Data frame with p-values for each interaction
#' @export
test_gee_interactions <- function(gee_model, visit_levels = 2:5, 
                                 group_var = "group_new", 
                                 group_value = "Intervention") {
  
  p_values <- numeric(length(visit_levels))
  
  for (i in seq_along(visit_levels)) {
    k <- visit_levels[i]
    # Create the term to test
    interaction_term <- paste0("factor(visit)", k, ":factor(", group_var, ")", group_value)
    
    # Get linearHypothesis test results for this term
    test_result <- car::linearHypothesis(gee_model, interaction_term)
    
    # Extract p-value
    p_values[i] <- test_result[["Pr(>Chisq)"]][2]
  }
  
  # Create dataframe with results
  result_df <- data.frame(
    visit = visit_levels,
    p_value = p_values,
    formatted_p = sprintf("%.3f", p_values)
  )
  
  return(result_df)
}

#' Create ROC curve analysis
#'
#' @param data Data frame
#' @param predictor Predictor variable (e.g., wilson_adh)
#' @param outcome Binary outcome variable (e.g., ut_result)
#' @param cutoffs Optional vector of cutoff values
#' @param find_best_cutoff Whether to calculate the best cutoff (default: TRUE)
#' @param find_sens90_cutoff Whether to find cutoff with sensitivity closest to 90% (default: FALSE)
#' @return List containing ROC object and dataframe of sensitivity/specificity at cutoffs
#' @export
create_roc_analysis <- function(data, predictor, outcome, cutoffs = NULL, 
                               find_best_cutoff = TRUE, 
                               find_sens90_cutoff = FALSE) {
  # Create ROC object
  roc_obj <- roc(data[[outcome]], data[[predictor]])
  
  # Calculate AUC
  auc_value <- auc(roc_obj)
  
  # Get coordinates for all thresholds
  coords_all <- coords(roc_obj, "all", ret = c("threshold", "sensitivity", "specificity"))
  
  # Find best cutoff using Youden's J statistic if requested
  best_cutoff <- NULL
  if (find_best_cutoff) {
    best_cutoff <- coords(roc_obj, "best", ret = "threshold")
  }
  
  # Find cutoff with sensitivity closest to 90% if requested
  sens90_cutoff <- NULL
  if (find_sens90_cutoff) {
    coords_df <- as.data.frame(coords(roc_obj, "all"))
    closest_to_90 <- coords_df[which.min(abs(coords_df$sensitivity - 0.9)),]
    sens90_cutoff <- closest_to_90$threshold
  }
  
  # If cutoffs are provided, calculate sensitivity/specificity at those cutoffs
  results_df <- NULL
  if (!is.null(cutoffs)) {
    sensitivity_values <- numeric(length(cutoffs))
    specificity_values <- numeric(length(cutoffs))
    
    for(i in 1:length(cutoffs)) {
      coords_at_cutoff <- coords(roc_obj, cutoffs[i], input = "threshold")
      sensitivity_values[i] <- coords_at_cutoff$sensitivity
      specificity_values[i] <- coords_at_cutoff$specificity
    }
    
    # Create final dataframe
    results_df <- data.frame(
      cutoff = round(cutoffs, 3),
      sensitivity = round(sensitivity_values, 3),
      specificity = round(specificity_values, 3)
    )
  }
  
  return(list(
    roc = roc_obj,
    auc = auc_value,
    best_cutoff = best_cutoff,
    sens90_cutoff = sens90_cutoff,
    coords_all = coords_all,
    results = results_df
  ))
}

#' Plot ROC curves with AUC values
#'
#' @param roc_list List of ROC objects with names
#' @param colors Vector of colors for each ROC curve
#' @param main Plot title
#' @param print_auc Whether to print AUC on the plot
#' @return ggplot object
#' @export
plot_roc_curves <- function(roc_list, colors, main = "ROC Curves", print_auc = FALSE) {
  # Create initial plot with first ROC curve
  plot(roc_list[[1]]$roc, col = colors[1], main = main, lwd = 2, print.auc = print_auc)
  
  # Add additional ROC curves
  if (length(roc_list) > 1) {
    for (i in 2:length(roc_list)) {
      plot(roc_list[[i]]$roc, col = colors[i], add = TRUE, lwd = 2, print.auc = print_auc)
    }
  }
  
  # Create legend text with AUC values
  legend_text <- sapply(1:length(roc_list), function(i) {
    paste0(names(roc_list)[i], " (AUC = ", round(roc_list[[i]]$auc, 3), ")")
  })
  
  # Add legend
  legend("bottomright", legend = legend_text, col = colors, lwd = 2)
}

#' Calculate descriptive statistics for a numeric variable
#'
#' @param data Data frame
#' @param var Variable name
#' @param na.rm Whether to remove NAs (default: TRUE)
#' @return List of statistics
#' @export
describe_numeric <- function(data, var, na.rm = TRUE) {
  list(
    mean = mean(data[[var]], na.rm = na.rm),
    sd = sd(data[[var]], na.rm = na.rm),
    median = median(data[[var]], na.rm = na.rm),
    min = min(data[[var]], na.rm = na.rm),
    max = max(data[[var]], na.rm = na.rm),
    n = sum(!is.na(data[[var]]))
  )
}

#' Create comprehensive summary statistics for all variables in a dataset
#'
#' @param data Data frame to summarize
#' @param group_var Grouping variable (default: NULL for overall stats)
#' @param strata_vars Additional stratification variables (default: NULL)
#' @param show.all Show all variables (default: TRUE)
#' @param show.p.overall Show overall p-value (default: TRUE)
#' @return CompareGroups table
#' @export
create_complete_summary <- function(data, group_var = NULL, strata_vars = NULL,
                                  show.all = TRUE, show.p.overall = TRUE) {
  # Create formula based on parameters
  if (is.null(group_var)) {
    formula <- ~.
  } else {
    formula <- as.formula(paste(group_var, "~", "."))
  }
  
  # Create compareGroups object
  cg <- compareGroups(formula, data = data, subset = NULL)
  
  # Add stratification if provided
  if (!is.null(strata_vars)) {
    cg <- strataTable(cg, strata = strata_vars)
  }
  
  # Create table
  tab <- createTable(cg, show.all = show.all, show.p.overall = show.p.overall)
  
  return(tab)
}

#' Analyze BLQ (Below Limit of Quantification) data
#'
#' @param data Data frame with BLQ data
#' @param blq_var BLQ variable name
#' @param visit_var Visit variable name
#' @param group_var Group variable name
#' @return List of tables and test results
#' @export
analyze_blq_data <- function(data, blq_var = "blq", visit_var = "visit_coded", group_var = "group_new") {
  # Create tables by visit for overall data
  overall_table <- table(data[[visit_var]], data[[blq_var]])
  
  # Create tables by visit for each group
  group_tables <- list()
  groups <- unique(data[[group_var]])
  
  for (group in groups) {
    group_data <- data[data[[group_var]] == group, ]
    group_tables[[as.character(group)]] <- table(group_data[[visit_var]], group_data[[blq_var]])
  }
  
  # Run chi-square tests for each visit
  visits <- unique(data[[visit_var]])
  chisq_results <- list()
  
  for (visit in visits) {
    visit_data <- data[data[[visit_var]] == visit, ]
    if (nrow(visit_data) > 0) {
      test_table <- table(visit_data[[group_var]], visit_data[[blq_var]])
      if (sum(test_table) > 0 && min(dim(test_table)) > 1) {
        chisq_results[[as.character(visit)]] <- chisq.test(test_table)
      }
    }
  }
  
  return(list(
    overall_table = overall_table,
    group_tables = group_tables,
    chisq_results = chisq_results
  ))
}

#' Analyze longitudinal data by visit and group
#'
#' @param data Data frame in wide format
#' @param var_prefix Variable prefix (e.g., "wilson_adh")
#' @param group_var Grouping variable (default: "group_new")
#' @param visit_count Number of visits (default: 5)
#' @param na.rm Whether to remove NAs (default: TRUE)
#' @return List of summary statistics and test results
#' @export
analyze_longitudinal_var <- function(data, var_prefix, group_var = "group_new", 
                                   visit_count = 5, na.rm = TRUE) {
  results <- list()
  
  # Loop through each visit
  for (i in 1:visit_count) {
    var_name <- paste0(var_prefix, ".", i)
    
    # Check if variable exists
    if (var_name %in% names(data)) {
      # Overall statistics
      results[[paste0("visit", i, "_overall")]] <- list(
        mean = mean(data[[var_name]], na.rm = na.rm),
        sd = sd(data[[var_name]], na.rm = na.rm),
        n = sum(!is.na(data[[var_name]]))
      )
      
      # Run t-test between groups
      if (length(unique(data[[group_var]])) > 1) {
        formula <- as.formula(paste(var_name, "~", group_var))
        results[[paste0("visit", i, "_ttest")]] <- t.test(formula, data = data)
        
        # Group statistics
        group_stats <- data %>%
          group_by(across(all_of(group_var))) %>%
          summarize(
            mean = mean(get(var_name), na.rm = na.rm),
            sd = sd(get(var_name), na.rm = na.rm),
            n = sum(!is.na(get(var_name))),
            missing = sum(is.na(get(var_name)))
          )
        
        results[[paste0("visit", i, "_by_group")]] <- group_stats
      }
    }
  }
  
  return(results)
}

#' Create a combined model with multiple predictors
#'
#' @param data Data frame
#' @param outcome Outcome variable name
#' @param predictors Vector of predictor variable names
#' @param family Model family (default: binomial with logit link)
#' @return Model object
#' @export
create_combined_model <- function(data, outcome, predictors, 
                                family = binomial(link = "logit")) {
  # Create formula string
  formula_str <- paste(outcome, "~", paste(predictors, collapse = " + "))
  
  # Create model
  model <- glm(formula = as.formula(formula_str), family = family, data = data)
  
  return(model)
} 