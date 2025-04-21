# Utility functions for ROC curve analysis

#' Create a ROC curve object and calculate AUC
#' 
#' @param response The binary response variable
#' @param predictor The continuous predictor variable
#' @param direction Direction parameter for the ROC curve
#' @return A list containing the ROC object and AUC value
create_roc <- function(response, predictor, direction = "auto") {
  roc_obj <- pROC::roc(response, predictor, direction = direction)
  auc_value <- pROC::auc(roc_obj)
  
  return(list(
    roc = roc_obj,
    auc = auc_value,
    auc_formatted = round(auc_value, 3)
  ))
}

#' Create a legend string with AUC value
#' 
#' @param label The label for the curve
#' @param auc_value The AUC value to include
#' @return A formatted legend string
create_legend_with_auc <- function(label, auc_value) {
  paste0(label, " (AUC = ", round(auc_value, 3), ")")
}

#' Plot multiple ROC curves with AUC values in legend
#' 
#' @param roc_list A list of ROC objects
#' @param labels A vector of labels for the legend
#' @param colors A vector of colors for the curves
#' @param main Main title for the plot
#' @param lwd Line width
plot_multiple_roc <- function(roc_list, labels, colors, main = "ROC Curves", lwd = 2) {
  # Plot first ROC curve
  plot(roc_list[[1]], col = colors[1], main = main, lwd = lwd)
  
  # Add additional ROC curves
  if (length(roc_list) > 1) {
    for (i in 2:length(roc_list)) {
      plot(roc_list[[i]], col = colors[i], add = TRUE, lwd = lwd)
    }
  }
  
  # Add legend
  legend("bottomright", legend = labels, col = colors, lwd = lwd)
}

#' Calculate sensitivity and specificity for a sequence of cutoffs
#' 
#' @param roc_obj ROC curve object
#' @param cutoffs Vector of cutoff values to evaluate
#' @return A dataframe with cutoffs, sensitivity, and specificity
calculate_cutoff_metrics <- function(roc_obj, cutoffs) {
  sensitivity_values <- numeric(length(cutoffs))
  specificity_values <- numeric(length(cutoffs))
  
  for(i in 1:length(cutoffs)) {
    # Get coordinates for this specific threshold
    coords_at_cutoff <- pROC::coords(roc_obj, cutoffs[i], input = "threshold")
    
    # Store values
    sensitivity_values[i] <- coords_at_cutoff$sensitivity
    specificity_values[i] <- coords_at_cutoff$specificity
  }
  
  # Create dataframe
  results_df <- data.frame(
    cutoff = round(cutoffs, 3),
    sensitivity = round(sensitivity_values, 3),
    specificity = round(specificity_values, 3)
  )
  
  return(results_df)
}

#' Find the optimal cutoff point using Youden's J statistic
#' 
#' @param roc_obj ROC curve object
#' @return The optimal cutoff value
find_best_cutoff <- function(roc_obj) {
  coords(roc_obj, "best", ret = "threshold")
}

#' Get detailed coordinates for all threshold points with Youden's index
#' 
#' @param roc_obj ROC curve object
#' @param num_points Number of points to include in the output (subsampled)
#' @return A dataframe with threshold, sensitivity, specificity, and Youden's index
get_detailed_coords <- function(roc_obj, num_points = 20) {
  # Get all coordinates
  coords_all <- pROC::coords(roc_obj, "all", ret = c("threshold", "sensitivity", "specificity"))
  
  # Convert to a dataframe
  coords_df <- as.data.frame(coords_all)
  
  # Calculate Youden's Index
  coords_df$youden_index <- coords_df$sensitivity + coords_df$specificity - 1
  
  # Round values
  coords_df <- coords_df %>%
    dplyr::mutate(
      threshold = round(threshold, 3),
      sensitivity = round(sensitivity, 3),
      specificity = round(specificity, 3),
      youden_index = round(youden_index, 3)
    )
  
  # Subsample points if requested
  if (num_points < nrow(coords_df)) {
    indices <- seq(1, nrow(coords_df), length.out = num_points)
    coords_df <- coords_df[indices, ]
  }
  
  return(coords_df)
}

#' Create a combined ROC model from multiple predictors
#' 
#' @param response The binary response variable
#' @param predictors A dataframe of predictor variables
#' @param model_formula Model formula (if NULL, all predictors are used)
#' @param family Model family (default is binomial for logistic regression)
#' @return A list with the model, ROC object, and AUC
create_combined_roc <- function(response, predictors, model_formula = NULL, family = binomial(link = "logit")) {
  # Create data for modeling
  model_data <- cbind(response = response, predictors)
  
  # Create formula if not provided
  if (is.null(model_formula)) {
    predictors_names <- colnames(predictors)
    formula_str <- paste("response ~", paste(predictors_names, collapse = " + "))
    model_formula <- as.formula(formula_str)
  }
  
  # Fit the model
  model <- glm(model_formula, family = family, data = model_data)
  
  # Get predicted probabilities
  predicted_probs <- predict(model, type = "response")
  
  # Create ROC object
  roc_obj <- pROC::roc(response, predicted_probs)
  auc_value <- pROC::auc(roc_obj)
  
  return(list(
    model = model,
    roc = roc_obj,
    auc = auc_value,
    auc_formatted = round(auc_value, 3),
    predicted_probs = predicted_probs
  ))
} 

#' Create Combined ROC plot
#' @param roc_obj The roc object 
#' @param title Title for plot
#' @param color color of ROC curve line 
#' @param legend_text Text for subtitle; automatically includes AUC number at end of line
create_roc_plot <- function(roc_obj, title, color = "blue", legend_text = "AUC") {
  # Convert ROC object to dataframe for ggplot
  roc_df <- data.frame(
    specificity = 1 - roc_obj$specificities,  # Convert to 1-specificity for x-axis
    sensitivity = roc_obj$sensitivities
  )
  
  # Calculate AUC for legend
  auc_value <- round(auc(roc_obj), 3)
  
  # Create full subtitle text with the custom legend text
  subtitle_text <- paste(legend_text, "=", auc_value)
  
  # Create plot
  plot <- ggplot(roc_df, aes(x = specificity, y = sensitivity)) +
    geom_line(color = color, size = 1) +
    geom_abline(slope = 1, intercept = 0, linetype = "dashed", color = "gray") +
    labs(
      title = title,
      x = "1-Specificity",
      y = "Sensitivity",
      subtitle = subtitle_text
    ) +
    theme_minimal() +
    coord_equal() +  # Ensure square aspect ratio
    theme(
      plot.title = element_text(size = 12, face = "bold"),
      plot.subtitle = element_text(size = 10)
    )
  
  return(plot)
}