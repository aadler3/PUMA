# Visualization utilities for PUMA study
# This file contains functions for creating visualizations of PUMA data

# Load required packages for visualizations
library(ggplot2)
library(flextable)
library(ggalluvial)
library(gridExtra)
library(kableExtra)

#' Create a line plot for longitudinal data by group
#'
#' @param data Data frame in long format
#' @param x_var X-axis variable (usually visit)
#' @param y_var Y-axis variable
#' @param group_var Grouping variable
#' @param title Plot title
#' @param y_label Y-axis label
#' @param x_label X-axis label
#' @param colors Vector of colors for groups
#' @param ylim Y-axis limits (optional)
#' @return ggplot object
#' @export
create_line_plot <- function(data, x_var, y_var, group_var, title = NULL, 
                           y_label = NULL, x_label = NULL, 
                           colors = c("darkorchid4", "skyblue4"),
                           ylim = NULL) {
  
  # Create plot
  p <- ggplot(data, aes_string(x = x_var, y = y_var, color = group_var, group = group_var)) +
    geom_line(size = 1) +
    geom_point(size = 3) +
    scale_color_manual(values = colors) +
    theme_minimal() +
    theme(legend.position = "top", legend.title = element_blank())
  
  # Add title if provided
  if (!is.null(title)) {
    p <- p + ggtitle(title)
  }
  
  # Add axis labels if provided
  if (!is.null(y_label)) {
    p <- p + ylab(y_label)
  }
  
  if (!is.null(x_label)) {
    p <- p + xlab(x_label)
  }
  
  # Add y-axis limits if provided
  if (!is.null(ylim)) {
    p <- p + coord_cartesian(ylim = ylim)
  }
  
  return(p)
}

#' Add mean/SD labels to a line plot
#'
#' @param plot ggplot object
#' @param data Data frame with means and SDs
#' @param x_var X-axis variable
#' @param y_var Y-axis variable (means)
#' @param label_var Variable containing formatted labels
#' @param offset Position offset for labels
#' @return Updated ggplot object
#' @export
add_mean_labels <- function(plot, data, x_var, y_var, label_var, offset = 2) {
  # Create a new variable for label positioning
  data$label_y <- ifelse(data[[y_var]] > mean(data[[y_var]]), 
                        data[[y_var]] - offset, 
                        data[[y_var]] + offset)
  
  # Add text labels
  plot + 
    geom_text(
      data = data,
      aes_string(x = x_var, y = "label_y", label = label_var),
      size = 3,
      vjust = "middle"
    )
}

#' Create an alluvial plot for categorical data across visits
#'
#' @param data Data frame in long format
#' @param id_var ID variable
#' @param var_name Variable to track
#' @param visit_var Visit variable
#' @param fill_var Variable for fill color
#' @param title Plot title
#' @return ggplot object
#' @export
create_alluvial_plot <- function(data, id_var, var_name, visit_var, fill_var, title = NULL) {
  # Create alluvial plot
  p <- ggplot(data,
         aes_string(x = visit_var, stratum = var_name, alluvium = id_var, fill = fill_var)) +
    geom_stratum() +
    geom_flow() +
    theme_minimal() +
    theme(legend.position = "bottom")
  
  # Add title if provided
  if (!is.null(title)) {
    p <- p + ggtitle(title)
  }
  
  return(p)
}

#' Create a forest plot for multiple analysis results
#'
#' @param data Data frame with effect estimates, CIs, and p-values
#' @param x_var Variable with effect estimates
#' @param lower_var Variable with lower CI bounds
#' @param upper_var Variable with upper CI bounds
#' @param label_var Variable with row labels
#' @param p_var Variable with p-values
#' @param null_value Reference/null value (usually 1 for OR, 0 for mean diff)
#' @param xlab X-axis label
#' @param title Plot title
#' @return ggplot object
#' @export
create_forest_plot <- function(data, x_var, lower_var, upper_var, label_var, p_var = NULL,
                             null_value = 0, xlab = "Estimate", title = NULL) {
  
  # Create basic forest plot
  p <- ggplot(data, aes_string(y = label_var, x = x_var, xmin = lower_var, xmax = upper_var)) +
    geom_pointrange() +
    geom_vline(xintercept = null_value, linetype = "dashed") +
    xlab(xlab) +
    ylab("") +
    theme_minimal()
  
  # Add p-values if provided
  if (!is.null(p_var)) {
    p <- p + 
      geom_text(aes_string(x = max(data[[upper_var]]) * 1.1, label = p_var), hjust = 0)
  }
  
  # Add title if provided
  if (!is.null(title)) {
    p <- p + ggtitle(title)
  }
  
  return(p)
}

#' Create a formatted table using flextable
#'
#' @param data Data frame to display
#' @param caption Table caption/title
#' @param digits Number of digits to round numeric columns
#' @return Flextable object
#' @export
create_flex_table <- function(data, caption = NULL, digits = 2) {
  # Round numeric columns
  data <- data %>%
    mutate(across(where(is.numeric), ~ round(., digits)))
  
  # Create flextable
  ft <- flextable(data)
  
  # Add caption if provided
  if (!is.null(caption)) {
    ft <- set_caption(ft, caption = caption)
  }
  
  # Auto-adjust column widths
  ft <- autofit(ft)
  
  return(ft)
}

#' Create boxplot comparing groups
#'
#' @param data Data frame
#' @param y_var Y-axis variable
#' @param group_var Grouping variable
#' @param title Plot title
#' @param y_label Y-axis label
#' @param colors Vector of colors for boxes
#' @return ggplot object
#' @export
create_boxplot <- function(data, y_var, group_var, title = NULL, 
                          y_label = NULL, colors = NULL) {
  
  # Create boxplot
  p <- ggplot(data, aes_string(x = group_var, y = y_var, fill = group_var)) +
    geom_boxplot() +
    theme_minimal() +
    theme(legend.position = "none")
  
  # Add custom colors if provided
  if (!is.null(colors)) {
    p <- p + scale_fill_manual(values = colors)
  }
  
  # Add title if provided
  if (!is.null(title)) {
    p <- p + ggtitle(title)
  }
  
  # Add Y-axis label if provided
  if (!is.null(y_label)) {
    p <- p + ylab(y_label)
  }
  
  return(p)
}

#' Create a bar plot for frequencies
#'
#' @param data Data frame
#' @param x_var X-axis variable
#' @param fill_var Fill variable (optional)
#' @param title Plot title
#' @param x_label X-axis label
#' @param y_label Y-axis label
#' @param position Position adjustment (default: "dodge")
#' @param colors Vector of colors
#' @return ggplot object
#' @export
create_bar_plot <- function(data, x_var, fill_var = NULL, title = NULL,
                          x_label = NULL, y_label = "Count",
                          position = "dodge", colors = NULL) {
  
  # Create mapping based on whether fill variable is provided
  if (is.null(fill_var)) {
    p <- ggplot(data, aes_string(x = x_var))
  } else {
    p <- ggplot(data, aes_string(x = x_var, fill = fill_var))
  }
  
  # Add geom_bar
  p <- p + 
    geom_bar(position = position) +
    theme_minimal()
  
  # Add custom colors if provided
  if (!is.null(colors)) {
    p <- p + scale_fill_manual(values = colors)
  }
  
  # Add title if provided
  if (!is.null(title)) {
    p <- p + ggtitle(title)
  }
  
  # Add axis labels if provided
  if (!is.null(x_label)) {
    p <- p + xlab(x_label)
  }
  
  if (!is.null(y_label)) {
    p <- p + ylab(y_label)
  }
  
  return(p)
} 