# PUMA Analysis - Visualizations
# Author: PUMA Research Team
# Date: 2025

source("utils.R")

# Create adherence score line graph
create_adherence_line_plot <- function(line_graph_df, by_group = TRUE) {
  if (by_group) {
    # Group-specific line graph
    means_line_graph <- line_graph_df %>% 
      group_by(visit, group_new) %>% 
      summarise(
        mean_score = mean(wilson_adh, na.rm = TRUE),
        sd_score = sd(wilson_adh, na.rm = TRUE)
      ) %>%
      ungroup()
    
    # Create labels with mean and sd
    means_line_graph$label <- sprintf("%.1f\n(SD: %.1f)", 
                                    means_line_graph$mean_score, 
                                    means_line_graph$sd_score)
    
    # Add a small offset to position labels slightly above the points
    means_line_graph$label_y <- ifelse(
      means_line_graph$group_new == "Intervention",
      means_line_graph$mean_score + .5,
      ifelse(
        (means_line_graph$visit == 1 | means_line_graph$visit == 2 | means_line_graph$visit == 3) & 
          (means_line_graph$group_new == "Standard of Care"), 
        means_line_graph$mean_score - 2, 
        means_line_graph$mean_score + .5
      )
    )
    
    # Create the line plot
    plot <- ggplot(means_line_graph, 
                  aes(x = visit, y = mean_score, color = group_new, group = group_new)) +
      geom_line(size = 1) +
      geom_point(size = 3) +
      geom_text(
        aes(y = label_y, label = label),
        size = 3,
        vjust = "bottom"
      ) +
      scale_x_continuous(breaks = 1:5) +
      scale_color_manual(values = c("Intervention" = "darkorchid4", 
                                   "Standard of Care" = "skyblue4")) +
      labs(
        title = "Mean Adherence Score by Visit and Study Arm",
        x = "Visit",
        y = "Mean Wilson Adherence Score",
        color = "Study Arm"
      ) +
      theme_minimal() +
      theme(
        legend.position = "top",
        legend.title = element_blank()
      ) +
      coord_cartesian(ylim = c(70, 95))
  } else {
    # Overall line graph
    means_line_graph_t <- line_graph_df %>% 
      group_by(visit) %>% 
      summarise(
        mean_score = mean(wilson_adh, na.rm = TRUE),
        sd_score = sd(wilson_adh, na.rm = TRUE)
      ) %>%
      ungroup()
    
    # Create labels with mean and sd
    means_line_graph_t$label <- sprintf("%.1f\n(SD: %.1f)", 
                                      means_line_graph_t$mean_score, 
                                      means_line_graph_t$sd_score)
    
    # Add a small offset to position labels slightly above the points
    means_line_graph_t$label_y <- means_line_graph_t$mean_score + 2
    
    plot <- ggplot(means_line_graph_t, aes(x = visit, y = mean_score)) +
      geom_line(size = 1) +
      geom_point(size = 3) +
      geom_text(
        aes(y = label_y, label = label),
        size = 3,
        vjust = "bottom"
      ) +
      scale_x_continuous(breaks = 1:5) +
      labs(
        title = "Mean Adherence Score by Visit",
        x = "Visit",
        y = "Mean Wilson Adherence Score"
      ) +
      theme_minimal() +
      theme(
        legend.position = "top",
        legend.title = element_blank()
      ) +
      coord_cartesian(ylim = c(70, 95))
  }
  
  return(plot)
}

# Create stratified hair level plot
create_hair_level_plot <- function(data, visit_label) {
  aggregated_data <- data %>% 
    filter(!is.na(hair3code)) %>%
    group_by(hair3code) %>%
    summarise(
      mean_adh = mean(wilson_adh, na.rm = TRUE),
      sd_adh = sd(wilson_adh, na.rm = TRUE)
    )
  
  plot <- ggplot(aggregated_data, aes(x = hair3code, y = mean_adh)) +
    geom_bar(stat = "identity", fill = "slategray2") +
    labs(
      title = visit_label,
      x = "PrEP Doses per Week as Estimated by Hair TFV Levels",
      y = "Mean Wilson Adherence Score"
    ) +
    theme_minimal() +
    # Add mean values with SD in parentheses as text in the center of the bars
    geom_text(aes(label = sprintf("%.2f (%.2f)", mean_adh, sd_adh)), 
              vjust = -0.5,     
              color = "black") + 
    scale_y_continuous(expand = expansion(mult = c(0, 0.2)))
  
  return(plot)
}

# Create ROC curve plots
plot_roc_curves <- function(roc_urine, roc_hair) {
  # Get AUC values
  auc_urine <- auc(roc_urine)
  auc_hair <- auc(roc_hair)
  
  # Create the plot
  plot <- plot(roc_urine, col = "blue", 
             main = "ROC Curves for Adherence Prediction",
             xlab = "1 - Specificity", 
             ylab = "Sensitivity")
  
  # Add the second ROC curve
  plot <- plot(roc_hair, col = "red", add = TRUE)
  
  # Add legend with AUC values
  legend_text <- c(
    paste0("Urine Analysis (AUC = ", round(auc_urine, 3), ")"),
    paste0("Hair Analysis (AUC = ", round(auc_hair, 3), ")")
  )
  
  legend("bottomright", legend = legend_text, 
       col = c("blue", "red"), lwd = 2)
  
  return(plot)
}

# Create density plots for adherence by biomarker results
create_density_plots <- function(data, biomarker_type = "hair") {
  if (biomarker_type == "hair") {
    # First create a factor for visualization
    data <- data %>% 
      mutate(hair_uncoded = factor(
        ifelse(detect_tfv == 0, "Negative", "Positive"),
        levels = c("Negative", "Positive")
      ))
    
    # Filter to only include non-missing data
    viz_data <- data %>% filter(!is.na(detect_tfv))
    
    # Create the density plot
    plot <- ggplot(viz_data, 
                 aes(x = wilson_adh, fill = hair_uncoded, 
                     linetype = factor(group))) +
      geom_density(aes(color = hair_uncoded), alpha = 0.3) +
      scale_fill_manual(values = c("#3498db", "#e74c3c"), 
                       name = "PrEP Adherence (Hair)") +
      scale_color_manual(values = c("#3498db", "#e74c3c"), 
                        name = "PrEP Adherence (Hair)") +
      scale_linetype_manual(values = c("solid", "dashed"),
                           name = "Study Group") +
      labs(
        title = "Distribution of Wilson Adherence Scores",
        subtitle = "By PrEP adherence (hair) and study arm",
        x = "Wilson Adherence Score",
        y = "Density"
      ) +
      theme_minimal()
  } else {
    # First create a factor for visualization
    data <- data %>% 
      mutate(ut_uncoded = factor(
        ifelse(ut_result == 0, "Negative", "Positive"),
        levels = c("Negative", "Positive")
      ))
    
    # Filter to only include non-missing data
    viz_data <- data %>% filter(!is.na(ut_result))
    
    # Create the density plot
    plot <- ggplot(viz_data, 
                 aes(x = wilson_adh, fill = ut_uncoded, 
                     linetype = factor(group))) +
      geom_density(aes(color = ut_uncoded), alpha = 0.3) +
      scale_fill_manual(values = c("#3498db", "#e74c3c"), 
                       name = "Urine Assay Result") +
      scale_color_manual(values = c("#3498db", "#e74c3c"), 
                        name = "Urine Assay Result") +
      scale_linetype_manual(values = c("solid", "dashed"),
                           name = "Study Group") +
      labs(
        title = "Distribution of Wilson Adherence Scores",
        subtitle = "By urine assay result and study arm",
        x = "Wilson Adherence Score",
        y = "Density"
      ) +
      theme_minimal()
  }
  
  return(plot)
}

# Create histogram plots for adherence by biomarker results
create_histogram_plots <- function(data, biomarker_type = "hair") {
  if (biomarker_type == "hair") {
    # First create a factor for visualization
    data <- data %>% 
      mutate(hair_uncoded = factor(
        ifelse(detect_tfv == 0, "Negative", "Positive"),
        levels = c("Negative", "Positive")
      ))
    
    # Filter to only include non-missing data
    viz_data <- data %>% filter(!is.na(detect_tfv))
    
    # Create the histogram plot
    plot <- ggplot(viz_data, aes(x = wilson_adh, fill = hair_uncoded)) +
      geom_histogram(position = "identity", alpha = 0.6, bins = 10) +
      scale_fill_manual(values = c("#3498db", "#e74c3c"), 
                       name = "PrEP Adherence (Hair)") +
      facet_wrap(~group, ncol = 1) +
      labs(
        title = "Distribution of Wilson Adherence Scores",
        subtitle = "By PrEP adherence (hair) and study arm",
        x = "Wilson Adherence Score",
        y = "Count"
      ) +
      theme_minimal()
  } else {
    # First create a factor for visualization
    data <- data %>% 
      mutate(ut_uncoded = factor(
        ifelse(ut_result == 0, "Negative", "Positive"),
        levels = c("Negative", "Positive")
      ))
    
    # Filter to only include non-missing data
    viz_data <- data %>% filter(!is.na(ut_result))
    
    # Create the histogram plot
    plot <- ggplot(viz_data, aes(x = wilson_adh, fill = ut_uncoded)) +
      geom_histogram(position = "identity", alpha = 0.6, bins = 10) +
      scale_fill_manual(values = c("#3498db", "#e74c3c"), 
                       name = "Urine Assay Result") +
      facet_wrap(~group, ncol = 1) +
      labs(
        title = "Distribution of Wilson Adherence Scores",
        subtitle = "By urine assay result and study arm",
        x = "Wilson Adherence Score",
        y = "Count"
      ) +
      theme_minimal()
  }
  
  return(plot)
} 