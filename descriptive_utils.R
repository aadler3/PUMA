# Utility functions for descriptive statistics

#' Calculate and print summary statistics for a continuous variable by group
#' 
#' @param data The dataset
#' @param var_name The variable name
#' @param group_var The grouping variable name
#' @param na.rm Whether to remove NA values
#' @param round_digits Number of digits to round to
#' @return A data frame with summary statistics
summarize_by_group <- function(data, var_name, group_var = NULL, na.rm = TRUE, round_digits = 2) {
  if(is.null(group_var)) {
    # Overall summary
    result <- data.frame(
      Group = "Overall",
      N = sum(!is.na(data[[var_name]])),
      Missing = sum(is.na(data[[var_name]])),
      Mean = round(mean(data[[var_name]], na.rm = na.rm), round_digits),
      SD = round(sd(data[[var_name]], na.rm = na.rm), round_digits),
      Median = round(median(data[[var_name]], na.rm = na.rm), round_digits),
      Lower_Q = round(quantile(data[[var_name]], 0.25, na.rm = na.rm), round_digits),
      Upper_Q = round(quantile(data[[var_name]], 0.75, na.rm = na.rm), round_digits)
    )
    return(result)
  } else {
    # By group
    groups <- unique(data[[group_var]])
    results <- data.frame()
    
    # Overall summary
    overall <- summarize_by_group(data, var_name, NULL, na.rm, round_digits)
    results <- rbind(results, overall)
    
    # By group summary
    for(g in groups) {
      subset_data <- data[data[[group_var]] == g, ]
      group_result <- data.frame(
        Group = g,
        N = sum(!is.na(subset_data[[var_name]])),
        Missing = sum(is.na(subset_data[[var_name]])),
        Mean = round(mean(subset_data[[var_name]], na.rm = na.rm), round_digits),
        SD = round(sd(subset_data[[var_name]], na.rm = na.rm), round_digits),
        Median = round(median(subset_data[[var_name]], na.rm = na.rm), round_digits),
        Lower_Q = round(quantile(subset_data[[var_name]], 0.25, na.rm = na.rm), round_digits),
        Upper_Q = round(quantile(subset_data[[var_name]], 0.75, na.rm = na.rm), round_digits)
      )
      results <- rbind(results, group_result)
    }
    
    return(results)
  }
}

#' Run and print t-test results between two groups
#' 
#' @param data The dataset
#' @param var_name The variable name
#' @param group_var The grouping variable name
#' @return The t-test result
run_t_test <- function(data, var_name, group_var) {
  formula <- as.formula(paste(var_name, "~", group_var))
  test_result <- t.test(formula, data = data)
  return(test_result)
}

#' Calculate frequencies and percentages for a categorical variable by group
#' 
#' @param data The dataset
#' @param var_name The variable name
#' @param group_var The grouping variable name (optional)
#' @param print_results Whether to print the results
#' @return A list with count and percentage tables
categorical_summary <- function(data, var_name, group_var = NULL, print_results = TRUE) {
  if(is.null(group_var)) {
    # Overall counts
    count_table <- table(data[[var_name]])
    pct_table <- prop.table(count_table) * 100
    
    if(print_results) {
      cat("\nOverall:\n")
      print(count_table)
      cat("\nPercentages:\n")
      print(round(pct_table, 1))
    }
    
    return(list(counts = count_table, percentages = pct_table))
  } else {
    # By group
    groups <- unique(data[[group_var]])
    result_list <- list()
    
    # Overall summary
    overall <- categorical_summary(data, var_name, NULL, print_results)
    result_list[["overall"]] <- overall
    
    # By group summary
    for(g in groups) {
      if(print_results) cat("\nGroup:", g, "\n")
      subset_data <- data[data[[group_var]] == g, ]
      
      count_table <- table(subset_data[[var_name]])
      pct_table <- prop.table(count_table) * 100
      
      if(print_results) {
        print(count_table)
        cat("\nPercentages:\n")
        print(round(pct_table, 1))
      }
      
      result_list[[as.character(g)]] <- list(counts = count_table, percentages = pct_table)
    }
    
    return(result_list)
  }
}

#' Run chi-square test and print results
#' 
#' @param data The dataset
#' @param var1 First variable name
#' @param var2 Second variable name (optional for one-way tests)
#' @return The chi-square test result
run_chi_square <- function(data, var1, var2 = NULL) {
  if(is.null(var2)) {
    # One-way chi-square test
    observed <- table(data[[var1]])
    result <- chisq.test(observed)
  } else {
    # Two-way chi-square test
    observed <- table(data[[var1]], data[[var2]])
    result <- chisq.test(observed)
  }
  
  return(result)
}

#' Summarize variable across multiple visits
#' 
#' @param data The dataset in wide format
#' @param var_prefix The variable prefix (e.g., "ut_result" for "ut_result.1", "ut_result.2", etc.)
#' @param visits The visit numbers to process
#' @param group_var The grouping variable (optional)
#' @param stat_func The function to apply (e.g., mean, sum)
#' @param na.rm Whether to remove NA values
#' @param print_results Whether to print the results
#' @return A data frame with the results
summarize_across_visits <- function(data, var_prefix, visits = 1:5, 
                                   group_var = NULL, stat_func = mean, 
                                   na.rm = TRUE, print_results = TRUE) {
  results <- data.frame(Visit = visits)
  
  # Overall statistics
  overall_stats <- numeric(length(visits))
  overall_n <- numeric(length(visits))
  overall_missing <- numeric(length(visits))
  
  for(i in seq_along(visits)) {
    v <- visits[i]
    var_name <- paste0(var_prefix, ".", v)
    
    # Count non-missing and missing values
    overall_n[i] <- sum(!is.na(data[[var_name]]))
    overall_missing[i] <- sum(is.na(data[[var_name]]))
    
    # Calculate statistic
    overall_stats[i] <- stat_func(data[[var_name]], na.rm = na.rm)
  }
  
  results$Overall_Stat <- overall_stats
  results$Overall_N <- overall_n
  results$Overall_Missing <- overall_missing
  
  # By group statistics (if group_var is provided)
  if(!is.null(group_var)) {
    groups <- unique(data[[group_var]])
    
    for(g in groups) {
      subset_data <- data[data[[group_var]] == g, ]
      group_stats <- numeric(length(visits))
      group_n <- numeric(length(visits))
      group_missing <- numeric(length(visits))
      
      for(i in seq_along(visits)) {
        v <- visits[i]
        var_name <- paste0(var_prefix, ".", v)
        
        # Count non-missing and missing values
        group_n[i] <- sum(!is.na(subset_data[[var_name]]))
        group_missing[i] <- sum(is.na(subset_data[[var_name]]))
        
        # Calculate statistic
        group_stats[i] <- stat_func(subset_data[[var_name]], na.rm = na.rm)
      }
      
      results[[paste0(g, "_Stat")]] <- group_stats
      results[[paste0(g, "_N")]] <- group_n
      results[[paste0(g, "_Missing")]] <- group_missing
    }
  }
  
  if(print_results) {
    print(results)
  }
  
  return(results)
}

#' Run comparison tests across visits
#'
#' @param data The dataset in wide format
#' @param var_prefix The variable prefix
#' @param visits The visits to analyze
#' @param group_var The grouping variable
#' @param test_func The test function to use (e.g., t.test, chisq.test)
#' @param print_results Whether to print results
#' @return A data frame with test results
compare_across_visits <- function(data, var_prefix, visits = 1:5, 
                                 group_var, test_func = t.test, 
                                 print_results = TRUE) {
  results <- data.frame(Visit = visits)
  p_values <- numeric(length(visits))
  
  for(i in seq_along(visits)) {
    v <- visits[i]
    var_name <- paste0(var_prefix, ".", v)
    
    # Skip if all values are NA
    if(all(is.na(data[[var_name]]))) {
      p_values[i] <- NA
      next
    }
    
    # Run the test
    formula <- as.formula(paste(var_name, "~", group_var))
    test_result <- test_func(formula, data = data)
    
    # Extract p-value
    p_values[i] <- test_result$p.value
  }
  
  results$p_value <- p_values
  
  if(print_results) {
    print(results)
  }
  
  return(results)
} 

#' Analyze a categorical variable across multiple visits
#' 
#' @param data The dataset (should contain a 'visit' column)
#' @param var_names A character vector of variable names to analyze
#' @param visits The visit numbers to analyze
#' @param group_var The grouping variable name
#' @param print_header Whether to print section headers
#' @param filter_func Optional filtering function to apply before analysis
#' @return A list with results for each variable and visit
analyze_categorical_across_visits <- function(data, var_names, visits = 1:5, 
                                              group_var = "group_new", print_header = TRUE,
                                              filter_func = NULL) {
  results <- list()
  
  for (var_name in var_names) {
    if (print_header) {
      cat("\n----- ", var_name, " -----\n")
    }
    
    visit_results <- list()
    
    for (i in visits) {
      if (print_header) {
        cat("\n----- Visit", i, "-----\n")
      }
      
      # Filter data for current visit
      visit_data <- data[data$visit == i, ]
      
      # Apply optional filtering function if provided
      if (!is.null(filter_func)) {
        visit_data <- filter_func(visit_data)
      }
      
      # Run categorical summary
      result <- categorical_summary(
        data = visit_data,
        var_name = var_name,
        group_var = group_var,
        print_results = print_header
      )
      
      visit_results[[paste0("visit_", i)]] <- result
    }
    
    results[[var_name]] <- visit_results
  }
  
  return(results)
}