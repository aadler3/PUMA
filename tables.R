# PUMA Analysis - Table Generation
# Author: PUMA Research Team
# Date: 2025

source("utils.R")

# Generate descriptive tables for hair levels
generate_hair_table <- function(puma_all) {
  results <- list()
  
  # Loop through visits 1-5 and generate tables by study arm and overall
  for (i in 1:5) {
    # Create data structure for this visit
    visit_data <- list()
    
    # Get all unique values of group_new for this visit
    groups <- unique(puma_all$group_new[puma_all$visit == i])
    
    # Loop through each group
    for (g in groups) {
      # Count table for this visit and group
      visit_group_table <- table(puma_all$hair3code[puma_all$visit == i & 
                                                  puma_all$group_new == g])
      
      # Percentage table for this visit and group
      visit_group_prop <- prop.table(visit_group_table) * 100
      
      # Store results
      visit_data[[paste0("group_", g)]] <- list(
        counts = visit_group_table,
        percentages = round(visit_group_prop, 1)
      )
    }
    
    # Show overall distribution for this visit
    visit_table <- table(puma_all$hair3code[puma_all$visit == i])
    visit_prop <- prop.table(visit_table) * 100
    
    visit_data[["overall"]] <- list(
      counts = visit_table,
      percentages = round(visit_prop, 1)
    )
    
    # Add to results
    results[[paste0("visit_", i)]] <- visit_data
  }
  
  return(results)
}

# Generate descriptive tables for hair levels with BLQ
generate_hair_blq_table <- function(puma_all) {
  results <- list()
  
  # Loop through visits 1-5 and generate tables by study arm and overall
  for (i in 1:5) {
    # Create data structure for this visit
    visit_data <- list()
    
    # Get all unique values of group_new for this visit
    groups <- unique(puma_all$group_new[puma_all$visit == i])
    
    # Loop through each group
    for (g in groups) {
      # Count table for this visit and group
      visit_group_table <- table(puma_all$hair4blq[puma_all$visit == i & 
                                                  puma_all$group_new == g])
      
      # Percentage table for this visit and group
      visit_group_prop <- prop.table(visit_group_table) * 100
      
      # Store results
      visit_data[[paste0("group_", g)]] <- list(
        counts = visit_group_table,
        percentages = round(visit_group_prop, 1)
      )
    }
    
    # Show overall distribution for this visit
    visit_table <- table(puma_all$hair4blq[puma_all$visit == i])
    visit_prop <- prop.table(visit_table) * 100
    
    visit_data[["overall"]] <- list(
      counts = visit_table,
      percentages = round(visit_prop, 1)
    )
    
    # Add to results
    results[[paste0("visit_", i)]] <- visit_data
  }
  
  return(results)
}

# Generate continuous hair summaries
generate_hair_continuous_summary <- function(puma_all) {
  results <- list()
  
  # Loop through visits 1-5 to calculate mean and sd of tfvngml
  for (i in 1:5) {
    visit_data <- list()
    
    # Calculate overall mean and sd for this visit
    overall_mean <- mean(puma_all$tfvngml[puma_all$visit == i], na.rm = TRUE)
    overall_sd <- sd(puma_all$tfvngml[puma_all$visit == i], na.rm = TRUE)
    
    visit_data[["overall"]] <- list(
      mean = round(overall_mean, 4),
      sd = round(overall_sd, 4)
    )
    
    # Get all unique values of group_new for this visit
    groups <- unique(puma_all$group_new[puma_all$visit == i])
    
    # Loop through each group
    for (g in groups) {
      # Calculate mean and sd for this visit and group
      group_mean <- mean(puma_all$tfvngml[puma_all$visit == i & 
                                        puma_all$group_new == g], 
                         na.rm = TRUE)
      
      group_sd <- sd(puma_all$tfvngml[puma_all$visit == i & 
                                     puma_all$group_new == g], 
                     na.rm = TRUE)
      
      visit_data[[paste0("group_", g)]] <- list(
        mean = round(group_mean, 4),
        sd = round(group_sd, 4)
      )
    }
    
    # Add to results
    results[[paste0("visit_", i)]] <- visit_data
  }
  
  return(results)
}

# Generate Table 1
generate_table1 <- function(table1_data) {
  # Variables for main demographic table
  table1_char <- table1_data %>% 
    select("age", "edu", "income", "sex1", "prep_length", 
           "padh3", "padh6", "emp_coded", 
           "ppoz_coded", "marital", "mh14", "bc_any", "travel_coded",
           "incsource_coded", "group")
  
  # Create compareGroups object for demographics by arm
  table1_arms <- compareGroups(group ~ ., data = table1_data, na.omit(TRUE))
  export_table1_arms <- createTable(table1_arms, show.p.overall = FALSE)
  
  # Extract sexual behavior variables
  table1_sex_data <- table1_data %>% 
    filter(sex2 > 0) %>% 
    select("sex5", "sex2", "condomless_sex",
           "new_partner", "HIV_pos_partner", 
           "unknown_HIV", "HIV_neg_partner", "sex14", "group")
  
  # Create compareGroups object for sexual behaviors by arm
  table1_sex <- compareGroups(group ~ ., data = table1_sex_data, 
                             na.omit(TRUE), method = 2)
  export_table1_sex <- createTable(table1_sex, show.p.overall = FALSE)
  
  # Overall tables without group comparison
  table1_all <- compareGroups(~ ., data = table1_char, na.omit(TRUE))
  export_table1_all <- createTable(table1_all, show.p.overall = FALSE)
  
  table1_sex_all <- compareGroups(~ ., data = table1_sex_data, 
                                na.omit(TRUE), method = 2)
  export_table1_sex_all <- createTable(table1_sex_all, show.p.overall = FALSE)
  
  # Return all tables
  return(list(
    demographics_by_arm = export_table1_arms,
    sexual_behavior_by_arm = export_table1_sex,
    demographics_overall = export_table1_all,
    sexual_behavior_overall = export_table1_sex_all
  ))
}

# Generate adherence summary tables
generate_adherence_summary <- function(adherence_wide) {
  results <- list()
  
  # Wilson adherence score summary
  # List of variables to process
  wilson_vars <- paste0("wilson_adh.", 1:5)
  
  # List of datasets to process
  datasets <- list(
    all = adherence_wide,
    intervention = adherence_wide %>% filter(group_new == "Intervention"),
    soc = adherence_wide %>% filter(group_new == "Standard of Care")
  )
  
  dataset_names <- c("all", "intervention", "soc")
  
  # Create a dataframe to store results
  wilson_results <- data.frame(
    Dataset = character(),
    Variable = character(),
    NonNA_Count = numeric(),
    NA_Count = numeric(),
    Mean = numeric(),
    SD = numeric(),
    stringsAsFactors = FALSE
  )
  
  # Loop through each dataset
  for (d in 1:length(datasets)) {
    dataset <- datasets[[d]]
    dataset_name <- dataset_names[d]
    
    # Loop through each variable
    for (var in wilson_vars) {
      # Extract the variable
      values <- dataset[[var]]
      
      # Calculate metrics
      nonNA_count <- sum(!is.na(values))
      na_count <- sum(is.na(values))
      mean_val <- mean(values, na.rm = TRUE)
      sd_val <- sd(values, na.rm = TRUE)
      
      # Add results to dataframe
      wilson_results <- rbind(wilson_results, data.frame(
        Dataset = dataset_name,
        Variable = var,
        NonNA_Count = nonNA_count,
        NA_Count = na_count,
        Mean = round(mean_val, 2),
        SD = round(sd_val, 2)
      ))
    }
  }
  
  results[["wilson_adh"]] <- wilson_results
  
  # Repeat for padh7 scores
  padh7_vars <- paste0("padh7_recode.", 1:5)
  
  padh7_results <- data.frame(
    Dataset = character(),
    Variable = character(),
    NonNA_Count = numeric(),
    NA_Count = numeric(),
    Mean = numeric(),
    SD = numeric(),
    stringsAsFactors = FALSE
  )
  
  # Loop through each dataset
  for (d in 1:length(datasets)) {
    dataset <- datasets[[d]]
    dataset_name <- dataset_names[d]
    
    # Loop through each variable
    for (var in padh7_vars) {
      # Extract the variable
      values <- dataset[[var]]
      
      # Calculate metrics
      nonNA_count <- sum(!is.na(values))
      na_count <- sum(is.na(values))
      mean_val <- mean(values, na.rm = TRUE)
      sd_val <- sd(values, na.rm = TRUE)
      
      # Add results to dataframe
      padh7_results <- rbind(padh7_results, data.frame(
        Dataset = dataset_name,
        Variable = var,
        NonNA_Count = nonNA_count,
        NA_Count = na_count,
        Mean = round(mean_val, 2),
        SD = round(sd_val, 2)
      ))
    }
  }
  
  results[["padh7"]] <- padh7_results
  
  # Repeat for padh4 scores
  padh4_vars <- paste0("padh4_recode.", 1:5)
  
  padh4_results <- data.frame(
    Dataset = character(),
    Variable = character(),
    NonNA_Count = numeric(),
    NA_Count = numeric(),
    Mean = numeric(),
    SD = numeric(),
    stringsAsFactors = FALSE
  )
  
  # Loop through each dataset
  for (d in 1:length(datasets)) {
    dataset <- datasets[[d]]
    dataset_name <- dataset_names[d]
    
    # Loop through each variable
    for (var in padh4_vars) {
      # Extract the variable
      values <- dataset[[var]]
      
      # Calculate metrics
      nonNA_count <- sum(!is.na(values))
      na_count <- sum(is.na(values))
      mean_val <- mean(values, na.rm = TRUE)
      sd_val <- sd(values, na.rm = TRUE)
      
      # Add results to dataframe
      padh4_results <- rbind(padh4_results, data.frame(
        Dataset = dataset_name,
        Variable = var,
        NonNA_Count = nonNA_count,
        NA_Count = na_count,
        Mean = round(mean_val, 2),
        SD = round(sd_val, 2)
      ))
    }
  }
  
  results[["padh4"]] <- padh4_results
  
  # Add t-test results
  # Wilson adherence
  t_test_results <- list()
  for (i in 1:5) {
    t_test_results[[paste0("wilson_visit", i)]] <- 
      t.test(as.formula(paste0("wilson_adh.", i, " ~ group_new")), 
             data = adherence_wide)
  }
  
  # padh7
  for (i in 1:5) {
    t_test_results[[paste0("padh7_visit", i)]] <- 
      t.test(as.formula(paste0("padh7_recode.", i, " ~ group_new")), 
             data = adherence_wide)
  }
  
  # padh4
  for (i in 1:5) {
    t_test_results[[paste0("padh4_visit", i)]] <- 
      t.test(as.formula(paste0("padh4_recode.", i, " ~ group_new")), 
             data = adherence_wide)
  }
  
  results[["t_tests"]] <- t_test_results
  
  return(results)
}

# Generate study retention table
generate_retention_table <- function(primary_wide) {
  results <- list()
  
  # Variables to track
  ut_vars <- c("ut_result.1", "ut_result.2", "ut_result.3", "ut_result.4", "ut_result.5")
  
  # Total Sample
  total_sample <- sapply(ut_vars, function(var) {
    sum(!is.na(primary_wide[[var]]))
  })
  
  results[["total"]] <- total_sample
  
  # Intervention Arm
  poc <- primary_wide %>% filter(group_new == "Intervention")
  
  poc_counts <- sapply(ut_vars, function(var) {
    sum(!is.na(poc[[var]]))
  })
  
  poc_percentages <- (poc_counts / length(poc$ptid)) * 100
  
  results[["intervention"]] <- list(
    counts = poc_counts,
    percentages = round(poc_percentages, 1)
  )
  
  # Control Arm
  soc <- primary_wide %>% filter(group_new == "Standard of Care")
  
  soc_counts <- sapply(ut_vars, function(var) {
    sum(!is.na(soc[[var]]))
  })
  
  soc_percentages <- (soc_counts / length(soc$ptid)) * 100
  
  results[["control"]] <- list(
    counts = soc_counts,
    percentages = round(soc_percentages, 1)
  )
  
  # For hair and urine both available
  hair_ut_available <- list()
  
  # Non-cumulative retention
  month1 <- primary_wide %>% filter((!is.na(ut_result.1)) & (!is.na(detect_tfv.1)))
  month3 <- primary_wide %>% filter((!is.na(ut_result.2)) & (!is.na(detect_tfv.2)))
  month6 <- primary_wide %>% filter((!is.na(ut_result.3)) & (!is.na(detect_tfv.3)))
  month9 <- primary_wide %>% filter((!is.na(ut_result.4)) & (!is.na(detect_tfv.4)))
  month12 <- primary_wide %>% filter((!is.na(ut_result.5)) & (!is.na(detect_tfv.5)))
  
  hair_ut_available[["non_cumulative"]] <- list(
    month1 = table(month1$group_new),
    month3 = table(month3$group_new),
    month6 = table(month6$group_new),
    month9 = table(month9$group_new),
    month12 = table(month12$group_new)
  )
  
  # Cumulative retention
  month3_c <- month1 %>% filter((!is.na(ut_result.2)) & (!is.na(detect_tfv.2)))
  month6_c <- month3_c %>% filter((!is.na(ut_result.3)) & (!is.na(detect_tfv.3)))
  month9_c <- month6_c %>% filter((!is.na(ut_result.4)) & (!is.na(detect_tfv.4)))
  month12_c <- month9_c %>% filter((!is.na(ut_result.5)) & (!is.na(detect_tfv.5)))
  
  hair_ut_available[["cumulative"]] <- list(
    month3 = table(month3_c$group_new),
    month6 = table(month6_c$group_new),
    month9 = table(month9_c$group_new),
    month12 = table(month12_c$group_new)
  )
  
  results[["hair_ut_available"]] <- hair_ut_available
  
  return(results)
} 