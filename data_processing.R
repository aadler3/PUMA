# PUMA Analysis - Data Processing
# Author: PUMA Research Team
# Date: 2025

source("utils.R")

# Process and create new variables in the PUMA dataset
process_puma_data <- function(puma_all) {
  # Create any BC variable
  puma_all$bc_any <- ifelse(puma_all$par5_recode == 0, 0, 
                           ifelse(puma_all$par5_recode == 99, NA, 1))
  
  # Create partner variables
  puma_all <- puma_all %>% mutate(
    new_partner = ifelse(sex6 >= 1, 1, 0),
    HIV_pos_partner = ifelse(sex7 >= 1, 1, 0),
    unknown_HIV = ifelse(sex9 >= 1, 1, 0),
    HIV_neg_partner = ifelse(sex11 >= 1, 1, 0)
  )
  
  # Create condomless sex variable
  puma_all <- puma_all %>% 
    mutate(condomless_sex = ifelse(sex2 == sex3, 0, 1)) %>%
    relocate(condomless_sex, .after = sex3)
  
  # Recode demographic variables
  puma_all <- puma_all %>% mutate(
    emp_coded = case_when(
      emp == 1 ~ "Laborer/semi-skilled",
      emp == 2 ~ "Trade/sales",
      emp == 3 ~ "Student",
      emp == 4 ~ "Professional",
      emp == 5 ~ "Farming/animal raising",
      emp == 6 ~ "Housewife",
      emp == 99 ~ "Multiple/Other"
    ),
    incsource_coded = case_when(
      incsource == 1 ~ "My work",
      incsource == 2 ~ "My husband's work",
      incsource == 3 ~ "My family",
      incsource == 9 ~ "Other"
    ),
    ppoz_coded = case_when(
      ppoz == 0 ~ "No",
      ppoz == 1 ~ "Yes",
      ppoz == 2 ~ "Don't know"
    ),
    travel_coded = case_when(
      travel == 1 ~ "Less than 30 minutes",
      travel == 2 ~ "30-59 minutes",
      travel == 3 ~ "1-2 hours",
      travel == 4 ~ "More than 2 hours"
    )
  )
  
  # Process group variable
  val_labels(puma_all$group) <- c("Intervention" = 1, "Standard of Care" = 2)
  puma_all$group <- to_character(puma_all$group)
  
  puma_all$group_new <- case_when(
    puma_all$group == "Standard of Care" ~ 1,
    puma_all$group == "Intervention" ~ 2
  )
  
  puma_all <- puma_all %>% relocate(group_new, .after = group)
  
  val_labels(puma_all$group_new) <- c("Intervention" = 2, "Standard of Care" = 1)
  puma_all$group_new <- to_character(puma_all$group_new)
  puma_all$group_new <- factor(puma_all$group_new, 
                             levels = c("Standard of Care", "Intervention"))
  
  # Process adherence variables
  val_labels(puma_all$padh7) <- c("Very poor" = 0, "Poor" = 1, "Fair" = 2, 
                                 "Good" = 3, "Very good" = 4, "Excellent" = 5)
  puma_all$padh7 <- to_character(puma_all$padh7)
  
  puma_all <- puma_all %>% 
    mutate(padh7_recode = case_when(
      padh7 == "Very poor" ~ 0,
      padh7 == "Poor" ~ .2,
      padh7 == "Fair" ~ .4,
      padh7 == "Good" ~ .6,
      padh7 == "Very good" ~ .8,
      padh7 == "Excellent" ~ 1
    )) %>%
    relocate(padh7_recode, .after = padh7)
  
  puma_all <- puma_all %>% 
    mutate(padh4 = ifelse(padh3 == 1, 0, padh4))
  
  puma_all <- puma_all %>% 
    mutate(padh4_recode = (30-padh4)/30) %>%
    relocate(padh4_recode, .after = padh4)
  
  puma_all <- puma_all %>% 
    mutate(avg_padh = ((padh4_recode + padh7_recode)/2)*100) %>%
    relocate(avg_padh, .after = wilson_adh)
  
  puma_all <- puma_all %>% 
    mutate(wilson_adh_scaled_sd = wilson_adh/sd(wilson_adh, na.rm = TRUE)) %>%
    relocate(wilson_adh_scaled_sd, .after = wilson_adh)
  
  # Create visit factor variable
  puma_all$visit_coded <- factor(puma_all$visit, 
                               levels = c(1, 2, 3, 4, 5),
                               labels = c("Baseline", "Month3", "Month6",
                                          "Month9", "Month12"))
  
  # Convert ID to character
  puma_all$ptid <- as.character(puma_all$ptid)
  
  # Process hair data measurements
  
  # First keep a numeric copy
  puma_all$hair3numeric <- puma_all$hair3code
  
  # 3 levels 
  puma_all$hair3code <- factor(puma_all$hair3code, 
                              levels = c(1, 2, 3),
                              labels = c("<2", "2 to 3", ">=4"))
  
  # 4 levels
  puma_all$hair4code <- factor(puma_all$hair4code, 
                              levels = c(1, 2, 3, 4),
                              labels = c("<2", "2 to 3", "4-6", "Daily"))
  
  # 3 with BLQ separated out
  puma_all$hair4blq <- factor(puma_all$hair4blq, 
                             levels = c(0, 1, 2, 3),
                             labels = c("BLQ", "<2", "2 to 3", ">=4"))
  
  # 4 with BLQ separated out
  puma_all$hair5blq <- factor(puma_all$hair5blq, 
                             levels = c(0, 1, 2, 3, 4),
                             labels = c("BLQ", "<2", "2 to 3", "4-6", "Daily"))
  
  # Add a wilson_adh10 variable (scaled by 10)
  puma_all$wilson_adh10 <- puma_all$wilson_adh/10
  
  return(puma_all)
}

# Create datasets for specific time points or groups
create_subset_data <- function(puma_all) {
  subsets <- list()
  
  # Baseline data for Table 1
  subsets$baseline <- puma_all %>% filter(visit == 1)
  
  # Month 12 data
  subsets$month12 <- puma_all %>% filter(visit == 5)
  
  # Intervention and SOC separately
  subsets$poc_all <- puma_all %>% filter(group == "Intervention")
  subsets$soc_all <- puma_all %>% filter(group == "Standard of Care")
  
  return(subsets)
}

# Create binary variables coded as "Yes"/"No" for table 1
create_table1_data <- function(baseline_data) {
  bin_variables <- c(
    "marital", "sex8", "sex10", "sex12", "sex14", "par5___1", "par5___2",
    "par5___3", "par5___4", "par5___5", "par5___6", "par5___7", 
    "par5___8", "par5___9", "padh3", "padh9___1", "padh9___2", 
    "padh9___3", "padh9___4", "padh9___5", "padh9___6", "padh9___7",
    "padh9___8", "padh9___9", "padh9___10", "padh9___11", "padh9___12",
    "padh9___13", "padh9___14", "padh9___15", "padh9____9", "padh10",
    "mh14", "mh16", "mh17", "ut_result", "detect_tfv", "condomless_sex",
    "new_partner", "HIV_pos_partner", "unknown_HIV", "HIV_neg_partner", 
    "bc_any"
  )
  
  table1_data <- baseline_data
  
  for (var in bin_variables) {
    table1_data[[var]] <- ifelse(table1_data[[var]] == 1, "Yes", "No")
  }
  
  return(table1_data)
}

# Create long and wide datasets for specific analyses
create_analysis_datasets <- function(puma_all) {
  datasets <- list()
  
  # Hair analysis datasets
  hair_long <- puma_all %>% select(group, detect_tfv, visit, ptid)
  hair_long <- as.data.frame(hair_long)
  
  hair_wide <- reshape(hair_long, v.names = "detect_tfv", 
                      idvar = c("ptid", "group"), 
                      timevar = "visit", direction = "wide")
  
  datasets$hair_long <- hair_long
  datasets$hair_wide <- hair_wide
  
  # Urine analysis datasets
  ut_long <- puma_all %>% select(group_new, ut_result, visit, ptid)
  ut_long <- as.data.frame(ut_long)
  
  ut_wide <- reshape(ut_long, v.names = "ut_result", 
                    idvar = c("ptid", "group_new"), 
                    timevar = "visit", direction = "wide")
  
  datasets$ut_long <- ut_long
  datasets$ut_wide <- ut_wide
  
  # Self-reported adherence datasets
  adherence <- puma_all %>% 
    select(group_new, wilson_adh, padh7_recode, padh4_recode, visit, ptid)
  
  adherence_wide <- reshape(adherence, 
                          v.names = c("wilson_adh", "padh7_recode", "padh4_recode"), 
                          idvar = "ptid", timevar = "visit", direction = "wide")
  
  datasets$adherence <- adherence
  datasets$adherence_wide <- adherence_wide
  
  # Primary analysis dataset in wide format
  primary_long <- puma_all %>% 
    select(ptid, group_new, visit, detect_tfv, ut_result,
           wilson_adh, wilson_adh_scaled_sd)
  
  primary_wide <- reshape(primary_long, 
                        idvar = c("ptid", "group_new"), 
                        timevar = "visit", direction = "wide")
  
  datasets$primary_long <- primary_long
  datasets$primary_wide <- primary_wide
  
  # GEE analysis datasets
  hair_gee_data <- puma_all %>%
    filter(complete.cases(detect_tfv, wilson_adh, group_new, ptid)) %>% 
    select(detect_tfv, wilson_adh_scaled_sd, wilson_adh10, 
           group_new, ptid, wilson_adh, mh14, hair3numeric)
  
  ut_gee_data <- puma_all %>%
    filter(complete.cases(ut_result, wilson_adh, group_new, ptid)) %>% 
    select(ut_result, wilson_adh, wilson_adh_scaled_sd, 
           wilson_adh10, group_new, ptid, mh14)
  
  datasets$hair_gee_data <- hair_gee_data
  datasets$ut_gee_data <- ut_gee_data
  
  return(datasets)
} 