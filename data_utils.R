# Data loading and preprocessing utilities for PUMA Analysis
# This file contains functions for loading and preprocessing data for the PUMA study

# Load required packages for data preprocessing
library(dplyr)
library(labelled)
library(collapse)

#' Load PUMA dataset
#' 
#' @param file_path Path to the CSV file or NULL to use the default path
#' @return Loaded dataframe
#' @export
load_puma_data <- function(file_path = NULL) {
  if (is.null(file_path)) {
    pathway <- "~/Library/CloudStorage/Box-Box/Puma Analysis/Data/"
    file_path <- file.path(pathway, "PUMAextract_13Mar2025.csv")
  }
  
  puma_all <- read.csv(file_path)
  return(puma_all)
}

#' Preprocess PUMA data with all recoding and variable creation
#' 
#' @param puma_all Raw PUMA dataset
#' @return Processed dataframe
#' @export
preprocess_puma_data <- function(puma_all) {
  # Create any BC variable
  puma_all$bc_any <- ifelse(puma_all$par5_recode == 0, 0, 
                           ifelse(puma_all$par5_recode == 99, NA, 1))
  
  # Create partner variables
  puma_all <- puma_all %>% mutate(
    new_partner = ifelse(sex6 >=1, 1, 0),
    HIV_pos_partner = ifelse(sex7 >=1, 1, 0),
    unknown_HIV = ifelse(sex9 >= 1, 1, 0),
    HIV_neg_partner = ifelse(sex11 >=1, 1, 0)
  )
  
  # Create condomless sex variable
  puma_all <- puma_all %>% 
    mutate(condomless_sex = ifelse(sex2 == sex3, 0, 1)) %>%
    relocate(condomless_sex, .after = sex3)
  
  # Recode categorical variables
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
  
  # Recode binary variables to Yes/No
  binary_vars <- c("marital", "sex8", "sex10", "sex12", "sex14", "par5___1", "par5___2",
                                    "par5___3", "par5___4", "par5___5", "par5___6", "par5___7", 
                                    "par5___8", "par5___9", "padh3", "padh9___1", "padh9___2", 
                                    "padh9___3", "padh9___4", "padh9___5", "padh9___6", "padh9___7",
                                    "padh9___8", "padh9___9", "padh9___10", "padh9___11", "padh9___12",
                                    "padh9___13", "padh9___14", "padh9___15", "padh9____9", "padh10",
                                    "mh14", "mh16", "mh17", "ut_result", "detect_tfv", "condomless_sex",
                                    "new_partner", "HIV_pos_partner","unknown_HIV", "HIV_neg_partner", 
                                    "bc_any")
  for (var in binary_vars) {
    if (var %in% names(puma_all)) {
      # Create a new variable with the _yn suffix
      new_var_name <- paste0(var, "_yn")
      puma_all[[new_var_name]] <- ifelse(is.na(puma_all[[var]]), NA, 
                                       ifelse(puma_all[[var]] == 1, "Yes", "No"))
      
      # Place the new variable after the original
      puma_all <- puma_all %>% relocate(all_of(new_var_name), .after = all_of(var))
    }
  }
  
  # Recode study group
  val_labels(puma_all$group) <- c("Intervention" = 1, "Standard of Care" = 2)
  puma_all$group <- to_character(puma_all$group)
  
  puma_all$group_new <- case_when(
    puma_all$group == "Standard of Care" ~ 1,
    puma_all$group == "Intervention" ~ 2
  )
  
  puma_all <- puma_all %>% relocate(group_new, .after = group)
  
  val_labels(puma_all$group_new) <- c("Intervention" = 2, "Standard of Care" = 1)
  puma_all$group_new <- to_character(puma_all$group_new)
  puma_all$group_new <- factor(puma_all$group_new, levels = c("Standard of Care", "Intervention"))
  
  # Adherence variables
  val_labels(puma_all$padh7) <- c("Very poor" = 0, "Poor" = 1, "Fair" = 2, 
                                 "Good" = 3, "Very good" = 4, "Excellent" = 5)
  puma_all$padh7 <- to_character(puma_all$padh7)
  
  puma_all <- puma_all %>% mutate(
    padh7_recode = case_when(
      padh7 == "Very poor" ~ 0,
      padh7 == "Poor" ~ .2,
      padh7 == "Fair" ~ .4,
      padh7 == "Good" ~ .6,
      padh7 == "Very good" ~ .8,
      padh7 == "Excellent" ~ 1
    )
  ) %>% relocate(padh7_recode, .after = padh7)
  
  puma_all <- puma_all %>% 
    mutate(padh4 = ifelse(padh3 == 1, 0, padh4)) %>%
    mutate(padh4_recode = (30-padh4)/30) %>%
    relocate(padh4_recode, .after = padh4) %>%
    mutate(avg_padh = ((padh4_recode + padh7_recode)/2)*100) %>%
    relocate(avg_padh, .after = wilson_adh) %>%
    mutate(wilson_adh_scaled_sd = wilson_adh/sd(wilson_adh, na.rm = TRUE)) %>%
    relocate(wilson_adh_scaled_sd, .after = wilson_adh)
  
  # Visit coding
  puma_all$visit_coded <- factor(
    puma_all$visit, 
    levels = c(1, 2, 3, 4, 5),
    labels = c("Baseline", "Month3", "Month6", "Month9", "Month12")
  )
  
  puma_all$ptid <- as.character(puma_all$ptid)
  
  # Hair data recoding
  puma_all$hair3numeric <- puma_all$hair3code
  
  puma_all$hair3code <- factor(
    puma_all$hair3code, 
    levels = c(1, 2, 3),
    labels = c("<2", "2 to 3", ">=4")
  )
  
  puma_all$hair4code <- factor(
    puma_all$hair4code, 
    levels = c(1, 2, 3, 4),
    labels = c("<2", "2 to 3", "4-6", "Daily")
  )
  
  puma_all$hair4blq <- factor(
    puma_all$hair4blq, 
    levels = c(0, 1, 2, 3),
    labels = c("BLQ", "<2", "2 to 3", ">=4")
  )
  
  puma_all$hair5blq <- factor(
    puma_all$hair5blq, 
    levels = c(0, 1, 2, 3, 4),
    labels = c("BLQ", "<2", "2 to 3", "4-6", "Daily")
  )
  
  # Low vs High Hair Variable
  puma_all <- puma_all %>% 
    mutate(
      hair2code = case_when(
        hair3code == "<2" ~ "<=3",
        hair3code == "2 to 3" ~ "<=3",
        hair3code == ">=4" ~ ">=4"
      ),
      hair2_coded = ifelse(hair2code == ">=4", 1, 0)
    ) %>%
    relocate(c(hair2code, hair2_coded), .after = hair3code)
  
  # Create categorical Wilson Adherence Score for Alluvial Plots
  puma_all <- puma_all %>% 
    mutate(
      wilson_cat = case_when(
        wilson_adh >= 0 & wilson_adh <= 20 ~ "0-20",
        wilson_adh >= 21 & wilson_adh <= 40 ~ "21-40",
        wilson_adh >= 41 & wilson_adh <= 60 ~ "41-60",
        wilson_adh >= 61 & wilson_adh <= 80 ~ "61-80",
        wilson_adh >= 81 & wilson_adh <= 100 ~ "81-100"
      )
    ) %>%
    relocate(wilson_cat, .after = wilson_adh)
  
  return(puma_all)
}

#' Create a dataframe containing only data from a specific visit 
#'
#' @param data Processed PUMA dataset
#' @param visit_num Visit number (1-5)
#' @return Filtered dataset
#' @export
get_visit_data <- function(data, visit_num) {
  if (!visit_num %in% 1:5) {
    stop("Visit number must be between 1 and 5")
  }
  
  visit_data <- data %>% filter(visit == visit_num)
  return(visit_data)
}

#' Reshape data from long to wide format for analyses
#'
#' @param data Processed PUMA dataset
#' @param vars Variables to reshape
#' @return Wide format dataset
#' @export
reshape_wide <- function(data, vars) {
  data_selected <- data %>% select(c("group_new", "ptid", "visit", vars))
  
  # Convert to wide format
  data_wide <- reshape(
    data_selected, 
    v.names = vars, 
    idvar = c("ptid", "group_new"),
    timevar = "visit", 
    direction = "wide"
  )
  
  return(data_wide)
} 