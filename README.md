# PUMA Analysis

This repository contains code for analyzing data from the PUMA (Point-of-care Urine Monitoring of Adherence) study.

## Project Structure

The project has been organized into modular R files for better maintainability and readability:

1. **Puma Analyses.Rmd**: The main analysis file that imports functions from the supporting files and presents the results.

2. **utils.R**: Utility functions including package loading, data path management, and common helper functions.

3. **data_processing.R**: Functions for data processing, variable creation, and dataset preparation.

4. **models.R**: Statistical modeling functions for GEE models, ROC analysis, and model selection.

5. **visualizations.R**: Functions for creating plots and visualizations.

6. **tables.R**: Functions for generating descriptive tables and summaries.

## Getting Started

To run the analysis:

1. Open `Puma Analyses.Rmd` in RStudio
2. Make sure all required packages are installed
3. Ensure the data file is in the correct location
4. Knit the document to produce the analysis report

## Required Packages

The analysis requires the following R packages:
- dplyr
- flextable
- ggplot2
- collapse
- labelled
- compareGroups
- gee
- geepack
- pROC
- ROCR
- plotROC
- gridExtra

These are automatically loaded when sourcing `utils.R`.

## Data

The analysis expects a data file named `PUMAextract_13Mar2025.csv` in the data directory.

## Analysis Overview

The analysis includes:
- Data processing and variable creation
- Descriptive statistics
- Statistical modeling with GEE
- ROC curve analysis
- Visualizations of adherence patterns 