# PUMA Study Analysis

This repository contains code for analyzing data from the PUMA (Point-of-care Urine Monitoring of Adherence) study, which evaluates the impact of point-of-care urine tenofovir assays on PrEP adherence. 
This repository specifically focuses on an analysis looking at the association between different modes of monitoring PrEP adhernece.

## Project Structure

The analysis code is organized into the following files:

### Utility Files

- `data_utils.R`: Functions for data loading and preprocessing
- `analysis_utils.R`: Functions for statistical analyses
- `roc_utils.R`: Functions for creating ROC curves

### Analysis Files

- `PUMA Data Cleaning.Rmd`: All data cleaning required to run subsequent analyses
- `PUMA Descriptive Statistics.Rmd`: All descriptive statistics for analysis
- `PUMA Modeling Analysis.Rmd`: All modeling of primary and secondary analyses
- `PUMA Presenation Statistics.Rmd`: Specific analyses used for a presentation of the material
- `PUMA ROC Curves.Rmd`: ROC curve code
- `PUMA Visualizations.Rmd`: All visualizations excluding ROC curves

### Output Directory

- `Output/` : Directory for saving generated tables and figures

## How to Use

1. Make sure all required packages are installed (see package loading sections in R markdown files)
2. Run the analyses in the following order:
   - First, run the main analysis: `PUMA Data Cleaning.Rmd`
   - Then, all other analyses can be run in any order

## Dataset

The analysis uses an extracted CSV dataset that is not publically available. 


## Contact

For questions about the analysis code, please contact the study team: aliza.adler@ucsf.edu. 