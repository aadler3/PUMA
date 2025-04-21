# PUMA Study Analysis

This repository contains code for analyzing data from the PUMA (Point-of-care Urine Monitoring of Adherence) study, which evaluates the impact of point-of-care urine tenofovir assays on PrEP adherence.

## Project Structure

The analysis code is organized into the following files:

### Utility Files

- `data_utils.R`: Functions for data loading and preprocessing
- `analysis_utils.R`: Functions for statistical analyses
- `viz_utils.R`: Functions for creating visualizations

### Analysis Files

- `PUMA_main.Rmd`: Main analysis document containing ALL original analyses
- `PUMA_supplementary.Rmd`: Additional complementary analyses
- `Puma Analyses FULL.Rmd`: The original R markdown file (kept for reference)

### Output Directory

- `figures/`: Directory for saving generated figures

## How to Use

1. Make sure all required packages are installed (see package loading sections in R markdown files)
2. Run the analyses in the following order:
   - First, run the main analysis: `PUMA_main.Rmd`
   - Then, run supplementary analyses: `PUMA_supplementary.Rmd`

## Complete Analysis Components

The main analysis file includes all original analyses:

1. **Descriptive Statistics**: Comprehensive baseline comparisons between study arms
2. **Month 12 Statistics**: Direct statistics from the final study visit
3. **BLQ Analysis**: Below Limit of Quantification analyses for all biomarkers
4. **Urine Test Results**: Complete analysis of urine test results across all visits
5. **Self-reported Adherence**: Analysis of all adherence measures
6. **Hair Analysis**: Complete analysis of all hair categorizations
7. **ROC Curve Analysis**: Receiver Operating Characteristic analyses for all biomarkers
8. **Combined Model**: Analysis combining multiple adherence measures

## Modifying the Analyses

The code is organized to make it easy to modify existing analyses while preserving ALL original variables and analyses:

- Functions in the utility files can be modified to adjust calculations
- The main file preserves all original analyses but in a more organized structure
- All variables from the original file are preserved and analyzed

## Adding New Analyses

To add a new analysis:

1. Add new analysis sections to the main or supplementary Rmd files
2. Use the utility functions to keep code concise and maintainable
3. Any additional utility functions can be added to the appropriate utility files

## Dataset

The analysis uses the following dataset:

- `PUMAextract_13Mar2025.csv`: Contains participant data including demographic information, adherence measures, and biomarker results

## Note on Completeness

Great care has been taken to ensure that NO variables or analyses from the original file are excluded. All analyses, including detailed explorations of specific variables and special cases (like individual participant analysis), have been preserved in the reorganized structure.

## Outputs

The analyses produce the following outputs:

- Tables of descriptive statistics
- Visualizations of adherence patterns over time
- Statistical test results comparing study arms
- ROC curve analyses for different adherence measures
- Model results for factors associated with adherence

## Contact

For questions about the analysis code, please contact the study team. 