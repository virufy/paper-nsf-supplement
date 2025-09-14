================================================================================
Replication Files for:
"The Double Pivot: A Mixed-Methods Case Study of the NSF I-Corps Customer Discovery Process in a University-Spun Health AI Venture"
================================================================================

Authors: Amil Khanzada, Takuji Takemoto
Contact: Amil Khanzada (kad23802@u-fukui.ac.jp)
Date: September 14th, 2025

This archive contains the necessary data and code to reproduce the quantitative analysis and figures presented in the manuscript.

--------------------------------------------------------------------------------
CONTENTS
--------------------------------------------------------------------------------

1. icorps_interviews_coded.csv  - The anonymized, coded dataset.
   - Rows: 148 interviews
   - Columns:
     - InterviewID: A unique identifier for each interview.
     - Phase: The program phase (1=Early, 2=Middle, 3=Late).
     - Theme Columns: Eight columns representing the coded themes (1=present, 0=absent).

2. icorps_reanalysis.R
   - The R script for analysis and figure generation.
   - This script takes the CSV file as input and generates all statistical tables and figures reported in the paper.

--------------------------------------------------------------------------------
INSTRUCTIONS FOR REPRODUCTION
--------------------------------------------------------------------------------

To reproduce the analysis, please ensure you have R installed along with the following packages:
- data.table
- MASS
- ggplot2
- tidyr
- scales
- stringr

You can install these packages by running the following command in R:
install.packages(c("data.table", "MASS", "ggplot2", "tidyr", "scales", "stringr"))

Once the packages are installed, place both the `icorps_reanalysis.R` script and the `icorps_interviews_coded.csv` file in the same directory.

Navigate to this directory in your terminal and run the following command:

   Rscript icorps_reanalysis.R icorps_interviews_coded.csv --plot --png

This command will:
1. Run all statistical analyses (Chi-Square tests, GLM models).
2. Create a new directory named `icorps_reanalysis_results/`.
3. Save all statistical output tables (e.g., `chi_square_test.csv`) and all figures (e.g., `Figure1_Funding_Spotlight.png`) into that new directory.

The results in the generated files will match those reported in the manuscript.
================================================================================