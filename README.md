# Read Me
This repository stores the R code for simulations presented in the article:

Alternative adjustment for seasonality and long-term time-trend in time-series analysis for long-term environmental exposures and disease counts
Honghyok Kim, Jong-Tae Lee, Kelvin C. Fong, Michelle L. Bell. (Revised)


///Files description (Simulation)

test_data_public.csv: this is a test dataset for simulation. All numbers for death and environmental measurements were generated from random sampling based on the original data for public use. So, simulation results based on this dataset will not provide the same numbers in the main text.

S1. Predict.R: load test_data_public.csv and predict average PM10 time series and average mortality time series

S2. Set PM10 effects.R: specify hypothetical lag-structures of PM10 effects (lag0-730).

S3. Simulation Run(DLM0730).R: simulate distributed lag models (lag0–730) to estimate specified PM10 effects, with different adjustment methods for seasonality and time-trend

S4. Simulation Run(MA01).R: simulate two-day moving average models to estimate specified PM10 effects at only lag0-1, with different adjustment methods for seasonality and time-trend

S5. Simulation Figure in Main Text.R: Get Table 2 and Figure 3 in the main text.

///Files description (Real-data analysis)

These files provide example codes for real-data analysis including additional adjustment using multiple dummy variables. The codes were simplified for easier illustrations. Datasets for real-data analysis cannot be provided for public use.

R1. Regression and figures.R: run city-specific distributed lag models with different adjustments, pool city-specific estimates using multivariate meta-analysis, and get figures.

R2. Additional adjustment.R: run city-specific distributed lag models with additional adjustment using multiple dummy variables.


Questions or error reporting: honghyok.kim@yale.edu or honghyok@korea.ac.kr




