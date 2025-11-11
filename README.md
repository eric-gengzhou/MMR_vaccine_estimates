# County- and ZIP-Code–Level Small Area Estimation of MMR Vaccine Uptake

This repository contains the analytic code and supplementary materials for the paper:

> **Assessing MMR Vaccination Coverage Gaps in US Children with Digital Participatory Surveillance**  
> *Eric G. Zhou, John Brownstein, Benjamin Rader, 2025 (accepted and forthcoming in Nature Health)*

The code implements a small area estimation (SAE) with post-stratification workflow to generate county- and ZIP Code Tabulation Area (ZCTA)–level estimates of measles-mumps-rubella (MMR) vaccination coverage among children under five. The analysis integrates participatory surveillance data with demographic and contextual covariates to produce geographically smoothed estimates with uncertainty relfected in confidence intervals.

## Repository Contents

- R scripts for data preparation, model estimation, and post-stratification  
- Code for spatial smoothing and Local Moran’s I (LISA) cluster detection  
- Figure generation  
- Instructions for reproducing all main and supplemental figures

## Data Sources

The following publicly available datasets were used in this analysis:

| Dataset Key | Source & Access Link | Description |
|--------------|----------------------|--------------|
| **NHGIS_AGE_RACE_GENDER** | [NHGIS – U.S. Census Bureau](https://www.nhgis.org) | County-level tables combining *Hispanic Origin*, *Age*, *Sex*, and *Race* from the ACS 2019–2023 5-Year Estimates. Used for post-stratification population counts. |
| **CDC_COUNTY_COVID_VAX** | [CDC COVID-19 Vaccinations in the United States, County Level](https://data.cdc.gov/Vaccinations/COVID-19-Vaccinations-in-the-United-States-County/8xkx-amqh/about_data) | Cumulative county-level COVID-19 vaccination coverage (series complete and completeness metrics) used as a contextual covariate. |
| **CT_TOWN_VAX** | [Connecticut DPH: COVID-19 Vaccinations by Town and Age Group (Archived)](https://data.ct.gov/Health-and-Human-Services/COVID-19-Vaccinations-by-Town-and-Age-Group-ARCHIV/gngw-ukpw/about_data) | Town-level vaccination data from Connecticut Department of Public Health, aggregated to counties to supplement CDC coverage gaps. |
| **CDC_STATE_MMR_36M** | [KFF / CDC: Percent of Children Aged 0–35 Months Who Are Immunized](https://www.kff.org/other-health/state-indicator/percent-of-children-aged-0-35-months-who-are-immunized/) | State-level MMR vaccine coverage by age 36 months, used for validation against SAE-derived estimates. |
| **MEASLES_CASES_STATE** | [CDC Measles Surveillance Data](https://www.cdc.gov/measles/data-research/index.html) | Annual state-level confirmed measles case counts used to compare estimated vaccine uptake with outbreak patterns. |

## Output Data

Two prediction datasets generated from the final SAE model are included in this repository:

| File | Description |
|-------|--------------|
| [**county_pred_final.csv**](https://github.com/eric-gengzhou/MMR_vaccine_estimates/blob/main/county_pred_final.csv) | County-level posterior mean estimates of MMR vaccine uptake (≥1 dose), along with standard errors and 95% credible intervals. Each record corresponds to a county FIPS code. |
| [**zcta_pred_final.csv**](https://github.com/eric-gengzhou/MMR_vaccine_estimates/blob/main/zcta_pred_final.csv) | ZIP Code Tabulation Area (ZCTA)–level predicted MMR vaccine uptake (≥1 dose) and corresponding risk classification. Each record corresponds to a ZCTA5 code. |

Both files are directly hosted in this repository for reproducibility and visualization purposes.   

