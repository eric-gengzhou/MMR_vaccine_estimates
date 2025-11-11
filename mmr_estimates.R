# SAE-Pipeline.R
# Small-Area Estimation of county-level MMR uptake with post-stratification
# Author: Eric G. Zhou, Benjamin Rader


suppressPackageStartupMessages({
  library(sf)
  library(dplyr)
  library(tidyr)
  library(stringr)
  library(ggplot2)
  library(spdep)
  library(glmnet)
  library(lme4)
  library(survey)
  library(MASS)       # mvrnorm
  library(readr)
  library(ggrepel)
  library(patchwork)
  library(janitor)
  library(grid)
  library(tigris)
})

options(scipen = 999)
sf::sf_use_s2(FALSE)
set.seed(123)

# 0) PARAMETERS & INPUTS
# Required data objects/files (document these in README):
# - us_counties: sf with columns geoid, statefp, n_res, covariates, geometry
# - mmr_data$variables: data.frame with respondent-level vars + weight column
# - town_region_xwalk: CT town->county crosswalk with county_fips
# - states_shp: sf state boundaries (for plots)
# - county_map_newshp: sf county geometries with geoid; est_mean will be joined if missing
# - ZCTA supplement (optional): zcta_clean (sf), states_conus (sf) 

paths <- list(
  covid_county_csv   = "covid_vax_county.csv",
  covid_ct_town_csv  = "covid_vax_ct.csv",
  nhgis_county_csv   = "nhgis0009_ds268_20235_county.csv",
  zcta_pred_csv      = "zcta_pred_final.csv",
  cdc_state_mmr_csv  = "cdc_state_mmr_36m.csv",      # columns: state_name, mmr_1d (0..1)
  measles_cases_state= "measles_cases_state.csv"     # columns: state_name, cases (count)
)

cfg <- list(
  min_county_obs        = 5L,
  lasso_nfolds          = 10L,
  n_sims                = 1000L,
  contiguous_state_fips = c(
    "01","04","05","06","08","09","10","11","12","13","16","17","18","19",
    "20","21","22","23","24","25","26","27","28","29","30","31","32","33",
    "34","35","36","37","38","39","40","41","42","44","45","46","47","48",
    "49","50","51","53","54","55","56"
  )
)

fips_map <- tibble::tibble(
  state = c("01","04","05","06","08","09","10","11","12","13","16","17","18","19","20","21","22","23","24","25",
            "26","27","28","29","30","31","32","33","34","35","36","37","38","39","40","41","42","44","45","46",
            "47","48","49","50","51","53","54","55","56"),
  state_name = c("Alabama","Arizona","Arkansas","California","Colorado","Connecticut","Delaware","District of Columbia",
                 "Florida","Georgia","Idaho","Illinois","Indiana","Iowa","Kansas","Kentucky","Louisiana","Maine",
                 "Maryland","Massachusetts","Michigan","Minnesota","Mississippi","Missouri","Montana","Nebraska",
                 "Nevada","New Hampshire","New Jersey","New Mexico","New York","North Carolina","North Dakota","Ohio",
                 "Oklahoma","Oregon","Pennsylvania","Rhode Island","South Carolina","South Dakota","Tennessee","Texas",
                 "Utah","Vermont","Virginia","Washington","West Virginia","Wisconsin","Wyoming"),
  state_abbr = c("AL","AZ","AR","CA","CO","CT","DE","DC","FL","GA","ID","IL","IN","IA","KS","KY","LA","ME","MD","MA",
                 "MI","MN","MS","MO","MT","NE","NV","NH","NJ","NM","NY","NC","ND","OH","OK","OR","PA","RI","SC","SD",
                 "TN","TX","UT","VT","VA","WA","WV","WI","WY")
)

DATA_SOURCES <- list(
  NHGIS_AGE_RACE_GENDER = "",     #https://www.nhgis.org/ select and download county-level 
                                  #	Hispanic Origin AND Age AND Sex AND Race table from the ACS 5-year 2023 data
  CDC_COUNTY_COVID_VAX  = "",     # https://data.cdc.gov/Vaccinations/COVID-19-Vaccinations-in-the-United-States-County/8xkx-amqh/about_data
  CT_TOWN_VAX           = "",     # https://data.ct.gov/Health-and-Human-Services/COVID-19-Vaccinations-by-Town-and-Age-Group-ARCHIV/gngw-ukpw/about_data
  CDC_STATE_MMR_36M     = "",     # https://www.kff.org/other-health/state-indicator/percent-of-children-aged-0-35-months-who-are-immunized/
  MEASLES_CASES_STATE   = ""      # https://www.cdc.gov/measles/data-research/index.html
)

legend_theme <- ggplot2::theme(
  legend.position = "bottom",
  legend.box = "horizontal",
  legend.key.size = grid::unit(0.3, "cm"),
  legend.text = element_text(size = 7),
  legend.title = element_text(size = 8),
  plot.title = element_text(hjust = -0.05, size = 10, face = "bold"),
  axis.text = element_blank(),
  axis.ticks = element_blank(),
  panel.grid = element_blank()
)

risk_colors <- c(
  "Very High Risk (<60%)" = "#DEEBF7",
  "High Risk (60–69%)"    = "#C6DBEF",
  "Medium Risk (70–79%)"  = "#6BAED6",
  "Low Risk (80–84%)"     = "#2171B5",
  "Lowest Risk (85%+)"    = "#08306B",
  "No Direct Estimates"   = "grey80"
)

# 1) DATA PREPARATION
prepare_covid_county <- function(path_csv, contiguous_state_fips) {
  stopifnot(file.exists(path_csv))
  raw <- read.csv(path_csv) |> janitor::clean_names()
  raw$date <- as.Date(raw$date, format = "%m/%d/%Y")
  base <- raw |>
    mutate(fips_chr = str_pad(as.character(fips), 5, pad = "0"),
           is_numeric = str_detect(fips_chr, "^\\d{5}$"),
           st_fips = substr(fips_chr, 1, 2)) |>
    filter(is_numeric, st_fips %in% contiguous_state_fips) |>
    mutate(fips = fips_chr) |>
    select(-fips_chr, -is_numeric, -st_fips)
  
  max_non_na <- function(x) if (all(is.na(x))) NA_real_ else max(x, na.rm = TRUE)
  out <- base |>
    group_by(fips) |>
    summarise(
      max_series   = max_non_na(series_complete_pop_pct),
      max_complete = max_non_na(completeness_pct),
      covid_vax    = pmax(max_series, max_complete, na.rm = TRUE),
      .groups = "drop"
    ) |>
    select(fips, covid_vax)
  
  out
}

prepare_covid_ct_county <- function(path_ct_town, town_region_xwalk) {
  stopifnot(file.exists(path_ct_town))
  stopifnot(all(c("town") %in% names(town_region_xwalk)))
  d <- read.csv(path_ct_town) |> janitor::clean_names()
  d <- d |>
    group_by(town) |>
    summarise(
      fully_vax = sum(fully_vaccinated_percent * population, na.rm = TRUE) / 
        sum(population, na.rm = TRUE),
      population = sum(population), .groups = "drop"
    ) |>
    filter(!is.na(fully_vax))
  
  xwalk <- town_region_xwalk |>
    mutate(town = str_replace_all(town, c("^W\\. " = "West ",
                                          "^E\\. " = "East ",
                                          "^N\\. " = "North ",
                                          "^S\\. " = "South ")))
  
  d2 <- left_join(d, xwalk, by = "town") |>
    group_by(county_fips) |>
    summarise(
      covid_vax = sum(fully_vax * population, na.rm = TRUE) / sum(population, na.rm = TRUE),
      .groups = "drop"
    ) |>
    mutate(fips = county_fips) |>
    select(fips, covid_vax)
  
  d2
}

merge_covid_vax <- function(cty_all, ct_ct) {
  cty_all2 <- cty_all[substr(cty_all$fips, 1, 2) != "09", ]
  bind_rows(cty_all2, ct_ct)
}

# 2) SPATIAL AGGREGATION (merge counties within-state until n_res >= min_county_obs)
aggregate_small_counties <- function(us_counties, min_threshold = 5) {
  uc <- us_counties |>
    mutate(centroid = st_centroid(geometry)) |>
    mutate(lon = st_coordinates(centroid)[, 1],
           lat = st_coordinates(centroid)[, 2])
  
  df <- uc |>
    st_drop_geometry() |>
    select(geoid, n_res, lon, lat, statefp) |>
    mutate(group = geoid) |>
    arrange(lon)
  
  group_n_res <- function(x) x |>
    group_by(group) |>
    summarise(agg_n_res = sum(n_res, na.rm = TRUE), .groups = "drop")
  
  iter <- 0; max_iter <- 1000
  repeat {
    iter <- iter + 1
    low <- group_n_res(df) |> filter(agg_n_res < min_threshold)
    if (nrow(low) == 0) break
    groups_df <- df |>
      group_by(group) |>
      summarise(state = first(statefp),
                agg_n_res = sum(n_res, na.rm = TRUE),
                centroid_lon = mean(lon, na.rm = TRUE),
                centroid_lat = mean(lat, na.rm = TRUE),
                .groups = "drop")
    for (g in low$group) {
      low_row <- groups_df |> filter(group == g)
      cand <- groups_df |> filter(group != g, state == low_row$state)
      if (nrow(cand) == 0) next
      distances <- sqrt((cand$centroid_lon - low_row$centroid_lon)^2 +
                          (cand$centroid_lat - low_row$centroid_lat)^2)
      nearest_group <- cand$group[which.min(distances)]
      df$group[df$group == g] <- nearest_group
    }
    if (iter > max_iter) { warning("Reached max iterations in county aggregation"); break }
  }
  df
}

# 3) VARIABLE SELECTION + UNIT-LEVEL ESTIMATION (LASSO -> GLMM)
select_vars_lasso <- function(dat, y_var, weight_var, candidate_vars, fac_vars = c("ruca_3")) {
  dd <- dat |>
    select(all_of(c(y_var, candidate_vars, fac_vars))) |>
    mutate(across(all_of(fac_vars), as.factor)) |>
    drop_na()
  
  x <- model.matrix(reformulate(c(setdiff(names(dd), y_var)), response = NULL), data = dd)[, -1, drop = FALSE]
  y <- dd[[y_var]]
  lasso_cv <- cv.glmnet(x, y, alpha = 1, family = "binomial", nfolds = cfg$lasso_nfolds)
  coefs <- coef(lasso_cv, s = "lambda.1se")
  keep <- rownames(coefs)[which(as.numeric(coefs) != 0)]
  setdiff(keep, "(Intercept)")
}

fit_unit_glmer <- function(dat, formula, weight_var) {
  dat$age_group <- recode(as.factor(dat$age_group), "over 50" = "50–59 years")
  dat$age_group <- factor(dat$age_group, levels = c("18-29 years","30-39 years","40-49 years","50–59 years"))
  dat$state <- substr(dat$county, 1, 2)
  
  glmer(
    formula,
    data = dat,
    family = binomial(link = "logit"),
    weights = dat[[weight_var]],
    na.action = na.exclude,
    control = glmerControl(optimizer = "bobyqa", optCtrl = list(maxfun = 200000))
  )
}

# 4) SPATIAL SMOOTHING OF COUNTY RANDOM EFFECTS (neighbor imputation for n_res==0)
impute_missing_re_by_neighbors <- function(us_counties, us_counties_mapped, nb = NULL) {
  if (is.null(nb)) nb <- spdep::poly2nb(us_counties, queen = TRUE)
  county_ids <- us_counties$geoid
  us_counties_mapped$state <- substr(us_counties_mapped$geoid, 1, 2)
  
  zero_geoids <- us_counties$geoid[us_counties$n_res == 0]
  us_counties_mapped$res0 <- as.integer(us_counties_mapped$geoid %in% zero_geoids)
  
  imputed_df <- lapply(zero_geoids, function(ct) {
    idx <- which(us_counties$geoid == ct)
    missing_state <- us_counties$state[idx]
    neighbor_indices <- nb[[idx]]
    neighbor_ids <- county_ids[neighbor_indices]
    neighbor_re <- us_counties_mapped |>
      filter(geoid %in% neighbor_ids, state == missing_state) |>
      pull(re_county)
    avg_effect <- if (length(neighbor_re) == 0) 0 else mean(neighbor_re, na.rm = TRUE)
    data.frame(geoid = ct, imputed_re = avg_effect, stringsAsFactors = FALSE)
  }) |> dplyr::bind_rows()
  
  us_counties_mapped |>
    left_join(imputed_df, by = "geoid") |>
    mutate(final_re = ifelse(res0 == 1, imputed_re, re_county))
}

# 5) POST-STRATIFICATION (NHGIS age x race x gender) & SIMULATION SEs
build_nhgis_cells <- function(path_csv) {
  stopifnot(file.exists(path_csv))
  nhgis <- read_csv(path_csv, show_col_types = FALSE)
  
  age_groups <- list(
    "18-29 years" = c("007","008","009","022","023","024"),
    "30-39 years" = c("010","011","025","026"),
    "40-49 years" = c("012","013","027","028"),
    "50–59 years" = c("014","015","029","030")
  )
  race_prefixes <- list(
    white    = "ASZHE",
    black    = "ASZBE",
    asian    = "ASZDE",
    hispanic = "ASZIE",
    aian     = "ASZCE",
    nhpi     = "ASZEE",
    other    = "ASZFE",
    multi    = "ASZGE"
  )
  extract_counts <- function(df, prefix, suffixes, gender) {
    gender_offset <- if (gender == "Male") 0 else 15
    cols <- paste0(prefix, sprintf("%03d", as.integer(suffixes) + gender_offset))
    cols <- cols[cols %in% names(df)]
    if (!length(cols)) return(rep(0, nrow(df)))
    rowSums(df[, cols, drop = FALSE], na.rm = TRUE)
  }
  
  nhgis <- nhgis |>
    mutate(county_fips = paste0(sprintf("%02d", as.integer(STATEA)),
                                sprintf("%03d", as.integer(COUNTYA))))
  
  results <- list()
  for (race in c("white","black","asian","hispanic","other")) {
    for (ag in names(age_groups)) {
      for (g in c("Male","Female")) {
        suf <- age_groups[[ag]]
        count <- dplyr::case_when(
          race %in% c("white","black","asian","hispanic") ~ extract_counts(nhgis, race_prefixes[[race]], suf, g),
          race == "other" ~ extract_counts(nhgis, race_prefixes$aian,  suf, g) +
            extract_counts(nhgis, race_prefixes$nhpi,  suf, g) +
            extract_counts(nhgis, race_prefixes$other, suf, g) +
            extract_counts(nhgis, race_prefixes$multi, suf, g)
        )
        results[[paste(race, ag, g, sep = "_")]] <- count
      }
    }
  }
  
  out <- nhgis |> select(county_fips) |> bind_cols(as.data.frame(results))
  long <- out |>
    pivot_longer(cols = -county_fips, names_to = c("race","age_raw","gender"),
                 names_pattern = "^(.*?)_(.*?)_(.*?)$", values_to = "N_pop") |>
    mutate(
      age_group = recode(age_raw,
                         "18.29.years" = "18-29 years",
                         "30.39.years" = "30-39 years",
                         "40.49.years" = "40-49 years",
                         "50.59.years" = "50–59 years"),
      race_5 = recode(race,
                      "white"="White", "black"="Black or African American",
                      "asian"="Asian", "hispanic"="Hispanic", "other"="Other"),
      age_group = factor(age_group, levels = c("18-29 years","30-39 years","40-49 years","50–59 years")),
      gender = factor(gender, levels = c("Female","Male"))
    ) |>
    select(county_fips, race_5, age_group, gender, N_pop)
  
  long
}

poststratify_county <- function(post_cells, unit_model, county_re, mmr_variables, covariate_vars) {
  post_cells <- post_cells |>
    filter(gender != "No answer") |>
    mutate(
      age_group = factor(age_group, levels = c("18-29 years","30-39 years","40-49 years","50–59 years")),
      race_5    = factor(race_5,    levels = levels(mmr_variables$race_5)),
      gender    = factor(gender,    levels = levels(mmr_variables$gender))
    )
  
  scale_means <- sapply(mmr_variables[, covariate_vars], mean, na.rm = TRUE)
  scale_sds   <- sapply(mmr_variables[, covariate_vars], sd,   na.rm = TRUE)
  for (v in covariate_vars) post_cells[[v]] <- (post_cells[[v]] - scale_means[[v]])/scale_sds[[v]]
  
  X <- model.matrix(
    ~ as.factor(age_group) + as.factor(race_5) + as.factor(gender) +
      scale(median_household_income) + scale(pct_white) + scale(covid_vax) +
      scale(pct_single_parent) + scale(pct_on_medicaid) + scale(dem_vote_share),
    data = post_cells
  )
  b <- fixef(unit_model)
  X <- X[, names(b), drop = FALSE]
  
  post_cells <- post_cells |>
    left_join(county_re |> select(county_fips = geoid, final_re), by = "county_fips") |>
    mutate(
      eta = as.vector(X %*% b + final_re),
      p   = plogis(eta)
    )
  
  post_cells |>
    group_by(county_fips) |>
    summarise(
      est_mean = sum(p * N_pop, na.rm = TRUE)/sum(N_pop, na.rm = TRUE),
      total_pop = sum(N_pop, na.rm = TRUE),
      .groups = "drop"
    )
}

simulate_se <- function(post_cells, unit_model, county_re, mmr_variables, covariate_vars, n_sims = 1000) {
  X <- model.matrix(
    ~ as.factor(age_group) + as.factor(race_5) + as.factor(gender) +
      scale(median_household_income) + scale(pct_white) + scale(covid_vax) +
      scale(pct_single_parent) + scale(pct_on_medicaid) + scale(dem_vote_share),
    data = post_cells
  )
  b_hat <- fixef(unit_model)
  Sigma <- vcov(unit_model)
  X <- X[, names(b_hat), drop = FALSE]
  
  pc <- post_cells |>
    left_join(county_re |> select(county_fips = geoid, final_re), by = "county_fips")
  
  county_ids <- unique(pc$county_fips)
  sim_mat <- matrix(NA_real_, nrow = length(county_ids), ncol = n_sims, dimnames = list(county_ids, NULL))
  
  for (s in seq_len(n_sims)) {
    beta_sim <- MASS::mvrnorm(1, mu = b_hat, Sigma = Sigma)
    eta      <- as.vector(X %*% beta_sim + pc$final_re)
    p_sim    <- plogis(eta)
    tmp <- pc |>
      mutate(p_sim = p_sim) |>
      group_by(county_fips) |>
      summarise(p_hat = sum(p_sim * N_pop, na.rm = TRUE)/sum(N_pop, na.rm = TRUE), .groups = "drop")
    sim_mat[tmp$county_fips, s] <- tmp$p_hat
  }
  
  tibble(county_fips = rownames(sim_mat)) |>
    mutate(
      se       = apply(sim_mat, 1, sd, na.rm = TRUE),
      est_mean = rowMeans(sim_mat, na.rm = TRUE),
      lower_95 = est_mean - 1.96 * se,
      upper_95 = est_mean + 1.96 * se
    )
}

make_state_cis <- function(sim_cell_draws, levels_map = NULL) {
  sim_cell_draws$state <- substr(sim_cell_draws$county_fips, 1, 2)
  state_sim <- sim_cell_draws |>
    group_by(state, sim) |>
    summarise(state_est = sum(p_sim * N_pop, na.rm = TRUE)/sum(N_pop, na.rm = TRUE), .groups = "drop")
  
  state_CI <- state_sim |>
    group_by(state) |>
    summarise(
      state_est_mean = mean(state_est, na.rm = TRUE),
      lower_95 = quantile(state_est, 0.025, na.rm = TRUE),
      upper_95 = quantile(state_est, 0.975, na.rm = TRUE),
      .groups = "drop"
    )
  
  if (!is.null(levels_map)) {
    state_CI <- state_CI |>
      mutate(state = sprintf("%02d", as.numeric(state))) |>
      left_join(levels_map, by = "state") |>
      rename(mean = state_est_mean, lower = lower_95, upper = upper_95) |>
      select(state_name, mean, lower, upper) |>
      arrange(state_name)
  }
  
  state_CI
}

plot_state_cis <- function(state_table, file = NULL) {
  p <- state_table |>
    mutate(state_name = reorder(state_name, mean)) |>
    ggplot(aes(x = mean, y = state_name)) +
    geom_point(size = 2) +
    geom_errorbarh(aes(xmin = lower, xmax = upper), height = 0.2) +
    scale_x_continuous(labels = scales::percent) +
    labs(x = "Estimated MMR Uptake (≥1 dose)", y = "", title = "") +
    theme_minimal(base_size = 10) +
    theme(panel.grid.major.y = element_blank())
  if (!is.null(file)) ggsave(file, p, width = 8.27, height = 11.69, units = "in", dpi = 600)
  p
}

lisa_county <- function(county_sf, value_col, states_shp) {
  nb <- poly2nb(county_sf, queen = TRUE)
  lw <- nb2listw(nb, style = "W", zero.policy = TRUE)
  vals <- county_sf[[value_col]]
  loc  <- localmoran(vals, lw, zero.policy = TRUE)
  county_sf$local_i <- loc[, "Ii"]
  county_sf$p_value <- loc[, "Pr(z != E(Ii))"]
  mm <- mean(vals, na.rm = TRUE)
  county_sf <- county_sf |>
    mutate(
      highlow_cat = case_when(
        vals > mm & lag.listw(lw, vals) > mm & p_value < 0.05 ~ "High-High",
        vals < mm & lag.listw(lw, vals) < mm & p_value < 0.05 ~ "Low-Low",
        vals < mm & lag.listw(lw, vals) > mm & p_value < 0.05 ~ "Low-High",
        vals > mm & lag.listw(lw, vals) < mm & p_value < 0.05 ~ "High-Low",
        TRUE ~ "Not Significant"
      ),
      highlow_cat = factor(highlow_cat, levels = c("Low-Low","Low-High","High-Low","High-High","Not Significant"))
    )
  
  p <- ggplot() +
    geom_sf(data = county_sf, aes(fill = highlow_cat), color = "gray90", size = 0.1) +
    geom_sf(data = states_shp, fill = NA, color = "black", size = 0.5) +
    scale_fill_manual(
      values = c("Low-Low"="#BD0026","Low-High"="#FDBB84","High-Low"="#B3CDE3",
                 "High-High"="#045A8D","Not Significant"="#FFFFFF"),
      name = "LISA Cluster"
    ) +
    labs(title = "", subtitle = "", caption = "") +
    theme_minimal(base_size = 10) +
    theme(axis.text = element_blank(), axis.ticks = element_blank(), panel.grid = element_blank())
  list(sf = county_sf, plot = p)
}

# 7) ZCTA-LEVEL MAPPING (Supplement; optional)
zcta_maps <- function(zcta_clean, states_conus, xlim2163, ylim2163) {
  risk_levels_5 <- names(risk_colors)[1:5]
  legend_levels <- c(risk_levels_5, "Uninhabited")
  risk_colors_legend <- c(risk_colors[risk_levels_5], "Uninhabited" = "white")
  
  zcta_plot <- zcta_clean |>
    mutate(
      risk_level_raw = cut(est_mean, breaks = c(-Inf, .60, .70, .80, .85, Inf),
                           labels = risk_levels_5, ordered_result = TRUE),
      risk_level = ifelse(is.na(risk_level_raw), "Uninhabited", as.character(risk_level_raw)),
      risk_level = factor(risk_level, levels = legend_levels, ordered = TRUE)
    ) |>
    select(-risk_level_raw)
  
  p <- ggplot() +
    annotate("rect", xmin = xlim2163[1], xmax = xlim2163[2],
             ymin = ylim2163[1], ymax = ylim2163[2], fill = "white") +
    geom_sf(data = zcta_plot, aes(fill = risk_level), color = NA) +
    geom_sf(data = st_boundary(states_conus), color = "white", linewidth = 0.15) +
    scale_fill_manual(values = risk_colors_legend, breaks = legend_levels, drop = FALSE,
                      name = "Measles Risk Level (vaccination): ") +
    coord_sf(crs = st_crs(2163), xlim = xlim2163, ylim = ylim2163, expand = FALSE) +
    theme_void(base_size = 11) +
    theme(legend.position = "right")
  
  p
}

# 8) EXHIBIT HELPERS
build_state_compare <- function(state_table, paths, fips_map) {
  stopifnot(file.exists(paths$cdc_state_mmr_csv),
            file.exists(paths$measles_cases_state))
  cdc_state <- readr::read_csv(paths$cdc_state_mmr_csv, show_col_types = FALSE) %>%
    dplyr::select(state_name, mmr_1d)
  measles   <- readr::read_csv(paths$measles_cases_state, show_col_types = FALSE) %>%
    dplyr::select(state_name, cases)
  
  out <- state_table %>%
    dplyr::inner_join(cdc_state, by = "state_name") %>%
    dplyr::left_join(measles, by = "state_name") %>%
    dplyr::left_join(fips_map, by = "state_name") %>%
    dplyr::mutate(
      cases = ifelse(is.na(cases), 0L, as.integer(cases)),
      case_group = dplyr::case_when(
        cases == 0 ~ "0 (No cases)",
        cases <= 10 ~ "1–10",
        cases <= 50 ~ "11–50",
        cases <= 99 ~ "51–99",
        cases >= 100 ~ "100+",
        TRUE ~ "0 (No cases)"
      ),
      case_group = factor(case_group, levels = c("0 (No cases)","1–10","11–50","51–99","100+"))
    ) %>%
    dplyr::rename(est_mean = mean)
  
  out
}

plot_fig1a_glmm <- function(county_map_newshp, states_shp, out_file = "Figure_1a_GLMM_Risk_A4_600dpi.png") {
  if (!"risk_level" %in% names(county_map_newshp) & "est_mean" %in% names(county_map_newshp)) {
    rl_levels <- names(risk_colors)[1:5]
    county_map_newshp <- county_map_newshp %>%
      dplyr::mutate(
        risk_level_raw = cut(est_mean, breaks = c(-Inf, .60, .70, .80, .85, Inf),
                             labels = rl_levels, ordered_result = TRUE),
        risk_level = ifelse(is.na(risk_level_raw), "No Direct Estimates", as.character(risk_level_raw)),
        risk_level = factor(risk_level, levels = c(rl_levels, "No Direct Estimates"))
      ) %>%
      dplyr::select(-risk_level_raw)
  }
  
  p <- ggplot(data = county_map_newshp) +
    geom_sf(aes(fill = risk_level), color = "gray90", size = 0.1) +
    geom_sf(data = states_shp, fill = NA, color = "white", size = 0.6) +
    scale_fill_manual(values = risk_colors, name = "Measles Risk Level (vaccination): ", na.value = "grey80") +
    labs(title = "a. County-Level Measles Risk Classification Based on Estimated MMR Vaccine Uptake (≥1 Dose)\n") +
    theme_minimal(base_size = 10) + legend_theme
  
  ggsave(out_file, p, width = 8.27, height = 11.69, units = "in", dpi = 600)
  p
}

save_fig1b_lisa <- function(lisa_obj, out_file = "Figure_1b_LISA_Counties_A4_600dpi.png") {
  ggsave(out_file, lisa_obj$plot, width = 8.27, height = 11.69, units = "in", dpi = 600)
}

plot_fig2_state <- function(state_table, out_file = "Figure_2_State_Dotplot_A4_600dpi.png") {
  p <- state_table %>%
    dplyr::mutate(state_name = reorder(state_name, mean)) %>%
    ggplot(aes(x = mean, y = state_name)) +
    geom_point(size = 2) +
    geom_errorbarh(aes(xmin = lower, xmax = upper), height = 0.2) +
    scale_x_continuous(labels = scales::percent_format(accuracy = 1)) +
    labs(x = "Estimated MMR Uptake (≥1 Dose)", y = "", title = "") +
    theme_minimal(base_size = 10) +
    theme(panel.grid.major.y = element_blank())
  ggsave(out_file, p, width = 8.27, height = 11.69, units = "in", dpi = 600)
  p
}

plot_fig3_state_compare <- function(state_compare, out_file = "Figure_3_State_Comparison_A4_600dpi.png") {
  sc <- state_compare
  if (!"state_abbr" %in% names(sc)) {
    sc$state_abbr <- sc$state_name
  }
  p <- ggplot(sc) +
    geom_point(aes(x = est_mean, y = mmr_1d, size = case_group, color = case_group), alpha = 0.7) +
    ggrepel::geom_text_repel(aes(x = est_mean, y = mmr_1d, label = state_abbr),
                             size = 4, box.padding = 0.3, max.overlaps = 100) +
    geom_abline(slope = 1, intercept = 0, linetype = "dashed", color = "black") +
    scale_size_manual(
      name = "Measles Cases",
      values = c("0 (No cases)" = 2, "1–10" = 3, "11–50" = 8, "51–99" = 12, "100+" = 16),
      drop = FALSE
    ) +
    scale_color_manual(
      name = "Measles Cases",
      values = c("0 (No cases)" = "blue", "1–10" = "red", "11–50" = "red", "51–99" = "red", "100+" = "red")
    ) +
    coord_cartesian(xlim = c(0.55, 1.00), ylim = c(0.55, 1.00), expand = FALSE) +
    theme_minimal(base_size = 13) +
    labs(
      x = "SAE-estimated State-level MMR Vaccine Rates",
      y = "CDC-Reported MMR Vaccine Rates (≥1 dose) by 36 months",
      title = ""
    ) +
    theme(legend.position = "right",
          plot.title = element_text(face = "bold"),
          plot.subtitle = element_text(margin = margin(b = 10)))
  ggsave(out_file, p, width = 8.27, height = 11.69, units = "in", dpi = 600)
  p
}

# 9) MAIN DRIVER
run_pipeline <- function(
    us_counties,                 # sf with county covariates + n_res
    mmr_data,                    # list with $variables (respondents)
    town_region_xwalk,           # CT helper table with county_fips
    states_shp,                  # sf state boundaries (for plots)
    county_map_newshp,           # sf county map; est_mean joined if missing
    paths = paths, cfg = cfg
) {
  covid_all  <- prepare_covid_county(paths$covid_county_csv, cfg$contiguous_state_fips)
  covid_ct   <- prepare_covid_ct_county(paths$covid_ct_town_csv, town_region_xwalk)
  covid_cty  <- merge_covid_vax(covid_all, covid_ct)
  
  groups_df  <- aggregate_small_counties(us_counties, cfg$min_county_obs)
  
  county_vote <- us_counties |> st_drop_geometry() |> select(geoid, dem_vote_share)
  mmr_data$variables <- left_join(mmr_data$variables, county_vote, by = c("county" = "geoid"))
  
  cand_vars <- c(
    "doctor_per_capita","totrate_2020","pop_density","median_age","total_population",
    "median_household_income","pct_white","pct_hispanic","pct_foreign_born","pct_children",
    "pct_high_school_only","pct_bachelors_degree","pct_single_parent","pct_snap","pct_unemployed",
    "pct_on_medicaid","pct_poverty_children","dem_vote_share","ruca_3"
  )
  keep <- select_vars_lasso(mmr_data$variables, y_var = "m_vax_u5",
                            weight_var = "weight_daily_national_18plus",
                            candidate_vars = cand_vars, fac_vars = "ruca_3")
  
  mmr_data$variables$age_group <- recode(as.factor(mmr_data$variables$age_group), "over 50" = "50–59 years")
  mmr_data$variables$age_group <- factor(mmr_data$variables$age_group,
                                         levels = c("18-29 years","30-39 years","40-49 years","50–59 years"))
  mmr_data$variables$state <- substr(mmr_data$variables$county, 1, 2)
  mmr_data$variables <- left_join(mmr_data$variables, covid_cty, by = c("county" = "fips"))
  
  formula_glmm <- as.formula(
    m_vax_u5 ~ as.factor(age_group) + as.factor(race_5) + as.factor(gender) +
      scale(median_household_income) + scale(pct_white) + scale(pct_single_parent) +
      scale(pct_on_medicaid) + scale(dem_vote_share) + scale(covid_vax) +
      (1 | state/group)
  )
  
  unit_model <- fit_unit_glmer(mmr_data$variables, formula_glmm, "weight_daily_national_18plus")
  
  re_list <- ranef(unit_model)$`group:state`
  re_df <- tibble::rownames_to_column(as.data.frame(re_list), var = "group_state") |>
    rename(re_county = `(Intercept)`) |>
    mutate(group_id = substr(group_state, 1, 5)) |>
    select(group_id, re_county)
  
  counties_map <- groups_df |>
    left_join(re_df, by = c("group" = "group_id")) |>
    left_join(us_counties |> st_drop_geometry() |> select(geoid, n_res, state = statefp), by = c("geoid"))
  
  nb_cty <- spdep::poly2nb(us_counties, queen = TRUE)
  counties_mapped <- left_join(us_counties |> st_drop_geometry() |> select(geoid), counties_map, by = "geoid")
  counties_imputed <- impute_missing_re_by_neighbors(us_counties, counties_mapped, nb_cty)
  county_re <- counties_imputed |> select(geoid, final_re)
  
  nhgis_long <- build_nhgis_cells(paths$nhgis_county_csv)
  
  age_levels   <- levels(mmr_data$variables$age_group)
  race_levels  <- levels(mmr_data$variables$race_5)
  gender_levels<- levels(mmr_data$variables$gender)
  post_cells <- expand.grid(
    county_fips = us_counties$geoid,
    age_group   = age_levels,
    race_5      = race_levels,
    gender      = gender_levels,
    stringsAsFactors = FALSE
  ) |>
    as_tibble() |>
    left_join(nhgis_long |> rename(N_pop = N_pop), by = c("county_fips","age_group","race_5","gender")) |>
    mutate(N_pop = ifelse(is.na(N_pop), 0, N_pop))
  
  ct_covs <- us_counties |> st_drop_geometry() |>
    select(county_fips = geoid, median_household_income, pct_white, dem_vote_share,
           pct_single_parent, pct_on_medicaid) |>
    left_join(county_re, by = c("county_fips" = "geoid")) |>
    left_join(covid_cty, by = c("county_fips" = "fips"))
  
  post_cells <- left_join(post_cells, ct_covs, by = "county_fips")
  
  covariate_vars <- c("median_household_income","pct_white","pct_single_parent","pct_on_medicaid","dem_vote_share","covid_vax")
  county_est <- poststratify_county(post_cells, unit_model, county_re, mmr_data$variables, covariate_vars)
  
  X <- model.matrix(
    ~ as.factor(age_group) + as.factor(race_5) + as.factor(gender) +
      scale(median_household_income) + scale(pct_white) + scale(covid_vax) +
      scale(pct_single_parent) + scale(pct_on_medicaid) + scale(dem_vote_share),
    data = post_cells
  )
  b_hat <- fixef(unit_model); Sigma <- vcov(unit_model)
  X <- X[, names(b_hat), drop = TRUE]
  pc_draw <- post_cells |>
    left_join(county_re |> select(county_fips = geoid, final_re), by = "county_fips") |>
    mutate(cell_id = row_number())
  
  sim_list <- vector("list", cfg$n_sims)
  for (s in seq_len(cfg$n_sims)) {
    beta_sim <- MASS::mvrnorm(1, mu = b_hat, Sigma = Sigma)
    eta      <- as.vector(X %*% beta_sim + pc_draw$final_re)
    sim_list[[s]] <- tibble(cell_id = pc_draw$cell_id, p_sim = plogis(eta), sim = s)
  }
  sim_long <- bind_rows(sim_list) |>
    left_join(pc_draw |> select(cell_id, county_fips, age_group, race_5, gender, N_pop), by = "cell_id")
  
  sim_mat <- sim_long |>
    group_by(sim, county_fips) |>
    summarise(p_hat = sum(p_sim * N_pop, na.rm = TRUE)/sum(N_pop, na.rm = TRUE), .groups = "drop") |>
    pivot_wider(names_from = sim, values_from = p_hat)
  
  se_tbl <- tibble(county_fips = sim_mat$county_fips) |>
    mutate(
      se       = apply(as.matrix(sim_mat[,-1]), 1, sd, na.rm = TRUE),
      est_mean = rowMeans(as.matrix(sim_mat[,-1]), na.rm = TRUE),
      lower_95 = est_mean - 1.96*se,
      upper_95 = est_mean + 1.96*se
    )
  
  county_out <- county_est |>
    select(county_fips, est_mean) |>
    left_join(se_tbl |> select(county_fips, se, lower_95, upper_95), by = "county_fips")
  
  state_table <- make_state_cis(sim_cell_draws = sim_long, levels_map = fips_map)
  
  if (!"est_mean" %in% names(county_map_newshp)) {
    county_map_newshp <- county_map_newshp %>%
      dplyr::left_join(county_out %>% dplyr::select(county_fips, est_mean),
                       by = c("geoid" = "county_fips"))
  }
  
  fig1a <- plot_fig1a_glmm(county_map_newshp, states_shp)
  
  lisa_obj <- lisa_county(county_map_newshp %>% dplyr::mutate(vals = est_mean), value_col = "vals", states_shp = states_shp)
  save_fig1b_lisa(lisa_obj)
  
  fig2 <- plot_fig2_state(state_table)
  
  state_compare <- build_state_compare(state_table, paths, fips_map)
  fig3 <- plot_fig3_state_compare(state_compare)
  
  readr::write_csv(county_out, "Table_County_Estimates.csv")
  readr::write_csv(state_table, "Table_State_Estimates.csv")
  
  list(
    unit_model     = unit_model,
    county_est     = county_out,
    state_table    = state_table,
    lisa_county_sf = lisa_obj$sf,
    state_compare  = state_compare,
    fig1a          = fig1a,
    fig2          = fig2,
    fig3          = fig3
  )
}

# 10) OPTIONAL: ZCTA SUPPLEMENT
# Dependencies for this section
suppressPackageStartupMessages({
  library(tigris)
  library(sf)
  library(ggplot2)
})

options(tigris_use_cache = TRUE)
sf::sf_use_s2(FALSE)

# Recommended bounding box for continental U.S. in EPSG:2163
xlim2163 <- c(-2400000, 2550000)
ylim2163 <- c(-2100000,  730000)

# State boundaries in EPSG:2163 and CONUS-only
states_all  <- tigris::states(cb = TRUE) |> sf::st_transform(2163)
states_conus <- states_all |>
  dplyr::filter(!STUSPS %in% c("AK", "HI", "PR", "VI", "GU", "MP", "AS"))

# Wrapper to export ZCTA risk map and ZCTA LISA map at A4, 600 dpi
# Assumes you have:
# - zcta_maps(zcta_clean, states_conus, xlim2163, ylim2163)   # returns a ggplot
# - zcta_lisa_map(zcta_clean, states_conus, xlim2163, ylim2163) # returns list(sf=..., plot=ggplot)
run_zcta_supplement <- function(zcta_clean,
                                states_conus = states_conus,
                                xlim2163 = xlim2163,
                                ylim2163 = ylim2163) {
  # Risk choropleth
  p_risk <- zcta_maps(zcta_clean, states_conus, xlim2163, ylim2163)
  ggsave("Figure_ZCTA_Risk_A4_600dpi.png", p_risk,
         width = 8.27, height = 11.69, units = "in", dpi = 600)
  
  # LISA map
  lisa_obj <- zcta_lisa_map(zcta_clean, states_conus, xlim2163, ylim2163)
  ggsave("Figure_ZCTA_LISA_A4_600dpi.png", lisa_obj$plot,
         width = 8.27, height = 11.69, units = "in", dpi = 600)
  
  invisible(list(risk_plot = p_risk, lisa = lisa_obj))
}

# Example usage 
# res <- run_pipeline(us_counties, mmr_data, town_region_xwalk, states_shp, county_map_newshp)
# run_zcta_supplement(zcta_clean, states_conus, xlim2163, ylim2163)
# Outputs created in working directory:
#   Figure_1a_GLMM_Risk_A4_600dpi.png
#   Figure_1b_LISA_Counties_A4_600dpi.png
#   Figure_2_State_Dotplot_A4_600dpi.png
#   Figure_3_State_Comparison_A4_600dpi.png
#   Table_County_Estimates.csv
#   Table_State_Estimates.csv