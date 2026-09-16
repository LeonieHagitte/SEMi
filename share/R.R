library(dplyr)
library(tidyr)
library(haven)
library(forcats)
library(lavaan)
library(OpenMx)
library(semtree)

df <- sharew9_rel9_0_0_dn %>%
  select(
    mergeid,
    country,
    dn003_,
    dn042_
  ) %>%
  
  # Add age at interview
  left_join(
    sharew9_rel9_0_0_cv_r %>%
      select(
        mergeid,
        age_int
      ),
    by = "mergeid"
  ) %>%
  
  # Add CASP-12 items
  left_join(
    sharew9_rel9_0_0_ac %>%
      select(
        mergeid,
        ac014_,
        ac015_,
        ac016_,
        ac017_,
        ac018_,
        ac019_,
        ac020_,
        ac021_,
        ac022_,
        ac023_,
        ac024_,
        ac025_
      ),
    by = "mergeid"
  )

df <- df %>%
  rename(
    dn003 = dn003_,
    dn042 = dn042_,
    ac014 = ac014_,
    ac015 = ac015_,
    ac016 = ac016_,
    ac017 = ac017_,
    ac018 = ac018_,
    ac019 = ac019_,
    ac020 = ac020_,
    ac021 = ac021_,
    ac022 = ac022_,
    ac023 = ac023_,
    ac024 = ac024_,
    ac025 = ac025_
  )

casp_negative <- c(
  "ac014", # age prevents activities
  "ac015", # out of control
  "ac016", # left out
  "ac018", # family responsibilities prevent
  "ac019"  # shortage of money prevents
)

casp_positive <- c(
  "ac017", # can do things wanted
  "ac020", # look forward to each day
  "ac021", # life has meaning
  "ac022", # look back with happiness
  "ac023", # full of energy
  "ac024", # life full of opportunities
  "ac025"  # future looks good
)

casp_raw <- sprintf("ac%03d", 14:25)

df_prep <- df %>%
  mutate(
    across(
      all_of(casp_raw),
      ~ {
        x <- as.numeric(.x)
        ifelse(x %in% 1:4, x, NA_real_)
      }
    )
  )

df_prep <- df_prep %>%
  mutate(
    # negatively worded items:
    # often = poor QoL, never = good QoL -> already correctly oriented
    across(
      all_of(casp_negative),
      as.numeric
    ),
    
    # positively worded items:
    # often = good QoL, never = poor QoL -> reverse
    across(
      all_of(casp_positive),
      ~ 5 - as.numeric(.x)
    )
  )

# ----------------------------------------------------------
df_prep <- df_prep %>%
  rename(
    C1 = ac014,
    C2 = ac015,
    C3 = ac016,
    
    A1 = ac017,
    A2 = ac018,
    A3 = ac019,
    
    P1 = ac020,
    P2 = ac021,
    P3 = ac022,
    
    S1 = ac023,
    S2 = ac024,
    S3 = ac025
  )

casp_items <- c(
  "C1", "C2", "C3",
  "A1", "A2", "A3",
  "P1", "P2", "P3",
  "S1", "S2", "S3"
)
# ----------------------------------------------------------
# Gender
#
# SHARE:
#   1 = male
#   2 = female
#
# Analysis coding:
#   0 = male
#   1 = female
# -------------------------------------------------------------------------

df_prep <- df_prep %>%
  mutate(
    female = case_when(
      as.numeric(dn042) == 1 ~ 0,
      as.numeric(dn042) == 2 ~ 1,
      TRUE ~ NA_real_
    )
  )

df_prep <- df_prep %>%
  mutate(
    age = as.numeric(age_int),
    age_z = (age - mean(age, na.rm = TRUE)) /
      sd(age, na.rm = TRUE)
  )

# -----------------------------------------------

df_prep <- df_prep %>%
  mutate(
    country_code = as.numeric(country),
    country_name = as.character(haven::as_factor(country))
  )

#df_prep %>%
#  distinct(country_code, country_name) %>%
#  arrange(country_code)

df_prep <- df_prep %>%
  mutate(
    country_name = factor(country_name)
  )

df_prep <- df_prep %>%
  mutate(
    cultural_cluster = case_when(
      
      country_name %in% c(
        "Denmark",
        "Finland",
        "Sweden"
      ) ~ "Nordic",
      
      country_name %in% c(
        "Austria",
        "Belgium",
        "France",
        "Germany",
        "Luxembourg",
        "Netherlands",
        "Switzerland"
      ) ~ "Western",
      
      country_name %in% c(
        "Cyprus",
        "Greece",
        "Italy",
        "Malta",
        "Portugal",
        "Spain"
      ) ~ "Southern",
      
      country_name %in% c(
        "Bulgaria",
        "Croatia",
        "Czech Republic",
        "Estonia",
        "Hungary",
        "Latvia",
        "Lithuania",
        "Poland",
        "Romania",
        "Slovakia",
        "Slovenia"
      ) ~ "Eastern",
      
      TRUE ~ NA_character_
    ),
    
    cultural_cluster = factor(
      cultural_cluster,
      levels = c(
        "Western",
        "Nordic",
        "Southern",
        "Eastern"
      )
    )
  )

df_prep <- df_prep %>%
  filter(!is.na(cultural_cluster))

df_prep <- df_prep %>%
  mutate(
    culture_nordic   = as.numeric(cultural_cluster == "Nordic"),
    culture_southern = as.numeric(cultural_cluster == "Southern"),
    culture_eastern  = as.numeric(cultural_cluster == "Eastern")
  )

#----------------------------------------------------------------
mnlfa_moderators <- c(
  "age_z",
  "female",
  "culture_nordic",
  "culture_southern",
  "culture_eastern"
)

tree_predictors <- c(
  "age_z",
  "female",
  "cultural_cluster"
)

# ----------------------------------------------------------------
df_analysis <- df_prep %>%
  filter(
    !is.na(cultural_cluster),
    !is.na(female),
    !is.na(age)
  ) %>%
  mutate(
    age_z = (age - mean(age)) / sd(age)
  )
#------------------------------------------------------------------
# Item-specific missingness
casp_missing <- df_analysis %>%
  summarise(
    across(
      all_of(casp_items),
      ~ mean(is.na(.x)) * 100
    )
  ) %>%
  pivot_longer(
    everything(),
    names_to = "item",
    values_to = "percent_missing"
  )

casp_missing

df_analysis <- df_analysis %>%
  mutate(
    n_casp_observed = rowSums(!is.na(across(all_of(casp_items)))),
    n_casp_missing  = 12 - n_casp_observed
  )

table(df_analysis$n_casp_missing, useNA = "ifany")

df_analysis %>%
  count(n_casp_missing) %>%
  mutate(
    percent = 100 * n / sum(n)
  )

df_analysis <- df_analysis %>%
  filter(n_casp_observed > 0)
#-------------------------------------------------------------------
df_analysis <- df_prep %>%
  filter(
    !is.na(cultural_cluster),
    !is.na(female),
    !is.na(age)
  ) %>%
  mutate(
    n_casp_observed = rowSums(
      !is.na(across(all_of(casp_items)))
    )
  ) %>%
  filter(n_casp_observed > 0) %>%
  mutate(
    age_z = (age - mean(age)) / sd(age)
  )
#-------------------------------------------------------------------
df_complete <- df_analysis %>%
  filter(n_casp_observed == 12)

nrow(df_complete)
# -------------------------------------------------------------------------
# PRELIMINARY DIAGNOSTIC CFA CHECKS
# Not part of the reported MNLFA / SEM-tree analyses
# -------------------------------------------------------------------------
casp_model <- '
  Control =~ C1 + C2 + C3
  Autonomy =~ A1 + A2 + A3
  Pleasure =~ P1 + P2 + P3
  SelfRealization =~ S1 + S2 + S3
'

fit_casp_wlsmv <- cfa(
  model = casp_model,
  data = df_analysis,
  ordered = casp_items,
  estimator = "WLSMV",
  std.lv = TRUE
)

summary(
  fit_casp_wlsmv,
  fit.measures = TRUE,
  standardized = TRUE
)

fitMeasures(
  fit_casp_wlsmv,
  c(
    "chisq", "df", "pvalue",
    "cfi", "tli",
    "rmsea", "rmsea.ci.lower", "rmsea.ci.upper",
    "srmr"
  )
)

standardizedSolution(fit_casp_wlsmv) %>%
  filter(op == "=~") %>%
  select(
    factor = lhs,
    item = rhs,
    loading = est.std,
    pvalue
  )

lavInspect(fit_casp_wlsmv, "cor.lv")
#------------
fit_casp_ml <- cfa(
  model = casp_model,
  data = df_analysis,
  estimator = "MLR",
  missing = "fiml",
  std.lv = TRUE
)

summary(
  fit_casp_ml,
  fit.measures = TRUE,
  standardized = TRUE
)

fitMeasures(
  fit_casp_ml,
  c(
    "cfi", "tli",
    "rmsea", "rmsea.ci.lower", "rmsea.ci.upper",
    "srmr"
  )
)

standardizedSolution(fit_casp_ml) %>%
  filter(op == "=~") %>%
  select(
    factor = lhs,
    item = rhs,
    loading = est.std
  )

lavInspect(fit_casp_ml, "cor.lv")

#---------------------------------------------------------------------------
df_prep <- df_prep %>%
  mutate(
    age = as.numeric(age_int),
    age = if_else(age < 0, NA_real_, age)
  )

df_analysis <- df_prep %>%
  filter(
    !is.na(cultural_cluster),
    !is.na(female),
    !is.na(age)
  ) %>%
  mutate(
    n_casp_observed = rowSums(
      !is.na(across(all_of(casp_items)))
    )
  ) %>%
  filter(n_casp_observed > 0)

# Save mean and SD used for standardizing age
age_mean <- mean(df_analysis$age)
age_sd   <- sd(df_analysis$age)

# Standardize age
df_analysis <- df_analysis %>%
  mutate(
    age_z = (age - age_mean) / age_sd
  )

df_complete <- df_analysis %>%
  filter(n_casp_observed == 12)

age_cutoffs <- tibble(
  split = c(
    "Metric: first age split",
    "Metric: lower age split",
    "Metric: upper age split",
    "Scalar: root age split",
    "Scalar: middle age split",
    "Scalar: lower age split",
    "Scalar: upper age split 1",
    "Scalar: upper age split 2"
  ),
  age_z = c(
    -0.0636058,
    -0.474243,
    0.962988,
    0.962988,
    0.0390536,
    -0.0636058,
    1.88692,
    1.78426
  )
) %>%
  mutate(
    age_years = age_mean + age_z * age_sd,
    age_years_rounded = round(age_years, 1)
  )

age_cutoffs
age_mean
age_sd
range(df_analysis$age)
