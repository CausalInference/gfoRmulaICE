source("../ICE.R")
source("../helper.R")
source("../plot.R")
source("../ICE_pool.R")
source("../ICE_stratify.R")
source("../bootstrap_ice.R")

data <- readRDS("survival_censor_comp.rds")

library(tidyverse)

fit_hazard_strat <- ice(data = data,
                           K = 4,
                           id = "id",
                           time_name = "t0",
                           outcome_name = "Y",
                           censor_name = "C",
                           competing_name = "D",
                           estimator = strat(hazard = T),
                           comp_effect = 1,
                           outcome_model = Y ~ L1 + L2,
                           competing_model = D ~ L1 + L2,
                           censor_model = C ~ L1 + L2,
                           ref_idx = 0,
                           int_descript = c("Static Intervention",
                                            "Dynamic Intervention"),
                           intervention1.A1 = list(static(0)),
                           intervention1.A2 = list(static(1)),
                           intervention2.A1 = list(dynamic("absorbing", "L1", "=", 1)),
                           intervention2.A2 = list(dynamic("compare", "L1", "=", 0)),
                        outcomeModel.1 = Y ~ L1,
                        compModel.2 = D ~ L1
)

fit_hazard_strat$summary

# NP Risk Classical Stratified ICE Risk Ratio Risk Difference
# Natural Course (Direct Effect)       0.26947                  0.26967    1.00000         0.00000
# Static Intervention (Direct Effect)       NA                  0.33918    1.25776         0.06951
# Dynamic Intervention (Direct Effect)      NA                  0.18494    0.68580        -0.08473

# NP Risk Classical Stratified ICE Risk Ratio Risk Difference
# Natural Course (Direct Effect)       0.26947                  0.26967    1.00000         0.00000
# Static Intervention (Direct Effect)       NA                  0.33731    1.25083         0.06764
# Dynamic Intervention (Direct Effect)      NA                  0.18494    0.68580        -0.08473

dynamic_cat <- case_when(data$L2 < 0 ~ 1,
                         data$L2 >= 0 & data$L2 < 2 ~ 2,
                         T ~ 3)

# Please see the attached three data sets
# survival_censor_comp.csv - survival outcome with censoring and competing event
# survival_censor.csv - survival outcome with censoring event
# survival_comp.csv - survival outcome with competing event
# Number of time points = 4
# Two treatment variables: A1 (categorical with 1, 2, 3) and A2 (binary)
# Covariates: L1 (binary) and L2 (continuous)
# Censoring event: C
# Competing event: D
# Please run through the following interventions on each data set
# Intervention 1: static always assign 3 for A1 and static always treat (always assign 1) for A2
# Intervention 2: dynamic treat when L2 < 0 assign 1, when L2 >= 0 and L2 < 2 assign 2, and assign 3 otherwise for A1,
# dynamic treat when L1 = 0 for A2
# Intervention 3 (pooled ICE only): threshold intervention with threshold value 1 for A1, and threshold intervention
# with threshold value 3 for A2
# Intervention 4 (pooled ICE only): dynamic treat when L1 = 1 with uniform grace period 2 period for A1
# dynamic treat when L2 = 2 with uniform grace period 2 period for A2
# Models:
# Pooled ICE: Y ~ L1 + L2 + A1 + A2, C ~ L1 + L2 + A1 + A2, D ~ L1 + L2 + A1 + A2
# Stratified ICE: Y ~ L1 + L2, C ~ L1 + L2, D ~ L1 + L2
# Weighted ICE: Y ~ L1 + L2, C ~ L1 + L2, D ~ L1 + L2, A1 ~ L1 + L2, A2 ~ L1 + L2

static_A1 <- case_when(data$t0 %in% 1:2 ~ 3,
                       T ~ data$A1)
static_A2 <- case_when(data$t0 %in% 1:2 ~ 1,
                       T ~ data$A2)

data$t0 <- case_when(data$t0 == 1 ~ "day1",
                     data$t0 == 2 ~ "day2",
                     data$t0 == 3 ~ "day3",
                     data$t0 == 4 ~ "day4")

data$t0 <- as.numeric(data$t0) + 1

dynamic_cat <- case_when(data$L2 < 0 ~ 1,
                         data$L2 >= 0 & data$L2 < 2 ~ 2,
                         T ~ 3)

fit_classical_pool <- ice(
  data = data,
  time_points = 4,
  id = "id",
  time_name = "t0",
  outcome_name = "Y",
  censor_name = "C",
  compevent_name = "D",
  estimator = pool(hazard = F),
  comp_effect = 0,
  outcome_model = Y ~ L1 + L2 + A1 + A2,
  competing_model = D ~ L1 + L2 + A1 + A2,
  censor_model = C ~ L1 + L2 + A1 + A2,
  ref_idx = 0,
  int_descript = c("Static Intervention",
                   "Dynamic Intervention"),
  intervention1.A1 = list(static(3)),
  intervention1.A2 = list(static(1)),
  intervention2.A1 = list(dynamic_cat),
  intervention2.A2 = list(dynamic("compare", "L1", "=", 0)),
  nsamples = 2, parallel = T, ci_method = "percentile",
)

# pooled - matched
## all time points
# NP Risk Classical Pooled ICE Risk Ratio Risk Difference
# Natural Course (Direct Effect)      0.26947              0.26874    1.00000         0.00000
# Static Intervention (Direct Effect)      NA              0.20222    0.75249        -0.06651

# NP Risk Classical Pooled ICE Risk Ratio Risk Difference
# Natural Course (Direct Effect)      0.26947              0.26874    1.00000         0.00000
# Static Intervention (Direct Effect)      NA              0.24935    0.92784        -0.01939

## 1:2 time points
# NP Risk Classical Pooled ICE Risk Ratio Risk Difference
# Natural Course (Direct Effect)      0.26947              0.26874    1.00000         0.00000
# Static Intervention (Direct Effect)      NA              0.24319    0.90494        -0.02555

# NP Risk Classical Pooled ICE Risk Ratio Risk Difference
# Natural Course (Direct Effect)      0.26883              0.26874    1.00000         0.00000
# Static Intervention (Direct Effect)      NA              0.26011    0.96788        -0.00863

# stratified
## all time points

# NP Risk Classical Stratified ICE Risk Ratio Risk Difference
# Natural Course (Direct Effect)      0.26947                  0.26967     1.0000          0.0000
# Static Intervention (Direct Effect)      NA                  0.19907     0.7382         -0.0706

# NP Risk Classical Stratified ICE Risk Ratio Risk Difference
# Natural Course (Direct Effect)      0.26883                  0.26967    1.00000         0.00000
# Static Intervention (Direct Effect)      NA                  0.22639    0.83951        -0.04328

## 1:2 time points
# NP Risk Classical Stratified ICE Risk Ratio Risk Difference
# Natural Course (Direct Effect)      0.26947                  0.26967    1.00000         0.00000
# Static Intervention (Direct Effect)      NA                  0.24365    0.90354        -0.02601

# NP Risk Classical Stratified ICE Risk Ratio Risk Difference
# Natural Course (Direct Effect)      0.26947                  0.26967    1.00000         0.00000
# Static Intervention (Direct Effect)      NA                  0.24441    0.90634        -0.02526

# NP Risk Classical Stratified ICE Risk Ratio Risk Difference
# Natural Course (Direct Effect)      0.26947                  0.26967    1.00000         0.00000
# Static Intervention (Direct Effect)      NA                  0.25865    0.95916        -0.01101

# pooled classical ICE

static_A1 <- case_when(data$t0 %in% 1:2 ~ 3,
                                T ~ data$A1)
static_A2 <- case_when(data$t0 %in% 1:2 ~ 1,
                                T ~ data$A2)

dynamic_cat_A1 <- case_when((data$L2 < 0) & (data$t0 %in% 2:3) ~ 1,
                                     data$L2 >= 0 & data$L2 < 2 & data$t0 %in% 2:3 ~ 2,
                                     data$t0 %in% 2:3 ~ 3,
                                     T ~ data$A1)

dynamic_A2 <- case_when((data$L1 == 0) & (data$t0 %in% 2:3) ~ 1,
                                     data$t0 %in% 2:3 ~ 0,
                                     T ~ data$A2)

## need to revise so the int_time option works

fit_classical_pool <- ice(
  data = data,
  K = 4,
  id = "id",
  time_name = "t0",
  outcome_name = "Y",
  censor_name = "C",
  competing_name = "D",
  estimator = pool(hazard = F),
  comp_effect = 0,
  outcome_model = Y ~ L1 + L2 + A1 + A2,
  censor_model = C ~ L1 + L2 + A1 + A2,
  ref_idx = 0,
  int_descript = c("Static Intervention",
                   "Dynamic Intervention"),
  intervention1.A1 = list(static_A1),
  intervention1.A2 = list(static_A2),
  intervention2.A1 = list(dynamic_cat_A1),
  intervention2.A2 = list(dynamic_A2)
)

fit_classical_pool <- ice(
  data = data,
  K = 4,
  id = "id",
  time_name = "t0",
  outcome_name = "Y",
  censor_name = "C",
  competing_name = "D",
  estimator = pool(hazard = F),
  comp_effect = 0,
  outcome_model = Y ~ L1 + L2 + A1 + A2,
  censor_model = C ~ L1 + L2 + A1 + A2,
  ref_idx = 0,
  int_descript = c("Static Intervention",
                   "Dynamic Intervention"),
  intervention1.A1 = list(static(3)),
  intervention1.A2 = list(static(1)),
  intervention2.A1 = list(dynamic_cat),
  intervention2.A2 = list(dynamic("compare", "L1", "=", 0)),
  bootstrap = T, nboot = 1000, parallel = T
)

plot_risk(fit_classical_pool)

View(fit_classical_pool$summary)

fit_classical_pool_total_effect <- ice(
  data = data,
  K = 4,
  id = "id",
  time_name = "t0",
  outcome_name = "Y",
  censor_name = "C",
  competing_name = "D",
  estimator = pool(hazard = F),
  comp_effect = 1,
  outcome_model = Y ~ L1 + L2 + A1 + A2,
  censor_model = C ~ L1 + L2 + A1 + A2,
  ref_idx = 0,
  int_descript = c("Static Intervention",
                   "Dynamic Intervention"),
  intervention1.A1 = list(static(3)),
  intervention1.A2 = list(static(1)),
  intervention2.A1 = list(dynamic_cat),
  intervention2.A2 = list(dynamic("compare", "L1", "=", 0))
)

summary_table(
  fit_classical_pool,
  fit_classical_pool_total_effect)

plot_risk(
  fit_classical_pool,
  fit_classical_pool_total_effect)

fit_classical_pool <- ice(data = data,
                          K = 4,
                          id = "id",
                          time_name = "t0",
                          outcome_name = "Y",
                          censor_name = "C",
                          competing_name = "D",
                          estimator = pool(hazard = F),
                          comp_effect = 0,
                          outcome_model = Y ~ L1 + L2 + A1 + A2,
                          censor_model = C ~ L1 + L2 + A1 + A2,
                          ref_idx = 0,
                          int_descript = c("Static Intervention",
                                           "Dynamic Intervention"),
                          intervention1.A1 = list(static(3)),
                          intervention1.A2 = list(static(1)),
                          intervention2.A1 = list(dynamic_cat),
                          intervention2.A2 = list(dynamic("compare", "L1", "=", 0))
)

plot_risk(ice_total_effect, fit_classical_ice)

plot_risk(fit_classical_ice)

ice(data = data,
    K = 4,
    id = "id",
    time_name = "t0",
    outcome_name = "Y",
    censor_name = "C",
    competing_name = "D",
    estimator = pool(hazard = F),
    comp_effect = 0,
    outcome_model = Y ~ L1 + L2 + A1 + A2,
    competing_model = D ~ L1 + L2 + A1 + A2,
    censor_model = C ~ L1 + L2 + A1 + A2,
    ref_idx = 0,
    int_descript = c("Static Intervention",
                     "Dynamic Intervention",
                     "Threshold Intervention"),
    intervention1.A1 = list(static(3), 1:2),
    intervention1.A2 = list(static(1)),
    intervention2.A1 = list(dynamic_cat),
    intervention2.A2 = list(dynamic("compare", "L1", "=", 0), 2:3),
    intervention3.A1 = list(threshold(2, Inf))
)

ice_test$summary

# result

# IPW Risk Classical Pooled ICE Risk Ratio Risk Difference
# Natural Course (Total Effect)           0.22949              0.22776    1.00000         0.00000
# Static Intervention (Total Effect)           NA              0.21449    0.94175        -0.01327
# Dynamic Intervention (Total Effect)          NA              0.27596    1.21164         0.04820
# Threshold Intervention (Total Effect)        NA              0.05842    0.25651        -0.16933
# Natural Course (Direct Effect)          0.26947              0.26874    1.00000         0.00000
# Static Intervention (Direct Effect)          NA              0.24935    0.92784        -0.01939
# Dynamic Intervention (Direct Effect)         NA              0.34078    1.26807         0.07204
# Threshold Intervention (Direct Effect)       NA              0.05389    0.20053        -0.21485

# pooled HB ICE

fit_hazard_pool <- ice(data = data,
                       K = 4,
                       id = "id",
                       time_name = "t0",
                       outcome_name = "Y",
                       censor_name = "C",
                       competing_name = "D",
                       estimator = pool(hazard = T),
                       comp_effect = 0,
                       outcome_model = Y ~ L1 + L2 + A1 + A2,
                       competing_model = D ~ L1 + L2 + A1 + A2,
                       censor_model = C ~ L1 + L2 + A1 + A2,
                       ref_idx = 0,
                       int_descript = c("Static Intervention",
                                        "Dynamic Intervention"),
                       intervention1.A1 = list(static(3)),
                       intervention1.A2 = list(static(1)),
                       intervention2.A1 = list(dynamic_cat),
                       intervention2.A2 = list(dynamic("compare", "L1", "=", 0))
)

fit_hazard_pool$summary



# Stratified Classical ICE

# new intervention format

fit_classical_strat <- ice(data = data,
                           K = 4,
                           id = "id",
                           time_name = "t0",
                           outcome_name = "Y",
                           censor_name = "C",
                           competing_name = "D",
                           estimator = strat(hazard = F),
                           comp_effect = 0,
                           outcome_model = Y ~ L1 + L2,
                           censor_model = C ~ L1 + L2,
                           ref_idx = 0,
                           int_descript = c("Static Intervention",
                                            "Dynamic Intervention"),
                           intervention1.A1 = list(static(3)),
                           intervention1.A2 = list(static(1)),
                           intervention2.A1 = list(dynamic_cat),
                           intervention2.A2 = list(dynamic("compare", "L1", "=", 0))
)

ice_test$summary

# NP Risk Classical Stratified ICE Risk Ratio Risk Difference
# Natural Course (Direct Effect)         0.26947                  0.26967    1.00000         0.00000
# Static Intervention (Direct Effect)         NA                  0.22639    0.83951        -0.04328
# Dynamic Intervention (Direct Effect)        NA                  0.29826    1.10602         0.02859
# Threshold Intervention (Direct Effect)      NA                  0.29928    1.10980         0.02961

# Stratified HB ICE

fit_hazard_strat <- ice(data = data,
                        K = 4,
                        id = "id",
                        time_name = "t0",
                        outcome_name = "Y",
                        censor_name = "C",
                        competing_name = "D",
                        estimator = strat(hazard = T),
                        comp_effect = 0,
                        outcome_model = Y ~ L1 + L2,
                        competing_model = D ~ L1 + L2,
                        censor_model = C ~ L1 + L2,
                        ref_idx = 0,
                        int_descript = c("Static Intervention",
                                         "Dynamic Intervention"),
                        intervention1.A1 = list(static(3)),
                        intervention1.A2 = list(static(1)),
                        intervention2.A1 = list(dynamic_cat),
                        intervention2.A2 = list(dynamic("compare", "L1", "=", 0))
)

ice_test$summary

# NP Risk HB Stratified ICE Risk Ratio Risk Difference
# Natural Course (Direct Effect)         0.26947           0.26936    1.00000         0.00000
# Static Intervention (Direct Effect)         NA           0.22562    0.83763        -0.04373
# Dynamic Intervention (Direct Effect)        NA           0.29913    1.11053         0.02977
# Threshold Intervention (Direct Effect)      NA           0.29969    1.11262         0.03034

# Weighted ICE

# test_data$A3 <- as.numeric(test_data$A1 == 1)
# test_data$A4 <- as.numeric(test_data$A2 == 2)
# test_data$A5 <- as.numeric(test_data$A3 == 3)

# old intervention format

fit_weighted_ice <- ice(data = data,
                        K = 4,
                        id = "id",
                        time_name = "t0",
                        outcome_name = "Y",
                        censor_name = "C",
                        competing_name = "D",
                        estimator = weight(list(A1 ~ L1 + L2, A2 ~ L1 + L2)),
                        comp_effect = 0,
                        outcome_model = Y ~ L1 + L2,
                        competing_model = D ~ L1 + L2,
                        censor_model = C ~ L1 + L2,
                        ref_idx = 0,
                        int_descript = c("Static Intervention",
                                         "Dynamic Intervention"),
                        intervention1.A1 = list(static(3)),
                        intervention1.A2 = list(static(1)),
                        intervention2.A1 = list(dynamic_cat),
                        intervention2.A2 = list(dynamic("compare", "L1", "=", 0))
)

ice_test$summary


plot_risk(fit_classical_pool,
          fit_classical_strat,
          fit_hazard_pool,
          fit_hazard_strat,
          fit_weighted_ice)

summary_table(fit_classical_pool,
              fit_classical_strat,
              fit_hazard_pool,
              fit_hazard_strat,
              fit_weighted_ice)

static(1)

static(0)

dynamic("compare", var, direction, value)

dynamic("compare", var = "L1", direction = ">", value = 1)

dynamic("absorbing", var, direction, value)

dynamic("absorbing", var = "L1", direction = ">", value = 1)

threshold(lower_bound, upper_bound)

threshold(lower_bound = 1, upper_bound = Inf)

grace_period("uniform", nperiod, var, value)

grace_period("uniform", nperiod = 2, var = "L1", value = 1)

