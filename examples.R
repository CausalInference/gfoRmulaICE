library(gfoRmulaICE)
library(tidyverse)

## Below are some examples

data <- gfoRmulaICE::data

# For the following examples, we consider two interventions.
# Intervention 1 is a static intervention, and
# intervention 2 is a dynamic intervention, where the intervention
# for A1 is user-defined.

dynamic_cat <- case_when(data$L2 < 0 ~ 1,
                         data$L2 >= 0 & data$L2 < 2 ~ 2,
                         T ~ 3)

# 1. Classical Pooled ICE:
# Competing event as direct effect for all interventions;
# Bootstrap with 1000 samples using empirical quantile;
# Natural course as the reference intervention.

ice_pool_classic <- ice(data = data, time_points = 4, id = "id", time_name = "t0",
                        censor_name = "C", outcome_name = "Y",
                        compevent_name = "D",
                        comp_effect = 0,
                        outcome_model = Y ~ L1 + L2 + A1 + A2, 
                        censor_model = C ~ L1 + L2 + A1 + A2,
                        competing_model = D ~ L1 + L2 + A1 + A2,
                        ref_idx = 0,
                        estimator = pool(hazard = F),
                        int_descript = c("Static Intervention",
                                         "Dynamic Intervention"),
                        intervention1.A1 = list(static(3)),
                        intervention1.A2 = list(static(1)),
                        intervention2.A1 = list(dynamic_cat),
                        intervention2.A2 = list(dynamic("compare", "L1", "=", 0)),
                        nsamples = 1000, ci_method = "percentile"
)

summary_table(ice_pool_classic)
plot_risk(ice_pool_classic)

# 2. Classical Stratified ICE:
# Competing event as direct effect for all interventions;
# Bootstrap with 1000 samples using normal quantile;
# Natural course as the reference intervention.

ice_strat_classic <- ice(data = data, time_points = 4, id = "id", time_name = "t0",
                        censor_name = "C", outcome_name = "Y",
                        compevent_name = "D",
                        comp_effect = 0,
                        outcome_model = Y ~ L1 + L2, 
                        censor_model = C ~ L1 + L2,
                        competing_model = D ~ L1 + L2,
                        ref_idx = 0,
                        estimator = strat(hazard = F),
                        int_descript = c("Static Intervention",
                                         "Dynamic Intervention"),
                        intervention1.A1 = list(static(3)),
                        intervention1.A2 = list(static(1)),
                        intervention2.A1 = list(dynamic_cat),
                        intervention2.A2 = list(dynamic("compare", "L1", "=", 0)),
                        nsamples = 1000, ci_method = "normal"
)

summary_table(ice_strat_classic)
plot_risk(ice_strat_classic)

# 3. Hazard based pooled ICE,
# Competing event as direct effect for all interventions;
# Always treat as the reference intervention.
ice_pool_haz <- ice(data = data, time_points = 4, id = "id", time_name = "t0",
                    censor_name = "C", outcome_name = "Y",
                    compevent_name = "D",
                    comp_effect = 0,
                    outcome_model = Y ~ L1 + L2 + A1 + A2, 
                    censor_model = C ~ L1 + L2 + A1 + A2,
                    competing_model = D ~ L1 + L2 + A1 + A2,
                    ref_idx = 0,
                    estimator = pool(hazard = T),
                    int_descript = c("Static Intervention",
                                     "Dynamic Intervention"),
                    intervention1.A1 = list(static(3)),
                    intervention1.A2 = list(static(1)),
                    intervention2.A1 = list(dynamic_cat),
                    intervention2.A2 = list(dynamic("compare", "L1", "=", 0))
)

summary_table(ice_pool_haz)
plot_risk(ice_pool_haz)

# 4a. Hazard based stratified ICE:
# Competing event as total effect for all interventions;
# Y ~ L1 + L2 as outcome model for both intervention 1 and intervention 2;
# D ~ L1 + L2 as competing model for both intervention 1 and intervention 2;
# Natural course as the reference intervention

ice_strat_haz <- ice(data = data, time_points = 4, id = "id", time_name = "t0",
                     censor_name = "C", outcome_name = "Y",
                     compevent_name = "D",
                     outcome_model = Y ~ L1 + L2, 
                     censor_model = C ~ L1 + L2,
                     competing_model = D ~ L1 + L2,
                     comp_effect = 0,
                     ref_idx = 0,
                     estimator = strat(hazard = T),
                     int_descript = c("Static Intervention",
                                      "Dynamic Intervention"),
                     intervention1.A1 = list(static(3)),
                     intervention1.A2 = list(static(1)),
                     intervention2.A1 = list(dynamic_cat),
                     intervention2.A2 = list(dynamic("compare", "L1", "=", 0))
)

summary_table(ice_strat_haz)
plot_risk(ice_strat_haz)

# 4b. Hazard based stratified ICE:
# Competing event as total effect for all interventions;
# Y ~ L1 + L2 as outcome model for intervention 1;
# Y ~ L1 as outcome model for intervention 2;
# D ~ L1 as competing model for intervention 1;
# D ~ L1 + L2 as competing model for intervention 2;
# Natural course as the reference intervention

ice_strat_haz <- ice(data = data, time_points = 4, id = "id", time_name = "t0",
                     censor_name = "C", outcome_name = "Y",
                     compevent_name = "D",
                     outcome_model = Y ~ L1, 
                     censor_model = C ~ L1,
                     competing_model = D ~ L1,
                     comp_effect = 1,
                     ref_idx = 0,
                     estimator = strat(hazard = T),
                     int_descript = c("Static Intervention",
                                      "Dynamic Intervention"),
                     intervention1.A1 = list(static(3)),
                     intervention1.A2 = list(static(1)),
                     intervention2.A1 = list(dynamic_cat),
                     intervention2.A2 = list(dynamic("compare", "L1", "=", 0)),
                     outcomeModel.1 = Y ~ L1 + L2,
                     compModel.2 = D ~ L1 + L2
)

summary_table(ice_strat_haz)
plot_risk(ice_strat_haz)


# 5a. Weighted ICE:
# Competing event as direct effect for all interventions;
# Y ~ L1 as outcome model for both intervention 1 and intervention 2;
# D ~ L1 as competing model for both intervention 1 and intervention 2;
# Natural course as the reference intervention.

ice_weight <- ice(data = data, time_points = 4, id = "id", time_name = "t0",
                  censor_name = "C", outcome_name = "Y",
                  compevent_name = "D",
                  comp_effect = 0,
                  outcome_model = Y ~ L1, 
                  censor_model = C ~ L1,
                  competing_model = D ~ L1,
                  ref_idx = 0,
                  estimator = weight(list(A1 ~ L1 + L2, A2 ~ L1 + L2)),
                  int_descript = c("Static Intervention",
                                   "Dynamic Intervention"),
                  intervention1.A1 = list(static(3)),
                  intervention1.A2 = list(static(1)),
                  intervention2.A1 = list(dynamic_cat),
                  intervention2.A2 = list(dynamic("compare", "L1", "=", 0))
)

summary_table(ice_weight)
plot_risk(ice_weight)

# 5b. Weighted ICE:
# Competing event as direct effect for all interventions;
# Y ~ L1 as outcome model for intervention 1;
# Y ~ L1 + L2 as outcome model for intervention 2;
# D ~ L1 + L2 as competing model for intervention 1;
# D ~ L1 as competing model for intervention 2;
# Bootstrap with 1000 samples using normal quantile;
# Natural course as the reference intervention.
ice_weight <- ice(data = data, time_points = 4, id = "id", time_name = "t0",
                  censor_name = "C", outcome_name = "Y",
                  compevent_name = "D",
                  comp_effect = 0,
                  outcome_model = Y ~ L1, censor_model = C ~ L1,
                  competing_model = D ~ L1,
                  ref_idx = 0,
                  estimator = weight(list(A1 ~ L1 + L2, A2 ~ L1 + L2)),
                  int_descript = c("Static Intervention",
                                   "Dynamic Intervention"),
                  intervention1.A1 = list(static(3)),
                  intervention1.A2 = list(static(1)),
                  intervention2.A1 = list(dynamic_cat),
                  intervention2.A2 = list(dynamic("compare", "L1", "=", 0)),
                  outcomeModel.2 = Y ~ L1 + L2,
                  compModel.1 = D ~ L1 + L2
)
summary_table(ice_weight)
plot_risk(ice_weight)

# For the following example, we consider two interventions:
# Intervention 1: threshold intervention - when the natural value of L2 at time t is lower than -3, set its value to -3.
# Intervention 2: dynamic intervention on L1 (treat when L1 = 0) with uniform grace period of 2 periods

# 6. Classical pooled ICE:
# Competing event as direct effect for all interventions,
# Natural course as the reference intervention.

ice_pool_grace_period <- ice(data = data, time_points = 4, id = "id", time_name = "t0",
                             censor_name = "C", outcome_name = "Y",
                             compevent_name = "D",
                             comp_effect = 0,
                             outcome_model = Y ~ L1 + L2 + A1 + A2, 
                             censor_model = C ~ L1 + L2 + A1 + A2,
                             competing_model = D ~ L1 + L2 + A1 + A2,
                             ref_idx = 0,
                             estimator = pool(hazard = F),
                             int_descript = c("Grace Period", "Threshold Intervention"),
                             intervention1.A2 = list(grace_period("uniform", 2, "L1", 0)),
                             intervention2.L2 = list(threshold(-3, Inf))
)

summary_table(ice_pool_grace_period)
plot_risk(ice_pool_grace_period)

