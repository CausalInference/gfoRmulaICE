# load(file = "../../data/compData.rda")
# library(tidyverse)
# library(stringr)
# library(nnet)
# library(doParallel)
data <- gfoRmulaICE::compData
library(data.table)
library(splines)
library(Hmisc)
set.seed(1)

test_that(
  "check classical pooled ICE direct effect - dynamic interventions",
  {
    ice_fit1 <- ice(data = data, time_points = 4,
                    id = "id", time_name = "t0",
                    censor_name = "C", outcome_name = "Y",
                    compevent_name = "D",
                    comp_effect = 0,
                    outcome_model = Y ~ L1 + L2 + A1 + A2,
                    censor_model = C ~ L1 + L2 + A1 + A2,
                    ref_idx = 0,
                    estimator = pool(hazard = F),
                    int_descript = c("Dynamic Intervention 1", "Dynamic Intervention 2",
                                     "Dynamic Intervention 3"),
                    nsamples = 10, ci_method = "percentile", parallel = T, ncores = 5,
                    intervention1.A2 = list(dynamic("L1 == 0", static(0), static(1))),
                    intervention2.A2 = list(dynamic("L1 == 0", static(0), static(1), absorb = T)),
                    intervention3.A2 = list(dynamic("L1 == 0", static(0), natural_course()))
    )
    out <- as.matrix(ice_fit1$summary)
    print(out)
    rownames(out) <- c()
    out <- as.vector(out)
    expect_equal(
      out,
      c(0.26851, NA, NA, NA, 0.26874, 0.30568, 0.25787, 0.32326,
        1.00000, 1.13746, 0.95957, 1.20286, 0.00000, 0.03694, -0.01087, 0.05452,
        0.00358, 0.00639, 0.00316, 0.00645, 1.00000, 0.01299, 0.00295, 0.01318,
        0.00000, 0.00371, 0.00087, 0.00385, 0.26485, 0.29583, 0.25506, 0.31347,
        0.27604, 0.31498, 0.26515, 0.33242, 1.00000, 1.11652, 0.95484, 1.18311,
        1.00000, 1.15747, 0.96365, 1.22291, 0.00000, 0.03088, -0.01235, 0.04852,
        0.00000, 0.04248, -0.00963, 0.06013)
    )})

test_that(
  "check classical pooled ICE direct effect - dynamic interventions (using data table)",
  {
    ice_fit1 <- ice(data = data.table(data), time_points = 4,
                    id = "id", time_name = "t0",
                    censor_name = "C", outcome_name = "Y",
                    compevent_name = "D",
                    comp_effect = 0,
                    outcome_model = Y ~ L1 + L2 + A1 + A2,
                    censor_model = C ~ L1 + L2 + A1 + A2,
                    ref_idx = 0,
                    estimator = pool(hazard = F),
                    int_descript = c("Dynamic Intervention 1", "Dynamic Intervention 2",
                                     "Dynamic Intervention 3"),
                    nsamples = 10, ci_method = "percentile", parallel = T, ncores = 3,
                    intervention1.A2 = list(dynamic("L1 == 0", static(0), static(1))),
                    intervention2.A2 = list(dynamic("L1 == 0", static(0), static(1), absorb = T)),
                    intervention3.A2 = list(dynamic("L1 == 0", static(0), natural_course()))
    )
    out <- as.matrix(ice_fit1$summary)
    print(out)
    rownames(out) <- c()
    out <- as.vector(out)
    expect_equal(
      out,
      c(0.26851, NA, NA, NA, 0.26874, 0.30568, 0.25787, 0.32326,
        1.00000, 1.13746, 0.95957, 1.20286, 0.00000, 0.03694, -0.01087, 0.05452,
        0.00358, 0.00639, 0.00316, 0.00645, 1.00000, 0.01299, 0.00295, 0.01318,
        0.00000, 0.00371, 0.00087, 0.00385, 0.26485, 0.29583, 0.25506, 0.31347,
        0.27604, 0.31498, 0.26515, 0.33242, 1.00000, 1.11652, 0.95484, 1.18311,
        1.00000, 1.15747, 0.96365, 1.22291, 0.00000, 0.03088, -0.01235, 0.04852,
        0.00000, 0.04248, -0.00963, 0.06013)
    )})

test_that(
  "check different interventions",
  {
    ice_fit2 <- ice(data = data, time_points = 4,
                    id = "id", time_name = "t0",
                    censor_name = "C", outcome_name = "Y",
                    compevent_name = "D",
                    comp_effect = 0,
                    outcome_model = Y ~ L1 + L2 + A1 + A2,
                    censor_model = C ~ L1 + L2 + A1 + A2,
                    ref_idx = 0,
                    estimator = pool(hazard = F),
                    nsamples = 10, ci_method = "percentile", parallel = T, ncores = 5,
                    int_descript = c("Static Intervention", "Threshold Intervention",
                                     "Dynamic Intervention with Grace Period"),
                    intervention1.A1 = list(static(3)),
                    intervention1.A2 = list(static(1)),
                    intervention2.L2 = list(threshold(-3, Inf)),
                    intervention3.A2 = list(grace_period("uniform", 2, "L1 == 0"))
    )
    out <- as.matrix(ice_fit2$summary)
    print(out)
    rownames(out) <- c()
    out <- as.vector(out)
    expect_equal(
      out,
      c(0.26851, NA, NA, NA, 0.26874, 0.24935, 0.26874, 0.28722, 1.00000, 0.92784, 1.00000, 1.06876,
        0.000000, -0.01939, 0.00000, 0.01848, 0.00358, 0.00992, 0.00358, 0.00350,
        1.00000, 0.03318, 0.00000, 0.00837, 0.00000, 0.00894, 0.00000, 0.00212,
        0.26485, 0.23801, 0.26485, 0.28350, 0.27604, 0.26802, 0.27604, 0.29357,
        1.00000, 0.88547, 1.00000, 1.05415, 1.00000, 0.98762, 1.00000, 1.08150,
        0.00000, -0.03078, 0.00000, 0.01469, 0.00000, -0.00336, 0.00000, 0.02150)
    )})

test_that(
  "check user-defined intervention",
  {
    ice_fit3 <- ice(data = data, time_points = 4,
                    id = "id", time_name = "t0",
                    censor_name = "C", outcome_name = "Y",
                    compevent_name = "D",
                    comp_effect = 0,
                    outcome_model = Y ~ L1 + L2 + A1 + A2,
                    censor_model = C ~ L1 + L2 + A1 + A2,
                    ref_idx = 0,
                    estimator = pool(hazard = F),
                    nsamples = 10, ci_method = "percentile", parallel = T, ncores = 5,
                    int_descript = c("Static Intervention", "Dynamic Intervention"),
                    intervention1.A1 = list(static(3)),
                    intervention1.A2 = list(static(1)),
                    intervention2.A1 = list(case_when(data$L2 < 0 ~ 1,
                                                      data$L2 >= 0 & data$L2 < 2 ~ 2, T ~ 3)),
                    intervention2.A2 = list(dynamic("L1 == 0", static(0), static(1)))
    )
    out <- as.matrix(ice_fit3$summary)
    print(out)
    rownames(out) <- c()
    out <- as.vector(out)
    expect_equal(
      out,
      c(0.26851, NA, NA, 0.26874, 0.24935, 0.34078, 1.00000, 0.92784, 1.26807,
        0.00000, -0.01939, 0.07204, 0.00358, 0.00992, 0.00854, 1.00000, 0.03318, 0.01820,
        0.00000, 0.00894, 0.00553, 0.26485, 0.23801, 0.32746, 0.27604, 0.26802, 0.35214,
        1.00000, 0.88547, 1.23583, 1.00000, 0.98762, 1.28880, 0.00000, -0.03078, 0.06251,
        0.00000, -0.00336, 0.07789)
    )})

test_that(
  "check hazard-based pooled ICE - time-specific hazard model same as the outcome model",
  {

    ice_fit4a <- ice(data = data, time_points = 4,
                     id = "id", time_name = "t0",
                     censor_name = "C", outcome_name = "Y",
                     compevent_name = "D",
                     comp_effect = 0,
                     outcome_model = Y ~ L1 + L2 + A1 + A2,
                     censor_model = C ~ L1 + L2 + A1 + A2,
                     competing_model = D ~ L1 + L2 + A1 + A2,
                     ref_idx = 0,
                     estimator = pool(hazard = T),
                     nsamples = 10, ci_method = "percentile", parallel = T, ncores = 5,
                     int_descript = c("Static Intervention", "Dynamic Intervention"),
                     intervention1.A1 = list(static(3)),
                     intervention1.A2 = list(static(1)),
                     intervention2.A1 = list(case_when(data$L2 < 0 ~ 1,
                                                       data$L2 >= 0 & data$L2 < 2 ~ 2, T ~ 3)),
                     intervention2.A2 = list(dynamic("L1 == 0", static(0), static(1)))
    )

    out <- as.matrix(ice_fit4a$summary)
    print(out)
    rownames(out) <- c()
    out <- as.vector(out)
    expect_equal(
      out,
      c(0.26851, NA, NA, 0.21482, 0.27801, 0.23550, 1.00000, 1.29413, 1.09626,
        0.00000, 0.06319, 0.02068, 0.00651, 0.01042, 0.01173, 1.00000, 0.04093, 0.02693,
        0.00000, 0.00847, 0.00632, 0.20731, 0.26289, 0.22275, 0.22439, 0.29616, 0.25829,
        1.00000, 1.24413, 1.06832, 1.00000, 1.35151, 1.15315, 0.00000, 0.05192, 0.01443,
        0.00000, 0.07481, 0.03421)
    )})

test_that(
  "check hazard-based pooled ICE - time-specific hazard model Y ~ L1 + L2",
  {

    ice_fit4b <- ice(data = data, time_points = 4,
                     id = "id", time_name = "t0",
                     censor_name = "C", outcome_name = "Y",
                     compevent_name = "D",
                     comp_effect = 0,
                     outcome_model = Y ~ L1 + L2 + A1 + A2,
                     censor_model = C ~ L1 + L2 + A1 + A2,
                     competing_model = D ~ L1 + L2 + A1 + A2,
                     hazard_model = Y ~ L1 + L2,
                     ref_idx = 0,
                     estimator = pool(hazard = T),
                     nsamples = 10, ci_method = "percentile", parallel = T, ncores = 5,
                     int_descript = c("Static Intervention", "Dynamic Intervention"),
                     intervention1.A1 = list(static(3)),
                     intervention1.A2 = list(static(1)),
                     intervention2.A1 = list(case_when(data$L2 < 0 ~ 1,
                                                       data$L2 >= 0 & data$L2 < 2 ~ 2, T ~ 3)),
                     intervention2.A2 = list(dynamic("L1 == 0", static(0), static(1)))
    )

    out <- as.matrix(ice_fit4b$summary)
    print(out)
    rownames(out) <- c()
    out <- as.vector(out)
    expect_equal(
      out,
      c(0.26851, NA, NA, 0.21573, 0.23591, 0.23651, 1.00000, 1.09353, 1.09630, 0.00000, 0.02018, 0.02077,
        0.00735, 0.00956, 0.00780, 1.00000, 0.0229, 0.0037, 0.00000, 0.00506, 0.00089, 0.20279, 0.21755, 0.21164,
        0.22491, 0.24662, 0.23498, 1.00000, 1.05853, 1.04028, 1.00000, 1.12489, 1.05111,
        0.00000, 0.01285, 0.00859, 0.00000, 0.02722, 0.01105)
    )})

test_that(
  "check hazard-based pooled ICE - pooled-over-time hazard model",
  {

    ice_fit4c <- ice(data = data, time_points = 4,
                     id = "id", time_name = "t0",
                     censor_name = "C", outcome_name = "Y",
                     compevent_name = "D",
                     comp_effect = 0,
                     outcome_model = Y ~ L1 + L2 + A1 + A2,
                     censor_model = C ~ L1 + L2 + A1 + A2,
                     competing_model = D ~ L1 + L2 + A1 + A2,
                     hazard_model = Y ~ L1 + L2 + A1 + A2 + ns(t0, df = 2),
                     global_hazard = T,
                     ref_idx = 0,
                     estimator = pool(hazard = T),
                     nsamples = 10, ci_method = "percentile", parallel = T, ncores = 5,
                     int_descript = c("Static Intervention", "Dynamic Intervention"),
                     intervention1.A1 = list(static(3)),
                     intervention1.A2 = list(static(1)),
                     intervention2.A1 = list(case_when(data$L2 < 0 ~ 1,
                                                       data$L2 >= 0 & data$L2 < 2 ~ 2, T ~ 3)),
                     intervention2.A2 = list(dynamic("L1 == 0", static(0), static(1)))
    )

    out <- as.matrix(ice_fit4c$summary)
    print(out)
    rownames(out) <- c()
    out <- as.vector(out)
    expect_equal(
      out,
      c(0.26851, NA, NA, 0.26811, 0.23677, 0.35556, 1.00000, 0.88312, 1.32618, 0.00000, -0.03134, 0.08745,
        0.00475, 0.00584, 0.00717, 1.00000, 0.02224, 0.00786, 0.00000, 0.00636, 0.00291, 0.26576, 0.23028, 0.33622,
        0.28035, 0.24731, 0.35660, 1.00000, 0.84884, 1.26230, 1.00000, 0.90787, 1.28456,
        0.00000, -0.04225, 0.07013, 0.00000, -0.02510, 0.07872)
    )})

test_that(
  "check classical stratified ICE",
  {

    ice_fit4d <- ice(data = data, time_points = 4,
                     id = "id", time_name = "t0",
                     censor_name = "C", outcome_name = "Y",
                     compevent_name = "D",
                     comp_effect = 0,
                     outcome_model = Y ~ L1 + L2,
                     censor_model = C ~ L1 + L2,
                     ref_idx = 0,
                     estimator = strat(hazard = F),
                     nsamples = 10, ci_method = "percentile", parallel = T, ncores = 5,
                     int_descript = c("Static Intervention", "Dynamic Intervention"),
                     intervention1.A1 = list(static(3)),
                     intervention1.A2 = list(static(1)),
                     intervention2.A1 = list(case_when(data$L2 < 0 ~ 1,
                                                       data$L2 >= 0 & data$L2 < 2 ~ 2, T ~ 3)),
                     intervention2.A2 = list(dynamic("L1 == 0", static(0), static(1)))
    )

    out <- as.matrix(ice_fit4d$summary)
    print(out)
    rownames(out) <- c()
    out <- as.vector(out)
    expect_equal(
      out,
      as.vector(c(0.26947, NA, NA, 0.26967, 0.22639, 0.29826, 1.00000, 0.83951, 1.10602, 0.00000, -0.04328, 0.02859,
                  0.05596, 0.01566, 0.04314, 1.00000, 0.18866, 0.20060, 0.00000, 0.06707, 0.05424, 0.21892, 0.21681, 0.30288,
                  0.38310, 0.26405, 0.43180, 1.00000, 0.57039, 0.91857, 1.00000, 1.11667, 1.50712,
                  0.00000, -0.16463, -0.02726, 0.00000, 0.02416, 0.13272))
    )})

test_that(
  "check hazard-extended stratified ICE - time-specific hazard model Y ~ L1",
  {

    ice_fit4e <- ice(data = data, time_points = 4,
                     id = "id", time_name = "t0",
                     censor_name = "C", outcome_name = "Y",
                     compevent_name = "D",
                     comp_effect = 0,
                     outcome_model = Y ~ L1 + L2,
                     censor_model = C ~ L1 + L2,
                     competing_model = D ~ L1 + L2,
                     hazard_model = Y ~ L1,
                     ref_idx = 0,
                     estimator = strat(hazard = T),
                     nsamples = 10, ci_method = "percentile", parallel = T, ncores = 5,
                     int_descript = c("Static Intervention", "Dynamic Intervention"),
                     intervention1.A1 = list(static(3)),
                     intervention1.A2 = list(static(1)),
                     intervention2.A1 = list(case_when(data$L2 < 0 ~ 1,
                                                       data$L2 >= 0 & data$L2 < 2 ~ 2, T ~ 3)),
                     intervention2.A2 = list(dynamic("L1 == 0", static(0), static(1)))
    )

    out <- as.matrix(ice_fit4e$summary)
    print(out)
    rownames(out) <- c()
    out <- as.vector(out)
    expect_equal(
      out,
      as.vector(c(0.26947, NA, NA, 0.26941, 0.21893, 0.30104, 1.00000, 0.81264, 1.11744,
                  0.00000, -0.05048, 0.03164, 0.05155, 0.01875, 0.05161, 1.00000, 0.12205, 0.19777,
                  0.00000, 0.05418, 0.05152, 0.25253, 0.20648, 0.29869, 0.40821, 0.26629, 0.44062,
                  1.00000, 0.55239, 0.95591, 1.00000, 0.92452, 1.54891, 0.00000, -0.18326, -0.01451,
                  0.00000, -0.02072, 0.13689))
    )})

test_that(
  "check doubly robust ICE",
  {

    ice_fit4f <- ice(data = data, time_points = 4,
                     id = "id", time_name = "t0",
                     censor_name = "C", outcome_name = "Y",
                     compevent_name = "D",
                     comp_effect = 0,
                     outcome_model = Y ~ L1 + L2,
                     censor_model = C ~ L1 + L2,
                     ref_idx = 0,
                     estimator = weight(list(A1 ~ L1 + L2, A2 ~ L1 + L2)),
                     nsamples = 10, ci_method = "percentile", parallel = T, ncores = 5,
                     int_descript = c("Static Intervention", "Dynamic Intervention"),
                     intervention1.A1 = list(static(3)),
                     intervention1.A2 = list(static(1)),
                     intervention2.A1 = list(case_when(data$L2 < 0 ~ 1,
                                                       data$L2 >= 0 & data$L2 < 2 ~ 2, T ~ 3)),
                     intervention2.A2 = list(dynamic("L1 == 0", static(0), static(1)))
    )

    out <- as.matrix(ice_fit4f$summary)
    print(out)
    rownames(out) <- c()
    out <- as.vector(out)
    expect_equal(
      out,
      as.vector(c(0.26947, NA, NA, 0.26118, 0.22097, 0.30014, 1.00000, 0.84602, 1.14917,
                  0.00000, -0.04022, 0.03896, 0.06102, 0.01655, 0.04078, 1.00000, 0.15695, 0.33725,
                  0.00000, 0.05179, 0.09196, 0.20320, 0.20418, 0.30312, 0.38399, 0.24966, 0.39734,
                  1.00000, 0.63963, 0.84017, 1.00000, 1.12959, 1.72028,
                  0.0000, -0.13911, -0.06411, 0.00000, 0.02622, 0.16387))
    )})

test_that(
  "check hazard-based stratified ICE with intervention-specific models",
  {

    ice_fit4h <- ice(data = data, time_points = 4,
                     id = "id", time_name = "t0",
                     censor_name = "C", outcome_name = "Y",
                     compevent_name = "D",
                     outcome_model = Y ~ L1, censor_model = C ~ L1,
                     competing_model = D ~ L1,
                     comp_effect = 1,
                     ref_idx = 0,
                     estimator = strat(hazard = T),
                     nsamples = 10, ci_method = "normal", parallel = T, ncores = 5,
                     int_descript = c("Static Intervention",
                                      "Dynamic Intervention"),
                     intervention1.A1 = list(static(3)),
                     intervention1.A2 = list(static(1)),
                     intervention2.A1 = list(case_when(data$L2 < 0 ~ 1,
                                                       data$L2 >= 0 & data$L2 < 2 ~ 2, T ~ 3)),
                     intervention2.A2 = list(dynamic("L1 == 0", static(0), static(1))),
                     outcomeModel.1 = Y ~ L1 + L2,
                     compModel.2 = D ~ L1 + L2
    )

    out <- as.matrix(ice_fit4h$summary)
    print(out)
    rownames(out) <- c()
    out <- as.vector(out)
    expect_equal(
      out,
      as.vector(c(0.22987, NA, NA, 0.24359, 0.20125, 0.25957, 1.00000, 0.82617, 1.06557,
                  0.00000, -0.04234, 0.01597, 0.03432, 0.01988, 0.03586, 1.00000, 0.11311, 0.09651,
                  0.00000, 0.03653, 0.02823, 0.17632, 0.16228, 0.18929, 0.31087, 0.24022, 0.32985,
                  1.00000, 0.60449, 0.87641, 1.00000, 1.04786, 1.25472, 0.00000, -0.11395, -0.03935,
                  0.00000, 0.02926, 0.07129))
    )})

test_that(
  "check flexible model specification - complicated terms",
  {

    ice_fit5a <- ice(data = data, time_points = 4,
                     id = "id", time_name = "t0",
                     censor_name = "C", outcome_name = "Y",
                     compevent_name = "D",
                     comp_effect = 0,
                     outcome_model = Y ~ I(L1^2) + rcspline.eval(lag1_L2, knots = 1:3) + A1 + A2,
                     censor_model = C ~ lag1_L1 + poly(L2, degree = 2) + A1 + A2,
                     ref_idx = 0,
                     estimator = pool(hazard = F),
                     nsamples = 10, ci_method = "percentile", parallel = T, ncores = 5,
                     int_descript = c("Static Intervention", "Dynamic Intervention"),
                     intervention1.A1 = list(static(3)),
                     intervention1.A2 = list(static(1)),
                     intervention2.A1 = list(case_when(data$L2 < 0 ~ 1,
                                                       data$L2 >= 0 & data$L2 < 2 ~ 2, T ~ 3)),
                     intervention2.A2 = list(dynamic("L1 == 0", static(0), static(1)))
    )

    out <- as.matrix(ice_fit5a$summary)
    print(out)
    rownames(out) <- c()
    out <- as.vector(out)
    expect_equal(
      out,
      as.vector(c(0.26831, NA, NA, 0.26654, 0.19984, 0.34456, 1.00000, 0.74976, 1.29275,
                  0.00000, -0.06670, 0.07803, 0.00433, 0.00781, 0.00671, 1.00000, 0.0236, 0.0072,
                  0.00000, 0.00621, 0.00273, 0.26286, 0.19210, 0.33251, 0.27645, 0.21332, 0.35273,
                  1.00000, 0.71545, 1.26322, 1.00000, 0.78523, 1.28678,
                  0.00000, -0.07641, 0.06945, 0.00000, -0.05835, 0.07809))
    )})

test_that(
  "check flexible model specification - using static intervention as reference",
  {

    ice_fit5b <- ice(data = data, time_points = 4,
                     id = "id", time_name = "t0",
                     censor_name = "C", outcome_name = "Y",
                     compevent_name = "D",
                     comp_effect = 0,
                     outcome_model = Y ~ I(L1^2) + rcspline.eval(lag1_L2, knots = 1:3) + A1 + A2,
                     censor_model = C ~ lag1_L1 + poly(L2, degree = 2) + A1 + A2,
                     ref_idx = 1,
                     estimator = pool(hazard = F),
                     nsamples = 10, ci_method = "percentile", parallel = T, ncores = 5,
                     int_descript = c("Static Intervention", "Dynamic Intervention"),
                     intervention1.A1 = list(static(3)),
                     intervention1.A2 = list(static(1)),
                     intervention2.A1 = list(case_when(data$L2 < 0 ~ 1,
                                                       data$L2 >= 0 & data$L2 < 2 ~ 2, T ~ 3)),
                     intervention2.A2 = list(dynamic("L1 == 0", static(0), static(1)))
    )

    out <- as.matrix(ice_fit5b$summary)
    print(out)
    rownames(out) <- c()
    out <- as.vector(out)
    expect_equal(
      out,
      as.vector(c(0.26831, NA, NA, 0.26654, 0.19984, 0.34456, 1.33376, 1.00000, 1.72422, 0.06670, 0.00000, 0.14473,
                  0.00433, 0.00781, 0.00671, 0.04202, 1.00000, 0.05251, 0.00621, 0.00000, 0.00644,
                  0.26286, 0.19564, 0.33251, 0.27645, 0.20504, 0.35273, 1.27352, 1.00000, 1.63394,
                  1.39792, 1.00000, 1.78442, 0.05835, 0.00000, 0.13287, 0.07641, 0.00000, 0.15065))
    )})

test_that(
  "check complicated scenario 1",
  {

    ice_fit6a <- ice(data = data, time_points = 4,
                     id = "id", time_name = "t0",
                     censor_name = "C", outcome_name = "Y",
                     compevent_name = "D",
                     comp_effect = 0,
                     outcome_model = Y ~ I(L1^2) + I(L2^3) + A1 + lag1_A2,
                     censor_model = C ~ L1 + rcspline.eval(lag1_L2, knots = 1:3) + lag2_A1 + A2,
                     competing_model = D ~ lag3_L1 + poly(L2, degree = 2) + A1 + A2,
                     hazard_model = Y ~ L1 + ns(lag1_L2, df = 3) + A1 + A2 + ns(t0, df = 2),
                     global_hazard = T,
                     ref_idx = 0,
                     estimator = pool(hazard = T),
                     nsamples = 10, ci_method = "percentile", parallel = T, ncores = 5,
                     int_descript = c("Static Intervention", "Dynamic Intervention"),
                     intervention1.A1 = list(static(3)),
                     intervention1.A2 = list(static(1)),
                     intervention2.A1 = list(case_when(data$L2 < 0 ~ 1,
                                                       data$L2 >= 0 & data$L2 < 2 ~ 2, T ~ 3)),
                     intervention2.A2 = list(dynamic("L1 == 0", static(0), static(1)))
    )
    out <- as.matrix(ice_fit6a$summary)
    rownames(out) <- c()
    print(out)
    out <- as.vector(out)
    expect_equal(
      out,
      as.vector(c(0.27178, NA, NA, 0.26748, 0.22316, 0.34947, 1.00000, 0.83431, 1.30650, 0.00000, -0.04432, 0.08199,
                  0.00497, 0.00518, 0.00621, 1.00000, 0.01790, 0.00738, 0.00000, 0.00518, 0.00214,
                  0.26197, 0.21431, 0.33886, 0.27808, 0.22955, 0.35855, 1.00000, 0.80161, 1.28506,
                  1.00000, 0.85402, 1.30689, 0.00000, -0.05395, 0.07621, 0.00000, -0.03916, 0.08246))
    )})

test_that(
  "check complicated scenario 2 - doubly robust ICE with intervention-specific models",
  {
    ice_fit6b <- ice(data = data, time_points = 4,
                     id = "id", time_name = "t0",
                     censor_name = "C", outcome_name = "Y",
                     compevent_name = "D",
                     outcome_model = Y ~ I(L1^2), censor_model = C ~ lag1_L1,
                     competing_model = D ~ L1,
                     comp_effect = 1,
                     ref_idx = 0,
                     estimator = weight(list(A1 ~ L1 + I(L2^2) + lag1_L2, A2 ~ lag2_L1 + L1 + ns(L2, df = 2))),
                     int_descript = c("Static Intervention",
                                      "Dynamic Intervention"),
                     intervention1.A1 = list(static(3)),
                     intervention1.A2 = list(static(1)),
                     intervention2.A1 = list(case_when(data$L2 < 0 ~ 1,
                                                       data$L2 >= 0 & data$L2 < 2 ~ 2, T ~ 3)),
                     intervention2.A2 = list(dynamic("L1 == 0", static(0), static(1))),
                     outcomeModel.1 = Y ~ I(L1^2) + ns(lag1_L2, df = 3),
                     compModel.2 = D ~ lag1_L1 + ns(L2, df = 2)
    )

    out <- as.matrix(ice_fit6b$summary)
    print(out)
    rownames(out) <- c()
    out <- as.vector(out)
    expect_equal(
      out,
      as.vector(c(0.22985, NA, NA, 0.22724, 0.19601, 0.24222,
                  1.00000, 0.86256, 1.06592, 0.00000, -0.03123, 0.01498))
    )})

test_that(
  "check complicated scenario 3 - intervention-specific time options",
  {
    
    ice_fit7a <- ice(data = data, time_points = 4, 
                    id = "id", time_name = "t0",
                    censor_name = "C", outcome_name = "Y",
                    compevent_name = "D",
                    comp_effect = 0,
                    outcome_model = Y ~ L1 + L2 + A1 + A2, 
                    censor_model = C ~ L1 + L2 + A1 + A2,
                    ref_idx = 0,
                    estimator = pool(hazard = F),
                    int_descript = c("Static Intervention", "Dynamic Intervention"),
                    intervention1.A1 = list(static(3), 0:2),
                    intervention1.A2 = list(static(1), 1:3),
                    intervention2.A1 = list(case_when(data$L2 < 0 ~ 1,
                                                      data$L2 >= 0 & data$L2 < 2 ~ 2, T ~ 3), 1:2),
                    intervention2.A2 = list(dynamic("L1 == 0", static(0), static(1)))
    )
    
    ice_fit7b <- ice(data = data, time_points = 4, 
                     id = "id", time_name = "t0",
                     censor_name = "C", outcome_name = "Y",
                     compevent_name = "D",
                     comp_effect = 0,
                     outcome_model = Y ~ L1 + L2 + A1 + A2, 
                     censor_model = C ~ L1 + L2 + A1 + A2,
                     competing_model = D ~ L1 + L2 + A1 + A2,
                     ref_idx = 0,
                     estimator = pool(hazard = T),
                     int_descript = c("Static Intervention", "Dynamic Intervention"),
                     intervention1.A1 = list(static(3), 0:1),
                     intervention1.A2 = list(static(1), 1:2),
                     intervention2.A1 = list(case_when(data$L2 < 0 ~ 1,
                                                       data$L2 >= 0 & data$L2 < 2 ~ 2, T ~ 3), 2:3),
                     intervention2.A2 = list(dynamic("L1 == 0", static(0), static(1)), 0:2)
    )
    
    ice_fit7c <- ice(data = data, time_points = 4, 
                     id = "id", time_name = "t0",
                     censor_name = "C", outcome_name = "Y",
                     compevent_name = "D",
                     comp_effect = 0,
                     outcome_model = Y ~ L1 + L2 + A1 + A2, 
                     censor_model = C ~ L1 + L2 + A1 + A2,
                     competing_model = D ~ L1 + L2 + A1 + A2,
                     hazard_model = Y ~ L1 + L2,
                     ref_idx = 0,
                     estimator = pool(hazard = T),
                     int_descript = c("Static Intervention", "Dynamic Intervention"),
                     intervention1.A1 = list(static(3), 0:1),
                     intervention1.A2 = list(static(1), 1:2),
                     intervention2.A1 = list(case_when(data$L2 < 0 ~ 1,
                                                       data$L2 >= 0 & data$L2 < 2 ~ 2, T ~ 3), 2:3),
                     intervention2.A2 = list(dynamic("L1 == 0", static(0), static(1)), 0:2)
    )
    
    ice_fit7d <- ice(data = data, time_points = 4, 
                     id = "id", time_name = "t0",
                     censor_name = "C", outcome_name = "Y",
                     compevent_name = "D",
                     comp_effect = 0,
                     outcome_model = Y ~ L1 + L2 + A1 + A2, 
                     censor_model = C ~ L1 + L2 + A1 + A2,
                     competing_model = D ~ L1 + L2 + A1 + A2,
                     hazard_model = Y ~ L1 + L2 + A1 + A2 + ns(t0, df = 2),
                     global_hazard = T,
                     ref_idx = 0,
                     estimator = pool(hazard = T),
                     int_descript = c("Static Intervention", "Dynamic Intervention"),
                     intervention1.A1 = list(static(3), 0:1),
                     intervention1.A2 = list(static(1), 1:2),
                     intervention2.A1 = list(case_when(data$L2 < 0 ~ 1,
                                                       data$L2 >= 0 & data$L2 < 2 ~ 2, T ~ 3), 2:3),
                     intervention2.A2 = list(dynamic("L1 == 0", static(0), static(1)), 0:2)
    )
    
    ice_fit7e <- ice(data = data, time_points = 4, 
                     id = "id", time_name = "t0",
                     censor_name = "C", outcome_name = "Y",
                     compevent_name = "D",
                     comp_effect = 0,
                     outcome_model = Y ~ L1 + L2, 
                     censor_model = C ~ L1 + L2,
                     ref_idx = 0,
                     estimator = strat(hazard = F),
                     int_descript = c("Static Intervention", "Dynamic Intervention"),
                     intervention1.A1 = list(static(3), 0:1),
                     intervention1.A2 = list(static(1), 1:2),
                     intervention2.A1 = list(case_when(data$L2 < 0 ~ 1,
                                                       data$L2 >= 0 & data$L2 < 2 ~ 2, T ~ 3), 2:3),
                     intervention2.A2 = list(dynamic("L1 == 0", static(0), static(1)), 0:2)
    )
    
    ice_fit7f <- ice(data = data, time_points = 4, 
                     id = "id", time_name = "t0",
                     censor_name = "C", outcome_name = "Y",
                     compevent_name = "D",
                     comp_effect = 0,
                     outcome_model = Y ~ L1 + L2, 
                     censor_model = C ~ L1 + L2,
                     competing_model = D ~ L1 + L2,
                     hazard_model = Y ~ L1,
                     ref_idx = 0,
                     estimator = strat(hazard = T),
                     int_descript = c("Static Intervention", "Dynamic Intervention"),
                     intervention1.A1 = list(static(3), 0:1),
                     intervention1.A2 = list(static(1), 1:2),
                     intervention2.A1 = list(case_when(data$L2 < 0 ~ 1,
                                                       data$L2 >= 0 & data$L2 < 2 ~ 2, T ~ 3), 2:3),
                     intervention2.A2 = list(dynamic("L1 == 0", static(0), static(1)), 0:2)
    )
    
    ice_fit7g <- ice(data = data, time_points = 4,  
                     id = "id", time_name = "t0",
                     censor_name = "C", outcome_name = "Y",
                     compevent_name = "D",
                     comp_effect = 0,
                     outcome_model = Y ~ L1 + L2, 
                     censor_model = C ~ L1 + L2,
                     ref_idx = 0,
                     estimator = weight(list(A1 ~ L1 + L2, A2 ~ L1 + L2)),
                     int_descript = c("Static Intervention", "Dynamic Intervention"),
                     intervention1.A1 = list(static(3), 0:1),
                     intervention1.A2 = list(static(1), 1:2),
                     intervention2.A1 = list(case_when(data$L2 < 0 ~ 1,
                                                       data$L2 >= 0 & data$L2 < 2 ~ 2, T ~ 3), 2:3),
                     intervention2.A2 = list(dynamic("L1 == 0", static(0), static(1)), 0:2)
    )
    
    ice_fit7h <- ice(data = data, time_points = 4, 
                     id = "id", time_name = "t0",
                     censor_name = "C", outcome_name = "Y",
                     compevent_name = "D",
                     outcome_model = Y ~ L1, censor_model = C ~ L1,
                     competing_model = D ~ L1,
                     comp_effect = 1,
                     ref_idx = 0,
                     estimator = strat(hazard = T),
                     int_descript = c("Static Intervention",
                                      "Dynamic Intervention"),
                     intervention1.A1 = list(static(3), 0:1),
                     intervention1.A2 = list(static(1), 1:2),
                     intervention2.A1 = list(case_when(data$L2 < 0 ~ 1,
                                                       data$L2 >= 0 & data$L2 < 2 ~ 2, T ~ 3), 2:3),
                     intervention2.A2 = list(dynamic("L1 == 0", static(0), static(1)), 0:2),
                     outcomeModel.1 = Y ~ L1 + L2,
                     compModel.2 = D ~ L1 + L2
    )
    
    out <- as.matrix(summary_table(ice_fit7a, ice_fit7b, ice_fit7c, ice_fit7d, 
                                   ice_fit7e, ice_fit7f, ice_fit7g, ice_fit7h)[, -(1:2)])
    print(out)
    rownames(out) <- c()
    colnames(out) <- c()
    out <- as.vector(out)
    expect_equal(
      out,
      as.vector(c(0.26851, NA, NA, 0.26851, NA, NA, 0.26851, NA, NA, 0.26851, NA, NA, 
                  0.26947, NA, NA, 0.26947, NA, NA, 0.26947, NA, NA, 0.22987, NA, NA,
                  0.26874, 0.27218, 0.33107, 0.21482, 0.26318, 0.20892, 0.21573, 0.23855, 0.21463, 
                  0.26811, 0.27992, 0.32020, 0.26967, 0.27646, 0.31554, 0.26941, 0.27523, 0.31377,
                  0.26118, 0.27455, 0.31287, 0.24359, 0.24675, 0.27024, 
                  1.00000, 1.01280, 1.23195, 1.00000, 1.22513, 0.97255, 1.00000, 1.10579, 0.99492, 
                  1.00000, 1.04407, 1.19430, 1.00000, 1.02521, 1.17010, 1.00000, 1.02162, 1.16467, 
                  1.00000, 1.05118, 1.19790, 1.00000, 1.01294, 1.10937, 
                  0.00000, 0.00344, 0.06234, 0.00000, 0.04836, -0.00590, 0.00000, 0.02282, -0.00110,
                  0.00000, 0.01181, 0.05209, 0.00000, 0.00680, 0.04587, 0.00000, 0.00582, 0.04436, 
                  0.00000, 0.01337, 0.05169, 0.00000, 0.00315, 0.02664))
    )})

