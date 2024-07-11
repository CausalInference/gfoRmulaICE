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
        0.00366, 0.00544, 0.00392, 0.00559, 1.00000, 0.00794, 0.00390, 0.00915,
        0.00000, 0.00240, 0.00101, 0.00276, 0.26430, 0.29854, 0.25272, 0.31572,
        0.27542, 0.31585, 0.26444, 0.33257, 1.00000, 1.12726, 0.95402, 1.19067,
        1.00000, 1.15175, 0.96576, 1.21538, 0.00000, 0.03381, -0.01219, 0.05142,
        0.00000, 0.04108, -0.00919, 0.05805)
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
        0.00366, 0.00544, 0.00392, 0.00559, 1.00000, 0.00794, 0.00390, 0.00915,
        0.00000, 0.00240, 0.00101, 0.00276, 0.26430, 0.29854, 0.25272, 0.31572,
        0.27542, 0.31585, 0.26444, 0.33257, 1.00000, 1.12726, 0.95402, 1.19067,
        1.00000, 1.15175, 0.96576, 1.21538, 0.00000, 0.03381, -0.01219, 0.05142,
        0.00000, 0.04108, -0.00919, 0.05805)
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
      c(0.26851, NA, NA, NA, 0.26874, 0.24935, 0.26874, 0.28686, 1.00000, 0.92784, 1.00000, 1.06744,
        0.000000, -0.01939, 0.00000, 0.01812, 0.00366, 0.00795, 0.00366, 0.00417,
        1.00000, 0.02274, 0.00000, 0.00472, 0.00000, 0.00614, 0.00000, 0.00131,
        0.26430, 0.24356, 0.26430, 0.27986, 0.27542, 0.26676, 0.27542, 0.29260,
        1.00000, 0.90431, 1.00000, 1.05708, 1.00000, 0.97245, 1.00000, 1.07099,
        0.00000, -0.02610, 0.00000, 0.01518, 0.00000, -0.00756, 0.00000, 0.01919)
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
        0.00000, -0.01939, 0.07204, 0.00366, 0.00795, 0.00751, 1.00000, 0.02274, 0.01282,
        0.00000, 0.00614, 0.00423, 0.26430, 0.24356, 0.33125, 0.27542, 0.26676, 0.35561,
        1.00000, 0.90431, 1.25330, 1.00000, 0.97245, 1.29152, 0.00000, -0.02610, 0.06695,
        0.00000, -0.00756, 0.08019)
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
        0.00000, 0.06319, 0.02068, 0.00450, 0.01142, 0.00767, 1.00000, 0.04317, 0.02536,
        0.00000, 0.00964, 0.00558, 0.21127, 0.26626, 0.23060, 0.22274, 0.30049, 0.25383,
        1.00000, 1.24358, 1.07596, 1.00000, 1.36736, 1.14217, 0.00000, 0.05360, 0.01633,
        0.00000, 0.08027, 0.03140)
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
        0.00448, 0.00908, 0.00530, 1.00000, 0.03240, 0.00661, 0.00000, 0.00719, 0.00156, 0.21209, 0.22501, 0.23136,
        0.22335, 0.24974, 0.24656, 1.00000, 1.03170, 1.08809, 1.00000, 1.12947, 1.10514,
        0.00000, 0.00698, 0.01885, 0.00000, 0.02862, 0.02329)
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
        0.00349, 0.00662, 0.00752, 1.00000, 0.02063, 0.01777, 0.00000, 0.00557, 0.00523, 0.26351, 0.22851, 0.34562,
        0.27348, 0.24832, 0.36629, 1.00000, 0.85777, 1.30300, 1.00000, 0.91154, 1.34794,
        0.00000, -0.03878, 0.08176, 0.00000, -0.02406, 0.09379)
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
                  0.00338, 0.02779, 0.03733, 1.00000, 0.09910, 0.13234, 0.00000, 0.02688, 0.03588, 0.26586, 0.19614, 0.24405,
                  0.27467, 0.28071, 0.34795, 1.00000, 0.72432, 0.90929, 1.00000, 1.02581, 1.28054,
                  0.00000, -0.07456, -0.02413, 0.00000, 0.00706, 0.07623))
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
                  0.00000, -0.05048, 0.03164, 0.00348, 0.02870, 0.03689, 1.00000, 0.10323, 0.13059,
                  0.00000, 0.02800, 0.03538, 0.26527, 0.18933, 0.24715, 0.27439, 0.27260, 0.35027,
                  1.00000, 0.70241, 0.92243, 1.00000, 0.99786, 1.28879, 0.00000, -0.08021, -0.02057,
                  0.00000, -0.00058, 0.07834))
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
                  0.00000, -0.04022, 0.03896, 0.00296, 0.02574, 0.03761, 1.00000, 0.09376, 0.13698,
                  0.00000, 0.02459, 0.03600, 0.25730, 0.19277, 0.24790, 0.26523, 0.27260, 0.35164,
                  1.00000, 0.73755, 0.95711, 1.00000, 1.03042, 1.33675,
                  0.0000, -0.0688, -0.0109, 0.00000, 0.00804, 0.08831))
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
                  0.00000, -0.04234, 0.01597, 0.00370, 0.02508, 0.03130, 1.00000, 0.09823, 0.12097,
                  0.00000, 0.02404, 0.02960, 0.23635, 0.15209, 0.19821, 0.25084, 0.25041, 0.32092,
                  1.00000, 0.63364, 0.82848, 1.00000, 1.01870, 1.30265, 0.00000, -0.08945, -0.04205,
                  0.00000, 0.00477, 0.07399))
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
                  0.00000, -0.06670, 0.07803, 0.00401, 0.00785, 0.00775, 1.00000, 0.02301, 0.01227,
                  0.00000, 0.00594, 0.00415, 0.26219, 0.19275, 0.33415, 0.27367, 0.21568, 0.35862,
                  1.00000, 0.71721, 1.27439, 1.00000, 0.78811, 1.31171,
                  0.00000, -0.07603, 0.07195, 0.00000, -0.05799, 0.08511))
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
                  0.00401, 0.00785, 0.00775, 0.04072, 1.00000, 0.04982, 0.00594, 0.00000, 0.00627,
                  0.26219, 0.19275, 0.33415, 0.27367, 0.21568, 0.35862, 1.26890, 1.00000, 1.66217,
                  1.39451, 1.00000, 1.80204, 0.05799, 0.00000, 0.13700, 0.07603, 0.00000, 0.15461))
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
      as.vector(c(0.27178, NA, NA, 0.26739, 0.22322, 0.34933, 1.00000, 0.83480, 1.30645, 0.00000, -0.04417, 0.08194,
                  0.00375, 0.00664, 0.00766, 1.00000, 0.01703, 0.01819, 0.00000, 0.00440, 0.00528,
                  0.26304, 0.21706, 0.33996, 0.27424, 0.23700, 0.36130, 1.0000, 0.81849, 1.27760,
                  1.00000, 0.86589, 1.32655, 0.00000, -0.04916, 0.07516, 0.0000, -0.0367, 0.0875))
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
      as.vector(c(0.22985, NA, NA, 0.22724, 0.19630, 0.24222,
                  1.00000, 0.86387, 1.06592, 0.00000, -0.03093, 0.01498))
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

