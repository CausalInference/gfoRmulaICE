data <- gfoRmulaICE::compData
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
        1.00000, 1.13746, 0.95957, 1.20286, 0.00000, 0.03694, -0.01087, 0.05452)
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
        1.00000, 1.13746, 0.95957, 1.20286, 0.00000, 0.03694, -0.01087, 0.05452)
    )})

test_that(
  "check different interventions",
  {
    set.seed(1)
    
    ice_fit2 <- ice(data = data, time_points = 4,
                    id = "id", time_name = "t0",
                    censor_name = "C", outcome_name = "Y",
                    compevent_name = "D",
                    comp_effect = 0,
                    outcome_model = Y ~ L1 + L2 + A1 + A2,
                    censor_model = C ~ L1 + L2 + A1 + A2,
                    ref_idx = 0,
                    estimator = pool(hazard = F),
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
      c(0.26851, NA, NA, NA, 0.26874, 0.24935, 0.26874, 0.28701, 1.00000, 0.92784, 1.00000, 1.06799,
        0.000000, -0.01939, 0.00000, 0.01827)
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
        0.00000, -0.01939, 0.07204)
    )})

test_that(
  "check complicated scenario 2 - doubly robust ICE with intervention-specific models",
  {
    library(splines)
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
    
    library(splines)
    set.seed(1)
    
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
                    # nsamples = 5, ci_method = "percentile", parallel = T, ncores = 5,
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
                     # nsamples = 5, ci_method = "percentile", parallel = T, ncores = 5,
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
                     # nsamples = 5, ci_method = "percentile", parallel = T, ncores = 5,
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
                     # nsamples = 5, ci_method = "percentile", parallel = T, ncores = 5,
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
                     # nsamples = 5, ci_method = "percentile", parallel = T, ncores = 5,
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
                     # nsamples = 5, ci_method = "percentile", parallel = T, ncores = 5,
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
                     # nsamples = 5, ci_method = "percentile", parallel = T, ncores = 5,
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
                     # nsamples = 5, ci_method = "percentile", parallel = T, ncores = 5,
                     intervention1.A1 = list(static(3), 0:1),
                     intervention1.A2 = list(static(1), 1:2),
                     intervention2.A1 = list(case_when(data$L2 < 0 ~ 1,
                                                       data$L2 >= 0 & data$L2 < 2 ~ 2, T ~ 3), 2:3),
                     intervention2.A2 = list(dynamic("L1 == 0", static(0), static(1)), 0:2),
                     outcomeModel.1 = Y ~ L1 + L2,
                     compModel.2 = D ~ L1 + L2
    )
    
    out <- as.matrix(summary(ice_fit7a, ice_fit7b, ice_fit7c, ice_fit7d, 
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

