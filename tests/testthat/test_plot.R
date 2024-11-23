data <- gfoRmulaICE::compData
set.seed(1)

## test plots for ICE

test_that(
  "plot classical pooled ICE direct effect - dynamic interventions",
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
                    nsamples = 5, ci_method = "percentile", parallel = F, ncores = 5,
                    intervention1.A2 = list(dynamic("L1 == 0", static(0), static(1))),
                    intervention2.A2 = list(dynamic("L1 == 0", static(0), static(1), absorb = T)),
                    intervention3.A2 = list(dynamic("L1 == 0", static(0), natural_course()))
    )

    plot1 <- plot(ice_fit1)
    expect_equal(class(plot1), c("gg", "ggplot"))
    })

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
                    nsamples = 5, ci_method = "percentile", parallel = F, ncores = 5,
                    int_descript = c("Static Intervention", "Threshold Intervention",
                                     "Dynamic Intervention with Grace Period"),
                    intervention1.A1 = list(static(3)),
                    intervention1.A2 = list(static(1)),
                    intervention2.L2 = list(threshold(-3, Inf)),
                    intervention3.A2 = list(grace_period("uniform", 2, "L1 == 0"))
    )
    
    plot2 <- plot(ice_fit2)
    expect_equal(class(plot2), c("gg", "ggplot"))
    })

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
                    nsamples = 5, ci_method = "percentile", parallel = F, ncores = 5,
                    int_descript = c("Static Intervention", "Dynamic Intervention"),
                    intervention1.A1 = list(static(3)),
                    intervention1.A2 = list(static(1)),
                    intervention2.A1 = list(case_when(data$L2 < 0 ~ 1,
                                                      data$L2 >= 0 & data$L2 < 2 ~ 2, T ~ 3)),
                    intervention2.A2 = list(dynamic("L1 == 0", static(0), static(1)))
    )
    
    plot3 <- plot(ice_fit3)
    expect_equal(class(plot3), c("gg", "ggplot"))
    })

test_that(
  "plot hazard-based pooled ICE - time-specific hazard model same as the outcome model",
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
                     nsamples = 5, ci_method = "percentile", parallel = F, ncores = 5,
                     int_descript = c("Static Intervention", "Dynamic Intervention"),
                     intervention1.A1 = list(static(3)),
                     intervention1.A2 = list(static(1)),
                     intervention2.A1 = list(case_when(data$L2 < 0 ~ 1,
                                                       data$L2 >= 0 & data$L2 < 2 ~ 2, T ~ 3)),
                     intervention2.A2 = list(dynamic("L1 == 0", static(0), static(1)))
    )
    
    plot4a <- plot(ice_fit4a)
    expect_equal(class(plot4a), c("gg", "ggplot"))
    })

test_that(
  "plot hazard-based pooled ICE - time-specific hazard model Y ~ L1 + L2",
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
                     nsamples = 5, ci_method = "percentile", parallel = F, ncores = 5,
                     int_descript = c("Static Intervention", "Dynamic Intervention"),
                     intervention1.A1 = list(static(3)),
                     intervention1.A2 = list(static(1)),
                     intervention2.A1 = list(case_when(data$L2 < 0 ~ 1,
                                                       data$L2 >= 0 & data$L2 < 2 ~ 2, T ~ 3)),
                     intervention2.A2 = list(dynamic("L1 == 0", static(0), static(1)))
    )
    
    plot4b <- plot(ice_fit4b)
    expect_equal(class(plot4b), c("gg", "ggplot"))
    })

test_that(
  "plot hazard-based pooled ICE - pooled-over-time hazard model",
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
                     nsamples = 5, ci_method = "percentile", parallel = F, ncores = 5,
                     int_descript = c("Static Intervention", "Dynamic Intervention"),
                     intervention1.A1 = list(static(3)),
                     intervention1.A2 = list(static(1)),
                     intervention2.A1 = list(case_when(data$L2 < 0 ~ 1,
                                                       data$L2 >= 0 & data$L2 < 2 ~ 2, T ~ 3)),
                     intervention2.A2 = list(dynamic("L1 == 0", static(0), static(1)))
    )
    
    plot4c <- plot(ice_fit4c)
    expect_equal(class(plot4c), c("gg", "ggplot"))
    })

test_that(
  "plot classical stratified ICE",
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
                     nsamples = 5, ci_method = "percentile", parallel = F, ncores = 5,
                     int_descript = c("Static Intervention", "Dynamic Intervention"),
                     intervention1.A1 = list(static(3)),
                     intervention1.A2 = list(static(1)),
                     intervention2.A1 = list(case_when(data$L2 < 0 ~ 1,
                                                       data$L2 >= 0 & data$L2 < 2 ~ 2, T ~ 3)),
                     intervention2.A2 = list(dynamic("L1 == 0", static(0), static(1)))
    )
    
    plot4d <- plot(ice_fit4d)
    expect_equal(class(plot4d), c("gg", "ggplot"))
    })

test_that(
  "plot hazard-extended stratified ICE - time-specific hazard model Y ~ L1",
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
                     nsamples = 5, ci_method = "percentile", parallel = F, ncores = 5,
                     int_descript = c("Static Intervention", "Dynamic Intervention"),
                     intervention1.A1 = list(static(3)),
                     intervention1.A2 = list(static(1)),
                     intervention2.A1 = list(case_when(data$L2 < 0 ~ 1,
                                                       data$L2 >= 0 & data$L2 < 2 ~ 2, T ~ 3)),
                     intervention2.A2 = list(dynamic("L1 == 0", static(0), static(1)))
    )
    
    plot4e <- plot(ice_fit4e)
    expect_equal(class(plot4e), c("gg", "ggplot"))
    })

test_that(
  "plot doubly robust ICE",
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
                     nsamples = 5, ci_method = "percentile", parallel = F, ncores = 5,
                     int_descript = c("Static Intervention", "Dynamic Intervention"),
                     intervention1.A1 = list(static(3)),
                     intervention1.A2 = list(static(1)),
                     intervention2.A1 = list(case_when(data$L2 < 0 ~ 1,
                                                       data$L2 >= 0 & data$L2 < 2 ~ 2, T ~ 3)),
                     intervention2.A2 = list(dynamic("L1 == 0", static(0), static(1)))
    )
    
    plot4f <- plot(ice_fit4f)
    expect_equal(class(plot4f), c("gg", "ggplot"))
    })

test_that(
  "plot hazard-based stratified ICE with intervention-specific models",
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
                     nsamples = 5, ci_method = "normal", parallel = F, ncores = 5,
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
    
    plot4h <- plot(ice_fit4h)
    expect_equal(class(plot4h), c("gg", "ggplot"))
    })

test_that(
  "plot flexible model specification - complicated terms",
  {
    library(Hmisc)
    ice_fit5a <- ice(data = data, time_points = 4,
                     id = "id", time_name = "t0",
                     censor_name = "C", outcome_name = "Y",
                     compevent_name = "D",
                     comp_effect = 0,
                     outcome_model = Y ~ I(L1^2) + rcspline.eval(lag1_L2, knots = 1:3) + A1 + A2,
                     censor_model = C ~ lag1_L1 + poly(L2, degree = 2) + A1 + A2,
                     ref_idx = 0,
                     estimator = pool(hazard = F),
                     nsamples = 5, ci_method = "percentile", parallel = F, ncores = 5,
                     int_descript = c("Static Intervention", "Dynamic Intervention"),
                     intervention1.A1 = list(static(3)),
                     intervention1.A2 = list(static(1)),
                     intervention2.A1 = list(case_when(data$L2 < 0 ~ 1,
                                                       data$L2 >= 0 & data$L2 < 2 ~ 2, T ~ 3)),
                     intervention2.A2 = list(dynamic("L1 == 0", static(0), static(1)))
    )
    
    plot5a <- plot(ice_fit5a)
    expect_equal(class(plot5a), c("gg", "ggplot"))
    })

test_that(
  "plot flexible model specification - using static intervention as reference",
  {
    library(Hmisc)
    ice_fit5b <- ice(data = data, time_points = 4,
                     id = "id", time_name = "t0",
                     censor_name = "C", outcome_name = "Y",
                     compevent_name = "D",
                     comp_effect = 0,
                     outcome_model = Y ~ I(L1^2) + rcspline.eval(lag1_L2, knots = 1:3) + A1 + A2,
                     censor_model = C ~ lag1_L1 + poly(L2, degree = 2) + A1 + A2,
                     ref_idx = 1,
                     estimator = pool(hazard = F),
                     nsamples = 5, ci_method = "percentile", parallel = F, ncores = 5,
                     int_descript = c("Static Intervention", "Dynamic Intervention"),
                     intervention1.A1 = list(static(3)),
                     intervention1.A2 = list(static(1)),
                     intervention2.A1 = list(case_when(data$L2 < 0 ~ 1,
                                                       data$L2 >= 0 & data$L2 < 2 ~ 2, T ~ 3)),
                     intervention2.A2 = list(dynamic("L1 == 0", static(0), static(1)))
    )
    
    plot5b <- plot(ice_fit5b)
    expect_equal(class(plot5b), c("gg", "ggplot"))
    })

test_that(
  "plot complicated scenario 1",
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
                     nsamples = 5, ci_method = "percentile", parallel = F, ncores = 5,
                     int_descript = c("Static Intervention", "Dynamic Intervention"),
                     intervention1.A1 = list(static(3)),
                     intervention1.A2 = list(static(1)),
                     intervention2.A1 = list(case_when(data$L2 < 0 ~ 1,
                                                       data$L2 >= 0 & data$L2 < 2 ~ 2, T ~ 3)),
                     intervention2.A2 = list(dynamic("L1 == 0", static(0), static(1)))
    )
    
    plot6a <- plot(ice_fit6a)
    expect_equal(class(plot6a), c("gg", "ggplot"))
    })

test_that(
  "plot complicated scenario 2 - doubly robust ICE with intervention-specific models",
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
    
    plot6b <- plot(ice_fit6b)
    expect_equal(class(plot6b), c("gg", "ggplot"))
    })

test_that(
  "plot complicated scenario 3 - intervention-specific time options",
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
    
    plots <- plot(ice_fit7a, ice_fit7b, ice_fit7c, ice_fit7d, 
                       ice_fit7e, ice_fit7f, ice_fit7g, ice_fit7h)
    expect_equal(class(plots), c("gg", "ggplot"))
 
    })

