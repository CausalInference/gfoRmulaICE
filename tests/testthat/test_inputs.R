load(file = "../../data/comp_and_censor_data.rda")
library(tidyverse)
library(stringr)
library(nnet)
library(doRNG)

test_that(
  "check pooled ICE inputs",
  {

    expect_error(ice(data = data, K = 4, id = "i",
                     time_name = "t0", outcome_name = "Y",
                     censor_name = NULL, compevent_name = "D",
                     estimator = pool(hazard = F),
                     comp_effect = 2,
                     outcome_model = Y ~ L1 + A2,
                     ref_idx = 0,
                     int_descript = c("Always Treat", "Never Treat"),
                     intervention1.A2 = list(static(1)),
                     intervention2.A2 = list(static(0))))

    expect_error(ice(data = data, K = 4, id = "id",
                     time_name = "t", outcome_name = "Y",
                     censor_name = NULL, compevent_name = "D",
                     estimator = pool(hazard = F),
                     comp_effect = 2,
                     outcome_model = Y ~ L1 + A2,
                     ref_idx = 0,
                     int_descript = c("Always Treat", "Never Treat"),
                     intervention1.A2 = list(static(1)),
                     intervention2.A2 = list(static(0))))

    expect_error(ice(data = data, K = 4, id = "id",
                     time_name = "t0",
                     censor_name = NULL, compevent_name = "d",
                     estimator = pool(hazard = F),
                     comp_effect = 2,
                     outcome_model = Y ~ L1 + A2,
                     ref_idx = 0,
                     int_descript = c("Always Treat", "Never Treat"),
                     intervention_1.2 = list(static(1)),
                     intervention_2.2 = list(static(0))))

    expect_error(ice(data = data, K = 4, id = "id",
                     time_name = "t0", outcome_name = "Y",
                     censor_name = "c", compevent_name = "D",
                     estimator = pool(hazard = F),
                     comp_effect = 2,
                     outcome_model = Y ~ L1 + A2,
                     ref_idx = 0,
                     int_descript = c("Always Treat", "Never Treat"),
                     intervention1.A2 = list(static(1)),
                     intervention2.A2 = list(static(0))))

    expect_error(ice(data = data, K = 4, id = "id",
                     time_name = "t0", outcome_name = "Y",
                     censor_name = NULL, compevent_name = "D",
                     estimator = pool(hazard = F),
                     comp_effect = 2,
                     outcome_model = Y ~ L1 + A2,
                     ref_idx = 0,
                     int_descript = c("Always Treat", "Never Treat"),
                     intervention_1.A2 = list(static(1)),
                     intervention_2.A2 = list(static(0))))

    expect_error(ice(data = data, K = 4, id = "id",
                     time_name = "t0", outcome_name = "Y",
                     censor_name = "C", compevent_name = "D",
                     estimator = pool(hazard = F),
                     comp_effect = 2,
                     outcome_model = Y ~ L1 + A2,
                     censor_model = NULL,
                     competing_model = D ~ L1 + A2,
                     ref_idx = 0,
                     int_descript = c("Always Treat", "Never Treat"),
                     intervention1.A2 = list(static(1)),
                     intervention2.A2 = list(static(0))))

    expect_error(ice(data = data, K = 4, id = "id",
                     time_name = "t0", outcome_name = "Y",
                     censor_name = NULL, compevent_name = "D",
                     estimator = pool(hazard = T),
                     comp_effect = 2,
                     outcome_model = Y ~ L1 + A2,
                     competing_model = NULL,
                     ref_idx = 0,
                     int_descript = c("Always Treat", "Never Treat"),
                     intervention1.A2 = list(static(1)),
                     intervention2.A2 = list(static(0))))

    expect_error(ice(data = data, K = 4, id = "id",
                     time_name = "t0", outcome_name = "Y",
                     censor_name = NULL, compevent_name = "D",
                     estimator = pool(hazard = F),
                     comp_effect = 2,
                     outcome_model = Y ~ L1,
                     competing_model = D ~ L1 + A2,
                     ref_idx = 0,
                     int_descript = c("Always Treat", "Never Treat"),
                     intervention1.A2 = list(static(1)),
                     intervention2.A2 = list(static(0))))

    expect_error(ice(data = test_data, K = 5, id = "id",
                     time_name = "t0", outcome_name = "Y",
                     censor_name = NULL, compevent_name = "D",
                     estimator = pool(hazard = T),
                     comp_effect = 2,
                     outcome_model = Y ~ L1 + A2,
                     competing_model = D ~ L1,
                     ref_idx = 0,
                     int_descript = c("Always Treat", "Never Treat"),
                     intervention1.A2 = list(static(1)),
                     intervention2.A2 = list(static(0))))


    })

test_that(
  "check stratified ICE inputs",
  {

    expect_error(ice(data = data, K = 4, id = "i",
                     time_name = "t0", outcome_name = "Y",
                     censor_name = NULL, compevent_name = "D",
                     estimator = strat(hazard = F),
                     comp_effect = 2,
                     outcome_model = Y ~ L1,
                     ref_idx = 0,
                     int_descript = c("Always Treat", "Never Treat"),
                     intervention1.A2 = list(static(1)),
                     intervention2.A2 = list(static(0))))

    expect_error(ice(data = data, K = 4, id = "id",
                     time_name = "t", outcome_name = "Y",
                     censor_name = NULL, compevent_name = "D",
                     estimator = strat(hazard = F),
                     comp_effect = 2,
                     outcome_model = Y ~ L1,
                     ref_idx = 0,
                     int_descript = c("Always Treat", "Never Treat"),
                     intervention1.A2 = list(static(1)),
                     intervention2.A2 = list(static(0))))

    expect_error(ice(data = data, K = 4, id = "id",
                     time_name = "t0", outcome_name = "Y",
                     censor_name = NULL, compevent_name = "d",
                     estimator = strat(hazard = F),
                     comp_effect = 2,
                     outcome_model = Y ~ L1,
                     ref_idx = 0,
                     int_descript = c("Always Treat", "Never Treat"),
                     intervention1.A2 = list(static(1)),
                     intervention2.A2 = list(static(0))))

    expect_error(ice(data = data, K = 4, id = "id",
                     time_name = "t0", outcome_name = "Y",
                     censor_name = "c", compevent_name = "D",
                     estimator = strat(hazard = F),
                     comp_effect = 2,
                     outcome_model = Y ~ L1,
                     ref_idx = 0,
                     int_descript = c("Always Treat", "Never Treat"),
                     intervention1.A2 = list(static(1)),
                     intervention2.A2 = list(static(0))))

    expect_error(ice(data = data, K = 4, id = "id",
                     time_name = "t0", outcome_name = "Y",
                     censor_name = NULL, compevent_name = "D",
                     estimator = strat(hazard = F),
                     comp_effect = 2,
                     outcome_model = Y ~ L1,
                     ref_idx = 0,
                     int_descript = c("Always Treat", "Never Treat"),
                     intervention1.a2 = list(static(1)),
                     intervention2.a2 = list(static(0))))

    expect_error(ice(data = data, K = 4, id = "id",
                     time_name = "t0", outcome_name = "Y",
                     censor_name = "C", compevent_name = "D",
                     estimator = pool(hazard = F),
                     comp_effect = 2,
                     outcome_model = Y ~ L1,
                     censor_model = NULL,
                     competing_model = D ~ L1,
                     ref_idx = 0,
                     int_descript = c("Always Treat", "Never Treat"),
                     intervention1.A2 = list(static(1)),
                     intervention2.A2 = list(static(0))))

    expect_error(ice(data = data, K = 4, id = "id",
                     time_name = "t0", outcome_name = "Y",
                     censor_name = NULL, compevent_name = "D",
                     estimator = strat(hazard = T),
                     comp_effect = 2,
                     outcome_model = Y ~ L1,
                     competing_model = NULL,
                     ref_idx = 0,
                     int_descript = c("Always Treat", "Never Treat"),
                     intervention1.A2 = list(static(1)),
                     intervention2.A2 = list(static(0))))

    expect_error(ice(data = data, K = 4, id = "id",
                     time_name = "t0", outcome_name = "Y",
                     censor_name = NULL, compevent_name = "D",
                     estimator = strat(hazard = T),
                     comp_effect = 2,
                     outcome_model = Y ~ L1,
                     competing_model = D ~ L1 + A,
                     ref_idx = 0,
                     int_descript = c("Always Treat", "Never Treat"),
                     intervention1.A2 = list(static(1)),
                     intervention2.A2 = list(static(0))))

    expect_error(ice(data = data, K = 4, id = "id",
                     time_name = "t0", outcome_name = "Y",
                     censor_name = NULL, compevent_name = "D",
                     estimator = strat(hazard = T),
                     comp_effect = 2,
                     outcome_model = Y ~ L1 + A2,
                     competing_model = D ~ L1,
                     ref_idx = 0,
                     int_descript = c("Always Treat", "Never Treat"),
                     intervention1.A2 = list(static(1)),
                     intervention2.A2 = list(static(0))))

  })

test_that(
  "check weighted ICE inputs",
  {

    expect_error(ice(data = data, K = 4, id = "i",
                     time_name = "t0", outcome_name = "Y",
                     censor_name = NULL, compevent_name = "D",
                     estimator = weight(list(A2 ~ L1)),
                     comp_effect = 2,
                     outcome_model = Y ~ L1,
                     ref_idx = 0,
                     int_descript = c("Always Treat", "Never Treat"),
                     intervention1.A2 = list(static(1)),
                     intervention2.A2 = list(static(0))))

    expect_error(ice(data = test_data, K = 4, id = "id",
                     time_name = "t", outcome_name = "Y",
                     censor_name = NULL, compevent_name = "D",
                     estimator = weight(list(A2 ~ L1)),
                     comp_effect = 2,
                     outcome_model = Y ~ L1,
                     ref_idx = 0,
                     int_descript = c("Always Treat", "Never Treat"),
                     intervention1.A2 = list(static(1)),
                     intervention2.A2 = list(static(0))))

    expect_error(ice(data = test_data, K = 4, id = "id",
                     time_name = "t0", outcome_name = "Y",
                     censor_name = NULL, compevent_name = "d",
                     estimator = weight(list(A2 ~ L1)),
                     comp_effect = 2,
                     outcome_model = Y ~ L1,
                     ref_idx = 0,
                     int_descript = c("Always Treat", "Never Treat"),
                     intervention1.A2 = list(static(1)),
                     intervention2.A2 = list(static(0))))

    expect_error(ice(data = test_data, K = 4, id = "id",
                     time_name = "t0", outcome_name = "Y",
                     censor_name = "c", compevent_name = "D",
                     estimator = weight(list(A2 ~ L1)),
                     comp_effect = 2,
                     outcome_model = Y ~ L1,
                     ref_idx = 0,
                     int_descript = c("Always Treat", "Never Treat"),
                     intervention1.A2 = list(static(1)),
                     intervention2.A2 = list(static(0))))

    expect_error(ice(data = test_data, K = 4, id = "id",
                     time_name = "t0", outcome_name = "Y",
                     censor_name = NULL, compevent_name = "D",
                     estimator = weight(list(A2 ~ L1)),
                     comp_effect = 2,
                     outcome_model = Y ~ L1,
                     ref_idx = 0,
                     int_descript = c("Always Treat", "Never Treat"),
                     intervention1.a2 = list(static(1)),
                     intervention2.a2 = list(static(0))))

    expect_error(ice(data = test_data, K = 4, id = "id",
                     time_name = "t0", outcome_name = "Y",
                     censor_name = "C", compevent_name = "D",
                     estimator = weight(list(A2 ~ L1)),
                     comp_effect = 2,
                     outcome_model = Y ~ L1,
                     censor_model = NULL,
                     competing_model = D ~ L1,
                     ref_idx = 0,
                     int_descript = c("Always Treat", "Never Treat"),
                     intervention1.A2 = list(static(1)),
                     intervention2.A2 = list(static(0))))

    expect_error(ice(data = test_data, K = 4, id = "id",
                     time_name = "t0", outcome_name = "Y",
                     censor_name = "C", compevent_name = "D",
                     estimator = weight(list(A2 ~ L1)),
                     comp_effect = 2,
                     outcome_model = Y ~ L1 + A2,
                     censor_model = NULL,
                     competing_model = D ~ L1,
                     ref_idx = 0,
                     int_descript = c("Always Treat", "Never Treat"),
                     intervention1.A2 = list(static(1)),
                     intervention2.A2 = list(static(0))))

    expect_error(ice(data = test_data, K = 4, id = "id",
                     time_name = "t0", outcome_name = "Y",
                     censor_name = "C", compevent_name = "D",
                     estimator = weight(list('A2 ~ L1')),
                     comp_effect = 2,
                     outcome_model = Y ~ L1 + A2,
                     censor_model = NULL,
                     competing_model = D ~ L1,
                     ref_idx = 0,
                     int_descript = c("Always Treat", "Never Treat"),
                     intervention1.A2 = list(static(1)),
                     intervention2.A2 = list(static(0))))

    expect_error(ice(data = test_data, K = 4, id = "id",
                     time_name = "t0", outcome_name = "Y",
                     censor_name = "C", compevent_name = "D",
                     estimator = weight(list()),
                     comp_effect = 2,
                     outcome_model = Y ~ L1 + A2,
                     censor_model = NULL,
                     competing_model = D ~ L1,
                     ref_idx = 0,
                     int_descript = c("Always Treat", "Never Treat"),
                     intervention1.A2 = list(static(1)),
                     intervention2.A2 = list(static(0))))


  })

test_that(
  "check general inputs",
  {

    expect_error(ice(data = data, K = 4, id = "id",
                     time_name = "t0", outcome_name = "Y",
                     censor_name = NULL, compevent_name = "D",
                     estimator = pool(hazard = F),
                     comp_effect = 2,
                     outcome_model = Y ~ L1 + A2,
                     ref_idx = 0,
                     int_descript = c("Always Treat", "Never Treat"),
                     intervention1.A2 = list(static(1)),
                     intervention2.A2 = list(static(0))))

    expect_error(ice(data = data, K = 4, id = "id",
                     time_name = "t0", outcome_name = "Y",
                     censor_name = NULL, compevent_name = "D",
                     estimator = pool(hazard = F),
                     comp_effect = 2,
                     outcome_model = Y ~ L1 + A2,
                     ref_idx = 0,
                     int_descript = c("Always Treat", "Never Treat"),
                     intervention1.A2 = list(static(1)),
                     intervention2.A2 = list(static(0))))

    expect_error(ice(data = data, K = 4, id = "id",
                     time_name = "t0", outcome_name = "Y",
                     censor_name = NULL, compevent_name = "D",
                     estimator = pool(hazard = F),
                     comp_effect = 2,
                     outcome_model = Y ~ L1 + A2,
                     ref_idx = 0,
                     int_descript = c("Always Treat", "Never Treat"),
                     intervention1.A2 = list(static(1)),
                     intervention2.A2 = list(static(0))))

    expect_error(ice(data = data, K = 4, id = "id",
                     time_name = "t0", outcome_name = "Y",
                     censor_name = NULL, compevent_name = "D",
                     estimator = strat(hazard = F),
                     comp_effect = 2,
                     outcome_model = Y ~ L1,
                     ref_idx = 0,
                     int_descript = c("Always Treat", "Never Treat"),
                     intervention1.A2 = list(static(1)),
                     intervention2.A2 = list(static(0))))

    expect_error(ice(data = data, K = 4, id = "id",
                     time_name = "t0", outcome_name = "Y",
                     censor_name = NULL, compevent_name = "D",
                     estimator = strat(hazard = F),
                     comp_effect = 2,
                     outcome_model = Y ~ L1,
                     ref_idx = 0,
                     int_descript = c("Always Treat", "Never Treat"),
                     intervention1.A2 = list(static(1)),
                     intervention2.A2 = list(static(0))))

    expect_error(ice(data = data, K = 4, id = "id",
                     time_name = "t0", outcome_name = "Y",
                     censor_name = NULL, compevent_name = "D",
                     estimator = strat(hazard = F),
                     comp_effect = 2,
                     outcome_model = Y ~ L1,
                     ref_idx = 0,
                     int_descript = c("Always Treat", "Never Treat"),
                     intervention1.A2 = list(static(1)),
                     intervention2.A2 = list(static(0))))


  })


