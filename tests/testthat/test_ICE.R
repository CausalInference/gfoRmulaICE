load(file = "../../data/comp_and_censor_data.rda")
library(tidyverse)
library(stringr)
library(nnet)
library(doParallel)

test_that(
  "check classical pooled ICE direct effect",
  {
    out <- as.matrix(ice(data = data,
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
                     intervention2.A1 = list(case_when(data$L2 < 0 ~ 1,
                                                       data$L2 >= 0 & data$L2 < 2 ~ 2,
                                                       T ~ 3)),
                     intervention2.A2 = list(dynamic("compare", "L1", "=", 0)))$summary)
    rownames(out) <- c()
    out <- as.vector(out)
    expect_equal(
      out,
      c(0.26947, NA, NA, 0.26874, 0.24935, 0.34078, 1.00000 ,
        0.92784, 1.26807, 0.00000, -0.01939, 0.07204)
    )})

test_that(
  "check classical stratified ICE direct effect",
  {
    out <- as.matrix(ice(
      data = data,
      time_points = 4,
      id = "id",
      time_name = "t0",
      outcome_name = "Y",
      censor_name = "C",
      compevent_name = "D",
      estimator = strat(hazard = F),
      comp_effect = 0,
      outcome_model = Y ~ L1 + L2,
      competing_model = D ~ L1 + L2,
      censor_model = C ~ L1 + L2,
      ref_idx = 0,
      int_descript = c("Static Intervention",
                       "Dynamic Intervention"),
      intervention1.A1 = list(static(3)),
      intervention1.A2 = list(static(1)),
      intervention2.A1 = list(case_when(data$L2 < 0 ~ 1,
                                        data$L2 >= 0 & data$L2 < 2 ~ 2,
                                        T ~ 3)),
      intervention2.A2 = list(dynamic("compare", "L1", "=", 0))
    )$summary)
    rownames(out) <- c()
    out <- as.vector(out)
    expect_equal(
      out,
      c(0.26947, NA, NA, 0.26967, 0.22639, 0.29826, 1.00000, 0.83951, 1.10602,
        0.00000, -0.04328, 0.02859)
    )})

test_that(
  "check hazard-extended pooled ICE direct effect",
  {
    out <- as.matrix(ice(
      data = data,
      time_points = 4,
      id = "id",
      time_name = "t0",
      outcome_name = "Y",
      censor_name = "C",
      compevent_name = "D",
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
      intervention2.A1 = list(case_when(data$L2 < 0 ~ 1,
                                        data$L2 >= 0 & data$L2 < 2 ~ 2,
                                        T ~ 3)),
      intervention2.A2 = list(dynamic("compare", "L1", "=", 0))
    )$summary)
    rownames(out) <- c()
    out <- as.vector(out)
    expect_equal(
      out,
      c(0.26947, NA, NA, 0.26809, 0.23528, 0.35548, 1.00000, 0.87763, 1.32597,
        0.00000, -0.03281, 0.08739)
    )})

test_that(
  "check hazard-extended stratified ICE direct effect",
  {
    out <- as.matrix(ice(
      data = data,
      time_points = 4,
      id = "id",
      time_name = "t0",
      outcome_name = "Y",
      censor_name = "C",
      compevent_name = "D",
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
      intervention2.A1 = list(case_when(data$L2 < 0 ~ 1,
                                        data$L2 >= 0 & data$L2 < 2 ~ 2,
                                        T ~ 3)),
      intervention2.A2 = list(dynamic("compare", "L1", "=", 0))
    )$summary)
    rownames(out) <- c()
    out <- as.vector(out)
    expect_equal(
      out,
      c(0.26947, NA, NA, 0.26936, 0.22562, 0.29913, 1.00000, 0.83763, 1.11053,
        0.00000, -0.04373, 0.02977)
    )})

test_that(
  "check weighted ICE direct effect",
  {
    out <- as.matrix(ice(
      data = data,
      time_points = 4,
      id = "id",
      time_name = "t0",
      outcome_name = "Y",
      censor_name = "C",
      compevent_name = "D",
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
      intervention2.A1 = list(case_when(data$L2 < 0 ~ 1,
                                        data$L2 >= 0 & data$L2 < 2 ~ 2,
                                        T ~ 3)),
      intervention2.A2 = list(dynamic("compare", "L1", "=", 0))
    )$summary)
    rownames(out) <- c()
    out <- as.vector(out)
    expect_equal(
      out,
      c(0.26947, NA, NA, 0.26162, 0.22232, 0.29975, 1.00000, 0.84976, 1.14572,
        0.00000, -0.03931, 0.03812)
    )})

test_that(
  "check classical pooled ICE total effect",
  {
    out <- as.matrix(ice(
      data = data,
      time_points = 4,
      id = "id",
      time_name = "t0",
      outcome_name = "Y",
      censor_name = "C",
      compevent_name = "D",
      estimator = pool(hazard = F),
      comp_effect = 1,
      outcome_model = Y ~ L1 + L2 + A1 + A2,
      competing_model = D ~ L1 + L2 + A1 + A2,
      censor_model = C ~ L1 + L2 + A1 + A2,
      ref_idx = 0,
      int_descript = c("Static Intervention",
                       "Dynamic Intervention"),
      intervention1.A1 = list(static(3)),
      intervention1.A2 = list(static(1)),
      intervention2.A1 = list(case_when(data$L2 < 0 ~ 1,
                                        data$L2 >= 0 & data$L2 < 2 ~ 2,
                                        T ~ 3)),
      intervention2.A2 = list(dynamic("compare", "L1", "=", 0))
    )$summary)
    rownames(out) <- c()
    out <- as.vector(out)
    expect_equal(
      out,
      c(0.22949, NA, NA, 0.22776, 0.21449, 0.27596, 1.00000, 0.94175, 1.21164,
        0.00000, -0.01327, 0.04820)
    )})

test_that(
  "check stratified pooled ICE total effect",
  {
    out <- as.matrix(ice(
      data = data,
      time_points = 4,
      id = "id",
      time_name = "t0",
      outcome_name = "Y",
      censor_name = "C",
      compevent_name = "D",
      estimator = strat(hazard = F),
      comp_effect = 1,
      outcome_model = Y ~ L1 + L2,
      competing_model = D ~ L1 + L2,
      censor_model = C ~ L1 + L2,
      ref_idx = 0,
      int_descript = c("Static Intervention",
                       "Dynamic Intervention"),
      intervention1.A1 = list(static(3)),
      intervention1.A2 = list(static(1)),
      intervention2.A1 = list(case_when(data$L2 < 0 ~ 1,
                                        data$L2 >= 0 & data$L2 < 2 ~ 2,
                                        T ~ 3)),
      intervention2.A2 = list(dynamic("compare", "L1", "=", 0))
    )$summary)
    rownames(out) <- c()
    out <- as.vector(out)
    expect_equal(
      out,
      as.vector(c(0.22949, NA, NA, 0.23119, 0.19974, 0.24189, 1.00000, 0.86397, 1.04629,
                  0.00000, -0.03145, 0.01070))
    )})

test_that(
  "check hazard-extended pooled ICE total effect",
  {
    out <- as.matrix(ice(
      data = data,
      time_points = 4,
      id = "id",
      time_name = "t0",
      outcome_name = "Y",
      censor_name = "C",
      compevent_name = "D",
      estimator = pool(hazard = T),
      comp_effect = 1,
      outcome_model = Y ~ L1 + L2 + A1 + A2,
      competing_model = D ~ L1 + L2 + A1 + A2,
      censor_model = C ~ L1 + L2 + A1 + A2,
      ref_idx = 0,
      int_descript = c("Static Intervention",
                       "Dynamic Intervention"),
      intervention1.A1 = list(static(3)),
      intervention1.A2 = list(static(1)),
      intervention2.A1 = list(case_when(data$L2 < 0 ~ 1,
                                        data$L2 >= 0 & data$L2 < 2 ~ 2,
                                        T ~ 3)),
      intervention2.A2 = list(dynamic("compare", "L1", "=", 0))
    )$summary)
    rownames(out) <- c()
    out <- as.vector(out)
    expect_equal(
      out,
      as.vector(c(0.22949, NA, NA, 0.22582, 0.20231, 0.28264, 1.00000, 0.89588, 1.25163,
                  0.00000, -0.02351, 0.05682))
    )})

test_that(
  "check hazard-extended stratified ICE total effect",
  {
    out <- as.matrix(ice(
      data = data,
      time_points = 4,
      id = "id",
      time_name = "t0",
      outcome_name = "Y",
      censor_name = "C",
      compevent_name = "D",
      estimator = strat(hazard = T),
      comp_effect = 1,
      outcome_model = Y ~ L1 + L2,
      competing_model = D ~ L1 + L2,
      censor_model = C ~ L1 + L2,
      ref_idx = 0,
      int_descript = c("Static Intervention",
                       "Dynamic Intervention"),
      intervention1.A1 = list(static(3)),
      intervention1.A2 = list(static(1)),
      intervention2.A1 = list(case_when(data$L2 < 0 ~ 1,
                                        data$L2 >= 0 & data$L2 < 2 ~ 2,
                                        T ~ 3)),
      intervention2.A2 = list(dynamic("compare", "L1", "=", 0))
    )$summary)
    rownames(out) <- c()
    out <- as.vector(out)
    expect_equal(
      out,
      as.vector(c(0.22949, NA, NA, 0.24340, 0.20589, 0.25813, 1.00000, 0.84587, 1.06052,
                  0.00000, -0.03751, 0.01473))
    )})

test_that(
  "check weighted ICE total effect",
  {
    out <- as.matrix(ice(
      data = data,
      time_points = 4,
      id = "id",
      time_name = "t0",
      outcome_name = "Y",
      censor_name = "C",
      compevent_name = "D",
      estimator = weight(list(A1 ~ L1 + L2, A2 ~ L1 + L2)),
      comp_effect = 1,
      outcome_model = Y ~ L1 + L2,
      competing_model = D ~ L1 + L2,
      censor_model = C ~ L1 + L2,
      ref_idx = 0,
      int_descript = c("Static Intervention",
                       "Dynamic Intervention"),
      intervention1.A1 = list(static(3)),
      intervention1.A2 = list(static(1)),
      intervention2.A1 = list(case_when(data$L2 < 0 ~ 1,
                                        data$L2 >= 0 & data$L2 < 2 ~ 2,
                                        T ~ 3)),
      intervention2.A2 = list(dynamic("compare", "L1", "=", 0))
    )$summary)
    rownames(out) <- c()
    out <- as.vector(out)
    expect_equal(
      out,
      as.vector(c(0.22949, NA, NA, 0.22683, 0.19755, 0.24260, 1.00000, 0.87094, 1.06955,
                  0.00000, -0.02927, 0.01578))
    )})

test_that(
  "check classical pooled ICE direct effect with intervention 1 as reference",
  {
    out <- as.matrix(ice(
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
      ref_idx = 1,
      int_descript = c("Static Intervention",
                       "Dynamic Intervention"),
      intervention1.A1 = list(static(3)),
      intervention1.A2 = list(static(1)),
      intervention2.A1 = list(case_when(data$L2 < 0 ~ 1,
                                        data$L2 >= 0 & data$L2 < 2 ~ 2,
                                        T ~ 3)),
      intervention2.A2 = list(dynamic("compare", "L1", "=", 0))
    )$summary)
    rownames(out) <- c()
    out <- as.vector(out)
    expect_equal(
      out,
      as.vector(c(0.26947, NA, NA, 0.26874, 0.24935, 0.34078, 1.07778, 1.00000, 1.36669,
                  0.01939, 0.00000, 0.09143))
    )})

test_that(
  "check classical pooled ICE direct effect with intervention 2 as reference",
  {
    out <- as.matrix(ice(
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
      ref_idx = 2,
      int_descript = c("Static Intervention",
                       "Dynamic Intervention"),
      intervention1.A1 = list(static(3)),
      intervention1.A2 = list(static(1)),
      intervention2.A1 = list(case_when(data$L2 < 0 ~ 1,
                                        data$L2 >= 0 & data$L2 < 2 ~ 2,
                                        T ~ 3)),
      intervention2.A2 = list(dynamic("compare", "L1", "=", 0))
    )$summary)
    rownames(out) <- c()
    out <- as.vector(out)
    expect_equal(
      out,
      as.vector(c(0.26947, NA, NA, 0.26874, 0.24935, 0.34078, 0.78860, 0.73169, 1.00000,
                  -0.07204, -0.09143, 0.00000))
    )})

test_that(
  "check classical stratified ICE direct effect with intervention 1 as reference",
  {
    out <- as.matrix(ice(
      data = data,
      time_points = 4,
      id = "id",
      time_name = "t0",
      outcome_name = "Y",
      censor_name = "C",
      compevent_name = "D",
      estimator = strat(hazard = F),
      comp_effect = 0,
      outcome_model = Y ~ L1 + L2,
      competing_model = D ~ L1 + L2,
      censor_model = C ~ L1 + L2,
      ref_idx = 1,
      int_descript = c("Static Intervention",
                       "Dynamic Intervention"),
      intervention1.A1 = list(static(3)),
      intervention1.A2 = list(static(1)),
      intervention2.A1 = list(case_when(data$L2 < 0 ~ 1,
                                        data$L2 >= 0 & data$L2 < 2 ~ 2,
                                        T ~ 3)),
      intervention2.A2 = list(dynamic("compare", "L1", "=", 0))
    )$summary)
    rownames(out) <- c()
    out <- as.vector(out)
    expect_equal(
      out,
      as.vector(c(0.26947, NA, NA, 0.26967, 0.22639, 0.29826, 1.19117, 1.00000, 1.31746,
                  0.04328, 0.00000, 0.07187))
    )})

test_that(
  "check classical stratified ICE direct effect with intervention 2 as reference",
  {
    out <- as.matrix(ice(
      data = data,
      time_points = 4,
      id = "id",
      time_name = "t0",
      outcome_name = "Y",
      censor_name = "C",
      compevent_name = "D",
      estimator = strat(hazard = F),
      comp_effect = 0,
      outcome_model = Y ~ L1 + L2,
      competing_model = D ~ L1 + L2,
      censor_model = C ~ L1 + L2,
      ref_idx = 2,
      int_descript = c("Static Intervention",
                       "Dynamic Intervention"),
      intervention1.A1 = list(static(3)),
      intervention1.A2 = list(static(1)),
      intervention2.A1 = list(case_when(data$L2 < 0 ~ 1,
                                        data$L2 >= 0 & data$L2 < 2 ~ 2,
                                        T ~ 3)),
      intervention2.A2 = list(dynamic("compare", "L1", "=", 0))
    )$summary)
    rownames(out) <- c()
    out <- as.vector(out)
    expect_equal(
      out,
      as.vector(c(0.26947, NA, NA, 0.26967, 0.22639, 0.29826, 0.90414, 0.75904, 1.00000,
                  -0.02859, -0.07187, 0.00000))
    )})

test_that(
  "check classical pooled ICE direct effect with intervention times options",
  {
    out <- as.matrix(ice(
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
      intervention1.A1 = list(static(3), 1:2),
      intervention1.A2 = list(static(1), 1:2),
      intervention2.A1 = list(case_when(data$L2 < 0 ~ 1,
                                        data$L2 >= 0 & data$L2 < 2 ~ 2,
                                        T ~ 3), 2:3),
      intervention2.A2 = list(dynamic("compare", "L1", "=", 0), 2:3)
    )$summary)
    rownames(out) <- c()
    out <- as.vector(out)
    expect_equal(
      out,
      as.vector(c(0.26947, NA, NA, 0.26874, 0.26011, 0.31459, 1.00000, 0.96788, 1.17063,
                  0.00000, -0.00863, 0.04585))
    )})

test_that(
  "check classical stratified ICE direct effect with intervention times options",
  {
    out <- as.matrix(ice(
      data = data,
      time_points = 4,
      id = "id",
      time_name = "t0",
      outcome_name = "Y",
      censor_name = "C",
      compevent_name = "D",
      estimator = strat(hazard = F),
      comp_effect = 0,
      outcome_model = Y ~ L1 + L2,
      competing_model = D ~ L1 + L2,
      censor_model = C ~ L1 + L2,
      ref_idx = 0,
      int_descript = c("Static Intervention",
                       "Dynamic Intervention"),
      intervention1.A1 = list(static(3), 1:2),
      intervention1.A2 = list(static(1), 1:2),
      intervention2.A1 = list(case_when(data$L2 < 0 ~ 1,
                                        data$L2 >= 0 & data$L2 < 2 ~ 2,
                                        T ~ 3), 2:3),
      intervention2.A2 = list(dynamic("compare", "L1", "=", 0), 2:3)
    )$summary)
    rownames(out) <- c()
    out <- as.vector(out)
    expect_equal(
      out,
      as.vector(c(0.26947, NA, NA, 0.26967, 0.25865, 0.32675, 1.00000, 0.95916, 1.21168,
                  0.00000, -0.01101, 0.05708))
    )})

test_that(
  "check classical pooled ICE direct effect with intervention times options and intervention 1 as reference",
  {
    out <- as.matrix(ice(
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
      ref_idx = 1,
      int_descript = c("Static Intervention",
                       "Dynamic Intervention"),
      intervention1.A1 = list(static(3), 1:2),
      intervention1.A2 = list(static(1), 1:2),
      intervention2.A1 = list(case_when(data$L2 < 0 ~ 1,
                                        data$L2 >= 0 & data$L2 < 2 ~ 2,
                                        T ~ 3), 2:3),
      intervention2.A2 = list(dynamic("compare", "L1", "=", 0), 2:3)
    )$summary)
    rownames(out) <- c()
    out <- as.vector(out)
    expect_equal(
      out,
      as.vector(c(0.26947, NA, NA, 0.26874, 0.26011, 0.31459, 1.03319, 1.00000, 1.20948,
                  0.00863, 0.00000, 0.05449))
    )})

test_that(
  "check classical pooled ICE direct effect with intervention times options and intervention 2 as reference",
  {
    out <- as.matrix(ice(
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
      ref_idx = 2,
      int_descript = c("Static Intervention",
                       "Dynamic Intervention"),
      intervention1.A1 = list(static(3), 1:2),
      intervention1.A2 = list(static(1), 1:2),
      intervention2.A1 = list(case_when(data$L2 < 0 ~ 1,
                                        data$L2 >= 0 & data$L2 < 2 ~ 2,
                                        T ~ 3), 2:3),
      intervention2.A2 = list(dynamic("compare", "L1", "=", 0), 2:3)
    )$summary)
    rownames(out) <- c()
    out <- as.vector(out)
    expect_equal(
      out,
      as.vector(c(0.26947, NA, NA, 0.26874, 0.26011, 0.31459, 0.85424, 0.82680, 1.00000,
                  -0.04585, -0.05449, 0.00000))
    )})

test_that(
  "check hazard extended pooled ICE direct effect with intervention times options and intervention 1 as reference",
  {
    out <- as.matrix(ice(
      data = data,
      time_points = 4,
      id = "id",
      time_name = "t0",
      outcome_name = "Y",
      censor_name = "C",
      compevent_name = "D",
      estimator = pool(hazard = T),
      comp_effect = 0,
      outcome_model = Y ~ L1 + L2 + A1 + A2,
      competing_model = D ~ L1 + L2 + A1 + A2,
      censor_model = C ~ L1 + L2 + A1 + A2,
      ref_idx = 1,
      int_descript = c("Static Intervention",
                       "Dynamic Intervention"),
      intervention1.A1 = list(static(3), 1:2),
      intervention1.A2 = list(static(1), 1:2),
      intervention2.A1 = list(case_when(data$L2 < 0 ~ 1,
                                        data$L2 >= 0 & data$L2 < 2 ~ 2,
                                        T ~ 3), 2:3),
      intervention2.A2 = list(dynamic("compare", "L1", "=", 0), 2:3)
    )$summary)
    rownames(out) <- c()
    out <- as.vector(out)
    expect_equal(
      out,
      as.vector(c(0.26947, NA, NA, 0.26809, 0.25660, 0.31748, 1.04476, 1.00000, 1.23723,
                  0.01149, 0.00000, 0.06087))
    )})

test_that(
  "check hazard extended pooled ICE direct effect with intervention times options and intervention 2 as reference",
  {
    out <- as.matrix(ice(
      data = data,
      time_points = 4,
      id = "id",
      time_name = "t0",
      outcome_name = "Y",
      censor_name = "C",
      compevent_name = "D",
      estimator = pool(hazard = T),
      comp_effect = 0,
      outcome_model = Y ~ L1 + L2 + A1 + A2,
      competing_model = D ~ L1 + L2 + A1 + A2,
      censor_model = C ~ L1 + L2 + A1 + A2,
      ref_idx = 2,
      int_descript = c("Static Intervention",
                       "Dynamic Intervention"),
      intervention1.A1 = list(static(3), 1:2),
      intervention1.A2 = list(static(1), 1:2),
      intervention2.A1 = list(case_when(data$L2 < 0 ~ 1,
                                        data$L2 >= 0 & data$L2 < 2 ~ 2,
                                        T ~ 3), 2:3),
      intervention2.A2 = list(dynamic("compare", "L1", "=", 0), 2:3)
    )$summary)
    rownames(out) <- c()
    out <- as.vector(out)
    expect_equal(
      out,
      as.vector(c(0.26947, NA, NA, 0.26809, 0.25660, 0.31748, 0.84444, 0.80826, 1.00000,
                  -0.04939, -0.06087, 0.00000))
    )})

test_that(
  "check classical stratified ICE direct effect with intervention times options",
  {
    out <- as.matrix(ice(
      data = data,
      time_points = 4,
      id = "id",
      time_name = "t0",
      outcome_name = "Y",
      censor_name = "C",
      compevent_name = "D",
      estimator = strat(hazard = F),
      comp_effect = 0,
      outcome_model = Y ~ L1 + L2,
      competing_model = D ~ L1 + L2,
      censor_model = C ~ L1 + L2,
      ref_idx = 0,
      int_descript = c("Static Intervention",
                       "Dynamic Intervention"),
      intervention1.A1 = list(static(3), 1:2),
      intervention1.A2 = list(static(1), 1:2),
      intervention2.A1 = list(case_when(data$L2 < 0 ~ 1,
                                        data$L2 >= 0 & data$L2 < 2 ~ 2,
                                        T ~ 3), 2:3),
      intervention2.A2 = list(dynamic("compare", "L1", "=", 0), 2:3)
    )$summary)
    rownames(out) <- c()
    out <- as.vector(out)
    expect_equal(
      out,
      as.vector(c(0.26947, NA, NA, 0.26967, 0.25865, 0.32675, 1.00000, 0.95916, 1.21168,
                  0.00000, -0.01101, 0.05708))
    )})

test_that(
  "check classical stratified ICE direct effect with intervention times options and intervention 1 as reference",
  {
    out <- as.matrix(ice(
      data = data,
      time_points = 4,
      id = "id",
      time_name = "t0",
      outcome_name = "Y",
      censor_name = "C",
      compevent_name = "D",
      estimator = strat(hazard = F),
      comp_effect = 0,
      outcome_model = Y ~ L1 + L2,
      competing_model = D ~ L1 + L2,
      censor_model = C ~ L1 + L2,
      ref_idx = 1,
      int_descript = c("Static Intervention",
                       "Dynamic Intervention"),
      intervention1.A1 = list(static(3), 1:2),
      intervention1.A2 = list(static(1), 1:2),
      intervention2.A1 = list(case_when(data$L2 < 0 ~ 1,
                                        data$L2 >= 0 & data$L2 < 2 ~ 2,
                                        T ~ 3), 2:3),
      intervention2.A2 = list(dynamic("compare", "L1", "=", 0), 2:3)
    )$summary)
    rownames(out) <- c()
    out <- as.vector(out)
    expect_equal(
      out,
      as.vector(c(0.26947, NA, NA, 0.26967, 0.25865, 0.32675, 1.04258, 1.00000, 1.26328,
                  0.01101, 0.00000, 0.06810))
    )})


test_that(
  "check classical stratified ICE direct effect with intervention times options and intervention 2 as reference",
  {
    out <- as.matrix(ice(
      data = data,
      time_points = 4,
      id = "id",
      time_name = "t0",
      outcome_name = "Y",
      censor_name = "C",
      compevent_name = "D",
      estimator = strat(hazard = F),
      comp_effect = 0,
      outcome_model = Y ~ L1 + L2,
      competing_model = D ~ L1 + L2,
      censor_model = C ~ L1 + L2,
      ref_idx = 2,
      int_descript = c("Static Intervention",
                       "Dynamic Intervention"),
      intervention1.A1 = list(static(3), 1:2),
      intervention1.A2 = list(static(1), 1:2),
      intervention2.A1 = list(case_when(data$L2 < 0 ~ 1,
                                        data$L2 >= 0 & data$L2 < 2 ~ 2,
                                        T ~ 3), 2:3),
      intervention2.A2 = list(dynamic("compare", "L1", "=", 0), 2:3)
    )$summary)
    rownames(out) <- c()
    out <- as.vector(out)
    expect_equal(
      out,
      as.vector(c(0.26947, NA, NA, 0.26967, 0.25865, 0.32675, 0.82530, 0.79159, 1.00000,
                  -0.05708, -0.06810, 0.00000))
    )})

test_that(
  "check hazard extended stratified ICE direct effect with intervention times options and intervention 1 as reference",
  {
    out <- as.matrix(ice(
      data = data,
      time_points = 4,
      id = "id",
      time_name = "t0",
      outcome_name = "Y",
      censor_name = "C",
      compevent_name = "D",
      estimator = strat(hazard = T),
      comp_effect = 0,
      outcome_model = Y ~ L1 + L2,
      competing_model = D ~ L1 + L2,
      censor_model = C ~ L1 + L2,
      ref_idx = 1,
      int_descript = c("Static Intervention",
                       "Dynamic Intervention"),
      intervention1.A1 = list(static(3), 1:2),
      intervention1.A2 = list(static(1), 1:2),
      intervention2.A1 = list(case_when(data$L2 < 0 ~ 1,
                                        data$L2 >= 0 & data$L2 < 2 ~ 2,
                                        T ~ 3), 2:3),
      intervention2.A2 = list(dynamic("compare", "L1", "=", 0), 2:3)
    )$summary)
    rownames(out) <- c()
    out <- as.vector(out)
    expect_equal(
      out,
      as.vector(c(0.26947, NA, NA, 0.26936, 0.25894, 0.32713, 1.04024, 1.00000, 1.26335,
                  0.01042, 0.00000, 0.06819))
    )})

test_that(
  "check hazard extended stratified ICE direct effect with intervention times options and intervention 2 as reference",
  {
    out <- as.matrix(ice(
      data = data,
      time_points = 4,
      id = "id",
      time_name = "t0",
      outcome_name = "Y",
      censor_name = "C",
      compevent_name = "D",
      estimator = strat(hazard = T),
      comp_effect = 0,
      outcome_model = Y ~ L1 + L2,
      competing_model = D ~ L1 + L2,
      censor_model = C ~ L1 + L2,
      ref_idx = 2,
      int_descript = c("Static Intervention",
                       "Dynamic Intervention"),
      intervention1.A1 = list(static(3), 1:2),
      intervention1.A2 = list(static(1), 1:2),
      intervention2.A1 = list(case_when(data$L2 < 0 ~ 1,
                                        data$L2 >= 0 & data$L2 < 2 ~ 2,
                                        T ~ 3), 2:3),
      intervention2.A2 = list(dynamic("compare", "L1", "=", 0), 2:3)
    )$summary)
    rownames(out) <- c()
    out <- as.vector(out)
    expect_equal(
      out,
      as.vector(c(0.26947, NA, NA, 0.26936, 0.25894, 0.32713, 0.82340, 0.79155, 1.00000,
                  -0.05777, -0.06819, 0.00000))
    )})

test_that(
  "check weighted ICE direct effect with intervention times options and intervention 1 as reference",
  {
    out <- as.matrix(ice(
      data = data,
      time_points = 4,
      id = "id",
      time_name = "t0",
      outcome_name = "Y",
      censor_name = "C",
      compevent_name = "D",
      estimator = weight(list(A1 ~ L1 + L2, A2 ~ L1 + L2)),
      comp_effect = 0,
      outcome_model = Y ~ L1 + L2,
      competing_model = D ~ L1 + L2,
      censor_model = C ~ L1 + L2,
      ref_idx = 1,
      int_descript = c("Static Intervention",
                       "Dynamic Intervention"),
      intervention1.A1 = list(static(3), 1:2),
      intervention1.A2 = list(static(1), 1:2),
      intervention2.A1 = list(case_when(data$L2 < 0 ~ 1,
                                        data$L2 >= 0 & data$L2 < 2 ~ 2,
                                        T ~ 3), 2:3),
      intervention2.A2 = list(dynamic("compare", "L1", "=", 0), 2:3)
    )$summary)
    rownames(out) <- c()
    out <- as.vector(out)
    expect_equal(
      out,
      as.vector(c(0.26947, NA, NA, 0.26162, 0.25775, 0.32653, 1.01505, 1.00000, 1.26685,
                  0.00388, 0.00000, 0.06878))
    )})

test_that(
  "check weighted ICE direct effect with intervention times options and intervention 2 as reference",
  {
    out <- as.matrix(ice(
      data = data,
      time_points = 4,
      id = "id",
      time_name = "t0",
      outcome_name = "Y",
      censor_name = "C",
      compevent_name = "D",
      estimator = weight(list(A1 ~ L1 + L2, A2 ~ L1 + L2)),
      comp_effect = 0,
      outcome_model = Y ~ L1 + L2,
      competing_model = D ~ L1 + L2,
      censor_model = C ~ L1 + L2,
      ref_idx = 2,
      int_descript = c("Static Intervention",
                       "Dynamic Intervention"),
      intervention1.A1 = list(static(3), 1:2),
      intervention1.A2 = list(static(1), 1:2),
      intervention2.A1 = list(case_when(data$L2 < 0 ~ 1,
                                        data$L2 >= 0 & data$L2 < 2 ~ 2,
                                        T ~ 3), 2:3),
      intervention2.A2 = list(dynamic("compare", "L1", "=", 0), 2:3)
    )$summary)
    rownames(out) <- c()
    out <- as.vector(out)
    expect_equal(
      out,
      as.vector(c(0.26947, NA, NA, 0.26162, 0.25775, 0.32653, 0.80124, 0.78936, 1.00000,
                  -0.06490, -0.06878, 0.00000))
    )})


test_that(
  "check classical pooled ICE direct effect bootstrap with empirical percentile",
  {
    out <- as.matrix(ice(
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
      nsamples = 10, parallel = F, ci_method = "percentile"
    )$summary)
    rownames(out) <- c()
    colnames(out) <- c()
    dimnames(out) <- c()
    out <- as.matrix(out)
    expect_equal(
      out, 
      matrix(c(0.26947, 0.26874, 1.00000, 0.00000, 0.00400, 0.00000, 0.0000, 0.26763, 0.26985, 1.00000, 1.00000, 0.0000, 0.0000,
               NA, 0.24935, 0.92784, -0.01939, 0.00884, 0.02187, 0.0058, 0.24697, 0.25172, 0.22815, 0.27054, 0.2494, 0.2493),
             byrow = T, nrow = 2)
    )})

test_that(
  "check hazard-extended pooled ICE direct effect bootstrap with empirical percentile",
  {
    out <- as.matrix(ice(
      data = data,
      time_points = 4,
      id = "id",
      time_name = "t0",
      outcome_name = "Y",
      censor_name = "C",
      compevent_name = "D",
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
      nsamples = 10, parallel = F, ci_method = "percentile"
    )$summary)
    rownames(out) <- c()
    colnames(out) <- c()
    dimnames(out) <- c()
    out <- as.matrix(out)
    expect_equal(
      out,
      matrix(c(0.26947, 0.26809, 1.00000, 0.00000, 0.00388, 0.00000, 0.0000, 0.26702, 0.26915, 1.00000, 1.00000, 0.00000, 0.00000,
               NA, 0.23528, 0.87763, -0.03281, 0.00663, 0.01926, 0.0052, 0.23363, 0.23693, 0.21781, 0.25276, 0.23542, 0.23515), byrow = T, nrow = 2)
    )})

test_that(
  "check classical pooled ICE total effect bootstrap with empirical percentile",
  {
    out <- as.matrix(ice(
      data = data,
      time_points = 4,
      id = "id",
      time_name = "t0",
      outcome_name = "Y",
      censor_name = "C",
      compevent_name = "D",
      estimator = pool(hazard = F),
      comp_effect = 1,
      outcome_model = Y ~ L1 + L2 + A1 + A2,
      competing_model = D ~ L1 + L2 + A1 + A2,
      censor_model = C ~ L1 + L2 + A1 + A2,
      ref_idx = 0,
      int_descript = c("Static Intervention",
                       "Dynamic Intervention"),
      intervention1.A1 = list(static(3)),
      intervention1.A2 = list(static(1)),
      nsamples = 10, parallel = F, ci_method = "percentile"
    )$summary)
    rownames(out) <- c()
    colnames(out) <- c()
    dimnames(out) <- c()
    out <- as.matrix(out)
    expect_equal(
      out,
      matrix(c(0.22949, 0.22776, 1.00000, 0.00000, 0.00374, 0.00000, 0.00000, 0.22688, 0.22863, 1.00000, 1.00000, 0.00000, 0.00000,
        NA, 0.21449, 0.94175, -0.01327, 0.00778, 0.02127, 0.00478, 0.21269, 0.21628, 0.19356, 0.23542, 0.21451, 0.21447), byrow = T,
        nrow = 2)
    )})

test_that(
  "check hazard-extended pooled ICE total effect bootstrap with empirical percentile",
  {
    out <- as.matrix(ice(
      data = data,
      time_points = 4,
      id = "id",
      time_name = "t0",
      outcome_name = "Y",
      censor_name = "C",
      compevent_name = "D",
      estimator = pool(hazard = T),
      comp_effect = 1,
      outcome_model = Y ~ L1 + L2 + A1 + A2,
      competing_model = D ~ L1 + L2 + A1 + A2,
      censor_model = C ~ L1 + L2 + A1 + A2,
      ref_idx = 0,
      int_descript = c("Static Intervention",
                       "Dynamic Intervention"),
      intervention1.A1 = list(static(3)),
      intervention1.A2 = list(static(1)),
      nsamples = 10, parallel = F, ci_method = "percentile"
    )$summary)
    rownames(out) <- c()
    colnames(out) <- c()
    dimnames(out) <- c()
    out <- as.matrix(out)
    expect_equal(
      out,
      matrix(c(0.22949, 0.22582, 1.00000, 0.00000, 0.00358, 0.00000, 0.00000, 0.22499, 0.22665, 1.00000, 1.00000, 0.00000, 0.00000,
               NA, 0.20231, 0.89588, -0.02351, 0.00574, 0.01896, 0.00431, 0.20108, 0.20354, 0.18473, 0.21989, 0.20238, 0.20224), byrow = T,
             nrow = 2)
    )})

test_that(
  "check classical stratified ICE total effect bootstrap with empirical percentile",
  {
    out <- as.matrix(ice(
      data = data,
      time_points = 4,
      id = "id",
      time_name = "t0",
      outcome_name = "Y",
      censor_name = "C",
      compevent_name = "D",
      estimator = strat(hazard = F),
      comp_effect = 1,
      outcome_model = Y ~ L1 + L2,
      competing_model = D ~ L1 + L2,
      censor_model = C ~ L1 + L2,
      ref_idx = 0,
      int_descript = c("Static Intervention",
                       "Dynamic Intervention"),
      intervention1.A1 = list(static(3)),
      intervention1.A2 = list(static(1)),
      nsamples = 10, parallel = F, ci_method = "percentile"
    )$summary)
    rownames(out) <- c()
    colnames(out) <- c()
    dimnames(out) <- c()
    out <- as.matrix(out)
    expect_equal(
      out,
      matrix(c(0.22949, 0.23119, 1.00000, 0.00000, 0.00383, 0.0000, 0.00000, 0.23028, 0.23210, 1.00000, 1.00000, 0.00000, 0.00000,
               NA, 0.19974, 0.86397, -0.03145, 0.02783, 0.1094, 0.02531, 0.19272, 0.20677, 0.08314, 0.31635, 0.19935, 0.20014), byrow = T,
             nrow = 2)
    )})

test_that(
  "check classical pooled ICE total effect bootstrap with normal percentile",
  {
    out <- as.matrix(ice(
      data = data,
      time_points = 4,
      id = "id",
      time_name = "t0",
      outcome_name = "Y",
      censor_name = "C",
      compevent_name = "D",
      estimator = pool(hazard = F),
      comp_effect = 1,
      outcome_model = Y ~ L1 + L2 + A1 + A2,
      competing_model = D ~ L1 + L2 + A1 + A2,
      censor_model = C ~ L1 + L2 + A1 + A2,
      ref_idx = 0,
      int_descript = c("Static Intervention",
                       "Dynamic Intervention"),
      intervention1.A1 = list(static(3)),
      intervention1.A2 = list(static(1)),
      nsamples = 10, parallel = F, ci_method = "normal"
    )$summary)
    rownames(out) <- c()
    colnames(out) <- c()
    dimnames(out) <- c()
    out <- as.matrix(out)
    expect_equal(
      out, 
      matrix(c(0.22949, 0.22776, 1.00000, 0.00000, 0.00374, 0.00000, 0.00000, 0.22042, 0.23509, 1.00000, 1.00000, 0.00000, 0.00000,
               NA, 0.21449, 0.94175, -0.01327, 0.00778, 0.02127, 0.00478, 0.19924, 0.22974, 0.17279, 0.25618, 0.20512, 0.22386), byrow = T,
             nrow = 2)
    )})

test_that(
  "check hazard-extended pooled ICE total effect bootstrap with normal percentile",
  {
    out <- as.matrix(ice(
      data = data,
      time_points = 4,
      id = "id",
      time_name = "t0",
      outcome_name = "Y",
      censor_name = "C",
      compevent_name = "D",
      estimator = pool(hazard = T),
      comp_effect = 1,
      outcome_model = Y ~ L1 + L2 + A1 + A2,
      competing_model = D ~ L1 + L2 + A1 + A2,
      censor_model = C ~ L1 + L2 + A1 + A2,
      ref_idx = 0,
      int_descript = c("Static Intervention",
                       "Dynamic Intervention"),
      intervention1.A1 = list(static(3)),
      intervention1.A2 = list(static(1)),
      nsamples = 10, parallel = F, ci_method = "normal"
    )$summary)
    rownames(out) <- c()
    colnames(out) <- c()
    dimnames(out) <- c()
    out <- as.matrix(out)
    expect_equal(
      out,
      matrix(c(0.22949, 0.22582, 1.00000, 0.00000, 0.00358, 0.00000, 0.00000, 0.21881, 0.23283, 1.00000, 1.00000, 0.00000, 0.00000,
               NA, 0.20231, 0.89588, -0.02351, 0.00574, 0.01896, 0.00431, 0.19106, 0.21356, 0.16515, 0.23947, 0.19386, 0.21076), byrow = T,
             nrow = 2)
    )})

test_that(
  "check hazard-extended stratified ICE total effect bootstrap with normal percentile",
  {
    out <- as.matrix(ice(
      data = data,
      time_points = 4,
      id = "id",
      time_name = "t0",
      outcome_name = "Y",
      censor_name = "C",
      compevent_name = "D",
      estimator = strat(hazard = T),
      comp_effect = 1,
      outcome_model = Y ~ L1 + L2,
      competing_model = D ~ L1 + L2,
      censor_model = C ~ L1 + L2,
      ref_idx = 0,
      int_descript = c("Static Intervention",
                       "Dynamic Intervention"),
      intervention1.A1 = list(static(3)),
      intervention1.A2 = list(static(1)),
      nsamples = 10, parallel = F, ci_method = "normal"
    )$summary)
    rownames(out) <- c()
    colnames(out) <- c()
    dimnames(out) <- c()
    out <- as.matrix(out)
    expect_equal(
      out,
      matrix(c(0.22949, 0.24340, 1.00000, 0.00000, 0.00398, 0.00000, 0.0000, 0.23561, 0.25119, 1.00000, 1.00000, 0.00000, 0.00000,
               NA, 0.20589, 0.84587, -0.03751, 0.02949, 0.1102, 0.02684, 0.14810, 0.26368, -0.01011, 0.42188, 0.15329, 0.25849), byrow = T,
             nrow = 2)
    )})

test_that(
  "check weighted ICE total effect bootstrap with normal percentile",
  {
    out <- as.matrix(ice(
      data = data,
      time_points = 4,
      id = "id",
      time_name = "t0",
      outcome_name = "Y",
      censor_name = "C",
      compevent_name = "D",
      estimator = weight(list(A1 ~ L1 + L2, A2 ~ L1 + L2)),
      comp_effect = 1,
      outcome_model = Y ~ L1 + L2,
      competing_model = D ~ L1 + L2,
      censor_model = C ~ L1 + L2,
      ref_idx = 0,
      int_descript = c("Static Intervention",
                       "Dynamic Intervention"),
      intervention1.A1 = list(static(3)),
      intervention1.A2 = list(static(1)),
      nsamples = 10, parallel = F, ci_method = "normal"
    )$summary)
    rownames(out) <- c()
    colnames(out) <- c()
    dimnames(out) <- c()
    out <- as.matrix(out)
    expect_equal(
      out,
      matrix(c(0.22949, 0.22683, 1.00000, 0.00000, 0.00379, 0.0000, 0.00000, 0.21939, 0.23426, 1.00000, 1.00000, 0.00000, 0.00000,
               NA, 0.19755, 0.87094, -0.02927, 0.02652, 0.1051, 0.02383, 0.14557, 0.24954, -0.00845, 0.40355, 0.15085, 0.24426), byrow = T,
             nrow = 2)
    )})

