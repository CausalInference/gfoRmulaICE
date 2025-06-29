# gfoRmulaICE: Parametric Iterative Conditional Expectation G-Formula
  <!-- badges: start -->
  ![CRAN Version](https://www.r-pkg.org/badges/version/gfoRmulaICE)
  [![Lifecycle: stable](https://img.shields.io/badge/lifecycle-stable-brightgreen.svg)](https://lifecycle.r-lib.org/articles/stages.html#stable)
  [![R-CMD-check](https://github.com/CausalInference/gfoRmulaICE/actions/workflows/R-CMD-check.yaml/badge.svg)](https://github.com/CausalInference/gfoRmulaICE/actions/workflows/R-CMD-check.yaml)
  ![Downloads](https://cranlogs.r-pkg.org/badges/grand-total/gfoRmulaICE)
  <!-- badges: end -->

## How to install 

Installation from CRAN:

``` r
install.packages("gfoRmulaICE")
```
Installation from GitHub repo:

``` r
devtools::install_github("CausalInference/gfoRmulaICE")
```

## How to use

Please see the following for a few examples.

### Example 1: Dynamic Intervention

In this example, we illustrate how to specify dynamic intervention. We consider the following interventions, which are applied to all time points.
* Intervention 1 on A2: at time t, if L1 = 0, then treat; otherwise, not treat. 
* Intervention 2 on A2: never treat upon until L1 = 0, after which follows always treat.
* Intervention 3 on A2: never treat upon until L1 = 0, after which follows natural course.

We use classical pooled ICE estimator and natural course as the reference intervention. We estimate variance using bootstrap with 1000 replicates, normal quantile, and parallel computing. 

``` r
library(gfoRmulaICE)

ice_fit1 <- ice(data = gfoRmulaICE::compData, 
                time_points = 4, 
                id = "id", 
                time_name = "t0",
                censor_name = "C", 
                outcome_name = "Y",
                compevent_name = "D",
                comp_effect = 0,
                outcome_model = Y ~ L1 + L2 + A1 + A2, 
                censor_model = C ~ L1 + L2 + A1 + A2,
                ref_idx = 0,
                estimator = pool(hazard = F),
                nsamples = 1000, ci_method = "percentile",
                parallel = T, ncores = 5,
                int_descript = c("Dynamic Intervention 1", "Dynamic Intervention 2", "Dynamic Intervention 3"),
                intervention1.A2 = list(dynamic("L1 == 0", static(0), static(1))),
                intervention2.A2 = list(dynamic("L1 == 0", static(0), static(1), absorb = T)),
                intervention3.A2 = list(dynamic("L1 == 0", static(0), natural_course()))
)

summary(ice_fit1)

plot(ice_fit1)
```

### Example 2: Built-in Interventions

In this example, we illustrate the available intervention functions within the package. We consider the following interventions and apply to all time points.
* Intervention 1 on A1: always treat with value 3.
* Intervention 1 on A2: always treat with value 1.
* Intervention 2 on L2: when the natural value of L2 at time t is lower than -3, set its value to -3. Otherwise, do not intervene.
* Intervention 3 on A2: dynamic intervention (treat when L1 = 0) with uniform grace period of 2 periods.

We use classical pooled ICE estimator and natural course as the reference intervention. We estimate variance using bootstrap with 1000 replicates, normal quantile, and parallel computing.

``` r
ice_fit2 <- ice(data = gfoRmulaICE::compData, 
                time_points = 4, 
                id = "id", 
                time_name = "t0",
                censor_name = "C", 
                outcome_name = "Y",
                compevent_name = "D",
                comp_effect = 0,
                outcome_model = Y ~ L1 + L2 + A1 + A2, 
                censor_model = C ~ L1 + L2 + A1 + A2,
                ref_idx = 0,
                estimator = pool(hazard = F),
                nsamples = 1000, 
                ci_method = "percentile",
                parallel = T, 
                ncores = 5,
                int_descript = c("Static Intervention", "Threshold Intervention", "Dynamic Intervention with Grace Period"),
                intervention1.A1 = list(static(3)),
                intervention1.A2 = list(static(1)),
                intervention2.L2 = list(threshold(-3, Inf)),
                intervention3.A2 = list(grace_period("uniform", 2, "L1 == 0"))
)

summary(ice_fit2)

plot(ice_fit2)
```

### Example 3: User-defined Intervention

In this example, we illustrate how to specify user-defined intervention. We consider the following interventions and apply to all time points.
* Intervention 1 on A1: at time t, if L2 < 0, then assign 1; if 0 <= L2 < 2, then assign 2; otherwise, assign 3.

We use classical pooled ICE estimator and natural course as the reference intervention. We estimate variance using bootstrap with 1000 replicates and percentile quantile.


``` r
dynamic_cat <- case_when(gfoRmulaICE::compData$L2 < 0 ~ 1,
                         gfoRmulaICE::compData$L2 >= 0 & gfoRmulaICE::compData$L2 < 2 ~ 2, 
                         T ~ 3)

ice_fit3 <- ice(data = gfoRmulaICE::compData, 
                time_points = 4, 
                id = "id", 
                time_name = "t0",
                censor_name = "C", 
                outcome_name = "Y",
                compevent_name = "D",
                comp_effect = 0,
                outcome_model = Y ~ L1 + L2 + A1 + A2, 
                censor_model = C ~ L1 + L2 + A1 + A2,
                competing_model = D ~ L1 + L2 + A1 + A2,
                ref_idx = 0,
                estimator = pool(hazard = F),
                nsamples = 1000, 
                ci_method = "percentile",
                parallel = T, 
                ncores = 5,
                int_descript = c("Dynamic Intervention"),
                intervention1.A1 = list(dynamic_cat)
)

summary(ice_fit3)

plot(ice_fit3)

```

### Example 4: Different ICE Estimators

In this example, we illustrate all available ICE estimators in the package. We consider the following interventions and apply to all time points.
* Intervention 1 on A1: always treat with value 3.
* Intervention 1 on A2: always treat with value 1.
* Intervention 2 on A1: at time t, if L2 < 0, then assign 1; if 0 <= L2 < 2, then assign 2; otherwise, assign 3.
* Intervention 2 on A2: at time t, if L1 = 0, then treat; otherwise, not treat. 

We use natural course as the reference intervention and estimate variance using bootstrap with 1000 replicates and percentile quantile.

#### Classical pooled ICE

``` r
dynamic_cat <- case_when(gfoRmulaICE::compData$L2 < 0 ~ 1,
                         gfoRmulaICE::compData$L2 >= 0 & gfoRmulaICE::compData$L2 < 2 ~ 2, 
                         T ~ 3)

ice_fit4a <- ice(data = gfoRmulaICE::compData, 
                time_points = 4, 
                id = "id", 
                time_name = "t0",
                censor_name = "C", 
                outcome_name = "Y",
                compevent_name = "D",
                comp_effect = 0,
                outcome_model = Y ~ L1 + L2 + A1 + A2, 
                censor_model = C ~ L1 + L2 + A1 + A2,
                competing_model = D ~ L1 + L2 + A1 + A2,
                ref_idx = 0,
                estimator = pool(hazard = F),
                nsamples = 1000, 
                ci_method = "percentile",
                parallel = T, 
                ncores = 5,
                int_descript = c("Static Intervention: classical pooled ICE", 
                                 "Dynamic Intervention: classical pooled ICE"),
                intervention1.A1 = list(static(3)),
                intervention1.A2 = list(static(1)),
                intervention2.A1 = list(dynamic_cat),
                intervention2.A2 = list(dynamic("L1 == 0", static(0), static(1)))
)

summary(ice_fit4a)

plot(ice_fit4a)

```

#### Hazard-based pooled ICE

In case where the hazard model is the same as the outcome model, there is no need to specify hazard_model argument.

``` r
dynamic_cat <- case_when(gfoRmulaICE::compData$L2 < 0 ~ 1,
                         gfoRmulaICE::compData$L2 >= 0 & gfoRmulaICE::compData$L2 < 2 ~ 2, 
                         T ~ 3)

ice_fit4b <- ice(data = gfoRmulaICE::compData, 
                time_points = 4,
                id = "id", 
                time_name = "t0",
                censor_name = "C", 
                outcome_name = "Y",
                compevent_name = "D",
                comp_effect = 0,
                outcome_model = Y ~ L1 + L2 + A1 + A2, 
                censor_model = C ~ L1 + L2 + A1 + A2,
                competing_model = D ~ L1 + L2 + A1 + A2,
                ref_idx = 0,
                estimator = pool(hazard = T),
                nsamples = 1000, 
                ci_method = "percentile",
                parallel = T, 
                ncores = 5,
                int_descript = c("Static Intervention: hazard-based pooled ICE", 
                                 "Dynamic Intervention: hazard-based pooled ICE"),
                intervention1.A1 = list(static(3)),
                intervention1.A2 = list(static(1)),
                intervention2.A1 = list(dynamic_cat),
                intervention2.A2 = list(dynamic("L1 == 0", static(0), static(1)))
)

summary(ice_fit4b)

plot(ice_fit4b)

```

#### Classical stratified ICE

``` r
dynamic_cat <- case_when(gfoRmulaICE::compData$L2 < 0 ~ 1,
                         gfoRmulaICE::compData$L2 >= 0 & gfoRmulaICE::compData$L2 < 2 ~ 2, 
                         T ~ 3)

ice_fit4c <- ice(data = gfoRmulaICE::compData, 
                time_points = 4,
                id = "id", 
                time_name = "t0",
                censor_name = "C", 
                outcome_name = "Y",
                compevent_name = "D",
                comp_effect = 0,
                outcome_model = Y ~ L1 + L2, 
                censor_model = C ~ L1 + L2,
                ref_idx = 0,
                estimator = strat(hazard = F),
                nsamples = 1000, 
                ci_method = "percentile",
                parallel = T, 
                ncores = 5,
                int_descript = c("Static Intervention: classical stratified ICE", 
                                 "Dynamic Intervention: classical stratified ICE"),
                intervention1.A1 = list(static(3)),
                intervention1.A2 = list(static(1)),
                intervention2.A1 = list(dynamic_cat),
                intervention2.A2 = list(dynamic("L1 == 0", static(0), static(1)))
)

summary(ice_fit4c)

plot(ice_fit4c)

```

#### Hazard-based stratified ICE

``` r
dynamic_cat <- case_when(gfoRmulaICE::compData$L2 < 0 ~ 1,
                         gfoRmulaICE::compData$L2 >= 0 & gfoRmulaICE::compData$L2 < 2 ~ 2, 
                         T ~ 3)

ice_fit4d <- ice(data = gfoRmulaICE::compData, 
                time_points = 4, 
                id = "id", 
                time_name = "t0",
                censor_name = "C", 
                outcome_name = "Y",
                compevent_name = "D",
                comp_effect = 0,
                outcome_model = Y ~ L1 + L2, 
                censor_model = C ~ L1 + L2,
                competing_model = D ~ L1 + L2,
                hazard_model = Y ~ L1,
                ref_idx = 0,
                estimator = strat(hazard = T),
                nsamples = 1000, 
                ci_method = "percentile",
                parallel = T, 
                ncores = 5,
                int_descript = c("Static Intervention: hazard-based stratified ICE", 
                                 "Dynamic Intervention: hazard-based stratified ICE"),
                intervention1.A1 = list(static(3)),
                intervention1.A2 = list(static(1)),
                intervention2.A1 = list(dynamic_cat),
                intervention2.A2 = list(dynamic("L1 == 0", static(0), static(1)))
)

summary(ice_fit4d)

plot(ice_fit4d)

```

#### Doubly robust ICE

``` r
dynamic_cat <- case_when(gfoRmulaICE::compData$L2 < 0 ~ 1,
                         gfoRmulaICE::compData$L2 >= 0 & gfoRmulaICE::compData$L2 < 2 ~ 2, 
                         T ~ 3)

ice_fit4e <- ice(data = gfoRmulaICE::compData, 
                time_points = 4, 
                id = "id", 
                time_name = "t0",
                censor_name = "C", 
                outcome_name = "Y",
                compevent_name = "D",
                comp_effect = 0,
                outcome_model = Y ~ L1 + L2, 
                censor_model = C ~ L1 + L2,
                ref_idx = 0,
                estimator = weight(list(A1 ~ L1 + L2 + t0, A2 ~ L1 + L2 + t0)),
                nsamples = 1000, 
                ci_method = "percentile",
                parallel = T, 
                ncores = 5,
                int_descript = c("Static Intervention: doubly robust ICE", 
                                 "Dynamic Intervention: doubly robust ICE"),
                intervention1.A1 = list(static(3)),
                intervention1.A2 = list(static(1)),
                intervention2.A1 = list(dynamic_cat),
                intervention2.A2 = list(dynamic("L1 == 0", static(0), static(1)))
)

summary(ice_fit4e)

plot(ice_fit4e)

```

### Example 5: Options for hazard model in hazard-based pooled ICE

In this example, we illustrate available options for specifying hazard model in hazard-based ICE estimators. We consider the same interventions in Example 4.

We allow pooled-over-time hazard model with flexible terms of time variable and time-specific hazard model in hazard-based pooled ICE. Note, in hazard-based stratified ICE, only time-specific hazard model is allowed.

To specify pooled-over-time hazard model in hazard-based pooled ICE, we set the global_hazard argument to ```r{TRUE}```. To specify time-specific hazard model in hazard-based pooled ICE, we set the global_hazard argument to ```r{FALSE}```, which is the default value. Example 4b illustrates the time-specific hazard model. We illustrate the pooled-over-time option in this example.

``` r
library(splines)

dynamic_cat <- case_when(gfoRmulaICE::compData$L2 < 0 ~ 1,
                         gfoRmulaICE::compData$L2 >= 0 & gfoRmulaICE::compData$L2 < 2 ~ 2, 
                         T ~ 3)

ice_fit5 <- ice(data = gfoRmulaICE::compData, 
                time_points = 4, 
                id = "id", 
                time_name = "t0",
                censor_name = "C", 
                outcome_name = "Y",
                compevent_name = "D",
                comp_effect = 0,
                outcome_model = Y ~ L1 + L2 + A1 + A2, 
                censor_model = C ~ L1 + L2 + A1 + A2,
                competing_model = D ~ L1 + L2 + A1 + A2,
                hazard_model = Y ~ L1 + L2 + A1 + A2 + ns(t0, df = 2),
                global_hazard = T,
                ref_idx = 0,
                estimator = pool(hazard = T),
                nsamples = 1000, 
                ci_method = "percentile",
                parallel = T, 
                ncores = 5,
                int_descript = c("Static Intervention: pooled-over-time", 
                                 "Dynamic Intervention: pooled-over-time"),
                intervention1.A1 = list(static(3)),
                intervention1.A2 = list(static(1)),
                intervention2.A1 = list(dynamic_cat),
                intervention2.A2 = list(dynamic("L1 == 0", static(0), static(1)))
)

summary(ice_fit4b, ice_fit5)

plot(ice_fit4b, ice_fit5)

```

### Example 6: Intervention-specific models for stratified ICE and doubly robust ICE

We allow users to specify intervention-specific models in stratified ICE and doubly robust ICE. We consider the same analysis in Example 4d but with intervention-specific models.

In this example, the outcome model for intervention 1 is Y ~ L1; The outcome model for intervention 2: Y ~ L1 + L2; The competing model for intervention 1: D ~ L1 + L2; The competing model for intervention 2: D ~ L1.


``` r

dynamic_cat <- case_when(gfoRmulaICE::compData$L2 < 0 ~ 1,
                         gfoRmulaICE::compData$L2 >= 0 & gfoRmulaICE::compData$L2 < 2 ~ 2, 
                         T ~ 3)

ice_fit6 <- ice(data = gfoRmulaICE::compData,
                time_points = 4, 
                id = "id", 
                time_name = "t0",
                censor_name = "C", 
                outcome_name = "Y",
                compevent_name = "D",
                outcome_model = Y ~ L1, 
                censor_model = C ~ L1,
                competing_model = D ~ L1,
                comp_effect = 1,
                ref_idx = 0,
                estimator = strat(hazard = T),
                nsamples = 1000, 
                ci_method = "normal",
                parallel = T, 
                ncores = 5,
                int_descript = c("Static Intervention: intervention-specific model", 
                                 "Dynamic Intervention: intervention-specific model"),
                intervention1.A1 = list(static(3)),
                intervention1.A2 = list(static(1)),
                intervention2.A1 = list(dynamic_cat),
                intervention2.A2 = list(dynamic("L1 == 0", static(0), static(1))),
                outcomeModel.1 = Y ~ L1 + L2,
                compModel.2 = D ~ L1 + L2
)

summary(ice_fit4d, ice_fit6)

plot(ice_fit4d, ice_fit6)

```

### Example 7: Flexible Model Specification

#### Complicated terms in model statement

In this example, we illustrate that including polynomial, spline, and lagged terms in model statements is supported. We consider the same interventions as in Example 4. We illustrate using the pooled ICE estimator.

``` r

dynamic_cat <- case_when(gfoRmulaICE::compData$L2 < 0 ~ 1,
                         gfoRmulaICE::compData$L2 >= 0 & gfoRmulaICE::compData$L2 < 2 ~ 2, 
                         T ~ 3)

ice_fit7a <- ice(data = gfoRmulaICE::compData,
                time_points = 4, 
                id = "id", time_name = "t0",
                censor_name = "C", 
                outcome_name = "Y",
                compevent_name = "D",
                comp_effect = 0,
                outcome_model = Y ~ I(L1^2) + rcspline.eval(lag1_L2, knots = 1:3) + A1 + A2,
                censor_model = C ~ lag1_L1 + poly(L2, degree = 2) + A1 + A2,
                ref_idx = 0,
                estimator = pool(hazard = F),
                nsamples = 1000, 
                ci_method = "percentile",
                parallel = T, 
                ncores = 5,
                int_descript = c("Static Intervention", "Dynamic Intervention"),
                intervention1.A1 = list(static(3)),
                intervention1.A2 = list(static(1)),
                intervention2.A1 = list(dynamic_cat),
                intervention2.A2 = list(dynamic("L1 == 0", static(0), static(1)))
)

summary(ice_fit7a)

plot(ice_fit7a)
```

#### Use user-defined intervention as reference intervention

In this example, we use static intervention as the reference intervention. We consider the same interventions as in Example 4. We illustrate using the pooled ICE estimator.

``` r

dynamic_cat <- case_when(gfoRmulaICE::compData$L2 < 0 ~ 1,
                         gfoRmulaICE::compData$L2 >= 0 & gfoRmulaICE::compData$L2 < 2 ~ 2, 
                         T ~ 3)

ice_fit7b <- ice(data = gfoRmulaICE::compData,
                time_points = 4, 
                id = "id", 
                time_name = "t0",
                censor_name = "C", 
                outcome_name = "Y",
                compevent_name = "D",
                comp_effect = 0,
                outcome_model = Y ~ I(L1^2) + rcspline.eval(lag1_L2, knots = 1:3) + A1 + A2,
                censor_model = C ~ lag1_L1 + poly(L2, degree = 2) + A1 + A2,
                ref_idx = 1,
                estimator = pool(hazard = F),
                nsamples = 1000, 
                ci_method = "percentile",
                parallel = T, 
                ncores = 5,
                int_descript = c("Static Intervention", "Dynamic Intervention"),
                intervention1.A1 = list(static(3)),
                intervention1.A2 = list(static(1)),
                intervention2.A1 = list(dynamic_cat),
                intervention2.A2 = list(dynamic("L1 == 0", static(0), static(1)))
)

summary(ice_fit7b)

plot(ice_fit7b)

```
