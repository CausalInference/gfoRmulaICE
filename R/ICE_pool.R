#' ICE Pooling Over Treatment History Estimator including Both Classical and Hazard-Based Versions for A Single
#' Intervention Strategy.
#'
#' This function estimates the risk over time for survival outcome using the given observed data set following
#' a single user-defined intervention strategy by the parametric g-formula iterative conditional expectation (ICE)
#' estimator. This function provides the classical ICE pooling over treatment history method
#' and the hazard based ICE pooling over treatment history method.
#'
#'
#' @param data a data frame containing the observed data in long format.
#' @param K a numerical value that indicates the total number of time points.
#' @param id a character string indicating the ID variable name in \code{data}.
#' @param time_name a character string indicating the time variable name in \code{data}.
#' @param outcome_name a character string indicating the outcome variable name in \code{data}.
#' @param censor_name a character string indicating the censor variable name in \code{data}.
#' Could be \code{NULL} if there is no censoring variable in \code{data}. Default is \code{NULL}.
#' @param competing_name a character string indicating the competing variable name in \code{data}.
#' Could be \code{NULL} if there is no competing variable in \code{data}. Default is \code{NULL}.
#' @param total_effect a logical indicating how the competing event is handled. TRUE for total effect. FALSE for direct effect.
#' @param outcome_model a formula specifying the model statement for the outcome model.
#' @param censor_model a formula specifying the model statement for the censoring model in IP weighted natural course risk.
#' Could be \code{NULL} if censoring variable \code{censor_name} is \code{NULL}. Default is \code{NULL}.
#' @param competing_model a formula specifying the model statement for the competing model in hazard based ICE estimator.
#' Could be \code{NULL} for classical ICE estimator. Default is \code{NULL}.
#' @param interventions a list of functions indicating intervention strategy for the intervention variables in \code{intervention_names}.
#' Interventions could be specified using the following functions:
#' \itemize{
#' \item{Always treat:} {\code{static(value = 1)} }
#' \item{Never treat:} {\code{static(value = 0)} }
#' \item{Natural course:} {\code{natural_course()} }
#' \item{Dynamic treat:} {\code{dynamic(type, var, direction, value)} }
#' \item{Threshold Intervention:} {\code{threshold(value, var)}}
#' \item{Grace period:} {\code{grace_period(type, nperiod, var, value)} }
#' }
#' For more details, please read the corresponding documentation of each intervention function.
#' If the length of the \code{interventions} list is larger than 1, then the function will process them as multiple treatments simultaneously.
#' The user defined intervention functions must return the intervened value as a vector in the length of the number of rows in \code{data}.
#' The user defined intervention functions in \code{interventions} must be in corresponding order of \code{intervention_varnames}.
#' The same logic applies for the built in intervention functions.
#' @param intervention_names a list of character strings indicating the intervention variable names in correspondence with the intervention strategies in \code{interventions}.
#' Length of the vector must be in the same size of \code{interventions}.
#' @param compute_nc_risk a logical to indicate whether to compute observed natural course risk. TRUE for observed natural course risk computation. FALSE for no natural course risk computation. Default to be TRUE.
#' @param hazard_based a logical indicating whether to use hazard-based ICE estimator or classical ICE estimator. TRUE for hazard-based estimator. FALSE for classical estimator.
#' @param intervention_description a character string specifying a description of the implemented intervention strategy.
#'
#' @return A list containing the following components:
#' \item{gformula_risk_last_time} {The estimated risk for the specified intervention(s) at the last time step.}
#' \item{gformula_risk} {A table containing the estimated risk for the specified intervention(s) at each time step.}
#' \item{weight_h} {A list containing the inverse probability weighted estimates of the natural course risk. Could be \code{NULL} if \code{compute_nc_risk} is FALSE.}
#' \item{ipw_model} {A model object for the probability of censoring used in the inverse probability weighted estimate of the natural course risk. Could be \code{NULL} if \code{compute_nc_risk} is FALSE.}
#' \item{fit_models}{A list containing the fitted models for the outcome, the treatment (if applicable), and the competing event (if applicable).}
#' \item{model_summary}{A list containing the summary of the fitted models.}
#' \item{model_stderr}{A list containing the standard errors of the coefficients of the fitted models.}
#' \item{model_vcov}{A list containing the variance-covariance matrices of the parameters of the fitted models.}
#' \item{model_rmse}{A list containing the root mean square error (RMSE) values of the fitted models.}
#' @import reshape2 tidyverse
#' @export
#'
#' @references Wen L, Young JG, Robins JM, Hernán MA. Parametric g-formula implementations for causal survival analyses. Biometrics. 2021;77(2):740-753.
#' @references McGrath S, Lin V, Zhang Z, Petito LC, Logan RW, Hernán MA, and JG Young. gfoRmula: An R package for estimating the effects of sustained treatment strategies via the parametric g-formula. Patterns. 2020;1:100008.
#' @references Young JG, Herńan MA, Robins JM. Identification, estimation and approximation of risk under interventions that depend on the natural value of treatment using observational data. Epidemiologic Methods. 2014;3(1):1-19.
#' @references Young JG, Vatsa R, Murray EJ, Hernán MA. Interval-cohort designs and bias in the estimation of per-protocol effects: a simulation study. Trials. 2019;20(1):552.
#' @references Díaz, I, Williams, N, Hoffman, KL, & Schenck, EJ. Nonparametric causal effects based on longitudinal modified treatment policies. Journal of the American Statistical Association. 2021;118(542), 846–857.
#' @references Young JG, Stensrud MJ, Tchetgen Tchetgen EJ, Hernán MA. A causal framework for classical statistical estimands in failure-time settings with competing events. Statistics in medicine. 2020;39(8):1199-1236.
#' @references Wen L, Hernán MA, Robins JM. Multiply robust estimators of causal effects for survival outcomes. Scandinavian journal of statistics, theory and applications. 2022;49(3):1304-1328.
#' @references Haneuse S, Rotnitzky A. Estimation of the effect of interventions that modify the received treatment. Statistics in medicine. 2013;32(30):5260-5277.
#' @references McGrath S, Young JG, Hernán MA. Revisiting the g-null Paradox. Epidemiology. 2022;33(1):114-120.
#' 
#' @keywords internal
#'
#' @examples
#'
#' # Import data set
#' test_data <- readRDS("test_data_competing.rds")
#'
#' # Estimate risks for the always treat intervention
#' (i.e. constantly treat over all time points)
#' # using hazard extended pooled ICE
#' # competing event as direct effect
#'
#' ice_static <- ice_pool(data = test_data, K = 5, id = "id",
#' time_name = "t0", outcome_name = "Y",
#' competing_name = "D", censor_name = "C",
#' total_effect = F,
#' outcome_model = Y ~ L1 + A,
#' censor_model = C ~ L1 + A,
#' competing_model = D ~ L1 + A,
#' interventions = list(static(1)),
#' intervention_names = list("A"),
#' hazard_based = T,
#' intervention_description = "Always Treat")
#'
#' ice_static
#'
#' # Estimate risks for both the observed inverse probability weighted
#' # natural course risk and natural course strategy
#' # Classical pooled ICE, competing event as total effect
#'
#' ice_natural_course <- ice_pool(data = test_data, K = 5, id = "id",
#' time_name = "t0", outcome_name = "Y",
#' competing_name = "D", censor_name = "C",
#' total_effect = T,
#' outcome_model = Y ~ L1 + A,
#' censor_model = C ~ L1 + A,
#' interventions = list(natural_course()),
#' intervention_names = list("A"),
#' compute_nc_risk = T,
#' hazard_based = F,
#' intervention_description = "Natural Course")
#'
#' ice_natural_course
#'
#' # Estimate risks for the dynamic treat intervention
#' # based on the covariate L1 > 0
#' # (i.e. treat when L1 > 0 and absorbing once one initiates treatment;
#' # not treat otherwise)
#' # Classical pooled ICE, competing event as direct effect
#'
#' ice_dynamic <- ice_pool(data = test_data, K = 5, id = "id",
#' time_name = "t0", outcome_name = "Y",
#' competing_name = "D", censor_name = "C",
#' total_effect = F,
#' outcome_model = Y ~ L1 + A,
#' censor_model = C ~ L1 + A,
#' interventions = list(dynamic("absorbing", "L1", ">", 0)),
#' intervention_names = list("A"),
#' hazard_based = F,
#' intervention_description = "Dynamic Treat")
#'
#' ice_dynamic
#'
#' # Estimate risks for the dynamic intervention where treat when L1 = 0
#' # with uniform grace period of 2 periods
#' # Hazard extended pooled ICE, competing event as total effect
#'
#' ice_grace_period <- ice_pool(data = test_data, K = 5, id = "id",
#' time_name = "t0", outcome_name = "Y",
#' competing_name = "D", censor_name = "C",
#' total_effect = T,
#' outcome_model = Y ~ L1 + A,
#' censor_model = C ~ L1 + A,
#' competing_model = D ~ L1 + A,
#' interventions = list(grace_period("uniform", 2, "L1", 0)),
#' intervention_names = list("A"),
#' hazard_based = T,
#' intervention_description = "Dynamic Treat Grace Period")
#'
#' ice_grace_period
#'
#' # Estimate risks for the threshold intervention where
#' # when the natural value of treatment A at time t is lower
#' # than -3, set its value to -3. Otherwise, do not intervene.
#' # Hazard extended pooled ICE, competing event as total effect
#'
#' ice_threshold <- ice_pool(data = test_data, K = 5, id = "id",
#' time_name = "t0", outcome_name = "Y",
#' competing_name = "D", censor_name = "C",
#' total_effect = T,
#' outcome_model = Y ~ L1 + A,
#' censor_model = C ~ L1 + A,
#' competing_model = D ~ L1 + A,
#' interventions = list(threshold(-3)),
#' intervention_names = list("A"),
#' hazard_based = T,
#' intervention_description = "Threshold Intervention")
#'
#' ice_threshold

# source("helper.R")
# library(reshape2)
# library(ggplot2)
# library(tidyverse)

ice_pool <- function(data, K, id, time_name, outcome_name,
                     censor_name = NULL, competing_name = NULL,
                     total_effect, outcome_model, censor_model = NULL,
                     competing_model = NULL, interventions,
                     intervention_names, intervention_times = NULL,
                     compute_nc_risk = T, hazard_based,
                     intervention_description)
{

  ## 0. some pre-processing

  outcome_varname <- outcome_name
  competing_varname <- competing_name
  censor_varname <- censor_name
  intervention_varnames <- intervention_names
  time0 <- time_name

  if (is.null(intervention_times) | (length(intervention_times[[1]]) == 0)) {
    intervention_times <- list()
    for (i in 1:length(intervention_varnames[[1]])) {
      intervention_times <- c(intervention_times, list(0:(K-1)))
    }
    intervention_times <- list(intervention_times)
  }

  outcome_covar <- str_remove_all(unlist(str_split(as.character(outcome_model)[3], "[+]")), " ")
  censor_covar <- str_remove_all(unlist(str_split(as.character(censor_model)[3], "[+]")), " ")
  competing_covar <- str_remove_all(unlist(str_split(as.character(competing_model)[3], "[+]")), " ")
  
  ## first create lag terms
  
  outcome_covar_lags <- outcome_covar[str_detect(outcome_covar, "lag[0-9999]_")]
  censor_covar_lags <- censor_covar[str_detect(censor_covar, "lag[0-9999]_")]
  competing_covar_lags <- competing_covar[str_detect(competing_covar, "lag[0-9999]_")]
  
  # all_lags <- as.set(unlist(union(outcome_covar_lags, censor_covar_lags)))
  # all_lags <- unlist(union(all_lags, competing_covar_lags))
  
  all_lags <- c(outcome_covar_lags, censor_covar_lags, competing_covar_lags)
  all_lags <- unique(all_lags)
  
  # detach("package:sets", unload = T)
  
  if (length(all_lags) > 0) {
  
  data$factor_id <- as.factor(data[, id])
  
  for (i in 1:length(all_lags)) {
    ilag <- all_lags[i]
    lag_components <- str_split(ilag, "_")
    lag_unit <- as.numeric(str_replace_all(lag_components[[1]][1], "lag", ""))
    lag_var <- lag_components[[1]][2]
    

    group_dta <- data %>% group_by(factor_id)
    group_dta[, "lag_var_new"] <- data[, lag_var]
    group_dta <- group_dta %>%
      mutate(lagged_var = dplyr::lag(lag_var_new, n = lag_unit, default = 0))
    
    data[, ilag] <- group_dta$lagged_var
    
  }
  }
  
  ## need to rebuild covars by calling the transform functions
  ## outcome covar
  outcome_covar_new <- c()
  data_add <- data
  for (i in 1:length(outcome_covar)) {
    
    icovar <- outcome_covar[i]
    
    column_name <- get_column_name_covar(icovar)
    
    if (!column_name %in% colnames(data)) {
      old_ncol <- ncol(data_add)
      data_add <- data %>% mutate(eval(parse(text = icovar)))
      new_ncol <- ncol(data_add)
      
      new_column <- data_add[, (old_ncol+ 1) : new_ncol]
      new_column <- as.matrix(new_column)
      
      ## create new columns
      
      if (ncol(new_column) > 1) {
        multiple_varname <- paste0(paste0(column_name, "."), 1:(ncol(new_column)))
        data[, multiple_varname] <- new_column
        outcome_covar_new <- c(outcome_covar_new, multiple_varname)
      } else {
        data[, column_name] <- as.vector(new_column)
        outcome_covar_new <- c(outcome_covar_new, column_name)
      }
    } else {
      outcome_covar_new <- c(outcome_covar_new, icovar)
    }
  }
  
  ## censor covar
  censor_covar_new <- c()
  data_add <- data
  for (i in 1:length(censor_covar)) {
    
    icovar <- censor_covar[i]
    
    column_name <- get_column_name_covar(icovar)
    
    if (!column_name %in% colnames(data)) {
      old_ncol <- ncol(data_add)
      data_add <- data %>% mutate(eval(parse(text = icovar)))
      new_ncol <- ncol(data_add)
      
      new_column <- data_add[, (old_ncol+ 1) : new_ncol]
      new_column <- as.matrix(new_column)
      
      ## create new columns
      
      if (ncol(new_column) > 1) {
        multiple_varname <- paste0(paste0(column_name, "."), 1:(ncol(new_column)))
        data[, multiple_varname] <- new_column
        censor_covar_new <- c(censor_covar_new, multiple_varname)
      } else {
        data[, column_name] <- as.vector(new_column)
        censor_covar_new <- c(censor_covar_new, column_name)
      }
    } else {
      censor_covar_new <- c(censor_covar_new, icovar)
    }
  }
  
  ## competing covar
  competing_covar_new <- c()
  data_add <- data
  for (i in 1:length(competing_covar)) {
    
    icovar <- competing_covar[i]
    
    column_name <- get_column_name_covar(icovar)
    
    if (!column_name %in% colnames(data)) {
      old_ncol <- ncol(data_add)
      data_add <- data %>% mutate(eval(parse(text = icovar)))
      new_ncol <- ncol(data_add)
      
      new_column <- data_add[, (old_ncol+ 1) : new_ncol]
      new_column <- as.matrix(new_column)
      
      ## create new columns
      
      if (ncol(new_column) > 1) {
        multiple_varname <- paste0(paste0(column_name, "."), 1:(ncol(new_column)))
        data[, multiple_varname] <- new_column
        competing_covar_new <- c(competing_covar_new, multiple_varname)
      } else {
        data[, column_name] <- as.vector(new_column)
        competing_covar_new <- c(competing_covar_new, column_name)
      }
    } else {
      competing_covar_new <- c(competing_covar_new, icovar)
    }
  }
  
  
  ## replace the old covar
  outcome_covar <- outcome_covar_new
  censor_covar <- censor_covar_new
  competing_covar <- competing_covar_new
  


  natural_course_type <- ifelse(is.null(censor_varname), "mean", "ipw")

  fit_all <- fit_summary <- fit_stderr <- fit_vcov <- fit_rmse <- c()

  if (!is.null(competing_varname) & total_effect & hazard_based) {

   data[, outcome_varname] <- ifelse(data[, competing_varname] == 1, 0, data[, outcome_varname])

  }


  ## 1. compute the IPW hazard for natural course

  if (compute_nc_risk) {

    obs_treatment_varname <- unlist(intervention_varnames[[1]])
    
    censor_covar_nc <- censor_covar

    if (!is.null(competing_varname)) {
      
      if (is.null(competing_model)) {
        competing_covar_nc <- outcome_covar
      } else {
        competing_covar_nc <- competing_covar
      }

      if (total_effect == F) {

        competing_formula <- as.formula(paste0(competing_varname, "~",
                                               paste0(competing_covar_nc, collapse = "+")))
        competing_fit <- glm(competing_formula, data = data, family = binomial)

        ## add in this competing fit

        fit_np <- fit_np_summary <- fit_np_stderr <- fit_np_vcov <- fit_np_rmse <- c()

        this_np_fit <- list(competing_fit)
        this_np_summary <- list(get_summary(competing_fit))
        this_np_stderr <- list(get_stderr(competing_fit))
        this_np_vcov <- list(get_vcov(competing_fit))
        this_np_rmse <- list(get_rmse(competing_fit))

        names(this_np_fit) <- names(this_np_summary) <- names(this_np_stderr) <- names(this_np_vcov) <- names(this_np_rmse) <- "NP"

        fit_np <- c(fit_np, this_np_fit)
        fit_np_summary <- c(fit_np_summary, this_np_summary)
        fit_np_stderr <- c(fit_np_stderr, this_np_stderr)
        fit_np_vcov <- c(fit_np_vcov, this_np_vcov)
        fit_np_rmse <- c(fit_np_rmse, this_np_rmse)

        fit_all <- c(fit_all, list(fit_np))
        fit_summary <- c(fit_summary, list(fit_np_summary))
        fit_stderr <- c(fit_summary, list(fit_np_stderr))
        fit_vcov <- c(fit_vcov, list(fit_np_vcov))
        fit_rmse <- c(fit_rmse, list(fit_np_rmse))
      }
    }

    nc_compute <- natural_course_ipweighted(data, id, censor_varname,
                                            K, time0, outcome_varname, censor_covar_nc,
                                            competing_varname, competing_fit,
                                            total_effect)
    risk_weighted <- nc_compute$risk_weighted
    logit_censor <- nc_compute$logit_censor
  } else {
    risk_weighted <- logit_censor <- NULL
  }

  ## 2. compute the intervention for each treatment variable

  if (is.null(censor_varname)) {
    data$C <- 0
    censor_varname <- "C"
  }

  abar_all <- list()
  data <- as.data.frame(data)

  gp_indicator <<- F

  interv_data <<- data

  for (i in 1:length(intervention_varnames[[1]])) {

    treat <<- i
    id_var <<- id

    treatment_varname <<- intervention_varnames[[1]][[treat]]
    intervention_f <- interventions[[1]][[treat]]
    interv_it <- intervention_f
    interv_data[, paste0("interv_it_", treatment_varname, "_", treat)] <- interv_it
    my.arrayofA <- paste0("interv_it_", treatment_varname, "_", treat)


    abar_all[[treat]] <- my.arrayofA
  }

  data <- interv_data


  ## 3. reshape data and fit outcome model
  dffullwide <- reshape(data, idvar = id, timevar = time0, direction = "wide", sep = "_")


  ## might need to change the newly created column names
  tmpdata = as.data.frame(dffullwide)
  formula_full <- as.formula(paste0(outcome_varname,"~", paste0(c(outcome_covar), collapse = "+")))
  yfitog = glm(formula_full, family = binomial(), data = data) #This is from the data generation mechanism

  paramtmp = (yfitog)$coef

  ## add in this outcome fit

  fit_outcome <- fit_outcome_summary <- fit_outcome_stderr <- fit_outcome_vcov <- fit_outcome_rmse <- c()

  this_outcome_fit <- list(yfitog)
  this_outcome_summary <- list(get_summary(yfitog))
  this_outcome_stderr <- list(get_stderr(yfitog))
  this_outcome_vcov <- list(get_vcov(yfitog))
  this_outcome_rmse <- list(get_rmse(yfitog))

  names(this_outcome_fit) <- names(this_outcome_summary) <- names(this_outcome_stderr) <- names(this_outcome_vcov) <- names(this_outcome_rmse) <- "outcome"

  fit_outcome <- c(fit_outcome, this_outcome_fit)
  fit_outcome_summary <- c(fit_outcome_summary, this_outcome_summary)
  fit_outcome_stderr <- c(fit_outcome_stderr, this_outcome_stderr)
  fit_outcome_vcov <- c(fit_outcome_vcov, this_outcome_vcov)
  fit_outcome_rmse <- c(fit_outcome_rmse, this_outcome_rmse)

  if (!is.null(competing_varname) & (total_effect == T) & hazard_based) {
    formula_full_comp <- as.formula(paste0(competing_varname,"~", paste0(c(competing_covar), collapse = "+")))
    yfitog_comp = glm(formula_full_comp, family = binomial(), data = data)
    paramcomp = (yfitog_comp)$coef

    ## add in this competing fit

    fit_comp <- fit_comp_summary <- fit_comp_stderr <- fit_comp_vcov <- fit_comp_rmse <- c()

    this_comp_fit <- list(yfitog_comp)
    this_comp_summary <- list(get_summary(yfitog_comp))
    this_comp_stderr <- list(get_stderr(yfitog_comp))
    this_comp_vcov <- list(get_vcov(yfitog_comp))
    this_comp_rmse <- list(get_rmse(yfitog_comp))

    names(this_comp_fit) <- names(this_comp_summary) <- names(this_comp_stderr) <- names(this_comp_vcov) <- names(this_comp_rmse) <- "competing"

    fit_comp <- c(fit_comp, this_comp_fit)
    fit_comp_summary <- c(fit_comp_summary, this_comp_summary)
    fit_comp_stderr <- c(fit_comp_stderr, this_comp_stderr)
    fit_comp_vcov <- c(fit_comp_vcov, this_comp_vcov)
    fit_comp_rmse <- c(fit_comp_rmse, this_comp_rmse)

    fit_all <- c(fit_all, list(fit_comp))
    fit_summary <- c(fit_summary, list(fit_comp_summary))
    fit_stderr <- c(fit_summary, list(fit_comp_stderr))
    fit_vcov <- c(fit_vcov, list(fit_comp_vcov))
    fit_rmse <- c(fit_rmse, list(fit_comp_rmse))
  }


  ## 4. prepare data for regression at each time step
  for (i in 1:(K - 1)) {

    C_lag <- paste0(censor_varname, "_", i - 1)
    Y_lag <- paste0(outcome_varname, "_", i - 1)

    C <- paste0(censor_varname, "_", i)
    Y <- paste0(outcome_varname, "_", i)

    if (i == 1) {

      tmpdata[, C] <- ifelse(is.na(tmpdata[, Y_lag]) & tmpdata[, C_lag] == 1, 1, tmpdata[, C])
      tmpdata[, Y] <- ifelse(tmpdata[, Y_lag] == 1, 1, tmpdata[, Y])

    } else {

      tmpdata[, C] <- ifelse(is.na(tmpdata[, Y_lag]) & tmpdata[, C_lag] == 1, 1, tmpdata[, C])
      tmpdata[, Y] <- ifelse(!is.na(tmpdata[, Y_lag]) & tmpdata[, Y_lag] == 1, 1, tmpdata[, Y])

    }

  }

  for (i in 1:K) {

    C_lag <- paste0(censor_varname, "_", K - i)
    Y_lag <- paste0(outcome_varname, "_", K - i)

    C <- paste0(censor_varname, "_", K - i + 1)
    Y <- paste0(outcome_varname, "_", K - i + 1)

    tmpdata[, Y] <- tmpdata[, Y_lag]
    tmpdata[, C] <- tmpdata[, C_lag]
  }

  tmpdata$Y_0 = tmpdata$C_0 = NULL

  meany <- matrix(NA, ncol = K + 1, nrow = 1)

  ## 5. regression at each time point
    for (i in 1:K) {

      t <- K - i + 1

      covar_t <- paste0(outcome_covar, sep = paste0("_", t - 1))


        if (t == 1) {
          data_pred_tmp <- data_fit <- tmpdata
        } else {
          if (!is.null(competing_varname) & total_effect == F) {
          data_pred_tmp <- data_fit <- tmpdata[!is.na(tmpdata[, paste0(censor_varname, "_", t - 1)]) & tmpdata[, paste0(censor_varname, "_", t - 1)] == 0 &
                                                 !is.na(tmpdata[, paste0(competing_varname, "_", t - 2)]) & tmpdata[, paste0(competing_varname, "_", t - 2)] == 0, ]
          } else {
            data_pred_tmp <- data_fit <- tmpdata[!is.na(tmpdata[, paste0(censor_varname, "_", t - 1)]) & tmpdata[, paste0(censor_varname, "_", t - 1)] == 0, ]
          }
        }

        for (treat in 1:length(intervention_varnames[[1]])) {
          intervention_variable <- intervention_varnames[[1]][[treat]]
          abar <- abar_all[[treat]][[1]]
          int_time <- intervention_times[[1]][[treat]]

          if ((t-1) %in% int_time) {

        if(is.numeric(abar)){
          data_fit[, paste0(intervention_variable, "_", t - 1)] <- abar
        } else if (is.character(abar)) {

          data_fit[, paste0(intervention_variable, "_", t - 1)] <- data_fit[, paste0(abar, "_", t - 1)]
        }
          data_pred_tmp[, intervention_variable] = data_fit[, paste0(intervention_variable, "_", t - 1)]
          } else {
            data_pred_tmp[, intervention_variable] = data_fit[, paste0(intervention_variable, "_", t - 1)]
          }
        }

          data_pred_tmp[, outcome_covar] = data_fit[, covar_t]

        data_fit[, paste0("y", t, "pred")] = predict(yfitog, newdata = data_pred_tmp, type = "response")

        ### here we implement competing event in different methods (total vs. direct)
        if (!hazard_based) {
          if (t != 1) {
            if (!is.null(competing_varname) & (total_effect == TRUE)) {
              data_fit[, paste0("y", t, "pred")] = case_when(data_fit[, paste0(competing_varname, "_", t - 2)] == 1 ~ 0,
                                                             (data_fit[, paste0(outcome_varname, "_", t - 1)] == 1) ~ 1,
                                                             (data_fit[, paste0(outcome_varname, "_", t - 1)] == 0) & (data_fit[, paste0(competing_varname, "_", t - 2)] == 0) ~ data_fit[, paste0("y", t, "pred")])

            } else {
              data_fit[, paste0("y", t, "pred")] = ifelse(data_fit[, paste0(outcome_varname, "_", t - 1)] == 1, 1, data_fit[, paste0("y", t, "pred")])
            }
          }
        }


        if (t != 1) {

          for (q in 1:(t - 1)) {


            iter <- t - q

            covar_iter <- paste0(outcome_covar, sep = paste0("_", iter -1))

              interact_covar <- NULL

            fit_formula <- as.formula(paste0("y", t, "pred", "~",
                                             paste0(c(covar_iter, interact_covar), collapse = "+")))
            fit_temp = glm(fit_formula, family = quasibinomial(), data = data_fit)

            this_outcome_fit <- list(fit_temp)
            this_outcome_summary <- list(get_summary(fit_temp))
            this_outcome_stderr <- list(get_stderr(fit_temp))
            this_outcome_vcov <- list(get_vcov(fit_temp))
            this_outcome_rmse <- list(get_rmse(fit_temp))

            names(this_outcome_fit) <- names(this_outcome_summary) <- names(this_outcome_stderr) <- names(this_outcome_vcov) <- names(this_outcome_rmse) <- "outcome"

            fit_outcome <- c(fit_outcome, this_outcome_fit)
            fit_outcome_summary <- c(fit_outcome_summary, this_outcome_summary)
            fit_outcome_stderr <- c(fit_outcome_stderr, this_outcome_stderr)
            fit_outcome_vcov <- c(fit_outcome_vcov, this_outcome_vcov)
            fit_outcome_rmse <- c(fit_outcome_rmse, this_outcome_rmse)


            if ((iter - 1) == 0) {
              data_fit = tmpdata
            } else {

              if (!is.null(competing_varname) & total_effect == F) {
              data_fit = tmpdata[!is.na(tmpdata[, paste0(censor_varname, "_", iter - 1)]) & tmpdata[, paste0(censor_varname, "_", iter - 1)] == 0 &
                                   !is.na(tmpdata[, paste0(competing_varname, "_", iter - 2)]) & tmpdata[, paste0(competing_varname, "_", iter - 2)] == 0, ]
              } else {
                data_fit = tmpdata[!is.na(tmpdata[, paste0(censor_varname, "_", iter - 1)]) & tmpdata[, paste0(censor_varname, "_", iter - 1)] == 0, ]
              }
            }

            for (treat in 1:length(intervention_varnames[[1]])) {
              intervention_variable <- intervention_varnames[[1]][[treat]]
              abar <- abar_all[[treat]][[1]]
              int_time <- intervention_times[[1]][[treat]]

              if ((iter - 1) %in% int_time) {

            if(is.numeric(abar)){
              data_fit[, paste0(intervention_variable, "_", iter - 1)] = abar
            } else if (is.character(abar)){

              data_fit[, paste0(intervention_variable, "_", iter - 1)] = data_fit[, paste0(abar, "_", iter - 1)]
            }
              }
            }

            data_fit[, paste0("y", t, "pred")] = predict(fit_temp, newdata = data_fit, type = "response")

            ### total effect and direct effect have different calculation in this step
            if (!hazard_based) {
              if ((iter - 1) != 0) {
                if (!is.null(competing_varname) & (total_effect == TRUE)) {
                  data_fit[, paste0("y", t, "pred")] = case_when(data_fit[, paste0(competing_varname, "_", iter - 2)] == 1 ~ 0,
                                                                 (data_fit[, paste0(outcome_varname, "_", iter - 1)] == 1) & (data_fit[, paste0(competing_varname, "_", iter - 2)] == 0) ~ 1,
                                                                 (data_fit[, paste0(outcome_varname, "_", iter - 1)] == 0) & (data_fit[, paste0(competing_varname, "_", iter - 2)] == 0) ~ data_fit[, paste0("y", t, "pred")])

                } else {
                  data_fit[, paste0("y", t, "pred")] = ifelse(data_fit[, paste0(outcome_varname, "_", iter - 1)] == 1, 1, data_fit[, paste0("y", t, "pred")])
                }
              }
            } else {
                covar_mat <- cbind(rep(1, nrow(data_fit)),
                                   data_fit[, covar_iter])
              predict_tmp <- plogis(as.matrix(covar_mat)  %*% matrix(paramtmp, nrow = length(paramtmp)))

              ## need to calculate d_hat and multiply it with the hazard
              if (!is.null(competing_varname) & total_effect == T) {
                predict_comp <- plogis(as.matrix(covar_mat)  %*% matrix(paramcomp, nrow = length(paramcomp)))
                data_fit[, paste0("y", t, "pred")] <- predict(fit_temp, newdata = data_fit, type = "response") * (1 - predict_tmp) * (1 - predict_comp) + predict_tmp
              } else {
                data_fit[, paste0("y", t, "pred")] <- predict(fit_temp, newdata = data_fit, type = "response") * (1 - predict_tmp) + predict_tmp
              }
            }


          }

        }

        meany[1, c(1, t + 1)] <- c(0, mean(data_fit[, paste0("y", t, "pred")]))

      }


  meany <- as.data.frame(meany)

  time_name <- c()

  for (i in 0:K) {

    time_name <- c(time_name, paste0("time", i))
  }

  colnames(meany) <- time_name
  rownames(meany) <- intervention_description

  fit_all <- c(fit_all, list(fit_outcome))
  fit_summary <- c(fit_summary, list(fit_outcome_summary))
  fit_stderr <- c(fit_summary, list(fit_outcome_stderr))
  fit_vcov <- c(fit_vcov, list(fit_outcome_vcov))
  fit_rmse <- c(fit_rmse, list(fit_outcome_rmse))

  return(list(gformula_risk_last_time = meany[length(meany)], gformula_risk = meany,
              weight_h = risk_weighted, ipw_model = logit_censor,
              fit_models = fit_all, model_summary = fit_summary,
              model_stderr = fit_stderr, model_vcov = fit_vcov,
              model_rmse = fit_rmse))
}


