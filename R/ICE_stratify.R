#' ICE Stratifying Treatment History Estimator including Both Classical, Hazard-Based, and Doubly Robust Version for A Single
#' Intervention Strategy
#'
#' This function estimates the risk over time for survival outcome using the given observed data set following
#' a single user-defined intervention strategy by the parametric g-formula iterative conditional expectation (ICE)
#' estimator. This function provides the classical ICE stratifying treatment history method,
#' the hazard based ICE stratifying treatment history method, and doubly robust weighted ICE method.
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
#' \item{Dynamic treat:} {\code{dynamic(type, var, value, direction)} }
#' \item{Threshold Intervention:} {\code{threshold(value, var)}}
#' \item{Grace period:} {\code{grace_period(type, nperiod, var, value)} }
#' }
#' For more details, please read the corresponding documentation of each intervention function.
#' If the length of the \code{interventions} list is larger than 1, then the function will process them as multiple treatments simultaneously.
#' The user defined intervention functions must return the intervened value as a vector in the length of the number of rows in \code{data}.
#' The user defined intervention functions in \code{interventions} must be in corresponding order of \code{intervention_names}.
#' The same logic applies for the built in intervention functions.
#' @param intervention_names a list of character strings indicating the intervention variable names in correspondence with the intervention strategies in \code{interventions}.
#' Length of the vector must be in the same size of \code{interventions}.
#' @param compute_nc_risk a logical to indicate whether to compute observed natural course risk. TRUE for observed natural course risk computation. FALSE for no natural course risk computation. Default is TRUE.
#' @param hazard_based a logical indicating whether to use hazard-based ICE estimator or classical ICE estimator. TRUE for hazard-based estimator. FALSE for classical estimator.
#' @param weighted a logical indicating whether to use weighted ICE estimtor. TRUE for doubly robust weighted ICE estimator. FALSE for singly robust ICE estimator. Note, the doubly robust version uses the
#' classical stratified ICE. Default is FALSE.
#' @param treat_model a list of formulae specifying the model statement for the corresponding treatment variable used in the doubly robust ICE estimator. The length of list must match with
#' the length of \code{obs_treatment_names}. Multiple treatments passed in \code{obs_treatment_names} are allowed and must follow the corresponding order as in \code{treat_model_covar}.
#' @param obs_treatment_names a list of character strings specifying the treatment variables to be used in the model for observed treatments of the weighted ICE estimator.
#' @param intervention_description a character string specifying the description of the specified intervention strategy.
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
#' @import tidyverse, data.table, reshape2, nnet
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
#' # (i.e. constantly treat over all time points)
#' # Hazard extended stratified ICE, competing event as direct effect
#'
#' ice_static <- ice_strat(data = test_data, K = 5, id = "id",
#' time_name = "t0", outcome_name = "Y",
#' competing_name = "D", censor_name = "C",
#' total_effect = F,
#' outcome_model = Y ~ L1,
#' censor_model = C ~ L1,
#' competing_model = D ~ L1,
#' interventions = list(static(1)),
#' intervention_names = list("A"),
#' hazard_based = T,
#' obs_treatment_names = list("A"),
#' intervention_description = "Always Treat")
#'
#' ice_static
#'
#' # Estimate risks for both the observed inverse probability weighted
#' # natural course risk and natural course strategy
#' # Classical stratified ICE, competing event as total effect
#'
#' ice_natural_course <- ice_strat(data = test_data, K = 5, id = "id",
#' time_name = "t0", outcome_name = "Y",
#' competing_name = "D", censor_name = "C",
#' total_effect = T,
#' outcome_model = Y ~ L1,
#' censor_model = C ~ L1,
#' interventions = list(natural_course()),
#' intervention_names = list("A"),
#' compute_nc_risk = T,
#' hazard_based = F,
#' obs_treatment_names = list("A"),
#' intervention_description = "Natural Course")
#'
#' ice_natural_course
#'
#' # Estimate risks for the dynamic treat intervention based on the covariate L1 > 0
#' # (i.e. treat when L1 > 0 and absorbing once one initiates treatment;
#' # not treat otherwise)
#' # Classical stratified ICE, competing event as direct effect
#'
#' ice_dynamic <- ice_strat(data = test_data, K = 5, id = "id",
#' time_name = "t0", outcome_name = "Y",
#' competing_name = "D", censor_name = "C",
#' total_effect = F,
#' outcome_model = Y ~ L1,
#' censor_model = C ~ L1,
#' interventions = list(dynamic("absorbing", "L1", ">", 0)),
#' intervention_names = list("A"),
#' hazard_based = F,
#' obs_treatment_names = list("A"),
#' intervention_description = "Dynamic Treat")
#'
#' ice_dynamic
#'
#' # Estimate risks for the always treat intervention
#' # (i.e. constantly treat over all time points)
#' # Weighted ICE, competing event as total effect,
#' # with L1 included in the treatment model for A
#' # treatment model: A ~ L1
#'
#' ice_weighted <- ice_strat(data = test_data, K = 5, id = "id",
#' time_name = "t0", outcome_name = "Y",
#' competing_name = "D", censor_name = "C",
#' total_effect = T,
#' outcome_model = Y ~ L1,
#' censor_model = C ~ L1,
#' competing_model = D ~ L1,
#' interventions = list(static(1)),
#' intervention_names = list("A"),
#' weighted = T, treat_model = list(A ~ L1),
#' obs_treatment_names = list("A"),
#' intervention_description = "Always Treat")
#'
#' ice_weighted
#'

ice_strat <- function(data, K, id, time_name, outcome_name,
                      censor_name = NULL, competing_name = NULL, total_effect,
                      outcome_model, censor_model = NULL, competing_model = NULL,
                      interventions, intervention_names, intervention_times = NULL,
                      compute_nc_risk = T, hazard_based, weighted = F,
                      treat_model, obs_treatment_names,
                      intervention_description) {

  ## 0. some pre-processing

  outcome_varname <- outcome_name
  competing_varname <- competing_name
  censor_varname <- censor_name
  intervention_varnames <- intervention_names
  time0 <- time_name
  obs_treatment_varnames <- obs_treatment_names
  all_treat_vars <- unlist(obs_treatment_varnames)

  fit_all <- fit_summary <- fit_stderr <- fit_vcov <- fit_rmse <- c()

  if (is.null(intervention_times) | (length(intervention_times[[1]]) == 0)) {
    intervention_times <- list()
    for (i in 1:length(intervention_varnames[[1]])) {
      intervention_times <- c(intervention_times, list(0:(K-1)))
    }
    intervention_times <- list(intervention_times)
  }

  if (weighted) {
    hazard_based <- F
    treat_model_covar <- list()
    for (l in 1:length(treat_model)) {
      treat_model_covar <- c(treat_model_covar, list(str_remove_all(unlist(str_split(as.character(treat_model[[l]])[3], "[+]")), " ")))
    }
  }

  outcome_covar <- str_remove_all(unlist(str_split(as.character(outcome_model)[3], "[+]")), " ")
  censor_covar <- str_remove_all(unlist(str_split(as.character(censor_model)[3], "[+]")), " ")
  competing_covar <- str_remove_all(unlist(str_split(as.character(competing_model)[3], "[+]")), " ")


  natural_course_type <- ifelse(is.null(censor_varname), "mean", "ipw")

  if (!is.null(competing_varname) & total_effect & hazard_based) {

    data[, outcome_varname] <- ifelse(data[, competing_varname] == 1, 0, data[, outcome_varname])

  }


  ## 1. compute the IPW hazard for natural course
  if (compute_nc_risk) {

    obs_treatment_varname <- intervention_varnames[[1]]

    if (length(which(censor_covar == obs_treatment_varname)) > 0) {
      censor_covar_nc <- censor_covar[-which(censor_covar == obs_treatment_varname)]
    } else {
      censor_covar_nc <- censor_covar
    }

    if (!is.null(competing_varname)) {
      if (length(which(competing_covar == obs_treatment_varname)) > 0) {
        competing_covar_nc <- competing_covar[-which(competing_covar == obs_treatment_varname)]
      } else if (is.null(competing_model)) {
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
    null_censor <- T
  } else {
    null_censor <- F
  }

  abar_all <- list()
  data <- as.data.frame(data)

  interv_data <<- data

  for (i in 1:length(intervention_varnames[[1]])) {

    treat <<- i
    id_var <<- id

    treatment_varname <<- intervention_varnames[[1]][[treat]]
    intervention_f <- interventions[[1]][[treat]]

    if (any(str_detect(as.character(substitute(interventions[[treat]])), "grace_period"))) {
      my.arrayofA <- paste0("interv_it_", gp_treatment_var)
    } else {
      interv_it <- intervention_f
      interv_data[, paste0("interv_it_", treatment_varname, "_", treat)] <- interv_it
      my.arrayofA <- paste0("interv_it_", treatment_varname, "_", treat)
    }

    abar_all[[treat]] <- my.arrayofA
  }

  data <- interv_data

  ## 3a. if weighted ICE, then

  if (weighted & hazard_based == F) {
    ## weighted ICE
    ntreatment = length(obs_treatment_varnames)
    data[, "pred_obs_prod"] <- 1

    for (i in 1:ntreatment) {
      treatment_i <- obs_treatment_varnames[[i]]
      nA_level <- length(unique(data[, treatment_i]))
      acovar <- treat_model_covar[[i]]

      aformula <- as.formula(paste0(treatment_i, "~", paste0(acovar, collapse = "+")))

      cformula <- as.formula(paste0(censor_varname, "~", paste0(c(censor_covar, treatment_i), collapse = "+")))

      if (!is.null(competing_varname) & total_effect == F) {
        dformula <- as.formula(paste0(competing_varname, "~", paste0(c(competing_covar, treatment_i), collapse = "+")))
      }

      if (null_censor) {
        data$pred_c = 0
      } else {
        cfit = glm(cformula, family = binomial(), data = data)
        data$pred_c = predict(cfit, newdata = data, type="response")
      }

      if (!is.null(competing_varname) & total_effect == F) {
        dfit = glm(dformula, family = binomial(), data = data)
        data$pred_d = predict(dfit, newdata = data, type="response")

        data$pred_d = (1- data$pred_d) * (1-data$pred_c)

      }


      if (nA_level > 2) {
        afit = multinom(aformula, data = data, trace = F)
        predicted_df <- fitted(afit)
        A_names <- sort(unique(data[, treatment_i]))
        for (i in 1:nA_level) {
          A_name <- A_names[i]
          data[, paste0("pred_obs_", i)] <- predicted_df[, i]
          match_rows <- which(data[, treatment_i] == A_name)
          data[match_rows, "pred_obs_all"] <- data[match_rows, paste0("pred_obs_", i)]
          data <- data %>% dplyr::select(-c(paste0("pred_obs_", i)))

        }

        pred_obsa <- 1- data$pred_c

        if (!is.null(competing_varname) & total_effect == F) {
          pred_da = ifelse(data$pred_d != 0, (data$pred_d) * pred_obsa, pred_obsa)
          data[, paste0("pred_obs_", treatment_i)] = pred_da
        } else {
          data[, paste0("pred_obs_", treatment_i)] = pred_obsa
        }

        data <- data %>% dplyr::select(-c("pred_obs_all"))

        } else {
          afit = glm(aformula, family = binomial(), data = data)
          


          pred_obs = predict(afit, newdata = data, type="response")
          pred_obs = ifelse(data[, treatment_i]==1, pred_obs, 1-pred_obs)

          pred_obsa <- 1- data$pred_c

          if (!is.null(competing_varname) & total_effect == F) {
            pred_da = ifelse(data$pred_d != 0, (data$pred_d) * pred_obsa, pred_obsa)
            data[, paste0("pred_obs_", treatment_i)] = pred_da
          } else {
            data[, paste0("pred_obs_", treatment_i)] = pred_obsa
          }

        }


      data[, "pred_obs_prod"] <- data[, "pred_obs_prod"] * data[, paste0("pred_obs_", treatment_i)]

      data <- data %>% dplyr::select(-c(paste0("pred_obs_", treatment_i)))
    }



    dffullwide <- reshape(data, idvar = id, timevar = time0, direction = "wide", sep = "_")


    for (i in 1:(K-1)) {
      if (i == 1) {
        dffullwide[, paste0(outcome_varname, "_", i)] <- ifelse(dffullwide[, paste0(outcome_varname, "_", i-1)] == 1, 1, dffullwide[, paste0(outcome_varname, "_", i)])
      } else {
        dffullwide[, paste0(outcome_varname, "_", i)] <- ifelse(!is.na(dffullwide[, paste0(outcome_varname, "_", i-1)]) & dffullwide[, paste0(outcome_varname, "_", i-1)] == 1,
                                                                    1, dffullwide[, paste0(outcome_varname, "_", i)])
      }

    }

      for (t in 0:(K-1)) {
        dffullwide[, paste0("pred_obs_", t)] <- dffullwide[ ,paste0("pred_obs_", "prod", "_", t)]
      }



    for (t in 0:(K-1)) {
      if (t == 0) {
        dffullwide$pi0 = dffullwide$pred_obs_0
      } else {
        dffullwide[, paste0("pi", t)] <- dffullwide[, paste0("pi", t-1)] * dffullwide[, paste0("pred_obs", "_", t)]
      }

    }
  } else {
    dffullwide <- reshape(data, idvar = id, timevar = time0, direction = "wide", sep = "_")
  }

  ## 4. reshape data and fit outcome model
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

  ## 5. prepare data for regression at each time step
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

  if (hazard_based) {
    #### get outcome model for each time point ######

    outcome_pred_times <- list()
    fit_haz <- fit_haz_summary <- fit_haz_stderr <- fit_haz_vcov <- fit_haz_rmse <- c()

    for (i in 1:(K)) {

      stratify_dta_all <- c()
      no_interv_treat <- c()
      interv_treat <- c()
      for (treat in 1:length(intervention_varnames[[1]])) {
        intervention_variable <- intervention_varnames[[1]][[treat]]
        abar <- abar_all[[treat]][[1]]
        int_time <- intervention_times[[1]][[treat]]

        ## change here the index of competing variable
        if (!is.null(competing_varname) & total_effect == F) {
          stratify_dta <- c(paste0(censor_varname, "_", i), paste0(competing_varname, "_", i-1),
                            sapply(0:(i-1), function(x){paste0(intervention_variable, "_", x)}))
        } else {
          stratify_dta <- c(paste0(censor_varname, "_", i),
                            sapply(0:(i-1), function(x){paste0(intervention_variable, "_", x)}))
        }



        if (is.character(abar)) {
          abar_names_extra <- sapply(0:(i-1), function(x){paste0(abar, "_", x)})

          if (!is.null(competing_varname) & total_effect == F) {
            stratify_dta <- sapply(stratify_dta,
                                   function(x){
                                     idx <- which(stratify_dta == x) - 2
                                     if (grepl(censor_varname, x)) {
                                       paste0("(!is.na(", x,"))&", x, "==", 0)
                                     } else if (grepl(competing_varname, x)) {
                                       paste0("(!is.na(", x,"))&", x, "==", 0)
                                     } else if (grepl(outcome_varname, x)) {
                                       paste0("(!is.na(", x,"))&", x, "==", 0)
                                     } else {paste0(x, "==", abar_names_extra[idx])}})
          } else {
            stratify_dta <- sapply(stratify_dta,
                                   function(x){
                                     idx <- which(stratify_dta == x) - 1
                                     if (grepl(censor_varname, x)) {
                                       paste0("(!is.na(", x,"))&", x, "==", 0)
                                     } else {paste0(x, "==", abar_names_extra[idx])}})
          }


        } else {
          if (!is.null(competing_varname) & total_effect == F) {
            stratify_dta <- sapply(stratify_dta,
                                   function(x){
                                     if (grepl(censor_varname, x)) {
                                       paste0("(!is.na(", x,"))&", x, "==", 0)
                                     } else if (grepl(competing_varname, x)) {
                                       paste0("(!is.na(", x,"))&", x, "==", 0)
                                     } else if (grepl(outcome_varname, x)) {
                                       paste0("(!is.na(", x,"))&", x, "==", 0)
                                     } else {paste0(x, "==", abar)}})
          } else {
            stratify_dta <- sapply(stratify_dta,
                                   function(x){
                                     if (grepl(censor_varname, x)) {
                                       paste0("(!is.na(", x,"))&", x, "==", 0)
                                     } else {paste0(x, "==", abar)}})
          }

        }
          stratify_dta_all <- c(stratify_dta_all, stratify_dta)
      if ((i-1) %in% int_time) {
          interv_treat <- c(interv_treat, intervention_variable)
      } else {
          no_interv_treat <- c(no_interv_treat, intervention_variable)
        }
      }

      stratify_dta_combine <- paste0(stratify_dta_all, collapse = "&")

      pred_data <- tmpdata %>% filter(eval(parse(text = stratify_dta_combine)))

      treat_as_covar <- all_treat_vars[!all_treat_vars %in% interv_treat]

      outcome_covar_t <- c(outcome_covar, treat_as_covar)

      covar_outcome_t <- paste0(outcome_covar_t, "_", i-1)



      outcome_formula <- as.formula(paste0(outcome_varname, "_", i, "~",
                                           paste0(covar_outcome_t, collapse = "+")))

      tmp_fit <- glm(outcome_formula, family = binomial(), data = pred_data)

      outcome_pred_times[[i]] <- tmp_fit

      this_haz_fit <- list(tmp_fit)
      this_haz_summary <- list(get_summary(tmp_fit))
      this_haz_stderr <- list(get_stderr(tmp_fit))
      this_haz_vcov <- list(get_vcov(tmp_fit))
      this_haz_rmse <- list(get_rmse(tmp_fit))

      names(this_haz_fit) <- names(this_haz_summary) <- names(this_haz_stderr) <- names(this_haz_vcov) <- names(this_haz_rmse) <- "hazard"

      fit_haz <- c(fit_haz, this_haz_fit)
      fit_haz_summary <- c(fit_haz_summary, this_haz_summary)
      fit_haz_stderr <- c(fit_haz_stderr, this_haz_stderr)
      fit_haz_vcov <- c(fit_haz_vcov, this_haz_vcov)
      fit_haz_rmse <- c(fit_haz_rmse, this_haz_rmse)

    }

    fit_all <- c(fit_all, list(fit_haz))
    fit_summary <- c(fit_summary, list(fit_haz_summary))
    fit_stderr <- c(fit_summary, list(fit_haz_stderr))
    fit_vcov <- c(fit_vcov, list(fit_haz_vcov))
    fit_rmse <- c(fit_rmse, list(fit_haz_rmse))

    if (!is.null(competing_varname) & total_effect == T) {
      comp_pred_times <- list()
      fit_comp <- fit_comp_summary <- fit_comp_stderr <- fit_comp_vcov <- fit_comp_rmse <- c()

      for (i in 1:K) {

        stratify_dta_all <- c()
        no_interv_treat <- c()
        interv_treat <- c()
        for (treat in 1:length(intervention_varnames[[1]])) {
          intervention_variable <- intervention_varnames[[1]][[treat]]
          abar <- abar_all[[treat]][[1]]
          int_time <- intervention_times[[1]][[treat]]

          if (i == 1) {
            stratify_dta <- c(paste0(censor_varname, "_", i), paste0(competing_varname, "_", i-1),
                              sapply(0:(i-1), function(x){paste0(intervention_variable, "_", x)}))
          } else {
            stratify_dta <- c(paste0(censor_varname, "_", i), paste0(competing_varname, "_", i-2),
                              sapply(0:(i-1), function(x){paste0(intervention_variable, "_", x)}))
          }

          if (is.character(abar)) {
            abar_names_extra <- sapply(0:(i-1), function(x){paste0(abar, "_", x)})

            stratify_dta <- sapply(stratify_dta,
                                   function(x){
                                     idx <- which(stratify_dta == x) - 2
                                     if (grepl(censor_varname, x)) {
                                       paste0("(!is.na(", x,"))&", x, "==", 0)
                                     } else if (grepl(competing_varname, x)) {
                                       paste0("(!is.na(", x,"))&", x, "==", 0)
                                     } else if (grepl(outcome_varname, x)) {
                                       paste0("(!is.na(", x,"))&", x, "==", 0)
                                     } else {paste0(x, "==", abar_names_extra[idx])}})


          } else {
            stratify_dta <- sapply(stratify_dta,
                                   function(x){
                                     if (grepl(censor_varname, x)) {
                                       paste0("(!is.na(", x,"))&", x, "==", 0)
                                     } else if (grepl(competing_varname, x)) {
                                       paste0("(!is.na(", x,"))&", x, "==", 0)
                                     } else if (grepl(outcome_varname, x)) {
                                       paste0("(!is.na(", x,"))&", x, "==", 0)
                                     } else {paste0(x, "==", abar)}})

          }
          stratify_dta_all <- c(stratify_dta_all, stratify_dta)

          if ((i-1) %in% int_time) {
          interv_treat <- c(interv_treat, intervention_variable)
          } else {
            no_interv_treat <- c(no_interv_treat, intervention_variable)
          }
        }


        stratify_dta_combine <- paste0(stratify_dta_all, collapse = "&")
        pred_data <- tmpdata %>% filter(eval(parse(text = stratify_dta_combine)))

        treat_as_covar <- all_treat_vars[!all_treat_vars %in% interv_treat]

        competing_covar_t <- c(outcome_covar, treat_as_covar)

        covar_competing_t <- paste0(competing_covar_t, "_", i-1)

        comp_formula <- as.formula(paste0(competing_varname, "_", i-1, "~",
                                          paste0(covar_competing_t, collapse = "+")))

        tmp_fit <- glm(comp_formula, family = binomial(), data = pred_data)

        comp_pred_times[[i]] <- tmp_fit

        this_comp_fit <- list(tmp_fit)
        this_comp_summary <- list(get_summary(tmp_fit))
        this_comp_stderr <- list(get_stderr(tmp_fit))
        this_comp_vcov <- list(get_vcov(tmp_fit))
        this_comp_rmse <- list(get_rmse(tmp_fit))

        names(this_comp_fit) <- names(this_comp_summary) <- names(this_comp_stderr) <- names(this_comp_vcov) <- names(this_comp_rmse) <- "competing"

        fit_comp <- c(fit_comp, this_comp_fit)
        fit_comp_summary <- c(fit_comp_summary, this_comp_summary)
        fit_comp_stderr <- c(fit_comp_stderr, this_comp_stderr)
        fit_comp_vcov <- c(fit_comp_vcov, this_comp_vcov)
        fit_comp_rmse <- c(fit_comp_rmse, this_comp_rmse)

      }

      fit_all <- c(fit_all, list(fit_comp))
      fit_summary <- c(fit_summary, list(fit_comp_summary))
      fit_stderr <- c(fit_summary, list(fit_comp_stderr))
      fit_vcov <- c(fit_vcov, list(fit_comp_vcov))
      fit_rmse <- c(fit_rmse, list(fit_comp_rmse))
    }
  }


  ## 6. regression at each time point
  meany <- matrix(NA, ncol = K + 1, nrow = length(my.arrayofA))

  for (i in 1:K) {
    t <- K - i + 1
    it_run <- T
    it <- t

    if (hazard_based) {
      this_fit <- outcome_pred_times[[t]]
      if (t == 1) {
        fitdata <- tmpdata
        fitdata[, paste0("y", t, "pred")] <- predict(this_fit, newdata = fitdata, type="response")
        it_run <- F
      }

      it <- t - 1
    }

    if (it_run) {
      for (j in 1:it) {

        q <- it - j

        interv_treat <- c()

        for (treat in 1:length(intervention_varnames[[1]])) {
          intervention_variable <- intervention_varnames[[1]][[treat]]
          int_time <- intervention_times[[1]][[treat]]

          times_extra <- 0:q

          if (q %in% int_time) {
            interv_treat <- c(interv_treat, intervention_variable)
          }
        }


        if (q == it-1) {
          stratify_names_extra_all <- c()
          for (treat in 1:length(intervention_varnames[[1]])) {
            intervention_variable <- intervention_varnames[[1]][[treat]]
            abar <- abar_all[[treat]][[1]]

            times_extra <- 0:q

            if (!is.null(competing_varname) & total_effect == F) {
              stratify_names_extra <- c(paste0(censor_varname, "_", q+1), paste0(competing_varname, "_", q),
                                        sapply(times_extra, function(x){paste0(intervention_variable, "_", x)}))
            } else {
              stratify_names_extra <- c(paste0(censor_varname, "_", q+1),
                                        sapply(times_extra, function(x){paste0(intervention_variable, "_", x)}))
            }

            if (is.character(abar)) {
              abar_names_extra <- sapply(times_extra, function(x){paste0(abar, "_", x)})

              if (!is.null(competing_varname) & total_effect == F) {
                stratify_names_extra <- sapply(stratify_names_extra,
                                               function(x){
                                                 idx <- which(stratify_names_extra == x) - 2
                                                 if (grepl(censor_varname, x)) {
                                                   paste0("(!is.na(", x,"))&", x, "==", 0)
                                                 } else if (grepl(competing_varname, x)) {
                                                   paste0("(!is.na(", x,"))&", x, "==", 0)
                                                 } else if (grepl(outcome_varname, x)) {
                                                   paste0("(!is.na(", x,"))&", x, "==", 0)
                                                 } else {
                                                   paste0(x, "==", abar_names_extra[idx])}})
              } else {

                stratify_names_extra <- sapply(stratify_names_extra,
                                               function(x){
                                                 idx <- which(stratify_names_extra == x) - 1
                                                 if (grepl(censor_varname, x)) {
                                                   paste0("(!is.na(", x,"))&", x, "==", 0)
                                                 } else {
                                                     paste0(x, "==", abar_names_extra[idx])}})
              }

            } else {
              if (!is.null(competing_varname) & total_effect == F) {
                stratify_names_extra <- sapply(stratify_names_extra,
                                               function(x){
                                                 if (grepl(censor_varname, x)) {
                                                   paste0("(!is.na(", x,"))&", x, "==", 0)
                                                 } else if (grepl(competing_varname, x)) {
                                                   paste0("(!is.na(", x,"))&", x, "==", 0)
                                                 } else if (grepl(outcome_varname, x)) {
                                                   paste0("(!is.na(", x,"))&", x, "==", 0)
                                                 } else {
                                                     paste0(x, "==", abar)}})
              } else {

                stratify_names_extra <- sapply(stratify_names_extra,
                                               function(x){
                                                 if (grepl(censor_varname, x)) {
                                                   paste0("(!is.na(", x,"))&", x, "==", 0)
                                                 } else {
                                                     paste0(x, "==", abar)}})
              }
            }

            stratify_names_extra_all <- c(stratify_names_extra_all, stratify_names_extra)
          }
          if (q == t - 1) {
            stratify_names_extra_all <- c(stratify_names_extra_all, paste0("(!is.na(", outcome_varname, "_", t,"))"))
          }

          stratify_names_extra_combine <- paste0(stratify_names_extra_all, collapse = "&")
          fitdata <- tmpdata %>% filter(eval(parse(text = stratify_names_extra_combine)))
          if (hazard_based) {
            fitdata[, paste0("y", t, "pred")] <- predict(this_fit, newdata = fitdata, type="response")
          }

        }
        
        treat_as_covar <- all_treat_vars[!all_treat_vars %in% interv_treat]

        if ((length(treat_as_covar) > 0)) {
          warn_mssg <- paste(paste0(intervention_description, ":"), "Add", paste0(treat_as_covar, collapse = ", "), "to model at time", q, ".")
          warning(warn_mssg)
        }
        outcome_covar_t <- c(outcome_covar, treat_as_covar)
        covar_temp <- paste0(outcome_covar_t, "_", q)


        if (hazard_based) {
          temp_formula <- as.formula(paste0("y", t, "pred", "~",
                                            paste0(covar_temp, collapse = "+")))
        } else {
          if (q == t - 1) {
            temp_formula <- as.formula(paste0(outcome_varname, "_", t, "~",
                                              paste0(covar_temp, collapse = "+")))

          } else {
            temp_formula <- as.formula(paste0("y", t, "pred", "~",
                                              paste0(covar_temp, collapse = "+")))

          }
        }

        if (weighted & hazard_based == F) {
          fitdata$w <- 1/fitdata[, paste0("pi", q)]

          if (q == t - 1) {
            fit <- glm(temp_formula, family = quasibinomial(), weight = fitdata$w, data = fitdata)
          } else {
            fit <- glm(temp_formula, family = quasibinomial(), weight = fitdata$w, data = fitdata)
          }

        } else {

          if (q == t - 1) {
            fit <- glm(temp_formula, family = binomial(), data = fitdata)
          } else {
            fit <- glm(temp_formula, family = quasibinomial(), data = fitdata)
          }

        }

        # add fit into list here

        this_outcome_fit <- list(fit)
        this_outcome_summary <- list(get_summary(fit))
        this_outcome_stderr <- list(get_stderr(fit))
        this_outcome_vcov <- list(get_vcov(fit))
        this_outcome_rmse <- list(get_rmse(fit))

        names(this_outcome_fit) <- names(this_outcome_summary) <- names(this_outcome_stderr) <- names(this_outcome_vcov) <- names(this_outcome_rmse) <- "outcome"

        fit_outcome <- c(fit_outcome, this_outcome_fit)
        fit_outcome_summary <- c(fit_outcome_summary, this_outcome_summary)
        fit_outcome_stderr <- c(fit_outcome_stderr, this_outcome_stderr)
        fit_outcome_vcov <- c(fit_outcome_vcov, this_outcome_vcov)
        fit_outcome_rmse <- c(fit_outcome_rmse, this_outcome_rmse)

        stratify_names_all <- c()
        interv_treat <- c()
        for (treat in 1:length(intervention_varnames[[1]])) {
          intervention_variable <- intervention_varnames[[1]][[treat]]
          abar <- abar_all[[treat]][[1]]
          int_time <- intervention_times[[1]][[treat]]

          times <- 0:(q-1)

          if (q %in% int_time) {
            interv_treat <- c(interv_treat, intervention_variable)
          }

          if (hazard_based) {
            if (!is.null(competing_varname) & total_effect == F) {
              stratify_names <- c(paste0(censor_varname, "_", q+1), paste0(competing_varname, "_", q),
                                  sapply(times, function(x){paste0(intervention_variable, "_", x)}))
            } else {
              stratify_names <- c(paste0(censor_varname, "_", q+1),
                                  sapply(times, function(x){paste0(intervention_variable, "_", x)}))
            }
          } else {
            if (!is.null(competing_varname) & total_effect == F) {
              stratify_names <- c(paste0(censor_varname, "_", q), paste0(competing_varname, "_", q-1),
                                  sapply(times, function(x){paste0(intervention_variable, "_", x)}))
            } else {
              stratify_names <- c(paste0(censor_varname, "_", q),
                                  sapply(times, function(x){paste0(intervention_variable, "_", x)}))
            }
          }

          if (is.character(abar)){
            abar_names <- sapply(times, function(x){paste0(abar, "_", x)})

            if (!is.null(competing_varname) & total_effect == F) {
              stratify_names <- sapply(stratify_names,
                                       function(x){
                                         idx <- which(stratify_names == x) - 2
                                         if (grepl(censor_varname, x)) {
                                           paste0("(!is.na(", x,"))&", x, "==", 0)
                                         } else if (grepl(competing_varname, x)) {
                                           paste0("(!is.na(", x,"))&", x, "==", 0)
                                         } else {paste0(x, "==", abar_names[idx])}})
            } else {
              stratify_names <- sapply(stratify_names,
                                       function(x){
                                         idx <- which(stratify_names == x) - 1
                                         if (grepl(censor_varname, x)) {
                                           paste0("(!is.na(", x,"))&", x, "==", 0)
                                         } else {paste0(x, "==", abar_names[idx])}})
            }
          } else {

            if (!is.null(competing_varname) & total_effect == F) {
              stratify_names <- sapply(stratify_names,
                                       function(x){
                                         if (grepl(censor_varname, x)) {
                                           paste0("(!is.na(", x,"))&", x, "==", 0)
                                         } else if (grepl(competing_varname, x)) {
                                           paste0("(!is.na(", x,"))&", x, "==", 0)
                                         } else {paste0(x, "==", abar)}})
            } else {
              stratify_names <- sapply(stratify_names,
                                       function(x){
                                         if (grepl(censor_varname, x)) {
                                           paste0("(!is.na(", x,"))&", x, "==", 0)
                                         } else {paste0(x, "==", abar)}})
            }

          }
          stratify_names_all <- c(stratify_names_all, stratify_names)
        }
        stratify_names_combine <- paste0(stratify_names_all, collapse = "&")

        if (q == 0) {
          fitdata <- tmpdata
        } else {
          fitdata <- tmpdata %>% filter(eval(parse(text = stratify_names_combine)))
        }

        if (weighted & hazard_based == F) {
          fitdata$w <- 1/fitdata[, paste0("pi", q)]
          fitdata[, paste0("y", t, "pred")] = predict(fit, newdata = fitdata, type="response", weight = w)

        } else {
          fitdata[, paste0("y", t, "pred")] = predict(fit, newdata = fitdata, type="response")
        }

        if (hazard_based) {
          weight_haz <- predict(outcome_pred_times[[q+1]], newdata = fitdata, type="response")

          if (!is.null(competing_varname) & total_effect == T) {
            weight_comp <- predict(comp_pred_times[[q+1]], newdata = fitdata, type="response")
            fitdata[, paste0("y", t, "pred")] <- fitdata[, paste0("y", t, "pred")] * (1 - weight_haz) * (1 - weight_comp) + weight_haz
          } else {
            fitdata[, paste0("y", t, "pred")] = fitdata[, paste0("y", t, "pred")] *
              (1 - weight_haz) + weight_haz
          }
        } else {
          if (q != 0) {
            if (!is.null(competing_varname) & total_effect == T) {
              # if (q - 2 == -1) {
              fitdata[, paste0("y", t, "pred")] = case_when(fitdata[, paste0(competing_varname, "_", q-1)] == 1 ~ 0,
                                                            (fitdata[, paste0(outcome_varname, "_", q)] == 1) & (fitdata[, paste0(competing_varname, "_", q-1)] == 0) ~ 1,
                                                            (fitdata[, paste0(outcome_varname, "_", q)] == 0) & (fitdata[, paste0(competing_varname, "_", q-1)] == 0) ~ fitdata[, paste0("y", t, "pred")])

            } else {
              fitdata[, paste0("y", t, "pred")] = ifelse(fitdata[, paste0(outcome_varname, "_", q)] == 1,
                                                         1, fitdata[, paste0("y", t, "pred")])
            }
          }
        }

      }
    }

    meany[, t + 1] <- mean(fitdata[, paste0("y", t, "pred")])

  }
  meany[, 1] <- 0
  meany <- as.data.frame(meany)

  fit_all <- c(fit_all, list(fit_outcome))
  fit_summary <- c(fit_summary, list(fit_outcome_summary))
  fit_stderr <- c(fit_summary, list(fit_outcome_stderr))
  fit_vcov <- c(fit_vcov, list(fit_outcome_vcov))
  fit_rmse <- c(fit_rmse, list(fit_outcome_rmse))

  time_name <- c()

  for (i in 0:K) {

    time_name <- c(time_name, paste0("time", i))
  }

  colnames(meany) <- time_name
  rownames(meany) <- intervention_description

  return(list(gformula_risk_last_time = meany[length(meany)], gformula_risk = meany,
              weight_h = risk_weighted, ipw_model = logit_censor,
              fit_models = fit_all, model_summary = fit_summary,
              model_stderr = fit_stderr, model_vcov = fit_vcov,
              model_rmse = fit_rmse))
}
