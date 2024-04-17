#' Static Intervention Function
#'
#' This function implements the static intervention, specifically always treat and never treat.
#' Always treat refers to the constant treatment across all time points. Never treat refers to no treatment across all time points.
#'
#'
#' @param value a numeric specifying the intervention value. 1 for always treat and 0 for never treat.
#' @param data a data frame containing the observed data.
#'
#' @return a vector containing the intervened value in the same size as the number of rows in \code{data}.
#' @export
#'
#' @examples
#' data <- readRDS("test_data_competing.rds")
#' always_treat <- static(value = 1, data = data)
#' always_treat
static <- function(value, data = interv_data) {

  interv_it <- rep(value, nrow(data))

  return(interv_it)
}

#' Natural Course Intervention Function
#'
#' This function implements the natural course intervention under the observed treatment.
#' This function is to be used in the main ICE estimator as an input in the \code{interventions} argument.
#'
#' @param data a data frame containing the observed data.
#' @param treat_var a character string specifying the observed treatment variable in \code{data}.
#'
#' @return a vector containing the intervened value in the same size as the number of rows in \code{data}.
#' @export
#'
#' @examples
#' data <- readRDS("test_data_competing.rds")
#' natural_course <- natural_course(data = data, treat_var = "A")
#' natural_course
natural_course <- function(data = interv_data, treat_var = treatment_varname) {
  interv_it <- data[, treat_var]

  return(interv_it)
}

#' IPW Estimator for Natural Course Risk
#'
#' This is the internal helper function for the ICE estimator for the inverse probability weighted estimate of the natural course risk
#'
#' @param data a data frame containing the observed data.
#' @param id a character string indicating the ID variable name in \code{data}.
#' @param censor_varname a character string indicating the censoring event variable name in \code{data}.
#' @param K a numeric that indicates the total number of time points.
#' @param time0 a character string indicating the time variable name in \code{data}.
#' @param outcome_varname a character string indicating the outcome variable name in \code{data}.
#' @param covar a vector of character strings specifying the covariates used
#' in the probability of censoring model for the inverse probability weighted estimate of the natural course risk.
#' @param competing_varname a character string indicating the competing event variable name in \code{data}.
#' @param competing_fit a fitted model for the competing event model.
#' @param total_effect a logical indicating whether to treat competing event as censoring or total effect.
#'
#' @return A list containing the inverse probability weighted estimates of the natural course risk and
#' a model object for the probability of censoring used in the inverse probability weighted estimate of the natural course risk.
#' @keywords internal

natural_course_ipweighted <- function(data, id, censor_varname,
                                      K, time0, outcome_varname, covar,
                                      competing_varname, competing_fit,
                                      total_effect) {
  risk_weighted <- NULL
  logit_censor <- NULL

    if (!is.null(censor_varname)) {
      censor_formula <- as.formula(paste0(censor_varname, "~",paste0(c(covar), collapse = "+")))
      logit_censor <- glm(censor_formula, data = data, family = binomial)
      prob_censor <- predict(logit_censor, newdata = data, type = "response")
      risk_weighted <- compute_weighted_hazard(prob_censor, data, id, censor_varname,
                                               K, time0, outcome_varname,
                                               competing_varname, competing_fit,
                                               total_effect)
    } else {
      if (!is.null(competing_varname) & total_effect == F) {
        prob_censor <- predict(competing_fit, newdata = data, type = "response")
        risk_weighted <- compute_weighted_hazard(prob_censor, data, id, competing_varname,
                                                 K, time0, outcome_varname,
                                                 NULL, NULL,
                                                 total_effect)
      } else if (!is.null(competing_varname) & total_effect == T) {
        risk_weighted <- compute_weighted_hazard(rep(0, nrow(data)), data, id, censor_varname,
                                                 K, time0, outcome_varname,
                                                 competing_varname, competing_fit,
                                                 total_effect)
      } else {
        risk_weighted <- compute_weighted_hazard(rep(0, nrow(data)), data, id, censor_varname,
                                                 K, time0, outcome_varname,
                                                 competing_varname, competing_fit,
                                                 total_effect)
      }
    }


  return(list(risk_weighted = risk_weighted, logit_censor = logit_censor))
}

#' Dynamic Intervention Function
#'
#' This function implements the dynamic intervention strategy. Dynamic intervention strategy refers to the class of treatment mechanisms that depend on the covariate values. This function provides an example of dynamic intervention strategy. The specific mechanism is described as following:
#' \itemize{
#' \item For \code{type} being \code{compare}, the treatment is determined based on the condition of \code{direction} and \code{value}.
#' For example, if \code{direction} is >, then treat if \code{var} > \code{value}, and not treat otherwise.
#' \item For \code{type} being \code{absorbing}, the initiation of treatment is based on the same mechanism as \code{type} being \code{compare}, and
#' adopts always treat once the treatment is initiated for each individual.
#' }
#' Possible \code{direction} includes \code{>, <, >=, <=, =, !=}.
#'
#' @param type a character string specifying the type of dynamic intervention. Possible inputs are \code{compare} and \code{absorbing}.
#' @param var a character string specifying the intervention variable name in \code{data}.
#' @param direction a character string specifying the comparison of \code{var} and \code{value}.
#' Possible inputs are \code{>, <, >=, <=, =, !=}.
#' @param value a numeric specifying the comparison value on \code{var}.
#' @param id a character string specifying the id variable in \code{data}.
#' @param data a data frame containing the observed data.
#'
#' @return a vector containing the intervened value in the same size as the number of rows in \code{data}.
#' @export
#'
#' @examples
#' data <- readRDS("test_data_competing.rds")
#' dynamic <- dynamic(type = "absorbing", var = "L1", value = 0, direction = ">", id = "id", data = data)
#' dynamic
dynamic <- function(type, var, direction, value, id = id_var, data = interv_data) {

  if (type == "compare" | type == "absorbing") {

    if (direction == ">") {

      interv_it <- ifelse(data[, var] > value, 1, 0)

    } else if (direction == "<") {

      interv_it <- ifelse(data[, var] < value, 1, 0)

    } else if (direction == "<=") {

      interv_it <- ifelse(data[, var] <= value, 1, 0)

    } else if (direction == ">=") {

      interv_it <- ifelse(data[, var] >= value, 1, 0)

    } else if (direction == "=") {

      interv_it <- ifelse(data[, var] == value, 1, 0)

    } else if (direction == "!=") {

      interv_it <- ifelse(data[, var] != value, 1, 0)
    }
    data[, paste0("interv_it")] <- interv_it
    data$interv_it_compare <- interv_it
  }

  if (type == "absorbing") {

    data[, paste0("interv_it")] <- NA

    for (i in 1:length(unique(data[, id]))) {

      interv_index <- which(data[data[, id] == unique(data[, id])[i], ]$interv_it_compare == 1)[1]
      interv_it <- data[data[, id] == unique(data[, id])[i], ]$interv_it_compare
      if (is.na(interv_index)) {
        data[data[, id] == unique(data[, id])[i], paste0("interv_it")] <- interv_it
      } else {
        interv_it[interv_index:length(interv_it)] <- 1
        data[data[, id] == unique(data[, id])[i], paste0("interv_it")] <- interv_it
      }
    }
  }

  interv_it <- data[, paste0("interv_it")]

  return(interv_it)
}

#' Intervention with Grace Period Strategy
#'
#' This function implements a particular intervention with grace period. There are two choices of grace type: "natural" or "uniform," specified by the argument \code{type}.
#' \itemize{
#' \item{Natural grace period describes the following mechanism:
#' When the condition \code{var} = \code{value} is met, initiate treatment in \code{n} time units, where \code{n} is specified by \code{nperiod}.
#' If there is no intervention, follow the observed treatment initiation distribution.}
#' \item{Uniform grace period describes the following mechanism:
#' When the condition \code{var} = \code{value} is met, initiate treatment in \code{n} time units, where \code{n} is specified by \code{nperiod}.
#' If there is no intervention, follow the uniform distribution of treatment initiation.}
#' }
#'
#' @param type a character string specifying the type of grace period strategy. "Uniform" for uniform grace period, and "natural" for natural grace period.
#' @param nperiod a numeric specifying the length of grace period.
#' @param var a character string specifying the variable that the grace period strategy is based on.
#' @param value a numeric specifying the value that \code{var} is compared to.
#' @param data a data frame containing the observed data.
#' @param id a character string indicating the ID variable name in \code{data}.
#' @param time_name a character string indicating the time variable name in \code{data}.
#' @param outcome_name a character string indicating the outcome variable name in \code{data}.
#'
#' @return a vector containing the intervened value in the same size as the number of rows in \code{data}.
#' @export
#'
#' @examples
#' data <- readRDS("test_data_competing.rds")
#' grace_period <- grace_period(type = "uniform", nperiod = 2, var = "L1", value = 1, data = data,
#'                              id = "id", time_name = "t0", outcome_name = "Y")
#' grace_period
grace_period <- function(type, nperiod, var, value, data = interv_data, id = idvar, time_name = time0var, outcome_name = outcomevar) {
  
  my.arrayofA <- 0
  gp_indicator <<- T
  gp_interv_type <<- type
  gp_treatment_var <<- var
  gp_threshold_value <<- value
  ngrace_period <<- nperiod

  grace_period_var <<- gp_treatment_var
  treatment_varname_gp <<- paste0("interv_it_", gp_treatment_var)
  abar <- 0

  data$id_var <- data[, id]
  data$outcome_var <- data[, outcome_name]
  data$time_var <- data[, time_name]

  dffullwide <- reshape(data, idvar = id, timevar = time_name, direction = "wide", sep = "_")
  tmpdata = as.data.frame(dffullwide)
  K <- length(unique(data[, time_name]))

  if (gp_interv_type == "uniform") {

    unique_ids <- unique(tmpdata[, id])

    for (t in 1:K) {

      tmpdata[, paste0("cond_indicator", "_", t - 1)] <- 0
      tmpdata[, paste0("R", "_", t - 1)] <- 0
      tmpdata[, paste0(treatment_varname_gp, "_", t - 1)] <- 0

      for (idx in 1:length(unique_ids)) {
        id_i <- unique_ids[idx]

        if (t == 1) {
          tmpdata[tmpdata$id == id_i, paste0("cond_indicator", "_", t - 1)] <- as.numeric(ifelse(tmpdata[tmpdata$id == id_i, paste0(grace_period_var, "_", t - 1)] != gp_threshold_value, 0, 1))
          tmpdata[tmpdata$id == id_i, paste0("R", "_", t - 1)] <- 0

          tmpdata[tmpdata$id == id_i, paste0(treatment_varname_gp, "_", t - 1)] <- as.numeric(ifelse(tmpdata[tmpdata$id == id_i, paste0(grace_period_var, "_", t - 1)] != gp_threshold_value, 0,
                                                                                                     uniform_sample(tmpdata[tmpdata$id == id_i, paste0("R", "_", t - 1)], ngrace_period)))

        } else if (t > 1){
          if (is.na(tmpdata[tmpdata$id == id_i, paste0(grace_period_var, "_", t - 1)])) {
            tmpdata[tmpdata$id == id_i, paste0(treatment_varname_gp, "_", t - 1)] <- NA
          } else {
            if (tmpdata[tmpdata$id == id_i, paste0(grace_period_var, "_", t - 1)] != gp_threshold_value & tmpdata[tmpdata$id == id_i, paste0("cond_indicator", "_", t - 2)] == 0) {
              tmpdata[tmpdata$id == id_i, paste0("cond_indicator", "_", t - 1)] = 0
            } else {
              tmpdata[tmpdata$id == id_i, paste0("cond_indicator", "_", t - 1)] = 1
            }

            if (tmpdata[tmpdata$id == id_i, paste0(treatment_varname_gp, "_", t - 2)] == 0 & tmpdata[tmpdata$id == id_i, paste0("cond_indicator", "_", t - 2)] == 1) {
              tmpdata[tmpdata$id == id_i, paste0("R", "_", t - 1)] = tmpdata[tmpdata$id == id_i, paste0("R", "_", t - 2)] + 1
            } else {
              tmpdata[tmpdata$id == id_i, paste0("R", "_", t - 1)] = tmpdata[tmpdata$id == id_i, paste0("R", "_", t - 2)]
            }

            tmpdata[tmpdata$id == id_i, paste0(treatment_varname_gp, "_", t - 1)] <- as.numeric(ifelse(tmpdata[tmpdata$id == id_i, paste0("cond_indicator", "_", t - 1)] == 1, 0,
                                                                                                       uniform_sample(tmpdata[tmpdata$id == id_i, paste0("R", "_", t - 1)], ngrace_period)))

            tmpdata[tmpdata$id == id_i, paste0(treatment_varname_gp, "_", t - 1)] <- as.numeric(ifelse(tmpdata[tmpdata$id == id_i, paste0(treatment_varname_gp, "_", t - 2)] == 1, 1,
                                                                                                       tmpdata[tmpdata$id == id_i, paste0(treatment_varname_gp, "_", t - 1)]))
          }

          if (t > ngrace_period + 1) {
            tmpdata[tmpdata$id == id_i, paste0(treatment_varname_gp, "_", t - 1)] <- as.numeric(ifelse(tmpdata[tmpdata$id == id_i, paste0(grace_period_var, "_", t - 1 - ngrace_period)] == gp_threshold_value, 1,
                                                                                                       tmpdata[tmpdata$id == id_i, paste0(treatment_varname_gp, "_", t - 1)]))
          }

        }
      }
    }
  }

  if (gp_interv_type == "natural") {

    for (i in 1:K) {

      t <- K - i + 1

    if (t == 1) {
      tmpdata[, paste0(treatment_varname_gp, "_", t - 1)] <- ifelse(tmpdata[, paste0(grace_period_var, "_", t - 1)] != gp_threshold_value, 0, tmpdata[, paste0(treatment_varname_gp, "_", t - 1)])
    } else if (t > 1 & t <= ngrace_period){

      tmpdata[, paste0(treatment_varname_gp, "_", t - 1)] <- ifelse(tmpdata[, paste0(grace_period_var, "_", t - 1)] != gp_threshold_value, 0, tmpdata[, paste0(treatment_varname_gp, "_", t - 1)])

      tmpdata[, paste0(treatment_varname_gp, "_", t - 1)] <- ifelse(tmpdata[, paste0(treatment_varname_gp, "_", t - 2)] == 1, 1, tmpdata[, paste0(treatment_varname_gp, "_", t - 1)])

    } else if (t > 1 & t > ngrace_period) {

      tmpdata[, paste0(treatment_varname_gp, "_", t - 1)] <- ifelse(tmpdata[, paste0(grace_period_var, "_", t - 1)] != gp_threshold_value, 0, tmpdata[, paste0(treatment_varname_gp, "_", t - 1)])

      tmpdata[, paste0(treatment_varname_gp, "_", t - 1)] <- ifelse(tmpdata[, paste0(grace_period_var, "_", t - 1 - ngrace_period)] == gp_threshold_value, 1, tmpdata[, paste0(treatment_varname_gp, "_", t - 1)])
    }
    }

  }

  tmpdata <- tmpdata[, c(id, paste0(var, "_", 0:(K-1)), paste0(outcome_name, "_", 0:(K-1)))]

  data_long <- reshape(tmpdata, direction = "long",
                       varying = list(paste0(var, "_", 0:(K-1)), paste0(outcome_name, "_", 0:(K-1))),
                       v.names = c("interv_it", "outcome"), sep = "_")

  data_long$time <- data_long$time-1

  append_data <- left_join(data, data_long, by = c("id_var" = "id", "outcome_var" = "outcome", "time_var" = "time"))

  return(append_data$interv_it)
}

#' Uniform Sampling
#'
#' This function is an internal function that helps with the implementation of uniform
#' grace period intervention function.
#'
#' @param r a vector of numerics specifying the time until initiate the treatment.
#' @param duration a numeric specifying the duration of grace period.
#'
#' @return the randomly generated treatment.
#' @keywords internal

uniform_sample <- function(r, duration) {
  p <- 1 / (duration + 1 - r)
  treat <- rbinom(1, 1, p)
  return(treat)

}

#' Calculate Observed Natural Course Risk
#'
#' This functions calculates the inverse probability weighted observed natural course risk.
#'
#' @param prob_censor a vector of numerics specifying the estimated probability of censoring for each individual.
#' @param data a data frame containing the observed data.
#' @param id a character string indicating the ID variable name in \code{data}.
#' @param censor_varname a character string indicating the censor variable name in \code{data}.
#' @param time_points a numeric that indicates the total number of time points.
#' @param time_name a character string indicating the time variable name in \code{data}.
#' @param outcome_name a character string indicating the outcome variable name in \code{data}.
#' @param competing_varname a character string indicating the competing event variable name in \code{data}.
#' @param competing_fit a fitted model for the competing event model.
#' @param total_effect a logical indicating whether to treat competing event as censoring or total effect.
#'
#' @return A list with the first entry as a vector of the mean observed risk.
#' Its second entry is a vector of mean observed survival. Its third entry is a vector of inverse probability weight.
#' 

compute_weighted_hazard <- function(prob_censor, data, id, censor_varname,
                                    time_points, time_name, outcome_name,
                                    competing_varname, competing_fit,
                                    total_effect) {

  censor0_weight <- 1/ (1 - prob_censor)
  censor0_weight_cum <- unlist(tapply(censor0_weight, data[[id]], FUN = cumprod))
  if (!is.null(censor_varname)) {
    ipw_censor <- ifelse(data[[censor_varname]] == 1, 0, censor0_weight_cum)
  } else {
    ipw_censor <- censor0_weight_cum
  }

  if (!is.null(competing_varname) & total_effect == F){
    comprisk_p0_inv <- rep(0, length = nrow(data))
    row_ind <- !is.na(data[[competing_varname]])
    comprisk_p0_inv[row_ind] <- 1 / (1 - stats::predict(competing_fit, type = 'response', newdata = data[row_ind, ]))
    comprisk_inv_cum <- unlist(tapply(comprisk_p0_inv, data[[id]], FUN = cumprod))
    w_d <- ifelse(data[[competing_varname]] == 1 | is.na(data[[competing_varname]]), 0, comprisk_inv_cum)
    ipw_censor <- ipw_censor * w_d
  }

  h_k <- obs_risk_temp <- obs_survival <- rep(NA, times = time_points)
  compevent_risk_temp <- h_k2 <- rep(NA, times = time_points)


  for (i in 0:(time_points - 1)){
    cur_time_ind <- data[[time_name]] == i
    w_cur <- ipw_censor[cur_time_ind]

    if (!is.null(competing_varname) & total_effect == T) {
      w_cur_elimD <- ifelse(data[[competing_varname]][cur_time_ind] == 1, 0, w_cur)
      h_k[i + 1] <- weighted.mean(x = data[[outcome_name]][cur_time_ind], w = w_cur_elimD, na.rm = TRUE)
      h_k2[i + 1] <- weighted.mean(x = data[[competing_varname]][cur_time_ind], w = w_cur, na.rm = TRUE)
      if (i == 0){
        obs_risk_temp[i + 1] <- h_k[i + 1] * (1 - h_k2[i + 1])
        compevent_risk_temp[i + 1] <- h_k2[i + 1]
      } else {
        obs_risk_temp[i + 1] <- h_k[i + 1] * (1 - h_k2[i + 1]) * prod((1 - h_k[1:i]) * (1 - h_k2[1:i]))
        compevent_risk_temp[i + 1] <- h_k2[i + 1] * prod((1 - h_k[1:i]) * (1 - h_k2[1:i]))
      }
    } else {
      h_k[i + 1] <- weighted.mean(x = data[[outcome_name]][cur_time_ind], w = w_cur, na.rm = TRUE)
      if (i == 0){
        obs_risk_temp[i + 1] <- h_k[i + 1]
      } else {
        obs_risk_temp[i + 1] <- h_k[i + 1] * prod(1 - h_k[1:i])
      }
    }


    obs_risk <- cumsum(obs_risk_temp)
    obs_survival <- 1 - obs_risk
    res <- list(risk = obs_risk, surv = obs_survival, w = ipw_censor)

  }

  return(res)
}

#' Calculate product of multiple columns in a data frame
#'
#' This function is an internal function that helps calculate the product of multiple columns in the data frame.
#'
#' @param data a data frame containing the observed data.
#' @param columns a vector of character strings specifying the variables in \code{data} to be used to calculate the product.
#'
#' @return a vector of the product of the variables specified in \code{columns}.
#' @export
#' 
#' @keywords internal
#'
product <- function(data, columns){

  ncol <- length(columns)

  if (ncol == 1) {
    vector <- data[, columns[1]]
  } else {
    vector <- data[, columns[1]]

    for (i in 1:(ncol - 1)) {
      idx <- i + 1
      vector <- vector * data[, columns[idx]]
    }
  }

  return(vector)
}

#' Indicator for the pooling over treatment history ICE estimator in the main function \code{ice}
#'
#' This function identifies the pooling over treatment history ICE estimator. The classical pooling over
#' treatment history ICE estimator is specified by \code{hazard = F}. The hazard based pooling over treatment
#' history ICE estimator is specified by \code{hazard = T}.
#'
#' @param hazard a logical indicating whether to use hazard-based ICE estimator or classical ICE estimator.
#' TRUE for hazard-based estimator. FALSE for classical estimator.
#'
#' @return a logical on whether to use hazard-based ICE estimator.
#' @export
#'
pool <- function(hazard) {
  return(hazard)
}

#' Indicator for the stratifying treatment history ICE estimator in the main function \code{ice}
#'
#' This function identifies the stratifying treatment history ICE estimator. The classical stratifying
#' treatment history ICE estimator is specified by \code{hazard = F}. The hazard based stratifying treatment
#' history ICE estimator is specified by \code{hazard = T}.
#'
#' @param hazard a logical indicating whether to use hazard-based ICE estimator or classical ICE estimator.
#' TRUE for hazard-based estimator. FALSE for classical estimator.
#'
#' @return a logical on whether to use hazard-based ICE estimator.
#' @export
#'
strat <- function(hazard) {
  return(hazard)
}

#' Indicator for the weighted ICE estimator in the main function \code{ice}
#'
#' This function identifies the doubly robust weighted ICE estimator. The treatment models could be specified by
#' \code{treat_model} and the treatment variables could be specified by \code{obs_treatment_varnames}.
#'
#' @param treat_model a list of formulas specifying the model statement for the corresponding treatment variable used in the doubly robust ICE estimator. The length of list must match with
#' the length of \code{obs_treatment_varnames}. Multiple treatments passed in \code{obs_treatment_varnames} are allowed and must follow the corresponding order as in \code{treat_model_covar}.
#' @param obs_treatment_varnames a list of character strings specifying the treatment variables to be used in the model for observed treatments of the weighted ICE estimator.
#'
#' @return the model specification for each treatment model and the treatment variables and
#' the observed treatment variable names extracted from the specified model.
#' @export
#'
weight <- function(treat_model = list()) {
  split_treat <- str_split(as.character(treat_model), " ~ ")
  obs_treat <- lapply(split_treat, function(x) {x[1]})
  obs_treat_unique <- as.list(unique(unlist(obs_treat)))
  return(list(treat_model = treat_model, obs_treat = obs_treat_unique))
}

#' Threshold Intervention
#'
#' This function implements the threshold intervention.
#' At each time, if an individual's treatment value is within the user-specified range inclusively,
#' then follow the natural value of the treatment. Otherwise, if the treatment value is below the
#' lower bound specified by the user, then set to the lower bound of the range,
#' and if the treatment value is above the upper bound specified by the user, then set to the upper bound of the range.
#' For more details, please see Young et al 2014.
#'
#' @param lower_bound a numeric specifying the lower bound of the threshold.
#' @param upper_bound a numeric specifying the upper bound of the threshold.
#' @param var a character string specifying the treatment variable for the intervention.
#' @param data a data frame containing the observed data.
#'
#' @return the intervened values on the intervention variable.
#' @export
#'
#' @references Young JG, Herńan MA, Robins JM. Identification, estimation and approximation of risk under interventions that depend on the natural value of treatment using observational data. Epidemiologic Methods. 2014;3(1):1-19.
#'
#' @examples
#' data <- readRDS("test_data_competing.rds")
#' threshold_treat <- threshold(threshold = -3, var = "A", data = data)
#' threshold_treat
threshold <- function(lower_bound, upper_bound, var = threshold_treatment, data = interv_data){
  interv_it <- case_when(data[, var] >= lower_bound & data[, var] <= upper_bound ~ data[, var],
                         data[, var] > upper_bound ~ upper_bound,
                         data[, var] < lower_bound ~ lower_bound)

  return(interv_it)
}

#' Get Standard Deviation from Model Object
#'
#' @param fit a glm model object.
#'
#' @return the standard errors of the coefficients of the fitted model.
#' @export
#'
#' @keywords internal
get_stderr <- function(fit) {
  return(sqrt(diag(vcov(fit))))
}

#' Get Variance-Covariance Matrix from Model Object
#'
#' @param fit a glm model object.
#'
#' @return the variance-covariance matrices of the parameters of the fitted model.
#' @export
#' 
#' @keywords internal
#'
get_vcov <- function(fit) {
  return(vcov(fit))
}

#' Get Summary Table from Model Object
#'
#' @param fit a glm model object.
#'
#' @return the summary table of the fitted model.
#' @export
#' 
#' @keywords internal
#'
get_summary <- function(fit) {
  return(summary(fit))
}

#' Get Root Mean Square Error from the Model Object
#'
#' @param fit a glm model object.
#'
#' @return the root mean square error (RMSE) values of the fitted model.
#' @export
#' 
#' @keywords internal
#'
get_rmse <- function(fit) {
  return(sqrt(mean(fit$residuals^2)))
}


#' Output the Summary Table from ICE Estimator Object
#'
#' (This function is going to be converted to the S3 method in R so it is named "summary_table" for now.
#' After conversion, users could use the base function summary().)
#'
#' @param ... the ICE estimator objects.
#'
#' @return a data frame containing the summary table for all specified ICE estimator objects.
#' @export
#'
#' @examples
#'
#' fit_classical_pool <- ice(
#' data = data,
#' K = 4,
#' id = "id",
#' time_name = "t0",
#' outcome_name = "Y",
#' censor_name = "C",
#' competing_name = "D",
#' estimator = pool(hazard = F),
#' comp_effect = 0,
#' outcome_model = Y ~ L1 + L2 + A1 + A2,
#' censor_model = C ~ L1 + L2 + A1 + A2,
#' ref_idx = 0,
#' int_descript = c("Static Intervention"),
#' intervention1.A1 = list(static_A1),
#' intervention1.A2 = list(static_A2)
#' )
#'
#' fit_hazard_pool <- ice(
#' data = data,
#' K = 4,
#' id = "id",
#' time_name = "t0",
#' outcome_name = "Y",
#' censor_name = "C",
#' competing_name = "D",
#' estimator = pool(hazard = T),
#' comp_effect = 0,
#' outcome_model = Y ~ L1 + L2 + A1 + A2,
#' censor_model = C ~ L1 + L2 + A1 + A2,
#' ref_idx = 0,
#' int_descript = c("Static Intervention"),
#' intervention1.A1 = list(static_A1),
#' intervention1.A2 = list(static_A2)
#' )
#'
#' summary_table(fit_classical_pool, fit_hazard_pool)
summary_table <- function(...) {
  fit_ice <- list(...)

  summary_all <- data.frame()

  if (length(fit_ice) == 1) {
    return(fit_ice[[1]]$summary)
  } else {

  for (i in 1:length(fit_ice)) {
    ifit <- fit_ice[[i]]
    summary_ice <- ifit$summary
    estimator_name <- ifit$estimator.type
    colnames(summary_ice)[2] <- "ICE Risk"
    summary_ice$Estimator <- estimator_name
    summary_ice$Intervention <- rownames(summary_ice)
    rownames(summary_ice) <- NULL

    summary_ice <- summary_ice[, c("Intervention", "Estimator", head(colnames(summary_ice), -2))]

    summary_all <- rbind(summary_all, summary_ice)
  }

  return(summary_all)
  }
}

#' Create the name of transformed variable name
#'
#' @param icovar the name of covariate to be transformed
#'
#' @return the name of the transformed variable name
#' @internal
get_column_name_covar <- function(icovar) {
  
  if (str_detect(icovar, "I()")) {
    covar_name <- str_replace_all(icovar, "[I()]", "")
    covar_name <- str_replace_all(covar_name, "\\^", "_degree")
  } else if (str_detect(icovar, "poly()")) {
    covar_name <- paste0("poly_", str_split(str_replace_all(icovar, "[poly()]", ""), ",")[[1]][1])
  } else if (str_detect(icovar, "ns()")) {
    covar_name <- paste0("ns_", str_split(str_replace_all(icovar, "[ns()]", ""), ",")[[1]][1])
  } else if (str_detect(icovar, "rcspline.eval()")) {
    covar_name <- paste0("rcspline_", str_split(str_replace_all(icovar, "[rcspline.eval()]", ""), ",")[[1]][1])
  } else {
    covar_name <- icovar
  }
  
  return(covar_name)
}
