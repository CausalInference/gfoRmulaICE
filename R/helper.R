#' Static
#'
#' This function specifies the static intervention, either treat with a constant value or never treat, 
#' on the treatment variable in \code{data}.
#'
#'
#' @param value a number specifying the intervention value. 0 for never treat.
#' @param data a data frame containing the observed data.
#'
#' @return a vector containing the intervened values of the same size as the number of rows in \code{data}.
#' @export
#'
#' @examples
#' data <- gfoRmulaICE::compData
#' always_treat <- static(value = 1, data = data)
static <- function(value, data = interv_data) {

  interv_it <- rep(value, nrow(data))

  return(interv_it)
}

#' Natural Course
#'
#' This function specifies the natural course intervention on the treatment variable in \code{data}.
#'
#' @param data a data frame containing the observed data.
#' @param treat_var a string specifying the treatment variable in \code{data}.
#'
#' @return a vector containing the intervened values of the same size as the number of rows in \code{data}.
#' @export
#'
#' @examples
#' data <- gfoRmulaICE::compData
#' natural_course <- natural_course(data = data, treat_var = "A")
#' 
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
#' @internal

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


#' Strategy with Grace Period
#' 
#' This function specifies an intervention in which treatment is initiated within the grace period of \code{nperiod} time units. 
#' During the grace period, the treatment variable follows its natural value or initiate intervention with a uniform distribution at each time point.
#'
#' @param type a string specifying the type of grace period strategy. Possible values are "uniform" and "natural".
#' @param nperiod a number indicating the length of grace period.
#' @param condition a string specifying the logical expression, upon which is met, the treatment is initiated within \code{nperiod} time units.
#' @param data a data frame containing the observed data.
#' @param id a string specifying the ID variable name in \code{data}.
#' @param time_name a string specifying the time variable name in \code{data}.
#' @param outcome_name a string specifying the outcome variable name in \code{data}.
#'
#' @return a vector containing the intervened value of the same size as the number of rows in \code{data}.
#' @export
#'
#' @examples
#' data <- gfoRmulaICE::compData
#' grace_period <- grace_period(type = "uniform", nperiod = 2, condition = "L1 == 0", data = data,
#'                              id = "id", time_name = "t0", outcome_name = "Y")
grace_period <- function(type, nperiod, condition,
                         # condition = dynamic_cond, 
                         data = interv_data, id = idvar, time_name = time0var, 
                         outcome_name = outcomevar) {
  
  
  ## separate var and logical from condition
  
  ## remove spaces from condition
  
  condition <- str_remove_all(condition, " ")
  
  ## get every character in condition
  
  condition_elements <- unlist(str_split(condition, ""))
  
  ## get the variable
  
  var_idx <- which(str_detect(condition_elements, "[^a-zA-Z0-9]"))[1]
  if (var_idx - 1 <= 0) {
    stop("Please enter a valid condition.")
  } else {
    var <- paste0(condition_elements[1:(var_idx-1)], collapse = "")
    condition_sep <- paste0(condition_elements[(var_idx) : length(condition_elements)], collapse = "")
  }
  
  
  my.arrayofA <- 0
  gp_indicator <<- T
  gp_interv_type <<- type
  gp_treatment_var <<- var
  # gp_threshold_value <<- value
  ngrace_period <<- nperiod

  grace_period_var <<- gp_treatment_var
  treatment_varname_gp <<- paste0("interv_it_", gp_treatment_var)
  abar <- 0

  data$id_var <- data[, id]
  data$outcome_var <- data[, outcome_name]
  data$time_var <- data[, time_name]
  
  data[, "cond_indicator"] <- NA
  data[, "R"] <- NA
  data[, treatment_varname_gp] <- NA

  dffullwide <- reshape(data, idvar = id, timevar = time_name, direction = "wide", sep = "_")
  tmpdata = as.data.frame(dffullwide)
  K <- length(unique(data[, time_name]))

  if (gp_interv_type == "uniform") {

    for (i in 1:K) {
      
      t <- i - 1
      
      int_var_value <- tmpdata[, paste0(grace_period_var, "_", t)]
      opposite_formula <- paste0("!int_var_value", condition_sep)
      
      if (t == 0) {
        tmpdata[,paste0("cond_indicator", "_", t)] <- as.numeric(ifelse(eval(parse(text=opposite_formula)), 0, 1))
        tmpdata[,paste0("R", "_", t)] <- 0
        
        tmpdata[,paste0(treatment_varname_gp, "_", t)] <- as.numeric(ifelse(eval(parse(text=opposite_formula)), 0,
                                                                            uniform_sample(tmpdata[,paste0("R", "_", t)], ngrace_period)))
        
      } else if (t >= 1){
        
            
            tmpdata[, paste0("cond_indicator", "_", t)] <- ifelse(eval(parse(text=opposite_formula)) & tmpdata[, paste0("cond_indicator", "_", t - 1)] == 0, 0, 1)
         
          
          tmpdata[, paste0("R", "_", t)] <- ifelse(tmpdata[, paste0(treatment_varname_gp, "_", t - 1)] == 0 & 
                                                                             tmpdata[, paste0("cond_indicator", "_", t - 1)] == 1, 
                                                                           tmpdata[, paste0("R", "_", t - 1)] + 1, 
                                                                           tmpdata[, paste0("R", "_", t - 1)])
          
          
          
          tmpdata[, paste0(treatment_varname_gp, "_", t)] <- as.numeric(ifelse(tmpdata[, paste0("cond_indicator", "_", t)] == 1,
                                                                               tmpdata[, paste0(treatment_varname_gp, "_", t)],
                                                                               uniform_sample(tmpdata[, paste0("R", "_", t)], ngrace_period)))
          
          tmpdata[, paste0(treatment_varname_gp, "_", t)] <- as.numeric(ifelse(tmpdata[, paste0("cond_indicator", "_", t - 1)] == 1, 
                                                                               1,
                                                                               tmpdata[, paste0(treatment_varname_gp, "_", t)]))
         
          
          if (t > ngrace_period + 1) {
            gp_value <- tmpdata[, paste0(grace_period_var, "_", t - ngrace_period)]
            eval_formula <- paste0("gp_value", condition_sep)
            tmpdata[, paste0(treatment_varname_gp, "_", t)] <- as.numeric(ifelse(eval(parse(text = eval_formula)), 
                                                                                 1,
                                                                                 tmpdata[, paste0(treatment_varname_gp, "_", t)]))
          }
          
        }
      }
    
    }
  

  if (gp_interv_type == "natural") {

    for (i in 1:K) {

      # t <- K - i + 1
      t <- i - 1
      
      int_var_value <- tmpdata[, paste0(grace_period_var, "_", t)]
      opposite_formula <- paste0("!int_var_value", condition_sep)

      tmpdata[, paste0(treatment_varname_gp, "_", t)] <- ifelse(eval(parse(text = opposite_formula)), 0, tmpdata[, paste0(grace_period_var, "_", t)])
        
        if (t >= 1) {
          
          tmpdata[, paste0(treatment_varname_gp, "_", t)] <- 
            ifelse(tmpdata[,paste0(treatment_varname_gp, "_", t-1)] == 1, 
                   1, tmpdata[,paste0(treatment_varname_gp, "_", t)])
        }
        
        if (t >= ngrace_period) {
          
          gp_value <- tmpdata[, paste0(grace_period_var, "_", t - ngrace_period)]
          eval_formula <- paste0("gp_value", condition_sep)
          
          tmpdata[, paste0(treatment_varname_gp, "_", t)] <- ifelse(eval(parse(text = eval_formula)), 
                                                                       1,
                                                                       tmpdata[, paste0(treatment_varname_gp, "_", t)])
          
        }
      }

    
    }

  tmpdata <- tmpdata[, c(id, paste0(treatment_varname_gp, "_", 0:(K-1)), paste0(outcome_name, "_", 0:(K-1)))]
  
  data_long <- reshape(tmpdata, direction = "long",
                       varying = list(paste0(treatment_varname_gp, "_", 0:(K-1)), paste0(outcome_name, "_", 0:(K-1))),
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
#' @internal

uniform_sample <- function(r, duration) {
  p <- 1 / (duration + 1 - r)
  na_p <- p[-which(is.na(r))]
  if (anyNA(na_p)){
    
    warning("NA produced in uniform sampling probabilities in uniform grace period.")
    
  } else if ((any(na_p > 1 | na_p < 0))) {
    
    warning("Uniform sampling probabilities not between 0 and 1 in uniform grace period.")
  } 
  
  treat <- suppressWarnings(rbinom(length(p), 1, p))
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
#' 
#' @internal
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

#' Indicator for the pooling over treatment history ICE estimator
#'
#' This function identifies the pooling over treatment history ICE estimator. The classical pooling over
#' treatment history ICE estimator is specified by \code{pool(hazard = F)}. The hazard based pooling over treatment
#' history ICE estimator is specified by \code{pool(hazard = T)}.
#'
#' @param hazard a logical value indicating whether to use hazard-based ICE estimator.
#'
#' @return a logical value on whether to use hazard-based ICE estimator.
#' @export
#'
pool <- function(hazard) {
  return(hazard)
}

#' Indicator for the stratifying on treatment history ICE estimator
#'
#' This function identifies the stratifying on treatment history ICE estimator. The classical stratifying on 
#' treatment history ICE estimator is specified by \code{strat(hazard = F)}. The hazard based stratifying on treatment
#' history ICE estimator is specified by \code{strat(hazard = T)}.
#'
#' @param hazard a logical value indicating whether to use hazard-based ICE estimator.
#'
#' @return a logical value on whether to use hazard-based ICE estimator.
#' @export
#'
strat <- function(hazard) {
  return(hazard)
}

#' Indicator for the doubly robust ICE estimator
#'
#' This function identifies the doubly robust ICE estimator. The treatment models could be specified by
#' \code{treat_model}.
#'
#' @param treat_model a list of formulas specifying the treatment model for the corresponding treatment variable. The length of list must match
#' the number of treatment variables. 
#'
#' @return treatment model specifications and treatment variable names.
#' @export
#'
weight <- function(treat_model = list()) {
  split_treat <- str_split(as.character(treat_model), " ~ ")
  obs_treat <- lapply(split_treat, function(x) {x[1]})
  obs_treat_unique <- as.list(unique(unlist(obs_treat)))
  return(list(treat_model = treat_model, obs_treat = obs_treat_unique))
}

#' Threshold
#'
#' This function specifies the threshold intervention on the treatment variable in \code{data}.
#' If treatment value is between the lower bound and the upper bound, it follows the natural value of the treatment. 
#' If treatment value is either below the lower bound or above the upper bound, it is set to the lower bound or the upper bound, correspondingly.
#' See Young et al. (2014) for more details.
#'
#' @param lower_bound a number indicating the lower bound of the threshold.
#' @param upper_bound a number indicating the upper bound of the threshold.
#' @param var a string specifying the treatment variable for the intervention.
#' @param data a data frame containing the observed data.
#'
#' @return a vector containing the intervened values of the same size as the number of rows in \code{data}.
#' @export
#'
#' @references Young JG, Herńan MA, Robins JM. Identification, estimation and approximation of risk under interventions that depend on the natural value of treatment using observational data. Epidemiologic Methods. 2014;3(1):1-19.
#'
#' @examples
#' data <- gfoRmulaICE::compData
#' threshold_treat <- threshold(lower_bound = 0, upper_bound = 2, var = "A", data = data)
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
#'
#' @internal
get_stderr <- function(fit) {
  return(sqrt(diag(vcov(fit))))
}

#' Get Variance-Covariance Matrix from Model Object
#'
#' @param fit a glm model object.
#'
#' @return the variance-covariance matrices of the parameters of the fitted model.
#' 
#' @internal
#'
get_vcov <- function(fit) {
  return(vcov(fit))
}

#' Get Summary Table from Model Object
#'
#' @param fit a glm model object.
#'
#' @return the summary table of the fitted model.
#' 
#' @internal
#'
get_summary <- function(fit) {
  return(summary(fit))
}

#' Get Root Mean Square Error from the Model Object
#'
#' @param fit a glm model object.
#'
#' @return the root mean square error (RMSE) values of the fitted model.
#' 
#' @internal
#'
get_rmse <- function(fit) {
  return(sqrt(mean(fit$residuals^2)))
}


#' Summary Table for ICE Estimator Objects
#'
#' This function returns a summary table for ICE estimator objects. (This function is going to be converted to the S3 method in R so it is named "summary_table" for now.
#' After conversion, users could use the base function summary().)
#'
#' @param ... the ICE estimator objects.
#'
#' @return a data frame containing the summary table for all specified ICE estimator objects.
#' @export
#'
#' @examples
#' 
#' data <- gfoRmulaICE::compData
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
#' intervention1.A1 = list(static(3)),
#' intervention1.A2 = list(static(1))
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
#' competing_model = D ~ L1 + L2 + A1 + A2,
#' ref_idx = 0,
#' int_descript = c("Static Intervention"),
#' intervention1.A1 = list(static(3)),
#' intervention1.A2 = list(static(1))
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
    covar_name <- str_replace_all(str_replace_all(icovar, "I", ""), "[()]", "")
    covar_name <- str_replace_all(covar_name, "\\^", "_degree")
  } else if (str_detect(icovar, "poly()")) {
    covar_name <- paste0("poly_", str_split(str_replace_all(str_replace_all(icovar, "poly", ""), "[()]", ""), ",")[[1]][1])
  } else if (str_detect(icovar, "ns()")) {
    covar_name <- paste0("ns_", str_split(str_replace_all(str_replace_all(icovar, "ns", ""), "[()]", ""), ",")[[1]][1])
  } else if (str_detect(icovar, "rcspline.eval()")) {
    covar_name <- paste0("rcspline_", str_split(str_replace_all(str_replace_all(icovar, "rcspline.eval", ""), "[()]", ""), ",")[[1]][1])
  } else {
    covar_name <- icovar
  }
  
  return(covar_name)
}

#' Dynamic
#' 
#' This function specifies a dynamic intervention on the treatment variable specified in \code{data}. 
#' This function follows the treatment strategy specified in \code{strategy_before} until a user-defined condition that depends on covariate values is met. Upon the condition is met, the strategy specified in \code{strategy_after} is followed.
#'
#' @param condition a string that specifies a logical expression, upon which is met, the strategy specified in \code{strategy_after} is followed.
#' @param strategy_before a function or vector of intervened values that specifies the strategy followed after \code{condition} is met. 
#' The vector of intervened values should be the same length as the number of rows in the data frame \code{data}.
#' @param strategy_after a function or vector of intervened values that specifies the strategy followed before \code{condition} is met.
#' The vector of intervened values should be the same length as the number of rows in the data frame \code{data}.
#' @param absorb a logical value indicating whether the strategy specified in \code{strategy_after} becomes absorbing upon the first time when \code{condition} is met.
#' @param id a string specifying the ID variable name in \code{data}.
#' @param time a string specifying the time variable name in \code{data}.
#' @param data a data frame containing the observed data.
#'
#' @return a vector containing the intervened value of the same size as the number of rows in \code{data}.
#' @export
#'
#' @examples
#' data <- gfoRmulaICE::compData
#' # Dynamic intervention example 1: treat when L1 = 0, and not treat otherwise.
#' dynamic1 <- dynamic(condition = "L1 == 0", strategy_before = static(0), strategy_after = static(1), 
#' absorb = FALSE, id = "id", time_name = "t0", data = data)
#' 
#' # Dynamic intervention example 2: never treat upon until L1 = 0, after which follows always treat.
#' dynamic2 <- dynamic(condition = "L1 == 0", strategy_before = static(0), strategy_after = static(1), 
#' absorb = TRUE, id = "id", time_name = "t0", data = data)
#' 
#' # Dynamic intervention example 3: never treat upon until L1 = 0, after which follows natural course.
#' dynamic3 <- dynamic(condition = "L1 == 0", strategy_before = static(0), strategy_after = natural_course(), 
#' absorb = FALSE, id = "id", time_name = "t0", data = data)
dynamic <- function(condition, strategy_before, strategy_after, absorb = FALSE, 
                    id = id_var, time = time0var, data = interv_data) {
  
  first <- absorb
  
  # if (any(str_detect(as.character(substitute(strategy_after)), "grace_period"))) {
  #   
  #   dynamic_cond <<- condition
  #   interv_values <- strategy_after
  #   
  # } else {
    
    strategy_before_values <- strategy_before
    strategy_after_values <- strategy_after
    
    interv_values <- get_dynamic_interv_values(condition, strategy_before_values, 
                                               strategy_after_values, first, id, time, data)
    
  # }
  
  
  return(interv_values)
}

#' Calculate intervened values for dynamic strategy.
#'
#' @param condition a character string that specifies a logical expression, upon which is met, the strategy specified in \code{strategy_after} is followed.
#' @param strategy_before_values a function or vector of intervened values in the same length as the number of rows in the observed data \code{data} that specifies the strategy followed after the condition in \code{condition} is met.
#' @param strategy_after_values a function or vector of intervened values in the same length as the number of rows in the observed data \code{data} that specifies the strategy followed before the condition in \code{condition} is met.
#' @param first a logical indicating whether the strategy specified in \code{strategy_after} starts upon the first time when the condition specified in \code{condition} is met.
#' @param id a character string indicating the ID variable name in \code{data}.
#' @param time a character string indicating the time variable name in \code{data}.
#' @param data a data frame containing the observed data.
#'
#' @return a vector containing the intervened value in the same size as the number of rows in \code{data}.
#' @internal
get_dynamic_interv_values <- function(condition, strategy_before_values, strategy_after_values, first, id, time, data) {
  
  unique_times <- unique(data[, time])
  
  if (first) {
  init_info <- data %>% group_by(id) %>%
    mutate(first_init = unique_times[which(eval(parse(text = condition)))[1]]) %>%
    ungroup()
  
  init_info[, "is_init"] <- init_info[, time] >= init_info[, "first_init"]
  } else {
    init_info <- data %>% 
      mutate(is_init = eval(parse(text = condition)))
  }
  
  init_info[is.na(init_info$is_init), "is_init"] <- FALSE
  init_info[, "strategy_before"] <- strategy_before_values
  init_info[, "strategy_after"] <- strategy_after_values
  
  init_info <- init_info %>% 
    mutate(interv_values = ifelse(is_init, strategy_after_values, strategy_before_values))
  
  return(init_info[, "interv_values"]) 
  
}

#' Replicate lists
#'
#' @param v a vector that needs to be replicated into list(s).
#' @param n a numerical indicating the number of replicates.
#'
#' @return a list containing the replicates of the specified vector.
#' @internal
rep_list <- function(v, n) {
  
  out <- list()
  
  for (i in 1:n) {
    out <- c(out, list(v))
  }
  
  return(out)
}

#' Get model fits from ICE object
#'
#' @param this_fit an ICE object.
#' @param model_name a character string specifying which model to be extracted.
#' @param descript a character string specifying the description of the extracted models.
#'
#' @return a list containing the model information of the user-specified model name.
#' @internal
get_models <- function(this_fit, model_name, descript) {
  this_by_step <- this_fit[[model_name]]
  
  this_model <- this_by_step$fit
  this_summary <- this_by_step$summary
  this_stderr <- this_by_step$stderr
  this_vcov <- this_by_step$vcov
  this_rmse <- this_by_step$rmse
  
  this_fit_all <- list(list(fit = this_model, 
                            summary = this_summary, 
                            stderr = this_stderr, 
                            vcov = this_vcov, 
                            rmse = this_rmse))
  
  names(this_fit_all) <- descript
  
  return(this_fit_all)
}

#' Print out warning messages
#'
#' @param err a list containing whether there is any error from each replicate of bootstrap process.
#' @param err_mssg a list containing the error message of any error from each replicate of bootstrap process.
#' @param model a logical indicating whether the error processed is related to model error or not. Default to be \code{FALSE}.
#'
#' @return None.
#' @internal
give_warning <- function(err, err_mssg, model = F) {
  
  if (model) {
    err_explain <- ". The analysis should likely be repeated
                     with a more parsimonious model."
  } else {
    err_explain <- ". This is likely due to 
                       no data satisfies the defined treatment strategy."
  }
  
  if (length(err) > 0) {
    warning(paste0("NA value occurs in bootstrap replicate ", paste(err, collapse = ","), err_explain, 
                   "\n", paste0(err_mssg, collapse = "")))
  }
}

