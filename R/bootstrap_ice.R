#' Bootstrap for ICE estimator
#'
#' This function estimates the variance of the point estimates for all user-defined interventions from
#' the ICE estimator object. The variance is calculated by bootstrapping and is used to calculate the confidence
#' interval for each point estimate of each user-defined intervention specified in the ICE estimator object.
#'
#' @param f a function specifying which ICE estimator to use for bootstrap.
#' @param K a numerical value that indicates the total number of time points.
#' @param nboot a numeric indicating number of bootstrap samples.
#' @param sig_level a numeric indicating the significance level to be used for confidence interval.
#' @param parallel a logical indicating whether to parallelize the bootstrap process.
#' @param ncores a numeric indicating the number of CPU cores to be used in parallel computing.
#' @param ref_description a string specifying a description for the reference intervention.
#' @param ref_intervention_varname a string specifying the intervention variable used in the reference intervention.
#' @param ref_total_effect a logical indicating how the competing event is handled for the reference intervention.
#' TRUE for total effect. FALSE for direct effect.
#' @param total_effect a logical indicating how the competing event is handled for the defined intervention.
#' TRUE for total effect. FALSE for direct effect.
#' @param natural_course_type a character string indicating how the natural course strategy is computed for the defined intervention.
#' "ipw" for IPW natural course risk estimation.
#' @param ref_intervention a list of function specifying the intervention to be used as reference.
#' for the inverse probability weighted estimate of the natural course risk in the reference intervention.
#' @param interventions a list of functions defining the interventions to be used in this bootstrap.
#' @param intervention_varnames a list of character strings corresponding the interventions to be used in this bootstrap.
#' @param intervention_description a character string specifying a description for the define intervention in this bootstrap.
#' @param data a data frame containing the observed data in long format.
#' @param id a character string indicating the ID variable name in \code{data}.
#' @param set_seed a numeric indicating the starting seed for bootstrap.
#' @param ... any keyword arguments to be passed in \code{f}.
#'
#' @return A list containing the following components:
#' \item{ice_se}{The standard error for the ICE risk of the defined intervention.}
#' \item{ref_se}{The standard error for the ICE risk of the reference intervention.}
#' \item{rr_se}{The standard error for the ICE risk of the risk ratio.}
#' \item{rd_se}{The standard error for the ICE risk of the risk difference.}
#' \item{ice_cv}{The empirical quantile for the ICE risk of the defined intervention.}
#' \item{ref_cv}{The empirical quantile for the ICE risk of the reference intervention.}
#' \item{rr_cv}{The empirical quantile for the ICE risk of the risk ratio.}
#' \item{rd_cv}{The empirical quantile for the ICE risk of the risk difference.}
#' \item{boot_data}{A list containing all the bootstrapped data samples.}
#' \item{boot_models}{A list containing the fitted models for the outcome, the treatment (if applicable), and the competing event (if applicable) on the bootstrapped samples.}
#' \item{boot_summary}{A list containing the summary of the fitted models on the bootstrapped samples.}
#' \item{boot_stderr}{A list containing the standard errors of the coefficients of the fitted models on the bootstrapped samples.}
#' \item{boot_vcov}{A list containing the variance-covariance matrices of the parameters of the fitted models on the bootstrapped samples.}
#' \item{boot_rmse}{A list containing the root mean square error (RMSE) values of the fitted models on the bootstrapped samples.}
#'
#' @export
#' @import doParallel
#'
bootstrap_ice <- function(f, K, nboot, sig_level, parallel, ncores, ref_description,
                          ref_intervention_varname, ref_total_effect, total_effect,
                          natural_course_type,
                          ref_intervention,
                          interventions, intervention_varnames, intervention_description,
                          data, id, set_seed, ...) {

  set.seed(set_seed)

  n <- length(unique(data[, id]))

  ice_boot <- c()
  ref_boot <- c()
  rr_boot <- c()
  rd_boot <- c()

  ice_boot_all <- matrix(NA, ncol = K + 1, nrow = nboot)
  ref_boot_all <- matrix(NA, ncol = K + 1, nrow = nboot)
  ref_ipw_boot_all <- matrix(NA, ncol = K + 1, nrow = nboot)
  
  if (nboot < 10) {
    nshow <- 1
  } else {
  nshow <- floor(nboot / 10)
  }
  
  if (parallel == T) {
    
    old_dir <- getwd()

    cl <- makeCluster(ncores)
    registerDoParallel(cl)

    result <- matrix(NA, ncol = 4, nrow = nboot)
    result <- as.data.frame(result)
    

    boot_result <- foreach(i = 1:nboot, .combine = "c", .packages = c("tidyverse", "nnet", "stringr", "base"),
                           .errorhandling='pass') %do% {
      
      unique_idx <- unique(data[, id])

      select_index <- sample(1:length(unique_idx), replace = T)
      boot_index <- unique_idx[select_index]
      nindx <- length(boot_index)
      
      row_idx <- lapply(1:nindx, function(i) {
        select_idx <- which(data[, id] %in% boot_index[i])
        return(select_idx)
      }
      )
      
      row_idx <- unlist(row_idx)
      
      boot_data <- lapply(1:nindx, function(i) {
        select_data <- data[which(data[, id] %in% boot_index[i]), ] 
        select_data[, paste0("original_", id)] <- select_data[, id]
        select_data[, id] <- i
        return(select_data)
      }
      )
      
      boot_data <- dplyr::bind_rows(boot_data)
      
      # row_idx <- which(data[, id] %in% boot_index)

      # boot_data <- data[which(data[, id] %in% boot_index), ]
      

      ice <- try(
        f(data = boot_data, total_effect = total_effect, id = id, K = K,
               interventions = interv_boot(interventions, row_idx), intervention_names = intervention_varnames,
               intervention_description = intervention_description, ...),
        silent = T
      )
      

      ref <- f(data = boot_data, total_effect = ref_total_effect, id = id, K = K,
               interventions = interv_boot(ref_intervention, row_idx), intervention_names = list(ref_intervention_varname),
               intervention_description = ref_description, ...)
      
      if (class(ice) == "try-error") {
        
        if (str_detect(ice[1], "Argument mu must be a nonempty numeric vector")) {
          data_issue_err <- i
          data_issue_mssg <- ice[1]
          
          model_issue_err <- NA
          model_issue_mssg <- NA
        } else {
          data_issue_err <- NA
          data_issue_mssg <- NA
          
          model_issue_err <- i
          model_issue_mssg <- ice[1]
        }

        ice_value <- this_model <- this_stderr <- this_vcov <- this_rmse <- this_summary <- NA
        ice_all <- rep(NA, K+1)
        rr <- NA
        rd <- NA
        
        ref_all <- ref$gformula_risk
        ref_ipw_all <- c(0, ref$weight_h$risk)
        ref_value <- ref$gformula_risk_last_time

      } else {
        
        data_issue_err <- model_issue_err <- data_issue_mssg <- model_issue_mssg <- NA
        
        ice_value <- ice$gformula_risk_last_time
        ice_all <- ice$gformula_risk
        
        this_outcome_init <- list(ice$outcome_init)
        this_comp_init <- list(ice$comp_init)
        this_np_model <- list(ice$np_model)
        
        names(this_outcome_init) <- names(this_np_model) <- names(this_comp_init) <- paste0("Bootstrap Replicate ", i)
        
        this_outcome_by_step <- get_models(ice, "outcome_by_step", paste0("Bootstrap Replicate ", i))
        this_hazard_by_step <- get_models(ice, "hazard_by_step", paste0("Bootstrap Replicate ", i))
        this_comp_by_step <- get_models(ice, "comp_by_step", paste0("Bootstrap Replicate ", i))

        # this_model <- list(ice$fit_models)
        # this_summary <- list(ice$model_summary)
        # this_stderr <- list(ice$model_stderr)
        # this_vcov <- list(ice$model_vcov)
        # this_rmse <- list(ice$model_rmse)
        
        ref_all <- ref$gformula_risk
        ref_ipw_all <- c(0, ref$weight_h$risk)
        ref_value <- ref$gformula_risk_last_time
        rr <- ice_value / ref_value
        rd <- ice_value - ref_value

      }

      numerical_out <- rbind(c(ice_value, ref_value, rr, rd, ice_all, ref_all, ref_ipw_all))

      output <- list(boot_out = numerical_out, data = boot_data, data_issue = data_issue_err, 
                     model_issue = model_issue_err,
                     data_issue_mssg = data_issue_mssg, model_issue_mssg = model_issue_mssg, 
                     outcome_init_boot = this_outcome_init, comp_init_boot = this_comp_init, 
                     np_model_boot = this_np_model, 
                     outcome_by_step_boot = this_outcome_by_step, hazard_by_step_boot = this_hazard_by_step,
                     comp_by_step_boot = this_comp_by_step)
      
      }

    
    data_err <- na.omit(unlist(boot_result[names(boot_result) == "data_issue"]))
    model_err <- na.omit(unlist(boot_result[names(boot_result) == "model_issue"]))
    
    data_err_mssg <- unlist(boot_result[names(boot_result) == "data_issue_mssg"])
    model_err_mssg <- unlist(boot_result[names(boot_result) == "model_issue_mssg"])
    
    ## construct error message
    data_err_mssg <- na.omit(data_err_mssg)
    model_err_mssg <- na.omit(model_err_mssg)

    data_err_mssg_combine <- paste0("The error message in bootstrap sample ", data_err)
    data_err_mssg_combine <- paste0(data_err_mssg_combine, ": ", data_err_mssg, "\n")

    model_err_mssg_combine <- paste0("The error message in bootstrap sample ", model_err)
    model_err_mssg_combine <- paste0(model_err_mssg_combine, ": ", model_err_mssg, "\n")
    
    if (length(data_err) > 0) {
      warning(paste0("NA value occurs in bootstrap replicate ", paste(data_err, collapse = ","), ". This is likely due to 
                       no data satisfies the defined treatment strategy.", "\n", paste0(data_err_mssg_combine, collapse = "")))
    }
    
    if (length(model_err) > 0) {
      warning(paste0("NA value occurs in bootstrap replicate ", paste(model_err, collapse = ","), ". The analysis should likely be repeated
                     with a more parsimonious model.", "\n", paste0(model_err_mssg_combine, collapse = "")))
    }
    
    outcome_init <- get_boot_models("outcome_init_boot", boot_result)
    comp_init <- get_boot_models("comp_init_boot", boot_result)
    np_model <- get_boot_models("np_model_boot", boot_result)
    outcome_by_step <- get_boot_models("outcome_by_step_boot", boot_result)
    comp_by_step <- get_boot_models("comp_by_step_boot", boot_result)
    hazard_by_step <- get_boot_models("hazard_by_step_boot", boot_result)
    # outcome_init <- boot_result[names(boot_result) == "outcome_init_boot"]
    
    # comp_init <- boot_result[names(boot_result) == "comp_init_boot"]
    # np_model <- boot_result[names(boot_result) == "np_model_boot"]
    # outcome_by_step <- boot_result[names(boot_result) == "outcome_by_step_boot"]
    # comp_by_step <- boot_result[names(boot_result) == "comp_by_step_boot"]
    # hazard_by_step <- boot_result[names(boot_result) == "hazard_by_step_boot"]
    
    result <- matrix(unlist(boot_result[names(boot_result) == "boot_out"]), byrow = T, nrow = nboot)

    agg_result <- apply(result, 2, sd, na.rm = T)
    ice_se <- agg_result[1]
    ref_se <- agg_result[2]
    rr_se <- agg_result[3]
    rd_se <- agg_result[4]
    ice_se_all <- agg_result[5:(5+K)]
    ref_se_all <- agg_result[(6+K):(6+2*K)]
    ref_ipw_se_all <- agg_result[(7+2*K):(length(agg_result))]


    quantile_result_upper <- apply(result, 2, quantile, probs = 1 - sig_level/2, na.rm = T)
    ice_cv_upper <- quantile_result_upper[1]
    ref_cv_upper <- quantile_result_upper[2]
    rr_cv_upper <- quantile_result_upper[3]
    rd_cv_upper <- quantile_result_upper[4]
    ice_cv_all_upper <- quantile_result_upper[5:(5+K)]
    ref_cv_all_upper <- quantile_result_upper[(6+K):(6+2*K)]
    ref_ipw_cv_all_upper <- quantile_result_upper[(7+2*K):(length(quantile_result_upper))]
    
    quantile_result_lower <- apply(result, 2, quantile, probs = sig_level/2, na.rm = T)
    ice_cv_lower <- quantile_result_lower[1]
    ref_cv_lower <- quantile_result_lower[2]
    rr_cv_lower <- quantile_result_lower[3]
    rd_cv_lower <- quantile_result_lower[4]
    ice_cv_all_lower <- quantile_result_lower[5:(5+K)]
    ref_cv_all_lower <- quantile_result_lower[(6+K):(6+2*K)]
    ref_ipw_cv_all_lower <- quantile_result_lower[(7+2*K):(length(quantile_result_lower))]

    boot_data_all <- boot_result[names(boot_result) == "data"]
    # fit_all <- boot_result[names(boot_result) == "fit_all"]
    # fit_summary <- boot_result[names(boot_result) == "fit_summary"]
    # fit_stderr <- boot_result[names(boot_result) == "fit_stderr"]
    # fit_vcov <- boot_result[names(boot_result) == "fit_vcov"]
    # fit_rmse <- boot_result[names(boot_result) == "fit_rmse"]
    
    stopCluster(cl)

  } else {
    boot_data_all <- c()
    outcome_init <- comp_init <- np_model <- outcome_by_step <- comp_by_step <- hazard_by_step <- c()
    # fit_all <- fit_summary <- fit_stderr <- fit_vcov <- fit_rmse <- c()
    data_err <- model_err <- c()
    data_err_mssg <- model_err_mssg <- c()
    for (i in 1:nboot) {

      if (i %% nshow == 0) {
        cat(paste0(intervention_description, " Bootstrap Progress: ", i, "/", nboot, "\n"))
      }

      unique_idx <- unique(data[, id])

      select_index <- sample(1:length(unique_idx), replace = T)
      boot_index <- unique_idx[select_index]
      nindx <- length(boot_index)
      
      row_idx <- lapply(1:nindx, function(i) {
        select_idx <- which(data[, id] %in% boot_index[i])
        return(select_idx)
      }
      )
      
      row_idx <- unlist(row_idx)
      
      boot_data <- lapply(1:nindx, function(i) {
        select_data <- data[which(data[, id] %in% boot_index[i]), ] 
        select_data[, paste0("original_", id)] <- select_data[, id]
        select_data[, id] <- i
        return(select_data)
      }
      )
      
      boot_data <- dplyr::bind_rows(boot_data)

      # row_idx <- which(data[, id] %in% boot_index)
      # 
      # boot_data <- data[row_idx, ]

      boot_data_all <- c(boot_data_all, list(boot_data))

      ice <- try(
        f(data = boot_data, total_effect = total_effect, id = id, K = K,
               interventions = interv_boot(interventions, row_idx), intervention_names = intervention_varnames,
               intervention_description = intervention_description, ...),
        silent = T
      )

      ref <- f(data = boot_data, total_effect = ref_total_effect, id = id, K = K,
               interventions = interv_boot(ref_intervention, row_idx), intervention_names = list(ref_intervention_varname),
               intervention_description = ref_description, ...)
      
      if (class(ice) == "try-error") {
        
        if (str_detect(ice[1], "Argument mu must be a nonempty numeric vector")) {
          data_issue_err <- i
          data_issue_mssg <- ice[1]
          
          model_issue_err <- NA
          model_issue_mssg <- NA
          
        } else {
          
          data_issue_err <- NA
          data_issue_mssg <- NA
          
          model_issue_err <- i
          model_issue_mssg <- ice[1]
        }
        
        this_outcome_init <- this_comp_init <- this_np_model <- this_outcome_by_step <- 
          this_hazard_by_step <- this_comp_by_step <- NA 
        # ice_value <- this_model <- this_stderr <- this_vcov <- this_rmse <- this_summary <- NA
        ice_all <- rep(NA, K+1)
        rr <- NA
        rd <- NA
        
        ref_all <- ref$gformula_risk
        ref_ipw_all <- c(0, ref$weight_h$risk)
        ref_value <- ref$gformula_risk_last_time
        
      } else {
        
        data_issue_err <- data_issue_mssg <- model_issue_err <- model_issue_mssg <- NA
        
        ice_value <- ice$gformula_risk_last_time
        ice_all <- ice$gformula_risk
        
        this_outcome_init <- list(ice$outcome_init)
        this_comp_init <- list(ice$comp_init)
        this_np_model <- list(ice$np_model)
        
        names(this_outcome_init) <- names(this_np_model) <- names(this_comp_init) <- paste0("Bootstrap Replicate ", i)
        
        this_outcome_by_step <- get_models(ice, "outcome_by_step", paste0("Bootstrap Replicate ", i))
        this_hazard_by_step <- get_models(ice, "hazard_by_step", paste0("Bootstrap Replicate ", i))
        this_comp_by_step <- get_models(ice, "comp_by_step", paste0("Bootstrap Replicate ", i))
        
        # this_model <- list(ice$fit_models)
        # this_summary <- list(ice$model_summary)
        # this_stderr <- list(ice$model_stderr)
        # this_vcov <- list(ice$model_vcov)
        # this_rmse <- list(ice$model_rmse)
        
        ref_all <- ref$gformula_risk
        ref_ipw_all <- c(0, ref$weight_h$risk)
        ref_value <- ref$gformula_risk_last_time
        rr <- ice_value / ref_value
        rd <- ice_value - ref_value
        
      }

      # fit_all <- c(fit_all, this_model)
      # fit_summary <- c(fit_summary, this_summary)
      # fit_stderr <- c(fit_stderr, this_stderr)
      # fit_vcov <- c(fit_vcov, this_vcov)
      # fit_rmse <- c(fit_rmse, this_rmse)
      
      outcome_by_step <- c(outcome_by_step, this_outcome_by_step)
      hazard_by_step <- c(hazard_by_step, this_hazard_by_step)
      comp_by_step <- c(comp_by_step, this_comp_by_step)
      
      outcome_init <- c(outcome_init, this_outcome_init)
      comp_init <- c(comp_init, this_comp_init)
      np_model <- c(np_model, this_np_model)

      ice_boot <- c(ice_boot, ice_value)
      ref_boot <- c(ref_boot, ref_value)
      rr_boot <- c(rr_boot, rr)
      rd_boot <- c(rd_boot, rd)
      ice_boot_all[i, ] <- as.vector(as.matrix(ice_all)[1, ])
      ref_boot_all[i, ] <- as.vector(as.matrix(ref_all)[1, ])
      ref_ipw_boot_all[i, ] <- ref_ipw_all
      
      data_err <- c(data_err, data_issue_err)
      model_err <- c(model_err, model_issue_err)

      data_err_mssg <- c(data_issue_mssg, data_err_mssg)
      model_err_mssg <- c(model_issue_mssg, model_err_mssg)

    }
    
    ice_se <- sd(unlist(ice_boot), na.rm = T)
    ice_cv_upper <- quantile(unlist(ice_boot), probs = 1 - sig_level/2, na.rm = T)
    ice_cv_lower <- quantile(unlist(ice_boot), probs = sig_level/2, na.rm = T)
    ref_se <- sd(unlist(ref_boot), na.rm = T)
    ref_cv_upper <- quantile(unlist(ref_boot), probs = 1 - sig_level/2, na.rm = T)
    ref_cv_lower <- quantile(unlist(ref_boot), probs = sig_level/2, na.rm = T)
    rr_se <- sd(unlist(rr_boot), na.rm = T)
    rr_cv_upper <- quantile(unlist(rr_boot), probs = 1 - sig_level/2, na.rm = T)
    rr_cv_lower <- quantile(unlist(rr_boot), probs = sig_level/2, na.rm = T)
    rd_se <- sd(unlist(rd_boot), na.rm = T)
    rd_cv_upper <- quantile(unlist(rd_boot), probs = 1 - sig_level/2, na.rm = T)
    rd_cv_lower <- quantile(unlist(rd_boot), probs = sig_level/2, na.rm = T)
    ice_se_all <- apply(ice_boot_all, 2, sd, na.rm = T)
    ref_se_all <- apply(ref_boot_all, 2, sd, na.rm = T)
    ref_ipw_se_all <- apply(ref_ipw_boot_all, 2, sd, na.rm = T)
    
    ice_cv_all_upper <- apply(ice_boot_all, 2, quantile, probs = 1 - sig_level/2, na.rm = T)
    ref_cv_all_upper <- apply(ref_boot_all, 2, quantile, probs = 1 - sig_level/2, na.rm = T)
    ref_ipw_cv_all_upper <- apply(ref_ipw_boot_all, 2, quantile, probs = 1 - sig_level/2, na.rm = T)
    
    ice_cv_all_lower <- apply(ice_boot_all, 2, quantile, probs = sig_level/2, na.rm = T)
    ref_cv_all_lower <- apply(ref_boot_all, 2, quantile, probs = sig_level/2, na.rm = T)
    ref_ipw_cv_all_lower <- apply(ref_ipw_boot_all, 2, quantile, probs = sig_level/2, na.rm = T)
    
    data_err <- na.omit(data_err)
    model_err <- na.omit(model_err)
    
    data_err_mssg <- na.omit(data_err_mssg)
    model_err_mssg <- na.omit(model_err_mssg)
    
    ## construct error message
    data_err_mssg_combine <- paste0("The error message in bootstrap sample ", data_err)
    data_err_mssg_combine <- paste0(data_err_mssg_combine, ": ", data_err_mssg, "\n")

    model_err_mssg_combine <- paste0("The error message in bootstrap sample ", model_err)
    model_err_mssg_combine <- paste0(model_err_mssg_combine, ": ", model_err_mssg, "\n")

    
    if (length(data_err) > 0) {
      warning(paste0("NA value occurs in bootstrap replicate ", paste(data_err, collapse = ","), ". This is likely due to 
                       no data satisfies the defined treatment strategy.", "\n", paste0(data_err_mssg_combine, collapse = "")))
    }
    
    if (length(model_err) > 0) {
      warning(paste0("NA value occurs in bootstrap replicate ", paste(model_err, collapse = ","), ". The analysis should likely be repeated
                     with a more parsimonious model.", "\n", paste0(model_err_mssg_combine, collapse = "")))
    }
    

  }

  return(list(ice_se = ice_se_all, ref_se = ref_se_all, rr_se = rr_se, rd_se = rd_se,
              ice_cv_upper = ice_cv_upper, ref_cv_upper = ref_cv_upper, rr_cv_upper = rr_cv_upper, rd_cv_upper = rd_cv_upper,
              ice_cv_all_upper = ice_cv_all_upper, ref_cv_all_upper = ref_cv_all_upper,
              ice_cv_lower = ice_cv_lower, ref_cv_lower = ref_cv_lower, rr_cv_lower = rr_cv_lower, rd_cv_lower = rd_cv_lower,
              ice_cv_all_lower = ice_cv_all_lower, ref_cv_all_lower = ref_cv_all_lower,
              ref_ipw_se = ref_ipw_se_all, ref_ipw_cv_all_upper = ref_ipw_cv_all_upper, ref_ipw_cv_all_lower = ref_ipw_cv_all_lower,
              boot_data = boot_data_all, outcome_init = outcome_init, 
              comp_init = comp_init, np_model = np_model, outcome_by_step = outcome_by_step, 
              comp_by_step = comp_by_step, hazard_by_step = hazard_by_step))
}

#' Get the intervened values on bootstrap sample
#'
#' @param interv a list containing the intervened values corresponding to each intervention.
#' @param idx a vector of indices to be selected for bootstrap sample.
#'
#' @return a list containing the intervened values for bootstrap sample.
#' @export
#'
#' @keywords internal
interv_boot <- function(interv, idx) {
  interv_it <- c()
  for (itreat in 1:length(interv[[1]])) {
    interv_treat <- as.vector(unlist(interv[[1]][[itreat]]))[idx]
    # print(interv_treat)
    interv_it <- c(interv_it, list(interv_treat))
    # print(interv_it)
  }

  return(list(interv_it))
}

#' Extract bootstrap model information from each replicate
#'
#' @param name a character string specifying the name of the model to be extracted.
#' @param boot_result a list containing all the bootstrap results.
#'
#' @return a list containing the model information for each bootstrap replicate.
#' @internal
get_boot_models <- function(name, boot_result) {
  
  all_res <- boot_result[names(boot_result) == name]
  l <- length(all_res)
  
  out <- list()
  
  for (i in 1:l) {
    
    out <- c(out, all_res[i][[name]])
  }
  
  return(out)
  
}
