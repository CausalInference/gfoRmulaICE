#' Bootstrap for ICE estimator
#'
#' This function estimates the variance (and 95\% confidence intervals) of the point estimates via nonparametric bootstrapping.
#' 
#' @param f a function specifying which ICE estimator to use for bootstrap.
#' @param K a number indicating the total number of time points.
#' @param nboot a number indicating the number of bootstrap samples.
#' @param coverage a number greater than 0 and less than 100 indicating the coverage of the confidence interval. Default is 95.
#' @param parallel a logical value indicating whether to parallelize the bootstrap process.
#' @param ncores a number indicating the number of CPU cores to be used in parallel computing.
#' @param ref_description a string describing the reference intervention in this bootstrap.
#' @param ref_intervention_varnames a list of strings specifying the treatment variables to be used for the reference intervention in this bootstrap.
#' @param total_effect a logical value indicating how the competing event is handled for the defined intervention.
#' \code{TRUE} for total effect. \code{FALSE} for controlled direct effect.
#' @param ref_intervention a list of functions specifying the intervention to be used as reference.
#' @param interventions a list of functions defining the intervention to be used in this bootstrap.
#' @param intervention_varnames a list of strings specifying the treatment variables to be used for the defined intervention in this bootstrap.
#' @param intervention_description a string describing the defined intervention in this bootstrap.
#' @param intervention_times a list of numbers indicating the time points to which the defined intervention is applied.
#' @param ref_intervention_times a list of numbers indicating the time points to which the reference intervention is applied.
#' @param data a data frame containing the observed data in long format.
#' @param id a string indicating the ID variable name in \code{data}.
#' @param set_seed a number indicating the starting seed for bootstrap.
#' @param ... any keyword arguments to be passed in \code{f}.
#' 
#' @return A list containing the following components:
#' \item{ice_se}{Standard error for the risk of the defined intervention.}
#' \item{ref_se}{Standard error for the risk of the reference intervention.}
#' \item{rr_se}{Standard error for the risk ratio.}
#' \item{rd_se}{Standard error for the risk difference.}
#' \item{rr_cv_upper}{The (100 - (100 - \code{coverage}) / 2)th percentile for the risk ratio at the last time point. Default is the 97.5th percentile when \code{coverage} is set to its default value 95.}
#' \item{rd_cv_upper}{The (100 - (100 - \code{coverage}) / 2)th percentile for the risk difference at the last time point. Default is the 97.5th percentile when \code{coverage} is set to its default value 95.}
#' \item{ice_cv_all_upper}{The (100 - (100 - \code{coverage}) / 2)th percentile for the risk of the defined intervention at all time points. Default is the 97.5th percentile when \code{coverage} is set to its default value 95.}
#' \item{ref_cv_all_upper}{The (100 - (100 - \code{coverage}) / 2)th percentile for the risk of the reference intervention at all time points. Default is the 97.5th percentile when \code{coverage} is set to its default value 95.}
#' \item{rr_cv_lower}{The ((100 - \code{coverage}) / 2)th percentile for the risk ratio at the last time point. Default is the 2.5th percentile when \code{coverage} is set to its default value 95.}
#' \item{rd_cv_lower}{The ((100 - \code{coverage}) / 2)th percentile for the risk difference at the last time point. Default is the 2.5th percentile when \code{coverage} is set to its default value 95.}
#' \item{ice_cv_all_lower}{The ((100 - \code{coverage}) / 2)th percentile for the risk of the defined intervention at all time points. Default is the 2.5th percentile when \code{coverage} is set to its default value 95.}
#' \item{ref_cv_all_lower}{The ((100 - \code{coverage}) / 2)th percentile for the risk of the reference intervention at all time points. Default is the 2.5th percentile when \code{coverage} is set to its default value 95.}
#' \item{ref_ipw_se}{Standard error for the observed risk.}
#' \item{ref_ipw_cv_all_upper}{The (100 - (100 - \code{coverage}) / 2)th percentile for the observed risk at all time points. Default is the 97.5th percentile when \code{coverage} is set to its default value 95.}
#' \item{ref_ipw_cv_all_lower}{The ((100 - \code{coverage}) / 2)th percentile for the observed risk at all time points. Default is the 2.5th percentile when \code{coverage} is set to its default value 95.}
#' \item{boot_data}{A list of data samples from all bootstrap replicates.}
#' \item{outcome_init}{A list, where each sublist contains the fitted outcome models in the first step of algorithm for the defined intervention from each bootstrap replicate.}
#' \item{comp_init}{A list, where each sublist contains the fitted competing models in the first step of algorithm for the defined intervention from each bootstrap replicate (if applicable). }
#' \item{np_model}{A list, where each sublist contains the fitted censoring and/or competing models in estimating the observed risk from each bootstrap replicate (if applicable).}
#' \item{outcome_by_step}{A list, where each sublist contains the fitted outcome models in each iteration of algorithm for the defined intervention from each bootstrap replicate.}
#' \item{comp_by_step}{A list, where each sublist contains the fitted competing models in each iteration of algorithm for the defined intervention from each bootstrap replicate (if applicable).}
#' \item{hazard_by_step}{A list, where each sublist contains the fitted hazard model, either time-specific models at each time point or one pooled-over-time global model, for the defined intervention from each bootstrap replicate (if applicable).}
#' \item{ref_outcome_init}{A list, where each sublist contains the fitted outcome models in the first step of algorithm for the reference intervention from each bootstrap replicate.}
#' \item{ref_comp_init}{A list, where each sublist contains the fitted competing models in the first step of algorithm for the reference intervention from each bootstrap replicate (if applicable).}
#' \item{ref_outcome_by_step}{A list, where each sublist contains the fitted outcome models in each iteration of algorithm for the reference intervention from each bootstrap replicate.}
#' \item{ref_comp_by_step}{A list, where each sublist contains the fitted competing models in each iteration of algorithm for the reference intervention from each bootstrap replicate (if applicable).}
#' \item{ref_hazard_by_step}{A list, where each sublist contains the fitted hazard model, either time-specific models at each time point or one pooled-over-time global model, for the reference intervention from each bootstrap replicate (if applicable).}
#' \item{ref_data_err}{A list of logical values indicating whether there is any data error for the reference intervention from each bootstrap replicate.}
#' \item{ref_model_err}{A list of logical valules indicating whether there is any model error produced from each bootstrap replicate for the reference intervention.}
#' \item{ref_model_err_mssg}{A list of strings for error messages of any model error from each bootstrap replicate for the reference intervention.}
#' \item{ref_data_err_mssg}{A list of strings for error messages of any data error from each bootstrap replicate for the reference intervention.}
#'
#' @import doParallel parallel foreach magrittr dplyr stringr
#' @importFrom data.table SJ
#' @importFrom data.table as.data.table
#' @importFrom data.table setkey
#' @importFrom stats na.omit
#' @importFrom stats sd
#' @importFrom stats quantile
#' @noRd
bootstrap_ice <- function(f, K, nboot, coverage, parallel, ncores, ref_description,
                          ref_intervention_varnames, total_effect,
                          ref_intervention,
                          interventions, intervention_varnames, intervention_description,
                          intervention_times, ref_intervention_times,
                          data, id, set_seed, ...) {

  set.seed(set_seed)
  
  sig_level <- 100 - coverage

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
  
  if (parallel == TRUE) {
    
    old_dir <- getwd()

    cl <- makeCluster(ncores)
    registerDoParallel(cl)

    on.exit({
    try(parallel::stopCluster(cl), silent = TRUE)
    foreach::registerDoSEQ()
    }, add = TRUE)

    result <- matrix(NA, ncol = 4, nrow = nboot)
    result <- as.data.frame(result)

    boot_result <- foreach(i = 1:nboot, .combine = "c", .packages = c("tidyverse", "nnet", "stringr", "base"),
                           .errorhandling='pass') %do% {
                             
     unique_idx <- unique(data[, id])
     nidx <- length(unique_idx)
     select_ids <- as.data.table(sample(1:nidx, nidx, replace = TRUE))
     select_ids[, "new_id"] <- 1:nidx
     colnames(select_ids)[1] <- id
     boot_data <- data
     boot_data <- as.data.table(boot_data)
     setkey(boot_data, id)
     boot_data <- boot_data[SJ(select_ids), allow.cartesian = TRUE]
     boot_data <- as.data.frame(boot_data)
     boot_data[, id] <- boot_data[, "new_id"]
     row_idx <- match_ids(as.vector(select_ids[, id]), data[, id])


      ice <- try(
        f(data = boot_data, total_effect = total_effect, id = id, K = K,
               interventions = interv_boot(interventions, row_idx), intervention_names = intervention_varnames,
               intervention_description = intervention_description, intervention_times = intervention_times, ...),
        silent = TRUE
      )


      ref <- try(
        f(data = boot_data, total_effect = total_effect, id = id, K = K,
               interventions = interv_boot(ref_intervention, row_idx), intervention_names = ref_intervention_varnames,
               intervention_description = ref_description, intervention_times = ref_intervention_times, ...),
        silent = TRUE
      )

      if (inherits(ice, "try-error")) {
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
        fit_all <- rep(NA, K+1)
        rr <- NA
        rd <- NA
      } else {
        data_issue_err <- data_issue_mssg <- model_issue_err <- model_issue_mssg <- NA
        
        ice_value <- ice$gformula_risk_last_time
        ice_all <- ice$gformula_risk
        
        this_outcome_init <- list(ice$outcome_init)
        this_comp_init <- list(ice$comp_init)
        this_np_model <- list(ice$np_model)
        
        this_ref_outcome_init <- list(ref$outcome_init)
        this_ref_comp_init <- list(ref$comp_init)
        
        names(this_outcome_init) <- names(this_np_model) <- names(this_comp_init) <- paste0("Bootstrap Replicate ", i)
        names(this_ref_outcome_init) <- names(this_ref_comp_init) <- paste0("Bootstrap Replicate ", i)
        
        this_outcome_by_step <- get_models(ice, "outcome_by_step", paste0("Bootstrap Replicate ", i))
        this_hazard_by_step <- get_models(ice, "hazard_by_step", paste0("Bootstrap Replicate ", i))
        this_comp_by_step <- get_models(ice, "comp_by_step", paste0("Bootstrap Replicate ", i))
        
        this_ref_outcome_by_step <- get_models(ref, "outcome_by_step", paste0("Bootstrap Replicate ", i))
        this_ref_hazard_by_step <- get_models(ref, "hazard_by_step", paste0("Bootstrap Replicate ", i))
        this_ref_comp_by_step <- get_models(ref, "comp_by_step", paste0("Bootstrap Replicate ", i))
        
        ref_all <- ref$gformula_risk
        ref_ipw_all <- c(0, ref$weight_h$risk)
        ref_value <- ref$gformula_risk_last_time
        rr <- ice_value / ref_value
        rd <- ice_value - ref_value
      }
      
      if (inherits(ref, "try-error")) {
        
        if (str_detect(ref[1], "Argument mu must be a nonempty numeric vector")) {
          ref_data_issue_err <- i
          ref_data_issue_mssg <- ref[1]
          
          ref_model_issue_err <- NA
          ref_model_issue_mssg <- NA
          
        } else {
          
          ref_data_issue_err <- NA
          ref_data_issue_mssg <- NA
          
          ref_model_issue_err <- i
          ref_model_issue_mssg <- ref[1]
        }
        
        this_ref_outcome_init <- this_ref_comp_init <- this_ref_outcome_by_step <- 
          this_ref_hazard_by_step <- this_ref_comp_by_step <- NA 
        ref_all <- ref_ipw_all <- ref_value <- NA
        
      } else {
        
        ref_data_issue_err <- ref_data_issue_mssg <- ref_model_issue_err <- ref_model_issue_mssg <- NA
        
        this_ref_outcome_init <- list(ref$outcome_init)
        this_ref_comp_init <- list(ref$comp_init)
        
        names(this_ref_outcome_init) <- names(this_ref_comp_init) <- paste0("Bootstrap Replicate ", i)
        
        this_ref_outcome_by_step <- get_models(ref, "outcome_by_step", paste0("Bootstrap Replicate ", i))
        this_ref_hazard_by_step <- get_models(ref, "hazard_by_step", paste0("Bootstrap Replicate ", i))
        this_ref_comp_by_step <- get_models(ref, "comp_by_step", paste0("Bootstrap Replicate ", i))
        
        ref_all <- ref$gformula_risk
        ref_ipw_all <- c(0, ref$weight_h$risk)
        ref_value <- ref$gformula_risk_last_time
      }

      numerical_out <- rbind(c(ice_value, ref_value, rr, rd, ice_all, ref_all, ref_ipw_all))

      output <- list(boot_out = numerical_out, data = boot_data, data_issue = data_issue_err,
                     model_issue = model_issue_err, 
                     data_issue_mssg = data_issue_mssg, model_issue_mssg = model_issue_mssg,
                     ref_data_issue = ref_data_issue_err,
                     ref_model_issue = ref_model_issue_err, 
                     ref_data_issue_mssg = ref_data_issue_mssg, ref_model_issue_mssg = ref_model_issue_mssg,
                     outcome_init_boot = this_outcome_init, comp_init_boot = this_comp_init,
                     np_model_boot = this_np_model,
                     outcome_by_step_boot = this_outcome_by_step, hazard_by_step_boot = this_hazard_by_step,
                     comp_by_step_boot = this_comp_by_step, 
                     ref_outcome_init_boot = this_ref_outcome_init, ref_comp_init_boot = this_ref_comp_init,
                     ref_outcome_by_step_boot = this_ref_outcome_by_step, ref_hazard_by_step_boot = this_ref_hazard_by_step,
                     ref_comp_by_step_boot = this_ref_comp_by_step)

                           }
    
    message(paste0(" Bootstrap Progress: ", intervention_description, "\n"))

    
    data_err <- stats::na.omit(unlist(boot_result[names(boot_result) == "data_issue"]))
    model_err <- stats::na.omit(unlist(boot_result[names(boot_result) == "model_issue"]))
    
    data_err_mssg <- unlist(boot_result[names(boot_result) == "data_issue_mssg"])
    model_err_mssg <- unlist(boot_result[names(boot_result) == "model_issue_mssg"])
    
    ref_data_err <- stats::na.omit(unlist(boot_result[names(boot_result) == "ref_data_issue"]))
    ref_model_err <- stats::na.omit(unlist(boot_result[names(boot_result) == "ref_model_issue"]))
    
    ref_data_err_mssg <- unlist(boot_result[names(boot_result) == "ref_data_issue_mssg"])
    ref_model_err_mssg <- unlist(boot_result[names(boot_result) == "ref_model_issue_mssg"])
    
    ## construct error message
    data_err_mssg <- stats::na.omit(data_err_mssg)
    model_err_mssg <- stats::na.omit(model_err_mssg)

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
    
    ref_data_err_mssg <- stats::na.omit(ref_data_err_mssg)
    ref_model_err_mssg <- stats::na.omit(ref_model_err_mssg)
    
    ref_data_err_mssg_combine <- paste0("The error message in bootstrap sample ", ref_data_err)
    ref_data_err_mssg_combine <- paste0(ref_data_err_mssg_combine, ": ", ref_data_err_mssg, "\n")
    
    ref_model_err_mssg_combine <- paste0("The error message in bootstrap sample ", ref_model_err)
    ref_model_err_mssg_combine <- paste0(ref_model_err_mssg_combine, ": ", ref_model_err_mssg, "\n")
    
    outcome_init <- get_boot_models("outcome_init_boot", boot_result)
    comp_init <- get_boot_models("comp_init_boot", boot_result)
    np_model <- get_boot_models("np_model_boot", boot_result)
    outcome_by_step <- get_boot_models("outcome_by_step_boot", boot_result)
    comp_by_step <- get_boot_models("comp_by_step_boot", boot_result)
    hazard_by_step <- get_boot_models("hazard_by_step_boot", boot_result)
    
    ref_outcome_init <- get_boot_models("ref_outcome_init_boot", boot_result)
    ref_comp_init <- get_boot_models("ref_comp_init_boot", boot_result)
    ref_outcome_by_step <- get_boot_models("ref_outcome_by_step_boot", boot_result)
    ref_comp_by_step <- get_boot_models("ref_comp_by_step_boot", boot_result)
    ref_hazard_by_step <- get_boot_models("ref_hazard_by_step_boot", boot_result)
    
    result <- matrix(unlist(boot_result[names(boot_result) == "boot_out"]), byrow = TRUE, nrow = nboot)

    agg_result <- apply(result, 2, stats::sd, na.rm = TRUE)
    ice_se <- agg_result[1]
    ref_se <- agg_result[2]
    rr_se <- agg_result[3]
    rd_se <- agg_result[4]
    ice_se_all <- agg_result[5:(5+K)]
    ref_se_all <- agg_result[(6+K):(6+2*K)]
    ref_ipw_se_all <- agg_result[(7+2*K):(length(agg_result))]


    quantile_result_upper <- apply(result, 2, stats::quantile, probs = (100 - sig_level/2)/100, na.rm = TRUE)
    ice_cv_upper <- quantile_result_upper[1]
    ref_cv_upper <- quantile_result_upper[2]
    rr_cv_upper <- quantile_result_upper[3]
    rd_cv_upper <- quantile_result_upper[4]
    ice_cv_all_upper <- quantile_result_upper[5:(5+K)]
    ref_cv_all_upper <- quantile_result_upper[(6+K):(6+2*K)]
    ref_ipw_cv_all_upper <- quantile_result_upper[(7+2*K):(length(quantile_result_upper))]
    
    quantile_result_lower <- apply(result, 2, stats::quantile, probs = (sig_level/2)/100, na.rm = TRUE)
    ice_cv_lower <- quantile_result_lower[1]
    ref_cv_lower <- quantile_result_lower[2]
    rr_cv_lower <- quantile_result_lower[3]
    rd_cv_lower <- quantile_result_lower[4]
    ice_cv_all_lower <- quantile_result_lower[5:(5+K)]
    ref_cv_all_lower <- quantile_result_lower[(6+K):(6+2*K)]
    ref_ipw_cv_all_lower <- quantile_result_lower[(7+2*K):(length(quantile_result_lower))]

    boot_data_all <- boot_result[names(boot_result) == "data"]
    
    # stopCluster(cl)

  } else {
    boot_data_all <- c()
    outcome_init <- comp_init <- np_model <- outcome_by_step <- comp_by_step <- hazard_by_step <- c()
    ref_outcome_init <- ref_comp_init <- ref_outcome_by_step <- ref_comp_by_step <- ref_hazard_by_step <- c()
    data_err <- model_err <- c()
    data_err_mssg <- model_err_mssg <- c()
    
    ref_data_err <- ref_model_err <- c()
    ref_data_err_mssg <- ref_model_err_mssg <- c()
    for (i in 1:nboot) {

      unique_idx <- unique(data[, id])
      nidx <- length(unique_idx)
      select_ids <- as.data.table(sample(1:nidx, nidx, replace = TRUE))
      select_ids[, "new_id"] <- 1:nidx
      colnames(select_ids)[1] <- id
      boot_data <- data
      boot_data <- as.data.table(boot_data)
      setkey(boot_data, id)
      boot_data <- boot_data[SJ(select_ids), allow.cartesian = TRUE]
      boot_data <- as.data.frame(boot_data)
      boot_data[, id] <- boot_data[, "new_id"]
      
      row_idx <- match_ids(as.vector(select_ids[, id]), data[, id])

      boot_data_all <- c(boot_data_all, list(boot_data))

      ice <- try(
        f(data = boot_data, total_effect = total_effect, id = id, K = K,
               interventions = interv_boot(interventions, row_idx), intervention_names = intervention_varnames,
               intervention_description = intervention_description, intervention_times = intervention_times, ...),
        silent = TRUE
      )

      ref <- try(
        f(data = boot_data, total_effect = total_effect, id = id, K = K,
               interventions = interv_boot(ref_intervention, row_idx), intervention_names = ref_intervention_varnames,
               intervention_description = ref_description, intervention_times = ref_intervention_times, ...), 
        silent = TRUE
      )
      
        if (inherits(ice, "try-error")) {
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
          fit_all <- rep(NA, K+1)
          rr <- NA
          rd <- NA
        } else {
          data_issue_err <- data_issue_mssg <- model_issue_err <- model_issue_mssg <- NA
          
          ice_value <- ice$gformula_risk_last_time
          ice_all <- ice$gformula_risk
          
          this_outcome_init <- list(ice$outcome_init)
          this_comp_init <- list(ice$comp_init)
          this_np_model <- list(ice$np_model)
          
          this_ref_outcome_init <- list(ref$outcome_init)
          this_ref_comp_init <- list(ref$comp_init)
          
          names(this_outcome_init) <- names(this_np_model) <- names(this_comp_init) <- paste0("Bootstrap Replicate ", i)
          names(this_ref_outcome_init) <- names(this_ref_comp_init) <- paste0("Bootstrap Replicate ", i)
          
          this_outcome_by_step <- get_models(ice, "outcome_by_step", paste0("Bootstrap Replicate ", i))
          this_hazard_by_step <- get_models(ice, "hazard_by_step", paste0("Bootstrap Replicate ", i))
          this_comp_by_step <- get_models(ice, "comp_by_step", paste0("Bootstrap Replicate ", i))
          
          this_ref_outcome_by_step <- get_models(ref, "outcome_by_step", paste0("Bootstrap Replicate ", i))
          this_ref_hazard_by_step <- get_models(ref, "hazard_by_step", paste0("Bootstrap Replicate ", i))
          this_ref_comp_by_step <- get_models(ref, "comp_by_step", paste0("Bootstrap Replicate ", i))
          
          ref_all <- ref$gformula_risk
          ref_ipw_all <- c(0, ref$weight_h$risk)
          ref_value <- ref$gformula_risk_last_time
          rr <- ice_value / ref_value
          rd <- ice_value - ref_value
        }
        
        if (inherits(ref, "try-error")) {
          
          if (str_detect(ref[1], "Argument mu must be a nonempty numeric vector")) {
            ref_data_issue_err <- i
            ref_data_issue_mssg <- ref[1]
            
            ref_model_issue_err <- NA
            ref_model_issue_mssg <- NA
            
          } else {
            
            ref_data_issue_err <- NA
            ref_data_issue_mssg <- NA
            
            ref_model_issue_err <- i
            ref_model_issue_mssg <- ref[1]
          }
          
          this_ref_outcome_init <- this_ref_comp_init <- this_ref_outcome_by_step <- 
            this_ref_hazard_by_step <- this_ref_comp_by_step <- NA 
          ref_all <- ref_ipw_all <- ref_value <- NA
          
        } else {
          
          ref_data_issue_err <- ref_data_issue_mssg <- ref_model_issue_err <- ref_model_issue_mssg <- NA
          
        this_ref_outcome_init <- list(ref$outcome_init)
        this_ref_comp_init <- list(ref$comp_init)
        
        names(this_ref_outcome_init) <- names(this_ref_comp_init) <- paste0("Bootstrap Replicate ", i)
        
        this_ref_outcome_by_step <- get_models(ref, "outcome_by_step", paste0("Bootstrap Replicate ", i))
        this_ref_hazard_by_step <- get_models(ref, "hazard_by_step", paste0("Bootstrap Replicate ", i))
        this_ref_comp_by_step <- get_models(ref, "comp_by_step", paste0("Bootstrap Replicate ", i))
        
        ref_all <- ref$gformula_risk
        ref_ipw_all <- c(0, ref$weight_h$risk)
        ref_value <- ref$gformula_risk_last_time
        }
      
      
      outcome_by_step <- c(outcome_by_step, this_outcome_by_step)
      hazard_by_step <- c(hazard_by_step, this_hazard_by_step)
      comp_by_step <- c(comp_by_step, this_comp_by_step)
      
      outcome_init <- c(outcome_init, this_outcome_init)
      comp_init <- c(comp_init, this_comp_init)
      np_model <- c(np_model, this_np_model)
      
      ref_outcome_by_step <- c(ref_outcome_by_step, this_ref_outcome_by_step)
      ref_hazard_by_step <- c(ref_hazard_by_step, this_ref_hazard_by_step)
      ref_comp_by_step <- c(ref_comp_by_step, this_ref_comp_by_step)
      
      ref_outcome_init <- c(ref_outcome_init, this_ref_outcome_init)
      ref_comp_init <- c(ref_comp_init, this_ref_comp_init)

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
      
      if (i %% nshow == 0) {
        message(paste0(intervention_description, " Bootstrap Progress: ", i, "/", nboot, "\n"))
      }

    }
    
    ice_se <- stats::sd(unlist(ice_boot), na.rm = TRUE)
    ice_cv_upper <- stats::quantile(unlist(ice_boot), probs = (100 - sig_level/2)/100, na.rm = TRUE)
    ice_cv_lower <- stats::quantile(unlist(ice_boot), probs = (sig_level/2)/100, na.rm = TRUE)
    ref_se <- stats::sd(unlist(ref_boot), na.rm = TRUE)
    ref_cv_upper <- stats::quantile(unlist(ref_boot), probs = (100 - sig_level/2)/100, na.rm = TRUE)
    ref_cv_lower <- stats::quantile(unlist(ref_boot), probs = (sig_level/2)/100, na.rm = TRUE)
    rr_se <- stats::sd(unlist(rr_boot), na.rm = TRUE)
    rr_cv_upper <- stats::quantile(unlist(rr_boot), probs = (100 - sig_level/2)/100, na.rm = TRUE)
    rr_cv_lower <- stats::quantile(unlist(rr_boot), probs = (sig_level/2)/100, na.rm = TRUE)
    rd_se <- stats::sd(unlist(rd_boot), na.rm = TRUE)
    rd_cv_upper <- stats::quantile(unlist(rd_boot), probs = (100 - sig_level/2)/100, na.rm = TRUE)
    rd_cv_lower <- stats::quantile(unlist(rd_boot), probs = (sig_level/2)/100, na.rm = TRUE)
    ice_se_all <- apply(ice_boot_all, 2, stats::sd, na.rm = TRUE)
    ref_se_all <- apply(ref_boot_all, 2, stats::sd, na.rm = TRUE)
    ref_ipw_se_all <- apply(ref_ipw_boot_all, 2, stats::sd, na.rm = TRUE)
    
    ice_cv_all_upper <- apply(ice_boot_all, 2, stats::quantile, probs = (100 - sig_level/2)/100, na.rm = TRUE)
    ref_cv_all_upper <- apply(ref_boot_all, 2, stats::quantile, probs = (100 - sig_level/2)/100, na.rm = TRUE)
    ref_ipw_cv_all_upper <- apply(ref_ipw_boot_all, 2, stats::quantile, probs = (100 - sig_level/2)/100, na.rm = TRUE)
    
    ice_cv_all_lower <- apply(ice_boot_all, 2, stats::quantile, probs = (sig_level/2)/100, na.rm = TRUE)
    ref_cv_all_lower <- apply(ref_boot_all, 2, stats::quantile, probs = (sig_level/2)/100, na.rm = TRUE)
    ref_ipw_cv_all_lower <- apply(ref_ipw_boot_all, 2, stats::quantile, probs = (sig_level/2)/100, na.rm = TRUE)
    
    data_err <- stats::na.omit(data_err)
    model_err <- stats::na.omit(model_err)
    
    data_err_mssg <- stats::na.omit(data_err_mssg)
    model_err_mssg <- stats::na.omit(model_err_mssg)
    
    ref_data_err <- stats::na.omit(ref_data_err)
    ref_model_err <- stats::na.omit(ref_model_err)
    
    ref_data_err_mssg <- stats::na.omit(ref_data_err_mssg)
    ref_model_err_mssg <- stats::na.omit(ref_model_err_mssg)
    
    ## construct error message
    data_err_mssg_combine <- paste0("The error message in bootstrap sample ", data_err)
    data_err_mssg_combine <- paste0(data_err_mssg_combine, ": ", data_err_mssg, "\n")

    model_err_mssg_combine <- paste0("The error message in bootstrap sample ", model_err)
    model_err_mssg_combine <- paste0(model_err_mssg_combine, ": ", model_err_mssg, "\n")
    
    give_warning(data_err, data_err_mssg_combine)
    give_warning(model_err, model_err_mssg_combine)
    
    ## construct error message for reference intervention
    ref_data_err_mssg_combine <- paste0("The error message in bootstrap sample ", ref_data_err)
    ref_data_err_mssg_combine <- paste0(ref_data_err_mssg_combine, ": ", ref_data_err_mssg, "\n")
    
    ref_model_err_mssg_combine <- paste0("The error message in bootstrap sample ", ref_model_err)
    ref_model_err_mssg_combine <- paste0(ref_model_err_mssg_combine, ": ", ref_model_err_mssg, "\n")
    

  }

  return(list(ice_se = ice_se_all, ref_se = ref_se_all, rr_se = rr_se, rd_se = rd_se,
              rr_cv_upper = rr_cv_upper, rd_cv_upper = rd_cv_upper,
              ice_cv_all_upper = ice_cv_all_upper, ref_cv_all_upper = ref_cv_all_upper,
              rr_cv_lower = rr_cv_lower, rd_cv_lower = rd_cv_lower,
              ice_cv_all_lower = ice_cv_all_lower, ref_cv_all_lower = ref_cv_all_lower,
              ref_ipw_se = ref_ipw_se_all, ref_ipw_cv_all_upper = ref_ipw_cv_all_upper, ref_ipw_cv_all_lower = ref_ipw_cv_all_lower,
              boot_data = boot_data_all, outcome_init = outcome_init, 
              comp_init = comp_init, np_model = np_model, outcome_by_step = outcome_by_step, 
              comp_by_step = comp_by_step, hazard_by_step = hazard_by_step, 
              ref_outcome_init = ref_outcome_init, 
              ref_comp_init = ref_comp_init, 
              ref_outcome_by_step = ref_outcome_by_step, 
              ref_comp_by_step = ref_comp_by_step, ref_hazard_by_step = ref_hazard_by_step, 
              ref_data_err = ref_data_err, ref_model_err = ref_model_err, 
              ref_model_err_mssg = ref_model_err_mssg, ref_data_err_mssg = ref_data_err_mssg))
}
