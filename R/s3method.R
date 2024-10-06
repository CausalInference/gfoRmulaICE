#' Summary method for ICE Estimator Objects
#'
#' This function returns a summary table for ICE estimator objects. 
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
#' summary(fit_classical_pool, fit_hazard_pool)
summary.ICE <- function(...) {
  
  fit_ice <- list(...)
  
  if (length(fit_ice) > 1) {
  for (i in 1:length(fit_ice)) {
    if (!inherits(fit_ice[[i]], "ICE")){
      stop("Input arguments must be objects of ICE class.")
    }
  }
  } else {
    if (!inherits(fit_ice[[1]], "ICE")){
      stop("Input arguments must be objects of ICE class.")
    }
  }
  
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


#' Plot method for ICE estimator objects
#'
#' This function provides visualization of estimated risk 
#' for all specified interventions, estimated natural course risk, and observed risk at each time point.
#'
#' @param ... ICE estimator objects.
#' @param plot_obs a logical value indicating whether to plot the observed risk over time.
#' Default is \code{TRUE}.
#' @param label a number specifying which time label is used in x-axis. 0 represents using generic numerical time index, and
#' 1 represents using the original time label from the data. Default is 0.
#'
#' @return a plot for risks of all the interventions specified in \code{...}.
#' @export
#' @import ggplot2
#'
#' @examples
#' 
#' data <- gfoRmulaICE::compData
#'
#' ice_fit1 <- ice(
#' data = data, 
#' time_points = 4, 
#' id = "id", 
#' time_name = "t0",
#' censor_name = "C", 
#' outcome_name = "Y",
#' compevent_name = "D",
#' comp_effect = 0,
#' outcome_model = Y ~ L1 + L2 + A1 + A2, 
#' censor_model = C ~ L1 + L2 + A1 + A2,
#' ref_idx = 0,
#' estimator = pool(hazard = F),
#' int_descript = "Static Intervention",
#' intervention1.A1 = list(static(3)),
#' intervention1.A2 = list(static(1))
#' )
#'
#' ice_fit2 <- ice(
#' data = data, 
#' time_points = 4, 
#' id = "id", 
#' time_name = "t0",
#' censor_name = "C", 
#' outcome_name = "Y",
#' compevent_name = "D",
#' comp_effect = 0,
#' outcome_model = Y ~ L1 + L2 + A1 + A2, 
#' censor_model = C ~ L1 + L2 + A1 + A2,
#' competing_model = D ~ L1 + L2 + A1 + A2,
#' ref_idx = 0,
#' estimator = pool(hazard = T),
#' int_descript = "Static Intervention",
#' intervention1.A1 = list(static(3)),
#' intervention1.A2 = list(static(1))
#' )
#'
#' plot(ice_fit1, ice_fit2)
#'
plot.ICE <- function(..., plot_obs = T, label = 0) {
  
  plot_np <- plot_obs
  
  fit_ice <- list(...)
  
  if (length(fit_ice) > 1) {
    for (i in 1:length(fit_ice)) {
      if (!inherits(fit_ice[[i]], "ICE")){
        stop("Input arguments must be objects of ICE class.")
      }
    }
  } else {
    if (!inherits(fit_ice[[1]], "ICE")){
      stop("Input arguments must be objects of ICE class.")
    }
  }
  
  if (length(fit_ice) >= 1) {
    risk_df <- data.frame()
    for (i in 1:length(fit_ice)) {
      if (length(fit_ice[[i]]$estimator.type) > 0) {
        
        risk_fit <- fit_ice[[i]]$risk.over.time
        estimator_name <- fit_ice[[i]]$estimator.type
        risk_fit$estimator_name <- estimator_name
        risk_df <- rbind(risk_df, risk_fit)
        
      } else {
        stop("Please input a valid ICE model object")
      }
    }
  } else {
    stop("Please input a valid ICE model object.")
  }
  
  ## pre-process
  
  risk_df$Intervention <- ifelse(str_detect(risk_df$Intervention, "Natural Course") & !str_detect(risk_df$Intervention, "nonparametric"), 
                                 paste0("Natural Course (parametric ICE)"), #, str_split(risk_df$Intervention, "Natural Course")[[1]][2]), 
                                 risk_df$Intervention)
  
  if (!plot_np) {
    risk_df <- risk_df %>% filter(Intervention != "Natural Course (nonparametric)")
  }
  
  if (label == 0) {
    xlabels <- risk_df$Time
  } else if (label == 1) {
    xlabels <- risk_df$originTime
  }
  ## plot confidence interval if bootstrap
  
  if ("SE" %in% colnames(risk_df)) {
    risk_df$Risk <- as.numeric(risk_df$Risk)
    risk_df <- risk_df %>% mutate(lower = Risk - Critical_Value_Lower * SE,
                                  upper = Risk + Critical_Value_Upper * SE)
    risk_plot <- risk_df %>% ggplot(aes(as.numeric(Time), as.numeric(Risk), colour = as.factor(Intervention))) +
      geom_point() + geom_line() +
      geom_errorbar( aes(ymin = lower, ymax = upper),width = 0.2) +
      facet_wrap(~estimator_name) + scale_x_continuous(breaks = risk_df$Time, labels = xlabels) +
      xlab("Time") + ylab("Risk") + ggtitle("ICE Estimator Using Different Interventions") +
      guides(color = guide_legend(title = "Intervention"))
  } else {
    
    risk_plot <- risk_df %>% ggplot(aes(as.numeric(Time), as.numeric(Risk), colour = as.factor(Intervention))) +
      geom_point() + geom_line() +
      facet_wrap(~estimator_name) + scale_x_continuous(breaks = risk_df$Time, labels = xlabels) +
      xlab("Time") + ylab("Risk") + ggtitle("ICE Estimator Using Different Interventions") +
      guides(color = guide_legend(title = "Intervention"))
  }
  
  return(risk_plot)
  
}
