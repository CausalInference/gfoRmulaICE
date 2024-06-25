#' Plot risk estimated by ICE estimator over time
#'
#' This function provides visualization of estimated risk on all user-defined interventions
#' from the fitted ICE estimator object with natural course risk over time.
#' (This function is going to be converted to the S3 method in R so it is named "plot_risk" for now.
#' After conversion, users could use the base function plot().)
#'
#' @param ... ICE estimator objects.
#' @param plot_obs a logical indicating whether to add the observed risk over time to the plot.
#' Default is TRUE
#' @param label a numeric specifying which time label is used in x axis. 0 represents using generic numerical time index, and
#' 1 represents using the original time index in the data set. Default is 0.
#'
#' @return a plot for risks of all the interventions specified in \code{...}.
#' @export
#' @import ggplot2
#'
#' @examples
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
#' competing_model = D ~ L1 + L2 + A1 + A2,
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
#' plot_risk(ice_fit1, ice_fit2)
#'
plot_risk <- function(..., plot_np = T, label = 0) {
  
  plot_np <- plot_obs

  fit_ice <- list(...)

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
                                 paste0("Natural Course (parametric ICE)", str_split(risk_df$Intervention, "Natural Course")[[1]][2]), 
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
