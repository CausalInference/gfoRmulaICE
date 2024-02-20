#' Singly Robust and Doubly Robust Iterative Conditional Expectation (ICE) Estimator for the Specified Intervention
#'
#' This function estimates the risk over time for survival outcome using the given observed data set following
#' multiple user-defined intervention strategies by the parametric g-formula iterative conditional expectation (ICE)
#' estimator. This function allows users to access all singly robust and doubly robust ICE estimators:
#' classical pooling over treatment history ICE estimator, classical stratifying on treatment history ICE estimator,
#' hazard-based pooling over treatment history ICE estimator, hazard-based stratifying on treatment history ICE estimator,
#' and a doubly robust inverse probability weighted ICE estimator. Please see Wen et al. (2021) more details regarding the parameter g-formula
#' iterative conditional expectation estimator.
#'
#' Users could specify which version ICE estimator to use through \code{estimator}.
#' \itemize{
#' \item{\code{pool(hazard = F)} specifies the classical pooling over treatment history ICE estimator. }
#' \item{\code{pool(hazard = T)} specifies the hazard-based pooling over treatment history ICE estimator. }
#' \item{\code{strat(hazard = F)} specifies the classical stratifying on treatment history ICE estimator. }
#' \item{\code{strat(hazard = T)} specifies the hazard-based stratify on treatment history ICE estimator. }
#' \item{\code{weight(treat_model)} specifies the IP weighted ICE estimator
#' where \code{treat_model} specifies the treatment model.}
#' }
#' For stratified ICE, the estimator will automatically add in the treatment variables that are not intervened as covariates
#' into the model specification for each time point. To provide flexible choices on the model specification inputs, we also
#' provide keyword argument options for users to specify a different model statement, which is different from the model statement specified
#' in \code{outcome_model} or \code{competing_model}, for user-chosen interventions.
#' The input to each outcome or competing model keyword argument is a model statement formula.
#' To specify different outcome model for different intervention, please follow the keyword argument input convention below:
#'
#' \itemize{
#' \item{Each outcome model is specified using keyword argument name starting with \emph{outcomeModel} prefix.}
#' \item{Use \emph{.n} \emph{outcomeModel} prefix in keyword argument name to specify which intervention being applied to,
#' where \emph{n} represents the \emph{n}th intervention.}
#' }
#'
#' To specify different competing model for different intervention, please follow the keyword argument input convention below:
#'
#' \itemize{
#' \item{Each competing model is specified using keyword argument name starting with \emph{compModel} prefix.}
#' \item{Use \emph{.n} \emph{compModel} prefix in keyword argument name to represent applying to \emph{n}th intervention.}
#' }
#'
#' If an outcome or competing model is specified for an intervention using keyword argument, then the outcome or competing model will be used for the intervention.
#' If there is no outcome or competing model keyword argument specified, the outcome or competing model specified in \code{outcome_model} or \code{competing_model} will be
#' used for the intervention.
#'
#' For example, the following input specifies \code{Y ~ L1} as outcome model for intervention 1,
#' \code{D ~ L1} as competing model for intervention 2,
#' \code{D ~ L1 + L2} as competing model for intervention 1,
#' \code{Y ~ L1 + L2} as outcome model for intervention 2.
#'
#' \code{fit_hazard_strat <- ice(data = data, K = 4, id = "id", time_name = "t0", } \cr
#' \code{outcome_name = "Y", censor_name = "C", competing_name = "D", } \cr
#' \code{estimator = strat(hazard = T), comp_effect = 1, } \cr
#' \code{censor_model = C ~ L1 + L2, ref_idx = 0,} \cr
#' \code{int_descript = c("Static Intervention", "Dynamic Intervention"),} \cr
#' \code{outcome_model = Y ~ L1 + L2,} \cr
#' \code{competing_model = D ~ L1 + L2,} \cr
#' \code{intervention1.A1 = list(static(0)),} \cr
#' \code{intervention1.A2 = list(static(1)),} \cr
#' \code{intervention2.A1 = list(dynamic("absorbing", "L1", "=", 1)),} \cr
#' \code{intervention2.A2 = list(dynamic("compare", "L1", "=", 0)),} \cr
#' \code{outcomeModel.1 = Y ~ L1,} \cr
#' \code{compModel.2 = D ~ L1} \cr
#'
#' Because the keyword argument for outcome model is not specified for intervention 2, the outcome model for intervention 2 is
#' \code{Y ~ L1 + L2} as specified in \code{outcome_model}.
#' Similarly, because the keyword argument for competing model is not specified for intervention 1, the competing model for intervention 1 is
#' \code{D ~ L1 + L2} as specified in \code{competing_model}.
#' Please see more examples in the examples section.
#' Note that for stratified ICE in the case of direct effect, the keyword argument competing model statement inputs are ignored since 
#' the competing model is used for nonparametric risk estimation.
#'
#' Users could specify user-defined interventions or implement built-in interventions provided by the package
#' using the intervention input convention described in the parameter description section.
#'
#' For the package built-in intervention strategy functions, if one wants to implement the strategy within the ICE method,
#' please follow the instructions described below. Furthermore, users could assess the intervened values by
#' specifying additional arguments following the detailed documentation for each function.
#' Built-in interventions include following:
#' \itemize{
#' \item{Always Treat:} {\code{static(1)}}
#' \item{Never Treat:} {\code{static(0)}}
#' \item{Dynamic Treat Based on \code{var}:} {\code{dynamic(type, var, direction, value)}. \code{type} could be "compare" or "absorbing" and \code{direction} could be \code{>, >=, <, <=, =, !=}.
#' If \code{type} is "compare", then treat if \code{var} > or >= or < or <= or = or != \code{value} and not treat otherwise. If \code{type} is "absorbing", then implements always treat if treatment is initiated
#' based on the description of the "compare" type.}
#' \item{Threshold Intervention Based on \code{var}:} {\code{threshold(lower_bound, upper_bound, var)}. If treatment value is within the range between \code{lower_bound} and \code{upper_bound} inclusively,
#' then follow the natural value of the treatment. Otherwise, if the treatment value is below \code{lower_bound}, then set to \code{lower_bound}.
#' If the treatment value is above the upper bound, then set to \code{upper_bound}.}
#' \item{Grace Period:} {\code{grace_period(type, nperiod, var, value)}. \code{type} could be "uniform" or "natural". \code{nperiod} specifies the length of grace period.
#' When a defined intervention condition based on \code{var} = \code{value} is met, initiate treatment in \code{nperiod} time units. If there is no intervention,
#' follow the observed treatment initiation distribution (for natural grace period) or the uniform distribution of treatment initiation (for uniform grace period). }
#' \item{User-Defined Intervention:} {The output of the user-defined intervention should be a vector of intervened value of the intervention variable for each individual at each time point in the size as the number of rows in \code{data}.}
#' }
#' Some examples are provided in the example section.
#'
#' In order to obtain an inverse probability (IP) weighted natural course risk based on the observed data, users must specify
#' a censoring variable through \code{censor_name} and a corresponding censoring model through \code{censor_model}. Please see
#' Chiu et al. (2023) for more details regarding the IP weighted estimate of the natural course risk.
#'
#' If competing event exists in the data, users need to specify the name of the competing variable through \code{competing_name} and
#' the model specification through \code{competing_model} for hazard-based ICE estimator. Users need to specify whether to treat
#' the competing event as censoring or total effect through \code{total_effect}.
#'
#'
#'
#' @param data a data frame containing the observed data in long format.
#' @param time_points a numerical value that indicates the total number of time points.
#' @param id a character string indicating the ID variable name in \code{data}.
#' @param time_name a character string indicating the time variable name in \code{data}.
#' @param outcome_name a character string indicating the outcome variable name in \code{data}.
#' @param censor_name a character string indicating the censor variable name in \code{data}. Default is \code{NULL}.
#' @param compevent_name a character string indicating the competing variable name in \code{data}. Default is \code{NULL}.
#' @param comp_effect a numeric specifying how the competing event is handled for all the specified interventions. Default is 0.
#' \code{0} outputs direct effect for all the specified interventions.
#' \code{1} outputs total effect for all the specified interventions.
#' @param outcome_model a formula specifying the model statement for the outcome model.
#' @param censor_model a formula specifying the model statement for the censoring model for IP weighted natural course risk. Default is \code{NULL}.
#' @param competing_model a formula specifying the model statement for the competing model for hazard-based ICE estimator. Default is \code{NULL}.
#' @param ref_idx a numerical indicating which intervention to be used as the reference to calculate the risk ratio and risk difference. Default is 0.
#' 0 refers to using the natural course as the reference intervention.
#' Any other numbers refer to the corresponding intervention that users specify in the keywords argument.
#' @param estimator a function to specifying which ICE estimator to use for the estimation. Possible inputs are:
#' \itemize{
#' \item{Classical Pooling over Treatment History ICE Estimator: }{\code{pool(hazard = F)}}
#' \item{Hazard Based Pooling over Treatment History ICE Estimator: }{\code{pool(hazard = T)}}
#' \item{Classical Stratifying Treatment History ICE Estimator: }{\code{strat(hazard = F)}}
#' \item{Hazard Based Stratifying Treatment History ICE Estimator: }{\code{strat(hazard = T)}}
#' \item{Doubly Robust Weighted ICE Estimator: }{\cr
#' \code{weight(treat_model)}
#' where \code{treat_model} is a list specifying the treatment model.}
#' }
#' @param int_descript a vector of strings containing description for each specified intervention.
#' @param ci_method a character string specifying the method used for calculating the confidence interval, if \code{nsamples} is larger than 0.
#' Either "percentile" or "normal." Default is "percentile."
#' @param nsamples a numeric indicating number of bootstrap samples. If a number larger than 0 is specified, bootstrap samples will be used for
#' standard error estimate and confidence interval. Default is 0.
#' @param seed a numeric indicating the starting seed for bootstrap. Default is 1.
#' @param significance_level a numeric indicating the significance level to be used for confidence interval. Default is 0.05.
#' @param parallel a logical indicating whether to parallelize the bootstrap process. TRUE for using parallel computing. FALSE for not using parallel computing.
#' Default is FALSE.
#' @param ncores a numeric indicating the number of CPU cores to be used in parallel computing. Default is 2.
#' @param ... keywords arguments for intervention inputs, optional outcome models for stratified ICE, and optional competing models for stratified ICE.
#' The keyword argument for interventions should follow the below naming convention:
#' \itemize{
#' \item{Each intervention is specified using the keyword argument name with \emph{intervention} prefix.}
#' \item{Use \emph{i} after \emph{intervention} prefix in keyword argument name to represent the ith intervention strategy.}
#' \item{Use \emph{.} followed with \emph{treatment variable name} after \emph{interventioni} in keyword argument name to represent the treatment name within the ith intervention strategy.}
#' }
#' The input for each keyword argument must be a list encompassing two elements, which are:
#' a vector of intervened values and an optional vector of time points on which the corresponding intervention is applied.
#' If the second element is not specified, then the defined intervention will be applied to all time points.
#' A sample intervention input with two treatments of names A1 and A2 and two intervention strategies look like: \cr
#' \code{intervention1.A1 = list(static(1))} \cr
#' \code{intervention1.A2 = list(dynamic("compare", "L1", ">", 0))} \cr
#' \code{intervention2.A1 = list(static(0))} \cr
#' \code{intervention2.A2 = list(dynamic("absorbing", "L1", ">", 1))} \cr
#' A sample intervention input with one treatment of name A and two intervention strategies look like: \cr
#' \code{intervention1.A1 = list(static(1))} \cr
#' \code{intervention2.A1 = list(static(0))} \cr
#' For the above two intervention inputs, since there is no second element input, each specified intervention is applied on all time points.
#' If one wants to specify the custom time points on which each intervention is applied, a sample intervention input looks like:
#' \code{intervention1.A1 = list(static(1), 1:2)} \cr
#' \code{intervention1.A2 = list(static(0), 3:5)} \cr
#' where the intervention on A1 within intervention strategy 1 is applied on time 1 and 2,
#' and the intervention on A2 within intervention strategy 1 is applied on time 3 to 5.
#' If there is no intervention specified, only the natural course risk will be returned. \cr
#' To specify different outcome model for different intervention, please follow the keyword argument input convention below:
#'
#' \itemize{
#' \item{Each outcome model is specified using keyword argument name starting with \emph{outcomeModel} prefix.}
#' \item{Use \emph{.n} \emph{outcomeModel} prefix in keyword argument name to specify which intervention being applied to,
#' where \emph{n} represents the \emph{n}th intervention.}
#' }
#'
#' To specify different competing model for different intervention, please follow the keyword argument input convention below:
#'
#' \itemize{
#' \item{Each competing model is specified using keyword argument name starting with \emph{compModel} prefix.}
#' \item{Use \emph{.n} \emph{compModel} prefix in keyword argument name to represent applying to \emph{n}th intervention.}
#' }
#' The input to each outcome or competing model keyword argument is a model statement formula.
#' Please refer to the examples section for more examples.
#'
#'
#'
#' @return A list containing the following components:
#' \item{summary}{A summary table containing the estimated ICE risk, risk ratio, risk difference. If \code{bootstrap} is TRUE, then the table also includes standard error and confidence interval for ICE risk, risk ratio, and risk difference of each intervention.}
#' \item{risk.over.time}{A data frame containing the estimated ICE risk at each time point for each intervention.}
#' \item{models}{A list containing sublists whose names are the specified intervention descriptions, and each sublist contains the fitted models for the outcome, the treatment (if applicable), and the competing event (if applicable) of the corresponding intervention.}
#' \item{model.summary}{A list containing sublists whose names are the specified intervention descriptions, and each sublist contains the summary of the fitted models for the corresponding intervention.}
#' \item{model.stderr}{A list containing sublists whose names are the specified intervention descriptions, and each sublist contains the standard errors of the coefficients of the fitted models for the corresponding intervention.}
#' \item{model.vcov}{A list containing sublists whose names are the specified intervention descriptions, and each sublist contains the variance-covariance matrices of the parameters of the fitted models for the corresponding intervention.}
#' \item{model.rmse}{A list containing sublists whose names are the specified intervention descriptions, and each sublist contains root mean square error (RMSE) values of the fitted models for the corresponding intervention.}
#' \item{boot.data}{A list containing all the bootstrapped data if \code{boostrap} is set to \code{TRUE}.}
#' \item{boot.models}{A list containing the fitted models for the outcome, the treatment (if applicable), and the competing event (if applicable) on the bootstrapped samples.}
#' \item{boot.summary}{A list containing the summary of the fitted models on the bootstrapped samples.}
#' \item{boot.stderr}{A list containing the standard errors of the coefficients of the fitted models on the bootstrapped samples.}
#' \item{boot.vcov}{A list containing the variance-covariance matrices of the parameters of the fitted models on the bootstrapped samples.}
#' \item{boot.rmse}{A list containing the root mean square error (RMSE) values of the fitted models on the bootstrapped samples.}
#' \item{estimator.type}{A string specifying the type of ICE estimator.}
#' @export
#' @import tidyverse stringr data.table reshape2
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
#' @references Chiu YH, Wen L, McGrath S, Logan R, Dahabreh IJ, Hernán MA. Evaluating model specification when using the parametric g-formula in the presence of censoring. American Journal of Epidemiology. 2023;192:1887–1895.
#'
#' @examples
#'
#' data <- gfoRmulaICE::data
#' dynamic_cat <- case_when(data$L2 < 0 ~ 1,
#' data$L2 >= 0 & data$L2 < 2 ~ 2,
#' T ~ 3)
#'
#' # For the following examples, we consider two interventions.
#' # Intervention 1 is a static intervention, and
#' # intervention 2 is a dynamic intervention, where the intervention
#' # for A1 is user-defined.
#'
#' # Hazard based stratified ICE,
#' # competing event as total effect for all interventions,
#' # with Y ~ L1 + L2 as outcome model for both intervention 1 and intervention 2,
#' # and D ~ L1 + L2 as competing model for both intervention 1 and intervention 2,
#' # bootstrap with 1000 samples with normal quantile,
#' # natural course as the reference intervention
#'
#' ice_strat_haz <- ice(data = data, time_points = 4, id = "id", time_name = "t0",
#' censor_name = "C", outcome_name = "Y",
#' compevent_name = "D",
#' outcome_model = Y ~ L1 + L2, censor_model = C ~ L1 + L2,
#' competing_model = D ~ L1 + L2,
#' comp_effect = 1,
#' ref_idx = 0,
#' estimator = strat(hazard = T),
#' nsamples = 1000, ci_method = "normal",
#' int_descript = c("Static Intervention",
#' "Dynamic Intervention"),
#' intervention1.A1 = list(static(3)),
#' intervention1.A2 = list(static(1)),
#' intervention2.A1 = list(dynamic_cat),
#' intervention2.A2 = list(dynamic("compare", "L1", "=", 0))
#' )
#'
#' plot_risk(ice_strat_haz)
#'
#' # Hazard based stratified ICE,
#' # competing event as total effect for all interventions,
#' # with Y ~ L1 + L2 as outcome model for intervention 1,
#' # Y ~ L1 as outcome model for intervention 2,
#' # D ~ L1 as competing model for intervention 1,
#' # and D ~ L1 + L2 as competing model for intervention 2,
#' # bootstrap with 1000 samples with normal quantile,
#' # natural course as the reference intervention
#'
#' ice_strat_haz <- ice(data = data, time_points = 4, id = "id", time_name = "t0",
#' censor_name = "C", outcome_name = "Y",
#' compevent_name = "D",
#' outcome_model = Y ~ L1, censor_model = C ~ L1,
#' competing_model = D ~ L1,
#' comp_effect = 1,
#' ref_idx = 0,
#' estimator = strat(hazard = T),
#' nsamples = 1000, ci_method = "normal",
#' int_descript = c("Static Intervention",
#' "Dynamic Intervention"),
#' intervention1.A1 = list(static(3)),
#' intervention1.A2 = list(static(1)),
#' intervention2.A1 = list(dynamic_cat),
#' intervention2.A2 = list(dynamic("compare", "L1", "=", 0)),
#' outcomeModel.1 = Y ~ L1 + L2,
#' compModel.2 = D ~ L1 + L2
#' )
#'
#' plot_risk(ice_strat_haz)
#'
#' # Classical pooled ICE,
#' # competing event as direct effect for all interventions,
#' # bootstrap with 1000 samples using empirical quantile,
#' # natural course as the reference intervention.
#'
#' ice_pool_classic <- ice(data = data, time_points = 4, id = "id", time_name = "t0",
#' censor_name = "C", outcome_name = "Y",
#' compevent_name = "D",
#' comp_effect = 0,
#' outcome_model = Y ~ L1 + A1 + A2, censor_model = C ~ L1 + A1 + A2,
#' competing_model = D ~ L1 + A1 + A2,
#' ref_idx = 0,
#' estimator = pool(hazard = F),
#' nsamples = 1000, ci_method = "percentile",
#' int_descript = c("Static Intervention",
#' "Dynamic Intervention"),
#' intervention1.A1 = list(static(3)),
#' intervention1.A2 = list(static(1)),
#' intervention2.A1 = list(dynamic_cat),
#' intervention2.A2 = list(dynamic("compare", "L1", "=", 0))
#' )
#'
#' plot_risk(ice_pool_classic)
#'
#' # Hazard based pooled ICE,
#' # competing event as direct effect for all interventions,
#' # bootstrap with 1000 samples using empirical quantile,
#' # always treat as the reference intervention.
#'
#' ice_pool_haz <- ice(data = data, time_points = 4, id = "id", time_name = "t0",
#' censor_name = "C", outcome_name = "Y",
#' compevent_name = "D",
#' comp_effect = 0,
#' outcome_model = Y ~ L1 + A1 + A2, censor_model = C ~ L1 + A1 + A2,
#' competing_model = D ~ L1 + A1 + A2,
#' ref_idx = 1,
#' estimator = pool(hazard = T),
#' nsamples = 1000, ci_method = "percentile",
#' int_descript = c("Static Intervention",
#' "Dynamic Intervention"),
#' intervention1.A1 = list(static(3)),
#' intervention1.A2 = list(static(1)),
#' intervention2.A1 = list(dynamic_cat),
#' intervention2.A2 = list(dynamic("compare", "L1", "=", 0))
#' )
#'
#' plot_risk(ice_pool_haz)
#'
#' # Doubly robust ICE,
#' # competing event as direct effect for all interventions,
#' # with Y ~ L1 as outcome model for both intervention 1 and intervention 2,
#' # and D ~ L1 as competing model for both intervention 1 and intervention 2,
#' # bootstrap with 1000 samples using normal quantile,
#' # natural course as the reference intervention.
#'
#' ice_weight <- ice(data = data, time_points = 4, id = "id", time_name = "t0",
#' censor_name = "C", outcome_name = "Y",
#' compevent_name = "D",
#' comp_effect = 0,
#' outcome_model = Y ~ L1, censor_model = C ~ L1,
#' competing_model = D ~ L1,
#' ref_idx = 0,
#' estimator = weight(list(A1 ~ L1 + L2, A2 ~ L1 + L2)),
#' nsamples = 1000, ci_method = "normal",
#' int_descript = c("Static Intervention",
#' "Dynamic Intervention"),
#' intervention1.A1 = list(static(3)),
#' intervention1.A2 = list(static(1)),
#' intervention2.A1 = list(dynamic_cat),
#' intervention2.A2 = list(dynamic("compare", "L1", "=", 0))
#' )
#'
#' plot_risk(ice_weight)
#'
#' # Doubly robust ICE,
#' # competing event as direct effect for all interventions,
#' # with Y ~ L1 as outcome model for intervention 1,
#' # Y ~ L1 + L2 as outcome model for intervention 2,
#' # D ~ L1 + L2 as competing model for intervention 1,
#' # D ~ L1 as competing model for intervention 2,
#' # bootstrap with 1000 samples using normal quantile,
#' # natural course as the reference intervention.
#'
#' ice_weight <- ice(data = data, time_points = 4, id = "id", time_name = "t0",
#' censor_name = "C", outcome_name = "Y",
#' compevent_name = "D",
#' comp_effect = 0,
#' outcome_model = Y ~ L1, censor_model = C ~ L1,
#' competing_model = D ~ L1,
#' ref_idx = 0,
#' estimator = weight(list(A1 ~ L1 + L2, A2 ~ L1 + L2)),
#' nsamples = 1000, ci_method = "normal",
#' int_descript = c("Static Intervention",
#' "Dynamic Intervention"),
#' intervention1.A1 = list(static(3)),
#' intervention1.A2 = list(static(1)),
#' intervention2.A1 = list(dynamic_cat),
#' intervention2.A2 = list(dynamic("compare", "L1", "=", 0)),
#' outcomeModel.2 = Y ~ L1 + L2,
#' compModel.1 = D ~ L1 + L2
#' )
#'
#' plot_risk(ice_weight)
#'
#'
#' # The following example implements the natural value intervention on L1
#' # (i.e. if L1 > 0, then replace the observed value of L1 by 0,
#' # and keep the observed value of L1 otherwise),
#' # dynamic intervention on L1 (treat when L1 = 0)
#' # with uniform grace period of 2 periods, and
#' # threshold intervention when the natural value of L2 at time t is lower
#' # than -3, set its value to -3. Otherwise, do not intervene.
#'
#' # Classical pooled ICE,
#' # competing event as direct effect for all interventions,
#' # bootstrap with 1000 samples using empirical quantile,
#' # natural course as the reference intervention.
#'
#' ice_pool_grace_period <- ice(data = data, time_points = 4, id = "id", time_name = "t0",
#' censor_name = "C", outcome_name = "Y",
#' compevent_name = "D",
#' comp_effect = 0,
#' outcome_model = Y ~ L1 + L2 + A1 + A2, censor_model = C ~ L1 + L2 + A1 + A2,
#' competing_model = D ~ L1 + L2 + A1 + A2,
#' ref_idx = 0,
#' estimator = pool(hazard = F),
#' nsamples = 1000, ci_method = "percentile",
#' int_descript = c("Grace Period", "Threshold Intervention")
#' intervention1.A2 = list(grace_period("uniform", 2, "L1", 0)),
#' intervention2.L2 = list(threshold(-3))
#' )
#'
#' plot_risk(ice_pool_grace_period)
#'
#'

ice <- function(data, time_points, id, time_name,
                outcome_name, censor_name = NULL,
                compevent_name = NULL, comp_effect = 0,
                outcome_model, censor_model = NULL, competing_model = NULL,
                ref_idx = 0,
                estimator,
                int_descript,
                ci_method = "percentile",
                nsamples = 0, seed = 1,
                significance_level = 0.05, parallel = F, ncores = 2,
                ...) {

  ## pre-process user inputs

  set_seed <- seed
  K <- time_points
  competing_name <- compevent_name
  nboot <- nsamples

  if (nboot > 0) {
    bootstrap <- T
  } else {
    bootstrap <- F
  }

  if (ci_method == "normal") {
    normal_quantile <- T
  } else {
    normal_quantile <- F
  }

  interv_data <<- data
  idvar <<- id
  time0var <<- time_name
  outcomevar <<- outcome_name

  # preprocess the time index

  unique_time_names <- unique(data[, time_name])

  if (K > length(unique_time_names)) {
    stop("Please enter a number of total time points below the maximum time units in the data.")
  }

  map_time_column <- function(time_name, data, total_times) {
    data[, paste0("new_", time_name)] <- NA
    for (t in 0:(length(total_times)-1)) {
      map_idx <- which(data[, time_name] == total_times[t+1])
      data[map_idx, paste0("new_", time_name)] <- t
    }
    return(data)
  }

  if (class(data[, time_name]) != "integer") {
    data <- map_time_column(time_name, data, unique_time_names)
    time_name <- paste0("new_", time_name)
  } else {

    # if the time does not start with 0, then re-index
    if (min(unique_time_names) != 0) {
      data <- map_time_column(time_name, data, unique_time_names)
      time_name <- paste0("new_", time_name)
    }
  }
  
  ## preprocess threshold intervention

  sub_kwarg <- substitute(list(...))
  clean_kwarg_str <- paste0(str_remove_all(as.character(deparse(sub_kwarg)), " "), collapse = "")
  kwarg_list <- as.list(sub_kwarg)
  kwarg_name_list <- names(kwarg_list)
  threshold_idx <- str_which(str_split(as.character(substitute(list(...))), " = "), "threshold")

  if (length(threshold_idx) > 0) {
  threshold_kwarg <- kwarg_name_list[threshold_idx]
  threshold_interv_list_split <- str_split(threshold_kwarg, "[.]")
  threshold_interv_list <- lapply(threshold_interv_list_split, function(x) {x[1]})
  threshold_interv_list_all_names <- lapply(threshold_interv_list_split, function(x) {x[2]})

  for (i in 1:length(threshold_idx)) {
    ikwarg <- threshold_kwarg[i]
    ivar <- as.character(threshold_interv_list_all_names[i])
    int_name_length <- nchar(ikwarg) + 5
    origin_list <- as.character(kwarg_list[ikwarg])
    threshold_str_keep <- unlist(str_split(origin_list, "\\)"))[1]
    search_str <- paste0(ikwarg, "=", threshold_str_keep)
    search_str <- str_remove_all(search_str, " ")
    if (!str_detect(search_str, "var")) {
    start_idx <- as.vector(str_locate(clean_kwarg_str, ikwarg))[1]
    end_idx <- start_idx + nchar(search_str) - 1
    keep_str <- str_split(origin_list, "")[[1]][c(1:nchar(threshold_str_keep))]
    keep_str <- paste0(tail(keep_str, -5), collapse = "")
    add_str <- paste0(", var=\"", ivar, "\")")
    call_str <- paste0(keep_str, add_str)
    replace_str <- paste0(ikwarg, "=list(", call_str, ")")
    clean_arg_split <- str_split(clean_kwarg_str, "")[[1]]
    start_replace_idx <- start_idx + int_name_length
    end_replace_idx <- end_idx + 1
    clean_arg_split <- c(clean_arg_split[c(1:(start_replace_idx))], str_split(call_str, "")[[1]], clean_arg_split[c((end_replace_idx+1):length(clean_arg_split))])
    clean_kwarg_str <- paste0(clean_arg_split, collapse = "")
    }
  }

  }

  ## get argument list
  args <- c(as.list(environment()), eval(parse(text = clean_kwarg_str)))
  kwargs_len <- length(eval(parse(text = clean_kwarg_str)))
  if (kwargs_len == 0) {
    warning("No interventions are specified.")
  }

  ## get data info

  dta_len <- nrow(data)

  ## match internal args

  total_effect <- comp_effect

  ## get interventions

  arg_names <- names(args)
  arg_idx <- which(str_detect(arg_names, "intervention"))
  arg_interv <- args[arg_idx]

  ## get outcome models

  outcome_idx <- which(str_detect(arg_names, "outcomeModel"))
  outcome_interv <- args[outcome_idx]

  ## get competing models

  comp_idx <- which(str_detect(arg_names, "compModel"))
  comp_interv <- args[comp_idx]

  ## input checks

  if (kwargs_len != 0 & length(arg_idx) == 0) {
    stop("Please enter interventions through keywords arguments following the naming pattern:
    interventioni.treatment where treatment is the variable name in intervention i
    Example:
    intervention1.A1 (treatment 1 in intervention 1),
    intervention1.A2 (treatment 2 in intervention 1),
    intervention2.A1 (treatment 1 in intervention 2),
    intervention2.A2 (treatment 2 in intervention 2)")
  }

  if (!total_effect %in% 0:2) {

    stop("Please input 0, 1, or 2 for comp_effect.
         0 outputs direct effect for all interventions.
         1 outputs total effect for all interventions.
         2 outputs both direct effect and total effect for all interventions.")
  }

  if (!is.data.frame(data)) {

    stop("Please input data using a data frame.")
  }

  dta_cols <- colnames(data)

  if (!id %in% dta_cols) {
    stop("Please input the id column that is within the input data frame.")
  }

  if (!time_name %in% dta_cols) {
    stop("Please input the time column that is within the input data frame.")
  }

  if (!outcome_name %in% dta_cols) {
    stop("Please input the outcome column that is within the input data frame.")
  }

  if (!is.null(censor_name)) {
    if (!censor_name %in% dta_cols) {
    stop("Please input the censoring event column that is within the input data frame.")
    }
  }

  if (!is.null(competing_name)) {
    if (!competing_name %in% dta_cols) {
    stop("Please input the competing event column that is within the input data frame.")
    }
  }

  if (!is_formula(outcome_model)) {
    stop("Please input a formula for outcome_model.")
  }

  if (!is.null(censor_name)) {
    if (!is_formula(censor_model)) {
    stop("Please input a formula for censor_model.")
    }
  }

  if (any(str_detect(as.character(substitute(estimator)), "pool")) == F & any(str_detect(as.character(substitute(estimator)), "strat")) == F & any(str_detect(as.character(substitute(estimator)), "weight")) == F) {
    stop("Please input a valid estimator.
         pool(hazard = F) for classical pooling over treatment history ICE.
         pool(hazard = T) for hazard based pooling over treatment history ICE.
         strat(hazard = F) for classical stratifying on treatment history ICE.
         strat(hazard = T) for hazard based stratifying on treatment history ICE.
         weight() for weighted ICE and input the treatment models through the treat_model argument.")
  }

  ## get numbers of intervention

  if (kwargs_len == 0) {
    ninterv <- 0
    interventions <- list()
    intervention_varnames <- list()
    intervention_descriptions <- list()
  } else {
  interv_args_split <- split_args(arg_interv, "intervention")
  interv_list_origin <- interv_args_split$origin_list
  interv_list <- unlist(interv_args_split$prefix)
  interv_list_all_names <- interv_args_split$suffix
  ninterv <- length(unique(interv_list))

  # get outcome models
  outcomeModel_args_list <- split_args(outcome_interv, "outcomeModel")
  outcomeModel_interv <- outcomeModel_args_list$suffix

  # get competing models
  compModel_args_list <- split_args(comp_interv, "compModel")
  compModel_interv <- compModel_args_list$suffix

  # make outcome model list
  outcomeModels <- get_model_formula(interv_list, outcome_interv, outcomeModel_interv, ninterv, "outcomeModel")

  # make competing model list
  compModels <- get_model_formula(interv_list, comp_interv, compModel_interv, ninterv, "compModel")


  if (ref_idx > ninterv | ref_idx < 0) {
    stop("Please input a valid index for the reference intervention.
         0 represents natural course as reference.
         Other integer represents the corresponding specified intervention.")
  }

  ## convert to old intervention format

  interventions <- list()
  intervention_varnames <- list()
  intervention_descriptions <- list()
  intervention_times <- list()

  for (i in 1:ninterv) {

    # get all treatments for each intervention
    interv_idx <- which((unlist(interv_list)) == as.character(i))
    interv_all_names <- unlist(interv_list_all_names)[interv_idx]

    interv_single <- list()
    interv_var_single <- list()
    interv_time_single <- list()

    for (treat in 1:length(interv_idx)) {

      treat_idx <- interv_idx[treat]
      treat_name <- paste0("intervention", interv_list_origin[treat_idx])
      interv <- as.list(arg_interv[[which(names(arg_interv) == treat_name)]])[[1]]
      interv_var <- interv_all_names[treat]
      des <- as.list(int_descript[i])
      if (length(arg_interv[[which(names(arg_interv) == treat_name)]]) == 2) {
      times <- arg_interv[[which(names(arg_interv) == treat_name)]][[2]]
      } else {
        times <- 0:(K-1)
      }

      interv_single <- c(interv_single, list(interv))
      interv_var_single <- c(interv_var_single, list(interv_var))
      interv_time_single <- c(interv_time_single, list(times))

    }

    interventions <- c(interventions, list(interv_single))
    intervention_varnames <- c(intervention_varnames, list(interv_var_single))
    intervention_descriptions <- c(intervention_descriptions, list(des))
    intervention_times <- c(intervention_times, list(interv_time_single))

  }
  }

  obs_treatment_name <- unique(unlist(intervention_varnames))

  for (iobs in 1:length(obs_treatment_name)) {

  if (!obs_treatment_name[iobs] %in% dta_cols) {
    stop("Please input the observed treatment column that is within the input data frame.")
  }
  }

  if (total_effect == 0) {
    ref_total_effect_list <- list(F)
    total_effect_lists <- list(list(rep(F, ninterv)))

  } else if (total_effect == 1) {
    ref_total_effect_list <- list(T)
    total_effect_lists <- list(list(rep(T, ninterv)))

  } else if (total_effect == 2) {
    ref_total_effect_list <- list(T, F)
    total_effect_lists <- list(list(rep(T, ninterv)), list(rep(F, ninterv)))
  }

  if (kwargs_len == 0) {
    total_effect_lists <- ref_total_effect_list
  }

  if (any(str_detect(as.character(substitute(estimator)), "pool"))) {
    if (estimator) {
      ice_colname <- "HB Pooled ICE"
    } else {
      ice_colname <- "Classical Pooled ICE"
    }

  } else if (any(str_detect(as.character(substitute(estimator)), "strat"))) {
    if (estimator) {
      ice_colname <- "HB Stratified ICE"
    } else {
      ice_colname <- "Classical Stratified ICE"
    }

  } else if (any(str_detect(as.character(substitute(estimator)), "weight"))) {

      ice_colname <- "Weighted ICE"

  }

  risk_time_all <- data.frame()
  summary_all <- data.frame()

  for (ite in 1:length(ref_total_effect_list)) {

    ref_total_effect <- ref_total_effect_list[[ite]]
    total_effect <- as.list(total_effect_lists[[ite]][[1]])

  if (bootstrap) {
    summary <- as.data.frame(matrix(NA, nrow = (ninterv + 1), ncol = 13))
    rownames(summary) <- c('Natural Course', unlist(intervention_descriptions))
    colnames(summary) <- c('NP Risk', ice_colname, "Risk Ratio", "Risk Difference", "ICE Risk SE",
                           "RR SE", "RD SE", paste0("ICE Risk ", 100*(1-significance_level),"% Lower Bound"),
                           paste0("ICE Risk ", 100*(1-significance_level),"% Upper Bound"),
                           paste0("RR ", 100*(1-significance_level),"% Lower Bound"),
                           paste0("RR ", 100*(1-significance_level),"% Upper Bound"),
                           paste0("RD ", 100*(1-significance_level),"% Lower Bound"),
                           paste0("RD ", 100*(1-significance_level),"% Upper Bound"))

    risk_time <- as.data.frame(matrix(NA, ncol = 6, nrow = (K+1)*(ninterv+2)))
    colnames(risk_time) <- c("Intervention", "Risk", "Time", "SE", "Critical_Value_Upper", "Critical_Value_Lower")
    risk_time$Time <- rep(0:K, (ninterv+2))
  } else {
    summary <- as.data.frame(matrix(NA, nrow = (ninterv + 1), ncol = 4))
    rownames(summary) <- c('Natural Course', unlist(intervention_descriptions))
    colnames(summary) <- c('NP Risk', ice_colname, "Risk Ratio", "Risk Difference")

    risk_time <- as.data.frame(matrix(NA, ncol = 3, nrow = (K+1)*(ninterv+2)))
    colnames(risk_time) <- c("Intervention", "Risk", "Time")
    risk_time$Time <- rep(0:K, (ninterv+2))
  }


  if (any(str_detect(as.character(substitute(estimator)), "pool"))) {

    hazard <- estimator

    # check model input for pooled ICE


    if (!any(str_detect(as.character(outcome_model)[3], unlist(obs_treatment_name)))) {
      stop("Please include treatment variable as covariate in the outcome model for pooled ICE.")
    }

    if (!is.null(censor_name)) {
      if (!any(str_detect(as.character(censor_model)[3], unlist(obs_treatment_name)))) {
        stop("Please include treatment variable as covariate in the censoring event model for pooled ICE.")
      }
    }

    if (!is.null(competing_name)) {

      if (hazard) {
        if (!is_formula(competing_model)) {
          stop("Please input the competing model statement as a formula in competing_model for hazard based ICE.")
        }

        if (!any(str_detect(as.character(competing_model)[3], unlist(obs_treatment_name)))) {
          stop("Please include treatment variable as covariate in the competing model for pooled ICE.")
        }
      }
    }


    risk_descript <- c()
    risk_interv <- c()
    fit_all <- fit_summary <- fit_stderr <- fit_vcov <- fit_rmse <- c()

    nc_intervention_varlist <- list(obs_treatment_name)
    nc_descript <- "Natural Course"

    nc_interventions <- list()
    for (i_nc in 1:length(nc_intervention_varlist[[1]])) {
      i_nc_treat <- nc_intervention_varlist[[1]][[i_nc]]
      nc_interventions <- c(nc_interventions, list(natural_course(data = data, treat_var = i_nc_treat)))
    }

    nc_interventions <- list(nc_interventions)
    ref_int_times <- list()

    ref <- ice_pool(data = data, K = K, id = id, time_name = time_name, outcome_name = outcome_name,
                    censor_name = censor_name, competing_name = competing_name, total_effect = ref_total_effect,
                    outcome_model = outcome_model, censor_model = censor_model, competing_model = competing_model,
                    interventions = nc_interventions, intervention_names = nc_intervention_varlist,
                    compute_nc_risk = T,
                    hazard_based = hazard, intervention_description = nc_descript)

    summary[1, 1:2] <- c(ref$weight_h$risk[K], ref$gformula_risk_last_time)

    if (ref_idx == 0) {

      ref_intervention_varlist <- nc_intervention_varlist
      ref_description = nc_descript

      if (bootstrap) {
        summary[1, c(3:4, 6:7, 10:13)] <- c(1, 0, 1, 0, 1, 1, 0, 0)
      } else {
        summary[1, 3:4] <- c(1, 0)
      }

      boot_interv <- nc_interventions

    } else {

      ref_intervention_varlist <- list(intervention_varnames[[ref_idx]])
      ref_description <- intervention_descriptions[[ref_idx]]
      boot_interv <- list(interventions[[ref_idx]])
      ref_int_times <- list(intervention_times[[ref_idx]])

      if (!is.character(unlist(ref_intervention_varlist))) {
        stop("Please input the treatment variable at the second element of the input list for each intervention argument.")
      }

      if (is.character(unlist(ref_intervention_varlist))) {
        if (any(unlist(ref_intervention_varlist) %in% dta_cols) == F) {
          stop("The input treatment variable must be one of the columns in the data.")
        }
      }

      ref <- ice_pool(data = data, K = K, id = id, time_name = time_name, outcome_name = outcome_name,
                    censor_name = censor_name, competing_name = competing_name, total_effect = ref_total_effect,
                    outcome_model = outcome_model, censor_model = censor_model, competing_model = competing_model,
                    interventions = boot_interv, intervention_names = ref_intervention_varlist,
                    intervention_times = ref_int_times,
                    compute_nc_risk = T, hazard_based = hazard, intervention_description = ref_description)

      if (bootstrap) {
        summary[ref_idx + 1, c(2:4, 6:7, 10:13)] <- c(ref$gformula_risk_last_time, 1, 0, 1, 0, 1, 1, 0, 0)
      } else {
        summary[ref_idx + 1, c(2:4)] <- c(ref$gformula_risk_last_time, 1, 0)
      }
    }

    this_model <- list(ref$fit_models)
    names(this_model) <- ref_description

    this_summary <- list(ref$model_summary)
    names(this_summary) <- ref_description

    this_stderr <- list(ref$model_stderr)
    names(this_stderr) <- ref_description

    this_vcov <- list(ref$model_vcov)
    names(this_vcov) <- ref_description

    this_rmse <- list(ref$model_rmse)
    names(this_rmse) <- ref_description

    fit_all <- c(fit_all, this_model)
    fit_summary <- c(fit_summary, this_summary)
    fit_stderr <- c(fit_stderr, this_stderr)
    fit_vcov <- c(fit_vcov, this_vcov)
    fit_rmse <- c(fit_rmse, this_rmse)

    risk_descript <- c(risk_descript, c(rep(str_to_title(ref_description), K+1), rep("Natural Course NP", K+1)))
    risk_interv <- c(risk_interv, c(ref$gformula_risk[1, ], c(0, ref$weight_h$risk)))

    if (bootstrap) {
      intervention_descriptions <- c(intervention_descriptions, list("Natural Course"))
    }

    ninterv_new <- length(intervention_descriptions)

    ice_critical_value <- rr_critical_value <- rd_critical_value  <- matrix(NA, nrow = 2, ncol = ninterv + 1)

    if (bootstrap) {
      fit_all_boot <- fit_summary_boot <- fit_stderr_boot <- fit_vcov_boot <- fit_rmse_boot <- data_boot_all <- c()
    } else {
      fit_all_boot <- fit_summary_boot <- fit_stderr_boot <- fit_vcov_boot <- fit_rmse_boot <- data_boot_all <- NULL
    }

    if (ninterv_new > 0) {
      
      critical_value_all_lower <- critical_value_all_upper <- se_all <- c()

    for (int in 1:ninterv_new) {
      
      this_descript <- intervention_descriptions[[int]]
      
      if (this_descript == "Natural Course") {
        this_int_var <- list(obs_treatment_name)
        this_interv <- list()
        for (i_nc in 1:length(this_int_var[[1]])) {
          i_nc_treat <- this_int_var[[1]][[i_nc]]
          this_interv <- c(this_interv, list(natural_course(data = data, treat_var = i_nc_treat)))
        }
        
        this_interv <- list(this_interv)
        
        this_time <- NULL
        this_outcome_formula <- outcome_model
        this_comp_formula <- competing_model
        this_total_effect <- ref_total_effect
        
      } else {
        this_int_var <- list(intervention_varnames[[int]])
        this_interv <- list(interventions[[int]])
        this_total_effect <- total_effect[[int]]
        this_time <- list(intervention_times[[int]])
      }
      


      if (!is.character(unlist(this_int_var))) {
        stop("Please input the treatment variable at the second element of the input list for each intervention argument.")
      }

      if (is.character(unlist(this_int_var))) {
        if (any(unlist(this_int_var) %in% dta_cols) == F) {
          stop("The input treatment variable must be one of the columns in the data.")
        }
      }

      this_fit <- ice_pool(data = data, K = K, id = id, time_name = time_name, outcome_name = outcome_name,
                            censor_name = censor_name, competing_name = competing_name, total_effect = this_total_effect,
                           outcome_model = outcome_model, censor_model = censor_model, competing_model = competing_model,
                            interventions = this_interv, intervention_names = this_int_var, intervention_times = this_time,
                            hazard_based = hazard, intervention_description = this_descript )

      this_model <- list(this_fit$fit_models)
      names(this_model) <- this_descript

      this_summary <- list(this_fit$model_summary)
      names(this_summary) <- this_descript

      this_stderr <- list(this_fit$model_stderr)
      names(this_stderr) <- this_descript

      this_vcov <- list(this_fit$model_vcov)
      names(this_vcov) <- this_descript

      this_rmse <- list(this_fit$model_rmse)
      names(this_rmse) <- this_descript

      fit_all <- c(fit_all, this_model)
      fit_summary <- c(fit_summary, this_summary)
      fit_stderr <- c(fit_stderr, this_stderr)
      fit_vcov <- c(fit_vcov, this_vcov)
      fit_rmse <- c(fit_rmse, this_rmse)


      if (this_descript != "Natural Course") {
        summary[int+1, 2] <- this_fit$gformula_risk_last_time
        risk_descript <- c(risk_descript, rep(str_to_title(this_descript), K+1))
        risk_interv <- c(risk_interv, this_fit$gformula_risk[1, ])

      }


    if (bootstrap) {


      this_boot <- bootstrap_ice(ice_pool, K, nboot, significance_level, parallel, ncores, ref_description,
                                 ref_intervention_varlist[[1]], ref_total_effect, this_total_effect,
                                 "ipw", boot_interv,
                                 this_interv, this_int_var, this_descript, data, id, set_seed,
                                 time_name = time_name, outcome_name = outcome_name,
                                 censor_name = censor_name, competing_name = competing_name,
                                 hazard_based = hazard, outcome_model = outcome_model,
                                 censor_model = censor_model, competing_model = competing_model,
                                 intervention_times = this_time)

      this_se <- this_boot$ice_se[K+1]
      this_rr_se <- this_boot$rr_se
      this_rd_se <- this_boot$rd_se

      this_cv_upper <- this_boot$ice_cv_upper
      this_cv_lower <- this_boot$ice_cv_lower
      this_rr_cv_upper <- this_boot$rr_cv_upper
      this_rr_cv_lower <- this_boot$rr_cv_lower
      this_rd_cv_upper <- this_boot$rd_cv_upper
      this_rd_cv_lower <- this_boot$rd_cv_lower

      this_model_boot <- list(this_boot$boot_models)
      names(this_model_boot) <- this_descript

      this_summary_boot <- list(this_boot$boot_summary)
      names(this_summary_boot) <- this_descript

      this_stderr_boot <- list(this_boot$boot_stderr)
      names(this_stderr_boot) <- this_descript

      this_vcov_boot <- list(this_boot$boot_vcov)
      names(this_vcov_boot) <- this_descript

      this_rmse_boot <- list(this_boot$boot_rmse)
      names(this_rmse_boot) <- this_descript

      this_data_boot <- list(this_boot$boot_data)
      names(this_data_boot) <- this_descript

      fit_all_boot <- c(fit_all_boot, this_model_boot)
      fit_summary_boot <- c(fit_summary_boot, this_summary_boot)
      fit_stderr_boot <- c(fit_stderr_boot, this_stderr_boot)
      fit_vcov_boot <- c(fit_vcov_boot, this_vcov_boot)
      fit_rmse_boot <- c(fit_rmse_boot, this_rmse_boot)
      data_boot_all <- c(data_boot_all, this_data_boot)
      
      critical_value_all_upper <- append_list(this_boot, str_to_title(this_descript), critical_value_all_upper, "ice_cv_all_upper")
      critical_value_all_lower <- append_list(this_boot, str_to_title(this_descript), critical_value_all_lower, "ice_cv_all_lower")
      
      se_all <- append_list(this_boot, str_to_title(this_descript), se_all, "ice_se")


      if (this_descript != "Natural Course") {
        summary[int+1, 5:7] <- c(this_se, this_rr_se, this_rd_se)
      } else {
        summary[1, 5:7] <- c(this_se, this_rr_se, this_rd_se)
      }
      
      
    
      if (ref_idx == 0) {
        
        if (this_descript != "Natural Course") {
          ice_critical_value[1, int+1] <- this_cv_lower
          rr_critical_value[1, int+1] <- this_rr_cv_lower
          rd_critical_value[1, int+1] <- this_rd_cv_lower
          ice_critical_value[2, int+1] <- this_cv_upper
          rr_critical_value[2, int+1] <- this_rr_cv_upper
          rd_critical_value[2, int+1] <- this_rd_cv_upper
        } else {
          critical_value_all_upper <- append_list(this_boot, "Natural Course NP", critical_value_all_upper, "ref_ipw_cv_all_upper")
          critical_value_all_lower <- append_list(this_boot, "Natural Course NP", critical_value_all_lower, "ref_ipw_cv_all_lower")
          se_all <- append_list(this_boot, "Natural Course NP", se_all, "ref_ipw_se")
        }
      } else {
        if (this_descript != "Natural Course") {
          ice_critical_value[1, int+1] <- this_cv_lower
          rr_critical_value[1, int+1] <- this_rr_cv_lower
          rd_critical_value[1, int+1] <- this_rd_cv_lower
          ice_critical_value[2, int+1] <- this_cv_upper
          rr_critical_value[2, int+1] <- this_rr_cv_upper
          rd_critical_value[2, int+1] <- this_rd_cv_upper
        } else {
          ice_critical_value[1, 1] <- this_cv_lower
          rr_critical_value[1, 1] <- this_rr_cv_lower
          rd_critical_value[1, 1] <- this_rd_cv_lower
          ice_critical_value[2, 1] <- this_cv_upper
          rr_critical_value[2, 1] <- this_rr_cv_upper
          rd_critical_value[2, 1] <- this_rd_cv_upper
          critical_value_all_upper <- append_list(this_boot, "Natural Course NP", critical_value_all_upper, "ref_ipw_cv_all_upper")
          critical_value_all_lower <- append_list(this_boot, "Natural Course NP", critical_value_all_lower, "ref_ipw_cv_all_lower")
          se_all <- append_list(this_boot, "Natural Course NP", se_all, "ref_ipw_se")
        }
      }

    }
      if (bootstrap) {
        summary[ref_idx + 1, 5] <- this_boot$ref_se[K+1]

        if (ref_idx == 0) {
          ice_critical_value[1, 1] <- this_boot$ref_cv_lower
          ice_critical_value[2, 1] <- this_boot$ref_cv_upper
        }
      }

    }
    }

    risk_time$Intervention <- risk_descript
    risk_time$Risk <- risk_interv
    
    if (bootstrap) {
      
      risk_time <- match_boot_values(risk_time, data.frame(se_all), "SE")
      
      if (normal_quantile == F) {
      risk_time <- match_boot_values(risk_time, data.frame(critical_value_all_upper), "Critical_Value_Upper")
      risk_time <- match_boot_values(risk_time, data.frame(critical_value_all_lower), "Critical_Value_Lower")
      }
    }




  } else if (any(str_detect(as.character(substitute(estimator)), "strat")) | any(str_detect(as.character(substitute(estimator)), "weight"))) {

    risk_descript <- c()
    risk_interv <- c()
    fit_all <- fit_summary <- fit_stderr <- fit_vcov <- fit_rmse <- c()

    if (any(str_detect(as.character(substitute(estimator)), "strat"))) {
      weight <- F
      this_treat_model <- list()
      this_obs_treatment_varnames <- as.list(unique(unlist(intervention_varnames)))
      hazard <- estimator
    } else {
      weight <- T
      this_treat_model <- estimator$treat_model
      obs_treatment_name <- estimator$obs_treat
      this_obs_treatment_varnames <- as.list(obs_treatment_name)
      hazard <- F

    }


    # check model input for strat ICE
    if (!is.null(competing_name)) {
      if (hazard) {
        if (!is_formula(competing_model)) {
          stop("Please input a formula of competing model statement in competing_model for hazard based ICE.")
        }
      }
    }

    if (weight) {
      if (any(unlist(lapply(this_treat_model, function(i) is_formula(i)))) == F | length(this_treat_model) != length(this_obs_treatment_varnames)) {
      stop("Please input a formula to specify the treatment model statement for each treatment variable.")
      }
    }

    nc_intervention_varlist <- list(obs_treatment_name)
    nc_descript <- "Natural Course"


    nc_interventions <- list()
    for (i_nc in 1:length(nc_intervention_varlist[[1]])) {
      i_nc_treat <- nc_intervention_varlist[[1]][[i_nc]]
      nc_interventions <- c(nc_interventions, list(natural_course(data = data, treat_var = i_nc_treat)))
    }

    nc_interventions <- list(nc_interventions)
    
    if (any(!is.na(unlist(compModels)))) {
      if (any(str_detect(as.character(substitute(estimator)), "strat")) & (hazard == F) & (ref_total_effect == F)) {
      warning("The competing model is used for nonparametric risk estimation for direct effect case in stratified ICE. The keyword argument competing model statments are ignored.")
      } 
    }

    ref <- ice_strat(data = data, K = K, id = id, time_name = time_name, outcome_name = outcome_name,
                     censor_name = censor_name, competing_name = competing_name, total_effect = ref_total_effect,
                     outcome_model = outcome_model, censor_model = censor_model, competing_model = competing_model,
                     interventions = nc_interventions, intervention_names = nc_intervention_varlist,
                     compute_nc_risk = T, hazard_based = hazard, weighted = weight,
                     intervention_description = nc_descript,
                     treat_model = this_treat_model,
                     obs_treatment_names = this_obs_treatment_varnames)

    summary[1, c(1:2)] <- c(ref$weight_h$risk[K], ref$gformula_risk_last_time)

    if (ref_idx == 0) {

      ref_intervention_varlist <- nc_intervention_varlist
      ref_description <- nc_descript

      if (bootstrap) {
        summary[1, c(3:4, 6:7, 10:13)] <- c(1, 0, 1, 0, 1, 1, 0, 0)
      } else {
        summary[1, c(3:4)] <- c(1, 0)
      }

      boot_interv <- nc_interventions
    } else {

      ref_intervention_varlist <- list(intervention_varnames[[ref_idx]])
      ref_description <- intervention_descriptions[[ref_idx]]
      boot_interv <- list(interventions[[ref_idx]])
      ref_int_times <- list(intervention_times[[ref_idx]])

      if (!is.character(unlist(ref_intervention_varlist))) {
        stop("Please input the treatment variable at the second element of the input list for each intervention argument.")
      }

      if (is.character(unlist(ref_intervention_varlist))) {
        if (any(unlist(ref_intervention_varlist) %in% dta_cols) == F) {
          stop("The input treatment variable must be one of the columns in the data.")
        }
      }
      
      if (any(!is.na(unlist(compModels)))) {
        if (any(str_detect(as.character(substitute(estimator)), "strat")) & (hazard == F) & (ref_total_effect == F)) {
          warning("The competing model is used for nonparametric risk estimation for direct effect case in stratified ICE. The keyword argument competing model statments are ignored.")
        } 
      }

      ref <- ice_strat(data = data, K = K, id = id, time_name = time_name, outcome_name = outcome_name,
                      censor_name = censor_name, competing_name = competing_name, total_effect = ref_total_effect,
                      outcome_model = outcome_model, censor_model = censor_model, competing_model = competing_model,
                      interventions = boot_interv, intervention_names = ref_intervention_varlist,
                      compute_nc_risk = T, hazard_based = hazard, weighted = weight,
                      intervention_description = ref_description, intervention_times = ref_int_times,
                      treat_model = this_treat_model,
                      obs_treatment_names = this_obs_treatment_varnames)

      if (bootstrap) {
        summary[ref_idx + 1, c(2:4, 6:7, 10:13)] <- c(ref$gformula_risk_last_time, 1, 0, 1, 0, 1, 1, 0, 0)
      } else {
        summary[ref_idx + 1, c(2:4)] <- c(ref$gformula_risk_last_time, 1, 0)
      }
    }

    this_model <- list(ref$fit_models)
    names(this_model) <- ref_description

    this_summary <- list(ref$model_summary)
    names(this_summary) <- ref_description

    this_stderr <- list(ref$model_stderr)
    names(this_stderr) <- ref_description

    this_vcov <- list(ref$model_vcov)
    names(this_vcov) <- ref_description

    this_rmse <- list(ref$model_rmse)
    names(this_rmse) <- ref_description

    fit_all <- c(fit_all, this_model)
    fit_summary <- c(fit_summary, this_summary)
    fit_stderr <- c(fit_stderr, this_stderr)
    fit_vcov <- c(fit_vcov, this_vcov)
    fit_rmse <- c(fit_rmse, this_rmse)

    risk_descript <- c(risk_descript, c(rep(str_to_title(ref_description), K+1), rep("Natural Course NP", K+1)))
    risk_interv <- c(risk_interv, c(ref$gformula_risk[1, ], c(0, ref$weight_h$risk)))

    if (bootstrap) {
      intervention_descriptions <- c(intervention_descriptions, list("Natural Course"))
    }

    ninterv_new <- length(intervention_descriptions)

    ice_critical_value <- rr_critical_value <- rd_critical_value  <- matrix(NA, nrow = 2, ncol = ninterv + 1)

    if (bootstrap) {
      fit_all_boot <- fit_summary_boot <- fit_stderr_boot <- fit_vcov_boot <- fit_rmse_boot <- data_boot_all <- c()
    } else {
      fit_all_boot <- fit_summary_boot <- fit_stderr_boot <- fit_vcov_boot <- fit_rmse_boot <- data_boot_all <- NULL
    }

    if (ninterv_new > 0) {
      
      critical_value_all_upper <- critical_value_all_lower <- se_all <- c()

    for (int in 1:ninterv_new) {
      
      this_descript <- intervention_descriptions[[int]]
      
      if (this_descript == "Natural Course") {
        this_int_var <- list(obs_treatment_name)
        this_interv <- list()
        for (i_nc in 1:length(this_int_var[[1]])) {
          i_nc_treat <- this_int_var[[1]][[i_nc]]
          this_interv <- c(this_interv, list(natural_course(data = data, treat_var = i_nc_treat)))
        }
        
        this_interv <- list(this_interv)
        
        this_time <- NULL
        this_outcome_formula <- outcome_model
        this_comp_formula <- competing_model
        this_total_effect <- ref_total_effect
        
      } else {
        this_int_var <- list(intervention_varnames[[int]])
        this_interv <- list(interventions[[int]])
        this_total_effect <- total_effect[[int]]
        this_time <- list(intervention_times[[int]])
        this_interv <- construct_interv_value(data, time_name, this_interv[[1]], this_time[[1]], this_int_var[[1]])
        this_interv <- list(this_interv)
        
        this_outcome_formula <- outcomeModels[[int]]
        
        if (any(is.na(as.character(this_outcome_formula)))) {
          this_outcome_formula <- outcome_model
        }
        
        this_comp_formula <- compModels[[int]]
        
        if (any(is.na(as.character(this_comp_formula)))) {
          this_comp_formula <- competing_model
        }
      }


      if (!is.character(unlist(this_int_var))) {
        stop("Please input the treatment variable at the second element of the input list for each intervention argument.")
      }

      if (is.character(unlist(this_int_var))) {
        if (any(unlist(this_int_var) %in% dta_cols) == F) {
          stop("The input treatment variable must be one of the columns in the data.")
        }
      }

      this_fit <- ice_strat(data = data, K = K, id = id, time_name = time_name, outcome_name = outcome_name,
                            censor_name = censor_name, competing_name = competing_name, total_effect = this_total_effect,
                            outcome_model = this_outcome_formula, censor_model = censor_model, competing_model = this_comp_formula,
                            interventions = this_interv, intervention_names = this_int_var,
                            hazard_based = hazard, weighted = weight, treat_model = this_treat_model,
                            obs_treatment_names = this_obs_treatment_varnames, intervention_description = this_descript,
                            intervention_times = this_time)

      this_model <- list(this_fit$fit_models)
      names(this_model) <- this_descript

      this_summary <- list(this_fit$model_summary)
      names(this_summary) <- this_descript

      this_stderr <- list(this_fit$model_stderr)
      names(this_stderr) <- this_descript

      this_vcov <- list(this_fit$model_vcov)
      names(this_vcov) <- this_descript

      this_rmse <- list(this_fit$model_rmse)
      names(this_rmse) <- this_descript

      fit_all <- c(fit_all, this_model)
      fit_summary <- c(fit_summary, this_summary)
      fit_stderr <- c(fit_stderr, this_stderr)
      fit_vcov <- c(fit_vcov, this_vcov)
      fit_rmse <- c(fit_rmse, this_rmse)

      if (this_descript != "Natural Course") {

      summary[int+1, 2] <- this_fit$gformula_risk_last_time
      risk_descript <- c(risk_descript, rep(str_to_title(this_descript), K+1))
      risk_interv <- c(risk_interv, this_fit$gformula_risk[1, ])
      }

      if (bootstrap) {

        this_boot <- bootstrap_ice(ice_strat, K, nboot, significance_level, parallel, ncores, ref_description,
                                   ref_intervention_varlist[[1]], ref_total_effect, this_total_effect,
                                   "ipw", boot_interv,
                                   this_interv, this_int_var, this_descript, data, id, set_seed,
                                   time_name = time_name, outcome_name = outcome_name,
                                   censor_name = censor_name, competing_name = competing_name,
                                   hazard_based = hazard, weighted = weight, treat_model = this_treat_model,
                                   obs_treatment_names = this_obs_treatment_varnames,
                                   outcome_model = outcome_model, censor_model = censor_model,
                                   competing_model = competing_model,
                                   intervention_times = this_time)

        this_se <- this_boot$ice_se[K+1]
        this_rr_se <- this_boot$rr_se
        this_rd_se <- this_boot$rd_se

        this_cv_upper <- this_boot$ice_cv_upper
        this_cv_lower <- this_boot$ice_cv_lower
        this_rr_cv_upper <- this_boot$rr_cv_upper
        this_rr_cv_lower <- this_boot$rr_cv_lower
        this_rd_cv_upper <- this_boot$rd_cv_upper
        this_rd_cv_lower <- this_boot$rd_cv_lower

        this_model_boot <- list(this_boot$boot_models)
        names(this_model_boot) <- this_descript

        this_summary_boot <- list(this_boot$boot_summary)
        names(this_summary_boot) <- this_descript

        this_stderr_boot <- list(this_boot$boot_stderr)
        names(this_stderr_boot) <- this_descript

        this_vcov_boot <- list(this_boot$boot_vcov)
        names(this_vcov_boot) <- this_descript

        this_rmse_boot <- list(this_boot$boot_rmse)
        names(this_rmse_boot) <- this_descript

        this_data_boot <- list(this_boot$boot_data)
        names(this_data_boot) <- this_descript

        fit_all_boot <- c(fit_all_boot, this_model_boot)
        fit_summary_boot <- c(fit_summary_boot, this_summary_boot)
        fit_stderr_boot <- c(fit_stderr_boot, this_stderr_boot)
        fit_vcov_boot <- c(fit_vcov_boot, this_vcov_boot)
        fit_rmse_boot <- c(fit_rmse_boot, this_rmse_boot)
        data_boot_all <- c(data_boot_all, this_data_boot)
        
        critical_value_all_upper <- append_list(this_boot, str_to_title(this_descript), critical_value_all_upper, "ice_cv_all_upper")
        critical_value_all_lower <- append_list(this_boot, str_to_title(this_descript), critical_value_all_lower, "ice_cv_all_lower")
        
        se_all <- append_list(this_boot, str_to_title(this_descript), se_all, "ice_se")

        if (this_descript != "Natural Course") {
          summary[int+1, 5:7] <- c(this_se, this_rr_se, this_rd_se)
        } else {
          summary[1, 5:7] <- c(this_se, this_rr_se, this_rd_se)
        }

        if (ref_idx == 0) {
          
          if (this_descript != "Natural Course") {
          ice_critical_value[1, int+1] <- this_cv_lower
          rr_critical_value[1, int+1] <- this_rr_cv_lower
          rd_critical_value[1, int+1] <- this_rd_cv_lower
          ice_critical_value[2, int+1] <- this_cv_upper
          rr_critical_value[2, int+1] <- this_rr_cv_upper
          rd_critical_value[2, int+1] <- this_rd_cv_upper
          } else {
            critical_value_all_upper <- append_list(this_boot, "Natural Course NP", critical_value_all_upper, "ref_ipw_cv_all_upper")
            critical_value_all_lower <- append_list(this_boot, "Natural Course NP", critical_value_all_lower, "ref_ipw_cv_all_lower")
            se_all <- append_list(this_boot, "Natural Course NP", se_all, "ref_ipw_se")
          }
        } else {
          if (this_descript != "Natural Course") {
            ice_critical_value[1, int+1] <- this_cv_lower
            rr_critical_value[1, int+1] <- this_rr_cv_lower
            rd_critical_value[1, int+1] <- this_rd_cv_lower
            ice_critical_value[2, int+1] <- this_cv_upper
            rr_critical_value[2, int+1] <- this_rr_cv_upper
            rd_critical_value[2, int+1] <- this_rd_cv_upper
          } else {
            ice_critical_value[1, 1] <- this_cv_lower
            rr_critical_value[1, 1] <- this_rr_cv_lower
            rd_critical_value[1, 1] <- this_rd_cv_lower
            ice_critical_value[2, 1] <- this_cv_upper
            rr_critical_value[2, 1] <- this_rr_cv_upper
            rd_critical_value[2, 1] <- this_rd_cv_upper
            critical_value_all_upper <- append_list(this_boot, "Natural Course NP", critical_value_all_upper, "ref_ipw_cv_all_upper")
            critical_value_all_lower <- append_list(this_boot, "Natural Course NP", critical_value_all_lower, "ref_ipw_cv_all_lower")
            se_all <- append_list(this_boot, "Natural Course NP", se_all, "ref_ipw_se")
          }
        }
      }
      if (bootstrap) {
        summary[ref_idx + 1, 5] <- this_boot$ref_se[K+1]

        if (ref_idx == 0) {
          ice_critical_value[1, 1] <- this_boot$ref_cv_lower
          ice_critical_value[2, 1] <- this_boot$ref_cv_upper
        }
      }

    }
    }

    risk_time$Intervention <- risk_descript
    risk_time$Risk <- risk_interv
    
    if (bootstrap) {
      
      risk_time <- match_boot_values(risk_time, data.frame(se_all), "SE")
      
      if (normal_quantile == F) {
        risk_time <- match_boot_values(risk_time, data.frame(critical_value_all_upper), "Critical_Value_Upper")
        risk_time <- match_boot_values(risk_time, data.frame(critical_value_all_lower), "Critical_Value_Lower")
      }
    }


  }

  ## ref RR and RD
  summary[ref_idx + 1, 3:4] <- c(1, 0)
  ## RR
  summary[-(ref_idx + 1), 3] <- summary[-(ref_idx + 1), 2] / summary[ref_idx + 1, 2]
  ## RD
  summary[-(ref_idx + 1), 4] <- summary[-(ref_idx + 1), 2] - summary[ref_idx + 1, 2]

  if (bootstrap) {

    if (normal_quantile) {
      ice_critical_value_lower <- rr_critical_value_lower <- rd_critical_value_lower <- ice_critical_value_upper <- rr_critical_value_upper <- rd_critical_value_upper <- qnorm(1 - significance_level/2)
    } else {
      
      ice_critical_value_lower <- as.vector(ice_critical_value[1, ])
      ice_critical_value_upper <- as.vector(ice_critical_value[2, ])
      if (ref_idx != 0) {
        rr_critical_value_lower <- as.vector(rr_critical_value[1, ])[-ref_idx]
        rd_critical_value_lower <- as.vector(rd_critical_value[1, ])[-ref_idx]
        rr_critical_value_upper <- as.vector(rr_critical_value[2, ])[-ref_idx]
        rd_critical_value_upper <- as.vector(rd_critical_value[2, ])[-ref_idx]
      } else {
        rr_critical_value_lower <- as.vector(rr_critical_value[1, ])[-1]
        rd_critical_value_lower <- as.vector(rd_critical_value[1, ])[-1]
        rr_critical_value_upper <- as.vector(rr_critical_value[2, ])[-1]
        rd_critical_value_upper <- as.vector(rd_critical_value[2, ])[-1]
      }
    }

    summary[1:(ninterv + 1), 8] <- summary[1:(ninterv + 1), 2] - ice_critical_value_lower * summary[1:(ninterv + 1), 5]
    summary[1:(ninterv + 1), 9] <- summary[1:(ninterv + 1), 2] + ice_critical_value_upper * summary[1:(ninterv + 1), 5]
    summary[-(ref_idx + 1), 10] <- summary[-(ref_idx + 1), 2] - rr_critical_value_lower * summary[-(ref_idx + 1), 6]
    summary[-(ref_idx + 1), 11] <- summary[2:(ninterv + 1), 2] + rr_critical_value_upper * summary[-(ref_idx + 1), 6]
    summary[-(ref_idx + 1), 12] <- summary[2:(ninterv + 1), 2] - rd_critical_value_lower * summary[-(ref_idx + 1), 7]
    summary[-(ref_idx + 1), 13] <- summary[2:(ninterv + 1), 2] + rd_critical_value_upper * summary[-(ref_idx + 1), 7]
    
    if (normal_quantile) {
    risk_time$Critical_Value_Lower <- ice_critical_value_lower
    risk_time$Critical_Value_Upper <- ice_critical_value_upper
    }

  }

  if (any(unlist(total_effect))) {
    rownames(summary) <- paste0(rownames(summary), " (", "Total Effect" , ")")
    risk_time$Intervention <- paste0(risk_time$Intervention, " (", "Total Effect" , ")")
  } else {
    rownames(summary) <- paste0(rownames(summary), " (", "Direct Effect" , ")")
    risk_time$Intervention <- paste0(risk_time$Intervention, " (", "Direct Effect" , ")")
  }
  summary_all <- rbind(summary_all, summary)
  risk_time_all <- rbind(risk_time_all, risk_time)
  }

  summary_all <- round(summary_all, 5)

  ## match back the time points using original specification
  risk_time_all$originTime <- c(0, unique_time_names)[as.integer(risk_time_all$Time) + 1]

  risk_time_all$originTime[which(is.na(risk_time_all$originTime))] <- 0


  return(list(summary = summary_all, risk.over.time = risk_time_all,
              models = fit_all, model.summary = fit_summary,
              model.stderr = fit_stderr, model.vcov = fit_vcov,
              model.rmse = fit_rmse, boot.data = data_boot_all,
              boot.models = fit_all_boot, boot.summary = fit_summary_boot,
              boot.stderr = fit_stderr_boot, boot.vcov = fit_vcov_boot,
              boot.rmse = fit_rmse, estimator.type = ice_colname))

}

#' Split Keyword Arguments
#'
#' @param argument the keyword argument list.
#' @param target_string a string specifying the target keyword argument pattern.
#'
#' @return a list containing the information of the targeted keyword argument.
#' @export
#'
#' @keywords internal
split_args <- function(argument, target_string) {
  split_list <- str_split(names(argument), target_string)
  # print(split_list)
  origin_list <- lapply(split_list, function(x) {x[2]})
  split_by_dot <- str_split(origin_list, "[.]")
  prefix <- lapply(split_by_dot, function(x) {x[1]})
  suffix <- lapply(split_by_dot, function(x) {x[2]})

  return(list(prefix = prefix, suffix = suffix, origin_list = origin_list))
}

#' Get the model formula from keyword arguments
#'
#' @param interv_list the list of names of interventions.
#' @param arg_interv the list of arguments containing the target keyword argument.
#' @param model_interv the list of model inputs related to the target keyword argument.
#' @param ninterv the number of interventions.
#' @param arg_str the argument keyword.
#'
#' @return a list containing the model inputs corresponding to the interventions.
#' @export
#'
#' @keywords internal
get_model_formula <- function(interv_list, arg_interv, model_interv, ninterv, arg_str) {
  models <- as.list(rep(NA, ninterv))

  for (i in 1:ninterv) {

    interv_i <- interv_list[i]
    if (interv_i %in% model_interv) {
      arg_name_i <- paste0(arg_str, ".", interv_i)
      outcome_formula <- arg_interv[[which(names(arg_interv) == arg_name_i)]]
      models[[i]] <- outcome_formula
    }

  }

  return(models)
}


#' Construct intervened values for specified intervention time points
#'
#' @param data a dataframe containing the observed data.
#' @param timevar a character string specifying the time index column.
#' @param int_value a list containing the intervened values for each intervention on all time points.
#' @param int_time a list containing the intervention time points for each intervention.
#' @param int_var a list containing the intervention variable names for each intervention.
#'
#' @return the intervened values on specified intervention time points.
#' @export
#'
#' @keywords internal
construct_interv_value <- function(data, timevar, int_value, int_time, int_var){

  nint <- length(int_var)
  # print(int_var)

  # print(int_value)

  new_int_value <- list()

  for (i in 1:nint) {
    this_int_value <- int_value[[i]]
    not_int_idx <- which(!data[, timevar] %in% int_time[[i]])
    not_int_value <- data[, int_var[[i]]][not_int_idx]
    this_int_value[not_int_idx] <- not_int_value
    new_int_value <- c(new_int_value, list(this_int_value))
  }


  return(new_int_value)
}


#' Match bootstrap values to aggregate data frame
#'
#' @keywords internal
match_boot_values <- function(risk_df, cv_df, col) {
  risk_df[, col] <- NA
  names <- colnames(cv_df)
  for (icv in 1:ncol(cv_df)) {
    cv_name <- names[icv]
    cv_name <- str_replace_all(cv_name, "[.]", " ")
    risk_df[risk_df$Intervention == cv_name, col] <- as.vector(cv_df[, names[icv]])
  }
  
  return(risk_df)
}


#' Append item to list
#'
#' @keywords internal
append_list <- function(item, name, append_list, item_name) {
  item_list <- list(item[[item_name]])
  names(item_list) <- name
  append_list <- c(append_list, item_list)
  
  return(append_list)
}

