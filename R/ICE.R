#' Iterative Conditional Expectation Estimator
#'
#' This function implements iterative conditional expectation (ICE) estimators under user-defined treatment strategies.
#' Available ICE estimators are classical and hazard-based pooling over treatment history ICE, 
#' classical and hazard-based stratifying on treatment history ICE, and doubly robust ICE estimators. 
#' See Wen et al. (2021) for more details regarding the parametric g-formula iterative conditional expectation estimator.
#'
#' Users could specify which version ICE estimator to use through \code{estimator}.
#' \itemize{
#' \item{\code{pool(hazard = F)} specifies the classical pooling over treatment history ICE estimator. }
#' \item{\code{pool(hazard = T)} specifies the hazard-based pooling over treatment history ICE estimator. }
#' \item{\code{strat(hazard = F)} specifies the classical stratifying on treatment history ICE estimator. }
#' \item{\code{strat(hazard = T)} specifies the hazard-based stratifying on treatment history ICE estimator. }
#' \item{\code{weight(treat_model)} specifies the doubly robust weighted ICE estimator
#' where \code{treat_model} specifies the treatment model.}
#' }
#' To provide flexible choices on model inputs for stratified ICE and doubly robust ICE, we allow users to specify intervention-specific model statements through keyword arguments. 
#' In the case where intervention-specific model statements are specified, 
#' treatment variables that are not intervened under some strategies will be considered as a covariate and automatically added into the model specification at each time point.
#' Please see more details on how to specify intervention-specific model specifications in the "Arguments" section.
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
#' \code{intervention2.A1 = list(dynamic("L1 > 0", static(0), static(1), absorb = F)),} \cr
#' \code{intervention2.A2 = list(dynamic("L2 == 0", static(0), static(1), absorb = T)),} \cr
#' \code{outcomeModel.1 = Y ~ L1,} \cr
#' \code{compModel.2 = D ~ L1} \cr
#'
#' Because the keyword argument for outcome model is not specified for intervention 2, the outcome model for intervention 2 is
#' \code{Y ~ L1 + L2} as specified in \code{outcome_model}.
#' Similarly, because the keyword argument for competing model is not specified for intervention 1, the competing model for intervention 1 is
#' \code{D ~ L1 + L2} as specified in \code{competing_model}.
#' In the case of controlled direct effect, the keyword arguments for competing models are ignored. Please see more examples in the "Examples" section.
#'
#' Both built-in interventions and user-defined interventions are available. 
#' 
#' The following are the built-in intervention functions in the package:
#' \itemize{
#' \item{Static Intervention:} {\code{static(value)} specifies a constant intervention with \code{value}.}
#' \item{Dynamic Intervention:} {\cr \code{dynamic(condition, strategy_before, strategy_after, absorb)} 
#' specifies a dynamic intervention where the strategy in \code{strategy_before} is followed until \code{condition} is met. 
#' Upon \code{condition} is met, the strategy in \code{strategy_after} is followed. If absorb is \code{TRUE}, the intervention becomes absorbing once \code{condition} is met.}
#' \item{Threshold Intervention:} {\code{threshold(lower_bound, upper_bound)} specifies a threshold intervention. 
#' If the treatment value is between \code{lower_bound} and \code{upper_bound} inclusively, follow the natural value of the treatment. 
#' Otherwise, set to \code{lower_bound} or \code{upper_bound}, if the treatment value is below \code{lower_bound} or above \code{upper_bound}, correspondingly.}
#' \item{Grace Period:} {\code{grace_period(type, nperiod, condition)} specifies a dynamic intervention with grace period. 
#' Once \code{condition} is met, the intervention is initiated within \code{nperiod} time units. 
#' During the grace period, the treatment variable follows its natural value or initiate intervention with a uniform distribution at each time point.}
#' 
#' The following is the user-defined intervention:
#' 
#' \item{User-defined Interventions:} {The output of the user-defined intervention should contain the intervened value for each individual at each time point, and should be of the same size as the number of rows in \code{data}.
#' 
#' Please see examples in the "Examples" section.
#'
#' In order to obtain an inverse probability (IP) weighted natural course risk based on the observed data, users must specify
#' a censoring variable through \code{censor_name} and a corresponding censoring model through \code{censor_model}. Please see
#' Chiu et al. (2023) for more details regarding the IP weighted estimate of the natural course risk.
#'
#' If competing event exists in the data, users need to specify the name of the competing variable through \code{competing_name} and
#' the model specification through \code{competing_model} for hazard-based ICE estimator. Users need to specify whether to treat
#' the competing event as censoring or total effect through \code{total_effect}.
#' 
#' We provide flexible term options in model specification for the outcome, censoring, competing, and hazard model. 
#' Users could specify polynomial terms using functions \code{I} and \code{poly} and spline terms using \code{ns} from splines package
#' and \code{rcspline.eval} from Hmisc package. In addition, users could specify lagged terms using the format lag\code{n}_\code{var} 
#' to indicate lagging the variable \emph{var} with \emph{n} periods. 
#' If the lagged variable is a treatment variable, this variable is automatically intervened based on user-defined intervention.
#' The polynomial and spline terms could be used on lagged variables.
#'
#'
#'
#' @param data a data frame containing the observed data in long format.
#' @param time_points a number indicating the total number of time points.
#' @param id a string indicating the ID variable name in \code{data}.
#' @param time_name a string specifying the time variable name in \code{data}.
#' @param outcome_name a string specifying the outcome variable name in \code{data}.
#' @param censor_name a string specifying the censor variable name in \code{data}. Default is \code{NULL}.
#' @param compevent_name a string specifying the competing variable name in \code{data}. Default is \code{NULL}.
#' @param comp_effect a number indicating how the competing event is handled for all the specified interventions. Default is 0.
#' 0 for controlled direct effect. 1 for total effect.
#' @param outcome_model a formula specifying the model statement for the outcome. 
#' @param censor_model a formula specifying the model statement for the censoring event. Default is \code{NULL}.
#' @param competing_model a formula specifying the model statement for the competing event. Default is \code{NULL}.
#' @param hazard_model a formula specifying the model statement for the hazard, if hazard-based estimator is used. Default is \code{NULL}. 
#' If specified, the model in \code{hazard_model} will be used. If NULL, the model in \code{outcome_model} will be used.
#' @param global_hazard a logical value indicating whether to use global pooled-over-time hazard model or time-specific hazard models, for hazard-based pooled ICE only. 
#' If \code{TRUE}, use pooled-over-time hazard model. If \code{FALSE}, use time-specific hazard models. Default is \code{FALSE}.
#' @param ref_idx a number indicating which intervention to be used as the reference to calculate the risk ratio and risk difference. Default is 0.
#' 0 refers to the natural course as the reference intervention.
#' Any other numbers refer to the corresponding intervention that users specify in the keyword arguments.
#' @param estimator a function specifying which ICE estimator to use for the estimation. Possible inputs are:
#' \itemize{
#' \item{Classical pooling over treatment history ICE estimator (classical pooled ICE): }{ \code{pool(hazard = F)}}
#' \item{Hazard-Based pooling over treatment history ICE estimator (hazard-based pooled ICE): }{ \code{pool(hazard = T)}}
#' \item{Classical stratifying on treatment history ICE estimator (classical stratified ICE): }{ \code{strat(hazard = F)}}
#' \item{Hazard-Based stratifying on treatment history ICE estimator (hazard-based stratified ICE): }{ \code{strat(hazard = T)}}
#' \item{Doubly robust weighted ICE estimator (doubly robust ICE): }{\cr
#' \code{weight(treat_model)}
#' where \code{treat_model} is a list specifying the treatment model.}
#' }
#' @param int_descript a vector of strings containing descriptions for each specified intervention.
#' @param ci_method a string specifying the method for calculating the confidence interval, if \code{nsamples} is larger than 0.
#' Possible values are "percentile" and "normal." Default is "percentile."
#' @param nsamples a number larger than 0 indicating the number of bootstrap samples. Default is 0.
#' @param seed a number indicating the starting seed for bootstrapping. Default is 1.
#' @param coverage a number greater than 0 and less than 100 indicating the coverage of the confidence interval. Default is 95.
#' @param parallel a logical value indicating whether to parallelize the bootstrap process. Default is \code{FALSE}.
#' @param ncores a number indicating the number of CPU cores to use in parallel computing. Default is 2.
#' @param ... keyword arguments to specify intervention inputs. If stratified ICE is used, keyword arguments also allow intervention-specific outcome models and competing models. 
#' \cr 
#' To specify interventions, please follow the input convention below:
#' \itemize{
#' \item{Each intervention is specified using the keyword argument name with \emph{intervention} prefix.}
#' \item{Use \emph{i} after \emph{intervention} prefix in keyword argument name to represent the ith strategy.}
#' \item{Use \emph{.} followed with \emph{treatment variable name} after \emph{interventioni} in keyword argument name to represent the treatment name within the ith strategy.}
#' }
#' Each input of intervention keyword arguments is a list consisting of a vector of intervened values and an optional vector of time points on which the intervention is applied. 
#' If the intervention time points are not specified, the intervention is applied to all time points. For example, 
#' an input considers a simultaneous intervention with always treat on A1 and never treat on A2 at all time points looks like: \cr 
#' \cr
#' \code{intervention1.A1 = list(static(1))} \cr
#' \code{intervention1.A2 = list(static(0))} \cr \cr
#' The above intervention applies to all time points. The following is an example of custom intervention time points, 
#' with always treat on A1 at time point 1 and 2 and never treat on A2 at time point 3 to 5. \cr 
#' \cr
#' \code{intervention1.A1 = list(static(1), 1:2)} \cr
#' \code{intervention1.A2 = list(static(0), 3:5)} \cr
#' \cr If there is no intervention keyword argument specified, the function returns the natural course risk only. Please see the "Examples" section for more examples.
#' 
#' To specify different outcome model and/or competing model for different intervention, please follow the input convention below:
#'
#' \itemize{
#' \item{Each outcome model is specified using keyword argument name starting with \emph{outcomeModel} or \emph{compModel} prefix for outcome model or competing model correspondingly.}
#' \item{Use \emph{.n} after \emph{outcomeModel} or \emph{compModel} prefix in keyword argument name to specify which intervention being applied to,
#' where \emph{n} represents the \emph{n}th intervention.}
#' }
#' 
#' The input to each outcome or competing model keyword argument is a model statement formula. 
#' If no outcome model and competing model keyword argument is specified, the models specified in \code{outcome_model} and \code{comp_model} are used.
#' Please refer to the "Examples" section for more examples.
#'
#'
#'
#' @return A list containing the following components. Each component that contains the fitted models includes the model fits, the summary of the fitted model, 
#' standard errors of the coefficients, variance-covariance matrices of the parameters, and the root mean square error (RMSE) values.
#' \item{estimator.type}{A string describing the type of the estimator.}
#' \item{summary}{A summary table containing the estimated risk, risk ratio, and risk difference for user-defined interventions including estimated natural course risk and the observed risk.
#' If \code{nsamples} is greater than 0, the summary table includes standard error and confidence interval for the point estimates.}
#' \item{risk.over.time}{A data frame containing the estimated risk at each time point for each intervention.}
#' \item{initial.outcome}{A list, where the name of each sublist corresponds to each specified intervention description, and each sublist contains the fitted models for the outcome model in the first step of algorithm.}
#' \item{initial.comp}{A list, where the name of each sublist corresponds to each specified intervention description, and each sublist contains the fitted models for the competing model in the first step of algorithm (if applicable).}
#' \item{np.risk.model}{A list containing the fitted models for the censoring and/or competing model in estimating observed risk (if applicable).}
#' \item{outcome.models.by.step}{A list, where the name of each sublist corresponds to each specified intervention description, and each sublist contains the fitted models for the outcome model in each iteration of algorithm.}
#' \item{comp.models.by.step}{A list, where the name of each sublist corresponds to each specified intervention description, and each sublist contains the fitted models for the competing model in each iteration of algorithm (if applicable).}
#' \item{hazard.models.by.step}{A list, where the name of each sublist corresponds to each specified intervention description, and each sublist contains the fitted models for the hazard model (if applicable), either time-specific models at all time points or one pooled-over-time global model.} 
#' \item{boot.data}{A list of bootstrap samples. If \code{nsamples} is set to 0, a \code{NULL} value is returned.}
#' \item{boot.initial.outcome}{A list, where the name of each sublist corresponds to each specified intervention description, and each sublist contains the fitted models for the outcome model in the first step of algorithm on the bootstrap samples.
#' If \code{nsamples} is set to 0, a \code{NULL} value is returned.}
#' \item{boot.initial.comp}{A list, where the name of each sublist corresponds to each specified intervention description, and each sublist contains the fitted models for the competing model in the first step of algorithm on the bootstrap samples (if applicable).
#' If \code{nsamples} is set to 0, a \code{NULL} value is returned.}
#' \item{boot.np.risk.model}{A list containing the fitted models for the censoring and/or competing model in estimating observed risk on the bootstrap samples (if applicable).
#' If \code{nsamples} is set to 0, a \code{NULL} value is returned.}
#' \item{boot.outcome.models.by.step}{A list, where the name of each sublist corresponds to each specified intervention description, and each sublist contains the fitted models for the outcome model in each iteration of algorithm on the bootstrap samples (if applicable).
#' If \code{nsamples} is set to 0, a \code{NULL} value is returned.}
#' \item{boot.comp.models.by.step}{A list, where the name of each sublist corresponds to each specified intervention description, and each sublist contains the fitted models for the competing model in each iteration of algorithm on the bootstrap samples (if applicable).
#' If \code{nsamples} is set to 0, a \code{NULL} value is returned.}
#' \item{boot.hazard.models.by.step}{A list, where the name of each sublist corresponds to each specified intervention description, and each sublist contains the fitted models for the hazard model (if applicable), either time-specific models at all time points or one pooled-over-time global model, on the bootstrap samples.
#' If \code{nsamples} is set to 0, a \code{NULL} value is returned.}
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
#' data <- gfoRmulaICE::compData
#' 
#' # Example 1: Dynamic Intervention
#' 
#' # We consider the following interventions and intervened at all time points.
#' # Intervention 1 on A2: at time t, if L1 = 0, then treat; otherwise, not treat. 
#' # Intervention 2 on A2: never treat upon until L1 = 0, after which follows always treat.
#' # Intervention 3 on A2: never treat upon until L1 = 0, after which follows natural course.
#' 
#' # We use classical pooled ICE estimator, 
#' # natural course as the reference intervention, and the following models:
#' # a. outcome model: Y ~ L1 + L2 + A1 + A2
#' # b. censor model: C ~ L1 + L2 + A1 + A2
#' # c. competing model: D ~ L1 + L2 + A1 + A2.
#' # We estimate variance using bootstrap with 1000 replicates, normal quantile, and parallel computing.
#' 
#' ice_fit1 <- ice(data = data, time_points = 4, 
#' id = "id", time_name = "t0",
#' censor_name = "C", outcome_name = "Y",
#' compevent_name = "D",
#' comp_effect = 0,
#' outcome_model = Y ~ L1 + L2 + A1 + A2, 
#' censor_model = C ~ L1 + L2 + A1 + A2,
#' ref_idx = 0,
#' estimator = pool(hazard = F),
#' nsamples = 1000, ci_method = "percentile",
#' parallel = T, ncores = 5,
#' int_descript = c("Dynamic Intervention 1", "Dynamic Intervention 2", 
#' "Dynamic Intervention 3"),
#' intervention1.A2 = list(dynamic("L1 == 0", static(0), static(1))),
#' intervention2.A2 = list(dynamic("L1 == 0", static(0), static(1), absorb = T)),
#' intervention3.A2 = list(dynamic("L1 == 0", static(0), natural_course()))
#' )
#' 
#' plot_risk(ice_fit1)
#' 
#' # Example 2: Built-in Interventions
#' 
#' # We consider the following interventions and intervene at all time points.
#' # Intervention 1 on A1: always treat with value 3.
#' # Intervention 1 on A2: always treat with value 1.
#' # Intervention 2 on L2: when the natural value of L2 at time t is lower than -3, set its value to -3. 
#' # Otherwise, do not intervene.
#' # Intervention 3 on A2: dynamic intervention (treat when L1 = 0) with uniform grace period of 2 periods
#' 
#' # We use classical pooled ICE estimator, 
#' # natural course as the reference intervention, and the following models:
#' # a. outcome model: Y ~ L1 + L2 + A1 + A2
#' # b. censor model: C ~ L1 + L2 + A1 + A2
#' # c. competing model: D ~ L1 + L2 + A1 + A2.
#' # We estimate variance using bootstrap with 1000 replicates, normal quantile, and parallel computing.
#' 
#' ice_fit2 <- ice(data = data, time_points = 4, 
#' id = "id", time_name = "t0",
#' censor_name = "C", outcome_name = "Y",
#' compevent_name = "D",
#' comp_effect = 0,
#' outcome_model = Y ~ L1 + L2 + A1 + A2, 
#' censor_model = C ~ L1 + L2 + A1 + A2,
#' ref_idx = 0,
#' estimator = pool(hazard = F),
#' nsamples = 1000, ci_method = "percentile",
#' parallel = T, ncores = 5,
#' int_descript = c("Static Intervention", "Threshold Intervention", 
#' "Dynamic Intervention with Grace Period"),
#' intervention1.A1 = list(static(3)),
#' intervention1.A2 = list(static(1)),
#' intervention2.L2 = list(threshold(-3, Inf)),
#' intervention3.A2 = list(grace_period("uniform", 2, "L1 == 0"))
#' )
#' 
#' plot_risk(ice_fit2)
#' 
#' # Example 3: User-defined Intervention
#' 
#' # We consider the following interventions and intervene at all time points.
#' # Intervention 1 on A1: always treat with value 3.
#' # Intervention 1 on A2: always treat with value 1.
#' # Intervention 2 on A1: at time t, if L2 < 0, then assign 1; if 0 <= L2 < 2, then assign 2; otherwise, assign 3.
#' # Intervention 2 on A2: at time t, if L1 = 0, then treat; otherwise, not treat. 
#' 
#' # We use classical pooled ICE estimator, 
#' # natural course as the reference intervention, and the following models:
#' # a. outcome model: Y ~ L1 + L2 + A1 + A2
#' # b. censor model: C ~ L1 + L2 + A1 + A2
#' # c. competing model: D ~ L1 + L2 + A1 + A2.
#' # We estimate variance using bootstrap with 1000 replicates and percentile quantile.
#' 
#' dynamic_cat <- case_when(data$L2 < 0 ~ 1,
#' data$L2 >= 0 & data$L2 < 2 ~ 2, T ~ 3)
#' 
#' ice_fit3 <- ice(data = data, time_points = 4, 
#' id = "id", time_name = "t0",
#' censor_name = "C", outcome_name = "Y",
#' compevent_name = "D",
#' comp_effect = 0,
#' outcome_model = Y ~ L1 + L2 + A1 + A2, 
#' censor_model = C ~ L1 + L2 + A1 + A2,
#' ref_idx = 0,
#' estimator = pool(hazard = F),
#' nsamples = 1000, ci_method = "percentile",
#' parallel = T, ncores = 5,
#' int_descript = c("Static Intervention", "Dynamic Intervention"),
#' intervention1.A1 = list(static(3)),
#' intervention1.A2 = list(static(1)),
#' intervention2.A1 = list(dynamic_cat),
#' intervention2.A2 = list(dynamic("L1 == 0", static(0), static(1)))
#' )
#' 
#' plot_risk(ice_fit3)
#' 
#' # Example 4: Different ICE Estimators
#' 
#' # We use the interventions in Example 3 and implement each ICE estimator.
#' 
#' # a. hazard-based pooled ICE:
#' # hazard model is time-specific and shares the same model statement as the outcome model
#' 
#' ice_fit4a <- ice(data = data, time_points = 4, 
#' id = "id", time_name = "t0",
#' censor_name = "C", outcome_name = "Y",
#' compevent_name = "D",
#' comp_effect = 0,
#' outcome_model = Y ~ L1 + L2 + A1 + A2, 
#' censor_model = C ~ L1 + L2 + A1 + A2,
#' competing_model = D ~ L1 + L2 + A1 + A2,
#' ref_idx = 0,
#' estimator = pool(hazard = T),
#' nsamples = 1000, ci_method = "percentile",
#' parallel = T, ncores = 5,
#' int_descript = c("Static Intervention", "Dynamic Intervention"),
#' intervention1.A1 = list(static(3)),
#' intervention1.A2 = list(static(1)),
#' intervention2.A1 = list(dynamic_cat),
#' intervention2.A2 = list(dynamic("L1 == 0", static(0), static(1)))
#' )
#' 
#' plot_risk(ice_fit4a)
#' 
#' # b. hazard-based pooled ICE: 
#' # hazard model is time-specific and uses Y ~ L1 + L2
#' 
#' ice_fit4b <- ice(data = data, time_points = 4, 
#' id = "id", time_name = "t0",
#' censor_name = "C", outcome_name = "Y",
#' compevent_name = "D",
#' comp_effect = 0,
#' outcome_model = Y ~ L1 + L2 + A1 + A2, 
#' censor_model = C ~ L1 + L2 + A1 + A2,
#' competing_model = D ~ L1 + L2 + A1 + A2,
#' hazard_model = Y ~ L1 + L2,
#' ref_idx = 0,
#' estimator = pool(hazard = T),
#' nsamples = 1000, ci_method = "percentile",
#' parallel = T, ncores = 5,
#' int_descript = c("Static Intervention", "Dynamic Intervention"),
#' intervention1.A1 = list(static(3)),
#' intervention1.A2 = list(static(1)),
#' intervention2.A1 = list(dynamic_cat),
#' intervention2.A2 = list(dynamic("L1 == 0", static(0), static(1)))
#' )
#' 
#' plot_risk(ice_fit4b)
#' 
#' # c. hazard-based pooled ICE: 
#' # hazard model is pooled-over-time and includes flexible terms of time variable
#' 
#' library(splines)
#' 
#' ice_fit4c <- ice(data = data, time_points = 4, 
#' id = "id", time_name = "t0",
#' censor_name = "C", outcome_name = "Y",
#' compevent_name = "D",
#' comp_effect = 0,
#' outcome_model = Y ~ L1 + L2 + A1 + A2, 
#' censor_model = C ~ L1 + L2 + A1 + A2,
#' competing_model = D ~ L1 + L2 + A1 + A2,
#' hazard_model = Y ~ L1 + L2 + A1 + A2 + ns(t0, df = 2),
#' global_hazard = T,
#' ref_idx = 0,
#' estimator = pool(hazard = T),
#' nsamples = 1000, ci_method = "percentile",
#' parallel = T, ncores = 5,
#' int_descript = c("Static Intervention", "Dynamic Intervention"),
#' intervention1.A1 = list(static(3)),
#' intervention1.A2 = list(static(1)),
#' intervention2.A1 = list(dynamic_cat),
#' intervention2.A2 = list(dynamic("L1 == 0", static(0), static(1)))
#' )
#' 
#' plot_risk(ice_fit4c)
#' 
#' # d. classical stratified ICE:
#' 
#' ice_fit4d <- ice(data = data, time_points = 4, 
#' id = "id", time_name = "t0",
#' censor_name = "C", outcome_name = "Y",
#' compevent_name = "D",
#' comp_effect = 0,
#' outcome_model = Y ~ L1 + L2, 
#' censor_model = C ~ L1 + L2,
#' ref_idx = 0,
#' estimator = strat(hazard = F),
#' nsamples = 1000, ci_method = "percentile",
#' parallel = T, ncores = 5,
#' int_descript = c("Static Intervention", "Dynamic Intervention"),
#' intervention1.A1 = list(static(3)),
#' intervention1.A2 = list(static(1)),
#' intervention2.A1 = list(dynamic_cat),
#' intervention2.A2 = list(dynamic("L1 == 0", static(0), static(1)))
#' )
#' 
#' plot_risk(ice_fit4d)
#' 
#' # e. hazard-based stratified ICE:
#' # hazard model is time-specific and uses Y ~ L1 
#' # (Note: a pooled-over-time hazard model is not valid for stratified ICE.)
#' 
#' ice_fit4e <- ice(data = data, time_points = 4, 
#' id = "id", time_name = "t0",
#' censor_name = "C", outcome_name = "Y",
#' compevent_name = "D",
#' comp_effect = 0,
#' outcome_model = Y ~ L1 + L2, 
#' censor_model = C ~ L1 + L2,
#' competing_model = D ~ L1 + L2,
#' hazard_model = Y ~ L1,
#' ref_idx = 0,
#' estimator = strat(hazard = T),
#' nsamples = 1000, ci_method = "percentile",
#' parallel = T, ncores = 5,
#' int_descript = c("Static Intervention: Model 1", 
#' "Dynamic Intervention: Model 1"),
#' intervention1.A1 = list(static(3)),
#' intervention1.A2 = list(static(1)),
#' intervention2.A1 = list(dynamic_cat),
#' intervention2.A2 = list(dynamic("L1 == 0", static(0), static(1)))
#' )
#' 
#' plot_risk(ice_fit4e)
#' 
#' 
#' # f. doubly robust ICE:
#' 
#' ice_fit4f <- ice(data = data, time_points = 4, 
#' id = "id", time_name = "t0",
#' censor_name = "C", outcome_name = "Y",
#' compevent_name = "D",
#' comp_effect = 0,
#' outcome_model = Y ~ L1 + L2, 
#' censor_model = C ~ L1 + L2,
#' ref_idx = 0,
#' estimator = weight(list(A1 ~ L1 + L2, A2 ~ L1 + L2)),
#' nsamples = 1000, ci_method = "percentile",
#' parallel = T, ncores = 5,
#' int_descript = c("Static Intervention", 
#' "Dynamic Intervention"),
#' intervention1.A1 = list(static(3)),
#' intervention1.A2 = list(static(1)),
#' intervention2.A1 = list(dynamic_cat),
#' intervention2.A2 = list(dynamic("L1 == 0", static(0), static(1)))
#' )
#' 
#' plot_risk(ice_fit4f)
#' 
#' 
#' # g. hazard-based stratified ICE with intervention-specific models:
#' # hazard model is time-specific and same as the outcome model
#' # consider the total effect for competing event,
#' # using normal quantile for variance estimates,
#' # and the following outcome models and competing models:
#' # outcome model for intervention 1: Y ~ L1,
#' # outcome model for intervention 2: Y ~ L1 + L2,
#' # competing model for intervention 1: D ~ L1 + L2,
#' # competing model for intervention 2: D ~ L1
#' 
#' ice_fit4g <- ice(data = data, time_points = 4, 
#' id = "id", time_name = "t0",
#' censor_name = "C", outcome_name = "Y",
#' compevent_name = "D",
#' outcome_model = Y ~ L1, censor_model = C ~ L1,
#' competing_model = D ~ L1,
#' comp_effect = 1,
#' ref_idx = 0,
#' estimator = strat(hazard = T),
#' nsamples = 1000, ci_method = "normal",
#' parallel = T, ncores = 5,
#' int_descript = c("Static Intervention: Model 2",
#' "Dynamic Intervention: Model 2"),
#' intervention1.A1 = list(static(3)),
#' intervention1.A2 = list(static(1)),
#' intervention2.A1 = list(dynamic_cat),
#' intervention2.A2 = list(dynamic("L1 == 0", static(0), static(1))),
#' outcomeModel.1 = Y ~ L1 + L2,
#' compModel.2 = D ~ L1 + L2
#' )
#'
#' # Compare with the ICE estimates in Example 4e:
#' plot_risk(ice_fit4e, ice_fit4g)
#' summary_table(ice_fit4e, ice_fit4g)
#' 
#' # Example 5: Flexible Model Specification
#' 
#' # a. Complicated terms in model statement:
#' # We use the same interventions and ICE estimator in Example 3, 
#' # and include polynomial, spline, and lagged terms in models.
#' 
#' ice_fit5a <- ice(data = data, time_points = 4, 
#' id = "id", time_name = "t0",
#' censor_name = "C", outcome_name = "Y",
#' compevent_name = "D",
#' comp_effect = 0,
#' outcome_model = Y ~ I(L1^2) + rcspline.eval(lag1_L2, knots = 1:3) + A1 + A2,
#' censor_model = C ~ lag1_L1 + poly(L2, degree = 2) + A1 + A2,
#' ref_idx = 0,
#' estimator = pool(hazard = F),
#' nsamples = 1000, ci_method = "percentile",
#' parallel = T, ncores = 5,
#' int_descript = c("Static Intervention", "Dynamic Intervention"),
#' intervention1.A1 = list(static(3)),
#' intervention1.A2 = list(static(1)),
#' intervention2.A1 = list(dynamic_cat),
#' intervention2.A2 = list(dynamic("L1 == 0", static(0), static(1)))
#' )
#' 
#' plot_risk(ice_fit5a)
#' 
#' # b. Using static intervention as reference:
#' # We use the same interventions and ICE estimator in Example 3, 
#' # but use static intervention as the reference intervention.
#' 
#' ice_fit5b <- ice(data = data, time_points = 4, 
#' id = "id", time_name = "t0",
#' censor_name = "C", outcome_name = "Y",
#' compevent_name = "D",
#' comp_effect = 0,
#' outcome_model = Y ~ I(L1^2) + rcspline.eval(lag1_L2, knots = 1:3) + A1 + A2,
#' censor_model = C ~ lag1_L1 + poly(L2, degree = 2) + A1 + A2,
#' ref_idx = 1,
#' estimator = pool(hazard = F),
#' nsamples = 1000, ci_method = "percentile",
#' parallel = T, ncores = 5,
#' int_descript = c("Static Intervention", "Dynamic Intervention"),
#' intervention1.A1 = list(static(3)),
#' intervention1.A2 = list(static(1)),
#' intervention2.A1 = list(dynamic_cat),
#' intervention2.A2 = list(dynamic("L1 == 0", static(0), static(1)))
#' )
#' 
#' plot_risk(ice_fit5b)
#' 
#' 
#' @import tidyverse rlang dplyr stringr data.table reshape2 nnet
#' @export

ice <- function(data, time_points, id, time_name,
                outcome_name, censor_name = NULL,
                compevent_name = NULL, comp_effect = 0,
                outcome_model, censor_model = NULL, competing_model = NULL,
                hazard_model = NULL, global_hazard = F, 
                ref_idx = 0,
                estimator,
                int_descript,
                ci_method = "percentile",
                nsamples = 0, seed = 1,
                significance_level = 0.05, parallel = F, ncores = 2,
                ...) {
  
  interv_data <<- data
  idvar <<- id
  time0var <<- time_name
  outcomevar <<- outcome_name
  
  ## check if argument names are specified correctly
  input_args <- as.list(environment())
  input_arg_names <- names(input_args)
  named_args <- formalArgs(ice)
  
  for (i in 1:length(input_arg_names)) {
    
    iarg <- input_arg_names[i]
    
    if (!iarg %in% named_args) {
      stop(paste0("Please check the input argument name: ", iarg, ". An error occurs probably because of a typo in argument name."))
    }
  }
  
  ## pre-process user inputs
  
  data <- as.data.frame(data)

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
  
  # preprocess the hazard model
  
  if (any(str_detect(as.character(substitute(estimator)), "pool")) | 
      any(str_detect(as.character(substitute(estimator)), "strat")) & 
      any(str_detect(as.character(substitute(estimator)), "T")) & is.null(hazard_model)) {
    
    hazard_model <- outcome_model
  } else if (any(str_detect(as.character(substitute(estimator)), "pool")) | 
             any(str_detect(as.character(substitute(estimator)), "strat")) & 
             any(str_detect(as.character(substitute(estimator)), "F")) & !is.null(hazard_model)) {
    warning("Hazard model is ignored for Classical ICE estimators.")
  }
  

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
    
    # replace the time variable in the global hazard model with the new time name if any
    if (global_hazard) {

      global_hazard_model_str <- as.character(hazard_model)
      global_hazard_model_str[3] <- str_replace_all(global_hazard_model_str[3], time_name, paste0("new_", time_name))
      hazard_model <- as.formula(paste0(global_hazard_model_str[c(2, 1, 3)], collapse = ""))
    }
    
    time_name <- paste0("new_", time_name)
  } else {

    # if the time does not start with 0, then re-index
    if (min(unique_time_names) != 0) {
      data <- map_time_column(time_name, data, unique_time_names)
      
      # replace the time variable in the global hazard model with the new time name if any
      if (global_hazard) {
        
        global_hazard_model_str <- as.character(hazard_model)
        global_hazard_model_str[3] <- str_replace_all(global_hazard_model_str[3], time_name, paste0("new_", time_name))
        hazard_model <- as.formula(paste0(global_hazard_model_str[c(2, 1, 3)], collapse = ""))
      }
      
      time_name <- paste0("new_", time_name)
    }
  }
  
  
  ## preprocess threshold intervention

  sub_kwarg <- substitute(list(...))
  clean_kwarg_str <- paste0(str_remove_all(as.character(deparse(sub_kwarg)), " "), collapse = "")
  kwarg_list <- as.list(sub_kwarg)
  kwarg_name_list <- names(kwarg_list)
  threshold_idx <- str_which(str_split(as.character(substitute(list(...))), " = "), "threshold")
  
  ## preprocess natural course - need to add treatment variable into the argument
  
  nc_idx <- str_which(str_split(as.character(substitute(list(...))), " = "), "natural_course")
  
  if (length(nc_idx) > 0) {
    nc_kwarg <- kwarg_name_list[nc_idx]
    nc_interv_list_split <- str_split(nc_kwarg, "[.]")
    nc_interv_list <- lapply(nc_interv_list_split, function(x) {x[1]})
    nc_interv_list_all_names <- lapply(nc_interv_list_split, function(x) {x[2]})

    for (i in 1:length(nc_idx)) {
      ikwarg <- nc_kwarg[i]
      ivar <- as.character(nc_interv_list_all_names[i])
      raw_list <- str_remove_all(as.character(kwarg_list[ikwarg]), " ")
      ikwarg_start_idx <- str_locate(pattern = paste0(ikwarg, "="), clean_kwarg_str)[1, 2]
      ikwarg_end_idx <- ikwarg_start_idx + nchar(raw_list) + 1
      nc_arg_idx <- str_locate_all(pattern = "natural_course[\\(][\\)]", raw_list)[[1]]
      num_nc <- nrow(nc_arg_idx)
      add_string <- paste0("treat_var=\"", ivar, "\"")
      
      ## add the treatment variable in each call of natural course
      
      if (num_nc > 0) {
      for (j in 1:num_nc) {
        if (j == 1) {
        idx_replace <- nc_arg_idx[j, 2] - 1
        } else {
          idx_replace <- str_locate_all(pattern = "natural_course[\\(][\\)]", raw_list)[[1]][1, 2] - 1
        }
        
        raw_list <- paste0(substr(raw_list, 1, idx_replace), add_string, substr(raw_list, idx_replace + 1, nchar(raw_list)))
      }
      } 
      
      clean_kwarg_str <- paste0(substr(clean_kwarg_str, 1, ikwarg_start_idx), raw_list, 
                                substr(clean_kwarg_str, ikwarg_end_idx, nchar(clean_kwarg_str)))
    }
  }

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
  outcomeModels <- get_model_formula(interv_list, outcome_interv, outcomeModel_interv, ninterv, "outcomeModel", outcome_model)

  # make competing model list
  compModels <- get_model_formula(interv_list, comp_interv, compModel_interv, ninterv, "compModel", competing_model)


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
      ice_colname <- "Hazard-Based Pooled ICE"
    } else {
      ice_colname <- "Classical Pooled ICE"
    }

  } else if (any(str_detect(as.character(substitute(estimator)), "strat"))) {
    if (estimator) {
      ice_colname <- "Hazard-Based Stratified ICE"
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
    
    # if hazard based but no hazard model, then set it to outcome model
    if (hazard & is.null(hazard_model)) {
      
      hazard_model <- outcome_model
      
    }
    

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
    outcome_by_step <- outcome_init <- comp_init <- np_model <- c()

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
                    hazard_model = hazard_model, 
                    interventions = nc_interventions, intervention_names = nc_intervention_varlist,
                    compute_nc_risk = T, 
                    hazard_based = hazard, 
                    global_hazard = global_hazard,
                    intervention_description = nc_descript) 
                    #global_haz_model = global_hazard_model)

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
                    hazard_model = hazard_model, 
                    interventions = boot_interv, intervention_names = ref_intervention_varlist,
                    intervention_times = ref_int_times,
                    compute_nc_risk = T, hazard_based = hazard, 
                    global_hazard = global_hazard, 
                    intervention_description = ref_description)
                    #global_haz_model = global_hazard_model)

      if (bootstrap) {
        summary[ref_idx + 1, c(2:4, 6:7, 10:13)] <- c(ref$gformula_risk_last_time, 1, 0, 1, 0, 1, 1, 0, 0)
      } else {
        summary[ref_idx + 1, c(2:4)] <- c(ref$gformula_risk_last_time, 1, 0)
      }
    }

    this_outcome_init <- list(ref$outcome_init)
    this_comp_init <- list(ref$comp_init)
    this_np_model <- list(ref$np_model)
    
    names(this_outcome_init) <- names(this_np_model) <- names(this_comp_init) <- ref_description
    
    
    
    this_outcome_by_step <- ref$outcome_by_step
    
    this_model <- this_outcome_by_step$fit
    this_summary <- this_outcome_by_step$summary
    this_stderr <- this_outcome_by_step$stderr
    this_vcov <- this_outcome_by_step$vcov
    this_rmse <- this_outcome_by_step$rmse

    this_fit_all <- list(list(fit = this_model, 
                         summary = this_summary, 
                         stderr = this_stderr, 
                         vcov = this_vcov, 
                         rmse = this_rmse))
    
    names(this_fit_all) <- ref_description
    
    outcome_by_step <- c(outcome_by_step, this_fit_all)

    # fit_all <- c(fit_all, this_model)
    # fit_summary <- c(fit_summary, this_summary)
    # fit_stderr <- c(fit_stderr, this_stderr)
    # fit_vcov <- c(fit_vcov, this_vcov)
    # fit_rmse <- c(fit_rmse, this_rmse)
    outcome_init <- c(outcome_init, this_outcome_init)
    comp_init <- c(comp_init, this_comp_init)
    np_model <- c(np_model, this_np_model)

    risk_descript <- c(risk_descript, c(rep(str_to_title(ref_description), K+1), rep("Natural Course (nonparametric)", K+1)))
    risk_interv <- c(risk_interv, c(ref$gformula_risk[1, ], c(0, ref$weight_h$risk)))

    if (bootstrap) {
      intervention_descriptions <- c(intervention_descriptions, list("Natural Course"))
    }

    ninterv_new <- length(intervention_descriptions)

    ice_critical_value <- rr_critical_value <- rd_critical_value  <- matrix(NA, nrow = 2, ncol = ninterv + 1)

    if (bootstrap) {
      outcome_init_boot <- comp_init_boot <- np_model_boot <- outcome_by_step_boot <- 
        comp_by_step_boot <- hazard_by_step_boot <- data_boot_all <- c()
    } else {
      outcome_init_boot <- comp_init_boot <- np_model_boot <- outcome_by_step_boot <- 
        comp_by_step_boot <- hazard_by_step_boot <- data_boot_all <- NULL
    }

    if (ninterv_new > 0) {
      
      critical_value_all_lower <- critical_value_all_upper <- se_all <- c()
      
      if (ref_idx != 0) {
        idx_list <- rev(1:ninterv_new)
      } else {
        idx_list <- 1:ninterv_new
      }
      

    for (int in idx_list) {
      
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
                           hazard_model = hazard_model, 
                            interventions = this_interv, intervention_names = this_int_var, intervention_times = this_time,
                            hazard_based = hazard, intervention_description = this_descript, 
                           global_hazard = global_hazard) 
                           # global_haz_model = global_hazard_model)
      
      this_outcome_init <- list(this_fit$outcome_init)
      this_comp_init <- list(this_fit$comp_init)
      this_np_model <- list(this_fit$np_model)
      
      names(this_outcome_init) <- names(this_np_model) <- names(this_comp_init) <- this_descript

      this_fit_outcome <- this_fit$outcome_by_step
      
      this_model <- this_fit_outcome$fit
      this_summary <- this_fit_outcome$summary
      this_stderr <- this_fit_outcome$stderr
      this_vcov <- this_fit_outcome$vcov
      this_rmse <- this_fit_outcome$rmse
      
      this_fit_all <- list(list(fit = this_model, 
                           summary = this_summary, 
                           stderr = this_stderr, 
                           vcov = this_vcov, 
                           rmse = this_rmse))
      
      names(this_fit_all) <- this_descript
      
      outcome_by_step <- c(outcome_by_step, this_fit_all)
      
      # this_model <- list(this_fit_outcome$fit)
      # names(this_model) <- this_descript
      # 
      # this_summary <- list(this_fit_outcome$summary)
      # names(this_summary) <- this_descript
      # 
      # this_stderr <- list(this_fit_outcome$stderr)
      # names(this_stderr) <- this_descript
      # 
      # this_vcov <- list(this_fit_outcome$vcov)
      # names(this_vcov) <- this_descript
      # 
      # this_rmse <- list(this_fit_outcome$rmse)
      # names(this_rmse) <- this_descript
      # 
      # fit_all <- c(fit_all, this_model)
      # fit_summary <- c(fit_summary, this_summary)
      # fit_stderr <- c(fit_stderr, this_stderr)
      # fit_vcov <- c(fit_vcov, this_vcov)
      # fit_rmse <- c(fit_rmse, this_rmse)
      outcome_init <- c(outcome_init, this_outcome_init)
      comp_init <- c(comp_init, this_comp_init)
      np_model <- c(np_model, this_np_model)


      if (this_descript != "Natural Course") {
        summary[int+1, 2] <- this_fit$gformula_risk_last_time
        risk_descript <- c(risk_descript, rep(str_to_title(this_descript), K+1))
        risk_interv <- c(risk_interv, this_fit$gformula_risk[1, ])

      }


    if (bootstrap) {
      
      if (ref_idx == 0) {
        ref_time <- NULL
      } else {
        ref_time <- ref_int_times
      }
      
      # if this intervention is the reference intervention
      if (ref_idx == int & int != idx_list[1] | ref_idx == 0 & this_descript == "Natural Course") {
        
        
        this_boot <- list(ice_se = ref_first_boot$ref_se, 
                          ref_se = ref_first_boot$ref_se,
                          rr_se = 1,
                          rd_se = 0, 
                          ice_cv_all_upper = ref_first_boot$ref_cv_all_upper, 
                          ice_cv_all_lower = ref_first_boot$ref_cv_all_lower, 
                          # ice_cv_upper = tail(ref_first_boot$ref_cv_all_upper, 1), 
                          # ice_cv_lower = tail(ref_first_boot$ref_cv_all_lower, 1), 
                          ref_cv_all_upper = ref_first_boot$ref_cv_all_upper, 
                          ref_cv_all_lower = ref_first_boot$ref_cv_all_lower, 
                          # ref_cv_upper = tail(ref_first_boot$ref_cv_all_upper, 1), 
                          # ref_cv_lower = tail(ref_first_boot$ref_cv_all_lower, 1), 
                          ref_ipw_se = ref_first_boot$ref_ipw_se, 
                          ref_ipw_cv_all_upper = ref_first_boot$ref_ipw_cv_all_upper,
                          ref_ipw_cv_all_lower = ref_first_boot$ref_ipw_cv_all_lower,
                          rr_cv_upper = 1, 
                          rr_cv_lower = 1, 
                          rd_cv_upper = 0, 
                          rd_cv_lower = 0, 
                          outcome_init = ref_first_boot$ref_outcome_init, 
                          comp_init = ref_first_boot$ref_comp_init, 
                          np_model = ref_first_boot$np_model, 
                          outcome_by_step = ref_first_boot$ref_outcome_by_step, 
                          comp_by_step = ref_first_boot$ref_comp_by_step, 
                          hazard_by_step = ref_first_boot$ref_hazard_by_step)
        
        give_warning(ref_first_boot$ref_data_err, ref_first_boot$ref_data_err_mssg_combine)
        give_warning(ref_first_boot$ref_model_err, ref_first_boot$ref_model_err_mssg_combine)
        
      } else {

      this_boot <- bootstrap_ice(ice_pool, K, nboot, significance_level, parallel, ncores, ref_description,
                                 ref_intervention_varlist, this_total_effect, boot_interv,
                                 this_interv, this_int_var, this_descript, this_time, ref_time,
                                 data, id, set_seed,
                                 time_name = time_name, outcome_name = outcome_name,
                                 censor_name = censor_name, competing_name = competing_name,
                                 hazard_based = hazard, outcome_model = outcome_model,
                                 censor_model = censor_model, competing_model = competing_model,
                                 hazard_model = hazard_model, 
                                 global_hazard = global_hazard)
      }
      
      # record the reference intervention bootstrap
      if (int == idx_list[1]) {
        
        ref_first_boot <- this_boot
        
      } 

      this_se <- tail(this_boot$ice_se, 1)
      this_rr_se <- this_boot$rr_se
      this_rd_se <- this_boot$rd_se

      this_cv_upper <- tail(this_boot$ice_cv_all_upper, 1)
      this_cv_lower <- tail(this_boot$ice_cv_all_lower, 1)
      this_rr_cv_upper <- this_boot$rr_cv_upper
      this_rr_cv_lower <- this_boot$rr_cv_lower
      this_rd_cv_upper <- this_boot$rd_cv_upper
      this_rd_cv_lower <- this_boot$rd_cv_lower
      
      this_outcome_init_boot <- list(this_boot$outcome_init)
      this_comp_init_boot <- list(this_boot$comp_init)
      this_np_model_boot <- list(this_boot$np_model)
      this_outcome_by_step_boot <- list(this_boot$outcome_by_step)
      this_comp_by_step_boot <- list(this_boot$comp_by_step)
      this_hazard_by_step_boot <- list(this_boot$hazard_by_step)
      
      names(this_outcome_init_boot) <- names(this_np_model_boot) <- names(this_comp_init_boot) <- 
        names(this_outcome_by_step_boot) <- names(this_comp_by_step_boot) <- names(this_hazard_by_step_boot) <- this_descript
      
      outcome_init_boot <- c(outcome_init_boot, this_outcome_init_boot)
      comp_init_boot <- c(comp_init_boot, this_comp_init_boot)
      np_model_boot <- c(np_model_boot, this_np_model_boot)
      outcome_by_step_boot <- c(outcome_by_step_boot, this_outcome_by_step_boot)
      comp_by_step_boot <- c(comp_by_step_boot, this_comp_by_step_boot)
      hazard_by_step_boot <- c(hazard_by_step_boot, this_hazard_by_step_boot)


      # this_model_boot <- list(this_boot$boot_models)
      # names(this_model_boot) <- this_descript
      # 
      # this_summary_boot <- list(this_boot$boot_summary)
      # names(this_summary_boot) <- this_descript
      # 
      # this_stderr_boot <- list(this_boot$boot_stderr)
      # names(this_stderr_boot) <- this_descript
      # 
      # this_vcov_boot <- list(this_boot$boot_vcov)
      # names(this_vcov_boot) <- this_descript
      # 
      # this_rmse_boot <- list(this_boot$boot_rmse)
      # names(this_rmse_boot) <- this_descript
      # 
      this_data_boot <- list(this_boot$boot_data)
      names(this_data_boot) <- this_descript
      # 
      # fit_all_boot <- c(fit_all_boot, this_model_boot)
      # fit_summary_boot <- c(fit_summary_boot, this_summary_boot)
      # fit_stderr_boot <- c(fit_stderr_boot, this_stderr_boot)
      # fit_vcov_boot <- c(fit_vcov_boot, this_vcov_boot)
      # fit_rmse_boot <- c(fit_rmse_boot, this_rmse_boot)
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
          critical_value_all_upper <- append_list(this_boot, "Natural Course (nonparametric)", critical_value_all_upper, "ref_ipw_cv_all_upper")
          critical_value_all_lower <- append_list(this_boot, "Natural Course (nonparametric)", critical_value_all_lower, "ref_ipw_cv_all_lower")
          se_all <- append_list(this_boot, "Natural Course (nonparametric)", se_all, "ref_ipw_se")
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
          critical_value_all_upper <- append_list(this_boot, "Natural Course (nonparametric)", critical_value_all_upper, "ref_ipw_cv_all_upper")
          critical_value_all_lower <- append_list(this_boot, "Natural Course (nonparametric)", critical_value_all_lower, "ref_ipw_cv_all_lower")
          se_all <- append_list(this_boot, "Natural Course (nonparametric)", se_all, "ref_ipw_se")
        }
      }

    }
      if (bootstrap) {
        summary[ref_idx + 1, 5] <- this_boot$ref_se[K+1]

        if (ref_idx == 0) {
          ice_critical_value[1, 1] <- tail(this_boot$ref_cv_all_lower, 1)
          ice_critical_value[2, 1] <- tail(this_boot$ref_cv_all_upper, 1)
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



    comp_by_step <- hazard_by_step <- c()
    
    # outcome_by_step <- list(fit = fit_all, 
    #                            summary = fit_summary, 
    #                            stderr = fit_stderr, 
    #                            vcov = fit_vcov, 
    #                            rmse = fit_rmse)
    
    

  } else if (any(str_detect(as.character(substitute(estimator)), "strat")) | any(str_detect(as.character(substitute(estimator)), "weight"))) {

    risk_descript <- c()
    risk_interv <- c()
    comp_by_step <- outcome_by_step <- hazard_by_step <- np_model <- outcome_init <- comp_init <- c()

    if (any(str_detect(as.character(substitute(estimator)), "strat"))) {
      weight <- F
      this_treat_model <- list()
      this_obs_treatment_varnames <- as.list(unique(unlist(intervention_varnames)))
      hazard <- estimator
      
      if (global_hazard & hazard) {
        warning("Only time specific hazard model is valid for Stratified ICE.")
      }
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
      if ((any(str_detect(as.character(substitute(estimator)), "strat")) & (hazard == F) & (ref_total_effect == F)) | (any(str_detect(as.character(substitute(estimator)), "weight")) & (hazard == F) & (ref_total_effect == F))) {
      warning("The competing model is used for nonparametric risk estimation for direct effect case in stratified ICE. The keyword argument competing model statments are ignored.")
      } 
    }

    ref <- ice_strat(data = data, K = K, id = id, time_name = time_name, outcome_name = outcome_name,
                     censor_name = censor_name, competing_name = competing_name, total_effect = ref_total_effect,
                     outcome_model = outcome_model, censor_model = censor_model, competing_model = competing_model,
                     hazard_model = hazard_model,
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
      ref_outcome_model <- outcomeModels[[ref_idx]]
      ref_competing_model <- compModels[[ref_idx]]

      if (!is.character(unlist(ref_intervention_varlist))) {
        stop("Please input the treatment variable at the second element of the input list for each intervention argument.")
      }

      if (is.character(unlist(ref_intervention_varlist))) {
        if (any(unlist(ref_intervention_varlist) %in% dta_cols) == F) {
          stop("The input treatment variable must be one of the columns in the data.")
        }
      }
      
      if (any(!is.na(unlist(compModels)))) {
        if ((any(str_detect(as.character(substitute(estimator)), "strat")) & (hazard == F) & (ref_total_effect == F)) | (any(str_detect(as.character(substitute(estimator)), "weight")) & (hazard == F) & (ref_total_effect == F))) {
          warning("The competing model is used for nonparametric risk estimation for direct effect case in stratified ICE. The keyword argument competing model statments are ignored.")
        } 
      }

      ref <- ice_strat(data = data, K = K, id = id, time_name = time_name, outcome_name = outcome_name,
                      censor_name = censor_name, competing_name = competing_name, total_effect = ref_total_effect,
                      outcome_model = ref_outcome_model, censor_model = censor_model, 
                      competing_model = ref_competing_model,
                      hazard_model = hazard_model,
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
    
    this_outcome_init <- list(ref$outcome_init)
    # this_comp_init <- list(ref$comp_init)
    this_np_model <- list(ref$np_model)
    
    names(this_outcome_init) <- names(this_np_model) <- ref_description
    
    this_outcome_by_step <- get_models(ref, "outcome_by_step", ref_description)
    this_hazard_by_step <- get_models(ref, "hazard_by_step", ref_description)
    this_comp_by_step <- get_models(ref, "comp_by_step", ref_description)
    
    outcome_by_step <- c(outcome_by_step, this_outcome_by_step)
    hazard_by_step <- c(hazard_by_step, this_hazard_by_step)
    comp_by_step <- c(comp_by_step, this_comp_by_step)
    
    outcome_init <- c(outcome_init, this_outcome_init)
    # comp_init <- c(comp_init, this_comp_init)
    np_model <- c(np_model, this_np_model)

    # this_model <- list(ref$fit_models)
    # names(this_model) <- ref_description
    # 
    # this_summary <- list(ref$model_summary)
    # names(this_summary) <- ref_description
    # 
    # this_stderr <- list(ref$model_stderr)
    # names(this_stderr) <- ref_description
    # 
    # this_vcov <- list(ref$model_vcov)
    # names(this_vcov) <- ref_description
    # 
    # this_rmse <- list(ref$model_rmse)
    # names(this_rmse) <- ref_description
    # 
    # fit_all <- c(fit_all, this_model)
    # fit_summary <- c(fit_summary, this_summary)
    # fit_stderr <- c(fit_stderr, this_stderr)
    # fit_vcov <- c(fit_vcov, this_vcov)
    # fit_rmse <- c(fit_rmse, this_rmse)

    risk_descript <- c(risk_descript, c(rep(str_to_title(ref_description), K+1), rep("Natural Course (nonparametric)", K+1)))
    risk_interv <- c(risk_interv, c(ref$gformula_risk[1, ], c(0, ref$weight_h$risk)))

    if (bootstrap) {
      intervention_descriptions <- c(intervention_descriptions, list("Natural Course"))
    }

    ninterv_new <- length(intervention_descriptions)

    ice_critical_value <- rr_critical_value <- rd_critical_value  <- matrix(NA, nrow = 2, ncol = ninterv + 1)

    if (bootstrap) {
      outcome_init_boot <- comp_init_boot <- np_model_boot <- outcome_by_step_boot <- 
        comp_by_step_boot <- hazard_by_step_boot <- data_boot_all <- c()
    } else {
      outcome_init_boot <- comp_init_boot <- np_model_boot <- outcome_by_step_boot <- 
        comp_by_step_boot <- hazard_by_step_boot <- data_boot_all <- NULL
    }

    if (ninterv_new > 0) {
      
      critical_value_all_upper <- critical_value_all_lower <- se_all <- c()
      
      if (ref_idx != 0) {
        idx_list <- rev(1:ninterv_new)
      } else {
        idx_list <- 1:ninterv_new
      }

    for (int in idx_list) {
      
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
                            hazard_model = hazard_model,
                            interventions = this_interv, intervention_names = this_int_var,
                            hazard_based = hazard, weighted = weight, treat_model = this_treat_model,
                            obs_treatment_names = this_obs_treatment_varnames, intervention_description = this_descript,
                            intervention_times = this_time)
      
      this_outcome_init <- list(this_fit$outcome_init)
      # this_comp_init <- list(this_fit$comp_init)
      this_np_model <- list(this_fit$np_model)
      
      names(this_outcome_init) <- names(this_np_model) <- this_descript
      
      this_outcome_by_step <- get_models(this_fit, "outcome_by_step", this_descript)
      this_hazard_by_step <- get_models(this_fit, "hazard_by_step", this_descript)
      this_comp_by_step <- get_models(this_fit, "comp_by_step", this_descript)
      
      outcome_by_step <- c(outcome_by_step, this_outcome_by_step)
      hazard_by_step <- c(hazard_by_step, this_hazard_by_step)
      comp_by_step <- c(comp_by_step, this_comp_by_step)
      
      outcome_init <- c(outcome_init, this_outcome_init)
      # comp_init <- c(comp_init, this_comp_init)
      np_model <- c(np_model, this_np_model)
      
      # fit_all <- c(fit_all, this_model)
      # fit_summary <- c(fit_summary, this_summary)
      # fit_stderr <- c(fit_stderr, this_stderr)
      # fit_vcov <- c(fit_vcov, this_vcov)
      # fit_rmse <- c(fit_rmse, this_rmse)
      # outcome_init <- c(outcome_init, this_outcome_init)
      # comp_init <- c(comp_init, this_comp_init)
      # 
      # this_model <- list(this_fit$fit_models)
      # names(this_model) <- this_descript
      # 
      # this_summary <- list(this_fit$model_summary)
      # names(this_summary) <- this_descript
      # 
      # this_stderr <- list(this_fit$model_stderr)
      # names(this_stderr) <- this_descript
      # 
      # this_vcov <- list(this_fit$model_vcov)
      # names(this_vcov) <- this_descript
      # 
      # this_rmse <- list(this_fit$model_rmse)
      # names(this_rmse) <- this_descript
      # 
      # fit_all <- c(fit_all, this_model)
      # fit_summary <- c(fit_summary, this_summary)
      # fit_stderr <- c(fit_stderr, this_stderr)
      # fit_vcov <- c(fit_vcov, this_vcov)
      # fit_rmse <- c(fit_rmse, this_rmse)

      if (this_descript != "Natural Course") {

      summary[int+1, 2] <- this_fit$gformula_risk_last_time
      risk_descript <- c(risk_descript, rep(str_to_title(this_descript), K+1))
      risk_interv <- c(risk_interv, this_fit$gformula_risk[1, ])
      }

      if (bootstrap) {
        
        if (ref_idx == 0) {
          ref_time <- NULL
        } else {
          ref_time <- ref_int_times
        }
        
        # if this intervention is the reference intervention
        if (ref_idx == int & int != idx_list[1] | ref_idx == 0 & this_descript == "Natural Course") {
          
          
          this_boot <- list(ice_se = ref_first_boot$ref_se, 
                            ref_se = ref_first_boot$ref_se,
                            rr_se = 1,
                            rd_se = 0, 
                            ice_cv_all_upper = ref_first_boot$ref_cv_all_upper, 
                            ice_cv_all_lower = ref_first_boot$ref_cv_all_lower, 
                            # ice_cv_upper = tail(ref_first_boot$ref_cv_all_upper, 1), 
                            # ice_cv_lower = tail(ref_first_boot$ref_cv_all_lower, 1), 
                            ref_cv_all_upper = ref_first_boot$ref_cv_all_upper, 
                            ref_cv_all_lower = ref_first_boot$ref_cv_all_lower, 
                            # ref_cv_upper = tail(ref_first_boot$ref_cv_all_upper, 1), 
                            # ref_cv_lower = tail(ref_first_boot$ref_cv_all_lower, 1), 
                            ref_ipw_se = ref_first_boot$ref_ipw_se, 
                            ref_ipw_cv_all_upper = ref_first_boot$ref_ipw_cv_all_upper,
                            ref_ipw_cv_all_lower = ref_first_boot$ref_ipw_cv_all_lower,
                            rr_cv_upper = 1, 
                            rr_cv_lower = 1, 
                            rd_cv_upper = 0, 
                            rd_cv_lower = 0, 
                            outcome_init = ref_first_boot$ref_outcome_init, 
                            comp_init = ref_first_boot$ref_comp_init, 
                            np_model = ref_first_boot$np_model, 
                            outcome_by_step = ref_first_boot$ref_outcome_by_step, 
                            comp_by_step = ref_first_boot$ref_comp_by_step, 
                            hazard_by_step = ref_first_boot$ref_hazard_by_step)
          
          give_warning(ref_first_boot$ref_data_err, ref_first_boot$ref_data_err_mssg_combine)
          give_warning(ref_first_boot$ref_model_err, ref_first_boot$ref_model_err_mssg_combine)
          
        } else {

        this_boot <- bootstrap_ice(ice_strat, K, nboot, significance_level, parallel, ncores, ref_description,
                                   ref_intervention_varlist, this_total_effect, boot_interv,
                                   this_interv, this_int_var, this_descript, this_time, ref_time,
                                   data, id, set_seed,
                                   time_name = time_name, outcome_name = outcome_name,
                                   censor_name = censor_name, competing_name = competing_name,
                                   hazard_based = hazard, weighted = weight, treat_model = this_treat_model,
                                   obs_treatment_names = this_obs_treatment_varnames,
                                   outcome_model = outcome_model, censor_model = censor_model,
                                   competing_model = competing_model, hazard_model = hazard_model)
        }
        
        # record the reference intervention bootstrap
        if (int == idx_list[1]) {
          
          ref_first_boot <- this_boot
          
        } 
        
        

        this_se <- tail(this_boot$ice_se, 1)
        this_rr_se <- this_boot$rr_se
        this_rd_se <- this_boot$rd_se

        this_cv_upper <- tail(this_boot$ice_cv_all_upper, 1)
        this_cv_lower <- tail(this_boot$ice_cv_all_lower, 1)
        this_rr_cv_upper <- this_boot$rr_cv_upper
        this_rr_cv_lower <- this_boot$rr_cv_lower
        this_rd_cv_upper <- this_boot$rd_cv_upper
        this_rd_cv_lower <- this_boot$rd_cv_lower
        
        this_outcome_init_boot <- list(this_boot$outcome_init)
        this_comp_init_boot <- list(this_boot$comp_init)
        this_np_model_boot <- list(this_boot$np_model)
        this_outcome_by_step_boot <- list(this_boot$outcome_by_step)
        this_comp_by_step_boot <- list(this_boot$comp_by_step)
        this_hazard_by_step_boot <- list(this_boot$hazard_by_step)
        
        names(this_outcome_init_boot) <- names(this_np_model_boot) <- names(this_comp_init_boot) <- 
          names(this_outcome_by_step_boot) <- names(this_comp_by_step_boot) <- names(this_hazard_by_step_boot) <- this_descript
        
        outcome_init_boot <- c(outcome_init_boot, this_outcome_init_boot)
        comp_init_boot <- c(comp_init_boot, this_comp_init_boot)
        np_model_boot <- c(np_model_boot, this_np_model_boot)
        outcome_by_step_boot <- c(outcome_by_step_boot, this_outcome_by_step_boot)
        comp_by_step_boot <- c(comp_by_step_boot, this_comp_by_step_boot)
        hazard_by_step_boot <- c(hazard_by_step_boot, this_hazard_by_step_boot)

        # this_model_boot <- list(this_boot$boot_models)
        # names(this_model_boot) <- this_descript
        # 
        # this_summary_boot <- list(this_boot$boot_summary)
        # names(this_summary_boot) <- this_descript
        # 
        # this_stderr_boot <- list(this_boot$boot_stderr)
        # names(this_stderr_boot) <- this_descript
        # 
        # this_vcov_boot <- list(this_boot$boot_vcov)
        # names(this_vcov_boot) <- this_descript
        # 
        # this_rmse_boot <- list(this_boot$boot_rmse)
        # names(this_rmse_boot) <- this_descript

        this_data_boot <- list(this_boot$boot_data)
        names(this_data_boot) <- this_descript

        # fit_all_boot <- c(fit_all_boot, this_model_boot)
        # fit_summary_boot <- c(fit_summary_boot, this_summary_boot)
        # fit_stderr_boot <- c(fit_stderr_boot, this_stderr_boot)
        # fit_vcov_boot <- c(fit_vcov_boot, this_vcov_boot)
        # fit_rmse_boot <- c(fit_rmse_boot, this_rmse_boot)
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
            critical_value_all_upper <- append_list(this_boot, "Natural Course (nonparametric)", critical_value_all_upper, "ref_ipw_cv_all_upper")
            critical_value_all_lower <- append_list(this_boot, "Natural Course (nonparametric)", critical_value_all_lower, "ref_ipw_cv_all_lower")
            se_all <- append_list(this_boot, "Natural Course (nonparametric)", se_all, "ref_ipw_se")
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
            critical_value_all_upper <- append_list(this_boot, "Natural Course (nonparametric)", critical_value_all_upper, "ref_ipw_cv_all_upper")
            critical_value_all_lower <- append_list(this_boot, "Natural Course (nonparametric)", critical_value_all_lower, "ref_ipw_cv_all_lower")
            se_all <- append_list(this_boot, "Natural Course (nonparametric)", se_all, "ref_ipw_se")
          }
        }
      }
      if (bootstrap) {
        summary[ref_idx + 1, 5] <- this_boot$ref_se[K+1]

        if (ref_idx == 0) {
          ice_critical_value[1, 1] <- tail(this_boot$ref_cv_all_lower, 1)
          ice_critical_value[2, 1] <- tail(this_boot$ref_cv_all_upper, 1)
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
    summary[-(ref_idx + 1), 10] <- summary[-(ref_idx + 1), 3] - rr_critical_value_lower * summary[-(ref_idx + 1), 6]
    summary[-(ref_idx + 1), 11] <- summary[2:(ninterv + 1), 3] + rr_critical_value_upper * summary[-(ref_idx + 1), 6]
    summary[-(ref_idx + 1), 12] <- summary[2:(ninterv + 1), 4] - rd_critical_value_lower * summary[-(ref_idx + 1), 7]
    summary[-(ref_idx + 1), 13] <- summary[2:(ninterv + 1), 4] + rd_critical_value_upper * summary[-(ref_idx + 1), 7]
    
    if (normal_quantile) {
    risk_time$Critical_Value_Lower <- ice_critical_value_lower
    risk_time$Critical_Value_Upper <- ice_critical_value_upper
    }

  }

  if (!is.null(compevent_name)) {
  if (any(unlist(total_effect))) {
    rownames(summary) <- paste0(rownames(summary), " (", "Total Effect" , ")")
    risk_time$Intervention <- paste0(risk_time$Intervention, " (", "Total Effect" , ")")
  } else {
    rownames(summary) <- paste0(rownames(summary), " (", "Direct Effect" , ")")
    risk_time$Intervention <- paste0(risk_time$Intervention, " (", "Direct Effect" , ")")
  }
  }
  summary_all <- rbind(summary_all, summary)
  risk_time_all <- rbind(risk_time_all, risk_time)
  }

  summary_all <- round(summary_all, 5)

  ## match back the time points using original specification
  risk_time_all$originTime <- c(0, unique_time_names)[as.integer(risk_time_all$Time) + 1]

  risk_time_all$originTime[which(is.na(risk_time_all$originTime))] <- 0


  return(list(estimator.type = ice_colname, 
              summary = summary_all, risk.over.time = risk_time_all,
              initial.outcome = outcome_init, initial.comp = comp_init, 
              np.risk.model = np_model, outcome.models.by.step = outcome_by_step, 
              comp.models.by.step = comp_by_step, 
              hazard.models.by.step = hazard_by_step, 
              boot.data = data_boot_all,
              boot.initial.outcome = outcome_init_boot,
              boot.initial.comp = comp_init_boot,
              boot.np.risk.model = np_model_boot,
              boot.outcome.models.by.step = outcome_by_step_boot,
              boot.comp.models.by.step = comp_by_step_boot,
              boot.hazard.models.by.step = hazard_by_step_boot)
              # boot.models = fit_all_boot, boot.summary = fit_summary_boot,
              # boot.stderr = fit_stderr_boot, boot.vcov = fit_vcov_boot,
              # boot.rmse = fit_rmse, estimator.type = ice_colname)
  )

}

#' Split Keyword Arguments
#'
#' @param argument the keyword argument list.
#' @param target_string a string specifying the target keyword argument pattern.
#'
#' @return a list containing the information of the targeted keyword argument.
#'
#' @internal
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
#' @param original_model the global original model statements specified through regular argument input.
#'
#' @return a list containing the model inputs corresponding to the interventions.
#'
#' @internal
get_model_formula <- function(interv_list, arg_interv, model_interv, ninterv, arg_str, original_model) {
  models <- as.list(rep(NA, ninterv))
  
  interv_list <- unique(interv_list)
  
  for (i in 1:ninterv) {

    interv_i <- interv_list[i]
    if (interv_i %in% model_interv) {
      arg_name_i <- paste0(arg_str, ".", interv_i)
      outcome_formula <- arg_interv[[which(names(arg_interv) == arg_name_i)]]
      models[[i]] <- outcome_formula
    } else {
      if (!is.null(original_model)) {
        models[[i]] <- original_model
      }
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
#'
#' @internal
construct_interv_value <- function(data, timevar, int_value, int_time, int_var){

  nint <- length(int_var)
  # print(int_var)

  # print(int_value)

  new_int_value <- list()

  for (i in 1:nint) {
    this_int_value <- int_value[[i]]
    not_int_idx <- which(!data[, timevar] %in% int_time[[i]])
    
    if (length(not_int_idx) > 0) {
    not_int_value <- data[, int_var[[i]]][not_int_idx]
    this_int_value[not_int_idx] <- not_int_value
    }
    new_int_value <- c(new_int_value, list(this_int_value))
  }


  return(new_int_value)
}


#' Match bootstrap values to aggregate data frame
#'
#' @internal
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
#' @internal
append_list <- function(item, name, append_list, item_name) {
  item_list <- list(item[[item_name]])
  names(item_list) <- name
  append_list <- c(append_list, item_list)
  
  return(append_list)
}

