#' Estimate the emulated trial effect
#'
#' @param clones_weighted A named list of weighted clone data frames, or a
#'   single data frame containing an arm column.
#' @param method Analysis method: `"Cox"`, `"logistic"`, or `"KM"`.
#' @param cluster Column name used for robust clustering in Cox models.
#' @param weights Optional weight column name. Unquoted column names are also
#'   accepted when they exist in the data.
#' @param predictors Optional adjustment predictors. For KM, predictors define
#'   additional strata rather than covariate adjustment.
#' @param outcome Column name for the outcome indicator.
#' @param time_start Column name for interval start time.
#' @param time_stop Column name for interval stop time.
#' @param arm Column name for treatment arm after binding clone lists.
#'
#' @returns A fitted model object: `"coxph"` for Cox, `"glm"` for logistic, or
#'   `"survfit"` for KM.
#' @export
#' @examples
#' data(lungcancer)
#' arms <- c("Control", "Surgery")
#' clones <- clone_arms(lungcancer, arms)
#' policies <- create_policy_A(
#'   arms, "surgery", "timetosurgery", 182.62, "death", "fup_obs",
#'   clone_outcome = "outcome", clone_followup = "fup"
#' )
#' clones_policy <- apply_logics(clones, policies)
#' censoring_logics <- create_censoring_logics_A(
#'   arms, "surgery", "timetosurgery", 182.62, "fup_obs",
#'   clone_censoring = "censoring",
#'   clone_uncensored_followup = "fup_uncensored"
#' )
#' clones_censored <- apply_logics(clones_policy, censoring_logics)
#' clones_final <- create_final_data(
#'   clones_censored, "fup", "outcome", "censoring", "id"
#' )
#' clones_estimated <- estimate_censoring(
#'   clones_final,
#'   predictors = c("age", "sex"),
#'   method = "pooled_logit"
#' )
#' clones_weighted <- weight_cases(clones_estimated)
#' fit <- emul_estimate(
#'   clones_weighted,
#'   method = "Cox",
#'   weights = "weight_Cox",
#'   predictors = c("age", "sex")
#' )
emul_estimate <- function(
  clones_weighted,
  method = c("Cox", "logistic", "KM"),
  cluster = "id",
  weights = NULL,
  predictors = NULL,
  outcome = "outcome",
  time_start = "Tstart",
  time_stop = "Tstop",
  arm = "arms"
) {
  weights_expr <- substitute(weights)
  method <- match.arg(method)
  dat <- bind_clone_arms(clones_weighted, arm = arm)

  if (length(stats::na.omit(levels(dat[[arm]]))) < 2L) {
    stop("`", arm, "` must contain at least two levels.", call. = FALSE)
  }

  predictors <- normalize_predictors(predictors, exclude = arm)
  required_columns <- c(arm, outcome, predictors)
  if (method %in% c("Cox", "KM")) {
    required_columns <- c(required_columns, time_start, time_stop)
  }
  if (method == "Cox" && !is.null(cluster)) {
    required_columns <- c(required_columns, cluster)
  }
  .assert_required_columns(dat, required_columns)

  weight_col <- normalize_weights(weights, weights_expr, dat)
  weighted_data <- add_analysis_weights(dat, weight_col)
  dat <- weighted_data$data
  analysis_weight <- weighted_data$weight

  surv_response <- paste0(
    "survival::Surv(",
    backtick_name(time_start),
    ", ",
    backtick_name(time_stop),
    ", ",
    backtick_name(outcome),
    ")"
  )
  cox_formula <- emul_formula(
    surv_response,
    predictors = predictors,
    cluster = cluster,
    arm = arm
  )
  km_formula <- emul_formula(
    surv_response,
    predictors = predictors,
    arm = arm
  )
  logistic_formula <- emul_formula(
    backtick_name(outcome),
    predictors = predictors,
    arm = arm
  )

  if (method == "KM" && length(predictors) > 0L) {
    message(
      "`predictors` in KM create separate strata; they do not produce ",
      "covariate-adjusted survival curves."
    )
  }

  if (method == "Cox") {
    cox_args <- list(
      formula = cox_formula,
      data = quote(dat),
      robust = TRUE,
      ties = "efron"
    )
    if (!is.null(analysis_weight)) {
      cox_args$weights <- as.name(analysis_weight)
    }
    return(base::do.call(survival::coxph, cox_args))
  }

  if (method == "logistic") {
    logistic_args <- list(
      formula = logistic_formula,
      data = quote(dat),
      family = stats::binomial(link = "logit")
    )
    if (!is.null(analysis_weight)) {
      logistic_args$weights <- as.name(analysis_weight)
    }
    return(base::do.call(stats::glm, logistic_args))
  }

  survfit_args <- list(
    formula = km_formula,
    data = quote(dat)
  )
  if (!is.null(analysis_weight)) {
    survfit_args$weights <- as.name(analysis_weight)
  }
  base::do.call(survival::survfit, survfit_args)
}
