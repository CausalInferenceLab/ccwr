#' Bootstrap a complete clone-censor-weight analysis
#'
#' Resamples subjects from the original subject-level data and repeats the
#' complete analysis in every bootstrap replicate: cloning, strategy-specific
#' artificial censoring, person-time expansion, censoring-model estimation,
#' weighting, and outcome-model estimation.
#'
#' @param data A data frame with one row per subject.
#' @param arms A character vector of length two. The first element names the
#'   strategy of not receiving treatment during the grace period; the second
#'   names the strategy of receiving treatment during the grace period.
#' @param id Name of the column that uniquely identifies subjects.
#' @param treatment Name of the binary treatment column.
#' @param time_to_treatment Name of the numeric time-to-treatment column.
#' @param grace_period Positive numeric length of the treatment grace period.
#' @param outcome Name of the binary outcome column.
#' @param followup Name of the numeric follow-up-time column.
#' @param censoring_predictors Optional character vector of predictors for the
#'   denominator censoring model.
#' @param censoring_method Censoring model passed to [estimate_censoring()]:
#'   `"pooled_logit"`, `"stabilized_logit"`, or `"Cox"`.
#' @param numerator_predictors Optional character vector of predictors for the
#'   numerator model when `censoring_method = "stabilized_logit"`.
#' @param method Outcome analysis method: `"Cox"` or `"logistic"`.
#' @param predictors Optional adjustment predictors for the outcome model.
#' @param n_bootstrap Number of bootstrap resamples.
#' @param conf_level Confidence level for the percentile interval.
#' @param eps Probability floor passed to [estimate_censoring()] and
#'   [weight_cases()].
#' @param seed Optional random seed.
#'
#' @returns A list containing the full-sample effect estimate, percentile
#'   confidence limits, bootstrap estimates, and analysis settings. Effects are
#'   returned on the exponentiated coefficient scale: a hazard ratio for Cox
#'   models and an odds ratio for logistic models.
#' @export
#' @examples
#' data(lungcancer)
#' boot <- emul_estimate_bootstrap(
#'   lungcancer,
#'   arms = c("Control", "Surgery"),
#'   id = "id",
#'   treatment = "surgery",
#'   time_to_treatment = "timetosurgery",
#'   grace_period = 182.62,
#'   outcome = "death",
#'   followup = "fup_obs",
#'   censoring_predictors = c("age", "sex"),
#'   predictors = c("age", "sex"),
#'   n_bootstrap = 3,
#'   seed = 1
#' )
emul_estimate_bootstrap <- function(
  data,
  arms,
  id,
  treatment,
  time_to_treatment,
  grace_period,
  outcome,
  followup,
  censoring_predictors = NULL,
  censoring_method = c("pooled_logit", "stabilized_logit", "Cox"),
  numerator_predictors = NULL,
  method = c("Cox", "logistic"),
  predictors = NULL,
  n_bootstrap = 200,
  conf_level = 0.95,
  eps = 1e-6,
  seed = NULL
) {
  censoring_method <- match.arg(censoring_method)
  method <- match.arg(method)
  censoring_predictors <- normalize_predictors(censoring_predictors)
  if (!is.null(numerator_predictors)) {
    numerator_predictors <- normalize_predictors(numerator_predictors)
  }
  predictors <- normalize_predictors(predictors)

  validate_bootstrap_data(
    data = data,
    arms = arms,
    id = id,
    treatment = treatment,
    time_to_treatment = time_to_treatment,
    grace_period = grace_period,
    outcome = outcome,
    followup = followup,
    censoring_predictors = censoring_predictors,
    numerator_predictors = numerator_predictors,
    predictors = predictors
  )
  .assert_positive_integer(n_bootstrap, "n_bootstrap")
  .assert_conf_level(conf_level)
  .assert_probability_floor(eps)
  if (
    !is.null(seed) &&
      (!is.numeric(seed) || length(seed) != 1L || is.na(seed))
  ) {
    stop("`seed` must be NULL or a single non-missing number.", call. = FALSE)
  }

  full_fit <- run_ccw_analysis(
    data = data,
    arms = arms,
    subject_id = id,
    treatment = treatment,
    time_to_treatment = time_to_treatment,
    grace_period = grace_period,
    outcome = outcome,
    followup = followup,
    censoring_predictors = censoring_predictors,
    censoring_method = censoring_method,
    numerator_predictors = numerator_predictors,
    method = method,
    predictors = predictors,
    eps = eps
  )
  estimate <- extract_ccw_effect(full_fit, arms[[2]])

  if (!is.null(seed)) {
    set.seed(seed)
  }

  boot_estimates <- numeric(n_bootstrap)
  for (i in seq_len(n_bootstrap)) {
    sampled_rows <- sample.int(
      nrow(data),
      size = nrow(data),
      replace = TRUE
    )
    bootstrap_data <- data[sampled_rows, , drop = FALSE]
    bootstrap_data$.bootstrap_id <- seq_len(nrow(bootstrap_data))

    bootstrap_fit <- run_ccw_analysis(
      data = bootstrap_data,
      arms = arms,
      subject_id = ".bootstrap_id",
      treatment = treatment,
      time_to_treatment = time_to_treatment,
      grace_period = grace_period,
      outcome = outcome,
      followup = followup,
      censoring_predictors = censoring_predictors,
      censoring_method = censoring_method,
      numerator_predictors = numerator_predictors,
      method = method,
      predictors = predictors,
      eps = eps
    )
    boot_estimates[[i]] <- extract_ccw_effect(
      bootstrap_fit,
      arms[[2]]
    )
  }

  alpha <- 1 - conf_level
  interval <- stats::quantile(
    boot_estimates,
    probs = c(alpha / 2, 1 - alpha / 2),
    names = FALSE
  )

  list(
    estimate = estimate,
    ci_lower = interval[[1]],
    ci_upper = interval[[2]],
    estimates = boot_estimates,
    n_bootstrap = n_bootstrap,
    conf_level = conf_level,
    method = method,
    censoring_method = censoring_method
  )
}

#' Run one complete clone-censor-weight analysis
#'
#' @noRd
run_ccw_analysis <- function(
  data,
  arms,
  subject_id,
  treatment,
  time_to_treatment,
  grace_period,
  outcome,
  followup,
  censoring_predictors,
  censoring_method,
  numerator_predictors,
  method,
  predictors,
  eps
) {
  internal_columns <- c(
    ".ccw_outcome",
    ".ccw_followup",
    ".ccw_censoring",
    ".ccw_uncensored_followup",
    ".ccw_interval_id",
    ".ccw_tstart",
    ".ccw_tstop",
    ".ccw_weight",
    ".ccw_arm"
  )
  conflicting_columns <- intersect(internal_columns, names(data))
  if (length(conflicting_columns) > 0L) {
    stop(
      "Input data contains reserved columns: ",
      paste(conflicting_columns, collapse = ", "),
      call. = FALSE
    )
  }

  clones <- clone_arms(data, arms)
  policies <- create_policy_A(
    arms = arms,
    treatment = treatment,
    time_to_treatment = time_to_treatment,
    grace_period = grace_period,
    outcome = outcome,
    followup = followup,
    clone_outcome = ".ccw_outcome",
    clone_followup = ".ccw_followup"
  )
  clones <- apply_logics(clones, policies)

  censoring_logics <- create_censoring_logics_A(
    arms = arms,
    treatment = treatment,
    time_to_treatment = time_to_treatment,
    grace_period = grace_period,
    followup = followup,
    clone_censoring = ".ccw_censoring",
    clone_uncensored_followup = ".ccw_uncensored_followup"
  )
  clones <- apply_logics(clones, censoring_logics)

  clones <- create_final_data(
    clones = clones,
    clone_followup = ".ccw_followup",
    clone_outcome = ".ccw_outcome",
    clone_censoring = ".ccw_censoring",
    col_ids = subject_id,
    timestamp_start = ".ccw_tstart",
    id = ".ccw_interval_id",
    timestamp_stop = ".ccw_tstop"
  )
  clones <- estimate_censoring(
    clones = clones,
    predictors = censoring_predictors,
    method = censoring_method,
    numerator_predictors = numerator_predictors,
    censoring = ".ccw_censoring",
    id = subject_id,
    time_start = ".ccw_tstart",
    time_stop = ".ccw_tstop",
    eps = eps
  )
  clones <- weight_cases(
    clones = clones,
    uncensored_prob = "P_uncens",
    numerator_uncensored_prob = "P_uncens_num",
    weight = ".ccw_weight",
    eps = eps
  )

  emul_estimate(
    clones_weighted = clones,
    method = method,
    cluster = subject_id,
    weights = ".ccw_weight",
    predictors = predictors,
    outcome = ".ccw_outcome",
    time_start = ".ccw_tstart",
    time_stop = ".ccw_tstop",
    arm = ".ccw_arm"
  )
}

#' Extract the exponentiated treatment-arm coefficient
#'
#' @noRd
extract_ccw_effect <- function(fit, treated_arm) {
  coefficient <- find_arm_coefficient(
    stats::coef(fit),
    arm = ".ccw_arm",
    arm_level = treated_arm
  )
  unname(exp(stats::coef(fit)[[coefficient]]))
}

#' Validate original subject-level bootstrap data
#'
#' @noRd
validate_bootstrap_data <- function(
  data,
  arms,
  id,
  treatment,
  time_to_treatment,
  grace_period,
  outcome,
  followup,
  censoring_predictors,
  numerator_predictors,
  predictors
) {
  .assert_data_frame(data)
  if (nrow(data) < 2L) {
    stop("`data` must contain at least two subjects.", call. = FALSE)
  }
  if (
    !is.character(arms) ||
      length(arms) != 2L ||
      any(is.na(arms)) ||
      any(!nzchar(arms)) ||
      anyDuplicated(arms)
  ) {
    stop("`arms` must contain two distinct, non-empty names.", call. = FALSE)
  }

  column_arguments <- list(
    id = id,
    treatment = treatment,
    time_to_treatment = time_to_treatment,
    outcome = outcome,
    followup = followup
  )
  for (argument in names(column_arguments)) {
    value <- column_arguments[[argument]]
    if (
      !is.character(value) ||
        length(value) != 1L ||
        is.na(value) ||
        !nzchar(value)
    ) {
      stop(
        "`", argument, "` must be a single non-empty column name.",
        call. = FALSE
      )
    }
  }

  required_columns <- unique(c(
    id,
    treatment,
    time_to_treatment,
    outcome,
    followup,
    censoring_predictors,
    numerator_predictors,
    predictors
  ))
  .assert_required_columns(data, required_columns)
  if (".bootstrap_id" %in% names(data)) {
    stop("Input data contains reserved column: .bootstrap_id", call. = FALSE)
  }

  if (anyNA(data[[id]]) || anyDuplicated(data[[id]])) {
    stop("`id` must uniquely identify every subject.", call. = FALSE)
  }
  if (
    anyNA(data[[treatment]]) ||
      !all(data[[treatment]] %in% c(0, 1))
  ) {
    stop("`treatment` must contain only non-missing 0/1 values.", call. = FALSE)
  }
  if (anyNA(data[[outcome]]) || !all(data[[outcome]] %in% c(0, 1))) {
    stop("`outcome` must contain only non-missing 0/1 values.", call. = FALSE)
  }
  if (
    !is.numeric(data[[followup]]) ||
      anyNA(data[[followup]]) ||
      any(!is.finite(data[[followup]])) ||
      any(data[[followup]] < 0)
  ) {
    stop("`followup` must contain finite, non-negative times.", call. = FALSE)
  }
  if (
    !is.numeric(data[[time_to_treatment]]) ||
      any(
        !is.na(data[[time_to_treatment]]) &
          (
            !is.finite(data[[time_to_treatment]]) |
              data[[time_to_treatment]] < 0
          )
      )
  ) {
    stop(
      "`time_to_treatment` must contain non-negative times or NA.",
      call. = FALSE
    )
  }
  if (
    !is.numeric(grace_period) ||
      length(grace_period) != 1L ||
      is.na(grace_period) ||
      !is.finite(grace_period) ||
      grace_period <= 0
  ) {
    stop("`grace_period` must be a positive finite number.", call. = FALSE)
  }

  invisible(data)
}
