#' Bootstrap the emulated trial effect
#'
#' @param clones_weighted A named list of weighted clone data frames, or a
#'   single data frame containing an arm column.
#' @param method Analysis method: `"Cox"` or `"logistic"`.
#' @param cluster Column name identifying resampling clusters.
#' @param predictors Optional adjustment predictors.
#' @param weights Optional weight column name. Unquoted column names are also
#'   accepted when they exist in the data.
#' @param n_bootstrap Number of bootstrap resamples.
#' @param outcome Column name for the outcome indicator.
#' @param time_start Column name for interval start time.
#' @param time_stop Column name for interval stop time.
#' @param arm Column name for treatment arm after binding clone lists.
#' @param conf_level Confidence level for the percentile interval.
#' @param seed Optional random seed.
#'
#' @returns A list with `ci_lower`, `ci_upper`, and bootstrap `estimates`.
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
#' boot <- emul_estimate_bootstrap(
#'   clones_weighted,
#'   method = "Cox",
#'   weights = "weight_Cox",
#'   predictors = c("age", "sex"),
#'   n_bootstrap = 3,
#'   seed = 1
#' )
emul_estimate_bootstrap <- function(
  clones_weighted,
  method = c("Cox", "logistic"),
  cluster = "id",
  predictors = NULL,
  weights = NULL,
  n_bootstrap = 200,
  outcome = "outcome",
  time_start = "Tstart",
  time_stop = "Tstop",
  arm = "arms",
  conf_level = 0.95,
  seed = NULL
) {
  weights_expr <- substitute(weights)
  method <- match.arg(method)
  .assert_positive_integer(n_bootstrap, "n_bootstrap")
  .assert_conf_level(conf_level)

  dat <- bind_clone_arms(clones_weighted, arm = arm)
  arm_levels <- levels(dat[[arm]])
  if (length(arm_levels) != 2L) {
    stop("Bootstrap effect estimation currently requires exactly two arms.", call. = FALSE)
  }

  predictors <- normalize_predictors(predictors, exclude = arm)
  required_columns <- c(cluster, arm, outcome, predictors)
  if (method == "Cox") {
    required_columns <- c(required_columns, time_start, time_stop)
  }
  .assert_required_columns(dat, required_columns)

  weight_col <- normalize_weights(weights, weights_expr, dat)

  if (!is.null(seed)) {
    set.seed(seed)
  }

  boot_estimates <- vector("numeric", length = n_bootstrap)
  unique_clusters <- unique(dat[[cluster]])

  for (i in seq_len(n_bootstrap)) {
    bootstrap_clusters <- unique_clusters[
      sample.int(
        length(unique_clusters),
        size = length(unique_clusters),
        replace = TRUE
      )
    ]
    bootstrap_sample <- dplyr::bind_rows(
      lapply(seq_along(bootstrap_clusters), function(j) {
        cluster_rows <- dat[
          dat[[cluster]] == bootstrap_clusters[[j]],
          ,
          drop = FALSE
        ]
        cluster_rows$.bootstrap_id <- j
        cluster_rows
      })
    )
    bootstrap_sample[[arm]] <- factor(bootstrap_sample[[arm]], levels = arm_levels)

    fit <- emul_estimate(
      bootstrap_sample,
      method = method,
      cluster = ".bootstrap_id",
      weights = weight_col,
      predictors = predictors,
      outcome = outcome,
      time_start = time_start,
      time_stop = time_stop,
      arm = arm
    )
    coef_values <- stats::coef(fit)
    arm_coef <- find_arm_coefficient(coef_values, arm, arm_levels[[2]])
    boot_estimates[[i]] <- unname(exp(coef_values[[arm_coef]]))
  }

  alpha <- 1 - conf_level
  lower_ci <- stats::quantile(boot_estimates, probs = alpha / 2)
  upper_ci <- stats::quantile(boot_estimates, probs = 1 - alpha / 2)

  list(
    ci_lower = lower_ci,
    ci_upper = upper_ci,
    estimates = boot_estimates
  )
}
