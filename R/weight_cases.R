#' Add inverse probability of censoring weights
#'
#' @param clones A named list of clone data frames returned by
#'   `estimate_censoring()`.
#' @param uncensored_prob Column name for denominator uncensoring probability.
#' @param numerator_uncensored_prob Column name for numerator uncensoring
#'   probability. When this column exists, stabilized weights are created.
#' @param weight Column name for the output weight.
#' @param eps Small probability floor to avoid division by zero.
#'
#' @returns A named list of clone data frames with the weight column added.
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
weight_cases <- function(
  clones,
  uncensored_prob = "P_uncens",
  numerator_uncensored_prob = "P_uncens_num",
  weight = "weight_Cox",
  eps = 1e-6
) {
  .assert_clone_columns(clones, uncensored_prob)
  .assert_probability_floor(eps)

  arms <- names(clones)
  res <- vector("list", length = length(arms))
  names(res) <- arms

  for (arm in arms) {
    dat <- clones[[arm]]
    if (!is.numeric(dat[[uncensored_prob]])) {
      stop("`", uncensored_prob, "` must be numeric in clone arm `", arm, "`.", call. = FALSE)
    }

    denominator <- pmax(dat[[uncensored_prob]], eps)
    if (numerator_uncensored_prob %in% names(dat)) {
      if (!is.numeric(dat[[numerator_uncensored_prob]])) {
        stop(
          "`", numerator_uncensored_prob, "` must be numeric in clone arm `", arm, "`.",
          call. = FALSE
        )
      }
      weights <- dat[[numerator_uncensored_prob]] / denominator
    } else {
      weights <- 1 / denominator
    }

    res[[arm]] <-
      dat |>
      dplyr::mutate(!!weight := .env[["weights"]])
  }

  res
}
