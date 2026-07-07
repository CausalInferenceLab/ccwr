#' Estimate censoring probabilities
#'
#' @param clones A named list of long-form clone data frames.
#' @param predictors Optional character vector of denominator model predictors.
#' @param method Censoring model. `"Cox"` fits a Cox censoring model,
#'   `"pooled_logit"` fits a pooled logistic denominator model, and
#'   `"stabilized_logit"` additionally fits a numerator model.
#' @param numerator_predictors Optional character vector of numerator model
#'   predictors for stabilized pooled-logit weights. When `NULL`, `predictors`
#'   are used.
#' @param censoring Column name for the censoring indicator.
#' @param id Column name for the subject identifier.
#' @param time_start Column name for interval start time.
#' @param time_stop Column name for interval stop time.
#' @param eps Small probability floor to avoid division by zero.
#'
#' @returns A named list of clone data frames with censoring probability
#'   columns added. All methods add `P_uncens`; pooled-logit methods also add
#'   `p_cens_den`, and stabilized pooled logit adds `p_cens_num` and
#'   `P_uncens_num`.
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
estimate_censoring <- function(
  clones,
  predictors = NULL,
  method = c("Cox", "pooled_logit", "stabilized_logit"),
  numerator_predictors = NULL,
  censoring = "censoring",
  id = "id",
  time_start = "Tstart",
  time_stop = "Tstop",
  eps = 1e-6
) {
  method <- match.arg(method)
  predictors <- normalize_predictors(predictors)

  if (method == "stabilized_logit" && is.null(numerator_predictors)) {
    numerator_predictors <- predictors
  }
  numerator_predictors <- normalize_predictors(numerator_predictors)

  required_columns <- c(id, time_start, time_stop, censoring, predictors)
  if (method == "stabilized_logit") {
    required_columns <- c(required_columns, numerator_predictors)
  }
  .assert_clone_columns(clones, required_columns)
  .assert_probability_floor(eps)

  arms <- names(clones)
  res <- vector("list", length = length(arms))
  names(res) <- arms

  cox_response <- paste0(
    "survival::Surv(",
    backtick_name(time_start),
    ", ",
    backtick_name(time_stop),
    ", ",
    backtick_name(censoring),
    ")"
  )
  cox_formula <- make_censoring_formula(
    cox_response,
    predictors
  )
  pooled_formula <- make_censoring_formula(
    backtick_name(censoring),
    predictors,
    time_var = time_start
  )

  for (arm in arms) {
    dat <- clones[[arm]]

    if (method == "Cox") {
      ms_cens <- survival::coxph(
        cox_formula,
        ties = "efron",
        data = dat
      )

      lin_pred <- stats::predict(
        ms_cens,
        newdata = dat,
        type = "lp",
        reference = "zero"
      )

      base_hazard <- dplyr::as_tibble(
        survival::basehaz(ms_cens, centered = FALSE)
      )
      names(base_hazard) <- c("hazard", "t")

      res[[arm]] <-
        dat |>
        dplyr::mutate(lin_pred = .env[["lin_pred"]]) |>
        dplyr::left_join(
          base_hazard,
          by = stats::setNames("t", time_start)
        ) |>
        dplyr::mutate(
          hazard = dplyr::coalesce(.data$hazard, 0),
          P_uncens = exp(-.data$hazard * exp(.data$lin_pred))
        )
    } else {
      fit_den <- stats::glm(
        pooled_formula,
        data = dat,
        family = stats::binomial(link = "logit")
      )
      p_cens_den <- .clamp_probability(
        stats::predict(fit_den, type = "response"),
        eps = eps
      )
      p_uncens_den <- cumulative_uncensoring(
        dat,
        p_cens_den,
        id = id,
        time_start = time_start,
        time_stop = time_stop,
        eps = eps
      )

      res[[arm]] <-
        dat |>
        dplyr::mutate(
          p_cens_den = .env[["p_cens_den"]],
          P_uncens = .env[["p_uncens_den"]]
        )

      if (method == "stabilized_logit") {
        numerator_data <- add_baseline_predictors(
          dat,
          predictors = numerator_predictors,
          id = id,
          time_start = time_start,
          time_stop = time_stop
        )
        numerator_formula <- make_censoring_formula(
          backtick_name(censoring),
          numerator_data$predictors,
          time_var = time_start
        )
        fit_num <- stats::glm(
          numerator_formula,
          data = numerator_data$data,
          family = stats::binomial(link = "logit")
        )
        p_cens_num <- .clamp_probability(
          stats::predict(fit_num, type = "response"),
          eps = eps
        )
        p_uncens_num <- cumulative_uncensoring(
          dat,
          p_cens_num,
          id = id,
          time_start = time_start,
          time_stop = time_stop,
          eps = eps
        )

        res[[arm]] <-
          res[[arm]] |>
          dplyr::mutate(
            p_cens_num = .env[["p_cens_num"]],
            P_uncens_num = .env[["p_uncens_num"]]
          )
      }
    }
  }

  res
}

#' Create a censoring model formula
#'
#' @noRd
make_censoring_formula <- function(
  response,
  predictors = NULL,
  time_var = NULL
) {
  terms <- unique(c(
    normalize_predictors(time_var),
    normalize_predictors(predictors)
  ))
  terms <- backtick_name(terms)

  formula <- stats::reformulate(terms, response = response)
  environment(formula) <- parent.frame()
  formula
}

#' Estimate cumulative probability of remaining uncensored
#'
#' @noRd
cumulative_uncensoring <- function(
  data,
  p_censoring,
  id = "id",
  time_start = "Tstart",
  time_stop = "Tstop",
  eps = 1e-6
) {
  .assert_data_frame(data)
  .assert_required_columns(data, c(id, time_start, time_stop))
  .assert_probability_floor(eps)

  if (!is.numeric(p_censoring) || length(p_censoring) != nrow(data)) {
    stop("`p_censoring` must have one numeric value per row of `data`.", call. = FALSE)
  }

  p_uncensored_interval <- 1 - .clamp_probability(p_censoring, eps = eps)
  p_uncensored <- numeric(nrow(data))
  ordered_rows <- order(data[[id]], data[[time_start]], data[[time_stop]])
  rows_by_id <- split(ordered_rows, data[[id]][ordered_rows])

  for (rows in rows_by_id) {
    interval_prob <- p_uncensored_interval[rows]
    p_uncensored[rows] <- c(1, cumprod(interval_prob[-length(interval_prob)]))
  }

  pmax(p_uncensored, eps)
}

#' Add baseline values of predictors
#'
#' @noRd
add_baseline_predictors <- function(
  data,
  predictors = NULL,
  id = "id",
  time_start = "Tstart",
  time_stop = "Tstop",
  prefix = ".baseline_"
) {
  predictors <- normalize_predictors(predictors)
  if (length(predictors) == 0L) {
    return(list(data = data, predictors = character()))
  }

  .assert_data_frame(data)
  .assert_required_columns(data, c(id, time_start, time_stop, predictors))

  baseline_predictors <- paste0(prefix, predictors)
  conflicting_columns <- intersect(baseline_predictors, names(data))
  if (length(conflicting_columns) > 0L) {
    stop(
      "Baseline predictor columns already exist: ",
      paste(conflicting_columns, collapse = ", "),
      call. = FALSE
    )
  }

  ordered_rows <- order(data[[id]], data[[time_start]], data[[time_stop]])
  baseline_rows <- ordered_rows[!duplicated(data[[id]][ordered_rows])]
  baseline_values <- data[baseline_rows, c(id, predictors), drop = FALSE]
  names(baseline_values)[names(baseline_values) %in% predictors] <-
    baseline_predictors

  list(
    data = dplyr::left_join(data, baseline_values, by = id),
    predictors = baseline_predictors
  )
}
