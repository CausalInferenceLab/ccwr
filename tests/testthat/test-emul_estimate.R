make_weighted_clones <- function() {
  list(
    Control = tibble::tibble(
      id = 1:8,
      Tstart = 0,
      Tstop = 1,
      outcome = c(1, 0, 1, 0, 1, 0, 0, 1),
      weight_Cox = 1,
      age = c(50, 55, 60, 65, 52, 57, 62, 67)
    ),
    Surgery = tibble::tibble(
      id = 1:8,
      Tstart = 0,
      Tstop = 1,
      outcome = c(0, 1, 0, 0, 0, 1, 0, 0),
      weight_Cox = 1,
      age = c(50, 55, 60, 65, 52, 57, 62, 67)
    )
  )
}

test_that("emul_estimate fits Cox, logistic, and KM analyses", {
  weighted <- make_weighted_clones()

  cox_fit <- emul_estimate(weighted, method = "Cox", weights = "weight_Cox")
  logistic_fit <- suppressWarnings(
    emul_estimate(weighted, method = "logistic", weights = weight_Cox)
  )
  km_fit <- emul_estimate(weighted, method = "KM", weights = "weight_Cox")

  expect_s3_class(cox_fit, "coxph")
  expect_s3_class(logistic_fit, "glm")
  expect_s3_class(km_fit, "survfit")
})

test_that("emul_estimate accepts data frame input with an arm column", {
  dat <- dplyr::bind_rows(make_weighted_clones(), .id = "arms")

  fit <- emul_estimate(dat, method = "Cox", weights = "weight_Cox")

  expect_s3_class(fit, "coxph")
})

test_that("emul_estimate_bootstrap repeats the complete analysis", {
  data(lungcancer)
  estimate_calls <- 0L
  original_estimate_censoring <- estimate_censoring
  local_mocked_bindings(
    estimate_censoring = function(...) {
      estimate_calls <<- estimate_calls + 1L
      original_estimate_censoring(...)
    },
    .package = "ccwr"
  )

  result <- emul_estimate_bootstrap(
    lungcancer,
    arms = c("Control", "Surgery"),
    id = "id",
    treatment = "surgery",
    time_to_treatment = "timetosurgery",
    grace_period = 182.62,
    outcome = "death",
    followup = "fup_obs",
    censoring_predictors = c("age", "sex"),
    method = "Cox",
    predictors = c("age", "sex"),
    n_bootstrap = 3,
    seed = 1
  )

  expect_equal(estimate_calls, 4L)
  expect_named(
    result,
    c(
      "estimate",
      "ci_lower",
      "ci_upper",
      "estimates",
      "n_bootstrap",
      "conf_level",
      "method",
      "censoring_method",
      "censoring_time_spline_df"
    )
  )
  expect_null(result$censoring_time_spline_df)
  expect_equal(result$estimate, 0.4905093388, tolerance = 1e-7)
  expect_equal(
    result$estimates,
    c(0.5056717858, 0.6770176964, 0.6370150551),
    tolerance = 1e-7
  )
  expect_length(result$estimates, 3)
  expect_true(all(is.finite(result$estimates)))
  expect_true(is.finite(result$ci_lower))
  expect_true(is.finite(result$ci_upper))
})

test_that("emul_estimate_bootstrap requires unique subject-level rows", {
  data(lungcancer)
  duplicated_subjects <- lungcancer
  duplicated_subjects$id[[2]] <- duplicated_subjects$id[[1]]

  expect_error(
    emul_estimate_bootstrap(
      duplicated_subjects,
      arms = c("Control", "Surgery"),
      id = "id",
      treatment = "surgery",
      time_to_treatment = "timetosurgery",
      grace_period = 182.62,
      outcome = "death",
      followup = "fup_obs",
      n_bootstrap = 2
    ),
    "uniquely identify"
  )
})

test_that("estimation helper functions remain internal", {
  exports <- getNamespaceExports("ccwr")

  expect_false("backtick_name" %in% exports)
  expect_false("normalize_predictors" %in% exports)
  expect_false("normalize_weights" %in% exports)
  expect_false("emul_formula" %in% exports)
})
