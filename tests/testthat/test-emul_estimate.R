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

test_that("emul_estimate_bootstrap returns percentile intervals", {
  result <- emul_estimate_bootstrap(
    make_weighted_clones(),
    method = "Cox",
    weights = "weight_Cox",
    n_bootstrap = 3,
    seed = 1
  )

  expect_named(result, c("ci_lower", "ci_upper", "estimates"))
  expect_length(result$estimates, 3)
  expect_true(all(is.finite(result$estimates)))
  expect_true(is.finite(result$ci_lower))
  expect_true(is.finite(result$ci_upper))
})

test_that("estimation helper functions remain internal", {
  exports <- getNamespaceExports("clonecensorweighting")

  expect_false("backtick_name" %in% exports)
  expect_false("normalize_predictors" %in% exports)
  expect_false("normalize_weights" %in% exports)
  expect_false("emul_formula" %in% exports)
})
