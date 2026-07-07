make_censoring_clones <- function() {
  list(
    Control = tibble::tibble(
      id = rep(1:6, each = 2),
      Tstart = rep(c(0, 1), 6),
      Tstop = rep(c(1, 2), 6),
      outcome = c(0, 1, 0, 0, 0, 0, 0, 1, 0, 0, 0, 1),
      censoring = c(0, 1, 0, 0, 1, 0, 0, 1, 0, 0, 1, 0),
      age = rep(c(50, 60, 55, 65, 58, 62), each = 2)
    ),
    Surgery = tibble::tibble(
      id = rep(1:6, each = 2),
      Tstart = rep(c(0, 1), 6),
      Tstop = rep(c(1, 2), 6),
      outcome = c(0, 0, 0, 1, 0, 0, 0, 0, 0, 1, 0, 0),
      censoring = c(0, 0, 0, 1, 0, 0, 1, 0, 0, 0, 0, 1),
      age = rep(c(50, 60, 55, 65, 58, 62), each = 2)
    )
  )
}

test_that("estimate_censoring adds pooled-logit denominator probabilities", {
  result <- suppressWarnings(
    estimate_censoring(
      make_censoring_clones(),
      predictors = "age",
      method = "pooled_logit"
    )
  )

  expect_named(result, c("Control", "Surgery"))
  expect_true(all(c("p_cens_den", "P_uncens") %in% names(result$Control)))
  expect_true(all(result$Control$p_cens_den > 0))
  expect_true(all(result$Control$p_cens_den < 1))
  expect_true(all(result$Control$P_uncens > 0))
  expect_true(all(result$Control$P_uncens <= 1))
})

test_that("estimate_censoring supports stabilized pooled-logit weights", {
  result <- suppressWarnings(
    estimate_censoring(
      make_censoring_clones(),
      predictors = "age",
      method = "stabilized_logit"
    )
  )

  expect_true(all(c("p_cens_num", "P_uncens_num") %in% names(result$Surgery)))
  expect_true(all(result$Surgery$P_uncens_num > 0))
  expect_true(all(result$Surgery$P_uncens_num <= 1))
})

test_that("estimate_censoring supports Cox censoring models", {
  result <- suppressWarnings(
    estimate_censoring(
      make_censoring_clones(),
      predictors = "age",
      method = "Cox"
    )
  )

  expect_true(all(c("lin_pred", "hazard", "P_uncens") %in% names(result$Control)))
  expect_true(all(result$Control$P_uncens > 0))
  expect_true(all(result$Control$P_uncens <= 1))
})

test_that("weight_cases creates unstabilized and stabilized IPC weights", {
  unstabilized <- list(
    Control = tibble::tibble(P_uncens = c(1, 0.5)),
    Surgery = tibble::tibble(P_uncens = c(0.25, 0.8))
  )
  stabilized <- list(
    Control = tibble::tibble(P_uncens = c(1, 0.5), P_uncens_num = c(1, 0.75)),
    Surgery = tibble::tibble(P_uncens = c(0.25, 0.8), P_uncens_num = c(0.5, 0.4))
  )

  expect_equal(weight_cases(unstabilized)$Control$weight_Cox, c(1, 2))
  expect_equal(weight_cases(stabilized)$Surgery$weight_Cox, c(2, 0.5))
})

test_that("censoring helper functions remain internal", {
  exports <- getNamespaceExports("clonecensorweighting")

  expect_false("make_censoring_formula" %in% exports)
  expect_false("cumulative_uncensoring" %in% exports)
  expect_false("add_baseline_predictors" %in% exports)
})
