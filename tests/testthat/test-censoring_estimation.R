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

test_that("pooled-logit censoring uses a configurable natural time spline", {
  clones <- make_censoring_clones()
  result <- suppressWarnings(
    estimate_censoring(
      clones,
      predictors = "age",
      method = "pooled_logit",
      time_spline_df = 3
    )
  )
  expected_fit <- suppressWarnings(
    stats::glm(
      censoring ~ splines::ns(Tstart, df = 3) + age,
      data = clones$Control,
      family = stats::binomial(link = "logit")
    )
  )
  expected_probability <- .clamp_probability(
    stats::predict(expected_fit, type = "response")
  )

  expect_equal(result$Control$p_cens_den, expected_probability)
  expect_error(
    estimate_censoring(
      clones,
      method = "pooled_logit",
      time_spline_df = 1
    ),
    "at least 2"
  )
})

test_that("pooled-logit censoring retains a linear time option", {
  formula <- make_censoring_formula(
    "censoring",
    predictors = "age",
    time_var = "Tstart",
    time_spline_df = NULL
  )

  expect_equal(
    attr(stats::terms(formula), "term.labels"),
    c("Tstart", "age")
  )
})

test_that("censoring formula supports an intercept-only model", {
  formula <- make_censoring_formula(
    "survival::Surv(Tstart, Tstop, censoring)"
  )

  expect_equal(attr(stats::terms(formula), "term.labels"), character())
  expect_equal(attr(stats::terms(formula), "intercept"), 1L)
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

test_that("Cox baseline hazard is carried forward between event times", {
  dat <- tibble::tibble(
    id = c(1, 2, 3),
    Tstart = c(0, 1.5, 0),
    Tstop = c(1, 2, 2),
    censoring = c(1, 0, 0)
  )
  clones <- list(Control = dat, Surgery = dat)

  result <- estimate_censoring(
    clones,
    method = "Cox"
  )

  fit <- survival::coxph(
    survival::Surv(Tstart, Tstop, censoring) ~ 1,
    data = dat,
    ties = "efron"
  )
  base_hazard <- survival::basehaz(fit, centered = FALSE)
  hazard_index <- findInterval(dat$Tstart, base_hazard$time)
  expected_hazard <- numeric(nrow(dat))
  expected_hazard[hazard_index > 0L] <-
    base_hazard$hazard[hazard_index[hazard_index > 0L]]

  expect_equal(result$Control$hazard, expected_hazard)
  expect_gt(result$Control$hazard[dat$Tstart == 1.5], 0)
  expect_equal(
    result$Control$P_uncens,
    exp(-expected_hazard * exp(result$Control$lin_pred))
  )
})

test_that("cumulative uncensoring uses probabilities from prior intervals", {
  dat <- tibble::tibble(
    id = c(1, 1, 1, 2, 2),
    Tstart = c(0, 1, 2, 0, 1),
    Tstop = c(1, 2, 3, 1, 2)
  )

  result <- cumulative_uncensoring(
    dat,
    p_censoring = c(0.2, 0.5, 0.1, 0.25, 0.4)
  )

  expect_equal(result, c(1, 0.8, 0.4, 1, 0.75))
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
