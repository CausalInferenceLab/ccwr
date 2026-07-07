make_lungcancer_workflow <- function() {
  data(lungcancer, package = "clonecensorweighting")

  arms <- c("Control", "Surgery")
  clones <- clone_arms(lungcancer, arms)
  policies <- create_policy_A(
    arms,
    treatment = "surgery",
    time_to_treatment = "timetosurgery",
    grace_period = 182.62,
    outcome = "death",
    followup = "fup_obs",
    clone_outcome = "outcome",
    clone_followup = "fup"
  )
  clones_policy <- apply_logics(clones, policies)
  censoring_logics <- create_censoring_logics_A(
    arms,
    treatment = "surgery",
    time_to_treatment = "timetosurgery",
    grace_period = 182.62,
    followup = "fup_obs",
    clone_censoring = "censoring",
    clone_uncensored_followup = "fup_uncensored"
  )
  clones_censored <- apply_logics(clones_policy, censoring_logics)
  clones_final <- create_final_data(
    clones_censored,
    clone_followup = "fup",
    clone_outcome = "outcome",
    clone_censoring = "censoring",
    col_ids = "id"
  )
  clones_estimated <- estimate_censoring(
    clones_final,
    predictors = c("age", "sex"),
    method = "pooled_logit"
  )
  clones_weighted <- weight_cases(clones_estimated)

  list(
    clones = clones,
    clones_final = clones_final,
    clones_estimated = clones_estimated,
    clones_weighted = clones_weighted
  )
}

test_that("lungcancer data runs through clone-censor-weighting workflow", {
  workflow <- make_lungcancer_workflow()

  expect_named(workflow$clones, c("Control", "Surgery"))
  expect_equal(vapply(workflow$clones, nrow, integer(1)), c(Control = 200L, Surgery = 200L))
  expect_equal(
    vapply(workflow$clones_final, nrow, integer(1)),
    c(Control = 8792L, Surgery = 13980L)
  )

  required_final_columns <- c("id", "Tstart", "Tstop", "outcome", "censoring", "ID_t")
  expect_true(all(required_final_columns %in% names(workflow$clones_final$Control)))
  expect_true(all(c("p_cens_den", "P_uncens") %in% names(workflow$clones_estimated$Control)))
  expect_true(all(c("weight_Cox") %in% names(workflow$clones_weighted$Control)))
  expect_true(all(is.finite(workflow$clones_weighted$Control$weight_Cox)))
  expect_true(all(workflow$clones_weighted$Control$weight_Cox > 0))
})

test_that("lungcancer weighted data supports emulated trial estimation", {
  workflow <- make_lungcancer_workflow()

  cox_fit <- emul_estimate(
    workflow$clones_weighted,
    method = "Cox",
    weights = "weight_Cox",
    predictors = c("age", "sex")
  )
  km_fit <- emul_estimate(
    workflow$clones_weighted,
    method = "KM",
    weights = "weight_Cox"
  )
  boot <- emul_estimate_bootstrap(
    workflow$clones_weighted,
    method = "Cox",
    weights = "weight_Cox",
    predictors = c("age", "sex"),
    n_bootstrap = 3,
    seed = 1
  )

  expect_s3_class(cox_fit, "coxph")
  expect_true(all(is.finite(stats::coef(cox_fit))))
  expect_s3_class(km_fit, "survfit")
  expect_length(boot$estimates, 3)
  expect_true(all(is.finite(boot$estimates)))
  expect_true(is.finite(boot$ci_lower))
  expect_true(is.finite(boot$ci_upper))
})
