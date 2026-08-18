test_that("make_surv_response returns a Surv object", {
  trial_data <- tibble::tibble(
    follow_up = c(5, 8),
    event = c(1, 0)
  )

  surv_response <- make_surv_response(
    data = trial_data,
    follow_up = "follow_up",
    event = "event"
  )

  expect_s3_class(surv_response, "Surv")
  expect_equal(length(surv_response), 2)
})
