make_test_clones <- function() {
  list(
    Control = tibble::tibble(
      id = c(1, 2),
      fup = c(5, 10),
      outcome = c(1, 0),
      censoring = c(0, 1)
    ),
    Surgery = tibble::tibble(
      id = c(1, 2),
      fup = c(6, 8),
      outcome = c(0, 1),
      censoring = c(1, 0)
    )
  )
}

test_that("timestamp and split helpers are internal", {
  exports <- getNamespaceExports("clonecensorweighting")

  expect_false("create_timestamp_table" %in% exports)
  expect_false("split_at_timestamp" %in% exports)
})

test_that("create_timestamp_table creates ordered unique event times", {
  result <- create_timestamp_table(make_test_clones(), "fup")

  expect_equal(
    result,
    tibble::tibble(tevent = c(5, 6, 8, 10), ID_t = 1:4)
  )
})

test_that("split_at_timestamp splits each clone arm at common timestamps", {
  result <- split_at_timestamp(
    make_test_clones(),
    clone_followup = "fup",
    t_events = c(5, 6, 8, 10),
    event = "outcome"
  )

  expect_named(result, c("Control", "Surgery"))
  expect_equal(result$Control$Tstart, c(0, 0, 5, 6, 8))
  expect_equal(result$Control$fup, c(5, 5, 6, 8, 10))
  expect_equal(result$Control$outcome, c(1, 0, 0, 0, 0))
  expect_equal(result$Surgery$Tstart, c(0, 5, 0, 5, 6))
  expect_equal(result$Surgery$fup, c(5, 6, 5, 6, 8))
  expect_equal(result$Surgery$outcome, c(0, 0, 0, 0, 1))
})

test_that("create_final_data combines outcome and censoring intervals", {
  result <- create_final_data(
    make_test_clones(),
    clone_followup = "fup",
    clone_outcome = "outcome",
    clone_censoring = "censoring",
    col_ids = "id"
  )

  expect_named(result, c("Control", "Surgery"))
  expect_equal(result$Control$Tstop, result$Control$fup)
  expect_equal(result$Surgery$Tstop, result$Surgery$fup)
  expect_equal(result$Control$ID_t, c(0, 0, 1, 2, 3))
  expect_equal(result$Surgery$ID_t, c(0, 1, 0, 1, 2))
  expect_equal(result$Control$censoring, c(0, 0, 0, 0, 1))
  expect_equal(result$Surgery$censoring, c(0, 1, 0, 0, 0))
  expect_equal(result$Surgery$outcome, c(0, 0, 0, 0, 1))
})

test_that("create_final_data reports missing clone columns by arm", {
  clones <- make_test_clones()
  clones$Control$censoring <- NULL

  expect_error(
    create_final_data(
      clones,
      clone_followup = "fup",
      clone_outcome = "outcome",
      clone_censoring = "censoring",
      col_ids = "id"
    ),
    "In clone arm `Control`: Missing required columns: censoring"
  )
})
