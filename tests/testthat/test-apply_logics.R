test_that("definition of the outcome and survival time for each patient in each arm is correctly match to Maringe et al. (2020) for analysis model", {
  data(patients13)

  clones <- clone_arms(patients13, c("Control", "Surgery"))
  policy <- create_policy_A(
    arms = c("Control", "Surgery"),
    treatment = "surgery",
    time_to_treatment = "time_to_surgery",
    grace_period = 182,
    outcome = "death",
    followup = "followup",
    clone_outcome = "analysis_outcome",
    clone_followup = "analysis_followup"
  )
  res <- apply_logics(clones, policy)

  # check control arm results
  expect_equal(
    res$Control |>
      dplyr::select(id, analysis_outcome, analysis_followup),
    patients13 |>
      dplyr::transmute(
        id,
        analysis_outcome = dplyr::case_when(
          id %in% c("A", "B", "C", "D", "E") ~ 0,
          id %in% c("F", "G", "H", "I", "J", "K", "L", "M") ~ death,
        ),
        analysis_followup = dplyr::case_when(
          id %in% c("A", "B", "C", "D", "E") ~ time_to_surgery,
          id %in% c("F", "G", "H", "I", "J", "K", "L", "M") ~ followup,
        )
      )
  )

  # check surgery arm results
  expect_equal(
    res$Surgery |>
      dplyr::select(id, analysis_outcome, analysis_followup),
    patients13 |>
      dplyr::transmute(
        id,
        analysis_outcome = dplyr::case_when(
          id %in% c("A", "B", "C", "D", "E") ~ death,
          id %in% c("F", "G", "H", "I", "J", "L", "M") ~ 0,
          id %in% c("K") ~ 1,
        ),
        analysis_followup = dplyr::case_when(
          id %in% c("A", "B", "C", "D", "E") ~ followup,
          id %in% c("F", "G", "H", "I", "J", "K", "L", "M") ~ pmin(182, followup),
        )
      )
  )
})
