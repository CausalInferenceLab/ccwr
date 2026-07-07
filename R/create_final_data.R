
#' Create training data for censoring probability estimation
#'
#' @param clones A list of data frame. Each element of the list represents
#'   each treatment arm. This version of clones must contain a column that
#'   represents an emulated follow up time (corresponding to `clone_followup`
#'   argument), an emulated outcome (correspodning to `clone_outcome`), and a
#'   binary indicator variable that represents whether the observation violates
#'   arm's policy or not (corresponding to `clone_censoring` argument).
#' @param clone_followup A column name that represents the emulated follow up
#'   time in each arm of clones. The variable should exists in each element
#'   data frame of `clones` argument.
#' @param clone_outcome A column name that represents the emulated outcome
#'   in each arm of clones. The variable should exists in each element
#'   data frame of `clones` argument, and the variable value should be binary
#'   (0 or 1).
#' @param clone_censoring A column name that represent whether the observation
#'   violates arm's policy or not. The variable should exists in each element
#'   data frame of `clones` argument, and the variable value should be binary
#'   (0 or 1).
#' @param col_ids A vector of column names that a combination of their values
#'   uniquely identifies each observation.
#' @param timestamp_start A new variable name to denote start time of each
#'   subrecord of observations in a long-form data.
#' @param id A new variable name for a unique observation identifier, to
#'   represents that multiple rows in output data frame is associated with the
#'   same observation.
#' @param timestamp_stop A new variable name to denote end time of each
#'   subrecord of observations in a long-form data.
#'
#' @returns A list of long-form data frames. Each data frame represents each
#'   clone arm. Each row of the long-form data frame represents a subrecord of
#'   each observation associated with each specific time interval. The first
#'   subrecord starts with time 0, and the rows are expanded up to
#'   `clone_followup`, where cut times are determined by `t_events` argument.
#'
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
#'   clones_censored,
#'   clone_followup = "fup",
#'   clone_outcome = "outcome",
#'   clone_censoring = "censoring",
#'   col_ids = "id"
#' )
create_final_data <- function(
  clones,
  clone_followup,
  clone_outcome,
  clone_censoring,
  col_ids,
  timestamp_start = "Tstart",
  id = "ID",
  timestamp_stop = "Tstop"
) {
  .assert_clone_columns(
    clones,
    c(clone_followup, clone_outcome, clone_censoring, col_ids)
  )

  df_timestamp <- create_timestamp_table(clones, clone_followup)

  clones_splitted_by_outcome <- split_at_timestamp(
    clones,
    clone_followup,
    df_timestamp$tevent,
    clone_outcome,
    timestamp_start,
    id
  )

  clones_splitted_by_censoring <- split_at_timestamp(
    clones,
    clone_followup,
    df_timestamp$tevent,
    clone_censoring,
    timestamp_start,
    id
  )

  # merge two tables and create column to represent end of timestamp
  df_timestamp_with_time_zero <-
    dplyr::bind_rows(
      dplyr::tibble(tevent = 0, ID_t = 0),
      df_timestamp
    )
  df_timestamp_with_time_zero <-
    df_timestamp_with_time_zero[
      !duplicated(df_timestamp_with_time_zero$tevent),
      ,
      drop = FALSE
    ]

  n_clones <- length(clones)
  arms <- names(clones)
  res <- vector("list", length = n_clones)
  names(res) <- arms

  for (i in seq_len(n_clones)) {
    x <-
      clones_splitted_by_outcome[[i]] |>
      dplyr::select(-dplyr::all_of(clone_censoring))

    y <-
      clones_splitted_by_censoring[[i]] |>
      dplyr::select(
        dplyr::all_of(c(col_ids, clone_followup, clone_censoring))
      )

    res[[i]] <-
      dplyr::inner_join(
        x,
        y,
        by = c(col_ids, clone_followup)
      ) |>
      dplyr::mutate(!!timestamp_stop := .data[[clone_followup]])

    # Merge with timestamp table
    res[[i]] <-
      res[[i]] |>
      dplyr::left_join(
        df_timestamp_with_time_zero,
        by = stats::setNames("tevent", timestamp_start)
      )
  }

  res
}
