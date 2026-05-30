#############################################################################
# TRIAL EMULATION: SURGERY WITHIN 6 MONTHS AMONG OLDER LUNG CANCER PATIENTS
# Refactored Clemence Leyrat's R code published on https://academic.oup.com/ije/article/49/5/1719/5835351
#############################################################################

# Read data
tab <- read.csv(
  "example_data/ije-2019-08-1035-File008.csv",
  sep = ",",
  header = TRUE
)

# print(tab)

#########################
# Step: Cloning
#########################

#' Duplicate data frame for each treatment arm to emulate
#'
#' @param data Input data frame that contains all observations of interest.
#'   Each row represents an observation, and columns include
#'   observation identifiers, binary treatment variable (0/1),
#'   time to treatement (continuous), binary outcome variable (0/1),
#'   observed followup time (continuous), and covariates.
#' @param arms Character vector that each element represents each arm's name.
#'
#' @returns A list of data frame.
#'   Each element of list is associated with each arm.
#'
#' @export
#' @examples
#' tab <- read.csv(
#'   "example_data/ije-2019-08-1035-File008.csv",
#'   sep = ",",
#'   header = TRUE
#' )
#' clones <- clone_arms(tab, c("Control", "Surgery"))
clone_arms <- function(data, arms) {
  n <- length(arms)
  if (n <= 1) {
    stop("`arms` must have more than one value.")
  }

  res <- vector("list", length = n)
  names(res) <- arms
  for (i in seq_len(n)) {
    res[[arms[[i]]]] <- data
  }

  res
}

arms <- c("Control", "Surgery")

clones <- clone_arms(tab, arms)

apply_policy_old <- function(
  arms,
  clones,
  policies,
  clone_outcome,
  clone_followup
) {
  res <- vector("list", length = length(arms))
  names(res) <- arms

  for (arm in arms) {
    res[[arm]] <-
      clones[[arm]] |>
      dplyr::mutate(
        {{ clone_outcome }} := dplyr::case_when(
          !!!rlang::parse_exprs(policies[[arm]][[clone_outcome]]),
          TRUE ~ NA
        ),
        {{ clone_followup }} := dplyr::case_when(
          !!!rlang::parse_exprs(policies[[arm]][[clone_followup]]),
          TRUE ~ NA
        )
      )
  }

  res
}

#' Apply case_when logics (e.g. policy, censoring) to clones
#'
#' @param clones A list of data frame. Each element of the list represents
#'   each treatment arm.
#' @param logics A nested list. Each element of outer list represents each
#'   treatment arm. Each element of inner list represents each new variable
#'   to be created by applying logics. Each element of inner list containts a
#'   character vector that represents a sequence of logics to be passed into
#'   `case_when()` call to determine a value of the new variable.
#'
#' @returns A list of data frame that each data frame include new variables
#'   created by the provided logics.
#'
#' @export
#' @examples
apply_logics <- function(
  clones,
  logics
) {
  stopifnot(names(clones) == names(logics))

  arms <- names(clones)

  res <- vector("list", length = length(arms))
  names(res) <- arms

  for (arm in arms) {
    vars <- names(logics[[arm]])

    res[[arm]] <- clones[[arm]]

    for (var in vars) {
      res[[arm]] <-
        res[[arm]] |>
        dplyr::mutate(
          {{ var }} := dplyr::case_when(
            !!!rlang::parse_exprs(logics[[arm]][[var]]),
            TRUE ~ NA
          )
        )
    }
  }

  res
}


#' Generate clone policy for scenario A
#'
#' @param arms A character vector of length 2. The first element represents
#'   a name of the untreated arm, and the second element represents a name of
#'   the treated arm.
#' @param treatment A name of variable that represents whether each observation
#'   was treated or not in observational data. The treatment variable should
#'   exists in data frame that the return value of this function will be
#'   applied, and the treatment variable value in the data frame should be
#'   either 0 and 1, i.e. binary treatment.
#' @param time_to_treatment A name of variable that represent time to
#'   treatment in the observational data. The time-to-treatment variable should
#'   exists in data frame that the return value of this function will be
#'   applied, and the value in the data frame should be either numeric value
#'   or NA if the observation was untreated in the observational data.
#' @param grace_period A numeric value to represent grace period of treatment.
#'   Treatment policy is assumed to be "provide treatment within the grace
#'   period."
#' @param outcome A name of variable that represent outcome in the
#'   observational data. The outcome variable should exists in data frame that
#'   the return value of this function will be applied, and the outcome
#'   variable value in the data frame should be either 0 and 1, i.e. binary
#'   outcome.
#' @param followup A name of variable that represent follow up time in the
#'   observational data. The follow up time variable should exists in data
#'   frame that the return value of this function will be applied, and the
#'   follow up time variable value in the data frame should be numeric.
#' @param clone_outcome A name of variable to be newly created to represent
#'   emulated outcome in cloned data frame. The new variable name should not
#'   already exist in the data frame that the return value of this function
#'   will be applied, to avoid accidental overwriting.
#' @param clone_followup A name of variable to be newly created to represent
#'   emulated follow up time in cloned data frame. The new variable name should
#'   not already exist in the data frame that the return value of this function
#'   will be applied, to avoid accidental overwriting.
#'
#' @returns A nested list. The first element of the outer list represents
#'   untreated arm, while the second element of the outer list represents
#'   treated arm. For each element of outer list, the first element of the inner
#'   list represents emulated outcome, and the second element of the inner list
#'   represents emulated follow up time. Each element of the inner list
#'   represents a sequence of logics to be passed into `case_when()` when
#'   creating new variables for emulated outcome and follow up time.
#'
#' @export
#' @examples
create_policy_A <- function(
  arms,
  treatment,
  time_to_treatment,
  grace_period,
  outcome,
  followup,
  clone_outcome = ".outcome",
  clone_followup = ".fup"
) {
  res <- list(
    list(
      outcome = c(
        glue::glue(
          "{treatment} == 1 & {time_to_treatment} <= {grace_period} ~ 0"
        ),
        glue::glue(
          "{treatment} == 0 | ({treatment} == 1 & {time_to_treatment} > {grace_period}) ~ {outcome}"
        )
      ),
      followup = c(
        glue::glue(
          "{treatment} == 1 & {time_to_treatment} <= {grace_period} ~ {time_to_treatment}"
        ),
        glue::glue(
          "{treatment} == 0 | ({treatment} == 1 & {time_to_treatment} > {grace_period}) ~ {followup}"
        )
      )
    ),
    list(
      outcome = c(
        glue::glue(
          "{treatment} == 1 & {time_to_treatment} <= {grace_period} ~ {outcome}"
        ),
        glue::glue(
          "{treatment} == 0 & {followup} <= {grace_period} ~ {outcome}"
        ),
        glue::glue(
          "({treatment} == 0 & {followup} > {grace_period}) | ({treatment} == 0 & {time_to_treatment} > {grace_period}) ~ 0"
        )
      ),
      followup = c(
        glue::glue(
          "{treatment} == 1 & {time_to_treatment} <= {grace_period} ~ {followup}"
        ),
        glue::glue(
          "{treatment} == 0 & {followup} <= {grace_period} ~ {followup}"
        ),
        glue::glue(
          "({treatment} == 0 & {followup} > {grace_period}) | ({treatment} == 0 & {time_to_treatment} > {grace_period}) ~ {grace_period}"
        )
      )
    )
  )

  names(res) <- arms
  for (i in seq_along(res)) {
    names(res[[i]]) <- c(clone_outcome, clone_followup)
  }

  res
}

policies <- create_policy_A(
  arms,
  "surgery",
  "timetosurgery",
  182.62,
  "death",
  "fup_obs",
  "outcome",
  "fup"
)

testthat::expect_equal(
  policies,
  list(
    Control = list(
      outcome = c(
        "surgery == 1 & timetosurgery <= 182.62 ~ 0",
        "surgery == 0 | (surgery == 1 & timetosurgery > 182.62) ~ death"
      ),
      fup = c(
        "surgery == 1 & timetosurgery <= 182.62 ~ timetosurgery",
        "surgery == 0 | (surgery == 1 & timetosurgery > 182.62) ~ fup_obs"
      )
    ),
    Surgery = list(
      outcome = c(
        "surgery == 1 & timetosurgery <= 182.62 ~ death",
        "surgery == 0 & fup_obs <= 182.62 ~ death",
        "(surgery == 0 & fup_obs > 182.62) | (surgery == 0 & timetosurgery > 182.62) ~ 0"
      ),
      fup = c(
        "surgery == 1 & timetosurgery <= 182.62 ~ fup_obs",
        "surgery == 0 & fup_obs <= 182.62 ~ fup_obs",
        "(surgery == 0 & fup_obs > 182.62) | (surgery == 0 & timetosurgery > 182.62) ~ 182.62"
      )
    )
  )
)

clones_policy_applied <- apply_logics(clones, policies)

testthat::expect_equal(
  clones_policy_applied,
  apply_policy_old(arms, clones, policies, "outcome", "fup")
)

#########################
# Step: Censoring
#########################

censor_arms_old <- function(
  arms,
  clones,
  logics,
  censoring,
  followup_uncensored
) {
  res <- vector("list", length = length(arms))
  names(res) <- arms

  for (arm in arms) {
    res[[arm]] <-
      clones[[arm]] |>
      dplyr::mutate(
        {{ censoring }} := dplyr::case_when(
          !!!rlang::parse_exprs(logics[[arm]][["censoring"]]),
          TRUE ~ NA
        ),
        {{ followup_uncensored }} := dplyr::case_when(
          !!!rlang::parse_exprs(logics[[arm]][["fup_uncensored"]]),
          TRUE ~ NA
        )
      )
  }

  res
}

#' Generate censoring logic for scenario A
#'
#' @param arms A character vector of length 2. The first element represents
#'   a name of the untreated arm, and the second element represents a name of
#'   the treated arm.
#' @param treatment A name of variable that represents whether each observation
#'   was treated or not in observational data. The treatment variable should
#'   exists in data frame that the return value of this function will be
#'   applied, and the treatment variable value in the data frame should be
#'   either 0 and 1, i.e. binary treatment.
#' @param time_to_treatment A name of variable that represent time to
#'   treatment in the observational data. The time-to-treatment variable should
#'   exists in data frame that the return value of this function will be
#'   applied, and the value in the data frame should be either numeric value
#'   or NA if the observation was untreated in the observational data.
#' @param grace_period A numeric value to represent grace period of treatment.
#'   Treatment policy is assumed to be "provide treatment within the grace
#'   period."
#' @param followup A name of variable that represent follow up time in the
#'   observational data. The follow up time variable should exists in data
#'   frame that the return value of this function will be applied, and the
#'   follow up time variable value in the data frame should be numeric.
#' @param clone_censoring A name of binary indicator variable to be newly
#'   created to represent whether the observation violates arm's policy or not.
#'   The new variable name should not already exist in the data frame that the
#'   return value of this function will be applied, to avoid accidental
#'   overwriting.
#' @param clone_uncensored_followup A name of variable to be newly created to
#'   represent the earliest time that the value of the emulated censoring binary
#'   indicator (i.e. variable to be named according to `clone_censoring`
#'   argument) value can be determined for each observation. The new variable
#'   name should not already exist in the data frame that the return value of
#'   this function will be applied, to avoid accidental overwriting.
#'
#' @returns A nested list. The first element of the outer list represents
#'   untreated arm, while the second element of the outer list represents
#'   treated arm. For each element of outer list, the first element of the inner
#'   list represents emulated censoring binary indicator (0/1) that represents
#'   whether the observation violated the arm's policy or not within the grace
#'   period, and the second. The second element of the inner list represents
#'   the earliest time that the value of the emulated censoring binary
#'   indicator can be determined for each observation. Each element of the
#'   inner list represents a sequence of logics to be passed into `case_when()`
#'   when creating new variables for emulated censoring indicator and censoring
#'   time.
#'
#' @export
#' @examples
create_censoring_logics_A <- function(
  arms,
  treatment,
  time_to_treatment,
  grace_period,
  followup,
  clone_censoring = ".censoring",
  clone_uncensored_followup = ".fup_uncensored"
) {
  res <- list(
    list(
      censoring = c(
        glue::glue(
          "{treatment} == 1 & {time_to_treatment} <= {grace_period} ~ 1"
        ),
        glue::glue(
          "{treatment} == 0 & {followup} <= {grace_period} ~ 0"
        ),
        glue::glue(
          "({treatment} == 0 & {followup} > {grace_period}) | ({treatment} == 1 & {time_to_treatment} > {grace_period}) ~ 0"
        )
      ),
      fup_uncensored = c(
        glue::glue(
          "{treatment} == 1 & {time_to_treatment} <= {grace_period} ~ {time_to_treatment}"
        ),
        glue::glue(
          "{treatment} == 0 & {followup} <= {grace_period} ~ {followup}"
        ),
        glue::glue(
          "({treatment} == 0 & {followup} > {grace_period}) | ({treatment} == 1 & {time_to_treatment} > {grace_period}) ~ {grace_period}"
        )
      )
    ),
    list(
      censoring = c(
        glue::glue(
          "{treatment} == 1 & {time_to_treatment} <= {grace_period} ~ 0"
        ),
        glue::glue(
          "{treatment} == 0 & {followup} <= {grace_period} ~ 0"
        ),
        glue::glue(
          "({treatment} == 0 & {followup} > {grace_period}) | ({treatment} == 1 & {time_to_treatment} > {grace_period}) ~ 1"
        )
      ),
      fup_uncensored = c(
        glue::glue(
          "{treatment} == 1 & {time_to_treatment} <= {grace_period} ~ {time_to_treatment}"
        ),
        glue::glue(
          "{treatment} == 0 & {followup} <= {grace_period} ~ {followup}"
        ),
        glue::glue(
          "({treatment} == 0 & {followup} > {grace_period}) | ({treatment} == 1 & {time_to_treatment} > {grace_period}) ~ {grace_period}"
        )
      )
    )
  )

  names(res) <- arms
  for (i in seq_along(res)) {
    names(res[[i]]) <- c(clone_censoring, clone_uncensored_followup)
  }

  res
}

censoring_logics <- create_censoring_logics_A(
  arms,
  "surgery",
  "timetosurgery",
  182.62,
  "fup_obs",
  "censoring",
  "fup_uncensored"
)

testthat::expect_equal(
  censoring_logics,
  list(
    Control = list(
      censoring = c(
        "surgery == 1 & timetosurgery <= 182.62 ~ 1",
        "surgery == 0 & fup_obs <= 182.62 ~ 0",
        "(surgery == 0 & fup_obs > 182.62) | (surgery == 1 & timetosurgery > 182.62) ~ 0"
      ),
      fup_uncensored = c(
        "surgery == 1 & timetosurgery <= 182.62 ~ timetosurgery",
        "surgery == 0 & fup_obs <= 182.62 ~ fup_obs",
        "(surgery == 0 & fup_obs > 182.62) | (surgery == 1 & timetosurgery > 182.62) ~ 182.62"
      )
    ),
    Surgery = list(
      censoring = c(
        "surgery == 1 & timetosurgery <= 182.62 ~ 0",
        "surgery == 0 & fup_obs <= 182.62 ~ 0",
        "(surgery == 0 & fup_obs > 182.62) | (surgery == 1 & timetosurgery > 182.62) ~ 1"
      ),
      fup_uncensored = c(
        "surgery == 1 & timetosurgery <= 182.62 ~ timetosurgery",
        "surgery == 0 & fup_obs <= 182.62 ~ fup_obs",
        "(surgery == 0 & fup_obs > 182.62) | (surgery == 1 & timetosurgery > 182.62) ~ 182.62"
      )
    )
  )
)

clones_censored <- apply_logics(clones_policy_applied, censoring_logics)

testthat::expect_equal(
  clones_censored,
  censor_arms_old(
    arms,
    clones_policy_applied,
    censoring_logics,
    "censoring",
    "fup_uncensored"
  )
)


#########################
# Step: Weighting (No bootstrapping)
#########################

# create timestamp table
create_timestamp_table_old <- function(clones) {
  timestamps <- sapply(clones, `[[`, i = "fup", simplify = FALSE)
  t_events <- sort(unique(unlist(timestamps)))
  res <- dplyr::tibble(tevent = t_events, ID_t = seq_along(t_events))

  res
}

#' Create timestamp table
#'
#' @param clones A list of data frame. Each element of the list represents
#'   each treatment arm. This version of clones must contain a column that
#'   represents a emulated follow up time.
#' @param clone_followup A column name in emulcated clone (i.e. `clones`)
#'   that represents the emulated follow up time in cloned data frame.
#'
#' @returns A data frame with two columns: `tevent` and `ID_t`.
#'   `tevent` represents a timestamp that outcome event can occur based on
#'   observed data. `ID_t` represents an enumerated identifier of each
#'   timestamp, from 1 to n where n represents the number of unique `tevent`
#'   value.
#'
#' @export
#' @examples
create_timestamp_table <- function(clones, clone_followup) {
  timestamps <- sapply(clones, `[[`, i = clone_followup, simplify = FALSE)
  t_events <- sort(unique(unlist(timestamps)))
  res <- dplyr::tibble(tevent = t_events, ID_t = seq_along(t_events))

  res
}


# split data at each time event
split_at_timestamp_old <- function(clones, t_events, event) {
  arms <- names(clones)

  res <- vector("list", length = length(arms))
  names(res) <- arms
  for (arm in arms) {
    res[[arm]] <-
      clones[[arm]] |>
      survival::survSplit(
        cut = t_events,
        end = "fup",
        start = "Tstart",
        event = event,
        id = "ID"
      )
  }

  res
}

#' Split each observation into multiple subrecords at each time cut
#'
#' @param clones A list of data frame. Each element of the list represents
#'   each treatment arm. This version of clones must contain a column that
#'   represents a emulated follow up time (corresponding to `clone_followup`
#'   argument) and an event of interest (corresponding to `event` argument).
#' @param clone_followup A column name in emulcated clone (i.e. `clones`)
#'   that represents the emulated follow up time in cloned data frame.
#' @param t_events A vector of timestamp that outcome event can occur based on
#'   observed data.
#' @param event A variable name of an event of interest. The variable should
#'   exists in each data frame that is an element of `clones` argument, and
#'   the variable value should be a binary (0 or 1).
#' @param timestamp_start A new variable name to denote start time.
#' @param id A new variable name for a unique observation identifier, to
#'   represents that multiple rows in output data frame is associated with the
#'   same observation.
#'
#' @returns A list of long-form data frames. Each data frame represents each
#'   clone arm. Each row of the long-form data frame represents a subrecord of
#'   each observation associated with each specific time interval. The first
#'   subrecord starts with time 0, and the rows are expanded up to
#'   `clone_followup`, where cut times are determined by `t_events` argument.
#'
#' @export
#' @examples
split_at_timestamp <- function(
  clones,
  clone_followup,
  t_events,
  event,
  timestamp_start = "Tstart",
  id = "ID"
) {
  arms <- names(clones)

  res <- vector("list", length = length(arms))
  names(res) <- arms
  for (arm in arms) {
    res[[arm]] <-
      clones[[arm]] |>
      survival::survSplit(
        cut = t_events,
        end = clone_followup,
        start = timestamp_start,
        event = event,
        id = id
      )
  }

  res
}

df_timestamp <- create_timestamp_table(clones_censored, "fup")

clones_splitted_by_outcome <- split_at_timestamp(
  clones_censored,
  "fup",
  df_timestamp$tevent,
  "outcome",
  "Tstart",
  "ID"
)

clones_splitted_by_censoring <- split_at_timestamp(
  clones_censored,
  "fup",
  df_timestamp$tevent,
  "censoring",
  "Tstart",
  "ID"
)


# Function to create training data for censoring probability estimation
# TO DO: Eliminate dependency on {purrr}
create_final_data_old <- function(clones) {
  df_timestamp <- create_timestamp_table_old(clones)

  clones_splitted_by_outcome <- split_at_timestamp_old(
    clones,
    df_timestamp$tevent,
    "outcome"
  )

  clones_splitted_by_censoring <- split_at_timestamp_old(
    clones,
    df_timestamp$tevent,
    "censoring"
  )

  # merge two tables and create Tstop column
  clones_splitted <-
    purrr::map2(
      clones_splitted_by_outcome,
      clones_splitted_by_censoring,
      \(x, y) {
        dplyr::inner_join(
          x |> dplyr::select(!censoring),
          y |> dplyr::select(id, fup, censoring),
          by = dplyr::join_by(id, fup)
        ) |>
          dplyr::mutate(Tstop = fup)
      }
    )

  # Merge with timestamp table
  res <-
    purrr::map(
      clones_splitted,
      \(x, df_timestamp) {
        dplyr::left_join(x, df_timestamp, by = dplyr::join_by(Tstart == tevent))
      },
      df_timestamp = df_timestamp |>
        dplyr::bind_rows(dplyr::tibble(tevent = 0, ID_t = 0))
    )

  res
}

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
    df_timestamp |>
    dplyr::bind_rows(dplyr::tibble(tevent = 0, ID_t = 0))

  n_clones <- length(clones)
  arms <- names(clones)
  res <- vector("list", length = n_clones)
  names(res) <- arms

  for (i in seq_len(n_clones)) {
    x <-
      clones_splitted_by_outcome[[i]] |>
      dplyr::select(!{{ clone_censoring }})

    y <-
      clones_splitted_by_censoring[[i]] |>
      dplyr::select(
        all_of(col_ids),
        {{ clone_followup }},
        {{ clone_censoring }}
      )

    res[[i]] <-
      dplyr::inner_join(
        x,
        y,
        by = dplyr::join_by({{ col_ids }}, {{ clone_followup }})
      ) |>
      dplyr::mutate({{ timestamp_stop }} := .data[[clone_followup]])

    # Merge with timestamp table
    res[[i]] <-
      res[[i]] |>
      dplyr::left_join(
        df_timestamp_with_time_zero,
        by = dplyr::join_by({{ timestamp_start }} == tevent)
      )
  }

  res
}


clones_final_old <- create_final_data_old(clones_censored)
clones_final <- create_final_data(
  clones_censored,
  "fup",
  "outcome",
  "censoring",
  c("id"),
  "Tstart",
  "ID",
  "Tstop"
)

testthat::expect_equal(clones_final, clones_final_old)


# Function to estimate (un)censoring probability
make_censoring_formula <- function(
  response,
  predictors = NULL,
  time_var = NULL
) {
  terms <- unique(c(time_var, predictors))
  if (is.null(terms)) {
    terms <- character()
  }
  formula <- stats::reformulate(terms, response = response)
  environment(formula) <- parent.frame()
  formula
}

cumulative_uncensoring <- function(
  data,
  p_censoring,
  id = "id",
  time_start = "Tstart",
  time_stop = "Tstop",
  eps = 1e-6
) {
  missing_columns <- setdiff(c(id, time_start, time_stop), names(data))
  if (length(missing_columns) > 0L) {
    stop(
      "Missing required columns in clone data: ",
      paste(missing_columns, collapse = ", "),
      call. = FALSE
    )
  }
  if (length(p_censoring) != nrow(data)) {
    stop("`p_censoring` must have one value per row of `data`.", call. = FALSE)
  }

  p_uncensored_interval <- 1 - p_censoring
  p_uncensored <- numeric(nrow(data))
  ordered_rows <- order(data[[id]], data[[time_start]], data[[time_stop]])
  rows_by_id <- split(ordered_rows, data[[id]][ordered_rows])

  for (rows in rows_by_id) {
    interval_prob <- p_uncensored_interval[rows]
    p_uncensored[rows] <- c(1, cumprod(interval_prob[-length(interval_prob)]))
  }

  pmax(p_uncensored, eps)
}

add_baseline_predictors <- function(
  data,
  predictors = NULL,
  id = "id",
  time_start = "Tstart",
  time_stop = "Tstop",
  prefix = ".baseline_"
) {
  if (is.null(predictors) || length(predictors) == 0L) {
    return(list(data = data, predictors = character()))
  }
  predictors <- unique(predictors)

  missing_columns <- setdiff(
    c(id, time_start, time_stop, predictors),
    names(data)
  )
  if (length(missing_columns) > 0L) {
    stop(
      "Missing required columns in clone data: ",
      paste(missing_columns, collapse = ", "),
      call. = FALSE
    )
  }

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

estimate_censoring <- function(
  clones,
  predictors = NULL,
  method = c("Cox", "pooled_logit", "stabilized_logit"),
  numerator_predictors = NULL,
  id = "id",
  time_start = "Tstart",
  time_stop = "Tstop",
  eps = 1e-6
) {
  method <- match.arg(method)
  if (method == "stabilized_logit" && is.null(numerator_predictors)) {
    numerator_predictors <- predictors
  }
  arms <- names(clones)

  res <- vector("list", length = length(arms))
  names(res) <- arms

  surv_response <- paste0(
    "survival::Surv(",
    time_start,
    ", ",
    time_stop,
    ", censoring)"
  )
  cox_formula <- make_censoring_formula(
    surv_response,
    predictors
  )
  pooled_formula <- make_censoring_formula(
    "censoring",
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
      ) #We can also includes their interactions

      lin_pred <- stats::predict(
        ms_cens,
        newdata = dat,
        type = "lp",
        reference = "zero"
      )

      # Estimating the cumulative hazard (when covariates=0)
      base_hazard <- dplyr::as_tibble(
        survival::basehaz(ms_cens, centered = FALSE)
      )
      names(base_hazard) <- c("hazard", "t")

      # Merging and reordering the dataset
      # And estimating the probability of remaining uncensored at each time of event
      res[[arm]] <-
        dat |>
        dplyr::mutate(lin_pred = .env[["lin_pred"]]) |>
        dplyr::left_join(
          base_hazard,
          by = stats::setNames("t", time_start)
        ) |>
        dplyr::mutate(
          hazard = dplyr::coalesce(hazard, 0)
        ) |>
        dplyr::mutate(
          P_uncens = exp(-hazard * exp(lin_pred))
        )
    } else {
      fit_den <- stats::glm(
        pooled_formula,
        data = dat,
        family = stats::binomial(link = "logit")
      )
      p_cens_den <- pmin(
        pmax(stats::predict(fit_den, type = "response"), eps),
        1 - eps
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
          "censoring",
          numerator_data$predictors,
          time_var = time_start
        )
        fit_num <- stats::glm(
          numerator_formula,
          data = numerator_data$data,
          family = stats::binomial(link = "logit")
        )
        p_cens_num <- pmin(
          pmax(stats::predict(fit_num, type = "response"), eps),
          1 - eps
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

predictors <- c(
  "age",
  "sex",
  "emergency",
  "stage",
  "deprivation",
  "charlson",
  "perf"
)

clones_estimated <- estimate_censoring(clones_final, predictors)


weight_cases <- function(clones) {
  arms <- names(clones)

  res <- vector("list", length = length(arms))
  names(res) <- arms

  for (arm in arms) {
    if ("P_uncens_num" %in% names(clones[[arm]])) {
      res[[arm]] <-
        clones[[arm]] |>
        dplyr::mutate(
          weight_Cox = P_uncens_num / P_uncens # Stabilized IPCW
        )
    } else {
      res[[arm]] <-
        clones[[arm]] |>
        dplyr::mutate(
          weight_Cox = 1 / P_uncens # IPCW
        )
    }
  }

  res
}

clones_weighted <- weight_cases(clones_estimated)

#########################
# Step: Main analysis (Emulated trial with Cox weights (Cox regression, Kaplan-Meier, pooled logistic regression))
#########################

backtick_name <- function(x) {
  vapply(
    x,
    function(name) {
      if (make.names(name) == name) {
        name
      } else {
        paste0("`", gsub("`", "\\\\`", name), "`")
      }
    },
    character(1)
  )
}

normalize_predictors <- function(predictors = NULL) {
  if (is.null(predictors)) {
    return(character())
  }
  if (!is.character(predictors)) {
    stop(
      "`predictors` must be NULL or a character vector of column names.",
      call. = FALSE
    )
  }
  if (any(is.na(predictors)) || any(!nzchar(predictors))) {
    stop(
      "`predictors` must contain non-missing, non-empty column names.",
      call. = FALSE
    )
  }

  setdiff(unique(predictors), "arms")
}

normalize_weights <- function(weights, weights_expr, dat) {
  if (identical(weights_expr, quote(NULL))) {
    return(NULL)
  }

  if (is.symbol(weights_expr)) {
    weights_name <- as.character(weights_expr)
    if (weights_name %in% names(dat)) {
      weight_col <- weights_name
    } else {
      weights <- force(weights)
      if (is.null(weights)) {
        return(NULL)
      }
      if (!is.character(weights) || length(weights) != 1L) {
        stop(
          "`weights` must be NULL or a single column name.",
          call. = FALSE
        )
      }
      weight_col <- weights
    }
  } else {
    weights <- force(weights)
    if (is.null(weights)) {
      return(NULL)
    }
    if (!is.character(weights) || length(weights) != 1L) {
      stop(
        "`weights` must be NULL or a single column name.",
        call. = FALSE
      )
    }
    weight_col <- weights
  }

  if (is.na(weight_col) || !nzchar(weight_col)) {
    stop("`weights` must be a non-empty column name.", call. = FALSE)
  }
  if (!weight_col %in% names(dat)) {
    stop("`weights` column not found in data: ", weight_col, call. = FALSE)
  }
  if (!is.numeric(dat[[weight_col]])) {
    stop("`weights` column must be numeric: ", weight_col, call. = FALSE)
  }

  weight_col
}

emul_formula <- function(response, predictors = NULL, cluster = NULL) {
  predictors <- normalize_predictors(predictors)
  terms <- backtick_name(c("arms", predictors))
  if (!is.null(cluster)) {
    terms <- c(terms, paste0("cluster(", backtick_name(cluster), ")"))
  }

  stats::as.formula(
    paste(
      response,
      paste(terms, collapse = " + "),
      sep = " ~ "
    )
  )
}

emul_estimate <- function(
  clones_weighted,
  method = c("Cox", "logistic", "KM"),
  cluster = "id",
  weights = NULL,
  predictors = NULL
) {
  weights_expr <- substitute(weights)
  method <- match.arg(method)
  if (is.data.frame(clones_weighted)) {
    dat <- clones_weighted
    if (!"arms" %in% names(dat)) {
      stop("Data frame input must include an `arms` column.", call. = FALSE)
    }
  } else {
    dat <- dplyr::bind_rows(clones_weighted, .id = "arms")
  }
  if (length(unique(stats::na.omit(dat$arms))) < 2L) {
    stop("`arms` must contain at least two levels.", call. = FALSE)
  }
  predictors <- normalize_predictors(predictors)
  missing_predictors <- setdiff(predictors, names(dat))
  if (length(missing_predictors) > 0L) {
    stop(
      "`predictors` not found in data: ",
      paste(missing_predictors, collapse = ", "),
      call. = FALSE
    )
  }
  if (!cluster %in% names(dat)) {
    stop("`cluster` must be a column in the data.", call. = FALSE)
  }
  weight_col <- normalize_weights(weights, weights_expr, dat)

  cox_formula <- emul_formula(
    "survival::Surv(Tstart, Tstop, outcome)",
    predictors = predictors,
    cluster = cluster
  )
  km_formula <- emul_formula(
    "survival::Surv(Tstart, Tstop, outcome)",
    predictors = predictors
  )
  logistic_formula <- emul_formula(
    "outcome",
    predictors = predictors
  )

  if (method == "KM" && length(predictors) > 0L) {
    message(
      "`predictors` in KM create separate strata; they do not produce ",
      "covariate-adjusted survival curves."
    )
  }

  if (method == "Cox") {
    cox_args <- list(
      formula = cox_formula,
      data = quote(dat),
      robust = TRUE,
      ties = "efron"
    )
    if (!is.null(weight_col)) {
      cox_args$weights <- as.name(weight_col)
    }
    fit_hr <- base::do.call(survival::coxph, cox_args)
    return(fit_hr)
  } else if (method == "logistic") {
    logistic_args <- list(
      formula = logistic_formula,
      data = quote(dat),
      family = stats::binomial(link = "logit")
    )
    if (!is.null(weight_col)) {
      logistic_args$weights <- as.name(weight_col)
    }
    fit_logistic <- base::do.call(stats::glm, logistic_args)
    return(fit_logistic)
  } else if (method == "KM") {
    survfit_args <- list(
      formula = km_formula,
      data = quote(dat)
    )
    if (!is.null(weight_col)) {
      survfit_args$weights <- as.name(weight_col)
    }
    fit_km <- base::do.call(survival::survfit, survfit_args)
    return(fit_km)
  }
}

exp(
  emul_estimate(
    clones_weighted,
    method = "Cox",
    predictors = c("age", "sex")
  )$coefficients
)
emul_estimate(clones_weighted, method = "logistic")
emul_estimate(clones_weighted, method = "KM")


# 1 year survival curve in the surgery group
survminer::ggsurvplot(
  emul_estimate(clones_weighted, method = "KM"),
  data = dplyr::bind_rows(clones_weighted, .id = "arms")
)

#########################
# TO DO: Bootstrapping for confidence interval computation
#########################
emul_estimate_bootstrap <- function(
  clones_weighted,
  method = c("Cox", "logistic"),
  cluster = "id",
  predictors = NULL,
  weights = NULL,
  n_bootstrap = 200
) {
  weights_expr <- substitute(weights)
  method <- match.arg(method)
  dat <- dplyr::bind_rows(clones_weighted, .id = "arms")
  predictors <- normalize_predictors(predictors)
  missing_predictors <- setdiff(predictors, names(dat))
  if (length(missing_predictors) > 0L) {
    stop(
      "`predictors` not found in data: ",
      paste(missing_predictors, collapse = ", "),
      call. = FALSE
    )
  }
  if (!cluster %in% names(dat)) {
    stop("`cluster` must be a column in the data.", call. = FALSE)
  }
  weight_col <- normalize_weights(weights, weights_expr, dat)

  boot_estimates <- vector("numeric", length = n_bootstrap)
  for (i in seq_len(n_bootstrap)) {
    # Resample the data with replacement

    unique_clusters <- unique(dat[[cluster]])
    bootstrap_clusters <- unique_clusters[
      sample.int(
        length(unique_clusters),
        size = length(unique_clusters),
        replace = TRUE
      )
    ]
    bootstrap_sample <- dplyr::bind_rows(
      lapply(seq_along(bootstrap_clusters), function(j) {
        cluster_rows <- dat[
          dat[[cluster]] == bootstrap_clusters[[j]],
          ,
          drop = FALSE
        ]
        cluster_rows$.bootstrap_id <- j
        cluster_rows
      })
    )
    fit <- emul_estimate(
      bootstrap_sample,
      method = method,
      cluster = ".bootstrap_id",
      weights = weight_col,
      predictors = predictors
    )
    coef_values <- stats::coef(fit)
    arm_levels <- levels(factor(bootstrap_sample$arms))
    arm_coef <- intersect(paste0("arms", arm_levels[-1]), names(coef_values))
    if (length(arm_coef) != 1L) {
      stop(
        "Expected exactly one arm coefficient in bootstrap fit.",
        call. = FALSE
      )
    }
    boot_estimates[[i]] <- unname(exp(coef_values[[arm_coef]]))
  }
  lower_ci <- quantile(boot_estimates, probs = 0.025)
  upper_ci <- quantile(boot_estimates, probs = 0.975)
  list(
    ci_lower = lower_ci,
    ci_upper = upper_ci,
    estimates = boot_estimates
  )
}
emul_estimate_bootstrap(
  clones_weighted,
  method = "Cox",
  n_bootstrap = 10,
  predictors = c("age", "sex")
)
emul_estimate_bootstrap(clones_weighted, method = "logistic", n_bootstrap = 10)
