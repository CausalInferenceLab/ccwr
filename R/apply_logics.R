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
#' data(lungcancer)
#' arms <- c("Control", "Surgery")
#' clones <- clone_arms(lungcancer, arms)
#' policies <- create_policy_A(
#'   arms,
#'   treatment = "surgery",
#'   time_to_treatment = "timetosurgery",
#'   grace_period = 182.62,
#'   outcome = "death",
#'   followup = "fup_obs",
#'   clone_outcome = "outcome",
#'   clone_followup = "fup"
#' )
#' clones_policy <- apply_logics(clones, policies)
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
#' arms <- c("Control", "Surgery")
#' policies <- create_policy_A(
#'   arms,
#'   treatment = "surgery",
#'   time_to_treatment = "timetosurgery",
#'   grace_period = 182.62,
#'   outcome = "death",
#'   followup = "fup_obs",
#'   clone_outcome = "outcome",
#'   clone_followup = "fup"
#' )
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
#' arms <- c("Control", "Surgery")
#' censoring_logics <- create_censoring_logics_A(
#'   arms,
#'   treatment = "surgery",
#'   time_to_treatment = "timetosurgery",
#'   grace_period = 182.62,
#'   followup = "fup_obs",
#'   clone_censoring = "censoring",
#'   clone_uncensored_followup = "fup_uncensored"
#' )
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
