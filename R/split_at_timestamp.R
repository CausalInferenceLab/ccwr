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
#' @noRd
split_at_timestamp <- function(
  clones,
  clone_followup,
  t_events,
  event,
  timestamp_start = "Tstart",
  id = "ID"
) {
  .assert_clone_columns(clones, c(clone_followup, event))

  arms <- names(clones)

  res <- vector("list", length = length(arms))
  names(res) <- arms
  for (arm in arms) {
    res[[arm]] <-
      survival::survSplit(
        data = clones[[arm]],
        cut = t_events,
        end = clone_followup,
        start = timestamp_start,
        event = event,
        id = id
      )
  }

  res
}
