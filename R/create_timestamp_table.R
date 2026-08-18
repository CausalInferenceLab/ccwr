#' Create a timestamp lookup table
#'
#' @noRd
create_timestamp_table <- function(clones, clone_followup) {
  .assert_clone_columns(clones, clone_followup)

  timestamps <- unlist(
    lapply(clones, `[[`, clone_followup),
    use.names = FALSE
  )
  t_events <- sort(unique(timestamps))

  dplyr::tibble(tevent = t_events, ID_t = seq_along(t_events))
}
