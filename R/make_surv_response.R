#' Construct a survival response
#'
#' @param data A data frame with follow-up and event columns.
#' @param follow_up The name of the follow-up time column.
#' @param event The name of the event indicator column.
#'
#' @return An object of class `"Surv"`.
#' @export
#' @examples
#' data(lungcancer)
#' surv_response <- make_surv_response(
#'   lungcancer,
#'   follow_up = "fup_obs",
#'   event = "death"
#' )
make_surv_response <- function(data, follow_up, event) {
  .assert_data_frame(data)
  .assert_required_columns(data, c(follow_up, event))

  survival::Surv(time = data[[follow_up]], event = data[[event]])
}
