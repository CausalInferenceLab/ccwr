#' Read trial-style data from CSV
#'
#' @param file Path to a CSV file.
#' @param show_col_types Passed to [readr::read_csv()].
#'
#' @return A tibble.
#' @export
#' @examples
#' csv_file <- tempfile(fileext = ".csv")
#' writeLines(
#'   c(
#'     "id,follow_up,event,treatment",
#'     "1,10,1,A",
#'     "2,12,0,B"
#'   ),
#'   con = csv_file
#' )
#' trial_data <- read_trial_data(csv_file)
read_trial_data <- function(file, show_col_types = FALSE) {
  tibble::as_tibble(
    readr::read_csv(file, show_col_types = show_col_types)
  )
}
