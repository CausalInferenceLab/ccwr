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
#' data(lungcancer)
#' clones <- clone_arms(lungcancer, c("Control", "Surgery"))
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
