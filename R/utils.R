.assert_data_frame <- function(data) {
  if (!is.data.frame(data)) {
    stop("`data` must be a data frame.", call. = FALSE)
  }
}

.assert_required_columns <- function(data, columns) {
  missing_columns <- setdiff(columns, names(data))

  if (length(missing_columns) > 0) {
    stop(
      "Missing required columns: ",
      paste(missing_columns, collapse = ", "),
      call. = FALSE
    )
  }
}

.assert_named_data_frame_list <- function(data, arg = "clones") {
  if (!is.list(data) || is.null(names(data)) || any(names(data) == "")) {
    stop("`", arg, "` must be a named list of data frames.", call. = FALSE)
  }

  for (name in names(data)) {
    if (!is.data.frame(data[[name]])) {
      stop("Every element of `", arg, "` must be a data frame.", call. = FALSE)
    }
  }
}

.assert_clone_columns <- function(clones, columns) {
  .assert_named_data_frame_list(clones)

  for (arm in names(clones)) {
    tryCatch(
      .assert_required_columns(clones[[arm]], columns),
      error = function(cnd) {
        stop(
          "In clone arm `", arm, "`: ",
          conditionMessage(cnd),
          call. = FALSE
        )
      }
    )
  }
}

.assert_probability_floor <- function(eps) {
  if (!is.numeric(eps) || length(eps) != 1L || is.na(eps) || eps <= 0 || eps >= 1) {
    stop("`eps` must be a single number between 0 and 1.", call. = FALSE)
  }
}

.clamp_probability <- function(x, eps = 1e-6) {
  .assert_probability_floor(eps)
  if (!is.numeric(x)) {
    stop("Probabilities must be numeric.", call. = FALSE)
  }

  pmin(pmax(x, eps), 1 - eps)
}

.assert_positive_integer <- function(x, arg) {
  if (
    !is.numeric(x) ||
      length(x) != 1L ||
      is.na(x) ||
      x < 1 ||
      x != as.integer(x)
  ) {
    stop("`", arg, "` must be a positive integer.", call. = FALSE)
  }
}

.assert_conf_level <- function(conf_level) {
  if (
    !is.numeric(conf_level) ||
      length(conf_level) != 1L ||
      is.na(conf_level) ||
      conf_level <= 0 ||
      conf_level >= 1
  ) {
    stop("`conf_level` must be a single number between 0 and 1.", call. = FALSE)
  }
}
