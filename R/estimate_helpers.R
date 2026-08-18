#' Backtick a column name when needed for formulas
#'
#' @noRd
backtick_name <- function(x) {
  if (length(x) == 0L) {
    return(character())
  }

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

#' Normalize predictor column names
#'
#' @noRd
normalize_predictors <- function(predictors = NULL, exclude = NULL) {
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

  setdiff(unique(predictors), exclude)
}

#' Normalize a weights argument to a data column
#'
#' @noRd
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

#' Create an emulated-trial model formula
#'
#' @noRd
emul_formula <- function(
  response,
  predictors = NULL,
  cluster = NULL,
  arm = "arms"
) {
  predictors <- normalize_predictors(predictors, exclude = arm)
  terms <- backtick_name(c(arm, predictors))
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

#' Bind clone arms into one data frame
#'
#' @noRd
bind_clone_arms <- function(clones, arm = "arms") {
  if (is.data.frame(clones)) {
    dat <- tibble::as_tibble(clones)
    if (!arm %in% names(dat)) {
      stop("Data frame input must include an `", arm, "` column.", call. = FALSE)
    }
  } else {
    .assert_named_data_frame_list(clones)
    arm_conflicts <- vapply(clones, function(x) arm %in% names(x), logical(1))
    if (any(arm_conflicts)) {
      stop(
        "Clone data frames must not already include the `", arm, "` column.",
        call. = FALSE
      )
    }
    dat <- dplyr::bind_rows(clones, .id = arm)
  }

  dat[[arm]] <- factor(dat[[arm]], levels = unique(dat[[arm]]))
  dat
}

#' Add a temporary analysis weight column
#'
#' @noRd
add_analysis_weights <- function(dat, weight_col) {
  if (is.null(weight_col)) {
    return(list(data = dat, weight = NULL))
  }

  weight_name <- ".analysis_weight"
  while (weight_name %in% names(dat)) {
    weight_name <- paste0(".", weight_name)
  }
  dat[[weight_name]] <- dat[[weight_col]]

  list(data = dat, weight = weight_name)
}

#' Find the arm coefficient in a two-arm model
#'
#' @noRd
find_arm_coefficient <- function(coef_values, arm, arm_level) {
  candidates <- c(
    paste0(arm, arm_level),
    paste0(backtick_name(arm), arm_level)
  )
  arm_coef <- intersect(candidates, names(coef_values))

  if (length(arm_coef) != 1L) {
    stop(
      "Expected exactly one arm coefficient in the fitted model.",
      call. = FALSE
    )
  }

  arm_coef
}
