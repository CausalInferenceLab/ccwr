#' Contrast names by effect measure
#'
#' @noRd
contrast_names <- function(method = NULL) {
  names_by_method <- rbind(
    difference = c(
      risk = "RD", RMST = "dRMST", incidence = "IRD", quantile = "qdiff"
    ),
    ratio = c(
      risk = "RR", RMST = "rRMST", incidence = "IRR", quantile = "qratio"
    )
  )
  if (is.null(method)) {
    return(names_by_method)
  }

  names_by_method[, match(method, colnames(names_by_method))]
}

#' Split weighted clones by arm and fit a weighted Kaplan-Meier to each
#'
#' Accepts the same two shapes as [emul_estimate()], a named list of clone data
#' frames or one data frame carrying an arm column. A list is bound in the order
#' of its names, so the first clone is the reference arm; a data frame keeps the
#' levels its arm column already has.
#'
#' Each arm is fitted separately because a cloned patient appears in both arms,
#' and a single fit keyed on the patient id would see one subject with
#' overlapping intervals.
#'
#' @noRd
arm_curves <- function(
  data,
  cluster,
  weights,
  outcome,
  time_start,
  time_stop,
  arm
) {
  bound <- if (is.data.frame(data)) {
    data
  } else {
    do.call(rbind, Map(function(part, label) {
      part[[arm]] <- factor(label, levels = names(data))
      part
    }, data, names(data)))
  }

  dat <- data.frame(
    id = as.character(bound[[cluster]]),
    tstart = bound[[time_start]],
    tstop = bound[[time_stop]],
    event = as.integer(bound[[outcome]]),
    weight = if (is.null(weights)) 1 else bound[[weights]],
    arm = droplevels(as.factor(bound[[arm]]))
  )
  parts <- split(dat, dat$arm)

  curves <- lapply(parts, function(part) {
    fit <- survival::survfit(
      survival::Surv(tstart, tstop, event) ~ 1,
      data = part,
      weights = part$weight,
      id = part$id
    )
    data.frame(
      time = c(0, fit$time),
      surv = c(1, fit$surv),
      cumhaz = c(0, fit$cumhaz)
    )
  })

  list(data = dat, parts = parts, curves = curves, ids = unique(dat$id))
}

#' Warn when an effect measure is read at an unusable time
#'
#' Zero events on or before the time asked for gives an estimate of zero with a
#' zero standard error and a ratio of NaN, and a time past the end of follow-up
#' is read off a curve held flat. Both return something that looks like an
#' answer, which is why they are worth saying out loud. Everything else about
#' the input is left unchecked.
#'
#' @noRd
check_time <- function(dat, at) {
  if (is.null(at)) {
    return(invisible(NULL))
  }

  followup <- max(dat$tstop)
  if (sum(dat$event[dat$tstop <= at]) == 0L) {
    warning(
      "no events on or before ", at, ", so the estimate is 0 and the ratio ",
      "NaN. Follow-up runs to ", signif(followup, 4),
      ", so check the time unit.",
      call. = FALSE
    )
  }
  if (at > followup) {
    warning(
      "time ", at, " is past the end of follow-up at ", signif(followup, 4),
      ", so the estimate is extrapolated from a curve held flat.",
      call. = FALSE
    )
  }

  invisible(at)
}

#' Weighted number at risk at each of several times
#'
#' @noRd
weighted_at_risk <- function(tstart, tstop, weight, times) {
  tail_sum <- function(x, at) {
    ordered <- order(x)
    sorted <- x[ordered]
    running <- rev(cumsum(rev(weight[ordered])))
    k <- findInterval(at, sorted, left.open = TRUE) + 1L
    ifelse(k > length(sorted), 0, running[pmin(k, length(sorted))])
  }

  tail_sum(tstop, times) - tail_sum(tstart, times)
}

#' Read a survival curve at a time point
#'
#' @noRd
curve_value <- function(curve, column, at) {
  curve[[column]][pmax(findInterval(at, curve$time), 1L)]
}

#' Area under a survival curve up to a time point
#'
#' The curve is held flat past its last event, so a time beyond the end of
#' follow-up extrapolates rather than failing.
#'
#' @noRd
restricted_mean <- function(curve, at) {
  area <- c(0, cumsum(curve$surv[-nrow(curve)] * diff(curve$time)))
  i <- pmax(findInterval(at, curve$time), 1L)
  area[i] + curve$surv[i] * (at - curve$time[i])
}

#' First time a survival curve falls to a given level
#'
#' @noRd
curve_quantile <- function(curve, q) {
  hit <- which(curve$surv <= q)
  if (length(hit) == 0L) {
    return(NA_real_)
  }

  curve$time[min(hit)]
}

#' Density of the failure time distribution at a point on the curve
#'
#' A difference quotient of the curve, with Silverman's rule on the observed
#' event times setting the width. The width comes from the event times rather
#' than from the curve's own quartiles because those are missing whenever the
#' curve stops short of them.
#'
#' @noRd
curve_density <- function(curve, at, event_times) {
  bandwidth <- 1.06 * stats::IQR(event_times) / 1.349 *
    length(event_times)^(-1 / 5)
  lower <- max(at - bandwidth, 0)
  upper <- at + bandwidth

  (curve_value(curve, "surv", lower) - curve_value(curve, "surv", upper)) /
    (upper - lower)
}

#' Influence contribution of each subject to a weighted risk-set functional
#'
#' Integrates against the weighted risk set to give
#' `psi_i = sum_j (k_j / Y_j) * w_i(t_j) * [dN_i(t_j) - Y_i(t_j) * h_j]`, with
#' `h_j` the weighted hazard increment and the kernel `k_j` selecting the
#' effect measure. The result spans every subject, and is zero outside the arm
#' supplied, so that contrasts can be taken by subtracting columns.
#'
#' @noRd
influence_contributions <- function(dat, ids, mask, upto, kernel = NULL) {
  contributions <- stats::setNames(numeric(length(ids)), ids)
  times <- sort(unique(dat$tstop[mask & dat$tstop <= upto]))
  if (length(times) == 0L) {
    return(contributions)
  }

  at_risk <- weighted_at_risk(dat$tstart, dat$tstop, dat$weight, times)
  events <- as.numeric(tapply(
    dat$weight[mask],
    factor(dat$tstop[mask], levels = times),
    sum
  ))
  events[is.na(events)] <- 0
  weight_j <- if (is.null(kernel)) rep(1, length(times)) else kernel(times)

  cumulative <- cumsum(weight_j * (events / at_risk) / at_risk)
  cumulative_at <- function(at) {
    c(0, cumulative)[findInterval(at, times) + 1L]
  }

  counting <- numeric(nrow(dat))
  rows <- which(mask & dat$tstop <= upto)
  j <- match(dat$tstop[rows], times)
  counting[rows] <- weight_j[j] * dat$weight[rows] / at_risk[j]

  compensator <- dat$weight * (
    cumulative_at(pmin(dat$tstop, upto)) - cumulative_at(pmin(dat$tstart, upto))
  )
  by_subject <- rowsum(counting - compensator, dat$id, reorder = FALSE)
  contributions[rownames(by_subject)] <- by_subject[, 1L]
  contributions
}

#' Name per-arm estimates and append the two-arm contrasts
#'
#' Both estimation layers take their names from here, so an estimate and its
#' standard error cannot come out labelled differently.
#'
#' @noRd
label_estimates <- function(value, method, arms, q) {
  prefix <- c(
    survival = "S_",
    cumhaz = "H_",
    risk = "risk_",
    RMST = "RMST_",
    incidence = "rate_",
    quantile = paste0("q", q, "_")
  )[method]
  names(value) <- paste0(prefix, arms)

  contrast <- contrast_names(method)
  if (length(value) != 2L || is.na(contrast[["difference"]])) {
    return(value)
  }

  c(
    value,
    stats::setNames(
      c(value[2L] - value[1L], value[2L] / value[1L]),
      contrast
    )
  )
}

#' Append the influence contributions of the two-arm contrasts
#'
#' The difference is exact and the ratio is the delta method, both taken on the
#' subjects shared between the arms. Under cloning that is what nets a
#' patient's within-patient correlation out of the contrast.
#'
#' @noRd
contrast_influence <- function(psi, value, method) {
  contrast <- contrast_names(method)
  if (ncol(psi) != 2L || is.na(contrast[["difference"]])) {
    return(psi)
  }

  cbind(
    psi,
    psi[, 2L] - psi[, 1L],
    (psi[, 2L] - value[2L] / value[1L] * psi[, 1L]) / value[1L]
  )
}
