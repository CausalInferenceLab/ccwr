#' Estimate a marginal effect measure from weighted clones
#'
#' @param data A data frame of weighted clones in (start, stop] long format, as
#'   returned by [weight_cases()].
#' @param method Effect measure: `"survival"` for S(t) by weighted
#'   Kaplan-Meier, `"cumhaz"` for H(t) by weighted Nelson-Aalen, `"risk"` for
#'   1 - S(t), `"RMST"` for restricted mean survival time, `"incidence"` for the
#'   person-time rate, or `"quantile"` for the first time S(t) falls to `q`.
#' @param cluster Column name identifying the subject. Under cloning this is the
#'   patient id, not the clone id.
#' @param weights Weight column name, or `NULL` for an unweighted analysis.
#' @param outcome Column name for the 0/1 outcome indicator, 1 being the event.
#' @param time_start Column name for interval start time.
#' @param time_stop Column name for interval stop time.
#' @param arm Column name for treatment arm.
#' @param horizon Time point for `"survival"`, `"cumhaz"` and `"risk"`.
#' @param tau Restriction time for `"RMST"`.
#' @param q Survival level for `"quantile"`; `0.5` gives median survival.
#' @param per Person-time denominator for `"incidence"`.
#' @param conf_level Coverage of the confidence interval.
#'
#' @returns A data frame with columns `est`, `se`, `lcl`, `ucl` and `p`, one row
#'   per arm followed by the difference and ratio rows.
#'
#' @details
#' Complements [emul_estimate()]: that function fits a model and returns the
#' model object, this one returns the marginal effect measure itself with an
#' interval. Both take the output of [weight_cases()] and the same
#' column-naming arguments.
#'
#' The work runs in three layers, each callable on its own:
#' [emul_effect_point()] for the estimate, [emul_effect_variance()] for the
#' robust sandwich standard error, and [emul_confint()] for the interval and
#' test. One variance object can therefore serve several intervals, and either
#' layer can be replaced without touching the other.
#'
#' Little is validated, so three things are assumed and none of them raise
#' an error when they are wrong, only a plausible number. `outcome` must be
#' coded 0/1 with 1 the event, so survival's 1 = censored, 2 = event coding has
#' to be recoded first. `cluster` must be the patient id, since a clone id makes
#' a patient's two clones look independent and shrinks every contrast standard
#' error. And intervals of one subject within one arm must be disjoint.
#'
#' The one thing that is checked is the time the measure is read at: a
#' `horizon` or `tau` with no events before it, or past the end of follow-up,
#' warns rather than quietly returning zero or an extrapolation.
#'
#' Weights are treated as known: no standard error here accounts for having
#' estimated the weight model, which is conservative for a correctly specified
#' IPTW. A quantile exists only where S(t) reaches `q` inside follow-up; where
#' it does not the estimate is `NA`, and where it does but little of the risk
#' set is left beneath it, the interval covers less than it claims.
#'
#' @export
#' @examples
#' set.seed(1)
#' n <- 200
#' treated <- rbinom(n, 1, 0.5)
#' end <- pmin(rexp(n, ifelse(treated == 1, 0.10, 0.17)), runif(n, 2, 9))
#' died <- as.integer(end < 9)
#' rows <- pmax(ceiling(end), 1)
#' dat <- data.frame(
#'   id = rep(seq_len(n), rows),
#'   Tstart = sequence(rows) - 1,
#'   Tstop = pmin(sequence(rows), end[rep(seq_len(n), rows)]),
#'   arms = factor(treated[rep(seq_len(n), rows)], 0:1, c("Control", "Surgery")),
#'   weight_Cox = runif(sum(rows), 0.5, 2)
#' )
#' dat$outcome <- as.integer(dat$Tstop == end[dat$id] & died[dat$id] == 1)
#'
#' emul_effect(dat, "risk", weights = "weight_Cox", horizon = 5)
#' emul_effect(dat, "RMST", weights = "weight_Cox", tau = 6)
#' emul_effect(dat, "incidence", weights = "weight_Cox", per = 100)
#'
#' est <- emul_effect_point(dat, "risk", weights = "weight_Cox", horizon = 5)
#' v <- emul_effect_variance(dat, "risk", weights = "weight_Cox", horizon = 5)
#' emul_confint(est, v$se)
emul_effect <- function(
  data,
  method = c("survival", "cumhaz", "risk", "RMST", "incidence", "quantile"),
  cluster = "id",
  weights = NULL,
  outcome = "outcome",
  time_start = "Tstart",
  time_stop = "Tstop",
  arm = "arms",
  horizon = NULL,
  tau = NULL,
  q = 0.5,
  per = 1000,
  conf_level = 0.95
) {
  method <- match.arg(method)

  estimate <- emul_effect_point(
    data, method, cluster, weights, outcome, time_start, time_stop, arm,
    horizon, tau, q, per
  )
  variance <- emul_effect_variance(
    data, method, cluster, weights, outcome, time_start, time_stop, arm,
    horizon, tau, q, per
  )

  emul_confint(estimate, variance$se, conf_level)
}

#' Point estimate of a marginal effect measure
#'
#' @inheritParams emul_effect
#'
#' @returns A named numeric vector, one entry per arm. With exactly two arms the
#'   difference and the ratio are appended, always as arm 2 against arm 1.
#' @export
emul_effect_point <- function(
  data,
  method = c("survival", "cumhaz", "risk", "RMST", "incidence", "quantile"),
  cluster = "id",
  weights = NULL,
  outcome = "outcome",
  time_start = "Tstart",
  time_stop = "Tstop",
  arm = "arms",
  horizon = NULL,
  tau = NULL,
  q = 0.5,
  per = 1000
) {
  method <- match.arg(method)
  fitted <- arm_curves(
    data, cluster, weights, outcome, time_start, time_stop, arm
  )
  curves <- fitted$curves
  check_time(fitted$data, if (identical(method, "RMST")) tau else horizon)

  value <- switch(
    method,
    survival = vapply(curves, curve_value, numeric(1), "surv", horizon),
    cumhaz = vapply(curves, curve_value, numeric(1), "cumhaz", horizon),
    risk = 1 - vapply(curves, curve_value, numeric(1), "surv", horizon),
    RMST = vapply(curves, restricted_mean, numeric(1), tau),
    incidence = vapply(fitted$parts, function(part) {
      per * sum(part$weight * part$event) /
        sum(part$weight * (part$tstop - part$tstart))
    }, numeric(1)),
    quantile = vapply(curves, curve_quantile, numeric(1), q)
  )

  label_estimates(value, method, names(curves), q)
}

#' Robust sandwich variance of a marginal effect measure
#'
#' Standard errors in closed form from subject-level influence functions,
#' `Var(theta) = sum_i psi_i^2`. Contrasts are exact rather than assumed
#' independent, because the influence contributions of both arms live on the
#' same subject.
#'
#' @inheritParams emul_effect
#'
#' @returns A list with `se`, whose names match [emul_effect_point()], and
#'   `influence`, the subjects by estimands matrix of influence contributions.
#'   A custom contrast is column arithmetic on `influence`.
#' @export
emul_effect_variance <- function(
  data,
  method = c("survival", "cumhaz", "risk", "RMST", "incidence", "quantile"),
  cluster = "id",
  weights = NULL,
  outcome = "outcome",
  time_start = "Tstart",
  time_stop = "Tstop",
  arm = "arms",
  horizon = NULL,
  tau = NULL,
  q = 0.5,
  per = 1000
) {
  method <- match.arg(method)
  fitted <- arm_curves(
    data, cluster, weights, outcome, time_start, time_stop, arm
  )
  parts <- fitted$parts
  curves <- fitted$curves
  ids <- fitted$ids

  hazard <- function(at) {
    vapply(parts, function(part) {
      influence_contributions(part, ids, part$event == 1, at)
    }, numeric(length(ids)))
  }

  block <- switch(
    method,
    survival = {
      surv <- vapply(curves, curve_value, numeric(1), "surv", horizon)
      list(psi = sweep(hazard(horizon), 2L, -surv, "*"), value = surv)
    },
    cumhaz = list(
      psi = hazard(horizon),
      value = vapply(curves, curve_value, numeric(1), "cumhaz", horizon)
    ),
    risk = {
      surv <- vapply(curves, curve_value, numeric(1), "surv", horizon)
      list(psi = sweep(hazard(horizon), 2L, surv, "*"), value = 1 - surv)
    },
    RMST = {
      rmst <- vapply(curves, restricted_mean, numeric(1), tau)
      list(
        psi = vapply(seq_along(parts), function(i) {
          -influence_contributions(
            parts[[i]], ids, parts[[i]]$event == 1, tau,
            function(times) rmst[i] - restricted_mean(curves[[i]], times)
          )
        }, numeric(length(ids))),
        value = rmst
      )
    },
    incidence = list(
      psi = vapply(parts, function(part) {
        events <- rowsum(part$weight * part$event, part$id, reorder = FALSE)
        exposure <- rowsum(
          part$weight * (part$tstop - part$tstart), part$id, reorder = FALSE
        )
        total <- sum(exposure)
        contributions <- stats::setNames(numeric(length(ids)), ids)
        contributions[rownames(events)] <- per *
          (events[, 1L] - sum(events) / total * exposure[, 1L]) / total
        contributions
      }, numeric(length(ids))),
      value = vapply(parts, function(part) {
        per * sum(part$weight * part$event) /
          sum(part$weight * (part$tstop - part$tstart))
      }, numeric(1))
    ),
    quantile = {
      theta <- vapply(curves, curve_quantile, numeric(1), q)
      list(
        psi = vapply(seq_along(parts), function(i) {
          part <- parts[[i]]
          density <- curve_density(
            curves[[i]], theta[i], part$tstop[part$event == 1]
          )
          -curve_value(curves[[i]], "surv", theta[i]) *
            influence_contributions(part, ids, part$event == 1, theta[i]) /
            density
        }, numeric(length(ids))),
        value = theta
      )
    }
  )

  psi <- contrast_influence(block$psi, block$value, method)
  colnames(psi) <- names(
    label_estimates(block$value, method, names(curves), q)
  )

  list(se = sqrt(colSums(psi^2)), influence = psi)
}

#' Confidence interval and Wald test for a marginal effect measure
#'
#' Ratios are put on the log scale and tested against 1, differences are tested
#' against 0. A per-arm absolute quantity has no meaningful null value, so its
#' p is `NA` rather than a test against zero.
#'
#' @param estimate Named numeric vector from [emul_effect_point()].
#' @param se Named numeric vector of standard errors, the `se` element of
#'   [emul_effect_variance()].
#' @param conf_level Coverage of the confidence interval.
#'
#' @returns A data frame with columns `est`, `se`, `lcl`, `ucl` and `p`.
#' @export
emul_confint <- function(estimate, se, conf_level = 0.95) {
  terms <- names(estimate)
  se <- se[terms]
  contrast <- contrast_names()
  ratio <- terms %in% contrast["ratio", ]
  difference <- terms %in% contrast["difference", ]
  z <- stats::qnorm(1 - (1 - conf_level) / 2)

  lcl <- estimate - z * se
  ucl <- estimate + z * se
  lcl[ratio] <- exp(log(estimate[ratio]) - z * se[ratio] / estimate[ratio])
  ucl[ratio] <- exp(log(estimate[ratio]) + z * se[ratio] / estimate[ratio])

  null <- rep(NA_real_, length(terms))
  null[difference] <- 0
  null[ratio] <- 1
  statistic <- (estimate - null) / se
  statistic[ratio] <- log(estimate[ratio]) / (se[ratio] / estimate[ratio])

  data.frame(
    est = unname(estimate),
    se = unname(se),
    lcl = unname(lcl),
    ucl = unname(ucl),
    p = unname(2 * stats::pnorm(-abs(statistic))),
    row.names = terms
  )
}
