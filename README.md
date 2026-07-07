# clonecensorweighting

`clonecensorweighting` is an R package for building reproducible
clone-censor-weighting workflows for target trial emulation.

This repository is meant to be easy to use from a fresh GitHub clone and easy
to maintain as a shared collaboration project. The recommended setup uses:

- `rig` to install and switch to the project R version
- `renv` to restore the same package environment for every collaborator
- GitHub Actions to check reproducibility and package health

## Quick start

### 1. Clone the repository

```sh
git clone https://github.com/CausalInferenceLab/clonecensorweighting.git
cd clonecensorweighting
```

### 2. Use `rig` to install R 4.4.2

This project is pinned to R `4.4.2` for reproducibility.

```sh
rig add 4.4.2
rig default 4.4.2
```

If you already have R `4.4.2`, you can skip `rig add 4.4.2`.

### 3. Restore the project package library with `renv`

Open R in the project directory and run:

```r
install.packages("renv")
renv::restore()
```

`renv::restore()` installs the package versions recorded in `renv.lock`, so
everyone works with the same dependency set.

Note: this repository includes a `.Rprofile` file that activates `renv`
automatically when you open the project in R.

### 4. Install the package locally

From the project root:

```sh
R CMD INSTALL .
```

Then in R:

```r
library(clonecensorweighting)
```

## First example

```r
library(clonecensorweighting)

data(lungcancer)

arms <- c("Control", "Surgery")
clones <- clone_arms(lungcancer, arms)

policies <- create_policy_A(
  arms,
  treatment = "surgery",
  time_to_treatment = "timetosurgery",
  grace_period = 182.62,
  outcome = "death",
  followup = "fup_obs",
  clone_outcome = "outcome",
  clone_followup = "fup"
)
clones_policy <- apply_logics(clones, policies)

censoring_logics <- create_censoring_logics_A(
  arms,
  treatment = "surgery",
  time_to_treatment = "timetosurgery",
  grace_period = 182.62,
  followup = "fup_obs",
  clone_censoring = "censoring",
  clone_uncensored_followup = "fup_uncensored"
)
clones_censored <- apply_logics(clones_policy, censoring_logics)

clones_final <- create_final_data(
  clones_censored,
  clone_followup = "fup",
  clone_outcome = "outcome",
  clone_censoring = "censoring",
  col_ids = "id"
)

clones_estimated <- estimate_censoring(
  clones_final,
  predictors = c("age", "sex"),
  method = "pooled_logit"
)
clones_weighted <- weight_cases(clones_estimated)

fit <- emul_estimate(
  clones_weighted,
  method = "Cox",
  weights = "weight_Cox",
  predictors = c("age", "sex")
)

exp(stats::coef(fit))
```

The full process is:

1. Clone each patient into the target trial arms with `clone_arms()`.
2. Define and apply treatment-policy logic with `create_policy_A()` and
   `apply_logics()`.
3. Define and apply artificial censoring logic with
   `create_censoring_logics_A()` and `apply_logics()`.
4. Expand cloned observations into long-form interval data with
   `create_final_data()`.
5. Estimate censoring probabilities with `estimate_censoring()`.
6. Add inverse probability of censoring weights with `weight_cases()`.
7. Estimate the emulated trial effect with `emul_estimate()`.
8. Use `emul_estimate_bootstrap()` when bootstrap confidence intervals are
   needed.

## What the package currently provides

The package is still an early, lightweight foundation for future work. Right
now it includes:

- `read_trial_data()` to read trial-style CSV data into a tibble
- `clone_censor_weighting()` to create a starter cloned dataset across regimes
- `make_surv_response()` to build a `survival::Surv()` response object
- `clone_arms()` to duplicate observations across treatment strategies
- `create_policy_A()` and `create_censoring_logics_A()` to generate example
  treatment-policy and artificial-censoring logic for the lung cancer scenario
- `apply_logics()` to apply policy or censoring logic to cloned data
- `create_final_data()` to create long-form interval data for censoring models
- `estimate_censoring()` and `weight_cases()` to estimate censoring
  probabilities and add IPC weights
- `emul_estimate()` and `emul_estimate_bootstrap()` to estimate treatment
  effects and bootstrap confidence intervals

## Working together on this repository

For collaborative work, the safest pattern is:

1. Clone the repository.
2. Switch to R `4.4.2` with `rig`.
3. Run `renv::restore()`.
4. Make your changes in a branch.
5. Run checks before opening a pull request.

Useful commands:

```sh
R CMD INSTALL .
R CMD check --no-manual .
```

In R:

```r
testthat::test_local()
```

If you add, remove, or upgrade dependencies, update the lockfile from R:

```r
renv::settings$snapshot.type("explicit")
renv::status(dev = TRUE)
renv::snapshot(dev = TRUE)
```

Please commit both code changes and the updated `renv.lock` when dependency
changes are intentional.

## Continuous integration

This repository includes two GitHub Actions workflows:

- `check-reproducible.yaml` runs `R CMD check` with pinned R `4.4.2` and `renv`
- `check-latest.yaml` runs broader checks across operating systems and R versions

Together, these workflows help keep the project reproducible for collaborators
and stable for future users.
