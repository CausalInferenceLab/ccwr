# clonecensorweighting 0.0.0.9000

## Refactoring

* Split the final censoring-probability training data workflow into focused
  functions.
* Kept `create_final_data()` as the public interface for constructing long-form
  clone data.
* Moved `create_timestamp_table()` and `split_at_timestamp()` to internal helper
  status and removed them from the package exports.
* Removed duplicated reshape code and script-style execution leftovers from the
  package R files.
* Added clone-list and clone-column validation helpers for clearer errors when
  required variables are missing.
* Updated `create_final_data()` internals to use explicit character-vector
  column selection and joins.
* Added `estimate_censoring()` as the public interface for Cox, pooled-logit,
  and stabilized pooled-logit censoring probability estimation.
* Added `weight_cases()` as the public interface for unstabilized and stabilized
  inverse probability of censoring weights.
* Added `emul_estimate()` as the public interface for Cox, logistic, and
  Kaplan-Meier analyses of the emulated trial.
* Added `emul_estimate_bootstrap()` as the public interface for bootstrap
  confidence intervals for two-arm Cox and logistic analyses.
* Kept formula construction, cumulative uncensoring, baseline predictor
  extraction, predictor/weight normalization, arm binding, and coefficient
  extraction as internal helpers.

## Documentation

* Removed generated help pages for internal timestamp and splitting helpers.
* Regenerated the `create_final_data()` help source reference after moving the
  function to its own file.
* Added generated help pages for censoring estimation, case weighting,
  emulated-trial estimation, and bootstrap estimation.
* Added runnable examples for public workflow functions using the bundled
  `lungcancer` data where appropriate.
* Updated the README with a full `lungcancer` clone-censor-weighting workflow.

## Testing

* Added focused tests for timestamp table creation, interval splitting, final
  data construction, internal helper export status, and missing-column errors.
* Added focused tests for censoring probability estimation, IPC weighting,
  emulated-trial model fitting, bootstrap output, and internal helper export
  status.
* Added an integration test that runs the bundled `lungcancer` data through
  cloning, censoring, interval expansion, censoring probability estimation,
  weighting, emulated-trial estimation, and bootstrap estimation.
