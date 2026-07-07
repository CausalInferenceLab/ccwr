#' Simulated lung cancer patients data by Maringe et al. (2020)
#'
#' The dataset is from supplmentary files of [Maringe et al. (2020)](https://doi.org/10.1093/ije/dyaa057).
#' The dataset is a set of 200 simulated lung cancer patients. 
#' These patients are followed up for a year following their cancer diagnosis:
#' 106 of them received surgery within six months of their diagnosis and 
#' 48 died in the year.
#'
#' @format ## `lungcancer`
#' A data frame with 200 rows and 12 columns:
#' \describe{
#'   \item{id}{patient identifier}
#'   \item{fup_obs}{observed follow-up time (time to death or 1 year if censored alive)}
#'   \item{death}{observed event of interest (all-cause death) 1: dead, 0:alive}
#'   \item{timetosurgery}{time to surgery (NA if no surgery)}
#'   \item{surgery}{observed treatment 1 if the patient received surgery within 6 month, 0 otherwise}
#'   \item{age}{age at diagnosis}
#'   \item{sex}{patient's sex}
#'   \item{perf}{performance status at diagnosis}
#'   \item{stage}{stage at diagnosis}
#'   \item{deprivation}{deprivation score}
#'   \item{charlson}{Charlson's comorbidity index}
#'   \item{emergency}{route to diagnosis}
#' }
#' @source <https://doi.org/10.1093/ije/dyaa057>
#' @examples
#' data(lungcancer)
"lungcancer"


#' 13 types of patient records in Maringe et al. (2020)
#'
#' The dataset is from Figure 2 in [Maringe et al. (2020)](https://doi.org/10.1093/ije/dyaa057).
#' The dataset illustrates all possible censoring mechanisms with 13 types of 
#' patients records that could be seen in the cancer registry data, when 
#' allowing at most one treatment (surgery in this example) for each patient
#' and outcome event (death in this example) happens at most once for each
#' patient.
#'
#' @format ## `patients13`
#' A data frame with 13 rows and 5 columns:
#' \describe{
#'   \item{id}{patient identifier}
#'   \item{surgery}{observed treatment 1 if the patient received surgery within 6 month, 0 otherwise}
#'   \item{timetosurgery}{time to surgery (NA if no surgery)}
#'   \item{death}{observed event at the latest follow-up, 1: dead, 0: alive}
#'   \item{followup}{observed follow-up time (time to death or time to latest followup with 1 year (365 days) at maximum)}
#' }
#' @source <https://doi.org/10.1093/ije/dyaa057>
#' @examples
#' data(patients13)
"patients13"