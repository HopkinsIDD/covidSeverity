#'
#' US population by age category and GEOID
#'
#' A matrix containing the proportion of the population by 10 year age
#' categories in the US
#'
#' @format A matrix where rows are GEOIDs and columns are age groups
#' @docType data
#' @name US_age_geoid_pct
#' @usage data(US_age_geoid_pct)
NULL

#'
#' US population by age category and GEOID
#'
#' A matrix containing the total population by 10 year age categories in the US
#'
#' @format A matrix where rows are GEOIDs and columns are age groups
#' @docType data
#' @name US_age_geoid_pop
#' @usage data(US_age_geoid_pop)
NULL

#'
#' Raw age-specific estimates
#'
#' A data frame with age-specific parameter estimates from a variety of sources
#'
#' @format tbl with the following columns
#' \itemize{
#'   \item \code{param} character string of a conditional probability with the format p_x_y that is equivalent to P(X | Y)
#'   \item \code{USE_pipeline} logical, is the estimate used in the model pipeline?
#'   \item \code{ageL} numeric, lower bound of age category
#'   \item \code{ageR} numeric, upper bound of age category
#'   \item \code{X} integer, numerator value of conditional probability (e.g., the number of symptomatic cases in P(symptomatic | infection))
#'   \item \code{N} integer, denominator value of conditional probability (e.g., the number of infections in P(symptomatic | infection))
#'   \item \code{p} numeric, the proportion estimate for the conditional probability and age category.Is equal to X/N when both are known.
#'   \item \code{pL} numeric, lower bound of the conditional probability
#'   \item \code{pR} numeric, upper bound of the conditional probability
#'   \item \code{desc} character string, description of the observation
#'   \item \code{data_date} date, when was the data collected
#'   \item \code{publish_date} date, when was the source published
#'   \item \code{source} character string, where is the data from
#' }
#' @docType data
#' @name raw_age_estimates
#' @usage data(raw_age_estimates)
NULL

#'
#' Age-adjusted parameters for each US GEOID
#'
#' A data frame with the median conditional probabilities for each parameter for each GEOID
#'
#' @format tbl, where rows are GEOIDs and columns are parameters:
#' \itemize{
#'   \item \code{geoid} the GEOID of interest
#'   \item \code{p_symp_inf} probability of being symptomatic given infected
#'   \item \code{p_death_symp} probability of death given symptomatic
#'   \item \code{p_hosp_symp} probability of hospitalization given symptomatic
#'   \item \code{p_icu_hosp} probability of going to the ICU given hospitalized
#'   \item \code{p_vent_icu} probability of needing invasive mechanized ventilation given ICU
#'   \item \code{rr_death_symp} the relative risk of death given symptomatic (vs. US average)
#'   \item \code{death_symp_overall} US average probability of death given symptomatic
#'   \item \code{rr_hosp_symp} the relative risk of hospitalization given symptomatic (vs. US average)
#'   \item \code{hosp_symp_overall} US average probability of hospitalization given symptomatic
#'   \item \code{rr_death_inf} the relative risk of death given infection (vs. US average)
#'   \item \code{death_inf_overall} US average probability of death given infection
#'   \item \code{rr_hosp_inf} the relative risk of hospitalization given infection (vs. US average)
#'   \item \code{hosp_inf_overall} US average probability of hospitalization given infection
#' }
#' @docType data
#' @name US_geoid_params
#' @usage data(US_geoid_params)
NULL
