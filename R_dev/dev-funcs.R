##' Function to standardize a given probability by age for each
##' GEOID - assuming a single vector of age-specific probabilities
##'
##' @param age_county_pop age specific population by county
##' @param p_vec age specific probability of event
##' @param var_name variable name of standardized pop
##'
##' @return p, standardized probability for each GEOID
##'
##' @import tidyverse
##'
##' @export

get_p_standardized <- function(age_county_pop, p_vec, var_name){

  p_age <- age_county_pop %>%
    select(GEOID, cat_l, page) %>%
    pivot_wider(names_from = cat_l, values_from = page)
  GEOID <- p_age[,1]
  p_age <- as.matrix(p_age[,-1])

  p_tmp <- sweep(p_age, 2, p_vec, "*")
  p_stand <- rowSums(p_tmp)

  return(p_severe_)
}
