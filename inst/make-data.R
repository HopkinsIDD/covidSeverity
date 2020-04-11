## make data
library(covidSeverity)

## make geoid data
USpop_geoid_agecat <- load_age_county_pop(census_api_key = "c235e1b5620232fab506af060c5f8580604d89c1",
                                          modeled_states = state.abb, census_year = 2018) %>%
  format_age_county_pop() %>%
  make_age_matrix()

rownames(USpop_geoid_agecat$p_age) <- USpop_geoid_agecat$GEOID[,1, drop=T]

US_age_geoid_pct <- USpop_geoid_agecat$p_age

save(US_age_geoid_pct, file="data/US_age_geoid_pct.rda")

rownames(USpop_geoid_agecat$est_age) <- USpop_geoid_agecat$GEOID[,1, drop=T]

US_age_geoid_pop <- USpop_geoid_agecat$est_age

save(US_age_geoid_pop, file="data/US_age_geoid_pop.rda")

## pull data for age-specific parameters
raw_age_estimates <- googlesheets4::read_sheet("https://docs.google.com/spreadsheets/d/1pDfQ2SkO--F2WNfjZxx6t9V7K51IZW0bJKV5cwSCZfA/edit#gid=1769840547",
                                               sheet="age risk") %>%
  filter(publicly_available==T)
1
save(raw_age_estimates, file="data/raw_age_estimates.rda")
