## make data

USpop_geoid_agecat <- load_age_county_pop(census_api_key = "c235e1b5620232fab506af060c5f8580604d89c1",
                                          modeled_states = state.abb, census_year = 2018) %>%
  format_age_county_pop() %>%
  make_age_matrix()
save(USpop_geoid_agecat, file="data/USpop_geoid_agecat.Rdata")
