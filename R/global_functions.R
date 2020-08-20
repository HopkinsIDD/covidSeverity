

subset_geoid_params <- function(   .region = "Africa",
                                   .country = NULL,
                                   region_var = "region",
                                   geoid_params,
                                   geodata ){
  
  geoid_country_ <- geoid_region_ <- NULL
  
  # Filter to desired locations
  if (!is.null(.country)){
    geoid_country_ <- geodata %>% dplyr::filter(admin0 %in% .country) %>% pull(geoid)
  }
  
  if (!is.null(.region)){
    geoid_region_ <- geodata %>% 
      dplyr::mutate(reg = get(region_var)) %>%
      dplyr::filter(reg %in% .region) %>% pull(geoid)
  }
  geoids_ <- sort(unique(c(geoid_country_, geoid_region_)))
  
  geoid_params_ <- geoid_params %>% filter(geoid %in% geoids_)
  
  return(geoid_params_)
}