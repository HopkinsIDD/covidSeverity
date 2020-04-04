


##' Function to load age-specific population by county for states + year of interest
##'
##' @param census_api_key census API key required for tidycensus
##' @param modeled_states states for which to pull data
##' @param census_year year for which to generate data
##' 
##' @return a data frame with columns
##'          - GEOID
##'          - agecat
##'          - estimate
##'          - NAME
##'          - age_l
##'          - age_r
##'
##' @import tidyverse
##' @import tidycensus
##' 
##' @export

load_age_county_pop <- function(census_api_key,
                                modeled_states,
                                census_year){
  
  census_api_key(census_api_key)
  
  tmp <- get_acs(geography = "county",
                variables = c(age_0004_M="B01001_003",
                              age_0509_M="B01001_004",
                              age_1014_M="B01001_005",
                              age_1517_M="B01001_006",
                              age_1819_M="B01001_007",
                              age_2020_M="B01001_008",
                              age_2121_M="B01001_009",
                              age_2224_M="B01001_010",
                              age_2529_M="B01001_011",
                              age_3034_M="B01001_012",
                              age_3539_M="B01001_013",
                              age_4044_M="B01001_014",
                              age_4549_M="B01001_015",
                              age_5054_M="B01001_016",
                              age_5559_M="B01001_017",
                              age_6061_M="B01001_018",
                              age_6264_M="B01001_019",
                              age_6566_M="B01001_020",
                              age_6769_M="B01001_021",
                              age_7074_M="B01001_022",
                              age_7579_M="B01001_023",
                              age_8084_M="B01001_024",
                              age_8599_M="B01001_025",
                              age_0004_F="B01001_027",
                              age_0509_F="B01001_028",
                              age_1014_F="B01001_029",
                              age_1517_F="B01001_030",
                              age_1819_F="B01001_031",
                              age_2020_F="B01001_032",
                              age_2121_F="B01001_033",
                              age_2224_F="B01001_034",
                              age_2529_F="B01001_035",
                              age_3034_F="B01001_036",
                              age_3539_F="B01001_037",
                              age_4044_F="B01001_038",
                              age_4549_F="B01001_039",
                              age_5054_F="B01001_040",
                              age_5559_F="B01001_041",
                              age_6061_F="B01001_042",
                              age_6264_F="B01001_043",
                              age_6566_F="B01001_044",
                              age_6769_F="B01001_045",
                              age_7074_F="B01001_046",
                              age_7579_F="B01001_047",
                              age_8084_F="B01001_048",
                              age_8599_F="B01001_049"),
                state = modeled_states,
                year = census_year,
                keep_geo_vars = TRUE) %>%
                mutate(sex = substr(variable, 10, 10),
                       agecat = substr(variable, 1, 8)) %>%
                group_by(GEOID, agecat) %>%
                summarize(estimate = sum(estimate),
                          NAME = NAME[1]) %>%
                mutate(age_l = as.numeric(substr(agecat, 5, 6)),
                       age_r = as.numeric(substr(agecat, 7, 8)))
  return(tmp)
}



##' Function to format age-specific county population into desired bins
##'
##' @param age_county_pop ACS data from load_age_county_pop
##' @param age_lower lower limit of desired age categories
##' @param age_upper upper limit of desired age categories
##' @param agecat_labels if NULL, will be paste0(age_lower, "_", age_upper)
##' 
##' @return a data frame with columns
##'          - GEOID
##'          - cat_l (index of age category)
##'          - estimate
##'          - NAME
##'          - age_l
##'          - age_r
##'          - agecat
##'
##' @import tidyverse 
##'
##' @export

format_age_county_pop <- function(age_county_pop,
                                  age_lower = seq(0, 80, 10),
                                  age_upper = c(seq(9, 79, 10), 99),
                                  agecat_labels=NULL){
  
  if(is.null(agecat_labels)){agecat_labels = paste0(age_lower, "_", age_upper)}
  
  tmp <- age_county_pop
  tmp$cat_l <- apply(tmp, 1, function(x) max(which(age_lower<=as.numeric(x['age_l']))))
  tmp$cat_r <- apply(tmp, 1, function(x) which(age_upper>=as.numeric(x['age_r']))[1])
  
  test <- (tmp$cat_l==tmp$cat_r)
  if(sum(!test)>0){ stop("age categories are not inclusive of ACS categories") }
  
  tmp <- tmp %>% 
         group_by(GEOID, cat_l) %>%
         summarize(estimate = sum(estimate),
                   NAME = NAME[1]) %>%
         mutate(age_l = age_lower[cat_l],
                age_r = age_upper[cat_l],
                agecat = agecat_labels[cat_l]) %>%
         group_by(GEOID) %>%
         mutate(tot_pop = sum(estimate)) %>%
         ungroup() %>%
         mutate(page = estimate / tot_pop)
  
  return(tmp)
}


##' Function to make population data into matrix for multiplication
##'
##' @param age_county_pop_formatted formatted ACS data from format_age_county_pop
##' 
##' @return a listof GEOID vector, agecat labels, and p_age matrix
##' 
##' @import tidyverse 
##'
##' @export
make_age_matrix <- function(age_county_pop_formatted){
  p_age <- age_county_pop_formatted %>% 
           select(GEOID, agecat, page) %>%
           pivot_wider(names_from = agecat, values_from = page)
  GEOID <- p_age[,1]
  agecat <- names(p_age)[-1]
  p_age <- as.matrix(p_age[,-1])
  return(list(GEOID=GEOID, agecat=agecat, p_age=p_age))
}

## examples
## xx <- load_age_county_pop(census_api_key = "c235e1b5620232fab506af060c5f8580604d89c1", modeled_states = state.abb, census_year = 2018)
## xx <- format_age_county_pop(xx)
## tmp <- make_age_matrix(xx)
