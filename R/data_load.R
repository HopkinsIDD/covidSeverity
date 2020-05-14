


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
  p_age_dat <- age_county_pop_formatted %>%
           select(GEOID, agecat, page) %>%
           pivot_wider(names_from = agecat,
                       values_from = page)
  GEOID <- p_age_dat[,1]
  agecat <- names(p_age_dat)[-1]
  p_age <- as.matrix(p_age_dat[,-1])
  est_age <- age_county_pop_formatted %>%
    select(GEOID, agecat, estimate) %>%
    pivot_wider(names_from = agecat,
                values_from = estimate) %>%
    select(-GEOID) %>%
    as.matrix()

  return(list(GEOID=GEOID, agecat=agecat, p_age=p_age, est_age=est_age))
}









##' Function to get population data from a worldpop .tif and geounit shapefile
##' This function extracts raster values for overlayed polygons
##' 
##' @param wp_file WorldPop population file path for the country of interest
##' @param shp Shapefile file path for shapefile of geunits (i.e., admin2)
##'
##' @return 
##'
##' @import dplyr
##' @importFrom exactextractr exact_extract
##' @importFrom raster raster
##' @importFrom sf read_sf
##' @importFrom tibble as_tibble
##' 
##' @references WorldPop (https://www.worldpop.org/geodata)
##' 
load_worldpop <- function(wp_file, shp) {
  
  #load data
  pop <- raster::raster(wp_file) #raster from world pop
  adm2 <- sf::read_sf(shp)  #district level shapefile
  
  #extract raster values by summing across 100m grids within shapefile polygons
  loc_values <- adm2 %>% 
    dplyr::mutate(sum = exactextractr::exact_extract(pop, adm2, 'sum')) %>%
    dplyr::select(ADM2_EN,sum) %>%
    tibble::as_tibble()
  
}


# # EXAMPLE
# # Download the files and save them
# 
# library(doParallel)
# library(foreach)
# library(tidyverse)
# 
# wp_file <- "raw_data/BGD/BGD_ppp_2015_adj_v2.tif"
# shp <- "raw_data/BGD/bgd_admbnda_adm2_bbs_20180410/bgd_admbnda_adm2_bbs_20180410.shp"
# 
# wp_pop_BGD <- covidSeverity:::load_worldpop(wp_file, shp)
  






##' Function to get population data from a worldpop geotiff and geounit shapefile
##' This function extracts raster values for overlayed polygons
##' 
##' @param shp Shapefile file path for shapefile of geunits (i.e., admin2)
##' @param country ISO3 of country of interest
##' @param year Year of population data (2000 to 2020)
##' @param save_dir directory where to save geotiff files
##' @param cores number of cores to parallelize over
##'
##' @return file names of the geotiffs downloaded and saved.
##'
##' @import dplyr doParallel foreach
##' @importFrom RCurl getURL
##' 
##' @references WorldPop (https://www.worldpop.org/geodata)
##'
##' @export
##' 
download_worldpop_agetifs <- function(country="BGD", year="2020", save_dir="raw_data", cores=4){
  
  dir.create(save_dir, recursive = TRUE, showWarnings = FALSE)
  
  url <- paste0("ftp://ftp.worldpop.org.uk/GIS/AgeSex_structures/Global_2000_2020/", year, "/", country, "/")
  filenames = RCurl::getURL(url, ftp.use.epsv = FALSE, dirlistonly = TRUE)
  filenames <- unlist(strsplit(filenames, "\\r\\n"))
  
  doParallel::registerDoParallel(cores)
  foreach(f=seq_len(length(filenames))) %dopar% {
    url1 <- file.path(url, filenames[f])
    download.file(url1, destfile = file.path(save_dir, country, filenames[f]), mode="wb")
  }
  doParallel::stopImplicitCluster()
  
  print(paste0("Successfully downloaded age population files from Worldpop and save to ", save_dir,"/",country))
    
  return(filenames)
  
}









##' Function to get population data from a worldpop geotiff and geounit shapefile
##' This function extracts raster values for overlayed polygons
##' 
##' @param shp Shapefile file path for shapefile of geunits (i.e., admin2)
##' @param country ISO3 of country of interest
##' @param year Year of population data (2000 to 2020)
##' @param save_dir directory where to save geotiff files
##' @param cores number of cores to parallelize over
##'
##' @return long age population data by admin level 2
##'
##' @import dplyr tidyr doParallel foreach
##' @importFrom RCurl getURL
##' @importFrom exactextractr exact_extract
##' @importFrom raster raster
##' @importFrom sf read_sf
##' @importFrom tibble as_tibble
##' 
##' @references WorldPop (https://www.worldpop.org/geodata)
##'
##' @export
##' 
load_worldpop_age <- function(shp, country="BGD", year="2020", save_dir="raw_data", cores=4) {
  
  url <- paste0("ftp://ftp.worldpop.org.uk/GIS/AgeSex_structures/Global_2000_2020/", year, "/", country, "/")
  filenames = RCurl::getURL(url, ftp.use.epsv = FALSE, dirlistonly = TRUE)
  filenames <- unlist(strsplit(filenames, "\\r\\n"))
  
  age_grps <- sort(unique(as.integer(data.frame(matrix(unlist(strsplit(filenames, "_")), ncol=4, byrow=TRUE), stringsAsFactors = FALSE)[,3])))
  age_grps_full <- paste(age_grps, c(age_grps[-1],100), sep="_")
  
  
  #load data
  adm2 <- sf::read_sf(shp)  #district level shapefile
  
  doParallel::registerDoParallel(cores)
  age_pop_data <- foreach(f=seq_len(length(filenames)), .combine=rbind,
                          .packages = c("dplyr", "tibble", 
                                        "exactextractr", "raster", "tidyr")) %dopar% {
    
    wp_file = file.path(save_dir, country, filenames[f])
    pop <- raster::raster(wp_file) #raster from world pop
    
    male_female <- unlist(strsplit(filenames[f], "_"))[2]
    age_grp_ <- as.integer(unlist(strsplit(filenames[f], "_")))[3]
    age_grps_full_ <- age_grps_full[which(age_grp_==age_grps)]
    
    #extract raster values by summing across 100m grids within shapefile polygons
    loc_values <- adm2 %>% 
      dplyr::mutate(sum = exactextractr::exact_extract(pop, adm2, 'sum')) %>%
      tibble::as_tibble() %>% 
      dplyr::select(ADM2_EN, sum) %>%
      dplyr::mutate(sex = male_female,
                    age = age_grps_full_)
      
    loc_values <- loc_values %>%
      tidyr::separate(age, into=c("age_l","age_r"), sep="_", remove=FALSE) %>%
      dplyr::mutate(age_r = as.integer(age_r) - 1)
    
    loc_values
    
  }
  doParallel::stopImplicitCluster()    
  
  
  age_pop_data <- age_pop_data %>% 
    tibble::as_tibble() %>%
    dplyr::rename(adm2 = ADM2_EN, pop = sum) %>% 
    tidyr::pivot_wider(names_from="sex", values_from = pop) %>%
    dplyr::rename(pop_m = m, pop_f = f) %>%
    dplyr::mutate(pop = pop_m + pop_f)

  return(age_pop_data)
  
}




##' Function to transform population data from a worldpop age groups to 10-year age groups
##' 
##' @param age_pop_data data pulled from worldpop using `load_worldpop_age`
##'
##' @return long age population data by admin level 2, in 10-year age groups
##'
##' @import dplyr
##' @importFrom tibble as_tibble
##' 
##' @references WorldPop (https://www.worldpop.org/geodata)
##'
##' @export
##' 
convert_wp_10yr <- function(age_pop_data){
  
  age_pop_10yr <- age_pop_data %>% 
    dplyr::mutate(age10 = floor(as.integer(age_l)/10)*10) %>%
    dplyr::group_by(age10, adm2) %>%
    dplyr::summarise(pop = sum(pop), pop_m = sum(pop_m), pop_f = sum(pop_f)) %>%
    dplyr::mutate(age_l = age10, age_r = age10+9) %>% 
    tibble::as_tibble() %>%
    dplyr::mutate(age = paste(age_l, age_r, sep="_")) %>%
    dplyr::select(-age10)
    
  return(age_pop_10yr)
  
}



 